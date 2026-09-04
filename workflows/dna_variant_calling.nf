/*
 * DNA somatic variant calling (WGS/WES, tumor vs matched normal).
 *
 * Ported from Jayden's original Docker-based Python script:
 *   trim_and_qc -> align_and_index -> run_mutect2 -> filter_vcf_standard_chroms
 * VEP annotation + pVACseq moved to workflows/pvacseq_core.nf, since in the
 * per-patient neoantigen score design (see [[neoantigen-score]]) annotation
 * happens *after* the RNA-support check, on the smaller RNA-supported
 * variant subset, not on every raw DNA call.
 */

include { BWA_INDEX }                            from '../modules/local/dna/bwa_index'
include { SAMTOOLS_FAIDX; PICARD_CREATE_SEQUENCE_DICTIONARY } from '../modules/local/dna/prepare_genome'
include { FASTQC as FASTQC_PRE }                 from '../modules/local/dna/fastqc'
include { FASTQC as FASTQC_POST }                from '../modules/local/dna/fastqc'
include { TRIM_GALORE }                          from '../modules/local/dna/trim_galore'
include { BWA_MEM_ALIGN }                        from '../modules/local/dna/bwa_mem_align'
include { ADD_READ_GROUPS }                      from '../modules/local/dna/add_read_groups'
include { MARK_DUPLICATES }                      from '../modules/local/dna/mark_duplicates'
include { MUTECT2 }                              from '../modules/local/dna/mutect2'
include { FILTER_MUTECT_CALLS }                  from '../modules/local/dna/filter_mutect_calls'
include { FILTER_STANDARD_CHROMS }               from '../modules/local/dna/filter_standard_chroms'

workflow DNA_VARIANT_CALLING {

    take:
    dna_samplesheet_ch   // [patient_id, tumor_r1, tumor_r2, normal_r1, normal_r2, tumor_bam, normal_bam, hla_alleles]

    main:
    fasta = Channel.fromPath(params.genome_fasta).collect()

    // NOTE: BWA_INDEX always runs even for an all-BAM cohort (Nextflow
    // can't tell in advance that no row will need it) — harmless, just a
    // few wasted minutes on Gadi. Not worth the extra complexity to gate it.
    BWA_INDEX(fasta)
    SAMTOOLS_FAIDX(fasta)
    PICARD_CREATE_SEQUENCE_DICTIONARY(fasta)

    // fan each sample row out into a (meta, sample_type, ...) record per
    // tumor/normal, mirroring the original script's two ThreadPoolExecutor
    // branches — and, per-sample_type, take whichever of a FASTQ pair or an
    // already-aligned BAM the row supplies (tumor_bam/normal_bam columns;
    // see assets/samplesheet_schema.md). A row can mix both, e.g. an
    // already-aligned normal but a tumor still needing alignment.
    reads_ch = dna_samplesheet_ch.flatMap { row ->
        def meta = [id: row.patient_id]
        ['tumor', 'normal'].collect { sample_type ->
            def bam_path = row["${sample_type}_bam"]
            def r1_path  = row["${sample_type}_r1"]
            def r2_path  = row["${sample_type}_r2"]
            if (bam_path) {
                tuple(meta, sample_type, 'bam', file(bam_path))
            } else if (r1_path && r2_path) {
                tuple(meta, sample_type, 'fastq', file(r1_path), file(r2_path))
            } else {
                error "DNA_VARIANT_CALLING: patient ${row.patient_id} is missing both ${sample_type}_bam and a ${sample_type}_r1/${sample_type}_r2 pair"
            }
        }
    }

    fastq_rows_ch = reads_ch.filter { it[2] == 'fastq' }.map { meta, st, kind, r1, r2 -> tuple(meta, st, r1, r2) }
    // already-aligned BAMs skip QC/trim/align entirely and rejoin the
    // pipeline right where BWA_MEM_ALIGN's output would have been
    bam_rows_ch   = reads_ch.filter { it[2] == 'bam' }.map { meta, st, kind, bam -> tuple(meta, st, bam) }

    FASTQC_PRE(fastq_rows_ch, 'pre_trim')
    TRIM_GALORE(fastq_rows_ch)
    FASTQC_POST(TRIM_GALORE.out.trimmed_reads, 'post_trim')
    BWA_MEM_ALIGN(TRIM_GALORE.out.trimmed_reads, BWA_INDEX.out.fasta, BWA_INDEX.out.index)

    aligned_ch = BWA_MEM_ALIGN.out.bam.mix(bam_rows_ch)   // [meta, sample_type, bam], either source

    ADD_READ_GROUPS(aligned_ch)
    MARK_DUPLICATES(ADD_READ_GROUPS.out.bam)

    // regroup tumor+normal dedup BAMs back onto one row per patient for Mutect2
    tn_pairs_ch = MARK_DUPLICATES.out.bam
        .map { meta, sample_type, bam, bai -> tuple(meta.id, meta, sample_type, bam, bai) }
        .groupTuple(by: 0)
        .map { patient_id, metas, sample_types, bams, bais ->
            def idx_t = sample_types.indexOf('tumor')
            def idx_n = sample_types.indexOf('normal')
            if (idx_t < 0 || idx_n < 0) {
                error "DNA_VARIANT_CALLING: patient ${patient_id} is missing a tumor or normal BAM after dedup"
            }
            tuple(metas[idx_t], bams[idx_t], bais[idx_t], bams[idx_n], bais[idx_n])
        }

    MUTECT2(tn_pairs_ch, fasta, SAMTOOLS_FAIDX.out.fai, PICARD_CREATE_SEQUENCE_DICTIONARY.out.dict)
    FILTER_MUTECT_CALLS(MUTECT2.out.vcf, fasta, SAMTOOLS_FAIDX.out.fai, PICARD_CREATE_SEQUENCE_DICTIONARY.out.dict)
    FILTER_STANDARD_CHROMS(FILTER_MUTECT_CALLS.out.vcf)

    emit:
    filtered_vcf = FILTER_STANDARD_CHROMS.out.vcf   // [meta, vcf]
    tumor_bam    = MARK_DUPLICATES.out.bam.filter { it[1] == 'tumor' }.map { meta, st, bam, bai -> tuple(meta, bam, bai) }
    // normal (germline) BAM, not tumor — HLA typing runs on this to avoid
    // bias from tumor HLA-LOH, which [[loh-analysis]] is tracking separately
    normal_bam   = MARK_DUPLICATES.out.bam.filter { it[1] == 'normal' }.map { meta, st, bam, bai -> tuple(meta, bam, bai) }
}
