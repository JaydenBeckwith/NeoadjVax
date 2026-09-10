/*
 * RNA-seq realignment + variant calling, per patient-timepoint.
 *
 * Ported from Jayden's original RNA_variant_pipeline.sh (PBS, GATK Best
 * Practices for RNA-seq short variant discovery):
 *   BAM -> FASTQ -> STAR 2-pass -> AddRG -> MarkDup -> SplitNCigarReads
 *   -> BQSR -> HaplotypeCaller -> VariantFiltration -> VariantsToTable
 *
 * NOTE the original script's tail (a BRAF-locus allelic-count extraction
 * feeding a "PIPELINE SUMMARY" collector job) was cut off mid-heredoc in
 * what was pasted in: braf_allelic_counts.tsv is referenced by the
 * cleanup step but never actually produced. Flagged rather than guessed at;
 * see docs/ARCHITECTURE.md open-questions section.
 */

include { STAR_INDEX }                              from '../modules/local/rna/star_index'
include { SAMTOOLS_FAIDX; PICARD_CREATE_SEQUENCE_DICTIONARY } from '../modules/local/dna/prepare_genome'
include { TABIX_INDEX as TABIX_MILLS }               from '../modules/local/rna/index_known_sites'
include { TABIX_INDEX as TABIX_1000G }               from '../modules/local/rna/index_known_sites'
include { TABIX_INDEX as TABIX_DBSNP }                from '../modules/local/rna/index_known_sites'
include { BAM_TO_FASTQ }                             from '../modules/local/rna/bam_to_fastq'
include { STAR_ALIGN }                               from '../modules/local/rna/star_align'
include { ADD_READ_GROUPS_RNA }                      from '../modules/local/rna/add_read_groups_rna'
include { MARK_DUPLICATES_RNA }                      from '../modules/local/rna/mark_duplicates_rna'
include { SPLIT_NCIGAR_READS }                       from '../modules/local/rna/split_ncigar_reads'
include { BASE_RECALIBRATOR; APPLY_BQSR }            from '../modules/local/rna/base_recalibrator'
include { HAPLOTYPE_CALLER_RNA }                     from '../modules/local/rna/haplotype_caller_rna'
include { VARIANT_FILTRATION_RNA }                   from '../modules/local/rna/variant_filtration_rna'
include { VARIANTS_TO_TABLE }                        from '../modules/local/rna/variants_to_table'

workflow RNA_VARIANT_CALLING {

    take:
    rna_samplesheet_ch   // [patient_id, timepoint, rna_bam] or [..., rna_fastq_r1, rna_fastq_r2] or [..., rna_vcf]

    main:
    fasta = Channel.fromPath(params.genome_fasta).collect()
    gtf   = Channel.fromPath(params.gtf).collect()

    SAMTOOLS_FAIDX(fasta)
    PICARD_CREATE_SEQUENCE_DICTIONARY(fasta)
    // NOTE: like BWA_INDEX in the DNA branch, this always runs even when
    // --rna_skip_realignment leaves nothing that needs it: harmless.
    STAR_INDEX(fasta, gtf)

    // Shared references are reusable values, not single-use queue items.
    // Otherwise BQSR/calling stops after the first RNA sample.
    TABIX_MILLS(Channel.fromPath(params.known_mills).first())
    TABIX_1000G(Channel.fromPath(params.known_1000g).first())
    TABIX_DBSNP(Channel.fromPath(params.known_dbsnp).first())

    // patient-timepoints that already have a called+filtered RNA VCF skip
    // alignment AND calling entirely
    vcf_rows_ch = rna_samplesheet_ch
        .filter { row -> row.rna_vcf }
        .map { row -> tuple([id: row.patient_id, timepoint: row.timepoint], file(row.rna_vcf)) }

    // split the rest on whether the row supplies a BAM (realign, unless
    // --rna_skip_realignment) or FASTQs directly
    samples_ch = rna_samplesheet_ch.filter { row -> !row.rna_vcf }.map { row ->
        def meta = [id: row.patient_id, timepoint: row.timepoint]
        tuple(meta, row)
    }

    bam_rows_ch   = samples_ch.filter { meta, row -> row.rna_bam }
        .map { meta, row -> tuple(meta, file(row.rna_bam)) }
    fastq_rows_ch = samples_ch.filter { meta, row -> !row.rna_bam && row.rna_fastq_r1 }
        .map { meta, row -> tuple(meta, file(row.rna_fastq_r1), file(row.rna_fastq_r2)) }

    // --run_rna_skip_realignment / params.rna_skip_realignment (default
    // false, i.e. original behaviour: realign every rna_bam row through
    // STAR like the original PBS script always did). Set true when the
    // supplied BAMs are already aligned the way this pipeline wants
    // (STAR two-pass, --outFilterMultimapNmax 2: see star_align.nf) and
    // realigning them would just burn Gadi walltime for no benefit.
    // Mixing is fine: rna_bam rows skip straight to ADD_READ_GROUPS_RNA,
    // any rna_fastq_r1/r2 rows still align normally either way, since
    // there's no BAM to skip aligning for those.
    if (params.rna_skip_realignment) {
        prealigned_bam_ch      = bam_rows_ch
        bam_rows_to_realign_ch = Channel.empty()
    } else {
        prealigned_bam_ch      = Channel.empty()
        bam_rows_to_realign_ch = bam_rows_ch
    }

    BAM_TO_FASTQ(bam_rows_to_realign_ch)

    all_fastq_ch = BAM_TO_FASTQ.out.fastq.mix(fastq_rows_ch)

    STAR_ALIGN(all_fastq_ch, STAR_INDEX.out.index)

    aligned_bam_ch = STAR_ALIGN.out.bam.mix(prealigned_bam_ch)   // [meta, bam], either source

    ADD_READ_GROUPS_RNA(aligned_bam_ch)
    MARK_DUPLICATES_RNA(ADD_READ_GROUPS_RNA.out.bam)
    SPLIT_NCIGAR_READS(MARK_DUPLICATES_RNA.out.bam, fasta, SAMTOOLS_FAIDX.out.fai, PICARD_CREATE_SEQUENCE_DICTIONARY.out.dict)

    BASE_RECALIBRATOR(
        SPLIT_NCIGAR_READS.out.bam, fasta, SAMTOOLS_FAIDX.out.fai, PICARD_CREATE_SEQUENCE_DICTIONARY.out.dict,
        TABIX_MILLS.out.indexed.map { it[0] }, TABIX_MILLS.out.indexed.map { it[1] },
        TABIX_1000G.out.indexed.map { it[0] }, TABIX_1000G.out.indexed.map { it[1] },
        TABIX_DBSNP.out.indexed.map { it[0] }, TABIX_DBSNP.out.indexed.map { it[1] }
    )
    APPLY_BQSR(BASE_RECALIBRATOR.out.recal_table, fasta, SAMTOOLS_FAIDX.out.fai, PICARD_CREATE_SEQUENCE_DICTIONARY.out.dict)

    HAPLOTYPE_CALLER_RNA(
        APPLY_BQSR.out.bam, fasta, SAMTOOLS_FAIDX.out.fai, PICARD_CREATE_SEQUENCE_DICTIONARY.out.dict,
        TABIX_DBSNP.out.indexed.map { it[0] }, TABIX_DBSNP.out.indexed.map { it[1] }
    )
    VARIANT_FILTRATION_RNA(HAPLOTYPE_CALLER_RNA.out.vcf, fasta, SAMTOOLS_FAIDX.out.fai, PICARD_CREATE_SEQUENCE_DICTIONARY.out.dict)
    // VariantsToTable is a QC/reporting convenience, not consumed
    // downstream: skipped for rna_vcf rows since there's no reason to
    // regenerate it if a table wasn't already supplied alongside
    VARIANTS_TO_TABLE(VARIANT_FILTRATION_RNA.out.vcf)

    emit:
    filtered_vcf = VARIANT_FILTRATION_RNA.out.vcf.mix(vcf_rows_ch)   // [meta(id, timepoint), vcf]
    recal_bam    = APPLY_BQSR.out.bam               // needed downstream for RNA coverage/VAF annotation; empty for rna_vcf-only timepoints
}
