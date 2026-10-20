/*
 * Tumour and normal DNA from FASTQ (or an already-aligned BAM) to
 * analysis-ready BAMs:
 *
 *   FASTQ: FastQC -> Trim Galore -> FastQC -> BWA-MEM (+read group) ─┐
 *   BAM:   AddOrReplaceReadGroups ──────────────────────────────────┴─> MarkDuplicates -> BQSR
 *
 * Every BAM leaves with SM=tumor or SM=normal: Mutect2, bam-readcount,
 * WhatsHap and pVACseq all identify samples by those names.
 */

include { FASTQC as FASTQC_RAW     } from '../modules/local/dna/fastqc'
include { FASTQC as FASTQC_TRIMMED } from '../modules/local/dna/fastqc'
include { TRIM_GALORE              } from '../modules/local/dna/trim_galore'
include { BWA_MEM_ALIGN            } from '../modules/local/dna/bwa_mem_align'
include { ADD_READ_GROUPS          } from '../modules/local/dna/add_read_groups'
include { MARK_DUPLICATES          } from '../modules/local/dna/mark_duplicates'
include { BASE_RECALIBRATOR; APPLY_BQSR } from '../modules/local/dna/bqsr'

workflow PREPROCESS_DNA {

    take:
    samples_ch        // samplesheet rows (maps)
    fasta
    fasta_fai
    fasta_dict
    bwa_index
    known_sites       // [vcf.gz, ...]
    known_sites_tbi   // [vcf.gz.tbi, ...]

    main:
    // one record per patient x sample type, FASTQ or BAM as the row supplies
    inputs_ch = samples_ch.flatMap { row ->
        ['tumor', 'normal'].collect { sample_type ->
            def meta = [id: row.patient_id]
            row["${sample_type}_bam"]
                ? tuple(meta, sample_type, 'bam', [file(row["${sample_type}_bam"])])
                : tuple(meta, sample_type, 'fastq', [file(row["${sample_type}_r1"]), file(row["${sample_type}_r2"])])
        }
    }
    fastq_ch = inputs_ch.filter { it[2] == 'fastq' }.map { meta, st, kind, reads -> tuple(meta, st, reads[0], reads[1]) }
    bam_ch   = inputs_ch.filter { it[2] == 'bam' }.map { meta, st, kind, bams -> tuple(meta, st, bams[0]) }

    FASTQC_RAW(fastq_ch, 'raw')
    TRIM_GALORE(fastq_ch)
    FASTQC_TRIMMED(TRIM_GALORE.out.trimmed_reads, 'trimmed')
    BWA_MEM_ALIGN(TRIM_GALORE.out.trimmed_reads, fasta, bwa_index)

    ADD_READ_GROUPS(bam_ch)

    MARK_DUPLICATES(BWA_MEM_ALIGN.out.bam.mix(ADD_READ_GROUPS.out.bam))

    BASE_RECALIBRATOR(MARK_DUPLICATES.out.bam, fasta, fasta_fai, fasta_dict, known_sites, known_sites_tbi)
    APPLY_BQSR(BASE_RECALIBRATOR.out.table, fasta, fasta_fai, fasta_dict)

    emit:
    bam = APPLY_BQSR.out.bam   // [meta, sample_type, bam, bai]
}
