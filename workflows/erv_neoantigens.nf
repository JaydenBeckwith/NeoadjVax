/*
 * ERV (endogenous retrovirus) neoantigen discovery branch. New territory,
 * and: per the peptide-generation stub below: still the least
 * standardized of the three ("new territory" branches), but the calling
 * step itself is now real, not just sketched:
 *
 *   RNA FASTQ (direct, or BAM→FASTQ) ─▶ ERV_STAR_ALIGN ─▶ SAMTOOLS_COLLATE ─▶ TELESCOPE_QUANT ─▶ ERV_TO_PEPTIDE (stub)
 *
 * ERV_STAR_ALIGN runs its own dedicated STAR pass with relaxed multimapping
 * (--outFilterMultimapNmax 100 --winAnchorMultimapNmax 100): it can't reuse
 * RNA_VARIANT_CALLING's BAM, which discards multi-mapping reads at
 * --outFilterMultimapNmax 2, exactly the reads ERV/TE loci need since most
 * are repetitive. See modules/local/erv/erv_star_align.nf for the sourcing
 * of that parameter choice.
 *
 * SAMTOOLS_COLLATE re-collates STAR's BAM before Telescope: Telescope
 * requires collated (mate pairs adjacent), not coordinate-sorted, input;
 * see modules/local/erv/collate_bam.nf.
 *
 * Takes rna_bam OR rna_fastq_r1/r2 from the samplesheet, same pattern as
 * fusion_neoantigens.nf: rna_bam rows go through BAM_TO_FASTQ (reused from
 * the RNA variant-calling branch) first. rna_vcf-only rows have nothing for
 * this branch to align and are skipped.
 */

include { BAM_TO_FASTQ }    from '../modules/local/rna/bam_to_fastq'
include { ERV_STAR_ALIGN }  from '../modules/local/erv/erv_star_align'
include { SAMTOOLS_COLLATE } from '../modules/local/erv/collate_bam'
include { TELESCOPE_QUANT } from '../modules/local/erv/telescope_quant'
include { ERV_TO_PEPTIDE }  from '../modules/local/erv/erv_to_peptide'

workflow ERV_NEOANTIGENS {

    take:
    rna_samplesheet_ch   // [patient_id, timepoint, rna_bam, rna_fastq_r1, rna_fastq_r2, ...]

    main:
    if (!params.erv_annotation_gtf) {
        error "ERV_NEOANTIGENS needs --erv_annotation_gtf (a Telescope-format GTF from mlbendall/telescope_annotation_db, e.g. retro.hg38.v1: see modules/local/erv/telescope_quant.nf)"
    }
    erv_gtf = Channel.fromPath(params.erv_annotation_gtf).collect()
    star_index = Channel.fromPath(params.star_index_dir).collect()

    bam_rows_ch = rna_samplesheet_ch
        .filter { row -> row.rna_bam && !(row.rna_fastq_r1 && row.rna_fastq_r2) }
        .map { row ->
            def meta = [id: row.patient_id, timepoint: row.timepoint]
            tuple(meta, file(row.rna_bam))
        }
    fastq_rows_ch = rna_samplesheet_ch
        .filter { row -> row.rna_fastq_r1 && row.rna_fastq_r2 }
        .map { row ->
            def meta = [id: row.patient_id, timepoint: row.timepoint]
            tuple(meta, file(row.rna_fastq_r1), file(row.rna_fastq_r2))
        }

    BAM_TO_FASTQ(bam_rows_ch)
    fastq_ch = BAM_TO_FASTQ.out.fastq.mix(fastq_rows_ch)

    ERV_STAR_ALIGN(fastq_ch, star_index)
    SAMTOOLS_COLLATE(ERV_STAR_ALIGN.out.bam)
    TELESCOPE_QUANT(SAMTOOLS_COLLATE.out.bam, erv_gtf)
    peptide_ch = Channel.empty()
    if (params.erv_generate_peptides) {
        ERV_TO_PEPTIDE(TELESCOPE_QUANT.out.report)
        peptide_ch = ERV_TO_PEPTIDE.out.peptide_fasta
    }

    emit:
    report = TELESCOPE_QUANT.out.report
    peptide_fasta = peptide_ch
}
