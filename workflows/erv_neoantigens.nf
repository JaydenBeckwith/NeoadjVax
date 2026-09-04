/*
 * ERV (endogenous retrovirus) neoantigen discovery branch. New territory,
 * and the least standardized of the three — see
 * modules/local/erv/erv_to_peptide.nf for the open design decision.
 *
 *   RNA BAM (multi-mapping reads retained) ─▶ TELESCOPE_QUANT ─▶ ERV_TO_PEPTIDE (stub)
 *
 * NOTE: needs its own STAR alignment pass, not RNA_VARIANT_CALLING's BAM —
 * see modules/local/erv/telescope_quant.nf for why (multimapper filtering).
 * That re-alignment step isn't built yet; this workflow currently expects
 * the samplesheet to point at an already-suitable BAM via rna_bam (relaxed
 * multimapping) until we decide whether to add a dedicated STAR pass here.
 */

include { TELESCOPE_QUANT } from '../modules/local/erv/telescope_quant'
include { ERV_TO_PEPTIDE }  from '../modules/local/erv/erv_to_peptide'

workflow ERV_NEOANTIGENS {

    take:
    rna_samplesheet_ch   // [patient_id, timepoint, rna_bam, ...]

    main:
    if (!params.erv_annotation_gtf) {
        error "ERV_NEOANTIGENS needs --erv_annotation_gtf (Telescope-formatted retro.gtf)"
    }
    erv_gtf = Channel.fromPath(params.erv_annotation_gtf).collect()

    bam_ch = rna_samplesheet_ch
        .filter { row -> row.rna_bam }
        .map { row ->
            def meta = [id: row.patient_id, timepoint: row.timepoint]
            def bam  = file(row.rna_bam)
            def bai  = file("${row.rna_bam}.bai")
            tuple(meta, bam, bai)
        }

    TELESCOPE_QUANT(bam_ch, erv_gtf)
    ERV_TO_PEPTIDE(TELESCOPE_QUANT.out.report)

    emit:
    peptide_fasta = ERV_TO_PEPTIDE.out.peptide_fasta
}
