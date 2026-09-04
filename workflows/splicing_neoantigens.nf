/*
 * Splicing neoantigen discovery branch. New territory for this repo, but
 * the RNA analysis it leans on already exists and is running in production:
 *   - SpliceAI variant annotation ([[spliceai-variant-pipeline]], Gadi PBS)
 *   - IsoformSwitchAnalyzeR isoform switching ([[neoadjuvant-splicing]])
 *
 * This workflow deliberately does NOT re-run either of those — it ingests
 * SpliceAI's output and stops at a stub for the peptide-generation step,
 * which is the one part with no settled design yet (see
 * modules/local/splicing/splice_to_peptide.nf).
 */

include { SPLICEAI_FILTER }  from '../modules/local/splicing/spliceai_filter'
include { SPLICE_TO_PEPTIDE } from '../modules/local/splicing/splice_to_peptide'

workflow SPLICING_NEOANTIGENS {

    take:
    rna_samplesheet_ch   // [patient_id, timepoint, ...]

    main:
    if (!params.spliceai_annotated_vcf_dir) {
        error "SPLICING_NEOANTIGENS needs --spliceai_annotated_vcf_dir pointing at the existing SpliceAI pipeline's output"
    }

    spliceai_vcf_ch = rna_samplesheet_ch.map { row ->
        def meta = [id: row.patient_id, timepoint: row.timepoint]
        def vcf  = file("${params.spliceai_annotated_vcf_dir}/${row.patient_id}_${row.timepoint}*.vcf*")
        tuple(meta, vcf)
    }

    SPLICEAI_FILTER(spliceai_vcf_ch)
    SPLICE_TO_PEPTIDE(SPLICEAI_FILTER.out.vcf)

    // Once SPLICE_TO_PEPTIDE produces a real peptide FASTA, feed it through
    // pVACtools' generate_protein_fasta / pvacsplice path (pVACtools >=4.0
    // has an experimental `pvacsplice` command purpose-built for this) —
    // deferred until the peptide-generation design above is settled.

    emit:
    peptide_fasta = SPLICE_TO_PEPTIDE.out.peptide_fasta
}
