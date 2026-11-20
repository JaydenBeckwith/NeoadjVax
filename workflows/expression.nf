/*
 * Expression context for the ranked candidates:
 *
 *   RNA BAM (baseline, week 6) ─> featureCounts ─> TPM
 *     ─> pVACseq + CCF reports ─> GTEx normal tissue, TCGA cohort,
 *        sample TPM, expression flag
 *
 * Two independent halves, which is why the gates are separate. Counting
 * needs an RNA BAM and a GTF; annotation needs only a reference table, so
 * it still runs on a DNA-only patient and falls back to the TCGA cohort
 * median for expression. Every RNA timepoint given is quantified (the week
 * 6 table is published for comparison), but the candidate annotation uses
 * baseline, which is the sample the candidates were called from.
 */

include { FEATURECOUNTS; GENE_TPM } from '../modules/local/rna/quantify'
include { ANNOTATE_EXPRESSION } from '../modules/local/expression/annotate_expression'

workflow EXPRESSION {

    take:
    rna_ch        // [meta, [[sample_type, bam, bai], ...]] per patient, list may be empty
    reports_ch    // [meta, [report.tsv, ...]] aggregated reports to annotate
    gtf           // GTF for counting, or [] when params.rna_gtf is unset
    gtex_db       // GTEx table, or []
    tcga_db       // TCGA table, or []

    main:
    // Gated on the param, not on the `gtf` input: a workflow take: input is
    // a channel, so `if (gtf)` would be true even when it holds nothing.
    tpm_ch = Channel.empty()
    if (params.rna_gtf) {
        rna_bams_ch = rna_ch.flatMap { meta, rna -> rna.collect { st, bam, bai -> tuple(meta, st, bam, bai) } }
        FEATURECOUNTS(rna_bams_ch, gtf)
        GENE_TPM(FEATURECOUNTS.out.counts, gtf)
        tpm_ch = GENE_TPM.out.tpm
    }

    // Baseline is the sample the candidates were called from, so it is the
    // one whose TPM goes on the report; week 6 is published on its own.
    baseline_tpm_ch = tpm_ch
        .filter { meta, sample_type, tpm -> sample_type == 'rna_baseline' }
        .map { meta, sample_type, tpm -> tuple(meta, tpm) }

    // remainder: true keeps patients with no RNA, which then annotate from
    // the reference tables alone
    annotate_in_ch = reports_ch
        .join(baseline_tpm_ch, remainder: true)
        .filter { it[1] != null }
        .map { it[-1] == null ? it[0..-2] + [[]] : it }

    ANNOTATE_EXPRESSION(annotate_in_ch, gtex_db, tcga_db)

    emit:
    annotated = ANNOTATE_EXPRESSION.out.tsv
    tpm       = tpm_ch
}
