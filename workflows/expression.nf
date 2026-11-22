/*
 * Expression context for the ranked candidates:
 *
 *   RNA BAM (baseline, week 6) ─> featureCounts ─> TPM
 *     ─> pVACseq + CCF reports ─> GTEx normal tissue, TCGA cohort,
 *        sample TPM, expression flag
 *     ─> gene panel: presentation machinery, interferon route and
 *        melanoma drivers, with copy number and somatic mutations
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
include { GENE_PANEL } from '../modules/local/expression/gene_panel'
include { TAA_REPORT } from '../modules/local/expression/taa'

workflow EXPRESSION {

    take:
    rna_ch        // [meta, [[sample_type, bam, bai], ...]] per patient, list may be empty
    reports_ch    // [meta, [report.tsv, ...]] aggregated reports to annotate
    vep_vcf_ch    // [meta, vep.vcf] somatic, for the panel's mutation column
    clinical_ch   // [meta, purity, ploidy, sex, segments]
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

    // The panel is per patient, not per candidate, so it is driven by the
    // clinical row (always present) rather than by the reports, and runs
    // even when a patient has no RNA and no segments.
    panel_file = params.gene_panel ? file(params.gene_panel, checkIfExists: true) : []
    panel_in_ch = vep_vcf_ch
        .join(clinical_ch.map { meta, purity, ploidy, sex, segments -> tuple(meta, ploidy, segments) })
        .join(baseline_tpm_ch, remainder: true)
        .filter { it[1] != null }
        .map { it[-1] == null ? it[0..-2] + [[]] : it }
    GENE_PANEL(panel_in_ch, gtex_db, tcga_db, gtf, panel_file)

    // Shared antigens are a per-patient gene-level question, like the
    // panel, so it runs off the clinical row too and needs no candidates.
    taa_ch = Channel.empty()
    if (params.taa) {
        taa_panel_file = params.taa_panel
            ? file(params.taa_panel, checkIfExists: true)
            : file("${projectDir}/assets/taa_panel.tsv", checkIfExists: true)
        taa_in_ch = clinical_ch
            .map { meta, purity, ploidy, sex, segments -> tuple(meta, true) }
            .join(baseline_tpm_ch, remainder: true)
            .filter { it[1] != null }
            .map { meta, _flag, tpm -> tuple(meta, tpm == null ? [] : tpm) }
        TAA_REPORT(taa_in_ch, taa_panel_file, gtex_db, tcga_db)
        taa_ch = TAA_REPORT.out.tsv
    }

    emit:
    annotated = ANNOTATE_EXPRESSION.out.tsv
    panel     = GENE_PANEL.out.tsv
    taa       = taa_ch
    tpm       = tpm_ch
}
