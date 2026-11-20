/*
 * Per-patient read-outs:
 *
 *   readcounts VCF (baseline RNA) + week 6 RNA VCF + pVACseq CCF reports
 *     ─> RNA_SUPPORT     non-synonymous variants and clonal neoantigens seen
 *                        in RNA at each timepoint (patients with RNA only)
 *   everything above + HLA + gene panel + expression-annotated reports
 *     ─> PATIENT_REPORT  melanoma report (LaTeX)
 *     ─> REPORT_PDF      pdflatex (unless --report_pdf false)
 *
 * The EXPRESSION outputs are optional: when that step did not run, the
 * report falls back to the CCF reports and marks the panel and GTEx/TCGA
 * sections as not available.
 */

include { RNA_SUPPORT; PATIENT_REPORT; REPORT_PDF } from '../modules/local/report/report'

workflow REPORT {

    take:
    readcounts_vcf_ch   // [meta, vcf.gz, tbi]  tumour DNA + baseline RNA counts
    week6_vcf_ch        // [meta, vcf]          week 6 RNA counts (patients with week 6 only)
    rna_ch              // [meta, [[sample_type, bam, bai], ...]]
    ccf_ch              // [meta, [report.ccf.tsv, ...]]
    expression_ch       // [meta, [report.ccf.expression.tsv, ...]] or empty
    panel_ch            // [meta, gene_panel.tsv] or empty
    hla_ch              // [meta, allele string]
    clinical_ch         // [meta, purity, ploidy, sex, segments]
    fasta
    fasta_fai
    target_intervals    // capture BED/interval_list, or []

    main:
    // [meta, ['rna_baseline', 'rna_week6']] in time order, possibly empty
    timepoints_ch = rna_ch.map { meta, rna -> tuple(meta, rna.collect { it[0] }.sort()) }

    vcf_ch = readcounts_vcf_ch.map { meta, vcf, tbi -> tuple(meta, vcf) }

    support_in_ch = vcf_ch
        .join(week6_vcf_ch, remainder: true)
        .filter { it[1] != null }
        .map { meta, vcf, week6 -> tuple(meta, vcf, week6 ?: []) }
        .join(timepoints_ch)
        .filter { meta, vcf, week6, tps -> tps }
        .join(ccf_ch)
        .join(clinical_ch)
    RNA_SUPPORT(support_in_ch)

    // Optional inputs arrive as null after a remainder join; [] is "no file".
    report_in_ch = vcf_ch
        .join(hla_ch)
        .join(clinical_ch)
        .join(expression_ch, remainder: true)
        .filter { it[1] != null }
        .join(ccf_ch)
        .map { meta, vcf, hla, purity, ploidy, sex, segments, expression, ccf ->
            tuple(meta, vcf, hla, purity, ploidy, sex, segments, expression ?: ccf) }
        .join(panel_ch, remainder: true)
        .filter { it[1] != null }
        // one value per side, so an unmatched remainder is a single null
        .join(RNA_SUPPORT.out.summary.map { meta, summary, figure -> tuple(meta, [summary, figure]) }, remainder: true)
        .filter { it[1] != null }
        .join(timepoints_ch)
        .map { meta, vcf, hla, purity, ploidy, sex, segments, reports, panel, support, tps ->
            tuple(meta, vcf, hla, purity, ploidy, sex, segments, reports, panel ?: [],
                  support ? support[0] : [], support ? support[1] : [], tps) }
    PATIENT_REPORT(report_in_ch, fasta, fasta_fai, target_intervals)

    pdf_ch = Channel.empty()
    if (params.report_pdf) {
        REPORT_PDF(PATIENT_REPORT.out.tex
            .join(RNA_SUPPORT.out.summary.map { meta, summary, figure -> tuple(meta, figure) }, remainder: true)
            .map { meta, tex, figure -> tuple(meta, tex, figure ?: []) })
        pdf_ch = REPORT_PDF.out.pdf
    }

    emit:
    rna_support = RNA_SUPPORT.out.summary
    tex         = PATIENT_REPORT.out.tex
    pdf         = pdf_ch
}
