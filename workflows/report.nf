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
 * The EXPRESSION and FUSION outputs are optional: when a step did not run,
 * the report falls back to what it has and marks that section as not
 * available. A patient the fusion branch skipped (no RNA reads, no
 * STAR-Fusion table) is in the same position as a run with --fusion false.
 */

include { RNA_SUPPORT; REPORT_FIGURES; PATIENT_REPORT; REPORT_PDF } from '../modules/local/report/report'

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
    pvacfuse_ch         // [meta, pvacfuse/] or empty
    fusion_calls_ch     // [meta, star-fusion predictions.tsv] or empty
    prioritized_ch      // [meta, [prioritized.tsv, ...]] or empty
    hla_loh_ch          // [meta, hla_loh.tsv] or empty
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

    // The report draws nothing itself: every figure is made here first and
    // handed over as a file. Anchored on the clinical row so a patient with
    // neither candidates nor an LOH table still reaches PATIENT_REPORT.
    figures_ch = Channel.empty()
    if (params.report_figures) {
        REPORT_FIGURES(clinical_ch
            .map { meta, purity, ploidy, sex, segments -> tuple(meta, []) }
            .join(prioritized_ch, remainder: true)
            .filter { it[1] != null }
            .join(hla_loh_ch, remainder: true)
            .filter { it[1] != null }
            .map { meta, _anchor, prioritized, loh ->
                tuple(meta, prioritized == null ? [] : prioritized, loh == null ? [] : loh) })
        figures_ch = REPORT_FIGURES.out.figures.map { meta, pdfs, pngs -> tuple(meta, pdfs) }
    }

    // Anchored on the clinical row, which every patient has, so a run
    // without the somatic analysis still produces a report: the sections
    // that need variants then say they were not assessed. Optional inputs
    // arrive as null after a remainder join; [] is "no file". Each one is
    // wrapped in a single-element list first, so an unmatched remainder is
    // one null rather than a variable number of them.
    report_in_ch = clinical_ch
        .join(hla_ch)
        .join(vcf_ch.map { meta, vcf -> tuple(meta, [vcf]) }, remainder: true)
        .filter { it[1] != null }
        .join(expression_ch.map { meta, reports -> tuple(meta, [reports]) }, remainder: true)
        .filter { it[1] != null }
        .join(ccf_ch.map { meta, reports -> tuple(meta, [reports]) }, remainder: true)
        .filter { it[1] != null }
        .join(panel_ch, remainder: true)
        .filter { it[1] != null }
        .join(RNA_SUPPORT.out.summary.map { meta, summary, figure -> tuple(meta, [summary, figure]) }, remainder: true)
        .filter { it[1] != null }
        .join(timepoints_ch)
        .join(pvacfuse_ch.map { meta, dir -> tuple(meta, [dir]) }, remainder: true)
        .filter { it[1] != null }
        .join(fusion_calls_ch.map { meta, calls -> tuple(meta, [calls]) }, remainder: true)
        .filter { it[1] != null }
        .join(figures_ch.map { meta, pdfs -> tuple(meta, [pdfs]) }, remainder: true)
        .filter { it[1] != null }
        .join(hla_loh_ch.map { meta, loh -> tuple(meta, [loh]) }, remainder: true)
        .filter { it[1] != null }
        .map { meta, purity, ploidy, sex, segments, hla, vcf, expression, ccf, panel, support, tps,
               fusion, calls, figures, loh ->
            tuple(meta, vcf ? vcf[0] : [], hla, purity, ploidy, sex, segments,
                  (expression ?: ccf)?.getAt(0) ?: [], panel ?: [],
                  support ? support[0] : [], support ? support[1] : [], tps,
                  fusion ? fusion[0] : [], calls ? calls[0] : [],
                  figures ? figures[0] : [], loh ? loh[0] : []) }
    PATIENT_REPORT(report_in_ch, fasta, fasta_fai, target_intervals)

    pdf_ch = Channel.empty()
    if (params.report_pdf) {
        REPORT_PDF(PATIENT_REPORT.out.tex
            .join(RNA_SUPPORT.out.summary.map { meta, summary, figure -> tuple(meta, figure) }, remainder: true)
            .map { meta, tex, figure -> tuple(meta, tex, figure ?: []) }
            .join(figures_ch.map { meta, pdfs -> tuple(meta, [pdfs]) }, remainder: true)
            .filter { it[1] != null }
            .map { meta, tex, figure, figures -> tuple(meta, tex, figure, figures ? figures[0] : []) })
        pdf_ch = REPORT_PDF.out.pdf
    }

    emit:
    rna_support = RNA_SUPPORT.out.summary
    tex         = PATIENT_REPORT.out.tex
    pdf         = pdf_ch
}
