// Per-patient read-outs at the end of the run: the RNA-support figure and
// the melanoma report (bin/rna_support.py, bin/patient_report.py), and the
// report compiled to PDF.

// Non-synonymous variants and clonal neoantigens supported in each RNA
// timepoint (expressed both / week 6 only / baseline only / DNA only).
// Only patients with at least one RNA BAM reach this process.
process RNA_SUPPORT {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.plots
    publishDir "${params.outdir}/${meta.id}/report", mode: 'copy'

    input:
    tuple val(meta), path(readcounts_vcf), path(week6_vcf), val(timepoints), path(reports),
          val(purity), val(ploidy), val(sex), path(segments)

    output:
    tuple val(meta), path("${meta.id}.rna_support.summary.tsv"), path("${meta.id}.rna_support.pdf"), emit: summary
    tuple val(meta), path("${meta.id}.rna_support.variants.tsv"), emit: variants
    tuple val(meta), path("${meta.id}.rna_support.png"), emit: png

    script:
    // baseline counts live in the readcounts VCF pVACseq read; week 6 in its own copy
    def tp_args = timepoints.collect { st ->
        st == 'rna_baseline' ? "--timepoint 'Baseline=${readcounts_vcf}'" : "--timepoint 'Week 6=${week6_vcf}'"
    }.join(' ')
    def ploidy_arg   = ploidy   ? "--ploidy ${ploidy}" : ''
    def sex_arg      = sex      ? "--sex ${sex}" : ''
    def segments_arg = segments ? "--segments ${segments}" : ''
    """
    class_i=\$(for f in *.MHC_I.aggregated.ccf.tsv; do [ -e "\$f" ] && echo "\$f"; done | head -1)
    class_ii=\$(for f in *.MHC_II.aggregated.ccf.tsv; do [ -e "\$f" ] && echo "\$f"; done | head -1)
    rna_support.py \\
        --patient ${meta.id} \\
        --vcf ${readcounts_vcf} \\
        ${tp_args} \\
        \${class_i:+--class-i \$class_i} \\
        \${class_ii:+--class-ii \$class_ii} \\
        --purity ${purity} ${ploidy_arg} ${sex_arg} ${segments_arg} \\
        --min-alt-reads ${params.rna_support_min_alt_reads} \\
        --min-rna-vaf ${params.rna_support_min_vaf} \\
        --clonal-ccf ${params.clonal_ccf} \\
        --class-i-ic50 ${params.neoag_class_i_ic50} --class-i-percentile ${params.neoag_class_i_percentile} \\
        --class-ii-ic50 ${params.neoag_class_ii_ic50} --class-ii-percentile ${params.neoag_class_ii_percentile} \\
        -o ${meta.id}
    """

    stub:
    """
    printf 'Set\\tLabel\\tTotal\\n' > ${meta.id}.rna_support.summary.tsv
    printf 'ID\\tRNA Support\\n' > ${meta.id}.rna_support.variants.tsv
    touch ${meta.id}.rna_support.pdf ${meta.id}.rna_support.png
    """
}

// The report's prediction figures (bin/report_figures.py): the candidate
// landscape, the flag summary and the HLA-LOH overview.
//
// Separate from PATIENT_REPORT on purpose. The report writes LaTeX and
// draws nothing, so a figure can be regenerated from an existing results
// directory and looked at on its own, without pdflatex or the rest of the
// report, and the data behind each one is testable without matplotlib.
process REPORT_FIGURES {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.plots
    publishDir "${params.outdir}/${meta.id}/report", mode: 'copy'

    input:
    tuple val(meta), path(prioritized, stageAs: 'prioritized/*'), path(hla_loh)

    output:
    tuple val(meta), path("${meta.id}.*.pdf"), path("${meta.id}.*.png"), optional: true, emit: figures

    script:
    def loh_arg = hla_loh ? "--hla-loh ${hla_loh}" : ''
    """
    # Candidate figures come from the Combined table where there is one
    # (both MHC classes), otherwise whichever single class the run made.
    report=\$(ls prioritized/* 2>/dev/null | grep -i combined | head -n1 || true)
    if [ -z "\$report" ]; then
        report=\$(ls prioritized/* 2>/dev/null | head -n1 || true)
    fi

    report_figures.py \\
        --prefix ${meta.id} \\
        \${report:+--prioritized "\$report"} \\
        ${loh_arg} \\
        --top-labels ${params.report_figure_labels}
    """

    stub:
    """
    # Only the figures the inputs would actually have produced, so the
    # smoke test checks the wiring rather than the stub's optimism.
    if ls prioritized/* >/dev/null 2>&1; then
        touch ${meta.id}.candidates.pdf ${meta.id}.candidates.png
        touch ${meta.id}.flags.pdf ${meta.id}.flags.png
    fi
    if [ -n "${hla_loh}" ]; then
        touch ${meta.id}.hla_loh.pdf ${meta.id}.hla_loh.png
    fi
    """
}

// LaTeX report from everything the run produced for this patient. Each
// optional input arrives as [] when absent and its section says so.
process PATIENT_REPORT {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.pvactools
    publishDir "${params.outdir}/${meta.id}/report", mode: 'copy'

    input:
    tuple val(meta), path(readcounts_vcf), val(hla_alleles), val(purity), val(ploidy), val(sex), path(segments),
          path(reports, stageAs: 'reports/*'), path(panel), path(rna_summary), path(rna_figure), val(timepoints),
          path(pvacfuse_dir), path(fusion_predictions), path(report_figures), path(hla_loh)
    path fasta
    path fasta_fai
    path target_intervals

    output:
    tuple val(meta), path("${meta.id}.report.tex"), emit: tex

    script:
    // No VCF means the run did not include somatic calling; the variant,
    // TMB, UV and clonality sections then say they were not assessed.
    def vcf_arg      = readcounts_vcf ? "--vcf ${readcounts_vcf}" : ''
    def loh_arg      = hla_loh ? "--hla-loh ${hla_loh}" : ''
    def ploidy_arg   = ploidy      ? "--ploidy ${ploidy}" : ''
    def sex_arg      = sex         ? "--sex ${sex}" : ''
    def segments_arg = segments    ? "--segments ${segments}" : ''
    def panel_arg    = panel       ? "--gene-panel ${panel}" : ''
    def summary_arg  = rna_summary ? "--rna-support-summary ${rna_summary}" : ''
    def figure_arg   = rna_figure  ? "--rna-support-figure ${rna_figure}" : ''
    // absent means the fusion branch did not run for this patient
    def fusion_args  = (pvacfuse_dir ? "--pvacfuse-dir ${pvacfuse_dir}" : '') +
                       (fusion_predictions ? " --fusion-predictions ${fusion_predictions}" : '')
    def tp_args      = timepoints.collect { it == 'rna_baseline' ? "--timepoint-label Baseline" : "--timepoint-label 'Week 6'" }.join(' ')
    // TMB denominator: explicit Mb, else the capture targets, else the called contigs (WGS)
    def tmb_arg      = params.tmb_callable_mb ? "--callable-mb ${params.tmb_callable_mb}"
                     : (target_intervals ? "--target-bed ${target_intervals}" : "--contigs '${params.scatter_chromosomes}'")
    """
    # prefer the expression-annotated reports, which carry the CCF columns too
    pick() {
        for f in reports/*.\$1.aggregated.ccf.expression.tsv reports/*.\$1.aggregated.ccf.tsv; do
            if [ -e "\$f" ]; then echo "\$f"; return; fi
        done
    }
    class_i=\$(pick MHC_I)
    class_ii=\$(pick MHC_II)

    # Figures are drawn by REPORT_FIGURES and staged here beside the .tex,
    # which is where \\includegraphics looks for them. Each is optional:
    # a run without candidates has no candidate figure to embed.
    fig() { [ -e "${meta.id}.\$1.pdf" ] && echo "${meta.id}.\$1.pdf"; }
    candidate_figure=\$(fig candidates)
    flag_figure=\$(fig flags)
    loh_figure=\$(fig hla_loh)

    patient_report.py \\
        --patient ${meta.id} \\
        ${vcf_arg} \\
        ${loh_arg} \\
        \${candidate_figure:+--candidate-figure \$candidate_figure} \\
        \${flag_figure:+--flag-figure \$flag_figure} \\
        \${loh_figure:+--hla-loh-figure \$loh_figure} \\
        --hla '${hla_alleles}' \\
        --purity ${purity} ${ploidy_arg} ${sex_arg} ${segments_arg} \\
        --fasta ${fasta} --fai ${fasta_fai} \\
        ${tmb_arg} \\
        ${panel_arg} \\
        \${class_i:+--class-i \$class_i} \\
        \${class_ii:+--class-ii \$class_ii} \\
        ${summary_arg} ${figure_arg} ${tp_args} \\
        ${fusion_args} \\
        --project '${params.report_project}' \\
        --diagnosis '${params.report_diagnosis}' \\
        --tcga-cancer-type '${params.tcga_cancer_type}' \\
        --gtex-tissue '${params.gtex_tissue}' \\
        --pipeline-version '${workflow.manifest.version}' \\
        --clonal-ccf ${params.clonal_ccf} \\
        --min-alt-reads ${params.rna_support_min_alt_reads} \\
        --class-i-ic50 ${params.neoag_class_i_ic50} --class-i-percentile ${params.neoag_class_i_percentile} \\
        --class-ii-ic50 ${params.neoag_class_ii_ic50} --class-ii-percentile ${params.neoag_class_ii_percentile} \\
        --top-n ${params.report_top_n} \\
        -o ${meta.id}.report.tex
    """

    stub:
    """
    printf '\\\\documentclass{article}\\\\begin{document}stub\\\\end{document}\\n' > ${meta.id}.report.tex
    """
}

// pdflatex twice (longtable widths settle on the second pass). The figure
// is [] for a patient without RNA, whose report has no figure.
process REPORT_PDF {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.texlive
    publishDir "${params.outdir}/${meta.id}/report", mode: 'copy'

    input:
    // Every figure the .tex references has to be staged beside it, or
    // pdflatex stops on the missing file.
    tuple val(meta), path(tex), path(figure), path(report_figures)

    output:
    tuple val(meta), path("${meta.id}.report.pdf"), emit: pdf

    script:
    """
    pdflatex -interaction=nonstopmode -halt-on-error ${tex} > /dev/null
    pdflatex -interaction=nonstopmode -halt-on-error ${tex}
    """

    stub:
    """
    touch ${meta.id}.report.pdf
    """
}
