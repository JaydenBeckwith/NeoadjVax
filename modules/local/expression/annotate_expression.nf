// GTEx normal-tissue and TCGA cohort expression columns on each pVACseq
// aggregated report (bin/annotate_expression.py).
//
// Runs after ANNOTATE_CCF so one table carries clonality and expression
// together: those are the two axes that decide whether a well-predicted
// binder is worth putting in a vaccine. The patient's own baseline RNA TPM
// is used when there is RNA; otherwise the matched TCGA cohort's median
// stands in, which is why this process is still useful on a DNA-only run.
//
// The reference tables are built once per cohort by
// bin/build_expression_db.py from public GTEx/TCGA downloads; nothing is
// bundled with this repo.

process ANNOTATE_EXPRESSION {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.pvactools
    publishDir "${params.outdir}/${meta.id}/pvacseq_expression", mode: 'copy'

    input:
    tuple val(meta), path(reports), path(sample_tpm)
    path gtex_db
    path tcga_db

    output:
    tuple val(meta), path("*.expression.tsv"), emit: tsv

    script:
    def gtex_args = gtex_db ? "--gtex ${gtex_db} --gtex-tissue '${params.gtex_tissue}'" : ''
    def tcga_args = tcga_db ? "--tcga ${tcga_db} --tcga-cancer-type '${params.tcga_cancer_type}'" : ''
    def tpm_arg   = sample_tpm ? "--sample-tpm ${sample_tpm}" : ''
    """
    for report in ${reports}; do
        annotate_expression.py "\$report" \\
            ${gtex_args} \\
            ${tcga_args} \\
            ${tpm_arg} \\
            --quantile ${params.expression_quantile} \\
            --min-tpm ${params.expression_min_tpm} \\
            --expressed-tpm ${params.gtex_expressed_tpm} \\
            -o "\$(basename "\$report" .tsv).expression.tsv"
    done
    """

    stub:
    """
    printf 'ID\\tExpression Flag\\n' > ${meta.id}.Combined.aggregated.ccf.expression.tsv
    """
}
