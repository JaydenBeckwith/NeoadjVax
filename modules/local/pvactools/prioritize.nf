// One ranked candidate list per MHC class, with a flag per failed
// criterion and a reason per rejected candidate (bin/prioritize_candidates.py).
//
// Runs last on the candidate tables, because it integrates what every
// earlier step contributed: binding from pVACseq, clonality from
// ANNOTATE_CCF, expression and normal-tissue context from
// ANNOTATE_EXPRESSION, and HLA-LOH from PARSE_HLA_LOH. Whichever of those
// are present is used; the rest leave their flags as missing data rather
// than as a pass.

process PRIORITIZE {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.pvactools
    publishDir "${params.outdir}/${meta.id}/prioritized", mode: 'copy'

    input:
    tuple val(meta), path(reports), path(hla_loh)

    output:
    tuple val(meta), path("*.prioritized.tsv"), emit: tsv

    script:
    def loh_arg = hla_loh ? "--hla-loh ${hla_loh}" : ''
    """
    for report in ${reports}; do
        prioritize_candidates.py "\$report" \\
            ${loh_arg} \\
            --max-percentile ${params.priority_max_percentile} \\
            --max-ic50 ${params.priority_max_ic50} \\
            --max-agretopicity ${params.priority_max_agretopicity} \\
            --min-ccf ${params.priority_min_ccf} \\
            --min-dna-vaf ${params.priority_min_dna_vaf} \\
            --min-dna-depth ${params.priority_min_dna_depth} \\
            --min-rna-depth ${params.priority_min_rna_depth} \\
            -o "\$(basename "\$report" .tsv).prioritized.tsv"
    done
    """

    stub:
    """
    printf 'ID\\tGene\\tBest Peptide\\tAllele\\tHLA LOH\\tAgretopicity\\tFlags\\tReject Reason\\tSelection\\tPriority Rank\\n' \\
        > ${meta.id}.Combined.prioritized.tsv
    """
}
