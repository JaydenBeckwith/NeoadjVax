// New code: none of your Sequenza_tools scripts pick a single (purity,
// ploidy) "top solution" per sample; run-sequenza.R stops at
// sequenza.results(), which writes several candidate-solution files, not
// one answer. This step exists because HLA_LOH (workflows/hla_loh.nf)
// needs exactly one purity/ploidy value per patient.
//
// See bin/extract_sequenza_top_solution.py's own docstring for the full
// reasoning and caveat: it targets sequenza.results()' standard
// "<sample>_alternative_solutions.txt" output (documented sequenza R
// package behaviour, not something specific to your scripts) and picks
// the best-scoring row by a tolerantly-matched SLPP/LPP/score column,
// falling back to the first row with a loud warning if no such column is
// found. This is a defensible mechanical choice, not a confirmed match to
// NeoadjLOH's sequenza_top_solutions_summary.csv (which wasn't among what
// you sent): spot-check one real sample's output against this script's
// pick before trusting it in an actual HLA-LOH run.

process SEQUENZA_EXTRACT_TOP_SOLUTION {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/${meta.id}/purity_ploidy", mode: 'copy'

    input:
    tuple val(meta), val(sex), path(raw_results)

    output:
    tuple val(meta), env(PURITY), env(PLOIDY), emit: purity_ploidy
    path "${meta.id}_top_solution.tsv", emit: report

    script:
    """
    extract_sequenza_top_solution.py \\
        --sample ${meta.id} \\
        --results-dir ${raw_results} \\
        --out ${meta.id}_top_solution.tsv

    PURITY=\$(cut -f1 ${meta.id}_top_solution.tsv)
    PLOIDY=\$(cut -f2 ${meta.id}_top_solution.tsv)
    echo "[INFO] ${meta.id}: purity=\${PURITY} ploidy=\${PLOIDY}"
    """
}
