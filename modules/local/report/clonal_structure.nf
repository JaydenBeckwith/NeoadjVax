// Clonal clusters and week 6 persistence (bin/clonal_structure.py).
//
// The candidate tables carry a CCF per variant but no structure: whether
// this tumour is one clone or three, and whether the clones present at
// baseline are still there after treatment. For a neoadjuvant vaccine the
// second question is the point - a peptide against a clone the checkpoint
// inhibitor has already cleared is a wasted slot.
//
// Persistence comes from the RNA timepoints, because these patients have
// RNA at baseline and week 6 but tumour DNA at baseline only. A real
// tumour evolution tree needs several regions or DNA timepoints, so this
// stops at clusters and does not draw one.

process CLONAL_STRUCTURE {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.pvactools
    publishDir "${params.outdir}/${meta.id}/clonal_structure", mode: 'copy'

    input:
    tuple val(meta), path(reports)

    output:
    tuple val(meta), path("${meta.id}.clonal_structure.tsv"), emit: clusters
    tuple val(meta), path("${meta.id}.clonal_variants.tsv"), emit: variants
    tuple val(meta), path("${meta.id}.clonal_structure.p*"), optional: true, emit: figure

    script:
    def plot = params.clonal_plot ? '' : '--no-plot'
    """
    # Clonal structure is a property of variants, not of epitopes, so one
    # table is enough: the Combined report where there is one, otherwise
    # whichever single class the run produced.
    report=\$(ls ${reports} | grep -i combined | head -n1 || true)
    if [ -z "\$report" ]; then
        report=\$(ls ${reports} | head -n1)
    fi

    clonal_structure.py "\$report" \\
        --prefix ${meta.id} \\
        --gap ${params.clonal_gap} \\
        --min-cluster ${params.clonal_min_cluster} \\
        --clonal-ccf ${params.clonal_ccf} \\
        --min-rna-depth ${params.priority_min_rna_depth} \\
        ${plot}
    """

    stub:
    """
    printf 'Cluster\\tLabel\\tCCF Mean\\tVariants\\tWeek 6 Persistence\\n' > ${meta.id}.clonal_structure.tsv
    printf '1\\tclonal\\t0.92\\t3\\t1.00\\n' >> ${meta.id}.clonal_structure.tsv
    printf 'ID\\tGene\\tCCF\\tCluster\\tChange\\n' > ${meta.id}.clonal_variants.tsv
    printf 'chr1-100-100-A-T\\tGENE1\\t0.92\\t1\\tpersisted\\n' >> ${meta.id}.clonal_variants.tsv
    """
}
