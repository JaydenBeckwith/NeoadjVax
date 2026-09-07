// pvacfuse itself is a real, documented pvactools subcommand (same
// container/algorithm set as pvacseq): the uncertain part of this branch
// is AGFUSION_ANNOTATE upstream, not this step. Tagged/published per
// caller (starfusion/arriba run in parallel, see fusion_neoantigens.nf)
// so results from each don't collide.

process PVACFUSE_RUN {
    tag "${meta.id}:${meta.timepoint}:${caller}"
    label 'process_medium'
    container params.containers.pvactools
    publishDir "${params.outdir}/${meta.id}/pvacfuse/${meta.timepoint}/${caller}", mode: 'copy'

    input:
    tuple val(meta), val(caller), path(agfusion_dir), val(hla_alleles)

    output:
    tuple val(meta), val(caller), path("pvacfuse_output"), emit: results

    script:
    def hla_str = hla_alleles instanceof List ? hla_alleles.join(',') : hla_alleles
    """
    pvacfuse run \\
        ${agfusion_dir} \\
        ${meta.id}_${meta.timepoint}_${caller} \\
        ${hla_str} \\
        ${params.pvacseq_algorithms} \\
        pvacfuse_output \\
        -t ${task.cpus}
    """
}
