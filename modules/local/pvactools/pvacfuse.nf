// pvacfuse itself is a real, documented pvactools subcommand (same
// container/algorithm set as pvacseq) — the uncertain part of this branch
// is AGFUSION_ANNOTATE upstream, not this step.

process PVACFUSE_RUN {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_medium'
    container params.containers.pvactools
    publishDir "${params.outdir}/${meta.id}/pvacfuse/${meta.timepoint}", mode: 'copy'

    input:
    tuple val(meta), path(agfusion_dir), val(hla_alleles)

    output:
    tuple val(meta), path("pvacfuse_output"), emit: results

    script:
    def hla_str = hla_alleles instanceof List ? hla_alleles.join(',') : hla_alleles
    """
    pvacfuse run \\
        ${agfusion_dir} \\
        ${meta.id}_${meta.timepoint} \\
        ${hla_str} \\
        ${params.pvacseq_algorithms} \\
        pvacfuse_output \\
        -t ${task.cpus}
    """
}
