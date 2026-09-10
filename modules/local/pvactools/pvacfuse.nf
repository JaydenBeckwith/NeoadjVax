// Pinned separately from pVACseq so upgrading this branch cannot change core results.
process PVACFUSE_RUN {
    tag "${meta.id}:${meta.timepoint}:${caller}"
    label 'process_medium'
    cpus params.pvacfuse_threads
    container params.containers.pvacfuse
    publishDir "${params.outdir}/${meta.id}/pvacfuse/${meta.timepoint}/${caller}", mode: 'copy'

    input:
    tuple val(meta), val(caller), path(agfusion_dir), val(hla_alleles)
    path models, stageAs: 'mhcflurry_models'
    path iedb, stageAs: 'iedb'

    output:
    tuple val(meta), val(caller), path('pvacfuse_output'), emit: results

    script:
    def model_opt = models ? "--models '${models}'" : ''
    def iedb_opt = iedb ? "--iedb '${iedb}'" : ''
    """
    python '${projectDir}/bin/fusion_tools.py' predict \\
        --input '${agfusion_dir}' --sample '${meta.id}_${meta.timepoint}_${caller}' \\
        --caller '${caller}' --alleles '${hla_alleles}' \\
        --algorithms '${params.pvacfuse_algorithms}' --threads ${task.cpus} \\
        --lengths-i '${params.pvacfuse_epitope_lengths_i}' --lengths-ii '${params.pvacfuse_epitope_lengths_ii}' \\
        --min-reads ${params.fusion_min_read_support} --min-ffpm ${params.fusion_min_ffpm} \\
        --binding-threshold ${params.pvacfuse_binding_threshold} \\
        ${model_opt} ${iedb_opt} --output pvacfuse_output
    """

    stub:
    """
    mkdir -p pvacfuse_output
    printf '{"stub": true, "caller": "${caller}", "status": "no_predictable_fusion_proteins"}\\n' > pvacfuse_output/run_status.json
    """
}
