// AGFusion 1.5.0: native caller TSV -> transcript/CDS/protein/exon annotations.
// --middlestar is essential for pVACfuse to locate the fusion junction.
process AGFUSION_ANNOTATE {
    tag "${meta.id}:${meta.timepoint}:${caller}"
    label 'process_medium'
    container params.agfusion_container
    publishDir "${params.outdir}/${meta.id}/fusion/${meta.timepoint}/${caller}", mode: 'copy'
    // PyEnsembl lazily writes indexes/pickles. Use a task-local copy, never
    // mutate a shared cache concurrently across patients or timepoints.
    stageInMode 'copy'

    input:
    tuple val(meta), val(caller), path(fusion_file, stageAs: 'caller/*')
    path database, stageAs: 'reference/*'
    path cache, stageAs: 'pyensembl_cache'

    output:
    tuple val(meta), val(caller), path('agfusion_out'), emit: annotated

    script:
    def noncanonical = params.agfusion_noncanonical ? '--noncanonical' : ''
    """
    python ${projectDir}/bin/fusion_tools.py annotate \\
        --input '${fusion_file}' --caller '${caller}' \\
        --database '${database}' --cache '${cache}' \\
        --min-reads ${params.fusion_min_read_support} --min-ffpm ${params.fusion_min_ffpm} \\
        ${noncanonical} --output agfusion_out
    """

    stub:
    """
    mkdir -p agfusion_out
    printf '{"stub": true, "caller": "${caller}", "status": "no_fusions"}\\n' > agfusion_out/annotation_qc.json
    """
}
