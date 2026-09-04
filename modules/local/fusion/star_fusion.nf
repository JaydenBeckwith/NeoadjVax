// New — gene fusion detection. STAR-Fusion chosen over Arriba because it
// reuses the same STAR aligner already standing up the RNA branch (one
// aligner across the pipeline instead of two), at the cost of needing a
// CTAT genome resource lib (params.ctat_resource_lib) alongside the
// existing GRCh38/GENCODE reference. Runs on the same FASTQs the RNA
// variant-calling branch already produces from BAM (or takes directly).

process STAR_FUSION {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_high'
    container params.containers.starfusion
    publishDir "${params.outdir}/${meta.id}/fusion/${meta.timepoint}", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)
    path ctat_resource_lib

    output:
    tuple val(meta), path("star_fusion_out/star-fusion.fusion_predictions.abridged.tsv"), emit: fusions
    path "star_fusion_out/*", emit: all

    script:
    """
    STAR-Fusion \\
        --genome_lib_dir ${ctat_resource_lib} \\
        --left_fq ${r1} \\
        --right_fq ${r2} \\
        --CPU ${task.cpus} \\
        --output_dir star_fusion_out
    """
}
