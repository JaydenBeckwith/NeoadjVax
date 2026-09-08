// Gene fusion detection, caller 1 of 2 (run alongside Arriba: see
// arriba_align.nf/arriba.nf: rather than instead of it, per
// --fusion_callers). STAR-Fusion needs its own CTAT genome resource lib
// (params.ctat_resource_lib), separate from this pipeline's own
// GRCh38/GENCODE STAR index/genome_fasta/gtf that Arriba's STAR pass
// reuses: the two tools' STAR runs are NOT shared, each needs its own
// aligner invocation tuned to what it detects fusions from.

process STAR_FUSION {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_high'
    container params.containers.starfusion
    publishDir "${params.outdir}/${meta.id}/fusion/${meta.timepoint}", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)
    path ctat_resource_lib

    output:
    tuple val(meta), val('starfusion'), path("star_fusion_out/star-fusion.fusion_predictions.abridged.tsv"), emit: fusions
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

    stub:
    """
    mkdir -p star_fusion_out
    printf '#FusionName\\tLeftGene\\tRightGene\\tLeftBreakpoint\\tRightBreakpoint\\tJunctionReadCount\\tSpanningFragCount\\tFFPM\\n' > star_fusion_out/star-fusion.fusion_predictions.abridged.tsv
    printf 'Stub run: no biological analysis performed.\\n' > star_fusion_out/stub.log
    """
}
