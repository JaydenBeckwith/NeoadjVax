// Ported from: RNA_variant_pipeline.sh Step 5/12 (Picard MarkDuplicates)

process MARK_DUPLICATES_RNA {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_medium'
    container params.containers.picard

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}_${meta.timepoint}.dupMarked.bam"), emit: bam
    path "${meta.id}_${meta.timepoint}.dup.metrics", emit: metrics

    script:
    """
    picard MarkDuplicates \\
        I=${bam} \\
        O=${meta.id}_${meta.timepoint}.dupMarked.bam \\
        M=${meta.id}_${meta.timepoint}.dup.metrics \\
        CREATE_INDEX=true \\
        VALIDATION_STRINGENCY=SILENT
    """

    stub:
    """
    touch ${meta.id}_${meta.timepoint}.dupMarked.bam ${meta.id}_${meta.timepoint}.dup.metrics
    """
}
