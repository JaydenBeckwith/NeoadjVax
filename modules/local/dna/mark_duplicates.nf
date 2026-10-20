// Picard MarkDuplicates. The index is written by Picard itself (the picard
// image has no samtools) and renamed to .bam.bai, the name every later
// tool looks for first.

process MARK_DUPLICATES {
    tag "${meta.id}:${sample_type}"
    label 'process_medium'
    container params.containers.picard
    publishDir "${params.outdir}/${meta.id}/qc", mode: 'copy', pattern: '*.metrics.txt'

    input:
    tuple val(meta), val(sample_type), path(bam)

    output:
    tuple val(meta), val(sample_type), path("${meta.id}.${sample_type}.dedup.bam"), path("${meta.id}.${sample_type}.dedup.bam.bai"), emit: bam
    path "${meta.id}.${sample_type}.markdup.metrics.txt", emit: metrics

    script:
    """
    picard MarkDuplicates \\
        I=${bam} \\
        O=${meta.id}.${sample_type}.dedup.bam \\
        M=${meta.id}.${sample_type}.markdup.metrics.txt \\
        CREATE_INDEX=true \\
        VALIDATION_STRINGENCY=LENIENT
    mv ${meta.id}.${sample_type}.dedup.bai ${meta.id}.${sample_type}.dedup.bam.bai
    """

    stub:
    """
    touch ${meta.id}.${sample_type}.dedup.bam ${meta.id}.${sample_type}.dedup.bam.bai ${meta.id}.${sample_type}.markdup.metrics.txt
    """
}
