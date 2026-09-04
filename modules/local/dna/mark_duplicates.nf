// Ported from: DNA pVACseq script -> align_and_index() Picard MarkDuplicates + samtools index

process MARK_DUPLICATES {
    tag "${meta.id}:${sample_type}"
    label 'process_medium'
    container params.containers.picard
    publishDir "${params.outdir}/${meta.id}/dna_bam", mode: 'copy', pattern: '*.dedup.bam*'

    input:
    tuple val(meta), val(sample_type), path(bam)

    output:
    tuple val(meta), val(sample_type), path("${meta.id}.${sample_type}.dedup.bam"), path("${meta.id}.${sample_type}.dedup.bam.bai"), emit: bam
    path "${meta.id}.${sample_type}.metrics.txt", emit: metrics

    script:
    """
    picard MarkDuplicates \\
        I=${bam} \\
        O=${meta.id}.${sample_type}.dedup.bam \\
        M=${meta.id}.${sample_type}.metrics.txt \\
        CREATE_INDEX=false \\
        VALIDATION_STRINGENCY=LENIENT

    samtools index ${meta.id}.${sample_type}.dedup.bam
    """
}
