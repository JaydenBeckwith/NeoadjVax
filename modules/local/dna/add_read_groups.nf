// Only for already-aligned BAM input (tumor_bam/normal_bam columns): FASTQ
// input gets its read group from BWA_MEM_ALIGN. Forces SM to the literal
// sample type, whatever the supplied BAM used, so every downstream tool
// sees the same "tumor"/"normal" sample names.

process ADD_READ_GROUPS {
    tag "${meta.id}:${sample_type}"
    label 'process_medium'
    container params.containers.picard

    input:
    tuple val(meta), val(sample_type), path(bam)

    output:
    tuple val(meta), val(sample_type), path("${meta.id}.${sample_type}.rg.bam"), emit: bam

    script:
    """
    picard AddOrReplaceReadGroups \\
        I=${bam} \\
        O=${meta.id}.${sample_type}.rg.bam \\
        SORT_ORDER=coordinate \\
        RGID=${meta.id}_${sample_type} \\
        RGLB=${meta.id}_${sample_type} \\
        RGPL=ILLUMINA \\
        RGPU=unit1 \\
        RGSM=${sample_type} \\
        VALIDATION_STRINGENCY=LENIENT
    """

    stub:
    """
    touch ${meta.id}.${sample_type}.rg.bam
    """
}
