// Ported from: DNA pVACseq script -> align_and_index() Picard AddOrReplaceReadGroups call
// NOTE: original hardcoded RGID=1/RGLB=lib1/RGPU=unit1 for both tumor and
// normal — kept per-sample-type here (RGID includes sample_type) so tumor
// and normal read groups don't collide once merged in downstream QC.

process ADD_READ_GROUPS {
    tag "${meta.id}:${sample_type}"
    label 'process_low'
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
        RGID=${meta.id}_${sample_type} \\
        RGLB=lib1 \\
        RGPL=illumina \\
        RGPU=unit1 \\
        RGSM=${sample_type} \\
        VALIDATION_STRINGENCY=LENIENT
    """
}
