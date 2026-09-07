// Ported from: RNA_variant_pipeline.sh Step 4/12 (Picard AddOrReplaceReadGroups)
// RGSM = sample basename (bn), matching the original: distinct from the
// DNA branch's tumor/normal RGSM convention since RNA samples are one BAM
// per patient-timepoint, not a tumor/normal pair.

process ADD_READ_GROUPS_RNA {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_low'
    container params.containers.picard

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}_${meta.timepoint}.rg.bam"), emit: bam

    script:
    def sample_name = "${meta.id}_${meta.timepoint}"
    """
    picard AddOrReplaceReadGroups \\
        I=${bam} \\
        O=${sample_name}.rg.bam \\
        RGID=${sample_name} \\
        RGLB=lib1 \\
        RGPL=ILLUMINA \\
        RGPU=unit1 \\
        RGSM=${sample_name} \\
        CREATE_INDEX=true \\
        VALIDATION_STRINGENCY=SILENT
    """
}
