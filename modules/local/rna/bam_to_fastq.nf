// Ported from: RNA_variant_pipeline.sh Step 1-2/12
// (samtools sort -n -> Picard SamToFastq). Original checkpointed by
// checking for an existing FASTQ and skipping — Nextflow's own -resume
// cache handles that now, so the manual existence check is dropped.
// Only used when the RNA samplesheet supplies rna_bam rather than FASTQs
// directly (see RNA_VARIANT_CALLING workflow).

process BAM_TO_FASTQ {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_medium'
    container params.containers.picard

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}_${meta.timepoint}_R1.fastq.gz"), path("${meta.id}_${meta.timepoint}_R2.fastq.gz"), emit: fastq

    script:
    """
    samtools sort -n -@ ${task.cpus} -o namesorted.bam ${bam}
    picard SamToFastq \\
        I=namesorted.bam \\
        F=${meta.id}_${meta.timepoint}_R1.fastq.gz \\
        F2=${meta.id}_${meta.timepoint}_R2.fastq.gz \\
        VALIDATION_STRINGENCY=SILENT
    rm -f namesorted.bam
    """
}
