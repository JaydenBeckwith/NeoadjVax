// Ported from: RNA_variant_pipeline.sh Step 6/12 (GATK SplitNCigarReads)
// RNA-specific step (splits reads with N CIGAR ops at spliced junctions) —
// this is why the RNA and DNA branches diverge after MarkDuplicates.

process SPLIT_NCIGAR_READS {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_medium'
    container params.containers.gatk

    input:
    tuple val(meta), path(bam)
    path fasta
    path fasta_fai
    path fasta_dict

    output:
    tuple val(meta), path("${meta.id}_${meta.timepoint}.split.bam"), emit: bam

    script:
    """
    gatk SplitNCigarReads \\
        -R ${fasta} \\
        -I ${bam} \\
        -O ${meta.id}_${meta.timepoint}.split.bam
    """
}
