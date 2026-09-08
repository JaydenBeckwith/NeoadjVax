// Ported from: DNA pVACseq script -> align_and_index()
// Original did: bwa mem -> samtools view -> samtools sort as three separate
// docker calls. Piped together here (bwa mem | samtools sort) since that's
// standard practice and avoids writing an intermediate .sam to scratch.

process BWA_MEM_ALIGN {
    tag "${meta.id}:${sample_type}"
    label 'process_high'
    container params.containers.bwa

    input:
    tuple val(meta), val(sample_type), path(r1), path(r2)
    path fasta
    path fasta_index

    output:
    tuple val(meta), val(sample_type), path("${meta.id}.${sample_type}.sorted.bam"), emit: bam

    script:
    """
    bwa mem -t ${task.cpus} ${fasta} ${r1} ${r2} \\
        | samtools sort -@ ${task.cpus} -o ${meta.id}.${sample_type}.sorted.bam -
    """

    stub:
    """
    touch ${meta.id}.${sample_type}.sorted.bam
    """
}
