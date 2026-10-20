// BWA-MEM straight into a coordinate sort. The read group is set here
// (SM = tumor/normal) rather than with a later AddOrReplaceReadGroups
// pass, which would rewrite the whole WGS BAM a second time. SM is the
// literal sample type because Mutect2, bam-readcount, WhatsHap and
// pVACseq all key on it: the VCF columns end up named "tumor"/"normal".

process BWA_MEM_ALIGN {
    tag "${meta.id}:${sample_type}"
    label 'process_high'
    container params.containers.bwa

    input:
    tuple val(meta), val(sample_type), path(r1), path(r2)
    path fasta
    path bwa_index

    output:
    tuple val(meta), val(sample_type), path("${meta.id}.${sample_type}.sorted.bam"), emit: bam

    script:
    def rg = "@RG\\tID:${meta.id}_${sample_type}\\tSM:${sample_type}\\tLB:${meta.id}_${sample_type}\\tPL:ILLUMINA"
    def sort_threads = Math.max(1, (task.cpus / 4) as int)
    """
    bwa mem -t ${task.cpus} -R '${rg}' ${fasta} ${r1} ${r2} \\
        | samtools sort -@ ${sort_threads} -m 2G -o ${meta.id}.${sample_type}.sorted.bam -
    """

    stub:
    """
    touch ${meta.id}.${sample_type}.sorted.bam
    """
}
