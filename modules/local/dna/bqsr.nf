// GATK base quality score recalibration against the dbSNP, Mills and
// 1000G known-sites VCFs (each needs its .tbi beside it). The recalibrated
// BAMs are the analysis-ready BAMs every later step reads, so they are
// published.

process BASE_RECALIBRATOR {
    tag "${meta.id}:${sample_type}"
    label 'process_medium'
    container params.containers.gatk

    input:
    tuple val(meta), val(sample_type), path(bam), path(bai)
    path fasta
    path fasta_fai
    path fasta_dict
    path known_sites
    path known_sites_tbi

    output:
    tuple val(meta), val(sample_type), path(bam), path(bai), path("${meta.id}.${sample_type}.recal.table"), emit: table

    script:
    def known = known_sites.collect { "--known-sites ${it}" }.join(' ')
    """
    gatk --java-options "-Xmx${(task.memory.toGiga() * 0.8) as int}g" BaseRecalibrator \\
        -R ${fasta} \\
        -I ${bam} \\
        ${known} \\
        -O ${meta.id}.${sample_type}.recal.table
    """

    stub:
    """
    touch ${meta.id}.${sample_type}.recal.table
    """
}

process APPLY_BQSR {
    tag "${meta.id}:${sample_type}"
    label 'process_medium'
    container params.containers.gatk
    publishDir "${params.outdir}/${meta.id}/bam", mode: 'copy'

    input:
    tuple val(meta), val(sample_type), path(bam), path(bai), path(recal_table)
    path fasta
    path fasta_fai
    path fasta_dict

    output:
    tuple val(meta), val(sample_type), path("${meta.id}.${sample_type}.recal.bam"), path("${meta.id}.${sample_type}.recal.bam.bai"), emit: bam

    script:
    """
    gatk --java-options "-Xmx${(task.memory.toGiga() * 0.8) as int}g" ApplyBQSR \\
        -R ${fasta} \\
        -I ${bam} \\
        --bqsr-recal-file ${recal_table} \\
        -O ${meta.id}.${sample_type}.recal.bam
    mv ${meta.id}.${sample_type}.recal.bai ${meta.id}.${sample_type}.recal.bam.bai
    """

    stub:
    """
    touch ${meta.id}.${sample_type}.recal.bam ${meta.id}.${sample_type}.recal.bam.bai
    """
}
