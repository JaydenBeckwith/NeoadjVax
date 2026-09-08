// Ported from: RNA_variant_pipeline.sh Steps 7-8/12
// (GATK BaseRecalibrator + ApplyBQSR). Known-sites VCFs are the same
// Mills/1000G/dbSNP set the original script used.

process BASE_RECALIBRATOR {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_medium'
    container params.containers.gatk

    input:
    tuple val(meta), path(bam)
    path fasta
    path fasta_fai
    path fasta_dict
    path mills
    path mills_tbi
    path kg
    path kg_tbi
    path dbsnp
    path dbsnp_tbi

    output:
    tuple val(meta), path(bam), path("${meta.id}_${meta.timepoint}.recal.data.csv"), emit: recal_table

    script:
    """
    gatk BaseRecalibrator \\
        -R ${fasta} \\
        -I ${bam} \\
        -O ${meta.id}_${meta.timepoint}.recal.data.csv \\
        --known-sites ${kg} \\
        --known-sites ${mills} \\
        --known-sites ${dbsnp}
    """

    stub:
    """
    touch ${meta.id}_${meta.timepoint}.recal.data.csv
    """
}

process APPLY_BQSR {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_medium'
    container params.containers.gatk
    publishDir "${params.outdir}/${meta.id}/rna_bam/${meta.timepoint}", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(recal_table)
    path fasta
    path fasta_fai
    path fasta_dict

    output:
    tuple val(meta), path("${meta.id}_${meta.timepoint}.recal.bam"), emit: bam

    script:
    """
    gatk ApplyBQSR \\
        -R ${fasta} \\
        -I ${bam} \\
        -O ${meta.id}_${meta.timepoint}.recal.bam \\
        --bqsr-recal-file ${recal_table}
    """

    stub:
    """
    touch ${meta.id}_${meta.timepoint}.recal.bam
    """
}
