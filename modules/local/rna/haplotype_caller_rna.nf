// Ported from: RNA_variant_pipeline.sh Step 9/12 (GATK HaplotypeCaller)
// --dont-use-soft-clipped-bases is RNA-specific (avoids spurious calls at
// splice junctions): carried over verbatim.

process HAPLOTYPE_CALLER_RNA {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_high'
    container params.containers.gatk

    input:
    tuple val(meta), path(bam)
    path fasta
    path fasta_fai
    path fasta_dict
    path dbsnp
    path dbsnp_tbi

    output:
    tuple val(meta), path("${meta.id}_${meta.timepoint}.raw.vcf"), emit: vcf

    script:
    """
    gatk HaplotypeCaller \\
        -R ${fasta} \\
        -I ${bam} \\
        -O ${meta.id}_${meta.timepoint}.raw.vcf \\
        -D ${dbsnp} \\
        -stand-call-conf 20.0 \\
        --dont-use-soft-clipped-bases \\
        --native-pair-hmm-threads ${task.cpus}
    """

    stub:
    """
    touch ${meta.id}_${meta.timepoint}.raw.vcf
    """
}
