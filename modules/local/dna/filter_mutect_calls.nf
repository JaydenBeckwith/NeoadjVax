// Ported from: DNA pVACseq script -> run_mutect2() FilterMutectCalls call

process FILTER_MUTECT_CALLS {
    tag "${meta.id}"
    label 'process_medium'
    container params.containers.gatk
    publishDir "${params.outdir}/${meta.id}/dna_variants", mode: 'copy'

    input:
    tuple val(meta), path(raw_vcf), path(stats)
    path fasta
    path fasta_fai
    path fasta_dict

    output:
    tuple val(meta), path("${meta.id}_filtered.vcf.gz"), emit: vcf

    script:
    """
    gatk FilterMutectCalls \\
        -R ${fasta} \\
        -V ${raw_vcf} \\
        -O ${meta.id}_filtered.vcf.gz
    """

    stub:
    """
    touch ${meta.id}_filtered.vcf.gz
    """
}
