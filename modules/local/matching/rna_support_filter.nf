// New module: implements the RNA-support check from [[neoantigen-score]]:
// "mapping variant loci to matched RNA via variant calling on the RNA
// sample to see how many carry over". Did not exist in either original
// script. Runs once per (patient, timepoint) since RNA is longitudinal.

process RNA_SUPPORT_FILTER {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_low'
    container params.containers.gatk  // has python3
    publishDir "${params.outdir}/${meta.id}/rna_supported_variants/${meta.timepoint}", mode: 'copy'

    input:
    tuple val(meta), path(dna_vcf), path(rna_vcf)

    output:
    tuple val(meta), path("${meta.id}_${meta.timepoint}.rna_supported.vcf"), emit: vcf
    path "${meta.id}_${meta.timepoint}.rna_support_summary.tsv", emit: summary

    script:
    def min_reads = params.min_rna_alt_reads ?: 1
    """
    match_rna_support.py \\
        --dna-vcf ${dna_vcf} \\
        --rna-vcf ${rna_vcf} \\
        --out-vcf ${meta.id}_${meta.timepoint}.rna_supported.vcf \\
        --min-rna-alt-reads ${min_reads} \\
        --summary ${meta.id}_${meta.timepoint}.rna_support_summary.tsv
    """

    stub:
    """
    touch ${meta.id}_${meta.timepoint}.rna_supported.vcf ${meta.id}_${meta.timepoint}.rna_support_summary.tsv
    """
}
