// Ported from: RNA_variant_pipeline.sh Step 10/12 (GATK VariantFiltration)
// window/cluster/filter expression carried over verbatim (FS > 30.0 || QD < 2.0).

process VARIANT_FILTRATION_RNA {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_medium'
    container params.containers.gatk
    publishDir "${params.outdir}/${meta.id}/rna_variants/${meta.timepoint}", mode: 'copy'

    input:
    tuple val(meta), path(vcf)
    path fasta
    path fasta_fai
    path fasta_dict

    output:
    tuple val(meta), path("${meta.id}_${meta.timepoint}.filtered.vcf"), emit: vcf

    script:
    """
    gatk VariantFiltration \\
        -R ${fasta} \\
        -V ${vcf} \\
        -O ${meta.id}_${meta.timepoint}.filtered.vcf \\
        -window 35 \\
        -cluster 3 \\
        -filter "FS > 30.0 || QD < 2.0" \\
        --filter-name FSQD
    """
}
