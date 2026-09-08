// Ported from: RNA_variant_pipeline.sh Step 11/12 (GATK VariantsToTable)
// Same field selection as the original: CHROM POS TYPE REF ALT QUAL FILTER, per-sample AD.

process VARIANTS_TO_TABLE {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_low'
    container params.containers.gatk
    publishDir "${params.outdir}/${meta.id}/rna_variants/${meta.timepoint}", mode: 'copy'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}_${meta.timepoint}.filtered.tsv"), emit: tsv

    script:
    """
    gatk VariantsToTable \\
        -V ${vcf} \\
        -F CHROM -F POS -F TYPE -F REF -F ALT -F QUAL -F FILTER -GF AD \\
        -O ${meta.id}_${meta.timepoint}.filtered.tsv
    """

    stub:
    """
    touch ${meta.id}_${meta.timepoint}.filtered.tsv
    """
}
