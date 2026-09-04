// Ported from: DNA pVACseq script -> filter_vcf_standard_chroms()
// The original was a hand-rolled Python gzip/text loop; reimplemented
// verbatim as bin/filter_vcf_standard_chroms.py so behaviour is unchanged
// (keeps only chroms 1-22, X, Y, MT/M, with or without a "chr" prefix)
// but it now runs inside the pipeline rather than as an inline function.

process FILTER_STANDARD_CHROMS {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.gatk  // has python3; swap for a lighter python container if preferred

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}_filtered_for_vep.vcf"), emit: vcf

    script:
    """
    filter_vcf_standard_chroms.py ${vcf} ${meta.id}_filtered_for_vep.vcf
    """
}
