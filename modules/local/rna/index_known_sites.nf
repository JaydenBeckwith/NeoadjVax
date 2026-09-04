// Not present in the original script (it assumed .tbi indexes already sat
// next to the known-sites VCFs on Gadi) — added so the pipeline is
// self-contained if that's ever not true.

process TABIX_INDEX {
    tag "${vcf.simpleName}"
    label 'process_low'
    container params.containers.gatk

    input:
    path vcf

    output:
    tuple path(vcf), path("${vcf}.tbi"), emit: indexed

    script:
    """
    if [ -f "${vcf}.tbi" ]; then
        cp ${vcf}.tbi .
    else
        gatk IndexFeatureFile -I ${vcf}
    fi
    """
}
