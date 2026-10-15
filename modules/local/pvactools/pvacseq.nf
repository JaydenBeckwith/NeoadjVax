// Ported from: DNA pVACseq script -> run_pvacseq()
// Consumes a VEP-annotated VCF (Wildtype + Frameshift plugins required)
// plus typed/manual HLA alleles. Shared by the core DNA+RNA-matched branch
// and, longer-term, any of the splicing/fusion/ERV branches that end up
// producing a VEP-style annotated VCF rather than needing pVACfuse.

process PVACSEQ_RUN {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_medium'
    container params.containers.pvactools
    publishDir "${params.outdir}/${meta.id}/pvacseq/${meta.timepoint}", mode: 'copy'

    input:
    tuple val(meta), path(annotated_vcf), val(hla_alleles)

    output:
    tuple val(meta), path("pvacseq_output"), emit: results

    // pvacseq's sample_name must be a VCF sample column, not the patient id:
    // ADD_READ_GROUPS sets RGSM=tumor/normal, so Mutect2's VCF columns are
    // literally "tumor" and "normal". --pass-only drops FilterMutectCalls'
    // non-PASS records (germline, artefacts), which RNA_SUPPORT_FILTER keeps.
    script:
    def hla_str = hla_alleles instanceof List ? hla_alleles.join(',') : hla_alleles
    def normal_arg = params.pvacseq_normal_sample_name ? "--normal-sample-name ${params.pvacseq_normal_sample_name}" : ''
    """
    pvacseq run \\
        ${annotated_vcf} \\
        ${params.pvacseq_tumor_sample_name} \\
        ${hla_str} \\
        ${params.pvacseq_algorithms} \\
        pvacseq_output \\
        ${normal_arg} \\
        --pass-only \\
        -t ${task.cpus}
    """

    stub:
    """
    mkdir -p pvacseq_output
    printf 'Stub run: no biological analysis performed.\\n' > pvacseq_output/stub.log
    """
}
