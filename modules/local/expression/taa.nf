// Shared tumour antigens: the second channel beside the patient's own
// mutations (bin/taa_report.py, assets/taa_panel.tsv).
//
// A low-mutation-burden patient can run out of neoantigens. Cancer-testis
// antigens, the melanocyte differentiation antigens and the broadly
// overexpressed genes are not patient-specific, so they are not
// neoantigens, but they are real targets with characterised epitopes and
// whether they are usable in a patient is an expression question this
// pipeline already answers: expressed in this tumour, and silent in
// normal tissue outside testis.

process TAA_REPORT {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.pvactools
    publishDir "${params.outdir}/${meta.id}/taa", mode: 'copy'

    input:
    tuple val(meta), path(sample_tpm)
    path panel
    path gtex_db
    path tcga_db

    output:
    tuple val(meta), path("${meta.id}.taa.tsv"), emit: tsv

    script:
    def gtex_args = gtex_db ? "--gtex ${gtex_db} --gtex-tissue '${params.gtex_tissue}'" : ''
    def tcga_args = tcga_db ? "--tcga ${tcga_db} --tcga-cancer-type '${params.tcga_cancer_type}'" : ''
    def tpm_arg   = sample_tpm ? "--sample-tpm ${sample_tpm}" : ''
    """
    taa_report.py \\
        --panel ${panel} \\
        ${tpm_arg} \\
        ${gtex_args} \\
        ${tcga_args} \\
        --quantile ${params.expression_quantile} \\
        --min-tpm ${params.expression_min_tpm} \\
        --expressed-tpm ${params.gtex_expressed_tpm} \\
        -o ${meta.id}.taa.tsv
    """

    stub:
    """
    printf 'Class\\tGene\\tSample TPM\\tGTEx Specificity\\tCall\\n' > ${meta.id}.taa.tsv
    printf 'cancer_testis\\tMAGEA3\\t25\\ttestis_restricted\\tgood_target\\n' >> ${meta.id}.taa.tsv
    """
}

// Binding predictions for a peptide set that is not derived from this
// patient's variants: the published epitopes of the shared antigens
// above, or any peptide list of your own. pVACbind takes peptides
// directly, so nothing about it depends on the somatic calls - only on
// the patient's alleles.
//
// Needs --taa_peptides: no epitope set ships with this repository,
// because which epitopes are worth testing is a judgement about the
// literature rather than a pipeline default.
process PVACBIND {
    tag "${meta.id}"
    label 'process_pvacseq'
    container params.containers.pvactools
    publishDir "${params.outdir}/${meta.id}/taa", mode: 'copy'

    input:
    tuple val(meta), val(hla_alleles)
    path peptides

    output:
    tuple val(meta), path("pvacbind"), emit: results

    script:
    def iedb = params.pvacseq_iedb_install_dir ? "--iedb-install-directory ${params.pvacseq_iedb_install_dir}" : ''
    """
    pvacbind run \\
        ${peptides} \\
        ${meta.id} \\
        '${hla_alleles}' \\
        ${params.taa_algorithms} \\
        pvacbind \\
        -e1 ${params.pvacseq_class_i_lengths} \\
        -e2 ${params.pvacseq_class_ii_lengths} \\
        ${iedb} \\
        -t ${task.cpus} \\
        ${params.taa_extra_args}
    """

    stub:
    """
    mkdir -p pvacbind
    printf 'Stub run: no biological analysis performed.\\n' > pvacbind/stub.log
    """
}
