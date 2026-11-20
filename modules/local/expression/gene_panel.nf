// Status table for the antigen-presentation machinery, the interferon
// route that switches it on, and the melanoma drivers
// (bin/gene_panel_report.py).
//
// Separate from the per-candidate annotation because the question is
// different: not "is this peptide worth making?" but "can this tumour
// present anything at all, and what is driving it?". A B2M loss or a JAK1
// truncation invalidates the whole candidate list rather than any one
// entry, and nothing in a per-candidate report would show it.
//
// Takes whatever the run produced - sample TPM, GTEx/TCGA tables, the
// VEP-annotated somatic VCF, Sequenza segments - and leaves the rest NA.

process GENE_PANEL {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.pvactools
    publishDir "${params.outdir}/${meta.id}/gene_panel", mode: 'copy'

    input:
    tuple val(meta), path(vep_vcf), val(ploidy), path(segments), path(sample_tpm)
    path gtex_db
    path tcga_db
    path gtf
    path panel

    output:
    tuple val(meta), path("${meta.id}.gene_panel.tsv"), emit: tsv

    script:
    def gtex_args    = gtex_db  ? "--gtex ${gtex_db} --gtex-tissue '${params.gtex_tissue}'" : ''
    def tcga_args    = tcga_db  ? "--tcga ${tcga_db} --tcga-cancer-type '${params.tcga_cancer_type}'" : ''
    def tpm_arg      = sample_tpm ? "--sample-tpm ${sample_tpm}" : ''
    def vcf_arg      = vep_vcf ? "--vep-vcf ${vep_vcf}" : ''
    def gtf_arg      = gtf ? "--gtf ${gtf}" : ''
    // copy number needs gene coordinates, so segments only help with a GTF
    def segments_arg = (segments && gtf) ? "--segments ${segments}" : ''
    def ploidy_arg   = ploidy ? "--ploidy ${ploidy}" : ''
    def panel_arg    = panel ? "--panel ${panel}" : ''
    """
    gene_panel_report.py \\
        ${gtex_args} \\
        ${tcga_args} \\
        ${tpm_arg} \\
        ${vcf_arg} \\
        ${gtf_arg} \\
        ${segments_arg} \\
        ${ploidy_arg} \\
        ${panel_arg} \\
        --amp-gain ${params.panel_amp_gain} \\
        --fold ${params.panel_fold} \\
        -o ${meta.id}.gene_panel.tsv
    """

    stub:
    """
    printf 'Pathway\\tGene\\tPanel Flag\\n' > ${meta.id}.gene_panel.tsv
    printf 'apm_class_i\\tB2M\\tOK\\n' >> ${meta.id}.gene_panel.tsv
    """
}
