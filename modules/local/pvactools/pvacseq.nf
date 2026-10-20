// pVACseq on the VEP + readcount annotated somatic VCF.
//   sample "tumor" / --normal-sample-name normal: the read-group SM names
//     set at alignment, so the VCF columns always carry these names.
//   --phased-proximal-variants-vcf: WhatsHap-phased somatic+germline
//     calls, so in-phase neighbouring missense changes are built into the
//     mutant and wild-type peptides.
//   --tumor-purity: Sequenza cellularity from the samplesheet. pVACseq
//     uses it in the aggregate report to call each variant clonal or
//     subclonal from its VAF; without it pVACseq guesses purity from the
//     VAFs themselves. pVACseq has no ploidy input, so ploidy and local
//     copy number are applied afterwards in ANNOTATE_CCF.
//   --iedb-install-directory: the griffithlab/pvactools image ships
//     standalone IEDB at /opt/iedb, which NetMHCpan/NetMHCIIpan need.

process PVACSEQ {
    tag "${meta.id}"
    label 'process_pvacseq'
    container params.containers.pvactools
    publishDir "${params.outdir}/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi), path(phased_vcf), path(phased_tbi), val(hla_alleles), val(purity)

    output:
    tuple val(meta), path("pvacseq"), emit: results

    script:
    def iedb = params.pvacseq_iedb_install_dir ? "--iedb-install-directory ${params.pvacseq_iedb_install_dir}" : ''
    """
    pvacseq run \\
        ${vcf} \\
        tumor \\
        '${hla_alleles}' \\
        ${params.pvacseq_algorithms} \\
        pvacseq \\
        --normal-sample-name normal \\
        --phased-proximal-variants-vcf ${phased_vcf} \\
        --tumor-purity ${purity} \\
        -e1 ${params.pvacseq_class_i_lengths} \\
        -e2 ${params.pvacseq_class_ii_lengths} \\
        --pass-only \\
        ${iedb} \\
        -t ${task.cpus} \\
        ${params.pvacseq_extra_args}
    """

    stub:
    """
    mkdir -p pvacseq/MHC_Class_I pvacseq/MHC_Class_II pvacseq/combined
    printf 'Stub run: no biological analysis performed.\\n' > pvacseq/stub.log
    printf 'ID\\tGene\\tAA Change\\tDNA VAF\\n' > pvacseq/combined/tumor.Combined.all_epitopes.aggregated.tsv
    """
}

// Cancer-cell-fraction estimate per aggregated pVACseq candidate, which is
// where ploidy and copy number come in (bin/annotate_ccf.py). When the
// patient has a week 6 RNA BAM, its RNA depth and VAF are joined on first
// (bin/add_rna_timepoint.py); pVACseq's own RNA columns are the baseline.
process ANNOTATE_CCF {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.pvactools
    publishDir "${params.outdir}/${meta.id}/pvacseq_ccf", mode: 'copy'

    input:
    tuple val(meta), path(pvacseq_dir), val(purity), val(ploidy), val(sex), path(segments), path(rna_week6_vcf)

    output:
    tuple val(meta), path("*.ccf.tsv"), optional: true, emit: tsv

    script:
    def ploidy_arg   = ploidy   ? "--ploidy ${ploidy}" : ''
    def sex_arg      = sex      ? "--sex ${sex}" : ''
    def segments_arg = segments ? "--segments ${segments}" : ''
    """
    # tumor.MHC_I / tumor.MHC_II / tumor.Combined .all_epitopes.aggregated.tsv
    for report in ${pvacseq_dir}/*/tumor.*.all_epitopes.aggregated.tsv; do
        [ -e "\$report" ] || continue
        mhc_class=\$(basename "\$report" | cut -d. -f2)
        if [ -n "${rna_week6_vcf}" ]; then
            add_rna_timepoint.py "\$report" --vcf ${rna_week6_vcf} --label "week 6" -o \${mhc_class}.week6.tsv
            report=\${mhc_class}.week6.tsv
        fi
        annotate_ccf.py "\$report" \\
            --purity ${purity} ${ploidy_arg} ${sex_arg} ${segments_arg} \\
            -o ${meta.id}.\${mhc_class}.aggregated.ccf.tsv
    done
    """

    stub:
    def week6_columns = rna_week6_vcf ? '\\tRNA Depth (week 6)\\tRNA VAF (week 6)' : ''
    """
    printf 'ID${week6_columns}\\tCCF\\n' > ${meta.id}.Combined.aggregated.ccf.tsv
    """
}
