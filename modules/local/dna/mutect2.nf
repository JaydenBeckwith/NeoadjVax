// GATK somatic short-variant calling, following GATK's tumour/matched-normal
// best practice: Mutect2 with a germline resource and (optionally) a panel
// of normals, read-orientation (FFPE/OxoG) model, cross-sample
// contamination, then FilterMutectCalls using all three. Mutect2 itself
// is scattered by chromosome (one task per entry of
// params.scatter_chromosomes) because a single WGS Mutect2 runs for days.

process MUTECT2 {
    tag "${meta.id}:${interval}"
    label 'process_mutect2'
    container params.containers.gatk

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), val(interval)
    path fasta
    path fasta_fai
    path fasta_dict
    tuple path(germline_resource), path(germline_resource_tbi)
    path pon_files                 // [pon.vcf.gz, pon.vcf.gz.tbi], or [] without a PoN
    path target_intervals          // WES targets (BED/interval_list), or [] for WGS

    output:
    tuple val(meta), path("${meta.id}.${interval}.vcf.gz"), path("${meta.id}.${interval}.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("${meta.id}.${interval}.vcf.gz.stats"), emit: stats
    tuple val(meta), path("${meta.id}.${interval}.f1r2.tar.gz"), emit: f1r2

    script:
    def pon = pon_files ? "--panel-of-normals ${pon_files.find { it.name.endsWith('.vcf.gz') }}" : ''
    def targets = target_intervals ? "-L ${target_intervals} --interval-set-rule INTERSECTION" : ''
    """
    gatk --java-options "-Xmx${(task.memory.toGiga() * 0.8) as int}g" Mutect2 \\
        -R ${fasta} \\
        -I ${tumor_bam} \\
        -I ${normal_bam} \\
        -normal normal \\
        --germline-resource ${germline_resource} \\
        ${pon} \\
        -L ${interval} ${targets} \\
        --f1r2-tar-gz ${meta.id}.${interval}.f1r2.tar.gz \\
        --native-pair-hmm-threads ${task.cpus} \\
        -O ${meta.id}.${interval}.vcf.gz
    """

    stub:
    """
    touch ${meta.id}.${interval}.vcf.gz ${meta.id}.${interval}.vcf.gz.tbi ${meta.id}.${interval}.vcf.gz.stats ${meta.id}.${interval}.f1r2.tar.gz
    """
}

process MERGE_MUTECT2_SCATTERS {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.gatk

    input:
    tuple val(meta), path(vcfs), path(tbis), path(stats), path(f1r2s)

    output:
    tuple val(meta), path("${meta.id}.mutect2.unfiltered.vcf.gz"), path("${meta.id}.mutect2.unfiltered.vcf.gz.tbi"), path("${meta.id}.mutect2.unfiltered.vcf.gz.stats"), path("${meta.id}.read-orientation-model.tar.gz"), emit: merged

    script:
    def vcf_args   = vcfs.collect { "-I ${it}" }.join(' ')
    def stats_args = stats.collect { "-stats ${it}" }.join(' ')
    def f1r2_args  = f1r2s.collect { "-I ${it}" }.join(' ')
    """
    gatk MergeVcfs ${vcf_args} -O ${meta.id}.mutect2.unfiltered.vcf.gz
    gatk MergeMutectStats ${stats_args} -O ${meta.id}.mutect2.unfiltered.vcf.gz.stats
    gatk LearnReadOrientationModel ${f1r2_args} -O ${meta.id}.read-orientation-model.tar.gz
    """

    stub:
    """
    touch ${meta.id}.mutect2.unfiltered.vcf.gz ${meta.id}.mutect2.unfiltered.vcf.gz.tbi ${meta.id}.mutect2.unfiltered.vcf.gz.stats ${meta.id}.read-orientation-model.tar.gz
    """
}

// Pileups at common biallelic germline sites (GATK's small_exac_common_3),
// run on tumour and normal separately.
process GET_PILEUP_SUMMARIES {
    tag "${meta.id}:${sample_type}"
    label 'process_low'
    container params.containers.gatk

    input:
    tuple val(meta), val(sample_type), path(bam), path(bai)
    tuple path(common_sites), path(common_sites_tbi)
    path target_intervals

    output:
    tuple val(meta), val(sample_type), path("${meta.id}.${sample_type}.pileups.table"), emit: table

    script:
    def targets = target_intervals ? "-L ${target_intervals} --interval-set-rule INTERSECTION" : ''
    """
    gatk --java-options "-Xmx${(task.memory.toGiga() * 0.8) as int}g" GetPileupSummaries \\
        -I ${bam} \\
        -V ${common_sites} \\
        -L ${common_sites} ${targets} \\
        -O ${meta.id}.${sample_type}.pileups.table
    """

    stub:
    """
    touch ${meta.id}.${sample_type}.pileups.table
    """
}

process CALCULATE_CONTAMINATION {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.gatk
    publishDir "${params.outdir}/${meta.id}/somatic_variants", mode: 'copy'

    input:
    tuple val(meta), path(tumor_table), path(normal_table)

    output:
    tuple val(meta), path("${meta.id}.contamination.table"), path("${meta.id}.tumor_segmentation.table"), emit: tables

    script:
    """
    gatk CalculateContamination \\
        -I ${tumor_table} \\
        -matched ${normal_table} \\
        --tumor-segmentation ${meta.id}.tumor_segmentation.table \\
        -O ${meta.id}.contamination.table
    """

    stub:
    """
    touch ${meta.id}.contamination.table ${meta.id}.tumor_segmentation.table
    """
}

process FILTER_MUTECT_CALLS {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.gatk
    publishDir "${params.outdir}/${meta.id}/somatic_variants", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi), path(stats), path(orientation_model), path(contamination_files)
    path fasta
    path fasta_fai
    path fasta_dict

    output:
    tuple val(meta), path("${meta.id}.mutect2.filtered.vcf.gz"), path("${meta.id}.mutect2.filtered.vcf.gz.tbi"), emit: vcf
    path "${meta.id}.mutect2.filtering_stats.tsv", emit: filtering_stats

    script:
    def contamination = contamination_files.find { it.name.endsWith('.contamination.table') }
    def segmentation  = contamination_files.find { it.name.endsWith('.tumor_segmentation.table') }
    def contamination_args = contamination ? "--contamination-table ${contamination} --tumor-segmentation ${segmentation}" : ''
    """
    gatk FilterMutectCalls \\
        -R ${fasta} \\
        -V ${vcf} \\
        --stats ${stats} \\
        --ob-priors ${orientation_model} \\
        ${contamination_args} \\
        --filtering-stats ${meta.id}.mutect2.filtering_stats.tsv \\
        -O ${meta.id}.mutect2.filtered.vcf.gz
    """

    stub:
    """
    touch ${meta.id}.mutect2.filtered.vcf.gz ${meta.id}.mutect2.filtered.vcf.gz.tbi ${meta.id}.mutect2.filtering_stats.tsv
    """
}
