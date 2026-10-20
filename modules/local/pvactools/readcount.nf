// Read counts at every somatic site, as pVACtools' input-prep docs lay
// out ("Adding coverage data to your VCF"): bam-readcount via the
// mgibio/bam_readcount_helper-cwl container, then vcf-readcount-annotator
// from VAtools to write AD/DP/AF into the tumour and normal columns that
// pVACseq reads for its DNA VAF/depth filters and clonality call.
//
// Helper arguments (confirmed from the script inside the pinned 1.1.1
// image): <vcf> <sample> <fasta> <bam> <prefix|NOPREFIX> <outdir>
// [min_base_quality] [min_mapping_quality]. It writes
// <sample>_bam_readcount_snv.tsv and <sample>_bam_readcount_indel.tsv.

process BAM_READCOUNT {
    tag "${meta.id}:${sample_type}"
    label 'process_low'
    container params.containers.bam_readcount

    input:
    tuple val(meta), val(sample_type), path(bam), path(bai), path(vcf)
    path fasta
    path fasta_fai

    output:
    tuple val(meta), val(sample_type), path("${sample_type}_bam_readcount_snv.tsv"), path("${sample_type}_bam_readcount_indel.tsv"), emit: readcounts

    script:
    """
    /usr/bin/python /usr/bin/bam_readcount_helper.py \\
        ${vcf} ${sample_type} ${fasta} ${bam} NOPREFIX . \\
        ${params.readcount_min_base_quality} ${params.readcount_min_mapping_quality}
    """

    stub:
    """
    touch ${sample_type}_bam_readcount_snv.tsv ${sample_type}_bam_readcount_indel.tsv
    """
}

process VCF_READCOUNT_ANNOTATOR {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.vatools

    input:
    tuple val(meta), path(vcf), path(tumor_snv), path(tumor_indel), path(normal_snv), path(normal_indel)

    output:
    tuple val(meta), path("${meta.id}.readcounts.vcf"), emit: vcf

    script:
    """
    cp ${vcf} current.vcf
    annotate() {   # <readcount tsv> <sample> <snv|indel>
        if [ -s "\$1" ]; then
            vcf-readcount-annotator current.vcf "\$1" DNA -s "\$2" -t "\$3" -o next.vcf
            mv next.vcf current.vcf
        fi
    }
    annotate ${tumor_snv} tumor snv
    annotate ${tumor_indel} tumor indel
    annotate ${normal_snv} normal snv
    annotate ${normal_indel} normal indel
    mv current.vcf ${meta.id}.readcounts.vcf
    """

    stub:
    """
    touch ${meta.id}.readcounts.vcf
    """
}
