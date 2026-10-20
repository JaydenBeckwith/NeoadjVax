// Read counts at every somatic site, as pVACtools' input-prep docs lay
// out ("Adding coverage data to your VCF"): bam-readcount via the
// mgibio/bam_readcount_helper-cwl container, then vcf-readcount-annotator
// from VAtools to write the counts into the VCF columns pVACseq reads:
// AD/DP/AF (DNA) in the tumour and normal columns for its DNA VAF/depth
// filters and clonality call, and RAD/RDP/RAF (RNA) in the tumour column
// for its RNA VAF/depth filters.
//
// Helper arguments (confirmed from the script inside the pinned 1.1.1
// image): <vcf> <sample> <fasta> <bam> <prefix|NOPREFIX> <outdir>
// [min_base_quality] [min_mapping_quality]. <sample> must be a VCF column:
// the helper counts at the ALT alleles in that sample's genotype. It writes
// <sample>_bam_readcount_snv.tsv and <sample>_bam_readcount_indel.tsv,
// renamed here to <sample_type>_... so RNA counts (taken at the tumour's
// sites) do not collide with the DNA ones.

process BAM_READCOUNT {
    tag "${meta.id}:${sample_type}"
    label 'process_low'
    container params.containers.bam_readcount

    input:
    tuple val(meta), val(sample_type), val(vcf_sample), path(bam), path(bai), path(vcf)
    path fasta
    path fasta_fai

    output:
    tuple val(meta), val(sample_type), path("${sample_type}_bam_readcount_snv.tsv"), path("${sample_type}_bam_readcount_indel.tsv"), emit: readcounts

    script:
    """
    # bam-readcount looks for <bam>.bai; an index staged as <x>.bai is linked
    [ -e ${bam}.bai ] || ln -s ${bai} ${bam}.bai
    /usr/bin/python /usr/bin/bam_readcount_helper.py \\
        ${vcf} ${vcf_sample} ${fasta} ${bam} NOPREFIX . \\
        ${params.readcount_min_base_quality} ${params.readcount_min_mapping_quality}
    if [ "${vcf_sample}" != "${sample_type}" ]; then
        mv ${vcf_sample}_bam_readcount_snv.tsv ${sample_type}_bam_readcount_snv.tsv
        mv ${vcf_sample}_bam_readcount_indel.tsv ${sample_type}_bam_readcount_indel.tsv
    fi
    """

    stub:
    """
    touch ${sample_type}_bam_readcount_snv.tsv ${sample_type}_bam_readcount_indel.tsv
    """
}

// Writes each <sample_type>_bam_readcount_<snv|indel>.tsv it is given into
// the VCF: tumor/normal as DNA into that column, rna_* as RNA into the
// tumor column. Empty files (no sites of that type) are skipped.
process VCF_READCOUNT_ANNOTATOR {
    tag "${meta.id}:${name}"
    label 'process_low'
    container params.containers.vatools

    input:
    tuple val(meta), path(vcf), path(readcounts, stageAs: 'readcounts/*'), val(name)

    output:
    tuple val(meta), path("${meta.id}.${name}.vcf"), emit: vcf

    script:
    """
    cp ${vcf} current.vcf
    for tsv in readcounts/*_bam_readcount_*.tsv; do
        [ -s "\$tsv" ] || continue
        file=\$(basename "\$tsv" .tsv)
        sample_type=\${file%%_bam_readcount_*}
        variant_type=\${file##*_bam_readcount_}
        case "\$sample_type" in
            tumor|normal) sample=\$sample_type; data=DNA ;;
            rna_*)        sample=tumor;         data=RNA ;;
            *) echo "unexpected readcount file \$tsv" >&2; exit 1 ;;
        esac
        vcf-readcount-annotator current.vcf "\$tsv" \$data -s \$sample -t \$variant_type -o next.vcf
        mv next.vcf current.vcf
    done
    mv current.vcf ${meta.id}.${name}.vcf
    """

    stub:
    """
    touch ${meta.id}.${name}.vcf
    """
}
