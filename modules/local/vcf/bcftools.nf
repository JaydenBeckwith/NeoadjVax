// Small bcftools steps around the GATK/VEP/pVACtools calls.

// Somatic calls that go forward to pVACseq: FilterMutectCalls PASS only,
// left-aligned and split to one ALT per record. Splitting is what
// pVACtools asks for before bam-readcount (its docs use vt decompose;
// bcftools norm -m -any is the same operation), and doing it here means
// the main pVACseq VCF and the phasing VCF carry identical records.
process SELECT_PASS_SOMATIC {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.bcftools
    publishDir "${params.outdir}/${meta.id}/somatic_variants", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path fasta
    path fasta_fai

    output:
    tuple val(meta), path("${meta.id}.somatic.pass.vcf.gz"), path("${meta.id}.somatic.pass.vcf.gz.tbi"), emit: vcf

    script:
    """
    bcftools view -f PASS ${vcf} \\
        | bcftools norm -m -any -f ${fasta} -Oz -o ${meta.id}.somatic.pass.vcf.gz
    bcftools index -t ${meta.id}.somatic.pass.vcf.gz
    """

    stub:
    """
    touch ${meta.id}.somatic.pass.vcf.gz ${meta.id}.somatic.pass.vcf.gz.tbi
    """
}

// BED of windows around each somatic site: germline calling for phasing
// only needs the neighbourhood of somatic variants, not the whole genome.
process SOMATIC_WINDOWS_BED {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.bcftools

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.id}.phasing_windows.bed"), emit: bed

    script:
    def window = params.phasing_window_bp
    """
    bcftools query -f '%CHROM\\t%POS0\\t%END\\n' ${vcf} \\
        | awk -v OFS='\\t' -v w=${window} '{ s = \$2 - w; if (s < 0) s = 0; print \$1, s, \$3 + w }' \\
        > ${meta.id}.phasing_windows.bed
    """

    stub:
    """
    touch ${meta.id}.phasing_windows.bed
    """
}

// Builds the single-sample VCF WhatsHap phases against the tumour BAM, as
// pVACseq's proximal-variant docs describe: the somatic calls restricted
// to the tumour column, plus the normal's germline calls (quality-filtered
// and renamed to "tumor" so the sample names match the tumour BAM),
// combined and sorted.
process COMBINE_SOMATIC_GERMLINE {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.bcftools

    input:
    tuple val(meta), path(somatic_vcf), path(somatic_tbi), path(germline_vcf), path(germline_tbi)
    path fasta
    path fasta_fai

    output:
    tuple val(meta), path("${meta.id}.somatic_germline.vcf.gz"), path("${meta.id}.somatic_germline.vcf.gz.tbi"), emit: vcf

    script:
    """
    bcftools view -s tumor -c 1 ${somatic_vcf} -Oz -o somatic.tumor.vcf.gz
    bcftools index -t somatic.tumor.vcf.gz

    echo 'normal tumor' > rename.txt
    bcftools view -i 'QUAL>=${params.germline_min_qual} && GT="alt"' ${germline_vcf} \\
        | bcftools norm -m -any -f ${fasta} \\
        | bcftools reheader -s rename.txt \\
        | bcftools view -Oz -o germline.tumor.vcf.gz
    bcftools index -t germline.tumor.vcf.gz

    bcftools concat -a -D somatic.tumor.vcf.gz germline.tumor.vcf.gz \\
        | bcftools sort -Oz -o ${meta.id}.somatic_germline.vcf.gz
    bcftools index -t ${meta.id}.somatic_germline.vcf.gz
    """

    stub:
    """
    touch ${meta.id}.somatic_germline.vcf.gz ${meta.id}.somatic_germline.vcf.gz.tbi
    """
}

process BGZIP_INDEX {
    tag "${meta.id}:${name}"
    label 'process_low'
    container params.containers.bcftools
    publishDir "${params.outdir}/${meta.id}/pvacseq_input", mode: 'copy'

    input:
    tuple val(meta), path(vcf), val(name)

    output:
    tuple val(meta), path("${meta.id}.${name}.vcf.gz"), path("${meta.id}.${name}.vcf.gz.tbi"), emit: vcf

    script:
    """
    bcftools view ${vcf} -Oz -o ${meta.id}.${name}.vcf.gz
    bcftools index -t ${meta.id}.${name}.vcf.gz
    """

    stub:
    """
    touch ${meta.id}.${name}.vcf.gz ${meta.id}.${name}.vcf.gz.tbi
    """
}
