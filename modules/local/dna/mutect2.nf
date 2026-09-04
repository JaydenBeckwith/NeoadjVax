// Ported from: DNA pVACseq script -> run_mutect2()
// Original called Mutect2 with literal -tumor/-normal sample names "tumor"/
// "normal" (matching the RGSM values it set) — kept the same convention.

process MUTECT2 {
    tag "${meta.id}"
    label 'process_high'
    container params.containers.gatk

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path fasta
    path fasta_fai
    path fasta_dict

    output:
    tuple val(meta), path("${meta.id}_raw.vcf.gz"), path("${meta.id}_raw.vcf.gz.stats"), emit: vcf

    script:
    """
    gatk Mutect2 \\
        -R ${fasta} \\
        -I ${tumor_bam} -tumor tumor \\
        -I ${normal_bam} -normal normal \\
        -O ${meta.id}_raw.vcf.gz
    """
}
