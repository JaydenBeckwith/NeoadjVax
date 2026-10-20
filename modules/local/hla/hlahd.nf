// HLA typing with HLA-HD (Kawaguchi et al. 2017) on the matched NORMAL
// DNA, so tumour HLA loss cannot bias the call. HLA-HD types class I
// (A/B/C) and class II including both alpha and beta chains
// (DQA1/DQB1, DPA1/DPB1), which is what pVACseq needs to form class II
// heterodimers.
//
// HLA-HD only needs reads from the HLA region, so the normal BAM is cut
// down first: reads aligned to the MHC (params.hlahd_regions) plus the
// unplaced unmapped reads, written back out as paired FASTQ. Typing the
// whole WGS FASTQ would give the same answer for ~100x the bowtie2 time.
//
// HLA-HD is distributed under an academic licence and cannot be pulled
// from a public registry: build containers/hlahd from your own download
// and point params.hlahd_container at the resulting image.

process EXTRACT_HLA_READS {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.samtools

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.hla_R1.fastq"), path("${meta.id}.hla_R2.fastq"), emit: reads

    script:
    """
    samtools view -@ ${task.cpus} -b -o mhc.bam ${bam} ${params.hlahd_regions}
    samtools view -@ ${task.cpus} -b -o unplaced.bam ${bam} '*'
    samtools cat -o hla_candidates.bam mhc.bam unplaced.bam
    samtools sort -n -@ ${task.cpus} -o hla_candidates.qsorted.bam hla_candidates.bam
    samtools fastq -@ ${task.cpus} \\
        -1 ${meta.id}.hla_R1.fastq -2 ${meta.id}.hla_R2.fastq \\
        -0 /dev/null -s /dev/null \\
        hla_candidates.qsorted.bam
    rm -f mhc.bam unplaced.bam hla_candidates.bam hla_candidates.qsorted.bam
    """

    stub:
    """
    touch ${meta.id}.hla_R1.fastq ${meta.id}.hla_R2.fastq
    """
}

process HLAHD {
    tag "${meta.id}"
    label 'process_medium'
    container params.hlahd_container
    publishDir "${params.outdir}/${meta.id}/hla_typing", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}.pvacseq_alleles.txt"), emit: alleles
    path "${meta.id}_final.result.txt", emit: result
    path "${meta.id}.hla_alleles.tsv", emit: table

    script:
    def hlahd_dir = params.hlahd_install_dir
    """
    hlahd.sh \\
        -t ${task.cpus} \\
        -m ${params.hlahd_min_read_length} \\
        -f ${hlahd_dir}/freq_data/ \\
        ${r1} ${r2} \\
        ${hlahd_dir}/HLA_gene.split.txt \\
        ${hlahd_dir}/dictionary/ \\
        ${meta.id} \\
        .

    cp ${meta.id}/result/${meta.id}_final.result.txt .
    parse_hlahd_result.py ${meta.id}_final.result.txt \\
        --pvacseq-out ${meta.id}.pvacseq_alleles.txt \\
        --table-out ${meta.id}.hla_alleles.tsv
    """

    stub:
    """
    printf 'A\\tHLA-A*02:01:01:01\\tHLA-A*24:02:01:01\\n' > ${meta.id}_final.result.txt
    printf 'HLA-A*02:01,HLA-A*24:02\\n' > ${meta.id}.pvacseq_alleles.txt
    printf 'gene\\tallele_1\\tallele_2\\n' > ${meta.id}.hla_alleles.tsv
    """
}
