// Gene-level expression from an already-aligned RNA BAM.
//
// The samplesheet's RNA BAMs are on the same GRCh38 reference as the DNA
// (main.nf checks they have an index and they are never realigned), so
// expression needs counting only: featureCounts over the GTF's exons, then
// bin/gene_tpm.py to turn counts plus union-exon length into TPM. TPM is
// the unit the GTEx and TCGA reference tables use, so candidates can be
// compared against them without a second normalisation step.
//
// Strandedness matters: counting a reverse-stranded library as unstranded
// inflates antisense overlap, and counting it with the wrong strand loses
// most reads. params.rna_strandedness mirrors featureCounts -s
// (0 unstranded, 1 forward, 2 reverse) and has no safe default, so it is
// stated explicitly per run.

process FEATURECOUNTS {
    tag "${meta.id}:${sample_type}"
    label 'process_medium'
    container params.containers.subread
    publishDir "${params.outdir}/${meta.id}/expression", mode: 'copy', pattern: '*.summary'

    input:
    tuple val(meta), val(sample_type), path(bam), path(bai)
    path gtf

    output:
    tuple val(meta), val(sample_type), path("${meta.id}.${sample_type}.counts.txt"), emit: counts
    path "${meta.id}.${sample_type}.counts.txt.summary", emit: summary

    script:
    def paired = params.rna_paired_end ? '-p --countReadPairs' : ''
    """
    featureCounts \\
        -a ${gtf} \\
        -o ${meta.id}.${sample_type}.counts.txt \\
        -g gene_id \\
        -t exon \\
        -s ${params.rna_strandedness} \\
        ${paired} \\
        -T ${task.cpus} \\
        ${params.featurecounts_extra_args} \\
        ${bam}
    """

    stub:
    """
    printf '# Program:featureCounts stub\\n' > ${meta.id}.${sample_type}.counts.txt
    printf 'Geneid\\tChr\\tStart\\tEnd\\tStrand\\tLength\\t${bam}\\n' >> ${meta.id}.${sample_type}.counts.txt
    printf 'ENSG00000141510\\tchr17\\t1\\t100\\t-\\t2000\\t400\\n' >> ${meta.id}.${sample_type}.counts.txt
    printf 'Status\\t${bam}\\nAssigned\\t400\\n' > ${meta.id}.${sample_type}.counts.txt.summary
    """
}

process GENE_TPM {
    tag "${meta.id}:${sample_type}"
    label 'process_low'
    container params.containers.pvactools
    publishDir "${params.outdir}/${meta.id}/expression", mode: 'copy'

    input:
    tuple val(meta), val(sample_type), path(counts)
    path gtf

    output:
    tuple val(meta), val(sample_type), path("${meta.id}.${sample_type}.tpm.tsv"), emit: tpm

    script:
    def gtf_arg = gtf ? "--gtf ${gtf}" : ''
    """
    gene_tpm.py ${counts} ${gtf_arg} -o ${meta.id}.${sample_type}.tpm.tsv
    """

    stub:
    """
    printf 'ENSEMBL_ID\\tGENE_NAME\\tGENE_BIOTYPE\\tCOUNTS\\tLENGTH\\tTPM\\tEXPRESSION_QUARTILE\\n' > ${meta.id}.${sample_type}.tpm.tsv
    printf 'ENSG00000141510\\tTP53\\tprotein_coding\\t400\\t2000\\t120.0000\\t4\\n' >> ${meta.id}.${sample_type}.tpm.tsv
    """
}
