// Ported from: DNA pVACseq script -> trim_and_qc()
// Original ran FastQC before+after trimming as separate docker calls;
// kept as a distinct QC process (fastqc.nf) so it can be skipped/parallelised
// independently and doesn't block trimming.

process TRIM_GALORE {
    tag "${meta.id}:${sample_type}"
    label 'process_medium'
    container params.containers.trimgalore
    publishDir "${params.outdir}/${meta.id}/trimmed", mode: 'copy', pattern: '*_val_*.fq.gz'

    input:
    tuple val(meta), val(sample_type), path(r1), path(r2)

    output:
    tuple val(meta), val(sample_type), path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed_reads
    path "*trimming_report.txt", emit: report

    script:
    """
    trim_galore --paired --cores ${task.cpus} -o . ${r1} ${r2}
    """

    stub:
    """
    touch ${meta.id}.${sample_type}_val_1.fq.gz ${meta.id}.${sample_type}_val_2.fq.gz ${meta.id}.${sample_type}_trimming_report.txt
    """
}
