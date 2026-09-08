// Ported from: DNA pVACseq script -> trim_and_qc() (the FastQC docker calls,
// split out from Trim Galore so before/after QC can run independently)

process FASTQC {
    tag "${meta.id}:${sample_type}:${stage}"
    label 'process_low'
    container params.containers.fastqc
    publishDir "${params.outdir}/${meta.id}/fastqc/${stage}", mode: 'copy'

    input:
    tuple val(meta), val(sample_type), path(r1), path(r2)
    val stage   // 'pre_trim' | 'post_trim'

    output:
    path "*.{html,zip}", emit: qc

    script:
    """
    fastqc -o . ${r1} ${r2}
    """

    stub:
    """
    touch ${meta.id}.${sample_type}.${stage}_stub.html ${meta.id}.${sample_type}.${stage}_stub.zip
    """
}
