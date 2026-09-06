// v1.4 CLI verified against upstream, pinned in containers/ervcaller.
// The runner validates output because upstream can exit zero on errors.
process ERVCALLER_RUN {
    tag "${meta.id}:${sample_type}"
    label 'process_medium'
    container params.ervcaller_container
    publishDir "${params.outdir}/${meta.id}/erv_dna/${sample_type}", mode: 'copy'

    input:
    tuple val(meta), val(sample_type), path('input.bam'), path('input.bam.bai')
    path references, stageAs: 'refs'
    path runner, stageAs: 'run_ervcaller.py'

    output:
    tuple val(meta), val(sample_type), path('results/*.vcf'), emit: calls
    tuple val(meta), val(sample_type), path('results/*.evidence.tsv'), emit: evidence
    tuple val(meta), val(sample_type), path('results/*.status.json'), emit: status
    path 'results/*.log', emit: logs

    script:
    def bwa = params.ervcaller_bwa_mem ? '--bwa-mem' : ''
    def genotype = params.ervcaller_genotype ? '--genotype' : ''
    """
    python3 run_ervcaller.py \\
        --patient '${meta.id}' --role '${sample_type}' \\
        --bam input.bam --bai input.bam.bai --references refs \\
        --caller '${params.ervcaller_script}' \\
        --read-length ${params.ervcaller_read_length} \\
        --min-reads ${params.ervcaller_min_reads} \\
        --min-split ${params.ervcaller_min_split} \\
        --threads ${task.cpus} ${bwa} ${genotype} --outdir results
    """

    stub:
    """
    mkdir results
    printf '##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t${meta.id}.${sample_type}\\n' > results/${meta.id}.${sample_type}.vcf
    printf 'patient_id\\tsample_type\\n' > results/${meta.id}.${sample_type}.evidence.tsv
    printf '{"stub":true}\\n' > results/${meta.id}.${sample_type}.status.json
    printf 'Stub run: no biological analysis performed.\\n' > results/${meta.id}.${sample_type}.log
    """
}
