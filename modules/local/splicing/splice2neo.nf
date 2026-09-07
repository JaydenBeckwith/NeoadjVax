process SPLICE_REFERENCE {
    label 'process_high'
    container params.containers.splice2neo
    publishDir "${params.outdir}/splicing_reference", mode: 'copy', pattern: 'reference_summary.txt'
    input:
    path fasta, stageAs: 'genome.fa'
    path gtf
    path normal
    output:
    path 'reference.rds', emit: reference
    path 'genome.fa.fai', emit: fai
    path 'reference_summary.txt', emit: summary
    script:
    def panel = normal ? normal : '-'
    """
    Rscript '${projectDir}/bin/run_splice2neo.R' prepare '${gtf}' genome.fa '${panel}' reference.rds
    """
    stub:
    """
    touch reference.rds genome.fa.fai reference_summary.txt
    """
}

process SPLICEAI_ANNOTATE {
    tag "${meta.id}"
    label 'process_medium'
    container params.splicing_tools_container
    publishDir "${params.outdir}/${meta.id}/splicing/dna", mode: 'copy', pattern: '*.spliceai.vcf'
    input:
    tuple val(meta), path(vcf)
    path fasta, stageAs: 'genome.fa'
    path fai, stageAs: 'genome.fa.fai'
    path annotation
    output:
    tuple val(meta), path("${meta.id}.spliceai.vcf"), emit: vcf
    script:
    def annot = annotation ? annotation : params.spliceai_annotation
    """
    bcftools view -f PASS -Ou '${vcf}' | bcftools norm -f genome.fa -m -any -Ov -o normalized.vcf
    export TF_NUM_INTRAOP_THREADS=${task.cpus}
    export TF_NUM_INTEROP_THREADS=1
    spliceai -I normalized.vcf -O '${meta.id}.spliceai.vcf' -R genome.fa -A '${annot}' -D ${params.spliceai_distance} -M 0
    test -s '${meta.id}.spliceai.vcf'
    """
    stub:
    """
    touch '${meta.id}.spliceai.vcf'
    """
}

process SPLICEAI_FILTER {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.splicing_python
    publishDir "${params.outdir}/${meta.id}/splicing/dna", mode: 'copy'
    input:
    tuple val(meta), path(vcf)
    output:
    tuple val(meta), path('spliceai_effects.tsv'), emit: table
    path 'spliceai_qc.json', emit: qc
    script:
    """
    python '${projectDir}/bin/prepare_splicing_inputs.py' spliceai --input '${vcf}' \\
        --output spliceai_effects.tsv --qc spliceai_qc.json --threshold ${params.spliceai_min_delta_score}
    """
    stub:
    """
    touch spliceai_effects.tsv spliceai_qc.json
    """
}

process SPLICE_PREDICT {
    tag "${meta.id}"
    label 'process_medium'
    container params.containers.splice2neo
    publishDir "${params.outdir}/${meta.id}/splicing/dna", mode: 'copy', pattern: '*.tsv'
    input:
    tuple val(meta), path(effects)
    path reference
    output:
    tuple val(meta), path('predictions.rds'), emit: predictions
    path 'predicted_junctions.tsv', emit: junctions
    script:
    """
    Rscript '${projectDir}/bin/run_splice2neo.R' predict '${effects}' '${reference}' ${params.spliceai_min_delta_score} predictions.rds
    """
    stub:
    """
    touch predictions.rds predicted_junctions.tsv
    """
}

process SPLICE_TO_PEPTIDE {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_medium'
    container params.containers.splice2neo
    publishDir "${params.outdir}/${meta.id}/splicing/${meta.timepoint}", mode: 'copy'
    input:
    tuple val(meta), path(predictions), path(counts), path(evidence)
    path reference
    path fasta, stageAs: 'genome.fa'
    path fai, stageAs: 'genome.fa.fai'
    output:
    tuple val(meta), path('peptide_contexts.fasta'), emit: peptide_fasta
    tuple val(meta), path('peptide_candidates.tsv'), emit: candidates
    tuple val(meta), path('junction_evidence.tsv'), emit: evidence
    tuple val(meta), path('summary.tsv'), emit: summary
    path 'versions.txt', emit: versions
    script:
    """
    Rscript '${projectDir}/bin/run_splice2neo.R' integrate '${predictions}' '${reference}' genome.fa \\
        '${counts}' '${evidence}' '${meta.id}' '${meta.timepoint}' \\
        ${params.splicing_peptide_flank} ${params.splicing_min_peptide_length}
    """
    stub:
    """
    touch peptide_contexts.fasta peptide_candidates.tsv junction_evidence.tsv summary.tsv versions.txt
    """
}
