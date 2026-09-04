// Ported from: DNA pVACseq script -> run_vep()
// Wildtype + Frameshift plugins are required by pVACseq for correct
// wild-type/mutant peptide pairing and frameshift neoantigen calling —
// kept exactly as specified in the original script.

process VEP_ANNOTATE {
    tag "${meta.id}"
    label 'process_medium'
    container params.containers.vep
    publishDir "${params.outdir}/${meta.id}/dna_variants", mode: 'copy'

    input:
    tuple val(meta), path(vcf)
    path fasta
    path vep_cache
    path vep_plugins

    output:
    tuple val(meta), path("${meta.id}.annotated.vcf.gz"), emit: vcf

    script:
    """
    vep -i ${vcf} \\
        -o ${meta.id}.annotated.vcf.gz \\
        --format vcf --vcf --verbose --assembly GRCh38 \\
        --symbol --canonical --distance 5 \\
        --plugin Wildtype,${fasta} \\
        --plugin Frameshift,${fasta} \\
        --fasta ${fasta} \\
        --offline --cache \\
        --dir_plugins ${vep_plugins} \\
        --dir_cache ${vep_cache} \\
        --tsl --biotype --hgvs \\
        --terms SO \\
        --force_overwrite \\
        --compress_output bgzip \\
        --fork ${task.cpus}
    """
}
