// VEP annotation with the options pVACseq's input-prep docs list as
// required (--symbol --terms SO --mane_select --canonical --tsl --biotype
// --hgvs --fasta --offline --cache) plus pVACtools' Wildtype and
// Frameshift plugins, which must be in params.vep_plugins. Used twice per
// patient with the same settings: on the somatic VCF that pVACseq scores,
// and on the phased somatic+germline VCF it reads proximal variants from;
// both need matching CSQ transcripts. Output is plain VCF because the
// readcount steps after it rewrite the file anyway.

process VEP_ANNOTATE {
    tag "${meta.id}:${name}"
    label 'process_medium'
    container params.containers.vep

    input:
    tuple val(meta), path(vcf), val(name)
    path fasta
    path fasta_fai
    path vep_cache
    path vep_plugins

    output:
    tuple val(meta), path("${meta.id}.${name}.vep.vcf"), emit: vcf

    script:
    """
    vep \\
        --input_file ${vcf} \\
        --output_file ${meta.id}.${name}.vep.vcf \\
        --format vcf --vcf \\
        --assembly ${params.vep_assembly} \\
        --cache_version ${params.vep_cache_version} \\
        --offline --cache --dir_cache ${vep_cache} \\
        --dir_plugins ${vep_plugins} \\
        --fasta ${fasta} \\
        --symbol --terms SO --mane_select --canonical --tsl --biotype --hgvs \\
        --plugin Frameshift --plugin Wildtype \\
        --force_overwrite --no_stats \\
        --fork ${task.cpus}
    """

    stub:
    """
    touch ${meta.id}.${name}.vep.vcf
    """
}
