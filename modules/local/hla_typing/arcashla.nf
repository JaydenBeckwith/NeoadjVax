// TODO stub — HLA class II typing, not yet built or tested.
// [[neoantigen-score]] specifies "HLA alleles typed for MHC class I and II
// from the WGS" but neither original script did class II typing at all.
// arcasHLA is the natural fit since it works directly off an RNA (or DNA)
// BAM we already have — but this process is unverified: command syntax,
// container, and output parsing all need a real test run before use.
// Disabled by default (params.run_hla_typing_class_ii = false).

process HLA_TYPING_ARCASHLA {
    tag "${meta.id}"
    label 'process_medium'
    container params.containers.arcashla

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("arcashla_out/*.genotype.json"), emit: genotype

    script:
    """
    # TODO — unverified. Rough shape per arcasHLA docs:
    #   arcasHLA extract ${bam} -o extracted --paired -t ${task.cpus}
    #   arcasHLA genotype extracted/*_1.fastq.gz extracted/*_2.fastq.gz \\
    #       -g A,B,C,DPB1,DQB1,DQA1,DRB1 -o arcashla_out -t ${task.cpus}
    echo "HLA_TYPING_ARCASHLA is a stub — not implemented yet" >&2
    exit 1
    """
}
