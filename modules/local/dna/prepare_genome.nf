// Shared genome prep (not present in the original scripts, which assumed
// these existed already): needed once so GATK tools have a .fai/.dict.
// Used by both the DNA and RNA branches.

process SAMTOOLS_FAIDX {
    tag "${fasta.simpleName}"
    label 'process_low'
    container params.containers.samtools

    input:
    path fasta

    output:
    path "${fasta}.fai", emit: fai

    script:
    """
    if [ -f "${fasta}.fai" ]; then
        cp ${fasta}.fai .
    else
        samtools faidx ${fasta}
    fi
    """

    stub:
    """
    touch ${fasta}.fai
    """
}

process PICARD_CREATE_SEQUENCE_DICTIONARY {
    tag "${fasta.simpleName}"
    label 'process_low'
    container params.containers.picard

    input:
    path fasta

    output:
    path "${fasta.baseName}.dict", emit: dict

    script:
    """
    picard CreateSequenceDictionary R=${fasta} O=${fasta.baseName}.dict
    """

    stub:
    """
    touch ${fasta.baseName}.dict
    """
}
