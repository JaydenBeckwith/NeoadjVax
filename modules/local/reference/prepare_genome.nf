// Reference indexes the pipeline needs once per run. main.nf only calls
// these when the file is missing next to params.genome_fasta,
// so a pre-indexed Gadi reference is reused as-is
// rather than rebuilt: BWA_INDEX alone is several CPU-hours on GRCh38.

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
    samtools faidx ${fasta}
    """

    stub:
    """
    touch ${fasta}.fai
    """
}

process CREATE_SEQUENCE_DICTIONARY {
    tag "${fasta.simpleName}"
    label 'process_low'
    container params.containers.gatk

    input:
    path fasta

    output:
    path "${fasta.baseName}.dict", emit: dict

    script:
    """
    gatk CreateSequenceDictionary -R ${fasta} -O ${fasta.baseName}.dict
    """

    stub:
    """
    touch ${fasta.baseName}.dict
    """
}

process BWA_INDEX {
    tag "${fasta.simpleName}"
    label 'process_medium'
    container params.containers.bwa

    input:
    path fasta

    output:
    path "${fasta}.{amb,ann,bwt,pac,sa}", emit: index

    script:
    """
    bwa index ${fasta}
    """

    stub:
    """
    touch ${fasta}.amb ${fasta}.ann ${fasta}.bwt ${fasta}.pac ${fasta}.sa
    """
}
