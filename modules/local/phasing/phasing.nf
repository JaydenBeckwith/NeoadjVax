// Proximal-variant phasing for pVACseq's --phased-proximal-variants-vcf.
// Without it, pVACseq builds each peptide as if the somatic variant were
// the only change in the region, so an in-phase germline or somatic
// neighbour gives a wrong mutant (and wild-type) sequence.
//
// pVACtools' own recipe phases with GATK 3.6 ReadBackedPhasing, which no
// current GATK ships (and hangs on Gadi compute nodes). WhatsHap with
// --tag HP writes the same HP-style haplotype tags ReadBackedPhasing did.
// pVACseq 7.1.3 reads only HP (pvactools/lib/proximal_variant.py): with
// WhatsHap's default PS tags it would ignore the phasing and keep only
// homozygous neighbours.

// Germline calls in the matched normal, only inside the windows around
// somatic sites. A patient with no PASS somatic calls gets an empty
// window BED; HaplotypeCaller then runs over a single base so the step
// still produces a (header-only) VCF.
process GERMLINE_HAPLOTYPECALLER {
    tag "${meta.id}"
    label 'process_medium'
    container params.containers.gatk

    input:
    tuple val(meta), path(normal_bam), path(normal_bai), path(windows_bed)
    path fasta
    path fasta_fai
    path fasta_dict

    output:
    tuple val(meta), path("${meta.id}.germline.vcf.gz"), path("${meta.id}.germline.vcf.gz.tbi"), emit: vcf

    script:
    """
    if [ -s ${windows_bed} ]; then
        intervals=${windows_bed}
    else
        intervals="\$(head -1 ${fasta_fai} | cut -f1):1-1"
    fi
    gatk --java-options "-Xmx${(task.memory.toGiga() * 0.8) as int}g" HaplotypeCaller \\
        -R ${fasta} \\
        -I ${normal_bam} \\
        -L "\${intervals}" \\
        --native-pair-hmm-threads ${task.cpus} \\
        -O ${meta.id}.germline.vcf.gz
    """

    stub:
    """
    touch ${meta.id}.germline.vcf.gz ${meta.id}.germline.vcf.gz.tbi
    """
}

process WHATSHAP_PHASE {
    tag "${meta.id}"
    label 'process_medium'
    container params.containers.whatshap

    input:
    tuple val(meta), path(vcf), path(tbi), path(tumor_bam), path(tumor_bai)
    path fasta
    path fasta_fai

    output:
    tuple val(meta), path("${meta.id}.phased.vcf"), emit: vcf

    script:
    """
    # Same settings as the Gadi PBS phasing script: indels phased too (the
    # default from WhatsHap 2.0, needed explicitly before it), and every
    # read used whatever its read group, since the VCF has one sample.
    whatshap --version
    whatshap phase \\
        --tag HP \\
        --indels \\
        --ignore-read-groups \\
        --reference ${fasta} \\
        -o ${meta.id}.phased.vcf \\
        ${vcf} ${tumor_bam}
    """

    stub:
    """
    touch ${meta.id}.phased.vcf
    """
}
