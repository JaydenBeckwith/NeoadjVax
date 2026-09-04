// New module — HLA class I typing. Not present in either original script:
// the DNA script took --hla as a manually supplied list. Runs on the
// NORMAL (germline) DNA BAM, deliberately not tumor, to avoid bias from
// tumor HLA-LOH — which [[loh-analysis]] is already tracking as its own
// workstream, so typing off tumor reads here would double-count that
// effect into the neoantigen calls.

process HLA_TYPING_OPTITYPE {
    tag "${meta.id}"
    label 'process_medium'
    container params.containers.optitype
    publishDir "${params.outdir}/${meta.id}/hla_typing", mode: 'copy'

    input:
    tuple val(meta), path(normal_bam), path(normal_bai)

    output:
    tuple val(meta), env(HLA_ALLELES), emit: hla_alleles
    path "optitype_out/*", emit: results

    script:
    """
    samtools bam2fq -1 normal_R1.fastq -2 normal_R2.fastq -n ${normal_bam}

    OptiTypePipeline.py \\
        -i normal_R1.fastq normal_R2.fastq \\
        --dna \\
        -v \\
        -o optitype_out

    result_tsv=\$(find optitype_out -name '*_result.tsv' | head -n1)
    HLA_ALLELES=\$(parse_optitype_result.py "\$result_tsv")
    echo "[INFO] \${HLA_ALLELES}"
    """
}
