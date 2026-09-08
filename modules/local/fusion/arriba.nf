// New: Arriba fusion calling. CLI flags copied verbatim from Arriba's own
// documented quickstart (github.com/suhrig/arriba wiki, "Command-line
// options"): -x main alignment, -o/-O passed/discarded fusions, -a
// genome fasta, -g GTF, and the optional-but-strongly-recommended -b
// blacklist / -k known-fusions / -p protein-domains references.
//
// The blacklist/known-fusions/protein-domains files are genome+annotation-
// build-specific and distributed separately from the arriba binary itself
// (via its own bundled `download_references.sh ASSEMBLY+ANNOTATION`
// script, e.g. `download_references.sh hg38+GENCODE46` to match this
// pipeline's GENCODE v46 GTF): not something to fetch or fabricate here.
// params.arriba_blacklist/arriba_known_fusions/arriba_protein_domains
// default to null; -b/-k/-p are only added to the command when set, so
// Arriba still runs without them (Arriba documents this as supported, just
// lower-specificity: more false-positive fusions from read-through
// transcription, homology, etc.) rather than hard-erroring. Strongly
// recommended to fetch and set these before trusting real results.

process ARRIBA {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_medium'
    container params.containers.arriba
    publishDir "${params.outdir}/${meta.id}/fusion/${meta.timepoint}/arriba", mode: 'copy'

    input:
    tuple val(meta), path(bam)
    path fasta
    path gtf
    path blacklist, stageAs: 'blacklist/*'
    path known_fusions, stageAs: 'known_fusions/*'
    path protein_domains, stageAs: 'protein_domains/*'

    output:
    tuple val(meta), val('arriba'), path("${meta.id}_${meta.timepoint}.arriba_fusions.tsv"), emit: fusions
    path "${meta.id}_${meta.timepoint}.arriba_fusions.discarded.tsv", emit: discarded

    script:
    def known_fusions_opt = known_fusions ? "-k '${known_fusions}'" : ''
    def protein_domains_opt = protein_domains ? "-p '${protein_domains}'" : ''
    """
    arriba \\
        -x '${bam}' \\
        -o ${meta.id}_${meta.timepoint}.arriba_fusions.tsv \\
        -O ${meta.id}_${meta.timepoint}.arriba_fusions.discarded.tsv \\
        -a '${fasta}' \\
        -g '${gtf}' \\
        -b '${blacklist}' ${known_fusions_opt} ${protein_domains_opt}
    """

    stub:
    """
    printf '#gene1\\tgene2\\tgene_id1\\tgene_id2\\tbreakpoint1\\tbreakpoint2\\tsplit_reads1\\tsplit_reads2\\tdiscordant_mates\\n' > ${meta.id}_${meta.timepoint}.arriba_fusions.tsv
    printf '#gene1\\tgene2\\tgene_id1\\tgene_id2\\tbreakpoint1\\tbreakpoint2\\tsplit_reads1\\tsplit_reads2\\tdiscordant_mates\\n' > ${meta.id}_${meta.timepoint}.arriba_fusions.discarded.tsv
    """
}
