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

    output:
    tuple val(meta), val('arriba'), path("${meta.id}_${meta.timepoint}.arriba_fusions.tsv"), emit: fusions
    path "${meta.id}_${meta.timepoint}.arriba_fusions.discarded.tsv", emit: discarded

    script:
    def blacklist_opt = params.arriba_blacklist ? "-b ${params.arriba_blacklist}" : ''
    def known_fusions_opt = params.arriba_known_fusions ? "-k ${params.arriba_known_fusions}" : ''
    def protein_domains_opt = params.arriba_protein_domains ? "-p ${params.arriba_protein_domains}" : ''
    """
    if [ -z "${params.arriba_blacklist}" ]; then
        echo "WARNING: --arriba_blacklist not set: running Arriba without a blacklist, which Arriba's own docs say increases false-positive fusion calls. Run Arriba's download_references.sh for your assembly+annotation build and set --arriba_blacklist before trusting real results." >&2
    fi

    arriba \\
        -x ${bam} \\
        -o ${meta.id}_${meta.timepoint}.arriba_fusions.tsv \\
        -O ${meta.id}_${meta.timepoint}.arriba_fusions.discarded.tsv \\
        -a ${fasta} \\
        -g ${gtf} \\
        ${blacklist_opt} ${known_fusions_opt} ${protein_domains_opt}
    """
}
