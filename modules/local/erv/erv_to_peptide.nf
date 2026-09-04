// TODO stub — the open design question flagged when this pipeline was
// first scoped: pVACseq/pVACtools were built for point mutations and
// fusions, not retroelement expression, so there isn't a drop-in
// downstream tool here. Two paths to choose between before implementing:
//
//   1. Treat expressed ERV loci as novel ORFs — six-frame-translate the
//      Telescope-quantified locus (using its genomic coordinates) above an
//      expression threshold, filter for a plausible ORF, and feed the
//      resulting peptides through pVACtools' generate_protein_fasta as if
//      they were any other custom peptide source. Simple, reuses existing
//      infrastructure, but the expression-threshold and ORF-calling logic
//      would be hand-rolled and unvalidated.
//   2. A dedicated scorer, e.g. published ERV-neoantigen pipelines
//      (e.g. the approach in Smith et al./similar TE-neoantigen papers) —
//      more defensible biologically, more work to stand up, and less
//      actively maintained tooling to build on than option 1.
//
// Not picking one without your input — this is the least standardized of
// the three new branches by a wide margin.

process ERV_TO_PEPTIDE {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_medium'

    input:
    tuple val(meta), path(telescope_report)

    output:
    tuple val(meta), path("${meta.id}_${meta.timepoint}.erv_peptides.fasta"), emit: peptide_fasta

    script:
    """
    echo "ERV_TO_PEPTIDE is a design stub — see module comments for the novel-ORF vs dedicated-scorer decision" >&2
    exit 1
    """
}
