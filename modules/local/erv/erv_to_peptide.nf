// TODO stub — the open design question flagged when this pipeline was
// first scoped, now backed by actual published precedent rather than just
// two hypothetical options: pVACseq/pVACtools were built for point
// mutations and fusions, not retroelement expression, so there is still no
// drop-in downstream tool that takes a Telescope report and returns
// pVACtools-ready peptides. Two real published approaches to choose
// between, not equivalent in effort or rigor:
//
//   1. Six-frame ORF translation of expressed ERV loci above an expression
//      threshold, then feed candidate peptides through pVACtools'
//      generate_protein_fasta as a custom peptide source (same shape as the
//      splicing branch's open decision). This is close to what ObsERV does
//      (Frontiers/npj Vaccines, 2025 — glioblastoma ERV vaccine pipeline):
//      TPM > 1 as the expression cutoff, 8-11mers scored for MHC-I ligands
//      and 15mers for MHC-II, then 27mer fragments ranked by an
//      "ObsERV score" that prioritizes fragments carrying both. ObsERV
//      itself is built on gEVE + RSEM, not Telescope, and isn't published
//      as a standalone installable tool (no confirmed GitHub/PyPI/bioconda
//      package as of this research pass) — so this option means adapting
//      its published scoring *logic* on top of this pipeline's own
//      Telescope counts and existing pVACtools/NetMHCpan setup, not
//      installing ObsERV itself.
//   2. A dedicated de novo transcript-assembly pipeline in the style of
//      Attig et al. / TEProF2 (Nature Genetics, 2023,
//      doi.org/10.1038/s41588-023-01349-3) — STAR + Cufflinks/TopHat to
//      assemble TE-chimeric transcripts, CPC2 + Pfam for coding-potential
//      filtering, NetMHCpan-4.0 for binding, with mass-spec (MaxQuant)
//      validation in the original paper. TEProF2 itself IS a real,
//      published, installable pipeline (Zenodo DOIs in the paper), so this
//      option is closer to "install and wire in a tool" than option 1 — but
//      it's built around older assembly tools (Cufflinks/TopHat), targets
//      TE-driven chimeric transcripts specifically (a different, narrower
//      phenomenon than "any expressed ERV locus" that Telescope quantifies),
//      and is considerably more infrastructure to stand up than option 1.
//
// Not picking one without your input — this remains the least standardized
// of the three new branches by a wide margin, and unlike the fusion
// branch's AGFusion stub (blocked on running agfusion-build once, not on a
// design choice), this is a genuine fork in approach with no dominant
// convention to default to.

process ERV_TO_PEPTIDE {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_medium'

    input:
    tuple val(meta), path(telescope_report)

    output:
    tuple val(meta), path("${meta.id}_${meta.timepoint}.erv_peptides.fasta"), emit: peptide_fasta

    script:
    """
    echo "ERV_TO_PEPTIDE is a design stub — see module comments for the six-frame-ORF (ObsERV-style) vs TEProF2-style dedicated-pipeline decision" >&2
    exit 1
    """
}
