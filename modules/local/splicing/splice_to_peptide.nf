// TODO stub — the genuinely unsolved step in the splicing branch.
// Unlike DNA/RNA SNV->pVACseq, there's no single standard tool that turns
// a novel/altered splice junction directly into a candidate peptide FASTA.
// Two real options to evaluate before implementing this for real:
//
//   1. NeoSplice (Zhang lab) — purpose-built for splice-neoantigen
//      discovery from RNA-seq, but has a narrower/less-maintained toolchain.
//   2. Custom ORF translation — take novel junctions (from SpliceAI deltas
//      above + observed junctions in the STAR-aligned RNA BAM, e.g. via
//      SJ.out.tab or a tool like LeafCutter/MAJIQ for the isoform-switch
//      side), splice them into transcript sequences, translate in-frame
//      from the annotated start codon, and stop at the first novel stop or
//      frameshift — producing a mutant/wildtype peptide pair pVACseq's
//      generate_protein_fasta input format can accept directly.
//
// Given [[neoadjuvant-splicing]] already has IsoformSwitchAnalyzeR output
// (isoform switching across treatment arms) and PRADO scRNA integration in
// progress, option 2 built on top of that existing analysis is probably
// the more direct path — needs a design conversation before writing real
// code here, not a solo implementation guess.

process SPLICE_TO_PEPTIDE {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_medium'

    input:
    tuple val(meta), path(spliceai_filtered_vcf)

    output:
    tuple val(meta), path("${meta.id}_${meta.timepoint}.splice_peptides.fasta"), emit: peptide_fasta

    script:
    """
    echo "SPLICE_TO_PEPTIDE is a design stub — see module comments for the NeoSplice vs custom-ORF-translation decision" >&2
    exit 1
    """
}
