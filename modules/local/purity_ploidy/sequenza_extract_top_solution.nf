// TODO stub — the one genuinely-blocked piece left in the purity/ploidy
// branch. SEQUENZA_BAM2SEQZ_BINNED, SEQUENZA_MERGE_BINS and SEQUENZA_FIT
// are all real, faithful ports of Sequenza_tools' four PBS scripts; this
// step is new because Nextflow needs a (meta, purity, ploidy) value out of
// SEQUENZA_FIT's raw output to hand to HLA_LOH, and nothing in the four
// scripts you sent does that extraction — run-sequenza.R (not yet sent)
// presumably picks a "top solution" internally (that's what NeoadjLOH's
// sequenza_top_solutions_summary.csv sounds like), but its output
// filenames/columns aren't known here, so guessing a parse would risk
// silently picking the wrong solution or the wrong column.
//
// Once run-sequenza.R (or whatever separately builds the top-solution
// summary, if that's a distinct script) is available, replace this with
// real logic reading SEQUENZA_FIT's sequenza_raw_out/ directory — most
// likely something under the standard sequenza R package's own output
// naming (e.g. a *_confints_CP.txt / *_alternative_solutions.txt style
// file with cellularity + ploidy.estimate columns) if run-sequenza.R
// follows that convention, but that's exactly the kind of assumption not
// worth encoding without seeing the real output.
//
// Segment-level copy-number classification (HETEROZYGOUS/LOH_DELETION/
// CNN_LOH/LOH_AMPLIFIED/biallelic-loss) is a separate downstream step, not
// part of this extraction — bin/classify_cn_state.py already has a draft
// of that rule (reconstructed from conversation, not your real R
// function; see that file's own docstring) for whenever you're ready to
// wire it in, most likely against a `*_segments.txt`-shaped file that
// run-sequenza.R also produces.



process SEQUENZA_EXTRACT_TOP_SOLUTION {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), val(sex), path(raw_results)

    output:
    tuple val(meta), env(PURITY), env(PLOIDY), emit: purity_ploidy

    script:
    """
    echo "ERROR: SEQUENZA_EXTRACT_TOP_SOLUTION is a design stub for ${meta.id} — see modules/local/purity_ploidy/sequenza_extract_top_solution.nf header comment. Needs run-sequenza.R's real output format to parse purity/ploidy out of ${raw_results}." >&2
    exit 1
    """
}
