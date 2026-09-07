// Ported from: Sequenza_tools/sequenza_step2.sh, wrapping
// bin/merge-bin200-files.pl + bin/merge-header-bin200.pl (both now real:
// you sent their content; they're identical to what you sent, just
// living in this pipeline's bin/, which Nextflow auto-adds to PATH for
// every process, so they're invoked by bare name below exactly as
// sequenza_step2.sh does).
// https://github.com/JaydenBeckwith/Sequenza_tools
//
// Structural contract (unpacked from reading the two perl scripts, since
// they're only ever invoked with a sample id / sex, no explicit file
// list: they glob internally):
//   - merge-bin200-files.pl expects "./<sample>/<sample>.chrN_bin200.seqz.gz"
//     for each of chr1-22,X,(Y unless sex=='female': it skips Y itself
//     for female samples, so the WGS/female Y-exclusion logic already
//     lives in the script you sent, not duplicated here), builds
//     "./seqz.header" itself by taking the header line straight off
//     chr1's file (`zcat .../chr1..._bin200.seqz.gz | head -1 > seqz.header`
//    : it OVERWRITES seqz.header every time, so nothing needs to supply
//     one beforehand; the earlier draft of this module pre-copied a
//     reconstructed "standard" seqz.header in for exactly this step and
//     that was unnecessary once the real script turned out to build its
//     own), then concatenates all the per-chr files (minus their own
//     header lines, via `grep -v chromosome`) into
//     "./<sample>/<sample>_bin200_noheader.seqz".
//   - merge-header-bin200.pl then just `cat`s that shared "./seqz.header"
//     onto the noheader file and gzips it to
//     "./<sample>/<sample>_bin200.seqz.gz": the final merged file.
// Both scripts operate relative to the CURRENT directory ("."), so this
// module stages everything to match: the per-chromosome files go into a
// "<sample>/" subdirectory of the process work dir (mirroring
// sequenza_step2.sh's own SAMPLE_DIR), and the two perl scripts are
// invoked from the work dir itself, one level up.
//
// CONFIRMED (not a guess, now that the real script is in hand):
// merge-bin200-files.pl hardcodes "_bin200" in every filename it builds:
// it does not take or derive a bin size from anywhere. So the
// --sequenza_seq_type wes / 50bp-bin path this pipeline added as new
// plumbing on top of your WGS-only original is INCOMPATIBLE with this
// script as written: a WES run would produce "*_bin50.seqz.gz" files from
// SEQUENZA_BAM2SEQZ_BINNED that this script would never find (it always
// looks for "*_bin200.seqz.gz"). workflows/purity_ploidy.nf now errors
// out at launch if bin_size != 200 rather than letting that fail
// confusingly here: if you need WES support, this script would need a
// parametrized bin size, which isn't something to guess into your file.
//
// Same conda-env caveat as before: this process does NOT run in the
// sequenza Singularity container (Step 1's bam2seqz/binning does): it
// needs the r_sequenza conda env on whatever host executes it, carried
// over verbatim from sequenza_step2.sh even though the two scripts here
// are pure perl (system() calls to zcat/cat/gzip/grep) and arguably don't
// need it; kept as given rather than "optimised away" since it's not
// obvious from the shell script alone whether that matters.
// `params.sequenza_conda_sh` points at a path under a personal Gadi home
// directory (/home/562/jb1592/...) in the original: confirm it still
// resolves (or point it at wherever this conda env actually lives) before
// a real run.

process SEQUENZA_MERGE_BINS {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(binned_files), val(sex), val(bin_size)

    output:
    tuple val(meta), path("${meta.id}/${meta.id}_bin200.seqz.gz"), val(sex), emit: merged_seqz

    script:
    if (bin_size != 200)
        error "SEQUENZA_MERGE_BINS: bin_size=${bin_size} for ${meta.id}: merge-bin200-files.pl only ever looks for *_bin200.seqz.gz, see this module's header comment"
    """
    __conda_setup="\$('/home/562/jb1592/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
    if [ \$? -eq 0 ]; then
        eval "\$__conda_setup"
    else
        if [ -f "${params.sequenza_conda_sh}" ]; then
            source "${params.sequenza_conda_sh}"
        else
            export PATH="${params.sequenza_conda_bin_fallback}:\$PATH"
        fi
    fi
    unset __conda_setup
    conda activate ${params.sequenza_conda_env}
    module unload R 2>/dev/null || true
    export PATH="${params.sequenza_conda_bin_fallback}/../envs/${params.sequenza_conda_env}/bin:\$PATH"

    mkdir -p ${meta.id}
    cp ${binned_files.join(' ')} ${meta.id}/

    perl merge-bin200-files.pl ${meta.id} ${sex}
    perl merge-header-bin200.pl ${meta.id}

    if [ ! -f "${meta.id}/${meta.id}_bin${bin_size}.seqz.gz" ]; then
        echo "ERROR: expected merged file ${meta.id}/${meta.id}_bin${bin_size}.seqz.gz was not produced for ${meta.id}" >&2
        exit 1
    fi
    """
}
