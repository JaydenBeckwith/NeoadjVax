// Ported from: Sequenza_tools/sequenza_step2.sh
// https://github.com/JaydenBeckwith/Sequenza_tools
//
// Merges one sample's 24 per-chromosome binned seqz files (from
// SEQUENZA_BAM2SEQZ_BINNED) into a single file, using merge-bin200-files.pl
// + merge-header-bin200.pl. The original ran these from a shared prep
// directory, invoking each script with just (sample, sex) / (sample) and
// letting the script glob "${sample}.chr*_bin${size}.seqz.gz" itself in
// its own cwd — so the structural contract this module replicates is: put
// all of one sample's binned files, named exactly as they came out of
// bam2seqz/binning, into a "<sample>/" subdirectory, put seqz.header next
// to it (not inside it — matches the original's `cp seqz.header ${INDIR}/`
// one level up from the per-sample subdirs), then call the two perl
// scripts from that parent level.
//
// NOT YET REAL: merge-bin200-files.pl and merge-header-bin200.pl
// themselves haven't been ported — you sent the four PBS wrapper scripts,
// not these two perl scripts they call out to (originally read from
// `/scratch/jo11/neoadjuvant` on Gadi). Send them (paste or upload) and
// they drop straight into bin/ under those exact names — Nextflow adds a
// pipeline's bin/ directory to PATH for every process automatically, so
// no other change would be needed here. Until then this process fails
// fast with a clear message rather than a confusing "command not found".
//
// assets/seqz.header here is the STANDARD sequenza-utils seqz output
// header (chromosome/position/base.ref/depth.normal/depth.tumor/
// depth.ratio/Af/Bf/zygosity.normal/GC.percent/good.reads/AB.normal/
// AB.tumor/tumor.strand) reconstructed from sequenza-utils' own public
// format docs, NOT copied from your repo — you didn't send this file
// either. If your actual seqz.header differs (extra/reordered columns,
// a different tool version), replace assets/seqz.header with your real
// one; this is a much smaller guess than the perl scripts since it's a
// documented, stable public format, but still flagged rather than assumed.
//
// One more thing that can't be confirmed without the actual perl scripts:
// merge-bin200-files.pl's own name bakes in "200", which is Jayden's WGS
// bin size — it's not visible from sequenza_step2.sh's shell orchestration
// alone whether the script hardcodes a 200bp assumption internally or
// genuinely takes the bin size from the input filenames it globs. Treat
// --sequenza_bin_wes (50) as UNVERIFIED against this script until confirmed
// — the original was only ever run WGS-only.
//
// The env setup below (conda activate r_sequenza / module unload R) is
// copied verbatim from sequenza_step2.sh even though this step is mostly
// perl — kept exactly as given rather than "optimised away", since it's
// not obvious from the script alone whether the merge scripts shell out to
// R internally. This also means, unlike every other module in this repo,
// Step 2 and Step 4 (sequenza_fit.nf) do NOT run in the sequenza
// Singularity container — they need that r_sequenza conda env to actually
// exist on whatever host executes them, same real constraint the original
// scripts had. `params.sequenza_conda_sh` points at a path under
// `/home/562/jb1592/...` in the original — that's a personal home
// directory, not a shared/project location, so confirm it still resolves
// (or point it at wherever this conda env actually lives) before a real run.

process SEQUENZA_MERGE_BINS {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(binned_files), val(sex), val(bin_size)
    path seqz_header

    output:
    tuple val(meta), path("${meta.id}/${meta.id}_bin${bin_size}.seqz.gz"), val(sex), emit: merged_seqz

    script:
    """
    if ! command -v merge-bin200-files.pl >/dev/null 2>&1 && [ ! -f "\$(dirname "\$0")/merge-bin200-files.pl" ]; then
        echo "ERROR: merge-bin200-files.pl not found on PATH / in bin/ — this script hasn't been ported yet, see modules/local/purity_ploidy/sequenza_merge_bins.nf header comment. Send it and it drops straight into bin/." >&2
        exit 1
    fi
    if ! command -v merge-header-bin200.pl >/dev/null 2>&1 && [ ! -f "\$(dirname "\$0")/merge-header-bin200.pl" ]; then
        echo "ERROR: merge-header-bin200.pl not found on PATH / in bin/ — same as above, not ported yet." >&2
        exit 1
    fi

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
    cp ${seqz_header} seqz.header

    perl merge-bin200-files.pl ${meta.id} ${sex}
    perl merge-header-bin200.pl ${meta.id}

    if [ ! -f "${meta.id}/${meta.id}_bin${bin_size}.seqz.gz" ]; then
        echo "ERROR: expected merged file ${meta.id}/${meta.id}_bin${bin_size}.seqz.gz was not produced for ${meta.id}" >&2
        exit 1
    fi
    """
}
