// Ported from: Sequenza_tools/sequenza_step4.sh
// https://github.com/JaydenBeckwith/Sequenza_tools
//
// Decompresses the merged seqz file (R's vroom-based gzip reading was
// unreliable under R 4.0 per the original script's own comment — kept the
// same decompress-then-cleanup approach rather than trying to feed R the
// .gz directly) and runs run-sequenza.R sample sex, matching the original
// invocation exactly. The original cd'd into a per-sample SAMPLE_OUTDIR
// before running Rscript and let it write into cwd; a Nextflow process's
// own work directory already plays that role, so there's no separate cd —
// whatever run-sequenza.R writes lands directly in this task's outputs.
//
// NOT YET REAL: run-sequenza.R itself hasn't been ported — you sent the
// four PBS wrapper scripts, not this R script (originally read from
// `/scratch/jo11/neoadjuvant`). It presumably wraps the sequenza R
// package's sequenza.extract() / sequenza.fit() / sequenza.results() and
// picks a "top solution" feeding into something like NeoadjLOH's
// sequenza_top_solutions_summary.csv (sample,purity,ploidy) — but the
// exact output filenames/columns, gamma/kmin/other fit parameters, and
// whether the top-solution-per-sample CSV comes out of this script
// directly or a separate aggregation step are all unknown without seeing
// it. Send it (and whatever produces the aggregated top-solution CSV, if
// that's a different script) and it drops into bin/ as
// bin/run-sequenza.R; modules/local/purity_ploidy/sequenza_extract_top_solution.nf
// is the deliberately-stubbed next step that turns this process's raw
// output into the (meta, purity, ploidy) triple HLA_LOH needs — see that
// file for what's blocked on it.
//
// Same conda-env caveat as SEQUENZA_MERGE_BINS: this process does NOT run
// in the sequenza Singularity container, it needs the r_sequenza conda env
// on whatever host executes it — confirm params.sequenza_conda_sh still
// resolves (the original path is under a personal Gadi home directory).

process SEQUENZA_FIT {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(merged_seqz_gz), val(sex)

    output:
    tuple val(meta), val(sex), path("sequenza_raw_out"), emit: raw_results

    script:
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

    RUN_SEQUENZA_R="\$(command -v run-sequenza.R || true)"
    if [ -z "\${RUN_SEQUENZA_R}" ]; then
        echo "ERROR: run-sequenza.R not found on PATH / in bin/ — this script hasn't been ported yet, see modules/local/purity_ploidy/sequenza_fit.nf header comment. Send it and it drops straight into bin/." >&2
        exit 1
    fi

    mkdir -p sequenza_raw_out
    SEQZ_UNCOMPRESSED="sequenza_raw_out/${meta.id}_bin.seqz"
    gunzip -c ${merged_seqz_gz} > "\${SEQZ_UNCOMPRESSED}"

    cd sequenza_raw_out
    Rscript "\${RUN_SEQUENZA_R}" "${meta.id}_bin.seqz" ${meta.id} ${sex}
    rm -f "${meta.id}_bin.seqz"
    """
}
