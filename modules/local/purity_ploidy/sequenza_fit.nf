// Ported from: Sequenza_tools/sequenza_step4.sh
// https://github.com/JaydenBeckwith/Sequenza_tools
//
// Decompresses the merged seqz file (R's vroom-based gzip reading was
// unreliable under R 4.0 per the original script's own comment: kept the
// same decompress-then-cleanup approach rather than trying to feed R the
// .gz directly) and runs run-sequenza.R sample sex, matching the original
// invocation exactly. The original cd'd into a per-sample SAMPLE_OUTDIR
// before running Rscript and let it write into cwd; a Nextflow process's
// own work directory already plays that role, so there's no separate cd:
// whatever run-sequenza.R writes lands directly in this task's outputs.
//
// bin/run-sequenza.R is now your real script, unmodified: sequenza.extract()
// (normalization.method="median", kmin=500, gamma=100) -> sequenza.fit()
// (ratio.priority=FALSE) -> sequenza.results(), with the chromosome list
// and female=TRUE/FALSE + XY mapping branched on args[3] exactly as
// written. Nextflow auto-adds bin/ to PATH, so `command -v run-sequenza.R`
// below finds it the same way it would find any other bin/ script.
//
// IMPORTANT: the R script's own branch is `if (args[3] == "male") ... else
// <female path>`: literally any value other than the exact string "male"
// takes the FEMALE path (wrong chromosome list / X-Y handling for a male
// sample if the sex value isn't spelled exactly right). This pipeline's
// --sequenza_gender_csv lookup lowercases whatever's in column 7, so it
// works correctly if your CSV spells out "Male"/"Female" (any case):
// but if it instead uses single-letter codes ("M"/"F"), those lowercase to
// "m"/"f", neither of which equals "male", and every sample would
// silently run down the female path. Worth confirming your CSV's actual
// gender values before trusting a real run: see
// assets/samplesheet_schema.md.
//
// run-sequenza.R writes into sequenza.results()'s out.dir =
// "<sample_name>_OUTPUT" (a subdirectory it creates itself, using
// sample_name = args[2] verbatim, not the hyphen-sanitized job_name
// variable): relative to R's own working directory, which is why this
// process's script cd's into sequenza_raw_out/ before invoking Rscript, so
// that ends up at sequenza_raw_out/<meta.id>_OUTPUT/. See
// modules/local/purity_ploidy/sequenza_extract_top_solution.nf for what's
// read out of that directory next, and what's still unverified about it
// (the standard sequenza R package output convention it targets, not
// something confirmed against a real run of this exact script).
//
// Same conda-env caveat as SEQUENZA_MERGE_BINS: this process does NOT run
// in the sequenza Singularity container, it needs the r_sequenza conda env
// on whatever host executes it: confirm params.sequenza_conda_sh still
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
        echo "ERROR: run-sequenza.R not found on PATH / in bin/ for ${meta.id}: it should be at bin/run-sequenza.R in this pipeline and Nextflow auto-adds bin/ to PATH; check it's actually present/synced." >&2
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
