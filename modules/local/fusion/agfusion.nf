// AGFusion 1.4.3 on STAR-Fusion's predictions, as run_agfusion.pbs does:
//   agfusion batch -a starfusion -db agfusion.homo_sapiens.<release>.db
//       --middlestar --noncanonical
//   --middlestar marks the junction in each fusion protein; pVACfuse needs
//     it to find where fusion peptides differ from the partner proteins.
//   --noncanonical (params.agfusion_noncanonical) annotates every
//     transcript pair, not only the canonical ones.
//
// AGFusion reads Ensembl gene models through pyensembl, which must find an
// already-indexed cache for the same release as the DB (compute nodes are
// offline, so it cannot download one). main.nf checks this at launch.
// params.agfusion_pythonpath is the directory holding the numpy/datacache
// sitecustomize.py shim the image needs; Python imports sitecustomize from
// PYTHONPATH at start-up. Both are exported inside the container, which
// is why the script sets them rather than SINGULARITYENV_*.
//
// The cache is used in place (not staged): it is read-only once indexed,
// and copying it per patient would cost several GB each time.

process AGFUSION {
    tag "${meta.id}"
    label 'process_agfusion'
    container params.containers.agfusion
    publishDir "${params.outdir}/${meta.id}/fusion", mode: 'copy'

    input:
    tuple val(meta), path(predictions)
    path agfusion_db

    output:
    tuple val(meta), path("agfusion"), emit: annotated

    script:
    def noncanonical = params.agfusion_noncanonical ? '--noncanonical' : ''
    def pythonpath = params.agfusion_pythonpath ? "export PYTHONPATH='${params.agfusion_pythonpath}'\${PYTHONPATH:+:\$PYTHONPATH}" : ''
    """
    export PYENSEMBL_CACHE_DIR='${params.agfusion_pyensembl_cache}'
    ${pythonpath}

    # A STAR-Fusion table with no calls is a valid negative result
    if [ "\$(grep -vc '^#' ${predictions} || true)" -eq 0 ]; then
        mkdir -p agfusion
        echo "STAR-Fusion reported no fusions for ${meta.id}" > agfusion/no_fusions.txt
        exit 0
    fi

    agfusion batch \\
        -f ${predictions} \\
        -a starfusion \\
        -db ${agfusion_db} \\
        -o agfusion \\
        --middlestar \\
        ${noncanonical} \\
        ${params.agfusion_extra_args}
    """

    stub:
    """
    mkdir -p agfusion
    printf 'Stub run: no biological analysis performed.\\n' > agfusion/stub.log
    """
}
