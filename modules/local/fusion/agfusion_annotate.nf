// TODO stub — annotates STAR-Fusion predictions into the fusion-transcript
// format pVACfuse needs (a "AGFusion"-annotated directory per fusion event:
// context sequence + breakpoint + reading frame). AGFusion is the tool
// pVACtools' own docs point at for this, but its local annotation database
// (agfusion-build) needs to be built once against our GENCODE v46 GTF
// before this can run for real — not yet done, so left as a stub rather
// than guessing at untested paths/versions.

process AGFUSION_ANNOTATE {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_medium'

    input:
    tuple val(meta), path(star_fusion_tsv)

    output:
    tuple val(meta), path("agfusion_out"), emit: annotated

    script:
    """
    # TODO — unverified. Rough shape per pVACtools/AGFusion docs:
    #   agfusion batch -f ${star_fusion_tsv} -a starfusion \\
    #       -db agfusion.homo_sapiens.87.db -o agfusion_out
    echo "AGFUSION_ANNOTATE is a stub — needs agfusion-build run against our GTF first" >&2
    exit 1
    """
}
