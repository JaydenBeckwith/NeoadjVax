// TODO stub: annotates a caller's fusion predictions into the
// fusion-transcript format pVACfuse needs (a "AGFusion"-annotated
// directory per fusion event: context sequence + breakpoint + reading
// frame). AGFusion is the tool pVACtools' own docs point at for this, and
// its `agfusion batch -f <format>` command documents support for both
// callers now feeding it (`-f starfusion` / `-f arriba`), so this module
// is written caller-aware even though it's still a stub: the blocker
// hasn't changed: AGFusion's local annotation database (agfusion-build)
// needs to be built once against our GENCODE v46 GTF before this can run
// for real, and that hasn't been done. STAR-Fusion and Arriba themselves
// (the actual fusion-CALLING step upstream of this) are real as of this
// scaffold: see star_fusion.nf / arriba_align.nf / arriba.nf.

process AGFUSION_ANNOTATE {
    tag "${meta.id}:${meta.timepoint}:${caller}"
    label 'process_medium'

    input:
    tuple val(meta), val(caller), path(fusion_file)

    output:
    tuple val(meta), val(caller), path("agfusion_out"), emit: annotated

    script:
    """
    # TODO: unverified. Rough shape per pVACtools/AGFusion docs:
    #   agfusion batch -f ${fusion_file} -a ${caller} \\
    #       -db agfusion.homo_sapiens.87.db -o agfusion_out
    echo "AGFUSION_ANNOTATE is a stub for ${meta.id}:${meta.timepoint} (${caller}): needs agfusion-build run against our GTF first" >&2
    exit 1
    """
}
