// New — ERV/transposable-element expression quantification. Telescope
// chosen over ERVmap because it works directly off a standard genome BAM
// (no separate ERV-specific alignment needed) and has an actively
// maintained retro.gtf annotation for GRCh38. Needs the RNA branch's
// STAR-aligned BAM (before BQSR — Telescope wants the raw multi-mapped
// alignments, not a recalibrated/duplicate-marked BAM) with
// --outFilterMultimapNmax relaxed beyond the RNA variant-calling branch's
// current setting of 2, since ERV quantification specifically needs the
// multi-mapping reads that setting discards. That means this can't just
// reuse RNA_VARIANT_CALLING's BAM output as-is — flagged in the workflow.

process TELESCOPE_QUANT {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_medium'
    container params.containers.telescope
    publishDir "${params.outdir}/${meta.id}/erv/${meta.timepoint}", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path erv_annotation_gtf

    output:
    tuple val(meta), path("${meta.id}_${meta.timepoint}-telescope_report.tsv"), emit: report

    script:
    """
    telescope assign \\
        ${bam} \\
        ${erv_annotation_gtf} \\
        --outdir . \\
        --exp_tag ${meta.id}_${meta.timepoint}
    """
}
