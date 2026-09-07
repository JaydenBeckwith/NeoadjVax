// ERV/transposable-element expression quantification. Telescope chosen over
// alternatives (TEtranscripts/TElocal give TE-family- or locus-level counts
// for differential expression but aren't built with per-locus reassignment
// confidence in mind the way immunopeptidome work wants; ERVmap needs its
// own dedicated ERV-only alignment/reference rather than a standard genome
// BAM): confirmed still the standard choice for this specific use case by
// two independently published 2025/2026 HERV-in-cancer studies that both
// used Telescope on a standard STAR BAM with relaxed multimapping, not a
// guess. Fed by this branch's own ERV_STAR_ALIGN + SAMTOOLS_COLLATE steps
// (see erv_star_align.nf / collate_bam.nf): not RNA_VARIANT_CALLING's BAM,
// which discards the multi-mapping reads ERV loci need, and not a
// coordinate-sorted BAM either (see collate_bam.nf for why Telescope
// specifically rejects that).
//
// Annotation (--erv_annotation_gtf): use one of Telescope's own official
// GRCh38 builds from github.com/mlbendall/telescope_annotation_db/tree/master/builds
// rather than a generic RepeatMasker dump: `retro.hg38.v1` (curated
// retroelement build, what Telescope's own tutorial and most published
// Telescope analyses use) is the recommended default; `HERV_rmsk.hg38.v2`
// (broader, RepeatMasker-derived, HERV-only) is the alternative if you want
// wider HERV coverage at the cost of some annotation noise. Not fetched or
// bundled here: genome-build-specific, same category as the fusion
// branch's Arriba reference files.
//
// --max_iter 1000 (vs. Telescope's own default of 100, --telescope_max_iter
// param): matches the value used in a recently published cutaneous-melanoma
// HERV profiling study (Frontiers in Oncology, 2026) doing essentially the
// same STAR-relaxed-multimap + Telescope quantification as here: not the
// tool's default, but not invented either. --theta_prior is left at
// Telescope's own default (200000); that same melanoma study passed it
// explicitly but at its default value, so there's nothing to override.
//
// Container tag: the previous `1.0.3--py38h24c8ff8_1` tag did not exist on
// depot.galaxyproject.org (confirmed 404): fixed to `1.0.3--py38h8e05983_5`,
// confirmed against bioconda's own published build list for telescope 1.0.3
// (Linux-64) and cross-checked with a partial fetch against
// depot.galaxyproject.org that returned real multi-megabyte binary content
// for that tag (vs. a clean 404 for the old one): see nextflow.config.

process TELESCOPE_QUANT {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_medium'
    container params.containers.telescope
    publishDir "${params.outdir}/${meta.id}/erv/${meta.timepoint}", mode: 'copy'

    input:
    tuple val(meta), path(collated_bam)
    path erv_annotation_gtf

    output:
    tuple val(meta), path("${meta.id}_${meta.timepoint}-telescope_report.tsv"), emit: report

    script:
    """
    telescope assign \\
        ${collated_bam} \\
        ${erv_annotation_gtf} \\
        --outdir . \\
        --exp_tag ${meta.id}_${meta.timepoint} \\
        --max_iter ${params.telescope_max_iter}
    """
}
