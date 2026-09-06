// New — ERV/transposable-element quantification needs its own STAR pass,
// same shape of problem already flagged for Arriba (see
// modules/local/fusion/arriba_align.nf): this pipeline's RNA branch aligns
// with --outFilterMultimapNmax 2, which throws away exactly the
// multi-mapping reads Telescope needs (most TE/ERV loci are repetitive and
// reads from them multi-map by nature).
//
// --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 is the setting
// confirmed from two independent, converging sources rather than guessed:
// (1) TEtranscripts' own README (github.com/mhammell-laboratory/TEtranscripts,
// the standard companion aligner recipe for TE-aware quantification, tested
// by its authors at 100/100 and explicitly recommending the same value for
// both flags), and (2) STAR's own maintainer, in
// github.com/alexdobin/STAR/discussions/1410, giving the same 100/100
// pairing for retaining TE-derived multimappers. Telescope's own project
// (mlbendall/telescope_tutorial) documents its alignment recipe around
// Bowtie2 (-k 100 --very-sensitive-local), not STAR, so there's no
// STAR-specific number from Telescope itself to defer to instead — 100/100
// is the closest match to that same "keep up to ~100 mapping loci per read"
// intent, translated to STAR's own flags.
//
// Reuses this pipeline's existing GRCh38/GENCODE STAR index rather than a
// second reference set — same reasoning as Arriba's STAR pass.

process ERV_STAR_ALIGN {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_high'
    container params.containers.star

    input:
    tuple val(meta), path(r1), path(r2)
    path star_index_dir

    output:
    tuple val(meta), path("${meta.id}_${meta.timepoint}.erv.Aligned.out.bam"), emit: bam

    script:
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${star_index_dir} \\
        --genomeLoad NoSharedMemory \\
        --readFilesIn ${r1} ${r2} \\
        --readFilesCommand zcat \\
        --outSAMtype BAM Unsorted \\
        --outSAMunmapped Within \\
        --outFilterMultimapNmax 100 \\
        --winAnchorMultimapNmax 100 \\
        --outFileNamePrefix ${meta.id}_${meta.timepoint}.erv.
    """
}
