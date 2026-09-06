// Telescope requires its input SAM/BAM to be COLLATED (all alignments for a
// read pair adjacent in the file), not coordinate-sorted — confirmed from
// Telescope's own docs (github.com/mlbendall/telescope): "must be collated
// so that all alignments for a read pair appear sequentially in the
// file... coordinate-sorted BAMs do not work," recommending `samtools
// collate` as the faster alternative to `samtools sort -n`. STAR's own
// "Unsorted" BAM output is not documented to guarantee multi-mapped mate
// alignments stay strictly adjacent, so this runs collate explicitly rather
// than assuming ERV_STAR_ALIGN's output order is already good enough —
// cheap, and removes a way for Telescope to silently misassign or undercount.

process SAMTOOLS_COLLATE {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_medium'
    container params.containers.samtools

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}_${meta.timepoint}.erv.collated.bam"), emit: bam

    script:
    """
    samtools collate -@ ${task.cpus} -o ${meta.id}_${meta.timepoint}.erv.collated.bam ${bam}
    """
}
