// New: Arriba needs its own STAR pass, not RNA_VARIANT_CALLING's BAM.
// Same reasoning already flagged for the ERV/Telescope branch: Arriba's
// own quickstart (github.com/suhrig/arriba, "Quick start") recommends
// relaxed multimapping (--outFilterMultimapNmax 50, vs. this pipeline's
// RNA branch using 2) plus a specific set of chimeric-alignment flags
// STAR needs turned on for Arriba to see split/discordant reads at all
// (--chimSegmentMin etc.): flags copied verbatim from Arriba's own
// documented command, not invented. Reuses this pipeline's existing
// params.star_index_dir/genome_fasta/gtf (GRCh38/GENCODE) rather than a
// second reference set: Arriba just needs its `-a`/`-g` inputs to match
// whatever built the index, which they already do.

process ARRIBA_STAR_ALIGN {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_high'
    container params.containers.star

    input:
    tuple val(meta), path(r1), path(r2)
    path star_index_dir

    output:
    tuple val(meta), path("${meta.id}_${meta.timepoint}.arriba.Aligned.out.bam"), emit: bam

    script:
    def read_command = r1.name.endsWith('.gz') ? 'zcat' : 'cat'
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir '${star_index_dir}' \\
        --genomeLoad NoSharedMemory \\
        --readFilesIn '${r1}' '${r2}' \\
        --readFilesCommand ${read_command} \\
        --outSAMtype BAM Unsorted \\
        --outSAMunmapped Within \\
        --outBAMcompression 0 \\
        --outFilterMultimapNmax 50 \\
        --peOverlapNbasesMin 10 \\
        --alignSplicedMateMapLminOverLmate 0.5 \\
        --alignSJstitchMismatchNmax 5 -1 5 5 \\
        --chimSegmentMin 10 \\
        --chimOutType WithinBAM HardClip \\
        --chimJunctionOverhangMin 10 \\
        --chimScoreDropMax 30 \\
        --chimScoreJunctionNonGTAG 0 \\
        --chimScoreSeparation 1 \\
        --chimSegmentReadGapMax 3 \\
        --chimMultimapNmax 50 \\
        --outFileNamePrefix ${meta.id}_${meta.timepoint}.arriba.
    """

    stub:
    """
    touch ${meta.id}_${meta.timepoint}.arriba.Aligned.out.bam
    """
}
