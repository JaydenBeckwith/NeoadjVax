// Ported from: RNA_variant_pipeline.sh Step 3/12 (STAR two-pass alignment)
// outFilterMultimapNmax 2 and outSAMattributes NH HI AS NM MD carried over
// verbatim from the original PBS job.

process STAR_ALIGN {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_high'
    container params.containers.star
    publishDir "${params.outdir}/${meta.id}/rna_star/${meta.timepoint}", mode: 'copy', pattern: '*.Log.final.out'

    input:
    tuple val(meta), path(r1), path(r2)
    path star_index

    output:
    tuple val(meta), path("${meta.id}_${meta.timepoint}.Aligned.sortedByCoord.out.bam"), emit: bam
    path "${meta.id}_${meta.timepoint}.Log.final.out", emit: log

    script:
    """
    STAR \\
        --runThreadN ${task.cpus} \\
        --genomeDir ${star_index} \\
        --readFilesIn ${r1} ${r2} \\
        --readFilesCommand zcat \\
        --twopassMode Basic \\
        --outSAMtype BAM SortedByCoordinate \\
        --limitBAMsortRAM 40000000000 \\
        --outFileNamePrefix ${meta.id}_${meta.timepoint}. \\
        --outSAMattributes NH HI AS NM MD \\
        --outFilterMultimapNmax 2

    samtools index ${meta.id}_${meta.timepoint}.Aligned.sortedByCoord.out.bam
    """
}
