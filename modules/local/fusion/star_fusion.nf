// STAR-Fusion on the patient's baseline tumour RNA FASTQ pair.
//
// STAR-Fusion runs its own STAR alignment against the CTAT genome lib
// (params.ctat_genome_lib, the ctat_genome_lib_build_dir of a plug-n-play
// bundle), not against this pipeline's genome_fasta: the lib carries its
// own GRCh38 genome, GENCODE annotation and STAR index. That index must
// have been built by a STAR release the image's STAR can read (2.7.4a or
// later for STAR 2.7.11b), so use a lib released for STAR-Fusion 1.13+.
//
// Both prediction tables are kept: the abridged one is what people read;
// the full one is what AGFusion and pVACfuse are given, as in the Gadi
// PBS scripts (pVACfuse --starfusion-file takes the full table for read
// support and FFPM).

process STAR_FUSION {
    tag "${meta.id}"
    label 'process_starfusion'
    container params.containers.starfusion
    publishDir "${params.outdir}/${meta.id}/fusion/star_fusion", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)
    path ctat_genome_lib

    output:
    tuple val(meta), path("star-fusion.fusion_predictions.tsv"), emit: predictions
    tuple val(meta), path("star-fusion.fusion_predictions.abridged.tsv"), emit: abridged

    script:
    """
    # STAR-Fusion's temp files default to /tmp, which is small on Gadi nodes
    mkdir -p sf_tmp
    STAR-Fusion \\
        --genome_lib_dir ${ctat_genome_lib} \\
        --left_fq ${r1} \\
        --right_fq ${r2} \\
        --CPU ${task.cpus} \\
        --output_dir star_fusion \\
        --tmpdir "\$PWD/sf_tmp" \\
        ${params.starfusion_extra_args}
    cp star_fusion/star-fusion.fusion_predictions.tsv star_fusion/star-fusion.fusion_predictions.abridged.tsv .
    rm -rf sf_tmp
    """

    stub:
    """
    header='#FusionName\\tJunctionReadCount\\tSpanningFragCount\\test_J\\test_S\\tSpliceType\\tLeftGene\\tLeftBreakpoint\\tRightGene\\tRightBreakpoint\\tLargeAnchorSupport\\tFFPM\\tLeftBreakDinuc\\tLeftBreakEntropy\\tRightBreakDinuc\\tRightBreakEntropy\\tannots'
    printf "\$header\\tJunctionReads\\tSpanningFrags\\n" > star-fusion.fusion_predictions.tsv
    printf "\$header\\n" > star-fusion.fusion_predictions.abridged.tsv
    """
}
