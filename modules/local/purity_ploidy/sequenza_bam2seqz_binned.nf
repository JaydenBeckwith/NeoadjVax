// Ported from: Sequenza_tools/sequenza_step1.sh, Steps 1+2 (bam2seqz + seqz_binning)
// https://github.com/JaydenBeckwith/Sequenza_tools
//
// One CLI pair per chromosome, exactly as the original ran them: just
// re-homed onto Nextflow's own parallelism instead of the original's
// manual bash `&`/`wait` background-job loop. That's a deliberate
// deviation worth flagging: the original packed all 24 chromosomes for
// one sample into a single 24-cpu/64GB PBS allocation (`sequenza_step1.sh`'s
// own `#PBS -l ncpus=24 -l mem=64GB`); here, each (sample, chromosome) pair
// becomes its OWN Nextflow task/PBS job, so per-task resources are sized
// per-chromosome, not per-whole-sample: see the process_low label below
// (worth tuning after timing one real chromosome on Gadi, same as
// everywhere else this scaffold says "confirm before scaling up").
//
// -w 200 (WGS bin size) is hardcoded in the original script; wired here to
// params.sequenza_bin_wgs (default 200, matching) / params.sequenza_bin_wes
// (default 50) so a WES cohort can override it: the original only ever
// ran on WGS, so the WES path is new plumbing riding on your stated
// convention, not something the original script itself does.

process SEQUENZA_BAM2SEQZ_BINNED {
    tag "${meta.id}:${chrom}"
    label 'process_low'
    container params.containers.sequenza

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), val(chrom)
    path fasta
    path fasta_fai
    path gc_file
    val bin_size

    output:
    tuple val(meta), val(chrom), path("${meta.id}.${chrom}_bin${bin_size}.seqz.gz"), emit: binned_seqz

    script:
    """
    sequenza-utils bam2seqz \\
        -gc ${gc_file} \\
        --fasta ${fasta} \\
        -n ${normal_bam} \\
        -t ${tumor_bam} \\
        -C ${chrom} \\
        | gzip > ${meta.id}.${chrom}.seqz.gz

    sequenza-utils seqz_binning \\
        -w ${bin_size} \\
        -s ${meta.id}.${chrom}.seqz.gz \\
        | gzip > ${meta.id}.${chrom}_bin${bin_size}.seqz.gz
    """
}
