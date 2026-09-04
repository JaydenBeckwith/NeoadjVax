// HLA typing (class I + class II) via xHLA, replacing the earlier
// OptiType (class I only) + arcasHLA (class II, TODO stub) combination —
// xHLA types both from a single indexed BAM in one pass, no FASTQ
// conversion needed.
//
// Runs on the NORMAL (germline) DNA BAM, not tumor — same reasoning as
// before the swap: avoids bias from tumor HLA-LOH, which
// [[loh-analysis]] already tracks as its own workstream.
//
// BAM requirements (confirmed via the xHLA GitHub wiki, "BAMs compatible
// with xHLA"): an INDEXED BAM aligned with BWA-MEM against hg38 WITHOUT
// ALT contigs — GATK-bundle-style references with ALT contigs can trap
// reads that should land in chr6:29844528-33100696, breaking typing.
// params.genome_fasta here is GENCODE's "primary assembly" FASTA, which
// excludes ALT/patch/haplotype scaffolds by construction, so this should
// already be compatible — still worth a sanity check against a real BAM
// before trusting it in production.
//
// Container: humanlongevity/hla is a Docker Hub image (the upstream repo
// predates Biocontainers-style distribution); Singularity pulls it
// directly via the docker:// reference, same pattern as pvactools below.
// Pinned to the 0.0.0 tag — confirmed via a third-party Nextflow wrapper's
// Dockerfile (alexvpickering/nf-xhla), since the upstream repo/Docker Hub
// page don't document image tags directly.
//
// Entrypoint: Nextflow's `container` directive runs OUR script inside the
// container rather than the image's baked-in ENTRYPOINT, so the command
// below invokes xHLA's run script directly — `python /opt/bin/run.py`,
// also confirmed via nf-xhla (neither xHLA's repo nor its Docker Hub page
// documents the in-container script path). Worth a one-off check
// (`singularity exec <image> python /opt/bin/run.py --help`) before
// trusting this in production, since it's the one detail sourced from a
// third party rather than xHLA's own docs.
//
// IMPORTANT — Gadi compute nodes have no external network access, so this
// container must be pre-pulled into the shared Singularity cache from the
// login node (which does have network) before this ever runs as a PBS
// job — see docs/RUNNING_ON_GADI.md.

process HLA_TYPING_XHLA {
    tag "${meta.id}"
    label 'process_medium'
    container params.containers.xhla
    publishDir "${params.outdir}/${meta.id}/hla_typing", mode: 'copy'

    input:
    tuple val(meta), path(normal_bam), path(normal_bai)

    output:
    tuple val(meta), env(HLA_ALLELES), emit: hla_alleles
    path "xhla_out/report-${meta.id}-hla.json", emit: report
    path "excluded_dq_dp_beta_only.txt", emit: excluded_alleles, optional: true

    script:
    """
    mkdir -p xhla_out
    python /opt/bin/run.py \\
        --sample_id ${meta.id} \\
        --input_bam_path ${normal_bam} \\
        --output_path xhla_out

    parse_xhla_result.py xhla_out/report-${meta.id}-hla.json \\
        -o pvacseq_hla_alleles.txt \\
        --excluded-out excluded_dq_dp_beta_only.txt

    HLA_ALLELES=\$(cat pvacseq_hla_alleles.txt)
    echo "[INFO] \${HLA_ALLELES}"
    """
}
