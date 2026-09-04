// Ported from: NeoadjLOH/specHLA/run_sample.pbs, Step 1 ("extract HLA reads")
// https://github.com/JaydenBeckwith/NeoadjLOH
//
// Pulls reads mapping to the HLA loci (and ALT contigs) out of the TUMOR
// BAM — SpecHLA's LOH workflow is documented as tumor-only ("Detect HLA
// LOH in tumor samples"), so unlike HLA_TYPING_XHLA (which deliberately
// runs on the normal BAM to dodge tumor HLA-LOH bias), this branch runs on
// tumor by design: the LOH call *is* the tumor-vs-germline-baseline signal
// SpecHLA computes internally from purity/ploidy, not something we can get
// from a normal BAM.
//
// The fastq-pairing fallback below is copied verbatim from run_sample.pbs,
// gotcha comment included: this cohort's sample IDs are digit-heavy (e.g.
// "...0000141960..."), so a loose "*1*"/"*2*" glob spuriously matches
// digits inside the sample ID itself rather than the read-pair marker —
// this is what silently handed SpecHLA the wrong fastq and made every
// locus come back no_match on Jayden's first run. Anchor on the literal
// "_1."/"_2." suffix and exclude the unpaired-singleton file instead.

process SPECHLA_EXTRACT_HLA_READS {
    tag "${meta.id}"
    label 'process_medium'
    container params.containers.spechla
    publishDir "${params.outdir}/${meta.id}/hla_loh/extracted", mode: 'copy'

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai)

    output:
    tuple val(meta), path("${meta.id}_extract_1.fq.gz"), path("${meta.id}_extract_2.fq.gz"), emit: reads

    script:
    """
    mkdir -p extracted_out
    ${params.spechla_extract_cmd} -s ${meta.id} -b ${tumor_bam} -r ${params.spechla_ref_build} -o extracted_out

    FQ1="extracted_out/${meta.id}_extract_1.fq.gz"
    FQ2="extracted_out/${meta.id}_extract_2.fq.gz"
    if [[ ! -s "\$FQ1" || ! -s "\$FQ2" ]]; then
        # Fall back to a search, but anchor on the real "_1."/"_2." read-pair
        # suffix and explicitly exclude the unpaired-singleton file — see the
        # header comment above for why a bare "*1*"/"*2*" glob breaks on this
        # cohort's digit-heavy sample IDs.
        FQ1=\$(find extracted_out -iname "*_1.f*q.gz" ! -iname "*unpaired*" | sort | head -n1)
        FQ2=\$(find extracted_out -iname "*_2.f*q.gz" ! -iname "*unpaired*" | sort | head -n1)
    fi
    if [[ -z "\${FQ1:-}" || -z "\${FQ2:-}" || ! -s "\$FQ1" || ! -s "\$FQ2" ]]; then
        echo "ERROR: could not find paired extracted fastqs in extracted_out for ${meta.id}" >&2
        exit 1
    fi
    echo "[INFO] Extracted reads: \$FQ1 (\$(du -h "\$FQ1" | cut -f1)) / \$FQ2 (\$(du -h "\$FQ2" | cut -f1))"

    # standardize to the fixed filenames declared as this process's output,
    # whichever of the two paths above actually resolved
    cp "\$FQ1" ${meta.id}_extract_1.fq.gz
    cp "\$FQ2" ${meta.id}_extract_2.fq.gz
    """
}
