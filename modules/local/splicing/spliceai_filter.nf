// Reuses the OUTPUT of Jayden's existing SpliceAI variant-annotation
// pipeline on Gadi ([[spliceai-variant-pipeline]] — per-chromosome PBS,
// Singularity) rather than re-running SpliceAI here. This module just
// ingests that pipeline's annotated VCFs (per patient-timepoint) and keeps
// variants above a splice-altering delta-score threshold — the candidate
// set that then needs translating into novel-junction peptides.
//
// params.spliceai_annotated_vcf_dir must contain one VCF per
// patient-timepoint named <patient_id>_<timepoint>*.vcf(.gz) — adjust the
// glob below once we confirm the existing pipeline's actual naming.

process SPLICEAI_FILTER {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_low'
    container params.containers.gatk  // has bcftools/python; swap if a lighter container is preferred
    publishDir "${params.outdir}/${meta.id}/splicing/${meta.timepoint}", mode: 'copy'

    input:
    tuple val(meta), path(spliceai_vcf)

    output:
    tuple val(meta), path("${meta.id}_${meta.timepoint}.spliceai_filtered.vcf"), emit: vcf

    script:
    def min_delta = params.spliceai_min_delta_score ?: 0.5
    """
    # Keep records where any DS_* INFO field (SpliceAI delta score) exceeds
    # the threshold. bcftools filter expression carried over from
    # [[neoadjuvant-splicing]]'s existing delta-score cutoffs — confirm the
    # exact SpliceAI INFO field names match what the existing pipeline emits.
    bcftools view -i 'INFO/DS_AG>${min_delta} || INFO/DS_AL>${min_delta} || INFO/DS_DG>${min_delta} || INFO/DS_DL>${min_delta}' \\
        ${spliceai_vcf} -O v -o ${meta.id}_${meta.timepoint}.spliceai_filtered.vcf
    """
}
