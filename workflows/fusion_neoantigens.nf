/*
 * Gene fusion neoantigen discovery branch. New territory — no existing
 * fusion-calling infrastructure to reuse (unlike splicing).
 *
 *   RNA FASTQ ─▶ STAR_FUSION ─▶ AGFUSION_ANNOTATE (stub) ─▶ PVACFUSE_RUN
 *
 * HLA alleles reuse whatever PVACSEQ_CORE typed/was given for that
 * patient — this workflow takes an hla_by_patient channel rather than
 * retyping, to avoid running OptiType twice per patient.
 */

include { STAR_FUSION }       from '../modules/local/fusion/star_fusion'
include { AGFUSION_ANNOTATE } from '../modules/local/fusion/agfusion_annotate'
include { PVACFUSE_RUN }      from '../modules/local/pvactools/pvacfuse'

workflow FUSION_NEOANTIGENS {

    take:
    rna_samplesheet_ch   // [patient_id, timepoint, rna_fastq_r1, rna_fastq_r2, ...]

    main:
    if (!params.ctat_resource_lib) {
        error "FUSION_NEOANTIGENS needs --ctat_resource_lib (STAR-Fusion CTAT genome lib)"
    }
    ctat_lib = Channel.fromPath(params.ctat_resource_lib).collect()

    fastq_ch = rna_samplesheet_ch
        .filter { row -> row.rna_fastq_r1 && row.rna_fastq_r2 }
        .map { row ->
            def meta = [id: row.patient_id, timepoint: row.timepoint]
            tuple(meta, file(row.rna_fastq_r1), file(row.rna_fastq_r2))
        }
    // NOTE: rows that only supply rna_bam (the common case per
    // assets/samplesheet_schema.md) aren't picked up here yet — either
    // reuse RNA_VARIANT_CALLING's BAM_TO_FASTQ output (requires wiring this
    // workflow to run after/alongside it rather than standalone) or add a
    // bam->fastq step here too. Left as-is until we decide whether fusion
    // calling should share the RNA branch's FASTQs or run fully independently.

    STAR_FUSION(fastq_ch, ctat_lib)
    AGFUSION_ANNOTATE(STAR_FUSION.out.fusions)

    if (!params.hla_alleles_manual) {
        log.warn "FUSION_NEOANTIGENS: no --hla_alleles_manual given and this workflow doesn't yet consume PVACSEQ_CORE's typed HLA output — pvacfuse will need HLA alleles supplied another way for now"
    }
    hla_ch = AGFUSION_ANNOTATE.out.annotated.map { meta, dir -> tuple(meta, dir, params.hla_alleles_manual) }

    PVACFUSE_RUN(hla_ch)

    emit:
    results = PVACFUSE_RUN.out.results
}
