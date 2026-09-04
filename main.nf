#!/usr/bin/env nextflow
/*
 * ===========================================================================
 * NeoadjVax — neoadjuvant melanoma vaccine neoantigen discovery pipeline
 * ===========================================================================
 * See docs/ARCHITECTURE.md for the full stage-by-stage design and open
 * decisions on the splicing/fusion/ERV branches.
 *
 *   DNA WGS/WES (tumor+normal, baseline) ─┐
 *                                          ├─▶ variant matching ─▶ pVACseq ─▶ core neoantigens (per RNA timepoint)
 *   RNA-seq (per timepoint) ───────────────┘        │
 *                                                    ├─▶ splicing branch  ─▶ pVACseq   (--run_splicing_neoantigens)
 *                                                    ├─▶ fusion branch    ─▶ pVACfuse (--run_fusion_neoantigens)
 *                                                    └─▶ ERV branch       ─▶ TBD       (--run_erv_neoantigens)
 *
 * DNA and RNA are joined on patient_id, not row-for-row — see
 * assets/samplesheet_schema.md for why there are two separate samplesheets.
 */

nextflow.enable.dsl = 2

include { DNA_VARIANT_CALLING }  from './workflows/dna_variant_calling'
include { RNA_VARIANT_CALLING }  from './workflows/rna_variant_calling'
include { PVACSEQ_CORE }         from './workflows/pvacseq_core'
include { SPLICING_NEOANTIGENS } from './workflows/splicing_neoantigens'
include { FUSION_NEOANTIGENS }   from './workflows/fusion_neoantigens'
include { ERV_NEOANTIGENS }      from './workflows/erv_neoantigens'

workflow {

    dna_samplesheet_ch = Channel.empty()
    rna_samplesheet_ch = Channel.empty()

    if (params.run_dna_variant_calling || params.run_pvacseq_core) {
        if (!params.dna_samplesheet) error "Please provide --dna_samplesheet (see assets/samplesheet_schema.md)"
        dna_samplesheet_ch = Channel.fromPath(params.dna_samplesheet).splitCsv(header: true)
    }

    if (params.run_rna_variant_calling || params.run_pvacseq_core) {
        if (!params.rna_samplesheet) error "Please provide --rna_samplesheet (see assets/samplesheet_schema.md)"
        rna_samplesheet_ch = Channel.fromPath(params.rna_samplesheet).splitCsv(header: true)
    }

    dna_variants_ch   = Channel.empty()
    dna_normal_bam_ch = Channel.empty()
    rna_variants_ch   = Channel.empty()

    if (params.run_dna_variant_calling) {
        DNA_VARIANT_CALLING(dna_samplesheet_ch)
        dna_variants_ch   = DNA_VARIANT_CALLING.out.filtered_vcf
        dna_normal_bam_ch = DNA_VARIANT_CALLING.out.normal_bam
    }

    if (params.run_rna_variant_calling) {
        RNA_VARIANT_CALLING(rna_samplesheet_ch)
        rna_variants_ch = RNA_VARIANT_CALLING.out.filtered_vcf
    }

    if (params.run_pvacseq_core) {
        if (!params.run_dna_variant_calling || !params.run_rna_variant_calling) {
            error "--run_pvacseq_core needs both --run_dna_variant_calling and --run_rna_variant_calling (the core neoantigen score is baseline DNA calls filtered for per-timepoint RNA support)"
        }
        PVACSEQ_CORE(dna_samplesheet_ch, dna_variants_ch, dna_normal_bam_ch, rna_variants_ch)
    }

    if (params.run_splicing_neoantigens) {
        SPLICING_NEOANTIGENS(rna_samplesheet_ch)
    }

    if (params.run_fusion_neoantigens) {
        FUSION_NEOANTIGENS(rna_samplesheet_ch)
    }

    if (params.run_erv_neoantigens) {
        ERV_NEOANTIGENS(rna_samplesheet_ch)
    }
}

workflow.onComplete {
    log.info "NeoadjVax finished — status: ${workflow.success ? 'OK' : 'FAILED'} — results: ${params.outdir}"
}
