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

// ---------------------------------------------------------------------
// Samplesheet validation — runs synchronously at launch (splitCsv() on a
// plain `file()` reads it eagerly, not as a reactive Channel) so a missing
// or inconsistent row fails immediately with a clear message instead of
// surfacing as a confusing empty-channel gap deep into execution. Worth
// having now that each row can supply FASTQ, an already-aligned BAM, or an
// already-called VCF, in any combination.
// ---------------------------------------------------------------------
def validateDnaSamplesheet(path) {
    file(path).splitCsv(header: true).each { row ->
        def hasVcf  = row.dna_vcf as boolean
        def hasHla  = row.hla_alleles as boolean
        ['tumor', 'normal'].each { sample_type ->
            def hasBam = row["${sample_type}_bam"] as boolean
            def hasFq  = (row["${sample_type}_r1"] as boolean) && (row["${sample_type}_r2"] as boolean)
            if (!hasVcf && !hasBam && !hasFq) {
                error "dna_samplesheet: patient ${row.patient_id} has no ${sample_type}_bam and no complete ${sample_type}_r1/${sample_type}_r2 pair, and no dna_vcf to fall back on"
            }
        }
        def hasNormalSource = (row.normal_bam as boolean) || ((row.normal_r1 as boolean) && (row.normal_r2 as boolean))
        if (hasVcf && params.run_pvacseq_core && !hasHla && !hasNormalSource) {
            error "dna_samplesheet: patient ${row.patient_id} supplies dna_vcf with no normal sample and no hla_alleles — PVACSEQ_CORE has no normal BAM to type HLA from and no manual override. Add hla_alleles, or a normal_bam/normal_r1+r2 pair."
        }
    }
}

def validateRnaSamplesheet(path) {
    file(path).splitCsv(header: true).each { row ->
        def hasVcf = row.rna_vcf as boolean
        def hasBam = row.rna_bam as boolean
        def hasFq  = (row.rna_fastq_r1 as boolean) && (row.rna_fastq_r2 as boolean)
        if (!hasVcf && !hasBam && !hasFq) {
            error "rna_samplesheet: patient ${row.patient_id} timepoint ${row.timepoint} has none of rna_vcf, rna_bam, or a complete rna_fastq_r1/rna_fastq_r2 pair"
        }
    }
}

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
        validateDnaSamplesheet(params.dna_samplesheet)
        dna_samplesheet_ch = Channel.fromPath(params.dna_samplesheet).splitCsv(header: true)
    }

    if (params.run_rna_variant_calling || params.run_pvacseq_core) {
        if (!params.rna_samplesheet) error "Please provide --rna_samplesheet (see assets/samplesheet_schema.md)"
        validateRnaSamplesheet(params.rna_samplesheet)
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
