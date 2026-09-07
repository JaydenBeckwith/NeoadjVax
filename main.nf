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
 *                                                    ├─▶ splicing branch  ─▶ splice2neo peptide contexts (--run_splicing_neoantigens)
 *                                                    ├─▶ fusion branch    ─▶ pVACfuse (--run_fusion_neoantigens)
 *                                                    └─▶ ERV branch       ─▶ TBD       (--run_erv_neoantigens)
 *
 * DNA and RNA are joined on patient_id, not row-for-row — see
 * assets/samplesheet_schema.md for why there are two separate samplesheets.
 */

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------
// Pipeline selection (--pipelines) — a friendlier alternative to setting
// each --run_* flag individually. Maps a comma-separated list of names
// (case/spacing-insensitive: "geneFusion", "gene-fusion", "Gene Fusion"
// and "gene_fusion" all resolve to the same entry) onto the underlying
// run_* flags. Some names turn on more than one flag together —
// "somatic_neoantigen" needs DNA calling + RNA calling + pVACseq core all
// three, since the core neoantigen score is DNA calls filtered for RNA
// support (see PVACSEQ_CORE below), not any one of those alone.
//
// When --pipelines is set, it's authoritative for every run_* flag —
// unlisted branches are turned OFF even if their individual --run_*
// default is true. Leave --pipelines unset to keep the old behaviour
// (each --run_* flag controls its own branch, independently).
// ---------------------------------------------------------------------
def applyPipelineSelection(String pipelinesParam) {
    // key: normalized name (lowercase, punctuation/spaces stripped) ->
    // the run_* flags it turns on. Some names turn on more than one flag
    // together — "somatic_neoantigen" needs DNA calling + RNA calling +
    // pVACseq core all three, since the core neoantigen score is DNA
    // calls filtered for RNA support (see PVACSEQ_CORE below), not any
    // one of those alone.
    def registry = [
        'dnavariantcalling'  : ['run_dna_variant_calling'],
        'rnavariantcalling'  : ['run_rna_variant_calling'],
        'somaticneoantigen'  : ['run_dna_variant_calling', 'run_rna_variant_calling', 'run_pvacseq_core'],
        'coreneoantigen'     : ['run_dna_variant_calling', 'run_rna_variant_calling', 'run_pvacseq_core'],
        'pvacseqcore'        : ['run_dna_variant_calling', 'run_rna_variant_calling', 'run_pvacseq_core'],
        'splicingneoantigen' : ['run_splicing_neoantigens'],
        'splicing'           : ['run_splicing_neoantigens'],
        'genefusion'         : ['run_fusion_neoantigens'],
        'fusion'             : ['run_fusion_neoantigens'],
        'fusionneoantigen'   : ['run_fusion_neoantigens'],
        'ervneoantigen'      : ['run_erv_neoantigens'],
        'erv'                : ['run_erv_neoantigens'],
        'ervdna'             : ['run_erv_dna'],
        'ervcaller'          : ['run_erv_dna'],
        'purityploidy'       : ['run_purity_ploidy'],
        'sequenza'           : ['run_purity_ploidy'],
        'hlaloh'             : ['run_hla_loh'],
        'loh'                : ['run_hla_loh'],
    ]
    def allFlags = [
        'run_dna_variant_calling', 'run_rna_variant_calling', 'run_pvacseq_core',
        'run_splicing_neoantigens', 'run_fusion_neoantigens', 'run_erv_neoantigens',
        'run_purity_ploidy', 'run_hla_loh', 'run_erv_dna',
    ]
    def normalize = { String s -> s.toLowerCase().replaceAll(/[^a-z0-9]/, '') }

    def requested = pipelinesParam.split(',').collect { it.trim() }.findAll { it }
    def flagsToEnable = [] as Set
    def unknown = []
    requested.each { name ->
        def key = normalize(name)
        if (registry.containsKey(key)) {
            flagsToEnable.addAll(registry[key])
        } else {
            unknown << name
        }
    }
    if (unknown) {
        def validNames = registry.keySet().sort().join(', ')
        error "--pipelines: unrecognized name(s) ${unknown} — valid names (case/spacing/punctuation-insensitive): ${validNames}. Comma-separate to run more than one, e.g. --pipelines somatic_neoantigen,gene_fusion"
    }
    allFlags.each { flag -> params[flag] = flagsToEnable.contains(flag) }
    log.info "--pipelines '${pipelinesParam}' resolved to: ${allFlags.findAll { params[it] }.sort()}"
}

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

// Only needed when --run_hla_loh is used WITHOUT --run_dna_variant_calling
// (standalone mode) — see docs/ARCHITECTURE.md / assets/samplesheet_schema.md.
def validateHlaLohSamplesheet(path) {
    file(path).splitCsv(header: true).each { row ->
        if (!(row.patient_id as boolean)) {
            error "hla_loh_samplesheet: row missing patient_id"
        }
        if (!(row.tumor_bam as boolean)) {
            error "hla_loh_samplesheet: patient ${row.patient_id} has no tumor_bam"
        }
        if (!(row.purity as boolean) || !(row.ploidy as boolean)) {
            error "hla_loh_samplesheet: patient ${row.patient_id} needs both purity and ploidy (SpecHLA's LOH step requires both)"
        }
    }
}

include { DNA_VARIANT_CALLING }  from './workflows/dna_variant_calling'
include { RNA_VARIANT_CALLING }  from './workflows/rna_variant_calling'
include { PVACSEQ_CORE }         from './workflows/pvacseq_core'
include { SPLICING_NEOANTIGENS; validateSplicingSamplesheets } from './workflows/splicing_neoantigens'
include { FUSION_NEOANTIGENS }   from './workflows/fusion_neoantigens'
include { ERV_NEOANTIGENS }      from './workflows/erv_neoantigens'
include { ERV_DNA; validateErvDnaSamplesheet } from './workflows/erv_dna'
include { PURITY_PLOIDY }        from './workflows/purity_ploidy'
include { HLA_LOH }              from './workflows/hla_loh'

workflow {

    if (params.pipelines) {
        applyPipelineSelection(params.pipelines)
    }

    dna_samplesheet_ch = Channel.empty()
    rna_samplesheet_ch = Channel.empty()

    if (params.run_dna_variant_calling || params.run_pvacseq_core || params.run_erv_dna || params.run_splicing_neoantigens) {
        if (!params.dna_samplesheet) error "Please provide --dna_samplesheet (see assets/samplesheet_schema.md)"
        if (params.run_dna_variant_calling || params.run_pvacseq_core) {
            validateDnaSamplesheet(params.dna_samplesheet)
        }
        if (params.run_erv_dna) validateErvDnaSamplesheet(params.dna_samplesheet, params.run_dna_variant_calling)
        dna_samplesheet_ch = Channel.fromPath(params.dna_samplesheet).splitCsv(header: true)
    }

    // Every RNA-derived branch consumes the longitudinal samplesheet, even
    // when it runs independently of RNA variant calling.  In particular,
    // `--pipelines gene_fusion` and `--pipelines erv` must not silently
    // receive an empty input channel.
    if (params.run_rna_variant_calling || params.run_pvacseq_core ||
        params.run_splicing_neoantigens || params.run_fusion_neoantigens ||
        params.run_erv_neoantigens) {
        if (!params.rna_samplesheet) error "Please provide --rna_samplesheet (see assets/samplesheet_schema.md)"
        if (!params.run_splicing_neoantigens || params.run_rna_variant_calling || params.run_pvacseq_core ||
            params.run_fusion_neoantigens || params.run_erv_neoantigens) {
            validateRnaSamplesheet(params.rna_samplesheet)
        }
        rna_samplesheet_ch = Channel.fromPath(params.rna_samplesheet).splitCsv(header: true)
    }

    if (params.run_splicing_neoantigens) {
        validateSplicingSamplesheets(params.dna_samplesheet, params.rna_samplesheet, params.run_dna_variant_calling, params)
    }

    dna_variants_ch   = Channel.empty()
    dna_normal_bam_ch = Channel.empty()
    dna_tumor_bam_ch  = Channel.empty()
    rna_variants_ch   = Channel.empty()

    if (params.run_dna_variant_calling) {
        DNA_VARIANT_CALLING(dna_samplesheet_ch)
        dna_variants_ch   = DNA_VARIANT_CALLING.out.filtered_vcf
        dna_normal_bam_ch = DNA_VARIANT_CALLING.out.normal_bam
        dna_tumor_bam_ch  = DNA_VARIANT_CALLING.out.tumor_bam
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
        SPLICING_NEOANTIGENS(dna_samplesheet_ch, rna_samplesheet_ch, dna_variants_ch, params.run_dna_variant_calling)
    }

    if (params.run_fusion_neoantigens) {
        FUSION_NEOANTIGENS(rna_samplesheet_ch)
    }

    if (params.run_erv_neoantigens) {
        ERV_NEOANTIGENS(rna_samplesheet_ch)
    }

    if (params.run_erv_dna) {
        ERV_DNA(dna_samplesheet_ch, dna_tumor_bam_ch, dna_normal_bam_ch, params.run_dna_variant_calling)
    }

    // ---------------------------------------------------------------
    // Tumor purity/ploidy (Sequenza) — optional, --run_purity_ploidy.
    // Fully ported from https://github.com/JaydenBeckwith/Sequenza_tools —
    // WGS/200bp bins only, see workflows/purity_ploidy.nf.
    // ---------------------------------------------------------------
    purity_ploidy_ch = Channel.empty()
    if (params.run_purity_ploidy) {
        if (!params.run_dna_variant_calling) {
            error "--run_purity_ploidy needs --run_dna_variant_calling too (it runs on DNA_VARIANT_CALLING's tumor+normal BAM pairs by patient) — a from-samplesheet standalone mode isn't built yet, ask if you need that."
        }
        PURITY_PLOIDY(dna_tumor_bam_ch, dna_normal_bam_ch)
        purity_ploidy_ch = PURITY_PLOIDY.out.purity_ploidy
    }

    // ---------------------------------------------------------------
    // HLA loss-of-heterozygosity (SpecHLA) — optional, --run_hla_loh.
    // Independently runnable two ways:
    //   1. --hla_loh_samplesheet (patient_id,tumor_bam,purity,ploidy) —
    //      no DNA/RNA branches needed at all.
    //   2. Chained off --run_dna_variant_calling's tumor BAM, with
    //      purity/ploidy from either --purity_ploidy_csv (manual,
    //      NeoadjLOH PURITY_CSV format: sample,purity,ploidy) or
    //      --run_purity_ploidy (once that branch is real).
    // See assets/samplesheet_schema.md.
    // ---------------------------------------------------------------
    if (params.run_hla_loh) {
        def hla_loh_tumor_bam_ch
        def hla_loh_purity_ploidy_ch

        if (params.hla_loh_samplesheet) {
            validateHlaLohSamplesheet(params.hla_loh_samplesheet)
            def hla_loh_rows_ch = Channel.fromPath(params.hla_loh_samplesheet).splitCsv(header: true)
            hla_loh_tumor_bam_ch = hla_loh_rows_ch.map { row ->
                def meta = [id: row.patient_id]
                tuple(meta, file(row.tumor_bam), file("${row.tumor_bam}.bai"))
            }
            hla_loh_purity_ploidy_ch = hla_loh_rows_ch.map { row ->
                tuple([id: row.patient_id], row.purity, row.ploidy)
            }
        } else {
            if (!params.run_dna_variant_calling) {
                error "--run_hla_loh needs either --hla_loh_samplesheet (standalone) or --run_dna_variant_calling (to supply the tumor BAM)"
            }
            hla_loh_tumor_bam_ch = dna_tumor_bam_ch
            if (params.purity_ploidy_csv) {
                hla_loh_purity_ploidy_ch = Channel.fromPath(params.purity_ploidy_csv)
                    .splitCsv(header: true)
                    .map { row -> tuple([id: row.sample], row.purity, row.ploidy) }
            } else if (params.run_purity_ploidy) {
                hla_loh_purity_ploidy_ch = purity_ploidy_ch
            } else {
                error "--run_hla_loh needs purity/ploidy from somewhere: pass --purity_ploidy_csv, or also set --run_purity_ploidy true"
            }
        }

        HLA_LOH(hla_loh_tumor_bam_ch, hla_loh_purity_ploidy_ch)
    }
}

workflow.onComplete {
    log.info "NeoadjVax finished — status: ${workflow.success ? 'OK' : 'FAILED'} — results: ${params.outdir}"
}
