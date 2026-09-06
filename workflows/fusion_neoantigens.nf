/*
 * Gene fusion neoantigen discovery branch. Runs the requested fusion
 * caller(s) in parallel per patient-timepoint — default both:
 *
 *   RNA FASTQ (direct, or BAM→FASTQ) ─┬─▶ STAR_FUSION      ─┐
 *                                     └─▶ ARRIBA_STAR_ALIGN ┴─▶ ARRIBA ─┐
 *                                                                       ├─▶ AGFUSION_ANNOTATE (stub, per caller) ─▶ PVACFUSE_RUN (per caller)
 *                                              (STAR_FUSION fusions) ───┘
 *
 * --fusion_callers (default 'starfusion,arriba') picks which caller(s)
 * run; results from each are kept separate all the way through (tagged by
 * caller, published to separate subdirectories) rather than merged —
 * AGFusion/pVACfuse take one caller's format at a time, and there's no
 * settled reconciliation logic for two callers disagreeing on the same
 * fusion event (that's a real open decision if/when you want a combined
 * call set, not something to silently invent here).
 *
 * STAR_FUSION and ARRIBA_STAR_ALIGN each run their OWN STAR pass — they
 * are not sharing one alignment. STAR-Fusion's wrapper needs the CTAT
 * genome lib's bundled index; Arriba needs relaxed multimapping + specific
 * chimeric flags this pipeline's RNA branch doesn't set. See each
 * module's own comment.
 *
 * Takes rna_bam OR rna_fastq_r1/r2 from the samplesheet (previously only
 * rna_fastq_r1/r2 was consumed here — an open decision from the original
 * scaffold, now resolved by reusing RNA_VARIANT_CALLING's own BAM_TO_FASTQ
 * module for rna_bam rows).
 *
 * HLA alleles reuse whatever PVACSEQ_CORE typed/was given for that
 * patient — this workflow takes an hla_by_patient channel rather than
 * retyping, to avoid running xHLA twice per patient.
 */

include { BAM_TO_FASTQ }      from '../modules/local/rna/bam_to_fastq'
include { STAR_FUSION }       from '../modules/local/fusion/star_fusion'
include { ARRIBA_STAR_ALIGN } from '../modules/local/fusion/arriba_align'
include { ARRIBA }            from '../modules/local/fusion/arriba'
include { AGFUSION_ANNOTATE } from '../modules/local/fusion/agfusion_annotate'
include { PVACFUSE_RUN }      from '../modules/local/pvactools/pvacfuse'

workflow FUSION_NEOANTIGENS {

    take:
    rna_samplesheet_ch   // [patient_id, timepoint, rna_bam, rna_fastq_r1, rna_fastq_r2, ...]

    main:
    def requested_callers = params.fusion_callers.split(',').collect { it.trim().toLowerCase() }
    def valid_callers = ['starfusion', 'arriba']
    def unknown = requested_callers - valid_callers
    if (unknown) {
        error "FUSION_NEOANTIGENS: unknown --fusion_callers entr${unknown.size() > 1 ? 'ies' : 'y'} ${unknown} — valid values are ${valid_callers} (comma-separated to run more than one)"
    }
    def run_starfusion = 'starfusion' in requested_callers
    def run_arriba = 'arriba' in requested_callers

    if (run_starfusion && !params.ctat_resource_lib) {
        error "FUSION_NEOANTIGENS: 'starfusion' requested in --fusion_callers but --ctat_resource_lib is not set (STAR-Fusion CTAT genome lib)"
    }

    bam_rows_ch = rna_samplesheet_ch
        .filter { row -> row.rna_bam && !(row.rna_fastq_r1 && row.rna_fastq_r2) }
        .map { row ->
            def meta = [id: row.patient_id, timepoint: row.timepoint]
            tuple(meta, file(row.rna_bam))
        }
    fastq_rows_ch = rna_samplesheet_ch
        .filter { row -> row.rna_fastq_r1 && row.rna_fastq_r2 }
        .map { row ->
            def meta = [id: row.patient_id, timepoint: row.timepoint]
            tuple(meta, file(row.rna_fastq_r1), file(row.rna_fastq_r2))
        }

    BAM_TO_FASTQ(bam_rows_ch)
    fastq_ch = BAM_TO_FASTQ.out.fastq.mix(fastq_rows_ch)

    fusions_ch = Channel.empty()

    if (run_starfusion) {
        ctat_lib = Channel.fromPath(params.ctat_resource_lib).collect()
        STAR_FUSION(fastq_ch, ctat_lib)
        fusions_ch = fusions_ch.mix(STAR_FUSION.out.fusions)
    }

    if (run_arriba) {
        fasta = Channel.fromPath(params.genome_fasta).collect()
        gtf   = Channel.fromPath(params.gtf).collect()
        star_index = Channel.fromPath(params.star_index_dir).collect()

        ARRIBA_STAR_ALIGN(fastq_ch, star_index)
        ARRIBA(ARRIBA_STAR_ALIGN.out.bam, fasta, gtf)
        fusions_ch = fusions_ch.mix(ARRIBA.out.fusions)
    }

    // fusions_ch: [meta, caller, fusion_file] — one entry per
    // patient-timepoint-caller combination requested above
    AGFUSION_ANNOTATE(fusions_ch)

    if (!params.hla_alleles_manual) {
        log.warn "FUSION_NEOANTIGENS: no --hla_alleles_manual given and this workflow doesn't yet consume PVACSEQ_CORE's typed HLA output — pvacfuse will need HLA alleles supplied another way for now"
    }
    hla_ch = AGFUSION_ANNOTATE.out.annotated.map { meta, caller, dir -> tuple(meta, caller, dir, params.hla_alleles_manual) }

    PVACFUSE_RUN(hla_ch)

    emit:
    results = PVACFUSE_RUN.out.results   // [meta, caller, pvacfuse_output]
}
