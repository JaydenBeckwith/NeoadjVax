/*
 * Core neoantigen-score branch: the part of the pipeline that already had
 * a design ([[neoantigen-score]]) before this repo existed:
 *
 *   baseline DNA somatic calls ──┐
 *                                 ├─▶ RNA_SUPPORT_FILTER (per RNA timepoint) ─▶ VEP_ANNOTATE ─▶ PVACSEQ_RUN
 *   per-timepoint RNA calls ─────┘                                                    ▲
 *                                                                                       │
 *   normal (germline) DNA reads ─▶ HLA_TYPING_XHLA (class I + II, once per patient) ────┘
 *                                  (skipped per-patient if hla_alleles supplied in the
 *                                   DNA samplesheet)
 *
 * One DNA VCF per patient is checked against RNA support at every timepoint
 * that patient has: this is what produces the baseline-vs-week-6(+)
 * comparison [[ines]] asked about.
 *
 * HLA typing note: xHLA types class I (A/B/C) and class II DRB1/DQB1/DPB1,
 * but only the beta chain for DQ/DP: pVACtools needs those paired with an
 * alpha chain (DQA1/DPA1) it doesn't type, so HLA_TYPING_XHLA excludes
 * DQB1/DPB1 from the allele string it hands to PVACSEQ_RUN (see that
 * module's comments and bin/parse_xhla_result.py). Only class I + DRB1
 * reach pVACseq today. Also note: pvacseq_algorithms currently defaults to
 * MHCflurry (class I only): DRB1 alleles won't actually get predicted
 * against until a class II algorithm is added to that param too.
 */

include { RNA_SUPPORT_FILTER }   from '../modules/local/matching/rna_support_filter'
include { VEP_ANNOTATE }         from '../modules/local/dna/vep_annotate'
include { HLA_TYPING_XHLA }      from '../modules/local/hla_typing/xhla'
include { PVACSEQ_RUN }          from '../modules/local/pvactools/pvacseq'

workflow PVACSEQ_CORE {

    take:
    dna_samplesheet_ch   // for the optional manual hla_alleles column
    dna_filtered_vcf_ch  // [meta(id), vcf]        <- DNA_VARIANT_CALLING.out.filtered_vcf
    dna_normal_bam_ch    // [meta(id), bam, bai]   <- DNA_VARIANT_CALLING.out.normal_bam
    rna_filtered_vcf_ch  // [meta(id, timepoint), vcf] <- RNA_VARIANT_CALLING.out.filtered_vcf

    main:
    fasta       = Channel.fromPath(params.genome_fasta).collect()
    vep_cache   = Channel.fromPath(params.vep_cache).collect()
    vep_plugins = Channel.fromPath(params.vep_plugins).collect()

    // --- join baseline DNA calls onto every RNA timepoint for that patient ---
    dna_keyed = dna_filtered_vcf_ch.map { meta, vcf -> tuple(meta.id, vcf) }
    rna_keyed = rna_filtered_vcf_ch.map { meta, vcf -> tuple(meta.id, meta, vcf) }

    matched_ch = rna_keyed.join(dna_keyed)
        .map { patient_id, rna_meta, rna_vcf, dna_vcf -> tuple(rna_meta, dna_vcf, rna_vcf) }

    RNA_SUPPORT_FILTER(matched_ch)
    VEP_ANNOTATE(RNA_SUPPORT_FILTER.out.vcf, fasta, vep_cache, vep_plugins)

    // --- HLA alleles: manual override from the DNA samplesheet, else type it ---
    // samplesheet field uses '|' as the allele separator (a literal comma
    // would be parsed as another CSV column): converted to the ','-joined
    // format pvacseq run expects, matching HLA_TYPING_OPTITYPE's output shape
    manual_hla_ch = dna_samplesheet_ch
        .filter { row -> row.hla_alleles }
        .map { row -> tuple(row.patient_id, row.hla_alleles.replace('|', ',')) }

    patients_needing_typing_ch = dna_samplesheet_ch
        .filter { row -> !row.hla_alleles }
        .map { row -> row.patient_id }

    normal_bam_for_typing_ch = dna_normal_bam_ch
        .map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }
        .join(patients_needing_typing_ch.map { pid -> tuple(pid, true) })
        .map { patient_id, meta, bam, bai, _flag -> tuple(meta, bam, bai) }

    HLA_TYPING_XHLA(normal_bam_for_typing_ch)

    typed_hla_ch = HLA_TYPING_XHLA.out.hla_alleles
        .map { meta, alleles -> tuple(meta.id, alleles) }

    hla_by_patient_ch = manual_hla_ch.mix(typed_hla_ch)

    // --- attach the right HLA string to each patient-timepoint's annotated VCF ---
    pvacseq_input_ch = VEP_ANNOTATE.out.vcf
        .map { meta, vcf -> tuple(meta.id, meta, vcf) }
        .combine(hla_by_patient_ch, by: 0)
        .map { patient_id, meta, vcf, hla_alleles -> tuple(meta, vcf, hla_alleles) }

    PVACSEQ_RUN(pvacseq_input_ch)

    emit:
    pvacseq_results = PVACSEQ_RUN.out.results
    rna_support_summary = RNA_SUPPORT_FILTER.out.summary
    hla_by_patient = hla_by_patient_ch
}
