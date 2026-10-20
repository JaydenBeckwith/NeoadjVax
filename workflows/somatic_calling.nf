/*
 * Mutect2 tumour vs matched normal, GATK best practice:
 *
 *   Mutect2 (per chromosome) -> MergeVcfs + MergeMutectStats + LearnReadOrientationModel
 *   GetPileupSummaries (tumour, normal) -> CalculateContamination   [if params.contamination_sites]
 *   -> FilterMutectCalls -> PASS only, left-aligned, one ALT per record
 */

include { MUTECT2; MERGE_MUTECT2_SCATTERS; GET_PILEUP_SUMMARIES; CALCULATE_CONTAMINATION; FILTER_MUTECT_CALLS } from '../modules/local/dna/mutect2'
include { SELECT_PASS_SOMATIC } from '../modules/local/vcf/bcftools'

workflow SOMATIC_CALLING {

    take:
    pairs_ch          // [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
    bams_ch           // [meta, sample_type, bam, bai]
    fasta
    fasta_fai
    fasta_dict
    germline_resource // [vcf.gz, tbi]
    pon_files         // [vcf.gz, tbi] or []
    common_sites      // [vcf.gz, tbi] or null
    target_intervals  // file or []

    main:
    scatter_ch = Channel.fromList(params.scatter_chromosomes.tokenize(',')*.trim())
    MUTECT2(pairs_ch.combine(scatter_ch), fasta, fasta_fai, fasta_dict, germline_resource, pon_files, target_intervals)

    // regroup the per-chromosome outputs; size: lets each patient go on as
    // soon as its own scatters finish instead of waiting for the cohort
    n_scatter = params.scatter_chromosomes.tokenize(',').size()
    per_patient_ch = MUTECT2.out.vcf.groupTuple(size: n_scatter)
        .join(MUTECT2.out.stats.groupTuple(size: n_scatter))
        .join(MUTECT2.out.f1r2.groupTuple(size: n_scatter))
    MERGE_MUTECT2_SCATTERS(per_patient_ch)

    if (common_sites) {
        GET_PILEUP_SUMMARIES(bams_ch, common_sites, target_intervals)
        pileups_ch = GET_PILEUP_SUMMARIES.out.table
            .map { meta, st, table -> tuple(meta, [(st): table]) }
            .groupTuple(size: 2)
            .map { meta, maps -> def m = maps.sum(); tuple(meta, m.tumor, m.normal) }
        CALCULATE_CONTAMINATION(pileups_ch)
        contamination_ch = CALCULATE_CONTAMINATION.out.tables.map { meta, c, s -> tuple(meta, [c, s]) }
    } else {
        contamination_ch = pairs_ch.map { meta, tb, ti, nb, ni -> tuple(meta, []) }
    }

    FILTER_MUTECT_CALLS(MERGE_MUTECT2_SCATTERS.out.merged.join(contamination_ch), fasta, fasta_fai, fasta_dict)
    SELECT_PASS_SOMATIC(FILTER_MUTECT_CALLS.out.vcf, fasta, fasta_fai)

    emit:
    filtered_vcf = FILTER_MUTECT_CALLS.out.vcf     // [meta, vcf.gz, tbi], all records with FILTER set
    pass_vcf     = SELECT_PASS_SOMATIC.out.vcf     // [meta, vcf.gz, tbi], PASS, normalised
}
