/*
 * HLA loss-of-heterozygosity (LOH) branch: optional, independently
 * runnable (--run_hla_loh, default off). Ported from
 * https://github.com/JaydenBeckwith/NeoadjLOH ("specHLA" directory),
 * which wraps SpecHLA (https://github.com/deepomicslab/SpecHLA) in three
 * stages run on the TUMOR BAM only: extract HLA reads -> type HLA ->
 * compute per-locus allelic copy ratio and call LOH, using tumor
 * purity/ploidy to interpret allelic imbalance.
 *
 * This is a separate typing call from HLA_TYPING_XHLA (used by
 * PVACSEQ_CORE, on the NORMAL bam, for pVACseq's input alleles): the two
 * are not interchangeable and cannot share output. xHLA answers "what are
 * this patient's germline HLA alleles, for neoantigen prediction". SpecHLA
 * here answers "has this tumor lost an HLA allele", which is why it types
 * from the tumor and additionally needs purity/ploidy that xHLA never
 * touches. See docs/ARCHITECTURE.md.
 *
 * Purity/ploidy source: either the (currently stubbed) PURITY_PLOIDY
 * workflow's output, or a manually-supplied --purity_ploidy_csv /
 * --hla_loh_samplesheet: see assets/samplesheet_schema.md. Whichever
 * source, the format matches NeoadjLOH's PURITY_CSV: sample,purity,ploidy
 * (comma-delimited: see the note in docs/ARCHITECTURE.md about a stale
 * comment in the original config.sh implying tab-delimited).
 */

include { SPECHLA_EXTRACT_HLA_READS } from '../modules/local/hla_loh/spechla_extract'
include { SPECHLA_TYPING }            from '../modules/local/hla_loh/spechla_typing'
include { SPECHLA_LOH }               from '../modules/local/hla_loh/spechla_loh'

workflow HLA_LOH {

    take:
    tumor_bam_ch       // [meta, tumor_bam, tumor_bai]
    purity_ploidy_ch   // [meta, purity, ploidy]

    main:
    SPECHLA_EXTRACT_HLA_READS(tumor_bam_ch)
    SPECHLA_TYPING(SPECHLA_EXTRACT_HLA_READS.out.reads)

    // join typing output to purity/ploidy on patient_id (meta.id): a
    // patient with a tumor BAM but no matching purity/ploidy row is
    // silently dropped here rather than erroring, since .join() defaults
    // to inner-join; see docs/ARCHITECTURE.md for why that's the current
    // behaviour and not yet a hard validation error.
    typing_by_id_ch = SPECHLA_TYPING.out.typing
        .map { meta, typing_out, hla_result, freq_list -> tuple(meta.id, meta, typing_out, hla_result, freq_list) }
    purity_by_id_ch = purity_ploidy_ch
        .map { meta, purity, ploidy -> tuple(meta.id, purity, ploidy) }

    loh_input_ch = typing_by_id_ch
        .join(purity_by_id_ch, by: 0)
        .map { patient_id, meta, typing_out, hla_result, freq_list, purity, ploidy ->
            tuple(meta, typing_out, hla_result, freq_list, purity, ploidy)
        }

    SPECHLA_LOH(loh_input_ch)

    emit:
    loh_calls = SPECHLA_LOH.out.loh_calls   // [meta, merge.hla.copy.txt]
}
