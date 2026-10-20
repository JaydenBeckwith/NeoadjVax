/*
 * Somatic VCF to ranked neoantigen candidates:
 *
 *   PASS somatic ─> VEP ─> bam-readcount (tumour, normal) ─> vcf-readcount-annotator
 *     ─> bgzip + tabix ─> pVACseq (HLA, phased proximal VCF, purity) ─> CCF columns
 */

include { VEP_ANNOTATE } from '../modules/local/dna/vep_annotate'
include { BAM_READCOUNT; VCF_READCOUNT_ANNOTATOR } from '../modules/local/pvactools/readcount'
include { BGZIP_INDEX } from '../modules/local/vcf/bcftools'
include { PVACSEQ; ANNOTATE_CCF } from '../modules/local/pvactools/pvacseq'

workflow NEOANTIGENS {

    take:
    somatic_vcf_ch    // [meta, vcf.gz, tbi] PASS, normalised
    bams_ch           // [meta, sample_type, bam, bai]
    phased_vcf_ch     // [meta, vcf.gz, tbi]
    hla_ch            // [meta, allele string]
    clinical_ch       // [meta, purity, ploidy, sex, segments]
    fasta
    fasta_fai
    vep_cache
    vep_plugins

    main:
    VEP_ANNOTATE(somatic_vcf_ch.map { meta, vcf, tbi -> tuple(meta, vcf, 'somatic') }, fasta, fasta_fai, vep_cache, vep_plugins)

    readcount_in_ch = bams_ch
        .combine(VEP_ANNOTATE.out.vcf, by: 0)   // [meta, sample_type, bam, bai, vcf]
    BAM_READCOUNT(readcount_in_ch, fasta, fasta_fai)

    readcounts_ch = BAM_READCOUNT.out.readcounts
        .map { meta, st, snv, indel -> tuple(meta, [(st): [snv, indel]]) }
        .groupTuple(size: 2)
        .map { meta, maps -> def m = maps.sum(); tuple(meta, m.tumor[0], m.tumor[1], m.normal[0], m.normal[1]) }
    VCF_READCOUNT_ANNOTATOR(VEP_ANNOTATE.out.vcf.join(readcounts_ch))
    BGZIP_INDEX(VCF_READCOUNT_ANNOTATOR.out.vcf.map { meta, vcf -> tuple(meta, vcf, 'somatic.vep.readcounts') })

    pvacseq_in_ch = BGZIP_INDEX.out.vcf
        .join(phased_vcf_ch)
        .join(hla_ch)
        .join(clinical_ch.map { meta, purity, ploidy, sex, segments -> tuple(meta, purity) })
    PVACSEQ(pvacseq_in_ch)

    ANNOTATE_CCF(PVACSEQ.out.results.join(clinical_ch))

    emit:
    pvacseq = PVACSEQ.out.results
    ccf     = ANNOTATE_CCF.out.tsv
}
