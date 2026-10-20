/*
 * Phased VCF of proximal variants for pVACseq, following pVACtools'
 * "Creating a phased VCF of proximal variants" recipe with WhatsHap in
 * place of the retired GATK3 ReadBackedPhasing:
 *
 *   PASS somatic ─> windows ─> HaplotypeCaller on normal (windows only)
 *        │                              │
 *        └──────── tumour column + germline (renamed to tumour) ─> WhatsHap
 *                  on the tumour BAM ─> VEP ─> bgzip + tabix
 */

include { SOMATIC_WINDOWS_BED; COMBINE_SOMATIC_GERMLINE; BGZIP_INDEX } from '../modules/local/vcf/bcftools'
include { GERMLINE_HAPLOTYPECALLER; WHATSHAP_PHASE } from '../modules/local/phasing/phasing'
include { VEP_ANNOTATE } from '../modules/local/dna/vep_annotate'

workflow PROXIMAL_VARIANTS {

    take:
    somatic_vcf_ch    // [meta, vcf.gz, tbi] PASS, normalised
    tumor_bam_ch      // [meta, bam, bai]
    normal_bam_ch     // [meta, bam, bai]
    fasta
    fasta_fai
    fasta_dict
    vep_cache
    vep_plugins

    main:
    SOMATIC_WINDOWS_BED(somatic_vcf_ch)
    GERMLINE_HAPLOTYPECALLER(normal_bam_ch.join(SOMATIC_WINDOWS_BED.out.bed), fasta, fasta_fai, fasta_dict)

    COMBINE_SOMATIC_GERMLINE(somatic_vcf_ch.join(GERMLINE_HAPLOTYPECALLER.out.vcf), fasta, fasta_fai)
    WHATSHAP_PHASE(COMBINE_SOMATIC_GERMLINE.out.vcf.join(tumor_bam_ch), fasta, fasta_fai)

    VEP_ANNOTATE(WHATSHAP_PHASE.out.vcf.map { meta, vcf -> tuple(meta, vcf, 'phased') }, fasta, fasta_fai, vep_cache, vep_plugins)
    BGZIP_INDEX(VEP_ANNOTATE.out.vcf.map { meta, vcf -> tuple(meta, vcf, 'phased_proximal') })

    emit:
    phased_vcf = BGZIP_INDEX.out.vcf   // [meta, vcf.gz, tbi]
}
