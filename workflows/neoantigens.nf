/*
 * Somatic VCF to ranked neoantigen candidates:
 *
 *   PASS somatic ─> VEP ─> bam-readcount (tumour, normal DNA; RNA baseline
 *     and week 6 when given) ─> vcf-readcount-annotator ─> bgzip + tabix
 *     ─> pVACseq (HLA, phased proximal VCF, purity) ─> CCF + week 6 RNA columns
 *
 * pVACseq reads one RNA sample (RAD/RDP/RAF in the tumour column), so the
 * baseline RNA counts go into its input VCF and drive its RNA filters. The
 * week 6 counts are written into a separate copy of the VCF and joined onto
 * the aggregated reports as extra columns afterwards.
 */

include { VEP_ANNOTATE } from '../modules/local/dna/vep_annotate'
include { BAM_READCOUNT; VCF_READCOUNT_ANNOTATOR } from '../modules/local/pvactools/readcount'
include { VCF_READCOUNT_ANNOTATOR as RNA_WEEK6_ANNOTATOR } from '../modules/local/pvactools/readcount'
include { BGZIP_INDEX } from '../modules/local/vcf/bcftools'
include { PVACSEQ; ANNOTATE_CCF } from '../modules/local/pvactools/pvacseq'

workflow NEOANTIGENS {

    take:
    somatic_vcf_ch    // [meta, vcf.gz, tbi] PASS, normalised
    bams_ch           // [meta, sample_type, bam, bai]
    phased_vcf_ch     // [meta, vcf.gz, tbi]
    hla_ch            // [meta, allele string]
    clinical_ch       // [meta, purity, ploidy, sex, segments]
    rna_ch            // [meta, [[sample_type, bam, bai], ...]] per patient, list may be empty
    fasta
    fasta_fai
    vep_cache
    vep_plugins

    main:
    VEP_ANNOTATE(somatic_vcf_ch.map { meta, vcf, tbi -> tuple(meta, vcf, 'somatic') }, fasta, fasta_fai, vep_cache, vep_plugins)

    // RNA is counted at the sites in the tumour's genotype column
    readcount_in_ch = bams_ch
        .map { meta, st, bam, bai -> tuple(meta, st, st, bam, bai) }
        .mix(rna_ch.flatMap { meta, rna -> rna.collect { st, bam, bai -> tuple(meta, st, 'tumor', bam, bai) } })
        .combine(VEP_ANNOTATE.out.vcf, by: 0)   // [meta, sample_type, vcf_sample, bam, bai, vcf]
    BAM_READCOUNT(readcount_in_ch, fasta, fasta_fai)

    // [meta, [sample_type: [snv, indel]]] once every BAM of the patient is counted
    expected_ch = rna_ch.map { meta, rna -> tuple(meta, 2 + rna.size()) }
    readcounts_ch = BAM_READCOUNT.out.readcounts
        .combine(expected_ch, by: 0)
        .map { meta, st, snv, indel, n -> tuple(groupKey(meta, n), [(st): [snv, indel]]) }
        .groupTuple()
        .map { key, maps -> tuple(key.getGroupTarget(), maps.sum()) }

    annotator_in_ch = VEP_ANNOTATE.out.vcf.join(readcounts_ch)
    VCF_READCOUNT_ANNOTATOR(annotator_in_ch.map { meta, vcf, counts ->
        tuple(meta, vcf, ['tumor', 'normal', 'rna_baseline'].findAll { counts[it] }.collectMany { counts[it] }, 'readcounts')
    })
    BGZIP_INDEX(VCF_READCOUNT_ANNOTATOR.out.vcf.map { meta, vcf -> tuple(meta, vcf, 'somatic.vep.readcounts') })

    RNA_WEEK6_ANNOTATOR(annotator_in_ch
        .filter { meta, vcf, counts -> counts.rna_week6 }
        .map { meta, vcf, counts -> tuple(meta, vcf, counts.rna_week6, 'rna_week6.readcounts') })

    pvacseq_in_ch = BGZIP_INDEX.out.vcf
        .join(phased_vcf_ch)
        .join(hla_ch)
        .join(clinical_ch.map { meta, purity, ploidy, sex, segments -> tuple(meta, purity) })
    PVACSEQ(pvacseq_in_ch)

    ANNOTATE_CCF(PVACSEQ.out.results
        .join(clinical_ch)
        .join(RNA_WEEK6_ANNOTATOR.out.vcf, remainder: true)
        .map { it[-1] == null ? it[0..-2] + [[]] : it })

    emit:
    pvacseq = PVACSEQ.out.results
    ccf     = ANNOTATE_CCF.out.tsv
    rna_week6_vcf = RNA_WEEK6_ANNOTATOR.out.vcf
}
