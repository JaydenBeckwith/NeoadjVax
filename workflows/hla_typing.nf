/*
 * HLA alleles per patient: the samplesheet's hla_alleles column when
 * given, otherwise HLA-HD on the normal BAM (MHC reads only).
 */

include { EXTRACT_HLA_READS; HLAHD } from '../modules/local/hla/hlahd'

workflow HLA_TYPING {

    take:
    samples_ch      // samplesheet rows
    normal_bam_ch   // [meta, bam, bai]

    main:
    // '|' or ';' in the samplesheet (a comma would split the CSV column)
    manual_ch = samples_ch
        .filter { row -> row.hla_alleles }
        .map { row -> tuple([id: row.patient_id], row.hla_alleles.tokenize('|;')*.trim().join(',')) }

    to_type_ch = samples_ch
        .filter { row -> !row.hla_alleles }
        .map { row -> tuple([id: row.patient_id], true) }

    EXTRACT_HLA_READS(normal_bam_ch.join(to_type_ch).map { meta, bam, bai, _flag -> tuple(meta, bam, bai) })
    HLAHD(EXTRACT_HLA_READS.out.reads)

    typed_ch = HLAHD.out.alleles.map { meta, txt -> tuple(meta, txt.text.trim()) }

    emit:
    alleles = manual_ch.mix(typed_ch)   // [meta, 'HLA-A*02:01,...,DQA1*01:02-DQB1*06:02']
}
