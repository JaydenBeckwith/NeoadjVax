/*
 * Non-reference ERV/TE insertions in DNA. Tumour and matched normal are
 * called independently; these are not somatic calls or peptide predictions.
 */
include { ERVCALLER_INDEX } from '../modules/local/erv/ervcaller_index'
include { ERVCALLER_RUN }   from '../modules/local/erv/ervcaller'

def ervDnaFile(path, String description) {
    def result = file(path.toString(), checkIfExists: true)
    if (!(result instanceof java.nio.file.Path) || !java.nio.file.Files.isRegularFile(result)) {
        error "ERV_DNA: ${description} must name one existing file: ${path}"
    }
    result
}

def ervDnaIndex(String bam) {
    def candidates = [file("${bam}.bai"), file(bam.replaceFirst(/(?i)\.bam$/, '.bai'))]
    def index = candidates.find { java.nio.file.Files.isRegularFile(it) }
    if (!index) error "ERV_DNA: missing BAM index for ${bam}; provide .bam.bai or .bai"
    index
}

def validateErvDnaSamplesheet(path, boolean reuseDnaBams) {
    if (!params.ervcaller_container) error "ERV_DNA needs --ervcaller_container; see docs/ERV_DNA.md"
    if (!(params.ervcaller_script ==~ /\/[A-Za-z0-9_\/.-]+\.pl/)) {
        error "ERV_DNA: --ervcaller_script must be an absolute in-container Perl path without whitespace or shell characters"
    }
    ['ervcaller_read_length', 'ervcaller_min_reads', 'ervcaller_min_split'].each { key ->
        if (params[key] == null || !(params[key].toString() ==~ /[1-9][0-9]*/)) {
            error "ERV_DNA: --${key} must be a positive integer (read length must be supplied explicitly)"
        }
    }
    if (params.ervcaller_min_split.toInteger() > params.ervcaller_read_length.toInteger()) {
        error "ERV_DNA: --ervcaller_min_split cannot exceed --ervcaller_read_length"
    }
    if (params.ervcaller_reference_dir) {
        def refs = file(params.ervcaller_reference_dir.toString(), checkIfExists: true)
        ['genome.fa', 'te.fa'].each { name ->
            ['', '.amb', '.ann', '.bwt', '.pac', '.sa', '.fai'].each { suffix ->
                ervDnaFile("${refs}/${name}${suffix}", 'prepared reference component')
            }
        }
    } else {
        if (!params.genome_fasta || !params.ervcaller_te_fasta) {
            error "ERV_DNA needs --genome_fasta and --ervcaller_te_fasta, or --ervcaller_reference_dir"
        }
        ervDnaFile(params.genome_fasta, 'genome FASTA')
        ervDnaFile(params.ervcaller_te_fasta, 'TE FASTA')
    }
    def ids = [] as Set
    def count = 0
    file(path, checkIfExists: true).splitCsv(header: true).each { row ->
        count++
        if (!(row.patient_id ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
            error "ERV_DNA: patient_id is required; use only letters, digits, underscore, hyphen and dot"
        }
        if (!ids.add(row.patient_id)) error "ERV_DNA: duplicate patient_id ${row.patient_id}"
        ['tumor', 'normal'].each { role ->
            // VCF rows bypass DNA_VARIANT_CALLING, so use their BAMs directly.
            if (!reuseDnaBams || row.dna_vcf) {
                def bam = row["${role}_bam"]
                if (!bam || !bam.toLowerCase().endsWith('.bam')) {
                    error "ERV_DNA: ${row.patient_id} requires ${role}_bam; VCFs contain no insertion-supporting reads. For FASTQs also select dna_variant_calling."
                }
                ervDnaFile(bam, "${role} BAM")
                ervDnaIndex(bam)
            } else {
                def bam = row["${role}_bam"]
                if (bam) {
                    ervDnaFile(bam, "${role} BAM")
                } else {
                    if (!row["${role}_r1"] || !row["${role}_r2"]) {
                        error "ERV_DNA: ${row.patient_id} needs ${role}_bam or a complete FASTQ pair"
                    }
                    ervDnaFile(row["${role}_r1"], "${role} R1")
                    ervDnaFile(row["${role}_r2"], "${role} R2")
                }
            }
        }
    }
    if (!count) error "ERV_DNA: dna_samplesheet contains no patients"
}

workflow ERV_DNA {
    take:
    dna_samplesheet_ch
    tumor_bam_ch
    normal_bam_ch
    reuse_dna_bams

    main:
    direct_bams = dna_samplesheet_ch
        .filter { row -> !reuse_dna_bams || row.dna_vcf }
        .flatMap { row ->
            ['tumor', 'normal'].collect { role ->
                def bam = row["${role}_bam"]
                tuple([id: row.patient_id], role, file(bam), ervDnaIndex(bam))
            }
        }
    sample_bams = direct_bams
        .mix(tumor_bam_ch.map { meta, bam, bai -> tuple(meta, 'tumor', bam, bai) })
        .mix(normal_bam_ch.map { meta, bam, bai -> tuple(meta, 'normal', bam, bai) })

    if (params.ervcaller_reference_dir) {
        references = Channel.value(file(params.ervcaller_reference_dir, checkIfExists: true))
    } else {
        ERVCALLER_INDEX(
            Channel.value(file(params.genome_fasta, checkIfExists: true)),
            Channel.value(file(params.ervcaller_te_fasta, checkIfExists: true))
        )
        references = ERVCALLER_INDEX.out.references.first()
    }
    ERVCALLER_RUN(sample_bams, references, Channel.value(file("${moduleDir}/../bin/run_ervcaller.py")))

    emit:
    calls = ERVCALLER_RUN.out.calls
    evidence = ERVCALLER_RUN.out.evidence
    status = ERVCALLER_RUN.out.status
}
