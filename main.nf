#!/usr/bin/env nextflow
/*
 * ===========================================================================
 * NeoadjVax: somatic neoantigen prediction for neoadjuvant melanoma vaccines
 * ===========================================================================
 *
 *   tumour + normal DNA (FASTQ or BAM)
 *     ─> PREPROCESS_DNA     QC, trim, BWA-MEM, MarkDuplicates, BQSR
 *     ─> SOMATIC_CALLING    Mutect2 + orientation + contamination, FilterMutectCalls, PASS
 *     ─> HLA_TYPING         HLA-HD on the normal (or samplesheet hla_alleles)
 *     ─> PROXIMAL_VARIANTS  germline near somatic sites, WhatsHap phasing, VEP
 *     ─> NEOANTIGENS        VEP, bam-readcount, pVACseq (Sequenza purity), CCF (ploidy/CN)
 *
 * One samplesheet row per patient: see assets/samplesheet_schema.md, and
 * docs/ARCHITECTURE.md for why each step is there.
 */

nextflow.enable.dsl = 2

include { PREPROCESS_DNA    } from './workflows/preprocess_dna'
include { SOMATIC_CALLING   } from './workflows/somatic_calling'
include { HLA_TYPING        } from './workflows/hla_typing'
include { PROXIMAL_VARIANTS } from './workflows/proximal_variants'
include { NEOANTIGENS       } from './workflows/neoantigens'
include { SAMTOOLS_FAIDX; CREATE_SEQUENCE_DICTIONARY; BWA_INDEX } from './modules/local/reference/prepare_genome'

// ---------------------------------------------------------------------
// Launch-time checks: the samplesheet and every reference file are
// validated before any task is submitted, so a typo fails in seconds
// rather than after a day of alignment.
// ---------------------------------------------------------------------
def isNumber(value) {
    try { Double.parseDouble(value as String); return true } catch (Exception e) { return false }
}

def validateSamplesheet(path) {
    def rows = file(path, checkIfExists: true).splitCsv(header: true)
    if (!rows) error "samplesheet ${path} has no rows"
    def seen = [] as Set
    rows.each { row ->
        def id = row.patient_id
        if (!id) error "samplesheet: a row has no patient_id"
        if (!(id ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) error "samplesheet: patient_id '${id}' must start with a letter or digit and contain only letters, digits, '.', '_' or '-'"
        if (!seen.add(id)) error "samplesheet: patient_id '${id}' appears more than once"

        ['tumor', 'normal'].each { st ->
            def hasBam = row["${st}_bam"] as boolean
            def hasR1 = row["${st}_r1"] as boolean
            def hasR2 = row["${st}_r2"] as boolean
            if (hasBam && (hasR1 || hasR2)) error "samplesheet: ${id} gives both ${st}_bam and ${st} FASTQ: supply one or the other"
            if (!hasBam && !(hasR1 && hasR2)) error "samplesheet: ${id} needs ${st}_r1 + ${st}_r2, or ${st}_bam"
            [row["${st}_bam"], row["${st}_r1"], row["${st}_r2"]].findAll { it }.each { f ->
                if (!file(f).exists()) error "samplesheet: ${id} ${st} input not found: ${f}"
            }
        }

        if (!row.purity) error "samplesheet: ${id} has no purity (Sequenza cellularity, between 0 and 1)"
        if (!isNumber(row.purity) || !((row.purity as double) > 0 && (row.purity as double) <= 1)) {
            error "samplesheet: ${id} purity '${row.purity}' must be a number in (0, 1]"
        }
        if (row.ploidy && !(isNumber(row.ploidy) && (row.ploidy as double) > 0)) {
            error "samplesheet: ${id} ploidy '${row.ploidy}' must be a positive number"
        }
        if (row.sequenza_segments && !file(row.sequenza_segments).exists()) {
            error "samplesheet: ${id} sequenza_segments not found: ${row.sequenza_segments}"
        }
        if (row.sex && !(row.sex.toLowerCase() in ['male', 'female'])) {
            error "samplesheet: ${id} sex '${row.sex}' must be male or female (or left blank)"
        }
        if (!row.hla_alleles && !params.hlahd_container) {
            error "samplesheet: ${id} has no hla_alleles, so HLA-HD must type it, but --hlahd_container is not set: build containers/hlahd (see its README) or give hla_alleles"
        }
    }
    rows
}

// A VCF param that must exist with its .tbi beside it; returns [vcf, tbi].
def indexedVcf(name, path) {
    def vcf = file(path)
    def tbi = file("${path}.tbi")
    if (!vcf.exists()) error "--${name}: not found: ${path}"
    if (!tbi.exists()) error "--${name}: no index at ${path}.tbi (create one with: gatk IndexFeatureFile -I ${path})"
    [vcf, tbi]
}

def requireParam(name) {
    if (!params[name]) error "--${name} is required (see docs/DEPENDENCIES.md)"
    params[name]
}

workflow {

    rows = validateSamplesheet(requireParam('samplesheet'))
    samples_ch = Channel.fromList(rows)

    // ---- reference files: reuse indexes found beside the FASTA ----
    def fasta_path = requireParam('genome_fasta')
    fasta = file(fasta_path, checkIfExists: true)
    fai = file("${fasta_path}.fai").exists() ? file("${fasta_path}.fai") : SAMTOOLS_FAIDX(fasta).fai
    def dict_path = fasta_path.replaceAll(/\.(fa|fasta|fna)(\.gz)?$/, '.dict')
    dict = file(dict_path).exists() ? file(dict_path) : CREATE_SEQUENCE_DICTIONARY(fasta).dict
    def bwa_suffixes = ['amb', 'ann', 'bwt', 'pac', 'sa']
    bwa_index = bwa_suffixes.every { file("${fasta_path}.${it}").exists() }
        ? bwa_suffixes.collect { file("${fasta_path}.${it}") }
        : BWA_INDEX(fasta).index

    def known = ['known_dbsnp', 'known_mills', 'known_1000g'].collect { indexedVcf(it, requireParam(it)) }
    germline_resource = indexedVcf('germline_resource', requireParam('germline_resource'))
    pon_files = params.panel_of_normals ? indexedVcf('panel_of_normals', params.panel_of_normals) : []
    common_sites = params.contamination_sites ? indexedVcf('contamination_sites', params.contamination_sites) : null
    target_intervals = params.target_intervals ? file(params.target_intervals, checkIfExists: true) : []
    vep_cache = file(requireParam('vep_cache'), checkIfExists: true)
    vep_plugins = file(requireParam('vep_plugins'), checkIfExists: true)

    if (!params.panel_of_normals) log.warn "No --panel_of_normals: Mutect2 will run without a PoN, so recurrent artefacts are not filtered"
    if (!params.contamination_sites) log.warn "No --contamination_sites: skipping GetPileupSummaries/CalculateContamination"

    // ---- per-patient clinical inputs, kept out of meta so editing a
    //      purity value does not invalidate the cached alignment ----
    clinical_ch = samples_ch.map { row ->
        tuple([id: row.patient_id],
              row.purity as String,
              row.ploidy ?: '',
              row.sex ? row.sex.toLowerCase() : '',
              row.sequenza_segments ? file(row.sequenza_segments) : [])
    }

    PREPROCESS_DNA(samples_ch, fasta, fai, dict, bwa_index, known*.get(0), known*.get(1))
    bams_ch = PREPROCESS_DNA.out.bam
    tumor_bam_ch  = bams_ch.filter { it[1] == 'tumor' }.map { meta, st, bam, bai -> tuple(meta, bam, bai) }
    normal_bam_ch = bams_ch.filter { it[1] == 'normal' }.map { meta, st, bam, bai -> tuple(meta, bam, bai) }

    SOMATIC_CALLING(tumor_bam_ch.join(normal_bam_ch), bams_ch, fasta, fai, dict,
                    germline_resource, pon_files, common_sites, target_intervals)

    HLA_TYPING(samples_ch, normal_bam_ch)

    PROXIMAL_VARIANTS(SOMATIC_CALLING.out.pass_vcf, tumor_bam_ch, normal_bam_ch,
                      fasta, fai, dict, vep_cache, vep_plugins)

    NEOANTIGENS(SOMATIC_CALLING.out.pass_vcf, bams_ch, PROXIMAL_VARIANTS.out.phased_vcf,
                HLA_TYPING.out.alleles, clinical_ch, fasta, fai, vep_cache, vep_plugins)

    workflow.onComplete = {
        log.info "NeoadjVax finished: status: ${workflow.success ? 'OK' : 'FAILED'}: results: ${params.outdir}"
    }
}
