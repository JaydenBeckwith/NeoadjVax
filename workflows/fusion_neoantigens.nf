/*
 * Caller-native RNA fusion calls -> AGFusion -> pVACfuse, per patient/timepoint/caller.
 * Raw reads use independent STAR passes for STAR-Fusion and Arriba. Precomputed
 * caller TSVs bypass only their corresponding caller. No cross-caller consensus
 * or RNA timepoint restriction is imposed. See docs/FUSION.md.
 */
include { BAM_TO_FASTQ }       from '../modules/local/rna/bam_to_fastq'
include { STAR_FUSION }        from '../modules/local/fusion/star_fusion'
include { ARRIBA_STAR_ALIGN } from '../modules/local/fusion/arriba_align'
include { ARRIBA }             from '../modules/local/fusion/arriba'
include { AGFUSION_ANNOTATE }  from '../modules/local/fusion/agfusion_annotate'
include { PVACFUSE_RUN }       from '../modules/local/pvactools/pvacfuse'

def fusionAlleles(value) {
    def alleles = (value instanceof List ? value.join(',') : (value ?: '').toString())
        .split(/[|,]/).collect { it.trim() }.findAll { it }.unique().sort()
    if (!alleles || alleles.any { !(it ==~ /[A-Za-z0-9*:_-]+/) }) {
        error 'Fusion: HLA alleles must be non-empty pVACtools-format names, separated by | or commas'
    }
    alleles.join(',')
}

def fusionPath(value, label, directory = false) {
    if (!value || value.toString().contains("'") || (value.toString() =~ /[\r\n]/)) {
        error "Fusion: provide a valid ${label} path (no apostrophes/newlines)"
    }
    def path = file(value)
    if (!path.exists() || (directory ? !path.isDirectory() : !path.isFile())) {
        error "Fusion: ${label} does not exist or is not a ${directory ? 'directory' : 'file'}: ${value}"
    }
    path
}

// Eager launch validation, before expensive alignments or HLA typing.
def validateFusionInputs(rnaSheet, dnaSheet, useCoreHla, p) {
    def callers = (p.fusion_callers ?: '').toString().split(',').collect { it.trim().toLowerCase() }.findAll { it }
    if (!callers || callers.size() != callers.unique(false).size() || callers.any { !(it in ['starfusion', 'arriba']) }) {
        error "Fusion: --fusion_callers must contain starfusion and/or arriba without duplicates"
    }
    if (!p.agfusion_container) error 'Fusion: set --agfusion_container (containers/agfusion/Dockerfile)'
    def db = fusionPath(p.agfusion_database, '--agfusion_database')
    if (!(db.name ==~ /agfusion\.homo_sapiens\.[0-9]+\.db/)) {
        error 'Fusion: keep the reference filename agfusion.homo_sapiens.<release>.db'
    }
    fusionPath(p.agfusion_pyensembl_cache, '--agfusion_pyensembl_cache', true)
    def algorithms = (p.pvacfuse_algorithms ?: '').toString().split(',').collect { it.trim() }.findAll { it }
    if (!algorithms || algorithms.unique(false).size() != algorithms.size() ||
        algorithms.any { !(it in ['MHCflurry', 'NetMHCpan', 'NetMHCIIpan']) }) {
        error 'Fusion: --pvacfuse_algorithms supports MHCflurry,NetMHCpan,NetMHCIIpan without duplicates'
    }
    if ('MHCflurry' in algorithms) {
        def models = fusionPath(p.pvacfuse_mhcflurry_models, '--pvacfuse_mhcflurry_models', true)
        fusionPath(models.resolve('manifest.csv'), 'MHCflurry manifest.csv')
    }
    if (algorithms.any { it in ['NetMHCpan', 'NetMHCIIpan'] }) {
        def iedb = fusionPath(p.pvacfuse_iedb_install_dir, '--pvacfuse_iedb_install_dir', true)
        fusionPath(iedb.resolve('mhc_i/src/predict_binding.py'), 'IEDB MHC I executable')
        fusionPath(iedb.resolve('mhc_ii/mhc_II_binding.py'), 'IEDB MHC II executable')
    }
    ['fusion_min_read_support', 'pvacfuse_threads', 'pvacfuse_binding_threshold'].each { key ->
        if (!(p[key].toString() ==~ /[0-9]+/) || (key != 'fusion_min_read_support' && p[key].toString().toInteger() < 1)) {
            error "Fusion: invalid --${key}"
        }
    }
    def ffpm = p.fusion_min_ffpm.toString()
    if (!ffpm.isBigDecimal() || ffpm.toBigDecimal() < 0) error 'Fusion: --fusion_min_ffpm must be finite and non-negative'
    ['pvacfuse_epitope_lengths_i', 'pvacfuse_epitope_lengths_ii'].each { key ->
        if (!(p[key].toString() ==~ /[1-9][0-9]*(,[1-9][0-9]*)*/)) error "Fusion: invalid --${key}"
    }
    def rows = file(rnaSheet).splitCsv(header: true)
    if (!rows) error 'Fusion: RNA samplesheet is empty'
    def keys = [] as Set
    def patientHla = [:]
    def needsCalling = [] as Set
    rows.each { row ->
        ['patient_id', 'timepoint'].each { key ->
            if (!(row[key] ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
                error "Fusion: ${key} must use letters, digits, dots, underscores or hyphens"
            }
        }
        def key = [row.patient_id, row.timepoint]
        if (!keys.add(key)) error "Fusion: duplicate patient/timepoint ${key}"
        if (row.hla_alleles) {
            def hla = fusionAlleles(row.hla_alleles)
            if (patientHla[row.patient_id] && patientHla[row.patient_id] != hla) {
                error "Fusion: inconsistent HLA alleles across timepoints for ${row.patient_id}"
            }
            patientHla[row.patient_id] = hla
        }
        def missing = callers.findAll { !row["${it}_tsv"] }
        callers.findAll { row["${it}_tsv"] }.each { fusionPath(row["${it}_tsv"], "${it}_tsv for ${key}") }
        if (missing) {
            needsCalling.addAll(missing)
            if ((row.rna_fastq_r1 as boolean) != (row.rna_fastq_r2 as boolean)) {
                error "Fusion: incomplete FASTQ pair for ${key}"
            }
            if (row.rna_fastq_r1 && row.rna_fastq_r2) {
                fusionPath(row.rna_fastq_r1, "RNA R1 for ${key}")
                fusionPath(row.rna_fastq_r2, "RNA R2 for ${key}")
                if (row.rna_fastq_r1.endsWith('.gz') != row.rna_fastq_r2.endsWith('.gz')) {
                    error "Fusion: FASTQ mates must use the same compression for ${key}"
                }
            } else if (row.rna_bam) {
                fusionPath(row.rna_bam, "RNA BAM for ${key}")
            } else {
                error "Fusion: ${key} needs RNA BAM/paired FASTQ or ${missing} caller TSVs; an RNA VCF is not a fusion input"
            }
        }
    }
    def patients = rows.collect { it.patient_id }.unique()
    if (p.hla_alleles_manual) {
        if (patients.size() != 1) error 'Fusion: global --hla_alleles_manual is restricted to single-patient runs; use RNA-sheet hla_alleles for cohorts'
        def hla = fusionAlleles(p.hla_alleles_manual)
        if (patientHla[patients[0]] && patientHla[patients[0]] != hla) error 'Fusion: RNA-sheet HLA conflicts with --hla_alleles_manual'
        patientHla[patients[0]] = hla
    }
    def corePatients = useCoreHla ? file(dnaSheet).splitCsv(header: true).collect { it.patient_id } : []
    if (corePatients.unique(false).size() != corePatients.size()) error 'Fusion: duplicate DNA patients make core HLA reuse ambiguous'
    patients.each { patient ->
        if (!patientHla[patient] && !(patient in corePatients)) {
            error "Fusion: no HLA source for ${patient}; supply RNA-sheet hla_alleles or run the somatic core for this patient"
        }
    }
    if ('starfusion' in needsCalling) fusionPath(p.ctat_resource_lib, '--ctat_resource_lib', true)
    if ('arriba' in needsCalling) {
        fusionPath(p.genome_fasta, '--genome_fasta')
        fusionPath(p.gtf, '--gtf')
        fusionPath(p.star_index_dir, '--star_index_dir', true)
        // A missing blacklist is not an innocuous default for candidate discovery.
        fusionPath(p.arriba_blacklist, '--arriba_blacklist')
        if (p.arriba_known_fusions) fusionPath(p.arriba_known_fusions, '--arriba_known_fusions')
        if (p.arriba_protein_domains) fusionPath(p.arriba_protein_domains, '--arriba_protein_domains')
    }
    [callers: callers, needsCalling: needsCalling, patientHla: patientHla]
}

workflow FUSION_NEOANTIGENS {
    take:
    rna_samplesheet_ch
    core_hla_by_patient_ch
    settings

    main:
    rows_ch = rna_samplesheet_ch.map { row ->
        def meta = [id: row.patient_id, timepoint: row.timepoint,
                    callers: settings.callers.findAll { !row["${it}_tsv"] }]
        tuple(meta, row)
    }
    bam_ch = rows_ch.filter { meta, row -> meta.callers && row.rna_bam && !row.rna_fastq_r1 }
        .map { meta, row -> tuple(meta, file(row.rna_bam, checkIfExists: true)) }
    direct_fastq_ch = rows_ch.filter { meta, row -> meta.callers && row.rna_fastq_r1 && row.rna_fastq_r2 }
        .map { meta, row -> tuple(meta, file(row.rna_fastq_r1, checkIfExists: true), file(row.rna_fastq_r2, checkIfExists: true)) }
    BAM_TO_FASTQ(bam_ch)
    fastq_ch = BAM_TO_FASTQ.out.fastq.mix(direct_fastq_ch)
    fusions_ch = rows_ch.flatMap { meta, row ->
        settings.callers.findAll { row["${it}_tsv"] }.collect { caller ->
            tuple(meta, caller, file(row["${caller}_tsv"], checkIfExists: true))
        }
    }
    if ('starfusion' in settings.needsCalling) {
        STAR_FUSION(fastq_ch.filter { meta, r1, r2 -> 'starfusion' in meta.callers },
                    Channel.value(file(params.ctat_resource_lib, checkIfExists: true)))
        fusions_ch = fusions_ch.mix(STAR_FUSION.out.fusions)
    }
    if ('arriba' in settings.needsCalling) {
        ARRIBA_STAR_ALIGN(fastq_ch.filter { meta, r1, r2 -> 'arriba' in meta.callers },
                          Channel.value(file(params.star_index_dir, checkIfExists: true)))
        ARRIBA(ARRIBA_STAR_ALIGN.out.bam,
               Channel.value(file(params.genome_fasta, checkIfExists: true)),
               Channel.value(file(params.gtf, checkIfExists: true)),
               Channel.value(file(params.arriba_blacklist, checkIfExists: true)),
               Channel.value(params.arriba_known_fusions ? file(params.arriba_known_fusions, checkIfExists: true) : []),
               Channel.value(params.arriba_protein_domains ? file(params.arriba_protein_domains, checkIfExists: true) : []))
        fusions_ch = fusions_ch.mix(ARRIBA.out.fusions)
    }
    AGFUSION_ANNOTATE(fusions_ch,
                     Channel.value(file(params.agfusion_database, checkIfExists: true)),
                     Channel.value(file(params.agfusion_pyensembl_cache, checkIfExists: true)))

    // A single map broadcasts patient HLA to EVERY timepoint and caller.
    // Unlike an inner join, lookup fails if typing never emitted a patient.
    hla_map_ch = core_hla_by_patient_ch.toList().map { entries ->
        def result = [:]
        entries.each { patient, alleles ->
            if (result.containsKey(patient)) error "Fusion: duplicate core HLA output for ${patient}"
            result[patient] = fusionAlleles(alleles)
        }
        result
    }
    hla_ch = AGFUSION_ANNOTATE.out.annotated.combine(hla_map_ch)
        .map { meta, caller, dir, hlaMap ->
            def alleles = settings.patientHla[meta.id] ?: hlaMap[meta.id]
            if (!alleles) error "Fusion: HLA typing produced no alleles for ${meta.id}"
            tuple(meta, caller, dir, alleles)
        }
    PVACFUSE_RUN(hla_ch,
                 Channel.value(params.pvacfuse_mhcflurry_models ? file(params.pvacfuse_mhcflurry_models, checkIfExists: true) : []),
                 Channel.value(params.pvacfuse_iedb_install_dir ? file(params.pvacfuse_iedb_install_dir, checkIfExists: true) : []))
    emit:
    results = PVACFUSE_RUN.out.results
    annotations = AGFUSION_ANNOTATE.out.annotated
}
