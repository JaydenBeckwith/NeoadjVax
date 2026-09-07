include { SPLICE_REFERENCE; SPLICEAI_ANNOTATE; SPLICEAI_FILTER; SPLICE_PREDICT; SPLICE_TO_PEPTIDE } from '../modules/local/splicing/splice2neo'
include { SPLICING_STAR_INDEX; SPLICING_STAR_ALIGN; LEAFCUTTER; LEAFCUTTER_SELECT } from '../modules/local/splicing/leafcutter'

def validateSplicingSamplesheets(dnaPath, rnaPath, reuseDna, opts) {
    def dnaRows = file(dnaPath, checkIfExists: true).splitCsv(header: true)
    def rnaRows = file(rnaPath, checkIfExists: true).splitCsv(header: true)
    if (!dnaRows || !rnaRows) error "Splicing requires nonempty DNA and RNA samplesheets"
    def ids = [] as Set
    def rnaKeys = [] as Set
    def safeId = { value -> value && (value ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/) }
    def checkFile = { value ->
        if (value.toString().contains("'") || value.toString().contains('\n')) error "Splicing paths cannot contain single quotes or newlines"
        file(value, checkIfExists: true)
    }
    dnaRows.each { row ->
        if (!safeId(row.patient_id)) error "Splicing DNA rows need a safe patient_id (letters, digits, _, ., -)"
        if (!ids.add(row.patient_id)) error "Splicing DNA: duplicate patient_id ${row.patient_id}"
        if (!row.spliceai_vcf && !row.dna_vcf && !reuseDna) {
            error "Splicing DNA ${row.patient_id}: supply dna_vcf/spliceai_vcf or also select dna_variant_calling"
        }
        if (row.spliceai_vcf) checkFile(row.spliceai_vcf)
        else if (row.dna_vcf) checkFile(row.dna_vcf)
    }
    rnaRows.each { row ->
        if (!safeId(row.patient_id) || !safeId(row.timepoint)) error "Splicing RNA rows need safe patient_id and timepoint labels"
        if (!ids.contains(row.patient_id)) error "Splicing RNA patient ${row.patient_id} has no matching DNA row"
        if (!rnaKeys.add([row.patient_id, row.timepoint])) error "Duplicate splicing RNA patient/timepoint: ${row.patient_id}/${row.timepoint}"
        if (row.leafcutter_counts) {
            checkFile(row.leafcutter_counts)
            if (row.leafcutter_sample && !(row.leafcutter_sample ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) error "leafcutter_sample must be an exact, safe counts-header sample name"
        } else if (row.rna_bam) {
            checkFile(row.rna_bam)
            if (!((row.rna_strand ?: opts.splicing_rna_strand) in ['XS', 'RF', 'FR'])) {
                error "Splicing RNA BAM ${row.patient_id}/${row.timepoint}: set rna_strand or --splicing_rna_strand to XS, RF, or FR"
            }
        } else if (row.rna_fastq_r1 && row.rna_fastq_r2) {
            checkFile(row.rna_fastq_r1)
            checkFile(row.rna_fastq_r2)
            if (row.rna_fastq_r1.endsWith('.gz') != row.rna_fastq_r2.endsWith('.gz')) error "RNA FASTQ pair must use the same compression"
        } else {
            error "Splicing RNA ${row.patient_id}/${row.timepoint}: needs leafcutter_counts, an unsplit RNA BAM, or paired FASTQs; RNA VCF alone is not junction evidence"
        }
    }
    if (opts.spliceai_annotated_vcf_dir) error "Replace --spliceai_annotated_vcf_dir with explicit per-patient spliceai_vcf paths in the DNA sheet"
    checkFile(opts.genome_fasta)
    checkFile(opts.gtf)
    if (opts.splicing_normal_junctions) checkFile(opts.splicing_normal_junctions)
    else if (!opts.splicing_allow_missing_normal) error "Supply --splicing_normal_junctions (GTEx/normal junc_id TSV), or explicitly enable --splicing_allow_missing_normal for exploratory output"
    def needsAnnotation = dnaRows.any { !it.spliceai_vcf }
    def needsLeafcutter = rnaRows.any { !it.leafcutter_counts }
    if ((needsAnnotation || needsLeafcutter) && !opts.splicing_tools_container) error "Set --splicing_tools_container; build recipe: containers/splicing-tools/Dockerfile"
    if (needsAnnotation) {
        if (!opts.spliceai_annotation) error "DNA SpliceAI annotation requires --spliceai_annotation grch37/grch38 or a matching custom annotation file"
        if (!(opts.spliceai_annotation in ['grch37', 'grch38'])) checkFile(opts.spliceai_annotation)
    }
    if (rnaRows.any { !it.leafcutter_counts && !it.rna_bam } && opts.star_index_dir) checkFile(opts.star_index_dir)
    if (!(opts.spliceai_min_delta_score instanceof Number) || opts.spliceai_min_delta_score < 0 || opts.spliceai_min_delta_score > 1) error "spliceai_min_delta_score must be in [0,1]"
    if (!(opts.spliceai_distance instanceof Number) || opts.spliceai_distance < 1 || opts.spliceai_distance > 4999 || opts.spliceai_distance != opts.spliceai_distance.intValue()) error "spliceai_distance must be an integer in [1,4999]"
    ['splicing_min_rna_reads', 'splicing_peptide_flank', 'splicing_min_peptide_length', 'leafcutter_min_reads', 'splicing_star_overhang'].each { key ->
        if (!(opts[key] instanceof Number) || opts[key] < 1 || opts[key] != opts[key].intValue()) error "${key} must be a positive integer"
    }
    if (!(opts.leafcutter_min_ratio instanceof Number) || opts.leafcutter_min_ratio < 0 || opts.leafcutter_min_ratio > 1) error "leafcutter_min_ratio must be in [0,1]"
}

workflow SPLICING_NEOANTIGENS {
    take:
    dna_samplesheet_ch
    rna_samplesheet_ch
    called_dna_vcf_ch
    reuse_dna_calling

    main:
    fasta = Channel.value(file(params.genome_fasta, checkIfExists: true))
    gtf = Channel.value(file(params.gtf, checkIfExists: true))
    normal = Channel.value(params.splicing_normal_junctions ? file(params.splicing_normal_junctions, checkIfExists: true) : [])
    SPLICE_REFERENCE(fasta, gtf, normal)

    preannotated = dna_samplesheet_ch.filter { it.spliceai_vcf }
        .map { row -> tuple([id: row.patient_id], file(row.spliceai_vcf, checkIfExists: true)) }
    if (reuse_dna_calling) {
        wanted = dna_samplesheet_ch.filter { !it.spliceai_vcf }.map { tuple(it.patient_id, [id: it.patient_id]) }
        raw_dna = called_dna_vcf_ch.map { meta, vcf -> tuple(meta.id, vcf) }
            .combine(wanted, by: 0).map { id, vcf, meta -> tuple(meta, vcf) }
    } else {
        raw_dna = dna_samplesheet_ch.filter { !it.spliceai_vcf }
            .map { row -> tuple([id: row.patient_id], file(row.dna_vcf, checkIfExists: true)) }
    }
    annotation_file = Channel.value(params.spliceai_annotation && !(params.spliceai_annotation in ['grch37', 'grch38']) ? file(params.spliceai_annotation, checkIfExists: true) : [])
    SPLICEAI_ANNOTATE(raw_dna, fasta, SPLICE_REFERENCE.out.fai, annotation_file)
    SPLICEAI_FILTER(preannotated.mix(SPLICEAI_ANNOTATE.out.vcf))
    SPLICE_PREDICT(SPLICEAI_FILTER.out.table, SPLICE_REFERENCE.out.reference)

    supplied_counts = rna_samplesheet_ch.filter { it.leafcutter_counts }.map { row ->
        tuple([id: row.patient_id, timepoint: row.timepoint], file(row.leafcutter_counts, checkIfExists: true), row.leafcutter_sample ?: '')
    }
    supplied_bams = rna_samplesheet_ch.filter { !it.leafcutter_counts && it.rna_bam }.map { row ->
        tuple([id: row.patient_id, timepoint: row.timepoint, strand: row.rna_strand ?: params.splicing_rna_strand], file(row.rna_bam, checkIfExists: true))
    }
    fastqs = rna_samplesheet_ch.filter { !it.leafcutter_counts && !it.rna_bam }.map { row ->
        tuple([id: row.patient_id, timepoint: row.timepoint, strand: 'XS'], file(row.rna_fastq_r1, checkIfExists: true), file(row.rna_fastq_r2, checkIfExists: true))
    }
    // Alignment for this branch preserves N CIGARs and XS strand tags.
    // Never use GATK SplitNCigarReads/ApplyBQSR output for junction discovery.
    hasFastqs = file(params.rna_samplesheet).splitCsv(header: true).any { !it.leafcutter_counts && !it.rna_bam }
    aligned = Channel.empty()
    if (hasFastqs) {
        if (params.star_index_dir) {
            star_index = Channel.value(file(params.star_index_dir, checkIfExists: true))
        } else {
            SPLICING_STAR_INDEX(fasta, gtf)
            star_index = SPLICING_STAR_INDEX.out.index
        }
        SPLICING_STAR_ALIGN(fastqs, star_index)
        aligned = SPLICING_STAR_ALIGN.out.bam
    }
    LEAFCUTTER(supplied_bams.mix(aligned))
    LEAFCUTTER_SELECT(supplied_counts.mix(LEAFCUTTER.out.counts.map { meta, counts -> tuple(meta, counts, '') }))

    // One baseline DNA prediction must fan out to ALL arbitrary RNA timepoints.
    // combine(by:0), not join(), implements this deliberate one-to-many relation.
    per_timepoint = SPLICE_PREDICT.out.predictions.map { meta, predictions -> tuple(meta.id, predictions) }
        .combine(LEAFCUTTER_SELECT.out.selected.map { meta, counts, evidence -> tuple(meta.id, meta, counts, evidence) }, by: 0)
        .map { id, predictions, meta, counts, evidence -> tuple(meta, predictions, counts, evidence) }
    SPLICE_TO_PEPTIDE(per_timepoint, SPLICE_REFERENCE.out.reference, fasta, SPLICE_REFERENCE.out.fai)

    emit:
    peptide_fasta = SPLICE_TO_PEPTIDE.out.peptide_fasta
    candidates = SPLICE_TO_PEPTIDE.out.candidates
    evidence = SPLICE_TO_PEPTIDE.out.evidence
    summary = SPLICE_TO_PEPTIDE.out.summary
}
