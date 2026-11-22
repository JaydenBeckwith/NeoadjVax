#!/usr/bin/env nextflow
/*
 * ===========================================================================
 * NeoadjVax: somatic neoantigen prediction for neoadjuvant melanoma vaccines
 * ===========================================================================
 *
 *   tumour + normal DNA (FASTQ or BAM)
 *     ─> PREPROCESS_DNA     QC, trim, BWA-MEM, MarkDuplicates, BQSR
 *     ─> SOMATIC_CALLING    Mutect2 + orientation + contamination, FilterMutectCalls, PASS
 *     ─> HLA_TYPING         HLA-HD on the normal (or samplesheet hla_alleles), SpecHLA LOH read in
 *     ─> PROXIMAL_VARIANTS  germline near somatic sites, WhatsHap phasing, VEP
 *     ─> NEOANTIGENS        VEP, bam-readcount (DNA + optional RNA), pVACseq (Sequenza purity), CCF (ploidy/CN)
 *     ─> EXPRESSION         optional: featureCounts TPM, GTEx/TCGA context, presentation-machinery and driver panel, shared antigens
 *     ─> FUSION             optional: STAR-Fusion on baseline RNA FASTQ, AGFusion, pVACfuse
 *     ─> PRIORITIZE         one ranked list per MHC class: LOH, clonality, expression, agretopicity as flags with reasons
 *     ─> LONG_PEPTIDES      mutant protein context, 25mers to order, optional pVACvector construct
 *     ─> CLONAL_STRUCTURE   CCF clusters and week 6 persistence
 *     ─> REPORT             RNA support of variants/neoantigens per timepoint, melanoma patient report (LaTeX/PDF), per-patient workbook
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
include { EXPRESSION        } from './workflows/expression'
include { FUSION            } from './workflows/fusion'
include { REPORT            } from './workflows/report'
include { PRIORITIZE } from './modules/local/pvactools/prioritize'
include { GENERATE_PROTEIN_FASTA; DESIGN_LONG_PEPTIDES; PVACVECTOR } from './modules/local/pvactools/long_peptides'
include { WORKBOOK } from './modules/local/report/workbook'
include { PVACBIND } from './modules/local/expression/taa'
include { CLONAL_STRUCTURE } from './modules/local/report/clonal_structure'
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

        // RNA BAMs only feed bam-readcount; they are not realigned, so they
        // must already be on the same GRCh38 reference, with an index.
        ['rna_baseline_bam', 'rna_week6_bam'].findAll { row[it] }.each { col ->
            if (!file(row[col]).exists()) error "samplesheet: ${id} ${col} not found: ${row[col]}"
            if (!bamIndex(row[col])) error "samplesheet: ${id} ${col} has no index (${row[col]}.bai): run samtools index on it"
        }

        // Fusion inputs: a baseline tumour RNA FASTQ pair for STAR-Fusion,
        // or the full predictions table of an existing STAR-Fusion run.
        def hasFusionR1 = row.rna_baseline_r1 as boolean
        def hasFusionR2 = row.rna_baseline_r2 as boolean
        if (hasFusionR1 != hasFusionR2) error "samplesheet: ${id} needs both rna_baseline_r1 and rna_baseline_r2, or neither"
        if (hasFusionR1 && row.starfusion_tsv) error "samplesheet: ${id} gives both an RNA FASTQ pair and starfusion_tsv: supply one or the other"
        [row.rna_baseline_r1, row.rna_baseline_r2, row.starfusion_tsv].findAll { it }.each { f ->
            if (!file(f).exists()) error "samplesheet: ${id} fusion input not found: ${f}"
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
        // SpecHLA HLA-LOH result (merge.hla.copy.txt), computed outside
        // this pipeline from the tumour BAM and the Sequenza purity/ploidy.
        if (row.hla_loh && !file(row.hla_loh).exists()) {
            error "samplesheet: ${id} hla_loh not found: ${row.hla_loh}"
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

// <x>.bam.bai or <x>.bai beside a BAM, or null.
def bamIndex(path) {
    def candidates = [file("${path}.bai")]
    if (path.endsWith('.bam')) candidates << file(path[0..-5] + '.bai')
    candidates.find { it.exists() }
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

    // ---- optional expression context: each reference table needs the
    //      group name to report, and counting needs a GTF ----
    if (params.gtex_db && !params.gtex_tissue) error "--gtex_db needs --gtex_tissue (a GROUP in that table)"
    if (params.tcga_db && !params.tcga_cancer_type) error "--tcga_db needs --tcga_cancer_type (a GROUP in that table)"
    gtex_db = params.gtex_db ? file(params.gtex_db, checkIfExists: true) : []
    tcga_db = params.tcga_db ? file(params.tcga_db, checkIfExists: true) : []
    rna_gtf = params.rna_gtf ? file(params.rna_gtf, checkIfExists: true) : []
    if (!(params.rna_strandedness in [0, 1, 2])) error "--rna_strandedness must be 0 (unstranded), 1 (forward) or 2 (reverse)"

    // ---- optional fusion branch: the AGFusion DB, its pyensembl cache and
    //      the CTAT lib are checked here, because compute nodes are offline
    //      and AGFusion would otherwise try to download a missing cache ----
    def fusion_rows = params.fusion ? rows.findAll { it.rna_baseline_r1 || it.starfusion_tsv } : []
    agfusion_db = []
    ctat_genome_lib = []
    if (params.fusion) {
        if (!fusion_rows) log.warn "--fusion: no samplesheet row has rna_baseline_r1/_r2 or starfusion_tsv, so no fusions will be called"
        agfusion_db = file(requireParam('agfusion_db'), checkIfExists: true)
        def db_name = agfusion_db.name =~ /^agfusion\.homo_sapiens\.(\d+)\.db$/
        if (!db_name.matches()) error "--agfusion_db: expected a file named agfusion.homo_sapiens.<Ensembl release>.db, got ${agfusion_db.name}"
        def release = db_name.group(1)
        def cache_dir = file("${requireParam('agfusion_pyensembl_cache')}/pyensembl/GRCh38/ensembl${release}")
        if (!(cache_dir.isDirectory() && cache_dir.listFiles().any { it.name.endsWith('.db') })) {
            error "--agfusion_pyensembl_cache: no indexed pyensembl ${release} cache under ${cache_dir} (build it on a login node with: PYENSEMBL_CACHE_DIR=${params.agfusion_pyensembl_cache} pyensembl install --release ${release} --species homo_sapiens)"
        }
        if (params.agfusion_pythonpath && !file("${params.agfusion_pythonpath}/sitecustomize.py").exists()) {
            error "--agfusion_pythonpath: no sitecustomize.py in ${params.agfusion_pythonpath} (set it to '' if the AGFusion image needs no shim)"
        }
        if (release != params.vep_cache_version.toString()) {
            log.warn "--agfusion_db is Ensembl ${release} but the VEP cache is ${params.vep_cache_version}: fusion and SNV/indel candidates will use different gene models"
        }
        if (fusion_rows.any { it.rna_baseline_r1 }) {
            ctat_genome_lib = file(requireParam('ctat_genome_lib'), checkIfExists: true)
            if (!file("${ctat_genome_lib}/ref_genome.fa.star.idx").isDirectory()) {
                error "--ctat_genome_lib: ${ctat_genome_lib} has no ref_genome.fa.star.idx; point it at the ctat_genome_lib_build_dir inside the plug-n-play bundle"
            }
        }
    } else if (rows.any { it.rna_baseline_r1 || it.starfusion_tsv }) {
        log.warn "Samplesheet has fusion inputs (rna_baseline_r1/_r2 or starfusion_tsv) but --fusion is false: ignoring them"
    }

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

    // ---- optional tumour RNA BAMs, one list per patient (may be empty):
    //      [meta, [[sample_type, bam, bai], ...]] ----
    rna_ch = samples_ch.map { row ->
        tuple([id: row.patient_id],
              ['rna_baseline', 'rna_week6']
                  .findAll { row["${it}_bam"] }
                  .collect { st -> [st, file(row["${st}_bam"]), bamIndex(row["${st}_bam"])] })
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
                HLA_TYPING.out.alleles, clinical_ch, rna_ch, fasta, fai, vep_cache, vep_plugins)

    // Expression context runs when there is anything to compare against:
    // a reference table, or the patient's own RNA with a GTF to count it
    // against. With neither there is nothing to say, and the CCF reports
    // are the final candidate output.
    expression_ch = Channel.empty()
    panel_ch = Channel.empty()
    taa_ch = Channel.empty()
    if (gtex_db || tcga_db || rna_gtf) {
        if (!rna_gtf) log.warn "No --rna_gtf: skipping RNA quantification, so expression context comes from the reference tables only"
        if (!(gtex_db || tcga_db)) log.warn "No --gtex_db or --tcga_db: expression comes from this patient's RNA alone, with no normal-tissue or cohort context"
        EXPRESSION(rna_ch, NEOANTIGENS.out.ccf, NEOANTIGENS.out.vep_vcf, clinical_ch,
                   rna_gtf, gtex_db, tcga_db)
        expression_ch = EXPRESSION.out.annotated
        panel_ch = EXPRESSION.out.panel
        taa_ch = EXPRESSION.out.taa
    }

    // Binding predictions for a shared-antigen peptide set, which depends
    // only on the patient's alleles. Needs a peptide FASTA: none ships
    // here, because which epitopes are worth testing is a judgement.
    if (params.taa && params.taa_peptides) {
        PVACBIND(HLA_TYPING.out.alleles, file(params.taa_peptides, checkIfExists: true))
    }

    // One ranked list per MHC class, from whichever annotation the run
    // produced: the expression-annotated tables when EXPRESSION ran,
    // otherwise the CCF tables. HLA-LOH joins in for the patients that
    // have it (remainder: true keeps the rest, with no LOH column).
    prioritize_in_ch = (gtex_db || tcga_db || rna_gtf ? expression_ch : NEOANTIGENS.out.ccf)
        .join(HLA_TYPING.out.loh, remainder: true)
        .filter { it[1] != null }
        .map { it[-1] == null ? it[0..-2] + [[]] : it }
    PRIORITIZE(prioritize_in_ch)

    // Long peptides for the ranked candidates. The flanking sequence is in
    // no table, so pVACtools regenerates the mutant protein context from
    // the same VCF pVACseq read, with the same phasing.
    long_peptides_ch = Channel.empty()
    if (params.long_peptides) {
        GENERATE_PROTEIN_FASTA(NEOANTIGENS.out.readcounts_vcf.join(PROXIMAL_VARIANTS.out.phased_vcf))
        DESIGN_LONG_PEPTIDES(PRIORITIZE.out.tsv.join(GENERATE_PROTEIN_FASTA.out.fasta))
        long_peptides_ch = DESIGN_LONG_PEPTIDES.out.tsv

        // One string-of-beads construct per patient, from the Combined
        // peptide set where there is one (both MHC classes), else whatever
        // single class the run produced.
        if (params.vaccine_vector) {
            PVACVECTOR(DESIGN_LONG_PEPTIDES.out.fasta
                .map { meta, fastas ->
                    def list = fastas instanceof List ? fastas : [fastas]
                    tuple(meta, list.find { it.name.contains('Combined') } ?: list.sort { it.name }.first())
                }
                .join(HLA_TYPING.out.alleles))
        }
    }

    // Clonal clusters, and which survived to week 6 in the RNA.
    clonal_ch = Channel.empty()
    if (params.clonal_structure) {
        CLONAL_STRUCTURE(PRIORITIZE.out.tsv)
        clonal_ch = CLONAL_STRUCTURE.out.clusters
    }

    // One spreadsheet per patient, from whatever the run produced. Each
    // optional table joins with remainder: true and arrives as [] when
    // the patient has none, so a DNA-only patient still gets a workbook.
    if (params.workbook) {
        workbook_in_ch = PRIORITIZE.out.tsv
            .join(long_peptides_ch, remainder: true)
            .join(HLA_TYPING.out.loh, remainder: true)
            .join(panel_ch, remainder: true)
            .join(taa_ch, remainder: true)
            .join(clonal_ch, remainder: true)
            .filter { it[1] != null }
            .map { row -> row[0..0] + row[1..-1].collect { it == null ? [] : it } }
        WORKBOOK(workbook_in_ch)
    }

    // Fusion neoantigens are opt-in and need RNA reads (or an existing
    // STAR-Fusion table), so most patients may have none.
    pvacfuse_ch = Channel.empty()
    fusion_calls_ch = Channel.empty()
    if (params.fusion) {
        fusion_samples_ch = Channel.fromList(fusion_rows)
        FUSION(fusion_samples_ch
                   .filter { row -> row.rna_baseline_r1 }
                   .map { row -> tuple([id: row.patient_id], file(row.rna_baseline_r1), file(row.rna_baseline_r2)) },
               fusion_samples_ch
                   .filter { row -> row.starfusion_tsv }
                   .map { row -> tuple([id: row.patient_id], file(row.starfusion_tsv)) },
               HLA_TYPING.out.alleles, ctat_genome_lib, agfusion_db)
        pvacfuse_ch = FUSION.out.pvacfuse
        fusion_calls_ch = FUSION.out.predictions
    }

    // RNA-support figure for patients with RNA, and a report for every patient
    if (params.report) {
        REPORT(NEOANTIGENS.out.readcounts_vcf, NEOANTIGENS.out.rna_week6_vcf, rna_ch,
               NEOANTIGENS.out.ccf, expression_ch, panel_ch, HLA_TYPING.out.alleles, clinical_ch,
               pvacfuse_ch, fusion_calls_ch, fasta, fai, target_intervals)
    }

    workflow.onComplete = {
        log.info "NeoadjVax finished: status: ${workflow.success ? 'OK' : 'FAILED'}: results: ${params.outdir}"
    }
}
