/*
 * Tumor purity/ploidy branch — optional, independently runnable
 * (--run_purity_ploidy, default off). Ported from
 * https://github.com/JaydenBeckwith/Sequenza_tools's four PBS scripts
 * (sequenza_step1-4.sh): per-chromosome bam2seqz + seqz_binning, merge
 * into one per-sample file, then a Sequenza R fit. HLA_LOH
 * (workflows/hla_loh.nf) and the biallelic-loss/whole-genome-doubling
 * question from [[loh-analysis]] both need this branch's purity/ploidy
 * output.
 *
 * What changed vs. the original four scripts, and why:
 *   - Step 1's manual bash `&`/`wait` parallel loop (24 chromosomes inside
 *     one 24-cpu/64GB PBS job) becomes 24 separate Nextflow tasks — a
 *     better fit for this pipeline's per-step PBS-job-per-task model (see
 *     docs/RUNNING_ON_GADI.md) than replicating manual shell backgrounding.
 *     See modules/local/purity_ploidy/sequenza_bam2seqz_binned.nf.
 *   - Step 3 (sequenza_step3.sh, "collect merged files into one
 *     directory") is dropped — it exists purely to move files between two
 *     PBS jobs sharing a filesystem; Nextflow's own channel staging does
 *     that automatically between SEQUENZA_MERGE_BINS and SEQUENZA_FIT, so
 *     there's nothing left for it to do here.
 *   - Steps 2 and 4 still need the actual merge-bin200-files.pl /
 *     merge-header-bin200.pl / run-sequenza.R script bodies — you sent the
 *     PBS wrapper scripts that call them, not the scripts themselves. The
 *     orchestration around them (env setup, sex lookup, decompression,
 *     exact CLI invocations) is a faithful, complete port; the scripts
 *     themselves are checked-for-and-error-if-missing rather than guessed
 *     at. See modules/local/purity_ploidy/sequenza_merge_bins.nf and
 *     sequenza_fit.nf.
 *   - Extracting a single (purity, ploidy) pair out of whatever
 *     run-sequenza.R writes is a genuinely new step none of the four
 *     scripts do (NeoadjLOH's sequenza_top_solutions_summary.csv implies
 *     something aggregates a "top solution" per sample, but that isn't
 *     among the scripts you sent either) — left as an explicit stub, see
 *     modules/local/purity_ploidy/sequenza_extract_top_solution.nf.
 *
 * Sex lookup: ported exactly as sequenza_step2.sh/step4.sh do it — a
 * numeric "melpin" prefix pulled from the sample id, looked up against
 * column 2 (melpin) / column 7 (gender) of --sequenza_gender_csv (their
 * dna_neotrio_gender_metadata.csv). A sample with no match is skipped with
 * a warning, matching the original's behaviour exactly (not a hard error).
 *
 * Runs on tumor+normal DNA BAM pairs, matched by patient_id — either
 * chained from DNA_VARIANT_CALLING's tumor_bam/normal_bam outputs, or
 * (for standalone use, per "can be run separately") built directly from
 * the dna_samplesheet's BAM columns without running the rest of
 * DNA_VARIANT_CALLING. See main.nf for how the two are picked between.
 */

include { SEQUENZA_BAM2SEQZ_BINNED }        from '../modules/local/purity_ploidy/sequenza_bam2seqz_binned'
include { SEQUENZA_MERGE_BINS }             from '../modules/local/purity_ploidy/sequenza_merge_bins'
include { SEQUENZA_FIT }                    from '../modules/local/purity_ploidy/sequenza_fit'
include { SEQUENZA_EXTRACT_TOP_SOLUTION }   from '../modules/local/purity_ploidy/sequenza_extract_top_solution'

// Ported from sequenza_step2.sh/step4.sh's get_gender() — reads
// --sequenza_gender_csv (column 2 = numeric "melpin" id, column 7 =
// gender; header row skipped), matching the exact awk logic those scripts
// use. Returns a melpin -> lowercased-gender map.
def loadSequenzaGenderCsv(path) {
    def lookup = [:]
    def rows = file(path).splitCsv()
    rows.drop(1).each { cols ->
        if (cols.size() >= 7) {
            def melpin = cols[1]?.toString()?.replaceAll('\r', '')?.replaceAll('"', '')?.trim()
            def gender = cols[6]?.toString()?.replaceAll('\r', '')?.replaceAll('"', '')?.trim()?.toLowerCase()
            if (melpin) lookup[melpin] = gender
        }
    }
    lookup
}

workflow PURITY_PLOIDY {

    take:
    tumor_bam_ch    // [meta, tumor_bam, tumor_bai]
    normal_bam_ch   // [meta, normal_bam, normal_bai]

    main:
    if (!params.sequenza_gender_csv) {
        error "PURITY_PLOIDY needs --sequenza_gender_csv (sequenza_step2.sh/step4.sh's dna_neotrio_gender_metadata.csv format: column 2 = numeric melpin id matching the leading digits of the sample id, column 7 = gender, header row present)"
    }
    def gender_lookup = loadSequenzaGenderCsv(params.sequenza_gender_csv)

    bin_size = params.sequenza_seq_type == 'wes' ? params.sequenza_bin_wes : params.sequenza_bin_wgs

    fasta      = Channel.fromPath(params.genome_fasta).collect()
    fasta_fai  = Channel.fromPath("${params.genome_fasta}.fai").collect()
    gc_file    = Channel.fromPath(params.sequenza_gc_file).collect()
    seqz_header = Channel.fromPath(params.sequenza_seqz_header).collect()

    pairs_with_sex_ch = tumor_bam_ch
        .map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }
        .join(normal_bam_ch.map { meta, bam, bai -> tuple(meta.id, bam, bai) }, by: 0)
        .map { patient_id, meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            tuple(meta, tumor_bam, tumor_bai, normal_bam, normal_bai)
        }
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            def m = (meta.id =~ /^(\d+)/)
            def melpin = m.find() ? m.group(1) : null
            def sex = melpin ? gender_lookup[melpin] : null
            tuple(meta, tumor_bam, tumor_bai, normal_bam, normal_bai, sex)
        }

    // matches the original: no gender match -> skip that sample with a
    // warning, not a pipeline-wide error
    pairs_with_sex_ch
        .filter { meta, tb, ti, nb, ni, sex -> !sex }
        .subscribe { meta, tb, ti, nb, ni, sex ->
            log.warn "PURITY_PLOIDY: could not find gender for ${meta.id} in --sequenza_gender_csv (no melpin match) — skipping, matching sequenza_step2.sh/step4.sh's own behaviour"
        }

    pairs_ok_ch = pairs_with_sex_ch.filter { meta, tb, ti, nb, ni, sex -> sex }

    // ---- Step 1: per-chromosome bam2seqz + seqz_binning ----
    chrom_list = (1..22).collect { "chr${it}" } + ['chrX', 'chrY']
    chrom_ch = Channel.fromList(chrom_list)

    per_chr_ch = pairs_ok_ch
        .map { meta, tb, ti, nb, ni, sex -> tuple(meta, tb, ti, nb, ni) }
        .combine(chrom_ch)

    SEQUENZA_BAM2SEQZ_BINNED(per_chr_ch, fasta, fasta_fai, gc_file, bin_size)

    // ---- Step 2: merge the 24 per-chromosome files back into one ----
    binned_grouped_ch = SEQUENZA_BAM2SEQZ_BINNED.out.binned_seqz
        .map { meta, chrom, binned -> tuple(meta.id, meta, binned) }
        .groupTuple(by: 0)
        .join(pairs_ok_ch.map { meta, tb, ti, nb, ni, sex -> tuple(meta.id, sex) }, by: 0)
        .map { patient_id, metas, binned_files, sex -> tuple(metas[0], binned_files, sex, bin_size) }

    SEQUENZA_MERGE_BINS(binned_grouped_ch, seqz_header)

    // ---- Step 4: Sequenza R fit ----
    SEQUENZA_FIT(SEQUENZA_MERGE_BINS.out.merged_seqz)

    // ---- purity/ploidy extraction (STUB) ----
    SEQUENZA_EXTRACT_TOP_SOLUTION(SEQUENZA_FIT.out.raw_results)

    emit:
    purity_ploidy = SEQUENZA_EXTRACT_TOP_SOLUTION.out.purity_ploidy   // [meta, purity, ploidy]
}
