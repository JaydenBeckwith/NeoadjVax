#!/usr/bin/env Rscript
# splice2neo v0.6.14 adapter. Sequence always comes from the supplied FASTA,
# never the package's implicit hg19 default.
suppressPackageStartupMessages({
  library(splice2neo)
  library(GenomicRanges)
  library(dplyr)
})

write_tsv <- function(x, path) {
  x <- x[, !vapply(x, is.list, logical(1)), drop = FALSE]
  readr::write_tsv(x, path, na = "NA")
}

prepare_reference <- function(gtf, fasta, normal, output) {
  Rsamtools::indexFa(fasta)
  fai <- Rsamtools::scanFaIndex(fasta)
  gr <- rtracklayer::import(gtf)
  excluded <- setdiff(GenomeInfoDb::seqlevelsInUse(gr), as.character(seqnames(fai)))
  gr <- GenomeInfoDb::keepSeqlevels(gr, intersect(GenomeInfoDb::seqlevels(gr), as.character(seqnames(fai))), pruning.mode = "coarse")
  if (!length(gr)) stop("GTF and FASTA share no contigs; check genome build and chr naming")
  if (any(end(gr) > width(fai)[match(as.character(seqnames(gr)), as.character(seqnames(fai)))])) {
    stop("GTF coordinates exceed FASTA contig lengths")
  }
  txdb <- txdbmaker::makeTxDbFromGRanges(gr)
  tx <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
  tx_gr <- GenomicFeatures::transcripts(txdb, columns = c("gene_id", "tx_id", "tx_name"))
  # annotate_mut_effect internally indexes both objects with the same index.
  tx_gr <- tx_gr[match(names(tx), tx_gr$tx_name)]
  stopifnot(identical(names(tx), as.character(tx_gr$tx_name)))
  cds <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)
  if (!length(tx) || !length(cds)) stop("GTF must contain exons and coding sequences")
  gene_table <- as.data.frame(S4Vectors::mcols(gr))
  if (!all(c("gene_id", "gene_name") %in% names(gene_table))) stop("GTF needs gene_id and gene_name attributes")
  gene_table <- distinct(gene_table, gene_id, gene_name) %>% filter(!is.na(gene_id), !is.na(gene_name))
  # Also accept ENSEMBL IDs as SYMBOL in custom SpliceAI annotations.
  gene_table <- bind_rows(gene_table, transmute(gene_table, gene_id, gene_name = gene_id)) %>% distinct()
  normals <- character()
  if (normal != "-") {
    normal_df <- readr::read_tsv(normal, col_types = readr::cols(.default = "c"))
    if (!identical(names(normal_df), "junc_id")) stop("Normal panel requires exactly one column, junc_id")
    normals <- unique(normal_df$junc_id)
    if (!length(normals) || anyNA(normals) || any(!grepl("^[^:[:space:]]+:[0-9]+-[0-9]+:[+-]$", normals))) {
      stop("Normal panel must contain nonempty chr:joinedExonBase1-joinedExonBase2:strand IDs")
    }
    nr <- splice2neo::junc_to_gr(normals)
    if (any(start(nr) < 1 | width(nr) < 2)) stop("Invalid normal-panel coordinates")
    if (!any(as.character(seqnames(nr)) %in% as.character(seqnames(fai)))) stop("Normal panel contigs do not match FASTA")
  }
  saveRDS(list(transcripts = tx, transcripts_gr = tx_gr, cds = cds,
               gene_table = gene_table, canonical = splice2neo::canonical_junctions(tx),
               normals = normals, normal_filtered = normal != "-",
               contigs = as.character(seqnames(fai))), output)
  writeLines(c(paste0("splice2neo=", packageVersion("splice2neo")),
               paste0("normal_filter=", normal != "-"),
               paste0("excluded_gtf_contigs=", paste(excluded, collapse = ","))),
             "reference_summary.txt")
}

empty_predictions <- function() {
  tibble(mut_id = character(), junc_id = character(), tx_id = character(),
         gene_id = character(), effect = character(), score = double(),
         event_type = character(), is_canonical = logical(), is_in_normal = logical())
}

predict_junctions <- function(table, reference, threshold, output) {
  ref <- readRDS(reference)
  raw <- readr::read_tsv(table, na = ".", col_types = readr::cols(
    .default = "c", POS = "i", DS_AG = "d", DS_AL = "d", DS_DG = "d", DS_DL = "d",
    DP_AG = "i", DP_AL = "i", DP_DG = "i", DP_DL = "i"))
  result <- empty_predictions()
  if (nrow(raw)) {
    if (any(!raw$CHROM %in% ref$contigs)) stop("SpliceAI VCF contigs absent from FASTA; no automatic liftover/chr renaming")
    if (any(!raw$SYMBOL %in% ref$gene_table$gene_name)) {
      stop(paste("SpliceAI genes absent from GTF:", paste(head(setdiff(raw$SYMBOL, ref$gene_table$gene_name)), collapse = ", ")))
    }
    effects <- splice2neo::format_spliceai(raw, gene_table = ref$gene_table) %>%
      filter(!is.na(pos), score >= threshold)
    if (nrow(effects)) {
      result <- splice2neo::annotate_mut_effect(
        effects, ref$transcripts, ref$transcripts_gr,
        gene_mapping = TRUE, consider_intron_retention = FALSE
      ) %>% filter(!is.na(junc_id), !is.na(tx_id))
      if (nrow(result)) {
        result <- splice2neo::unique_mut_junc(result) %>%
          select(mut_id, junc_id, tx_id, gene_id, effect, score, event_type) %>%
          mutate(is_canonical = junc_id %in% ref$canonical,
                 is_in_normal = junc_id %in% ref$normals)
      } else result <- empty_predictions()
    }
  }
  result <- arrange(result, mut_id, junc_id, tx_id)
  saveRDS(result, output)
  write_tsv(result, "predicted_junctions.tsv")
}

integrate_rna <- function(predictions, reference, fasta, counts, evidence,
                          patient, timepoint, flank, min_length) {
  ref <- readRDS(reference)
  predicted <- readRDS(predictions)
  observed <- readr::read_tsv(evidence, col_types = readr::cols(
    junc_id = "c", rna_junction_reads = "i", rna_cluster_reads = "i",
    rna_cluster_ratio = "d", leafcutter_cluster = "c"))
  # The upstream directory importer cannot select a sample from a cohort
  # matrix. Use its actual format adapter on our selected nonzero counts.
  if (nrow(observed)) {
    imported <- splice2neo::import_leafcutter_counts(counts)
    formatted <- splice2neo::transform_leafcutter_counts(imported)
    stopifnot(setequal(formatted$junc_id, observed$junc_id))
  }
  annotated <- predicted %>%
    left_join(observed, by = "junc_id") %>%
    mutate(patient_id = patient, timepoint = timepoint,
           rna_junction_reads = coalesce(rna_junction_reads, 0L),
           is_in_rnaseq = junc_id %in% observed$junc_id,
           normal_filter_applied = ref$normal_filtered,
           passes_filters = !is_canonical & !is_in_normal & is_in_rnaseq)
  write_tsv(annotated, "junction_evidence.tsv")
  supported <- filter(annotated, passes_filters)
  # Avoid upstream empty-data edge cases, including no matching CDS.
  coding <- filter(supported, tx_id %in% names(ref$cds))
  if (nrow(coding)) {
    genome <- Rsamtools::FaFile(fasta)
    coding <- splice2neo::add_peptide(coding, cds = ref$cds, bsg = genome, flanking_size = flank)
    peptide_columns <- setdiff(names(coding), names(supported))
    supported <- left_join(supported, select(coding, mut_id, junc_id, tx_id, all_of(peptide_columns)),
                           by = c("mut_id", "junc_id", "tx_id"))
  } else {
    supported$peptide_context <- rep(NA_character_, nrow(supported))
    supported$junc_in_orf <- rep(NA, nrow(supported))
  }
  supported <- supported %>%
    mutate(peptide_exported = !is.na(peptide_context) & coalesce(junc_in_orf, FALSE) &
             nchar(peptide_context) >= min_length &
             grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", peptide_context),
           candidate_id = sprintf("SPLICE_%s_%s_%06d", patient, timepoint, seq_len(n())))
  write_tsv(supported, "peptide_candidates.tsv")
  chosen <- filter(supported, peptide_exported)
  # One-to-one IDs preserve mutation/transcript provenance even when
  # several events predict the same sequence. Contexts, not ranked epitopes.
  fasta_lines <- if (nrow(chosen)) as.vector(rbind(paste0(">", chosen$candidate_id), chosen$peptide_context)) else character()
  writeLines(fasta_lines, "peptide_contexts.fasta")
  write_tsv(tibble(patient_id = patient, timepoint = timepoint,
                   predicted_junction_transcripts = nrow(predicted),
                   leafcutter_supported_junctions = nrow(observed),
                   supported_noncanonical = nrow(supported),
                   exported_peptide_contexts = nrow(chosen),
                   normal_filter_applied = ref$normal_filtered), "summary.tsv")
  writeLines(capture.output(sessionInfo()), "versions.txt")
}

main <- function(args) {
  if (as.character(packageVersion("splice2neo")) != "0.6.14") stop("This adapter requires splice2neo 0.6.14")
  switch(args[1],
    prepare = prepare_reference(args[2], args[3], args[4], args[5]),
    predict = predict_junctions(args[2], args[3], as.numeric(args[4]), args[5]),
    integrate = integrate_rna(args[2], args[3], args[4], args[5], args[6],
                              args[7], args[8], as.integer(args[9]), as.integer(args[10])),
    stop("Expected prepare, predict, or integrate"))
}

if (sys.nframe() == 0L) main(commandArgs(trailingOnly = TRUE))
