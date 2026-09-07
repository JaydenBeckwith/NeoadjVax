#!/usr/bin/env Rscript
# Usage: Rscript tests/splicing/test_splice2neo.R /absolute/generated-fixture
args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(".", winslash = "/", mustWork = TRUE)
fixture <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
source(file.path(repo, "bin/run_splice2neo.R"))
stopifnot(as.character(packageVersion("splice2neo")) == "0.6.14")
out <- tempfile("splice2neo-tests-")
dir.create(out)
setwd(out)
file.copy(file.path(fixture, "genome.fa"), "genome.fa")
prepare_reference(file.path(fixture, "genes.gtf"), "genome.fa", file.path(fixture, "normal.tsv"), "reference.rds")
ref <- readRDS("reference.rds")
stopifnot(all(c("chr1:160-301:+", "chr1:760-901:-") %in% ref$canonical))
predict_junctions(file.path(fixture, "effects.tsv"), "reference.rds", 0.5, "predictions.rds")
pred <- readRDS("predictions.rds")
stopifnot(nrow(pred) == 4, all(pred$score >= 0.5), !any(pred$event_type == "intron retention"))
stopifnot(pred$is_canonical[pred$junc_id == "chr1:160-301:+"],
          pred$is_in_normal[pred$junc_id == "chr1:160-298:+"],
          "chr1:966-1101:-" %in% pred$junc_id)
integrate_rna("predictions.rds", "reference.rds", "genome.fa",
              file.path(fixture, "selected.gz"), file.path(fixture, "evidence.tsv"),
              "P1", "day0", 13L, 8L)
peptides <- readr::read_tsv("peptide_candidates.tsv", show_col_types = FALSE)
stopifnot(nrow(peptides) == 2, all(peptides$peptide_exported),
          any(grepl("M{13}KK", peptides$peptide_context)),
          any(grepl("M{13}GG", peptides$peptide_context)))
stopifnot(sum(grepl("^>", readLines("peptide_contexts.fasta"))) == 2)
message("PASS: real splice2neo translation, both strands, canonical/normal exclusions, per-effect threshold")

# No observed RNA and no DNA predictions must be successful, truly empty output.
empty_evidence <- file.path(out, "empty.tsv")
write_tsv(readr::read_tsv(file.path(fixture, "evidence.tsv"), show_col_types = FALSE)[0, ], empty_evidence)
integrate_rna("predictions.rds", "reference.rds", "genome.fa", "unused", empty_evidence, "P1", "day99", 13L, 8L)
stopifnot(file.size("peptide_contexts.fasta") == 0, nrow(readr::read_tsv("peptide_candidates.tsv", show_col_types = FALSE)) == 0)
saveRDS(empty_predictions(), "empty.rds")
integrate_rna("empty.rds", "reference.rds", "genome.fa",
              file.path(fixture, "selected.gz"), file.path(fixture, "evidence.tsv"), "P1", "emptyDNA", 13L, 8L)
stopifnot(file.size("peptide_contexts.fasta") == 0)
message("PASS: empty RNA and empty DNA outputs")

# An exploratory opt-out must propagate its warning flag, not claim a normal filter.
ref$normal_filtered <- FALSE
ref$cds <- ref$cds[0]
saveRDS(ref, "no_cds.rds")
integrate_rna("predictions.rds", "no_cds.rds", "genome.fa",
              file.path(fixture, "selected.gz"), file.path(fixture, "evidence.tsv"), "P1", "noCDS", 13L, 8L)
peptides <- readr::read_tsv("peptide_candidates.tsv", show_col_types = FALSE)
stopifnot(nrow(peptides) == 2, !any(peptides$peptide_exported),
          !any(peptides$normal_filter_applied), file.size("peptide_contexts.fasta") == 0)
message("PASS: no coding transcript and normal-filter provenance")
cat("All splice2neo integration assertions passed.\n")
