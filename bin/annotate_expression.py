#!/usr/bin/env python3
"""
Add normal-tissue (GTEx) and tumour-cohort (TCGA) expression context to a
pVACseq aggregated report, plus this sample's own TPM when RNA was run.

pVACseq filters on the patient's own RNA at the variant site; it says
nothing about whether the gene is expressed in healthy tissue, nor whether
the tumour type expresses it at all. Both matter when picking vaccine
targets:

  * A candidate in a gene that healthy tissues express is a worse target:
    the same peptide's wild-type context is presented by normal cells, so
    T cells against it are more likely to be tolerised, and on-target
    off-tumour toxicity is a real risk. The GTEx columns say how widely the
    gene is expressed in normal tissue and how high it goes in the one
    tissue this tumour arose from.
  * A candidate whose gene the matched TCGA cohort expresses highly is more
    plausibly expressed in this patient than one the cohort never
    expresses. That is weak evidence next to the patient's own RNA, but it
    is the only expression evidence available when RNA is missing, so the
    report falls back to it.

Reference tables come from bin/build_expression_db.py (GTEx and TCGA are
downloaded by the user, nothing ships here). The sample TPM table comes
from bin/gene_tpm.py. Genes are matched on the report's "Gene" symbol, or
its Ensembl ID where the report carries one.

Columns added:

    Sample TPM, Sample Expression Quartile   this sample's RNA (NA without RNA)
    GTEx <tissue> TPM                        the matched normal tissue
    GTEx Median TPM, GTEx Max TPM            across all tissues in the table
    GTEx Tissues Expressed                   how many tissues are above --expressed-tpm
    GTEx Specificity                         not_expressed / testis_restricted /
                                             tissue_restricted / broadly_expressed
    TCGA <cancer type> TPM                   the matched tumour cohort
    Expression Evidence                      which value Expression Flag used
    Expression Flag                          OK / LOW_EXPRESSION / NO_EXPRESSION / NO_DATA

The flag is a heuristic sort key, not a verdict: thresholds are arguments,
and a flagged candidate is ranked down, never dropped. "NO_EXPRESSION"
means this sample's RNA saw the gene at zero TPM, which with shallow RNA
is also what a lowly expressed gene looks like.

Usage:
    annotate_expression.py tumor.MHC_I.aggregated.tsv \\
        --gtex gtex_v8_median_tpm.tsv --gtex-tissue "Skin - Sun Exposed (Lower leg)" \\
        --tcga tcga_cancertype_tpm.tsv --tcga-cancer-type "TCGA-SKCM" \\
        [--sample-tpm patient_tpm.tsv] [--min-tpm 1.0] -o out.tsv
"""
import argparse
import csv
import math
from pathlib import Path

TESTIS_KEYS = ("testis",)


def parse_float(raw):
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(value) else value


def fmt(value, digits=3):
    return "NA" if value is None else f"{value:.{digits}f}"


def read_expression_db(path, quantile=0.5):
    """build_expression_db.py table -> {key: {group: tpm}} for one quantile.

    Indexed by both Ensembl ID and gene symbol, so a report carrying either
    one resolves. A symbol used by more than one Ensembl ID (a handful in
    GENCODE) keeps the entry with the higher value in the group, which keeps
    the annotation conservative: it never understates normal expression.
    """
    by_key = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"ENSEMBL_ID", "GENE_NAME", "GROUP", "QUANTILE", "TPM"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"{path}: missing column(s) {sorted(missing)}: rebuild it with bin/build_expression_db.py")
        for row in reader:
            if parse_float(row["QUANTILE"]) != quantile:
                continue
            tpm = parse_float(row["TPM"])
            if tpm is None:
                continue
            group = row["GROUP"]
            for key in (row["ENSEMBL_ID"], row["GENE_NAME"]):
                if not key or key == "NA":
                    continue
                groups = by_key.setdefault(key, {})
                if group not in groups or tpm > groups[group]:
                    groups[group] = tpm
    if not by_key:
        raise SystemExit(f"{path}: no rows at quantile {quantile:g}")
    return by_key


def read_sample_tpm(path):
    """gene_tpm.py table -> {key: (tpm, quartile)}, indexed by ID and symbol."""
    by_key = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            tpm = parse_float(row.get("TPM"))
            if tpm is None:
                continue
            entry = (tpm, row.get("EXPRESSION_QUARTILE", "NA"))
            for key in (row.get("ENSEMBL_ID"), row.get("GENE_NAME")):
                if not key or key == "NA":
                    continue
                if key not in by_key or tpm > by_key[key][0]:
                    by_key[key] = entry
    return by_key


def gene_keys(row):
    """Candidate lookup keys from a report row, in order of preference.

    pVACseq's aggregated report keys candidates by gene symbol; a few rows
    carry several symbols for one variant (overlapping genes), so each is
    tried in turn.
    """
    keys = []
    for column in ("Gene", "Ensembl Gene ID", "Gene Name", "gene"):
        raw = (row.get(column) or "").strip()
        if not raw or raw == "NA":
            continue
        for part in raw.replace(";", ",").split(","):
            part = part.strip()
            if part and part not in keys:
                keys.append(part)
    return keys


def lookup(db, keys, group):
    for key in keys:
        groups = db.get(key)
        if groups and group in groups:
            return groups[group]
    return None


def gtex_summary(db, keys, expressed_tpm):
    """(tissue tpm lookup table, median, max, n expressed, specificity call)."""
    for key in keys:
        groups = db.get(key)
        if groups:
            values = sorted(groups.values())
            middle = len(values) // 2
            median = values[middle] if len(values) % 2 else (values[middle - 1] + values[middle]) / 2
            expressed = [tissue for tissue, tpm in groups.items() if tpm > expressed_tpm]
            if not expressed:
                specificity = "not_expressed"
            elif all(any(token in tissue.lower() for token in TESTIS_KEYS) for tissue in expressed):
                specificity = "testis_restricted"
            elif len(expressed) <= max(1, round(0.1 * len(groups))):
                specificity = "tissue_restricted"
            else:
                specificity = "broadly_expressed"
            return groups, median, max(values), len(expressed), specificity
    return None, None, None, None, "NA"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("report", type=Path, help="pVACseq aggregated report (or the CCF-annotated copy of it)")
    parser.add_argument("--gtex", type=Path, help="GTEx table from build_expression_db.py")
    parser.add_argument("--gtex-tissue", help="GROUP to report, e.g. 'Skin - Sun Exposed (Lower leg)'")
    parser.add_argument("--tcga", type=Path, help="TCGA table from build_expression_db.py")
    parser.add_argument("--tcga-cancer-type", help="GROUP to report, e.g. 'TCGA-SKCM'")
    parser.add_argument("--sample-tpm", type=Path, help="this sample's TPM table from gene_tpm.py")
    parser.add_argument("--quantile", type=float, default=0.5, help="quantile to read from the tables (default: 0.5)")
    parser.add_argument("--min-tpm", type=float, default=1.0,
                        help="at or below this TPM a candidate is flagged LOW_EXPRESSION (default: 1.0)")
    parser.add_argument("--expressed-tpm", type=float, default=1.0,
                        help="TPM above which a GTEx tissue counts as expressing the gene (default: 1.0)")
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args(argv)

    if args.gtex and not args.gtex_tissue:
        parser.error("--gtex needs --gtex-tissue")
    if args.tcga and not args.tcga_cancer_type:
        parser.error("--tcga needs --tcga-cancer-type")

    gtex = read_expression_db(args.gtex, args.quantile) if args.gtex else {}
    tcga = read_expression_db(args.tcga, args.quantile) if args.tcga else {}
    sample = read_sample_tpm(args.sample_tpm) if args.sample_tpm else {}

    for name, db, group in (("--gtex-tissue", gtex, args.gtex_tissue), ("--tcga-cancer-type", tcga, args.tcga_cancer_type)):
        if db and group not in {g for groups in db.values() for g in groups}:
            available = sorted({g for groups in db.values() for g in groups})
            raise SystemExit(f"{name} '{group}' is not in the table. Available: " + ", ".join(available))

    new_columns = ["Sample TPM", "Sample Expression Quartile"]
    if gtex:
        new_columns += [f"GTEx {args.gtex_tissue} TPM", "GTEx Median TPM", "GTEx Max TPM",
                        "GTEx Tissues Expressed", "GTEx Specificity"]
    if tcga:
        new_columns += [f"TCGA {args.tcga_cancer_type} TPM"]
    new_columns += ["Expression Evidence", "Expression Flag"]

    with args.report.open(newline="") as src, args.output.open("w", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit(f"{args.report}: empty report")
        writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
        writer.writerow(reader.fieldnames + new_columns)

        for row in reader:
            keys = gene_keys(row)
            sample_tpm, sample_quartile = sample.get(next((k for k in keys if k in sample), ""), (None, "NA"))
            extra = [fmt(sample_tpm), sample_quartile]

            if gtex:
                groups, median, maximum, n_expressed, specificity = gtex_summary(gtex, keys, args.expressed_tpm)
                tissue_tpm = groups.get(args.gtex_tissue) if groups else None
                extra += [fmt(tissue_tpm), fmt(median), fmt(maximum),
                          "NA" if n_expressed is None else str(n_expressed), specificity]

            tcga_tpm = lookup(tcga, keys, args.tcga_cancer_type) if tcga else None
            if tcga:
                extra += [fmt(tcga_tpm)]

            # The patient's own RNA wins; the cohort stands in only without it.
            if sample:
                evidence, value = "sample_rna", sample_tpm
            elif tcga:
                evidence, value = "tcga_cohort", tcga_tpm
            else:
                evidence, value = "none", None

            if value is None:
                flag = "NO_DATA"
            elif value <= 0:
                flag = "NO_EXPRESSION"
            elif value <= args.min_tpm:
                flag = "LOW_EXPRESSION"
            else:
                flag = "OK"
            extra += [evidence, flag]

            writer.writerow([row.get(field, "") for field in reader.fieldnames] + extra)


if __name__ == "__main__":
    main()
