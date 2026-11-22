#!/usr/bin/env python3
"""
Shared tumour antigens as a second channel beside the patient's own
mutations.

A low-mutation-burden patient can run out of neoantigens. Cancer-testis
antigens, melanocyte differentiation antigens and the broadly
overexpressed genes are not patient-specific and so are not neoantigens,
but they are real targets with characterised epitopes, and whether they
are usable in a given patient is an expression question this pipeline can
already answer.

For each gene in the panel this reports the patient's own tumour TPM, how
much the matched normal tissue expresses it (GTEx), how much the cohort
does (TCGA), and a call:

    good_target        expressed in this tumour and not in normal tissue
                       outside testis
    target_with_risk   expressed in this tumour, but normal tissue
                       expresses it too, so an on-target off-tumour
                       response is possible
    not_expressed      not expressed in this tumour; nothing to target
    no_data            no RNA for this patient and no reference table

The call is about expression only. It says nothing about whether the
patient's alleles can present anything from the gene: that needs the
epitope set, which is the --taa_peptides half of this channel (pVACbind).

The reference-table readers are shared with gene_panel_report.py, so GTEx
and TCGA are read exactly as they are everywhere else in the pipeline.

Usage:
    taa_report.py --panel assets/taa_panel.tsv -o P1.taa.tsv \\
        [--sample-tpm P1.tpm.tsv] [--gtex gtex.tsv --gtex-tissue 'Skin - Sun Exposed (Lower leg)'] \\
        [--tcga tcga.tsv --tcga-cancer-type TCGA-SKCM]
"""
import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gene_panel_report as panel_report  # noqa: E402  (same directory, on PATH in the pipeline)

COLUMNS = ["Class", "Gene", "Note", "Sample TPM", "Sample Expressed",
           "GTEx Tissue TPM", "GTEx Max TPM", "GTEx Max Tissue", "GTEx Specificity",
           "TCGA Cohort TPM", "Evidence", "Call"]

# Testis is immune-privileged, so expression there does not make a target
# risky the way expression in any other tissue does.
PRIVILEGED = ("testis",)


def read_panel(path):
    """gene<TAB>class[<TAB>note] -> [(gene, class, note)], comments skipped."""
    rows = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = [f.strip() for f in line.split("\t")]
            if not fields[0]:
                continue
            rows.append((fields[0], fields[1] if len(fields) > 1 else "", fields[2] if len(fields) > 2 else ""))
    if not rows:
        raise SystemExit(f"{path}: no gene rows found")
    return rows


def group_values(db, gene, quantile):
    """{group: tpm at this quantile} for one gene, across every group."""
    entry = db.get(gene) or {}
    values = {}
    for group, quantiles in entry.items():
        if quantile in quantiles:
            values[group] = quantiles[quantile]
        elif quantiles:
            nearest = min(quantiles, key=lambda q: abs(q - quantile))
            values[group] = quantiles[nearest]
    return values


def specificity(gtex_values, expressed_tpm):
    """How widely normal tissue expresses the gene."""
    if not gtex_values:
        return "NA"
    expressing = [tissue for tissue, tpm in gtex_values.items() if tpm > expressed_tpm]
    if not expressing:
        return "not_expressed"
    non_privileged = [t for t in expressing
                      if not any(p in t.lower() for p in PRIVILEGED)]
    if not non_privileged:
        return "testis_restricted"
    if len(non_privileged) <= 3:
        return "tissue_restricted"
    return "broadly_expressed"


def call_target(sample_tpm, cohort_tpm, spec, min_tpm):
    """-> (evidence, call). The patient's own RNA wins where there is any."""
    if sample_tpm is not None:
        expressed, evidence = sample_tpm > min_tpm, "sample_rna"
    elif cohort_tpm is not None:
        expressed, evidence = cohort_tpm > min_tpm, "tcga_cohort"
    else:
        return "none", "no_data"
    if not expressed:
        return evidence, "not_expressed"
    # Expressed in the tumour. Whether that is clean depends on normal
    # tissue; with no GTEx table at all (spec "NA") the risk is unknown,
    # which is reported as risk rather than as a clean target.
    if spec in ("not_expressed", "testis_restricted"):
        return evidence, "good_target"
    return evidence, "target_with_risk"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--panel", type=Path, required=True, help="gene<TAB>class[<TAB>note] table")
    parser.add_argument("--sample-tpm", type=Path, help="this patient's TPM table from bin/gene_tpm.py")
    parser.add_argument("--gtex", type=Path, help="GTEx table from bin/build_expression_db.py")
    parser.add_argument("--gtex-tissue", help="the matched normal tissue GROUP to report")
    parser.add_argument("--tcga", type=Path, help="TCGA table from bin/build_expression_db.py")
    parser.add_argument("--tcga-cancer-type", help="the matched cohort GROUP to report")
    parser.add_argument("--quantile", type=float, default=0.5, help="quantile to read from the tables (default 0.5)")
    parser.add_argument("--min-tpm", type=float, default=1.0, help="above this TPM a gene counts as expressed (default 1.0)")
    parser.add_argument("--expressed-tpm", type=float, default=1.0,
                        help="above this TPM a normal tissue counts as expressing the gene (default 1.0)")
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.gtex and not args.gtex_tissue:
        parser.error("--gtex needs --gtex-tissue")
    if args.tcga and not args.tcga_cancer_type:
        parser.error("--tcga needs --tcga-cancer-type")

    panel = read_panel(args.panel)
    sample = panel_report.read_sample_tpm(args.sample_tpm) if args.sample_tpm else {}
    gtex = panel_report.read_expression_db(args.gtex) if args.gtex else {}
    tcga = panel_report.read_expression_db(args.tcga) if args.tcga else {}

    for table, group, flag in ((gtex, args.gtex_tissue, "--gtex-tissue"), (tcga, args.tcga_cancer_type, "--tcga-cancer-type")):
        if table and group:
            groups = {g for entry in table.values() for g in entry}
            if group not in groups:
                raise SystemExit(f"{flag}: '{group}' is not in the table. Available: {', '.join(sorted(groups)[:20])}")

    fmt = panel_report.fmt
    rows = []
    for gene, gene_class, note in panel:
        sample_tpm = sample.get(gene)
        gtex_values = group_values(gtex, gene, args.quantile)
        tcga_values = group_values(tcga, gene, args.quantile)
        tissue_tpm = gtex_values.get(args.gtex_tissue) if args.gtex_tissue else None
        cohort_tpm = tcga_values.get(args.tcga_cancer_type) if args.tcga_cancer_type else None
        max_tissue, max_tpm = ("NA", None)
        if gtex_values:
            max_tissue = max(gtex_values, key=lambda t: gtex_values[t])
            max_tpm = gtex_values[max_tissue]
        spec = specificity(gtex_values, args.expressed_tpm)
        evidence, call = call_target(sample_tpm, cohort_tpm, spec, args.min_tpm)
        rows.append([
            gene_class, gene, note, fmt(sample_tpm),
            "NA" if sample_tpm is None else ("yes" if sample_tpm > args.min_tpm else "no"),
            fmt(tissue_tpm), fmt(max_tpm), max_tissue, spec, fmt(cohort_tpm), evidence, call,
        ])

    # Usable targets first, then risky ones, then the rest; within a group
    # the highest expressed in this tumour.
    order = {"good_target": 0, "target_with_risk": 1, "not_expressed": 2, "no_data": 3}
    rows.sort(key=lambda r: (order.get(r[-1], 9), -(panel_report.parse_float(r[3]) or 0.0), r[1]))

    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        writer.writerows(rows)

    counts = {}
    for row in rows:
        counts[row[-1]] = counts.get(row[-1], 0) + 1
    print(f"{len(rows)} shared antigens: " + ", ".join(f"{v} {k}" for k, v in sorted(counts.items())), file=sys.stderr)


if __name__ == "__main__":
    main()
