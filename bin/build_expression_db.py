#!/usr/bin/env python3
"""
Build the GTEx / TCGA reference expression tables that
bin/annotate_expression.py reads.

Both references answer the same question for one gene: what is its normal
expression in healthy tissue (GTEx), and its typical expression in this
tumour type (TCGA)? A neoantigen in a gene that every healthy tissue
expresses is a worse vaccine target than one in a gene that is silent
outside the tumour, and a gene that the matching TCGA cohort expresses
highly is more plausibly really expressed in this patient than one the
cohort never expresses - which matters most when there is no RNA.

Everything here is computed from public data you download yourself; no
reference data ships with this repo. Two input shapes are supported:

  --from-gct     A GTEx median-TPM GCT: genes x tissues, one median TPM per
                 cell. This is the small public summary file
                 (GTEx_Analysis_*_gene_median_tpm.gct.gz, ~15 MB) and is
                 enough for the annotation columns, but it has medians only,
                 so --quantiles is ignored.

  --from-matrix  A gene x sample TPM matrix plus --groups, a two-column
                 sample-to-group table (sample<TAB>group). Per-group
                 quantiles are computed across that group's samples. Use
                 this for TCGA (e.g. a UCSC Xena / GDC TPM matrix with its
                 phenotype table mapping each barcode to a project) and for
                 GTEx when you want percentiles rather than medians.

Output (TSV, one row per gene x group x quantile), sorted by gene:

    ENSEMBL_ID  GENE_NAME  GROUP  QUANTILE  TPM

Values are TPM as given by the source; nothing is re-normalised, so the
matrix you feed in must already be TPM (not counts, not log2(TPM+1) -
pass --expm1 for a log2(x+1) matrix such as Xena's).

A cohort matrix usually labels genes by Ensembl ID only, while pVACseq
reports gene symbols, so pass --gene-names with the GTF used elsewhere in
the pipeline to fill the GENE_NAME column; the annotator can then match
candidates either way.

Usage:
    build_expression_db.py --from-gct GTEx_..._gene_median_tpm.gct.gz -o gtex_v8_median_tpm.tsv
    build_expression_db.py --from-matrix tcga_tpm.tsv --groups tcga_sample_project.tsv \\
        --quantiles 0.25 0.5 0.75 --expm1 -o tcga_cancertype_tpm.tsv
"""
import argparse
import csv
import gzip
import math
from pathlib import Path


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def strip_version(gene_id):
    return gene_id.split(".", 1)[0] if gene_id.startswith("ENSG") else gene_id


def quantile(sorted_values, fraction):
    """Nearest-rank quantile of a sorted, non-empty list."""
    if not sorted_values:
        return None
    index = max(0, min(len(sorted_values) - 1, int(math.ceil(fraction * len(sorted_values))) - 1))
    return sorted_values[index]


def read_gene_names(path):
    """GTF -> {gene_id (unversioned): gene_name}, from gene rows."""
    import re
    attribute_re = re.compile(r'(\S+) "([^"]*)"')
    names = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attributes = dict(attribute_re.findall(fields[8]))
            if attributes.get("gene_id") and attributes.get("gene_name"):
                names[strip_version(attributes["gene_id"])] = attributes["gene_name"]
    if not names:
        raise SystemExit(f"{path}: no gene rows with a gene_name attribute")
    return names


def read_groups(path):
    """sample<TAB>group table -> {sample: group}, header row tolerated."""
    mapping = {}
    with open_text(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or not fields[0] or not fields[1]:
                continue
            mapping[fields[0].strip()] = fields[1].strip()
    if not mapping:
        raise SystemExit(f"{path}: no sample-to-group rows found")
    # a header like "sample<TAB>project" maps a sample called "sample": harmless
    return mapping


def rows_from_gct(path):
    """GTEx median-TPM GCT -> (gene_id, gene_name, {tissue: [(0.5, tpm)]})."""
    with open_text(path) as handle:
        version = handle.readline()
        if not version.startswith("#1.2"):
            raise SystemExit(f"{path}: not a GCT file (first line is {version.strip()!r}, expected '#1.2')")
        handle.readline()  # dimensions
        header = handle.readline().rstrip("\n").split("\t")
        if header[:2] != ["Name", "Description"]:
            raise SystemExit(f"{path}: unexpected GCT header {header[:2]}")
        tissues = header[2:]
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(header):
                continue
            values = {}
            for tissue, raw in zip(tissues, fields[2:]):
                if raw in ("", "NA"):
                    continue
                values[tissue] = [(0.5, float(raw))]
            yield strip_version(fields[0]), fields[1] or "NA", values


def rows_from_matrix(path, groups, quantiles, expm1):
    """Gene x sample TPM matrix -> (gene_id, gene_name, {group: [(q, tpm)]}).

    Streams one gene at a time, so a full cohort matrix never has to fit in
    memory. Samples missing from --groups are skipped, which is what lets a
    phenotype table that covers several cohorts select just one of them.
    """
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        sample_groups = [groups.get(name.strip()) for name in header[1:]]
        if not any(sample_groups):
            raise SystemExit(f"{path}: none of the matrix's column names appear in the --groups table")
        for fields in reader:
            if not fields or not fields[0]:
                continue
            per_group = {}
            for group, raw in zip(sample_groups, fields[1:]):
                if group is None or raw in ("", "NA", "NaN"):
                    continue
                value = float(raw)
                if expm1:
                    value = (2.0 ** value) - 1.0
                per_group.setdefault(group, []).append(value)
            values = {}
            for group, group_values in per_group.items():
                group_values.sort()
                values[group] = [(q, quantile(group_values, q)) for q in quantiles]
            yield strip_version(fields[0]), "NA", values


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--from-gct", type=Path, help="GTEx median-TPM GCT")
    source.add_argument("--from-matrix", type=Path, help="gene x sample TPM matrix")
    parser.add_argument("--groups", type=Path, help="sample<TAB>group table (required with --from-matrix)")
    parser.add_argument("--quantiles", type=float, nargs="+", default=[0.5],
                        help="quantiles to store per group (default: 0.5)")
    parser.add_argument("--gene-names", type=Path,
                        help="GTF to take gene symbols from, for sources that label genes by Ensembl ID only")
    parser.add_argument("--expm1", action="store_true",
                        help="matrix holds log2(TPM+1) (e.g. UCSC Xena): convert back to TPM")
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args(argv)

    if args.from_matrix and not args.groups:
        parser.error("--from-matrix needs --groups")
    if any(not 0 < q <= 1 for q in args.quantiles):
        parser.error("--quantiles must be in (0, 1]")

    if args.from_gct:
        rows = rows_from_gct(args.from_gct)
    else:
        rows = rows_from_matrix(args.from_matrix, read_groups(args.groups), sorted(args.quantiles), args.expm1)

    gene_names = read_gene_names(args.gene_names) if args.gene_names else {}

    written = 0
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["ENSEMBL_ID", "GENE_NAME", "GROUP", "QUANTILE", "TPM"])
        for gene_id, gene_name, values in rows:
            gene_name = gene_names.get(gene_id, gene_name)
            for group in sorted(values):
                for q, tpm in values[group]:
                    if tpm is None:
                        continue
                    writer.writerow([gene_id, gene_name, group, f"{q:g}", f"{tpm:.4f}"])
                    written += 1
    if not written:
        raise SystemExit("no rows written: check the input really is a TPM table")
    print(f"wrote {written} rows to {args.output}")


if __name__ == "__main__":
    main()
