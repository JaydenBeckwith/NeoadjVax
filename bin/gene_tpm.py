#!/usr/bin/env python3
"""
Turn a featureCounts gene-level count table into TPM.

The pipeline's RNA BAMs are already aligned to GRCh38, so gene expression
only needs counting and normalising, not re-quantification from FASTQ.
featureCounts gives raw counts plus each gene's union-exon length, which is
all TPM needs:

    RPK_gene = counts / (exon_length_bp / 1000)
    TPM_gene = RPK_gene / (sum of all RPK / 1e6)

TPM (rather than FPKM/CPM) because the GTEx and TCGA reference tables this
feeds into are themselves TPM, and TPM is the unit that is comparable
between samples without a library-size correction term.

Each gene also gets an expression quartile, computed over the genes with
non-zero TPM in this sample (quartile 1 = lowest expressed quarter of the
expressed genes, 4 = highest; a gene with TPM 0 is reported as 0). It is a
within-sample rank, so it says "high for this tumour", not "high in
absolute terms" - useful for sorting candidates, not for cross-sample
comparison.

Gene names and biotypes are read from the GTF when one is given; without it
the output carries Ensembl IDs only and the GENE_NAME column repeats the ID.

Usage:
    gene_tpm.py counts.txt [--gtf gencode.v44.annotation.gtf.gz] -o tpm.tsv
"""
import argparse
import csv
import gzip
import re
from pathlib import Path

COLUMNS = ["ENSEMBL_ID", "GENE_NAME", "GENE_BIOTYPE", "COUNTS", "LENGTH", "TPM", "EXPRESSION_QUARTILE"]

ATTRIBUTE_RE = re.compile(r'(\S+) "([^"]*)"')


def strip_version(gene_id):
    """ENSG00000141510.16 -> ENSG00000141510 (GTEx/TCGA tables are unversioned)."""
    return gene_id.split(".", 1)[0] if gene_id.startswith("ENSG") else gene_id


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_gtf_genes(path):
    """GTF -> {gene_id (unversioned): (gene_name, gene_biotype)} from gene rows."""
    genes = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attributes = dict(ATTRIBUTE_RE.findall(fields[8]))
            gene_id = attributes.get("gene_id")
            if not gene_id:
                continue
            biotype = attributes.get("gene_type") or attributes.get("gene_biotype") or "NA"
            genes[strip_version(gene_id)] = (attributes.get("gene_name") or "NA", biotype)
    return genes


def read_featurecounts(path):
    """featureCounts output -> [(gene_id, length_bp, counts)].

    Column layout is Geneid, Chr, Start, End, Strand, Length, <one column per
    BAM>; the comment line featureCounts writes first is skipped. Counts from
    several BAM columns are summed, so passing more than one BAM of the same
    sample (e.g. lanes counted together) still gives one expression value.
    """
    rows = []
    with open_text(path) as handle:
        reader = csv.reader((l for l in handle if not l.startswith("#")), delimiter="\t")
        header = next(reader)
        if header[:1] != ["Geneid"]:
            raise SystemExit(f"{path}: not a featureCounts table (first column is '{header[:1]}', expected 'Geneid')")
        length_index = header.index("Length")
        for row in reader:
            if not row or not row[0]:
                continue
            counts = sum(float(value) for value in row[length_index + 1:] if value not in ("", "NA"))
            rows.append((row[0], float(row[length_index]), counts))
    if not rows:
        raise SystemExit(f"{path}: no gene rows")
    return rows


def expression_quartiles(tpms):
    """TPM values -> {tpm: quartile}, quartiles taken over non-zero TPMs only.

    Zero-TPM genes are quartile 0 rather than being lumped into the bottom
    quarter: in a typical tumour RNA sample most annotated genes are not
    expressed at all, and letting them set the quartile boundaries would push
    every real expression value into the top quartile.
    """
    expressed = sorted(value for value in tpms if value > 0)
    if not expressed:
        return {}
    bounds = [expressed[min(len(expressed) - 1, int(len(expressed) * fraction))] for fraction in (0.25, 0.5, 0.75)]
    cuts = {}
    for value in set(tpms):
        if value <= 0:
            cuts[value] = 0
            continue
        quartile = 1
        for bound in bounds:
            if value > bound:
                quartile += 1
        cuts[value] = quartile
    return cuts


def quantify(count_rows, gene_info):
    """[(gene_id, length, counts)] -> rows of COLUMNS, highest TPM first."""
    rpks = []
    for gene_id, length, counts in count_rows:
        rpks.append(counts / (length / 1000.0) if length > 0 else 0.0)
    rpk_per_million = sum(rpks) / 1e6
    if rpk_per_million <= 0:
        raise SystemExit("no reads assigned to any gene: check the GTF matches the BAM's reference and strandedness")

    tpms = [rpk / rpk_per_million for rpk in rpks]
    quartiles = expression_quartiles(tpms)

    rows = []
    for (gene_id, length, counts), tpm in zip(count_rows, tpms):
        unversioned = strip_version(gene_id)
        name, biotype = gene_info.get(unversioned, (unversioned, "NA"))
        rows.append([unversioned, name, biotype, f"{counts:.0f}", f"{length:.0f}",
                     f"{tpm:.4f}", str(quartiles.get(tpm, 0))])
    rows.sort(key=lambda row: float(row[5]), reverse=True)
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("counts", type=Path, help="featureCounts gene-level output")
    parser.add_argument("--gtf", type=Path, help="GTF used for counting, for gene names and biotypes")
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args(argv)

    gene_info = read_gtf_genes(args.gtf) if args.gtf else {}
    rows = quantify(read_featurecounts(args.counts), gene_info)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        writer.writerows(rows)


if __name__ == "__main__":
    main()
