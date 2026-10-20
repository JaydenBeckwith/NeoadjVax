#!/usr/bin/env python3
"""
Add one RNA timepoint's read counts to a pVACseq aggregated report.

pVACseq reads RNA counts for one RNA sample only (RAD/RDP/RAF in the tumour
column of its input VCF), so the pipeline gives it the baseline RNA. A
later timepoint (e.g. week 6 on treatment) is counted the same way into a
separate copy of the VCF; this joins that VCF's tumour RDP and RAF onto the
report as "RNA Depth (<label>)" and "RNA VAF (<label>)".

Report rows are matched by their pVACseq ID, "chrom-start-stop-ref-alt",
where ref and alt are the VCF REF and ALT and start is the VCF POS or one
less (pVACseq's start is vcfpy's affected_start: POS - 1 for an SNV, POS
for an insertion; for a deletion it adds 1 to POS - 1, giving POS). A variant
missing from the VCF, or without counts there, gets NA.

Usage:
    add_rna_timepoint.py tumor.Combined.all_epitopes.aggregated.tsv \\
        --vcf P1.rna_week6.readcounts.vcf --label "week 6" -o out.tsv
"""
import argparse
import csv
import gzip
from pathlib import Path


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_rna_counts(path, sample="tumor"):
    """{(chrom, ref, alt): [(pos, depth, vaf)]} from the sample's RDP/RAF."""
    counts = {}
    column = None
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                if sample not in fields:
                    raise ValueError(f"{path}: no '{sample}' sample column")
                column = fields.index(sample)
                continue
            chrom, pos, _, ref, alts = fields[:5]
            values = dict(zip(fields[8].split(":"), fields[column].split(":")))
            depth = values.get("RDP", ".")
            vafs = values.get("RAF", ".").split(",")
            for i, alt in enumerate(alts.split(",")):
                vaf = vafs[i] if i < len(vafs) else "."
                counts.setdefault((chrom, ref, alt), []).append((int(pos), depth, vaf))
    return counts


def lookup(counts, variant_id):
    """(depth, vaf) for a pVACseq ID, or (NA, NA)."""
    chrom, start, _stop, ref, alt = variant_id.split("-", 4)
    for pos, depth, vaf in counts.get((chrom, ref, alt), []):
        if int(start) in (pos - 1, pos):
            return ("NA" if depth == "." else depth), ("NA" if vaf == "." else vaf)
    return "NA", "NA"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("report", type=Path)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args(argv)

    counts = read_rna_counts(args.vcf)
    with args.report.open(newline="") as src, args.output.open("w", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
        writer.writerow(reader.fieldnames + [f"RNA Depth ({args.label})", f"RNA VAF ({args.label})"])
        for row in reader:
            writer.writerow([row[f] for f in reader.fieldnames] + list(lookup(counts, row["ID"])))


if __name__ == "__main__":
    main()
