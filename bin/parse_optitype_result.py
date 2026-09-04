#!/usr/bin/env python3
"""
Parse OptiType's <sample>_result.tsv into a pVACseq-formatted HLA class I
allele string, e.g.: HLA-A*02:01,HLA-A*24:02,HLA-B*07:02,HLA-B*15:01,
HLA-C*03:04,HLA-C*07:02

OptiType's result columns are: A1, A2, B1, B2, C1, C2, Reads, Objective
(no "HLA-" prefix on the allele values) — this just reformats them.
"""
import argparse
import csv
import sys


def main():
    p = argparse.ArgumentParser()
    p.add_argument("result_tsv")
    p.add_argument("-o", "--out", default=None, help="write to file instead of stdout")
    args = p.parse_args()

    with open(args.result_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        row = next(reader, None)
    if row is None:
        sys.exit(f"[ERROR] No rows found in {args.result_tsv}")

    alleles = []
    for col in ("A1", "A2", "B1", "B2", "C1", "C2"):
        val = row.get(col, "").strip()
        if val:
            alleles.append(f"HLA-{val}")

    out_str = ",".join(alleles)
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(out_str + "\n")
    else:
        print(out_str)


if __name__ == "__main__":
    main()
