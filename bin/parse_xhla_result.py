#!/usr/bin/env python3
"""
Parse xHLA's report-<sample>-hla.json into a pVACseq-formatted allele
string covering HLA class I (A/B/C) and class II DRB1.

Format confirmed from griffithlab/pVACtools discussion #1124: class I
alleles need an "HLA-" prefix (e.g. "HLA-A*02:01"); class II alleles do
NOT ("DRB1*15:01", not "HLA-DRB1*15:01"). DRB1 needs no alpha-chain
pairing (DRA is monomorphic) so a bare DRB1 call is usable as-is.

xHLA also types DQB1 and DPB1, but only the beta chain. pVACtools needs
DQ/DP alleles as an alpha-beta PAIR joined with a hyphen (e.g.
"DQA1*01:02-DQB1*06:02"): xHLA does not type DQA1/DPA1 at all, so there
is no correct way to build that pair from its output alone. Rather than
guess at a population-common alpha allele (which could silently produce
wrong binding predictions), those beta-only calls are written to a
separate file and excluded from the string handed to pVACseq.

xHLA reports each homozygous locus twice (e.g. "C*06:02","C*06:02") by
design: deduplicated here since pVACseq only needs each allele once.

Usage:
    parse_xhla_result.py report.json -o pvacseq_alleles.txt [--excluded-out excluded.txt]
"""
import argparse
import json
import sys

CLASS_I_LOCI = {"A", "B", "C"}
USABLE_CLASS_II_LOCI = {"DRB1"}               # DRA is monomorphic, no pairing needed
INCOMPLETE_CLASS_II_LOCI = {"DQB1", "DPB1"}   # beta chain only: needs an alpha chain xHLA doesn't type


def locus_of(allele):
    return allele.split("*")[0]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("result_json")
    p.add_argument("-o", "--out", required=True, help="pvacseq-ready comma-joined allele string")
    p.add_argument("--excluded-out", default=None, help="alleles typed but not usable as-is (DQB1/DPB1 with no alpha chain): only written if non-empty")
    args = p.parse_args()

    with open(args.result_json) as fh:
        data = json.load(fh)

    raw_alleles = data.get("hla", {}).get("alleles", [])
    if not raw_alleles:
        sys.exit(f"[ERROR] No alleles found in {args.result_json}")

    usable = set()
    excluded = set()
    for allele in raw_alleles:
        locus = locus_of(allele)
        if locus in CLASS_I_LOCI:
            usable.add(f"HLA-{allele}")
        elif locus in USABLE_CLASS_II_LOCI:
            usable.add(allele)
        elif locus in INCOMPLETE_CLASS_II_LOCI:
            excluded.add(allele)
        else:
            print(f"[WARN] Unrecognised locus in xHLA output: {allele}", file=sys.stderr)

    with open(args.out, "w") as fh:
        fh.write(",".join(sorted(usable)) + "\n")

    if excluded:
        print(f"[WARN] {len(excluded)} DQB1/DPB1 allele(s) typed but excluded from the "
              f"pVACseq allele string (no alpha chain: see this script's docstring): "
              f"{sorted(excluded)}", file=sys.stderr)
        if args.excluded_out:
            with open(args.excluded_out, "w") as fh:
                for a in sorted(excluded):
                    fh.write(a + "\n")


if __name__ == "__main__":
    main()
