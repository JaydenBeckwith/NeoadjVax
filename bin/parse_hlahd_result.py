#!/usr/bin/env python3
"""
Turn an HLA-HD <ID>_final.result.txt into the allele list pVACseq takes.

HLA-HD writes one line per gene: the gene name, then two alleles at up to
four-field resolution, e.g.

    A       HLA-A*24:02:01:01   HLA-A*02:07:01
    DQA1    HLA-DQA1*01:02:01   -
    DRB3    Not typed

A "-" second allele means HLA-HD reported only one (homozygous, or the
second could not be resolved); "Not typed" means no call for that gene.

pVACseq wants two-field names: class I as "HLA-A*24:02"; class II DR beta
chains as "DRB1*15:01" (DRA is effectively monomorphic, so pVACseq pairs
it implicitly); DQ and DP as explicit alpha-beta heterodimers,
"DQA1*01:02-DQB1*06:02". Phase between the alpha and beta loci is
unknown, so every alpha x beta combination is emitted (cis and trans),
the usual practice when typing from short reads.

Null alleles (expression suffix N) are dropped: no protein, no
presentation. Genes outside the set below (HLA-E/F/G, DRB2/6-9 pseudogenes,
non-classical class II) are ignored.

Usage:
    parse_hlahd_result.py <ID>_final.result.txt \\
        --pvacseq-out alleles.txt --table-out alleles.tsv
"""
import argparse
import re
import sys
from pathlib import Path

CLASS_I = ("A", "B", "C")
DR_BETA = ("DRB1", "DRB3", "DRB4", "DRB5")
HETERODIMERS = (("DQA1", "DQB1"), ("DPA1", "DPB1"))
ALLELE_RE = re.compile(r"^(?:HLA-)?([A-Z]+[0-9]*)\*(\d+):(\d+)(?::\d+)*([A-Z]?)$")


def two_field(raw):
    """'HLA-A*24:02:01:01' -> ('A', 'A*24:02'); None for null/blank/untyped."""
    raw = raw.strip()
    if not raw or raw == "-" or raw.lower().startswith("not"):
        return None
    match = ALLELE_RE.match(raw)
    if not match:
        raise ValueError(f"unrecognised HLA-HD allele: {raw!r}")
    gene, field1, field2, suffix = match.groups()
    if suffix == "N":
        return None
    return gene, f"{gene}*{field1}:{field2}"


def parse_result(text):
    """gene -> list of distinct two-field alleles, in HLA-HD's order."""
    calls = {}
    for line in text.splitlines():
        cols = line.split("\t")
        if len(cols) < 2 or not cols[0].strip():
            continue
        gene = cols[0].strip()
        alleles = []
        for raw in cols[1:]:
            parsed = two_field(raw)
            if parsed is None:
                continue
            allele_gene, allele = parsed
            if allele_gene != gene:
                raise ValueError(f"allele {raw!r} listed under gene {gene}")
            if allele not in alleles:
                alleles.append(allele)
        calls[gene] = alleles
    return calls


def pvacseq_alleles(calls):
    out = []
    for gene in CLASS_I:
        out += [f"HLA-{a}" for a in calls.get(gene, [])]
    for gene in DR_BETA:
        out += calls.get(gene, [])
    for alpha_gene, beta_gene in HETERODIMERS:
        for alpha in calls.get(alpha_gene, []):
            for beta in calls.get(beta_gene, []):
                out.append(f"{alpha}-{beta}")
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("result", type=Path)
    parser.add_argument("--pvacseq-out", type=Path, required=True)
    parser.add_argument("--table-out", type=Path, required=True)
    args = parser.parse_args(argv)

    calls = parse_result(args.result.read_text())
    missing = [g for g in CLASS_I if not calls.get(g)]
    if missing:
        sys.exit(f"HLA-HD gave no usable call for class I gene(s) {', '.join(missing)} in {args.result}: "
                 "check read depth over the MHC, or supply hla_alleles for this patient in the samplesheet")

    alleles = pvacseq_alleles(calls)
    args.pvacseq_out.write_text(",".join(alleles) + "\n")
    with args.table_out.open("w") as handle:
        handle.write("gene\tallele_1\tallele_2\n")
        for gene in CLASS_I + DR_BETA + tuple(g for pair in HETERODIMERS for g in pair):
            if gene in calls:
                found = calls[gene] + ["", ""]
                handle.write(f"{gene}\t{found[0]}\t{found[1]}\n")
    print(f"[INFO] {len(alleles)} pVACseq alleles: {','.join(alleles)}")


if __name__ == "__main__":
    main()
