#!/usr/bin/env python3
"""
Filter a VCF (plain or .gz) down to standard chromosomes (1-22, X, Y, MT/M),
with or without a "chr" prefix.

Ported verbatim from Jayden's original DNA pVACseq script
(filter_vcf_standard_chroms) so pipeline behaviour is unchanged: this now
runs as its own step (FILTER_STANDARD_CHROMS) instead of an inline function.

Usage: filter_vcf_standard_chroms.py <in.vcf[.gz]> <out.vcf>
"""
import gzip
import sys
from pathlib import Path


def filter_vcf_standard_chroms(vcf_path: Path, output_path: Path):
    standard_chroms = {str(i) for i in range(1, 23)} | {"X", "Y", "MT", "M"}
    open_func = gzip.open if vcf_path.suffix == ".gz" else open
    mode = "rt" if vcf_path.suffix == ".gz" else "r"
    with open_func(vcf_path, mode) as infile, output_path.open("w") as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
            else:
                chrom = line.split()[0].replace("chr", "")
                if chrom in standard_chroms:
                    outfile.write(line)
    print(f"[INFO] Filtered VCF written to: {output_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: filter_vcf_standard_chroms.py <in.vcf[.gz]> <out.vcf>")
    filter_vcf_standard_chroms(Path(sys.argv[1]), Path(sys.argv[2]))
