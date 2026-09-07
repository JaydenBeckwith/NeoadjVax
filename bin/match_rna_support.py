#!/usr/bin/env python3
"""
Filter a somatic DNA VCF down to loci that also carry RNA support, per the
per-patient neoantigen score design in [[neoantigen-score]]:

    "somatic variant calling from WGS (tumour + normal) ... then mapping
    those variant loci to matched RNA via variant calling on the RNA sample
    to see how many carry over. That RNA-supported VCF subset then run
    through pVACseq."

This is new code: it did not exist in either of Jayden's original scripts
(the DNA script went straight from Mutect2 to VEP/pVACseq; the RNA script
stopped at its own filtered.vcf and never fed back into the DNA calls).

Matching is by exact (CHROM, POS, REF, ALT); CHROM naming (chr-prefixed or
not) is normalized before comparison since both branches use the same
GRCh38/GENCODE reference so this should be a non-issue in practice, but is
handled defensively anyway. A locus is kept if it appears in the RNA VCF
AND (if --min-rna-alt-reads > 0) the RNA allelic depth for that ALT allele
meets the threshold, read from the record's AD FORMAT field.

Usage:
    match_rna_support.py --dna-vcf DNA.vcf --rna-vcf RNA.vcf \\
        --out-vcf out.vcf [--min-rna-alt-reads 3] [--summary summary.tsv]
"""
import argparse
import gzip
import sys
from pathlib import Path


def _open(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r")


def _norm_chrom(chrom):
    return chrom[3:] if chrom.startswith("chr") else chrom


def _parse_ad(format_field, sample_field, alt_index):
    """Pull the ALT-allele read count out of a FORMAT/sample AD field.
    Returns None if AD isn't present or can't be parsed (record is kept
    in that case: absence of AD shouldn't silently drop a variant)."""
    try:
        keys = format_field.split(":")
        vals = sample_field.split(":")
        ad_idx = keys.index("AD")
        ad = vals[ad_idx].split(",")
        # AD[0] is the REF count, AD[1..] are ALT counts in ALT order
        return int(ad[alt_index + 1])
    except (ValueError, IndexError):
        return None


def load_rna_loci(rna_vcf_path, min_rna_alt_reads):
    """Return {(chrom, pos, ref, alt): rna_alt_depth_or_None} for every RNA
    call that passes the depth threshold (or has no parseable AD, kept
    permissively)."""
    rna_loci = {}
    with _open(rna_vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, _id, ref, alt_str = fields[0:5]
            chrom = _norm_chrom(chrom)
            format_field = fields[8] if len(fields) > 8 else None
            sample_field = fields[9] if len(fields) > 9 else None
            for alt_index, alt in enumerate(alt_str.split(",")):
                key = (chrom, pos, ref, alt)
                depth = None
                if format_field and sample_field:
                    depth = _parse_ad(format_field, sample_field, alt_index)
                if min_rna_alt_reads > 0 and depth is not None and depth < min_rna_alt_reads:
                    continue
                rna_loci[key] = depth
    return rna_loci


def filter_dna_vcf(dna_vcf_path, rna_loci, out_path, summary_path=None):
    total = 0
    kept = 0
    with _open(dna_vcf_path) as infile, open(out_path, "w") as outfile:
        summary_fh = open(summary_path, "w") if summary_path else None
        if summary_fh:
            summary_fh.write("CHROM\tPOS\tREF\tALT\tRNA_SUPPORTED\tRNA_ALT_DEPTH\n")
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            chrom, pos, _id, ref, alt_str = fields[0:5]
            norm_chrom = _norm_chrom(chrom)
            alts = alt_str.split(",")
            total += 1
            supported = False
            depth_seen = None
            for alt in alts:
                key = (norm_chrom, pos, ref, alt)
                if key in rna_loci:
                    supported = True
                    depth_seen = rna_loci[key]
                    break
            if supported:
                kept += 1
                outfile.write(line)
            if summary_fh:
                summary_fh.write(f"{chrom}\t{pos}\t{ref}\t{alt_str}\t{int(supported)}\t{depth_seen if depth_seen is not None else 'NA'}\n")
        if summary_fh:
            summary_fh.close()
    print(f"[INFO] {kept}/{total} DNA somatic loci carried RNA support -> {out_path}", file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dna-vcf", required=True)
    p.add_argument("--rna-vcf", required=True)
    p.add_argument("--out-vcf", required=True)
    p.add_argument("--min-rna-alt-reads", type=int, default=1,
                    help="Minimum RNA ALT-allele read depth to count as supported (0 = presence-only). Default: 1")
    p.add_argument("--summary", default=None, help="Optional per-locus TSV summary (for QC/reporting, not fed into pVACseq)")
    args = p.parse_args()

    rna_loci = load_rna_loci(args.rna_vcf, args.min_rna_alt_reads)
    filter_dna_vcf(args.dna_vcf, rna_loci, args.out_vcf, args.summary)


if __name__ == "__main__":
    main()
