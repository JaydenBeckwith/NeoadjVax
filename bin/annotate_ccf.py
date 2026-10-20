#!/usr/bin/env python3
"""
Add copy-number and cancer-cell-fraction (CCF) columns to a pVACseq
aggregated report, from Sequenza purity, ploidy and (optionally) segments.

pVACseq takes purity (--tumor-purity) but has no notion of copy number, so
a variant on an amplified or lost segment gets the wrong clonality call
from VAF alone. For each candidate this looks up the tumour copy number at
the variant and applies the standard purity/copy-number correction
(e.g. McGranahan et al. 2016, Science):

    u   = VAF * (purity * CN_tumour + CN_normal * (1 - purity)) / purity
    m   = max(1, round(u)), capped at the major copy number when known
    CCF = min(1, u / m)

u is the number of mutated copies per tumour cell, m the estimated
mutation multiplicity. CN_tumour is, in order of preference: the Sequenza
segment total copy number (CNt) covering the variant; the sample ploidy
(when no segments file is given, or the variant falls between segments);
otherwise the normal copy number. CN_normal is 2, except chrX/chrY in a
male sample where it is 1 (pseudo-autosomal regions are not treated
separately). The "CN Source" column records which was used.

This is a point estimate with no confidence interval: low-depth variants
carry large VAF noise, so read it alongside the "DNA VAF" and depth columns.

Usage:
    annotate_ccf.py tumor.all_epitopes.aggregated.tsv --purity 0.62 \\
        [--ploidy 3.1] [--segments sample_segments.txt] [--sex male] -o out.tsv
"""
import argparse
import bisect
import csv
import math
from pathlib import Path

NEW_COLUMNS = ["Tumor CN", "Major CN", "Minor CN", "CN Source", "Normal CN", "Mutated Copies", "Multiplicity", "CCF"]


def norm_chrom(chrom):
    chrom = chrom.strip().strip('"')
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def read_segments(path):
    """Sequenza *_segments.txt -> {chrom: sorted [(start, end, CNt, A, B)]}."""
    segments = {}
    with open(path, newline="") as handle:
        rows = csv.reader(handle, delimiter="\t")
        header = [h.strip().strip('"') for h in next(rows)]
        col = {name: header.index(name) for name in ("chromosome", "start.pos", "end.pos", "CNt", "A", "B")}
        for row in rows:
            if not row or not row[0].strip():
                continue
            def value(name, cast):
                raw = row[col[name]].strip().strip('"')
                return cast(raw) if raw not in ("", "NA") else None
            segments.setdefault(norm_chrom(row[col["chromosome"]]), []).append((
                value("start.pos", int), value("end.pos", int),
                value("CNt", float), value("A", float), value("B", float)))
    for chrom in segments:
        segments[chrom].sort()
    return segments


def find_segment(segments, chrom, pos):
    rows = segments.get(norm_chrom(chrom), [])
    i = bisect.bisect_right([r[0] for r in rows], pos) - 1
    if i >= 0 and rows[i][0] <= pos <= rows[i][1] and rows[i][2] is not None:
        return rows[i]
    return None


def variant_position(variant_id):
    """pVACseq ID 'chrom-start-stop-ref-alt' -> (chrom, 1-based position)."""
    chrom, start, stop = variant_id.split("-", 3)[:3]
    return chrom, int(stop)


def ccf_columns(vaf, chrom, pos, purity, ploidy, sex, segments):
    normal_cn = 1 if sex == "male" and norm_chrom(chrom) in ("X", "Y") else 2
    major = minor = None
    segment = find_segment(segments, chrom, pos) if segments else None
    if segment:
        cn, major, minor, source = segment[2], segment[3], segment[4], "sequenza_segment"
    elif ploidy is not None:
        cn, source = ploidy, "ploidy"
    else:
        cn, source = float(normal_cn), "assumed_normal"
    fmt = lambda x: "NA" if x is None else f"{x:.3g}"
    base = [fmt(cn), fmt(major), fmt(minor), source, str(normal_cn)]
    if vaf is None:
        return base + ["NA", "NA", "NA"]
    mutated_copies = vaf * (purity * cn + normal_cn * (1 - purity)) / purity
    multiplicity = max(1, round(mutated_copies))
    if major is not None and major >= 1:
        multiplicity = min(multiplicity, int(major))
    ccf = min(1.0, mutated_copies / multiplicity)
    return base + [f"{mutated_copies:.3f}", str(multiplicity), f"{ccf:.3f}"]


def parse_float(raw):
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(value) else value


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("report", type=Path)
    parser.add_argument("--purity", type=float, required=True)
    parser.add_argument("--ploidy", type=float)
    parser.add_argument("--segments", type=Path)
    parser.add_argument("--sex", type=str.lower, choices=["male", "female"])
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if not 0 < args.purity <= 1:
        parser.error("--purity must be in (0, 1]")

    segments = read_segments(args.segments) if args.segments else None
    with args.report.open(newline="") as src, args.output.open("w", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
        writer.writerow(reader.fieldnames + NEW_COLUMNS)
        for row in reader:
            chrom, pos = variant_position(row["ID"])
            extra = ccf_columns(parse_float(row.get("DNA VAF")), chrom, pos,
                                args.purity, args.ploidy, args.sex, segments)
            writer.writerow([row[f] for f in reader.fieldnames] + extra)


if __name__ == "__main__":
    main()
