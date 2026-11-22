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
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import add_rna_timepoint  # noqa: E402  (same directory, on PATH in the pipeline)
import pvactools_output  # noqa: E402

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
    """pVACseq ID 'chrom-start-stop-ref-alt' -> (chrom, 1-based position).

    None for an ID that is not a variant position, which is what a
    pVACfuse report carries: its candidates are named after the fusion,
    not a locus, and they have no DNA VAF to turn into a CCF either.
    """
    fields = str(variant_id).split("-", 3)
    if len(fields) < 3:
        return None
    try:
        return fields[0], int(fields[2])
    except ValueError:
        return None


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


def parse_timepoint(spec):
    """'Week 6=path.vcf' -> ('Week 6', Path). Mirrors rna_support.py."""
    label, separator, path = str(spec).partition("=")
    if not separator or not label.strip() or not path.strip():
        raise SystemExit(f"--timepoint needs LABEL=VCF, got: {spec}")
    return label.strip(), Path(path.strip())


def annotate_report(source, destination, args, segments):
    """One report in, one report out. -> (rows written, rows with a CCF)."""
    written = scored = 0
    with open(source, newline="") as src, open(destination, "w", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        fieldnames = reader.fieldnames or []
        writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
        writer.writerow(fieldnames + NEW_COLUMNS)
        for row in reader:
            position = variant_position(row.get("ID", ""))
            if position is None:
                # A candidate with no locus (a fusion): copy number and CCF
                # are not missing values, they do not apply to it.
                extra = ["NA", "NA", "NA", "not_applicable", "NA", "NA", "NA", "NA"]
            else:
                chrom, pos = position
                extra = ccf_columns(parse_float(row.get("DNA VAF")), chrom, pos,
                                    args.purity, args.ploidy, args.sex, segments)
                scored += extra[-1] != "NA"
            writer.writerow([row.get(f, "") for f in fieldnames] + extra)
            written += 1
    return written, scored


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("report", type=Path,
                        help="an aggregated report, or a pvacseq/ or pvacfuse/ output directory")
    parser.add_argument("--purity", type=float, required=True)
    parser.add_argument("--ploidy", type=float)
    parser.add_argument("--segments", type=Path)
    parser.add_argument("--sex", type=str.lower, choices=["male", "female"])
    parser.add_argument("--timepoint", action="append", default=[], metavar="LABEL=VCF",
                        help="an RNA readcount VCF to join on, e.g. 'Week 6=P1.rna_week6.readcounts.vcf'")
    parser.add_argument("-o", "--output", type=Path, help="output file (single report only)")
    parser.add_argument("--outdir", type=Path, help="output directory (required for a directory input)")
    args = parser.parse_args(argv)
    if not 0 < args.purity <= 1:
        parser.error("--purity must be in (0, 1]")

    reports = pvactools_output.find_reports(args.report)
    if not reports:
        raise SystemExit(f"{args.report}: no aggregated report found "
                         f"(expected a .tsv, or a directory holding {pvactools_output.AGGREGATED_GLOB})")
    if len(reports) > 1 and not args.outdir:
        raise SystemExit(f"{args.report} holds {len(reports)} reports: give --outdir, not -o")
    if len(reports) == 1 and not (args.output or args.outdir):
        parser.error("give -o (one report) or --outdir")
    if args.outdir:
        args.outdir.mkdir(parents=True, exist_ok=True)

    timepoints = [parse_timepoint(spec) for spec in args.timepoint]
    for label, vcf in timepoints:
        if not vcf.exists():
            raise SystemExit(f"--timepoint {label}: no such file: {vcf}")

    segments = read_segments(args.segments) if args.segments else None
    for mhc, path in reports:
        if pvactools_output.kind(pvactools_output.header(path)) == "fusion":
            print(f"{path.name}: fusion candidates are called from RNA and carry no DNA VAF, so cancer "
                  "cell fraction does not apply; writing the copy-number source as not_applicable",
                  file=sys.stderr)
        destination = args.output if args.output and len(reports) == 1 else \
            args.outdir / f"{path.name.replace('.tsv', '')}.ccf.tsv"
        written, scored = annotate_report(path, destination, args, segments)
        for tp_label, vcf in timepoints:
            add_rna_timepoint.main([str(destination), "--vcf", str(vcf), "--label", tp_label,
                                    "-o", str(destination)])
        print(f"{mhc}: {written} candidates, {scored} with a CCF -> {destination}", file=sys.stderr)
        print(destination)


if __name__ == "__main__":
    main()
