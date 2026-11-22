#!/usr/bin/env python3
"""
Agretopicity and DAI for pVACseq candidates, straight from a pVACseq
output directory.

Agretopicity is the mutant peptide's IC50 divided by its wild-type
counterpart's. Below 1 the mutation is what creates the epitope, so the
T-cell repertoire against it is less likely to have been tolerised in the
thymus. Near or above 1 the wild-type peptide binds about as well and is
presented by normal cells too, which makes the mutant a poor target
however strong its own binding looks.

DAI is the same comparison as a difference, IC50 WT - IC50 MT. The ratio
is unstable when both peptides are weak binders (900 vs 1000 nM is a
ratio of 0.9 and a difference of 100), so both are reported and the
difference is the one to read at the weak end.

    Agretopicity Class   foreign        at or below --foreign (0.1)
                         intermediate   between the two thresholds
                         self_like      at or above --self-like (0.5)
                         NA             no wild-type IC50 to compare

Fusions have no wild-type counterpart at all: a fusion junction peptide
does not exist in the germline, so there is nothing to divide by. Pointed
at pVACfuse output this says so and writes the class as `not_applicable`
rather than filling the column with NA, which would read like a failed
calculation rather than one that does not apply.

RNA timepoints can be attached at the same time with --timepoint, so one
command gives a table carrying agretopicity and the baseline and week 6
RNA support together, ready for the report.

Usage:
    # straight from the pVACseq output directory, one file out per class
    agretopicity.py results/P1/pvacseq --outdir P1_agretopicity/

    # one report, with both RNA timepoints joined on
    agretopicity.py P1.Combined.aggregated.tsv -o P1.agretopicity.tsv \\
        --timepoint "Baseline=P1.somatic.vep.readcounts.vcf.gz" \\
        --timepoint "Week 6=P1.rna_week6.readcounts.vcf"
"""
import argparse
import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import add_rna_timepoint  # noqa: E402  (same directory, on PATH in the pipeline)
import pvactools_output  # noqa: E402

NEW_COLUMNS = ["Agretopicity", "DAI", "Agretopicity Class"]


def number(raw):
    if raw is None:
        return None
    text = str(raw).strip()
    if text in ("", "NA", ".", "None", "-"):
        return None
    try:
        value = float(text)
    except ValueError:
        return None
    if math.isnan(value) or math.isinf(value):
        return None
    return value


def agretopicity(mutant, wild_type):
    """-> (ratio, DAI). Either is None when it cannot be computed.

    A wild-type IC50 of 0 would make the ratio infinite; that is a
    degenerate prediction rather than an infinitely foreign peptide, so it
    is reported as missing.
    """
    if mutant is None or wild_type is None:
        return None, None
    difference = wild_type - mutant
    if wild_type <= 0:
        return None, difference
    return mutant / wild_type, difference


def classify(ratio, foreign, self_like):
    if ratio is None:
        return "NA"
    if ratio <= foreign:
        return "foreign"
    if ratio >= self_like:
        return "self_like"
    return "intermediate"


def annotate(rows, fieldnames, foreign, self_like, applicable=True):
    """-> [row + the three columns], as lists in fieldnames order."""
    out = []
    for row in rows:
        if applicable:
            ratio, dai = agretopicity(number(row.get("IC50 MT")), number(row.get("IC50 WT")))
            label = classify(ratio, foreign, self_like)
        else:
            ratio, dai, label = None, None, "not_applicable"
        fmt = lambda x: "NA" if x is None else f"{x:.4g}"
        out.append([row.get(name, "") for name in fieldnames] + [fmt(ratio), fmt(dai), label])
    return out


def read_report(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), (reader.fieldnames or [])


def write_report(path, fieldnames, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fieldnames + NEW_COLUMNS)
        writer.writerows(rows)


def parse_timepoint(spec):
    """'Week 6=path.vcf' -> ('Week 6', Path). Mirrors rna_support.py."""
    label, separator, path = str(spec).partition("=")
    if not separator or not label.strip() or not path.strip():
        raise SystemExit(f"--timepoint needs LABEL=VCF, got: {spec}")
    return label.strip(), Path(path.strip())


def attach_timepoints(path, timepoints):
    """Join each RNA timepoint's depth and VAF onto the report in place."""
    for label, vcf in timepoints:
        if not vcf.exists():
            raise SystemExit(f"--timepoint {label}: no such file: {vcf}")
        add_rna_timepoint.main([str(path), "--vcf", str(vcf), "--label", label, "-o", str(path)])


def summarise(rows, fieldnames):
    counts = {}
    index = len(fieldnames) + NEW_COLUMNS.index("Agretopicity Class")
    for row in rows:
        counts[row[index]] = counts.get(row[index], 0) + 1
    return ", ".join(f"{count} {name}" for name, count in sorted(counts.items()))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("report", type=Path,
                        help="an aggregated report, or a pvacseq/ or pvacfuse/ output directory")
    parser.add_argument("-o", "--output", type=Path, help="output file (single report only)")
    parser.add_argument("--outdir", type=Path, help="output directory (required for a directory input)")
    parser.add_argument("--foreign", type=float, default=0.1,
                        help="at or below this ratio the mutant is called foreign (default 0.1)")
    parser.add_argument("--self-like", type=float, default=0.5,
                        help="at or above this ratio the mutant is no better than self (default 0.5)")
    parser.add_argument("--timepoint", action="append", default=[], metavar="LABEL=VCF",
                        help="an RNA readcount VCF to join on, e.g. 'Week 6=P1.rna_week6.readcounts.vcf'")
    args = parser.parse_args(argv)
    if args.foreign > args.self_like:
        parser.error("--foreign cannot be above --self-like")

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
    written = []
    for label, path in reports:
        rows, fieldnames = read_report(path)
        kind = pvactools_output.kind(fieldnames)
        applicable = kind != "fusion"
        if not applicable:
            print(f"{path.name}: fusion candidates have no wild-type peptide, so agretopicity does not "
                  "apply; writing the class as not_applicable", file=sys.stderr)
        annotated = annotate(rows, fieldnames, args.foreign, args.self_like, applicable)

        destination = args.output if args.output and len(reports) == 1 else \
            args.outdir / f"{path.name.replace('.tsv', '')}.agretopicity.tsv"
        write_report(destination, fieldnames, annotated)
        attach_timepoints(destination, timepoints)
        written.append(destination)
        print(f"{label}: {len(annotated)} candidates -> {destination}"
              f"{' (' + summarise(annotated, fieldnames) + ')' if annotated else ''}", file=sys.stderr)

    for path in written:
        print(path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
