#!/usr/bin/env python3
"""
Pull a single (purity, ploidy) "top solution" out of bin/run-sequenza.R's
output for one sample.

run-sequenza.R calls the sequenza R package's sequenza.extract() /
sequenza.fit() / sequenza.results() with sample.id=<sample>,
out.dir=<sample>_OUTPUT — none of that is guessed, it's read straight out
of the script you sent (bin/run-sequenza.R). What IS a targeted-but-
unverified assumption is what happens in THIS script: none of your four
Sequenza_tools scripts (nor run-sequenza.R itself) pick a single top
solution out of sequenza.results()' output — that step doesn't exist yet
anywhere in what you've sent, so it's new code, not a port.

sequenza.results() is documented (sequenza R package docs, not anything
specific to your setup) to write, among other files,
"<sample>_alternative_solutions.txt" — one row per candidate
(cellularity, ploidy) solution the fit considered, with a column scoring
each one (SLPP - scaled log posterior probability - in current sequenza
versions; older versions may spell/case this differently, hence the
tolerant column matching below). This script takes the row with the
BEST score in that column as the "top solution" — a defensible, mechanical
choice, but still an assumption about what "top solution" should mean for
your analysis, not something confirmed against NeoadjLOH's actual
sequenza_top_solutions_summary.csv (which you also haven't sent — if a
script already builds that file with different logic than "best-scoring
row", use that instead of this one).

RECOMMENDATION: before trusting this in a real HLA-LOH run, manually open
one sample's <sample>_alternative_solutions.txt and confirm the row this
script picks matches the solution you'd pick by eye (e.g. against
<sample>_CP_contours.pdf).
"""
import argparse
import csv
import glob
import os
import re
import sys


CELLULARITY_COL_PATTERNS = [r'^cellularity$']
PLOIDY_COL_PATTERNS = [r'^ploidy$', r'^ploidy\.estimate$', r'^ploidy_estimate$']
SCORE_COL_PATTERNS = [r'^slpp$', r'^lpp$', r'^score$']


def find_column(fieldnames, patterns):
    for pat in patterns:
        for fn in fieldnames:
            if re.match(pat, fn.strip(), re.IGNORECASE):
                return fn
    return None


def sniff_delimiter(path):
    with open(path, newline='') as fh:
        first_line = fh.readline()
    if '\t' in first_line:
        return '\t'
    return ','


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sample', required=True)
    p.add_argument('--results-dir', required=True,
                    help='directory containing <sample>_OUTPUT/ (SEQUENZA_FIT\'s raw_results output)')
    p.add_argument('--out', required=True, help='tab-separated purity\\tploidy, one row, written here')
    args = p.parse_args()

    out_dir_candidates = glob.glob(os.path.join(args.results_dir, f'{args.sample}_OUTPUT'))
    if not out_dir_candidates:
        # sequenza.results() uses sample_name verbatim as out.dir's prefix,
        # but be a little tolerant of a trailing-slash/case mismatch rather
        # than failing on something cosmetic
        out_dir_candidates = glob.glob(os.path.join(args.results_dir, f'{args.sample}*_OUTPUT'))
    if not out_dir_candidates:
        sys.exit(f"[ERROR] No '{args.sample}_OUTPUT' directory found under {args.results_dir} — "
                  f"run-sequenza.R may not have completed successfully for {args.sample}.")
    out_dir = out_dir_candidates[0]

    sol_file = os.path.join(out_dir, f'{args.sample}_alternative_solutions.txt')
    if not os.path.isfile(sol_file):
        sys.exit(f"[ERROR] Expected {sol_file} (sequenza.results()' standard alternative-solutions output) "
                  f"not found in {out_dir} — contents: {os.listdir(out_dir) if os.path.isdir(out_dir) else '<missing>'}. "
                  f"This script targets the standard sequenza R package output convention; if your actual output "
                  f"is named/shaped differently, this needs adjusting rather than guessing further.")

    delim = sniff_delimiter(sol_file)
    with open(sol_file, newline='') as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        rows = list(reader)
        fieldnames = reader.fieldnames or []

    if not rows:
        sys.exit(f"[ERROR] {sol_file} has a header but no data rows for {args.sample}.")

    cellularity_col = find_column(fieldnames, CELLULARITY_COL_PATTERNS)
    ploidy_col = find_column(fieldnames, PLOIDY_COL_PATTERNS)
    score_col = find_column(fieldnames, SCORE_COL_PATTERNS)

    if not cellularity_col or not ploidy_col:
        sys.exit(f"[ERROR] Could not find cellularity/ploidy columns in {sol_file} — found columns: {fieldnames}. "
                  f"Column-name matching is deliberately strict here rather than guessing which column is which.")

    if score_col:
        best_row = max(rows, key=lambda r: float(r[score_col]))
        print(f"[INFO] {args.sample}: picked best-scoring row by '{score_col}' out of {len(rows)} candidate solutions", file=sys.stderr)
    else:
        best_row = rows[0]
        print(f"[WARN] {args.sample}: no recognizable score column in {fieldnames} — falling back to the FIRST "
              f"row of {len(rows)} candidate solutions. Verify this is actually the top solution before trusting it.", file=sys.stderr)

    purity = best_row[cellularity_col]
    ploidy = best_row[ploidy_col]

    with open(args.out, 'w') as out_fh:
        out_fh.write(f"{purity}\t{ploidy}\n")

    print(f"[INFO] {args.sample}: purity={purity} ploidy={ploidy} (from {sol_file})", file=sys.stderr)


if __name__ == '__main__':
    main()
