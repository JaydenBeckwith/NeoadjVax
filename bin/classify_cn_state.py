#!/usr/bin/env python3
"""
Classify each segment of a Sequenza *_segments.txt file into a copy-number
state, using Jayden's stated rule for the neoadjuvant LOH analysis
([[loh-analysis]]):

    HETEROZYGOUS   = B > 0
    LOH_DELETION   = B == 0 and CNt == 1   (1+0)
    CNN_LOH        = B == 0 and CNt == 2   (2+0)
    LOH_AMPLIFIED  = B == 0 and CNt > 2
    (CNt == 0, biallelic loss, is left UNCLASSIFIED and written to a
     separate file rather than dropped silently: this is exactly the
     "heavy/pervasive genome-wide biallelic loss in 3 of 98 melanoma
     samples" pattern already flagged as worth a second look, possibly
     whole-genome-doubling related)


Usage:
    classify_cn_state.py segments.txt -o classified_segments.tsv \\
        [--biallelic-loss-out biallelic_loss.tsv]
"""
import argparse
import csv
import sys


def classify(cnt, b):
    if cnt == 0:
        return "UNCLASSIFIED_BIALLELIC_LOSS"
    if b > 0:
        return "HETEROZYGOUS"
    if cnt == 1:
        return "LOH_DELETION"
    if cnt == 2:
        return "CNN_LOH"
    return "LOH_AMPLIFIED"  # b == 0 and cnt > 2


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("segments_file", help="Sequenza *_segments.txt (tab-separated, with a header row)")
    p.add_argument("-o", "--out", required=True)
    p.add_argument("--biallelic-loss-out", default=None,
                    help="separate file for CNt==0 segments, only written if non-empty")
    p.add_argument("--cnt-col", default="CNt")
    p.add_argument("--b-col", default="B")
    args = p.parse_args()

    with open(args.segments_file, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if args.cnt_col not in reader.fieldnames or args.b_col not in reader.fieldnames:
            sys.exit(f"[ERROR] Expected columns '{args.cnt_col}' and '{args.b_col}' in {args.segments_file}, "
                      f"found: {reader.fieldnames}. If your Sequenza output uses different column names, "
                      f"pass --cnt-col/--b-col.")
        rows = list(reader)

    biallelic_rows = []
    with open(args.out, "w", newline="") as out_fh:
        fieldnames = reader.fieldnames + ["cn_state"]
        writer = csv.DictWriter(out_fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            try:
                cnt = int(float(row[args.cnt_col]))
                b = int(float(row[args.b_col]))
            except (ValueError, TypeError):
                print(f"[WARN] Skipping unparseable segment: {row}", file=sys.stderr)
                continue
            state = classify(cnt, b)
            row["cn_state"] = state
            writer.writerow(row)
            if state == "UNCLASSIFIED_BIALLELIC_LOSS":
                biallelic_rows.append(row)

    print(f"[INFO] {len(rows)} segments classified -> {args.out} "
          f"({len(biallelic_rows)} biallelic-loss / unclassified)", file=sys.stderr)

    if biallelic_rows and args.biallelic_loss_out:
        with open(args.biallelic_loss_out, "w", newline="") as out_fh:
            writer = csv.DictWriter(out_fh, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(biallelic_rows)


if __name__ == "__main__":
    main()
