#!/usr/bin/env python3
"""
Group a patient's variants into clonal clusters by CCF, and say which ones
are still transcribed at week 6.

The candidate tables carry a CCF per variant but no structure: whether
this tumour is one clone or three, how much of it a given candidate
actually covers, and whether the clones that were there at baseline are
still there after treatment. For a neoadjuvant vaccine that last question
is the point - a peptide against a clone that the checkpoint inhibitor has
already cleared is a wasted slot, and one against a clone that survived is
the one worth making.

Clustering is deliberately plain: CCFs are sorted and cut wherever
consecutive values are more than --gap apart, then clusters below
--min-cluster variants are merged into their nearest neighbour. No mixture
model, no priors. With a few hundred exonic variants and a point-estimate
CCF per variant (no confidence interval), a fitted model would read as
more certain than the input supports; a gap-based cut can be checked by
eye against the histogram.

Persistence is read from the RNA columns rather than from DNA, because
these patients have RNA at baseline and week 6 but tumour DNA at baseline
only. That is weaker evidence: losing RNA support can mean the clone is
gone, or that the gene stopped being transcribed, or that the biopsy
caught less tumour. The status says which of those the data can and cannot
distinguish - "not_detected" means covered and absent, "no_coverage"
means the question was not asked.

A real phylogeny (a tumour evolution tree) needs several regions or
timepoints of tumour DNA. With one DNA timepoint there is nothing to build
a tree from, so this stops at clusters and does not draw one.

Outputs:
    <prefix>.clonal_structure.tsv   one row per cluster
    <prefix>.clonal_variants.tsv    one row per variant, with its cluster
    <prefix>.clonal_structure.png   CCF histogram and persistence (unless --no-plot)

Usage:
    clonal_structure.py P1.Combined.prioritized.tsv --prefix P1 [--no-plot]
"""
import argparse
import csv
import sys
from pathlib import Path

CLUSTER_COLUMNS = ["Cluster", "Label", "CCF Min", "CCF Max", "CCF Mean", "Variants",
                   "Candidates Kept", "Week 6 Detected", "Week 6 Not Detected",
                   "Week 6 No Coverage", "Week 6 Persistence"]
VARIANT_COLUMNS = ["ID", "Gene", "AA Change", "CCF", "Cluster", "Cluster Label",
                   "Baseline RNA", "Week 6 RNA", "Change", "Selection", "Priority Rank"]


def number(raw):
    if raw is None:
        return None
    text = str(raw).strip()
    if text in ("", "NA", ".", "None", "-"):
        return None
    try:
        return float(text)
    except ValueError:
        return None


def cluster_ccfs(values, gap, min_cluster):
    """Sorted CCFs -> [[values], ...], cut at gaps and small clusters merged."""
    if not values:
        return []
    ordered = sorted(values)
    clusters = [[ordered[0]]]
    for value in ordered[1:]:
        if value - clusters[-1][-1] > gap:
            clusters.append([value])
        else:
            clusters[-1].append(value)
    # Merge anything too small to be a clone into whichever neighbour is
    # closer in CCF, so a single stray variant does not become a "clone".
    while len(clusters) > 1:
        smallest = min(range(len(clusters)), key=lambda i: len(clusters[i]))
        if len(clusters[smallest]) >= min_cluster:
            break
        if smallest == 0:
            target = 1
        elif smallest == len(clusters) - 1:
            target = smallest - 1
        else:
            mean = sum(clusters[smallest]) / len(clusters[smallest])
            below = abs(mean - sum(clusters[smallest - 1]) / len(clusters[smallest - 1]))
            above = abs(mean - sum(clusters[smallest + 1]) / len(clusters[smallest + 1]))
            target = smallest - 1 if below <= above else smallest + 1
        clusters[target] = sorted(clusters[target] + clusters[smallest])
        clusters.pop(smallest)
    # Highest CCF first: cluster 1 is the one closest to being in every cell.
    return sorted(clusters, key=lambda c: -sum(c) / len(c))


def rna_status(depth, vaf, min_depth):
    """What the RNA at one timepoint can say about this variant."""
    if depth is None or depth < min_depth:
        return "no_coverage"
    if vaf is None or vaf <= 0:
        return "not_detected"
    return "detected"


def change(baseline, week6):
    """Baseline vs week 6, only where both were actually measured."""
    if baseline == "no_coverage" or week6 == "no_coverage":
        return "not_evaluable"
    if baseline == "detected" and week6 == "detected":
        return "persisted"
    if baseline == "detected" and week6 == "not_detected":
        return "lost"
    if baseline == "not_detected" and week6 == "detected":
        return "gained"
    return "absent_throughout"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("report", type=Path, help="an annotated candidate table carrying a CCF column")
    parser.add_argument("--prefix", required=True, help="output path prefix")
    parser.add_argument("--gap", type=float, default=0.15, help="CCF gap that starts a new cluster (default 0.15)")
    parser.add_argument("--min-cluster", type=int, default=3, help="variants a cluster needs to stand on its own (default 3)")
    parser.add_argument("--clonal-ccf", type=float, default=0.85, help="cluster mean CCF at or above which it is called clonal (default 0.85)")
    parser.add_argument("--min-rna-depth", type=float, default=5, help="RNA depth needed to call a variant present or absent (default 5)")
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args(argv)

    with args.report.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    # One row per variant, not per epitope: the same variant can appear
    # once per MHC class, and clonal structure is a property of variants.
    by_variant = {}
    for row in rows:
        key = (row.get("ID") or "").strip()
        if not key or key in by_variant:
            continue
        by_variant[key] = row

    ccfs = {key: number(row.get("CCF")) for key, row in by_variant.items()}
    clusters = cluster_ccfs([v for v in ccfs.values() if v is not None], args.gap, args.min_cluster)
    bounds = [(min(c), max(c)) for c in clusters]

    def cluster_of(ccf):
        if ccf is None:
            return None
        for index, (low, high) in enumerate(bounds, start=1):
            if low <= ccf <= high:
                return index
        return min(range(1, len(bounds) + 1),
                   key=lambda i: min(abs(ccf - bounds[i - 1][0]), abs(ccf - bounds[i - 1][1])),
                   default=None)

    labels = {}
    for index, values in enumerate(clusters, start=1):
        mean = sum(values) / len(values)
        labels[index] = "clonal" if mean >= args.clonal_ccf else "subclonal"

    variant_rows, per_cluster = [], {}
    for key, row in by_variant.items():
        ccf = ccfs[key]
        index = cluster_of(ccf)
        baseline = rna_status(number(row.get("RNA Depth")), number(row.get("RNA VAF")), args.min_rna_depth)
        week6 = rna_status(number(row.get("RNA Depth (week 6)")), number(row.get("RNA VAF (week 6)")), args.min_rna_depth)
        kept = (row.get("Selection") or "").strip() != "Rejected"
        variant_rows.append({
            "ID": key, "Gene": row.get("Gene", ""), "AA Change": row.get("AA Change", ""),
            "CCF": "NA" if ccf is None else f"{ccf:.3f}",
            "Cluster": "NA" if index is None else index,
            "Cluster Label": labels.get(index, "unclustered"),
            "Baseline RNA": baseline, "Week 6 RNA": week6,
            "Change": change(baseline, week6),
            "Selection": row.get("Selection", ""), "Priority Rank": row.get("Priority Rank", ""),
        })
        if index is not None:
            bucket = per_cluster.setdefault(index, {"variants": 0, "kept": 0, "detected": 0,
                                                    "not_detected": 0, "no_coverage": 0})
            bucket["variants"] += 1
            bucket["kept"] += 1 if kept else 0
            bucket[week6] += 1

    variant_rows.sort(key=lambda r: (r["Cluster"] == "NA", r["Cluster"], r["Gene"]))

    cluster_rows = []
    for index, values in enumerate(clusters, start=1):
        counts = per_cluster.get(index, {"variants": len(values), "kept": 0, "detected": 0,
                                         "not_detected": 0, "no_coverage": 0})
        evaluable = counts["detected"] + counts["not_detected"]
        cluster_rows.append({
            "Cluster": index, "Label": labels[index],
            "CCF Min": f"{min(values):.3f}", "CCF Max": f"{max(values):.3f}",
            "CCF Mean": f"{sum(values) / len(values):.3f}",
            "Variants": counts["variants"], "Candidates Kept": counts["kept"],
            "Week 6 Detected": counts["detected"],
            "Week 6 Not Detected": counts["not_detected"],
            "Week 6 No Coverage": counts["no_coverage"],
            "Week 6 Persistence": f"{counts['detected'] / evaluable:.2f}" if evaluable else "NA",
        })

    def write(path, columns, data):
        with open(path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                    lineterminator="\n", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(data)

    write(f"{args.prefix}.clonal_structure.tsv", CLUSTER_COLUMNS, cluster_rows)
    write(f"{args.prefix}.clonal_variants.tsv", VARIANT_COLUMNS, variant_rows)

    if not args.no_plot and clusters:
        plot(args.prefix, clusters, labels, cluster_rows)

    print(f"{len(by_variant)} variants in {len(clusters)} CCF cluster(s)", file=sys.stderr)
    for row in cluster_rows:
        print(f"  cluster {row['Cluster']} ({row['Label']}): {row['Variants']} variants, "
              f"CCF {row['CCF Min']}-{row['CCF Max']}, week 6 persistence {row['Week 6 Persistence']}",
              file=sys.stderr)


def plot(prefix, clusters, labels, cluster_rows):
    """CCF histogram by cluster, and week 6 persistence per cluster."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("note: matplotlib is not available, skipping the figure", file=sys.stderr)
        return

    figure, (top, bottom) = plt.subplots(2, 1, figsize=(7, 6), constrained_layout=True)
    colours = plt.get_cmap("tab10").colors
    for index, values in enumerate(clusters):
        top.hist(values, bins=20, range=(0, 1), alpha=0.75, color=colours[index % len(colours)],
                 label=f"cluster {index + 1} ({labels[index + 1]}, n={len(values)})")
    top.set_xlabel("cancer cell fraction")
    top.set_ylabel("variants")
    top.set_title("Clonal structure")
    top.legend(fontsize="small")

    evaluable = [row for row in cluster_rows if row["Week 6 Persistence"] != "NA"]
    if evaluable:
        positions = range(len(evaluable))
        bottom.bar(positions, [float(row["Week 6 Persistence"]) for row in evaluable], width=0.5,
                   color=[colours[(int(row["Cluster"]) - 1) % len(colours)] for row in evaluable])
        bottom.set_xlim(-0.75, len(evaluable) - 0.25)
        bottom.set_xticks(list(positions))
        bottom.set_xticklabels([f"cluster {row['Cluster']}\n{row['Label']}" for row in evaluable])
        bottom.set_ylim(0, 1)
        bottom.set_ylabel("fraction still detected in RNA")
        bottom.set_title("Week 6 persistence (RNA-evaluable variants only)")
    else:
        bottom.text(0.5, 0.5, "no variant had RNA coverage at both timepoints",
                    ha="center", va="center", transform=bottom.transAxes)
        bottom.set_axis_off()

    figure.savefig(f"{prefix}.clonal_structure.png", dpi=150)
    figure.savefig(f"{prefix}.clonal_structure.pdf")
    plt.close(figure)


if __name__ == "__main__":
    main()
