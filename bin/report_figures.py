#!/usr/bin/env python3
"""
Draw the report's prediction figures, separately from the report itself.

patient_report.py writes LaTeX. It embeds figures but draws none, so
every plot is made here, from tables the pipeline already wrote, and
handed to the report as a file. The split is deliberate: a figure can be
regenerated and looked at on its own, from an existing results directory,
without running pdflatex or the rest of the report; and the data behind
each one is a plain function that can be tested with no matplotlib
installed at all.

Three figures, each from one table and each skipped when its table is
absent:

    candidates   The decision space: cancer cell fraction against
                 agretopicity, coloured by what prioritisation decided and
                 sized by expression. A good target sits top right -
                 clonal, and binding far better than its wild type. This
                 is the one figure that shows why the ranking came out as
                 it did rather than only what it decided.
    flags        How many candidates carry each flag, rejections drawn
                 apart from cautions. Answers "what did we lose them to?":
                 a patient losing everything to HLA_LOST is a different
                 problem from one losing everything to NO_EXPRESSION.
    hla_loh      One row per locus, coloured by call, with the alleles
                 written out. Small, but it is what a reader checks first
                 when a candidate was rejected on its allele.

Each is written as .png (for looking at) and .pdf (for LaTeX, where the
text stays sharp at any size).

Usage:
    report_figures.py --prefix P1 \\
        --prioritized P1.Combined.prioritized.tsv \\
        [--hla-loh P1.hla_loh.tsv] [--only candidates,flags]
"""
import argparse
import csv
import sys
from pathlib import Path

# Colours by what prioritisation decided. Not a red/green pair alone: the
# three also run light to dark, so the figure survives greyscale printing
# and the common colour-vision differences.
SELECTION_COLOURS = {
    "Top ranking": "#2b7a4b",
    "Caution": "#d99000",
    "Rejected": "#9aa0a6",
}
REJECT_COLOUR = "#b3261e"
CAUTION_COLOUR = "#d99000"
LOH_COLOURS = {
    "loss": "#b3261e",
    "imbalance": "#d99000",
    "retained": "#2b7a4b",
    "not_evaluable": "#9aa0a6",
    "low_purity": "#9aa0a6",
}

FIGURES = ("candidates", "flags", "hla_loh")

# Kept in step with bin/prioritize_candidates.py: a flag that rejects a
# candidate is drawn apart from one that only cautions.
REJECT_FLAGS = ("HLA_LOST", "LOW_BINDING", "NO_EXPRESSION", "NO_RNA_SUPPORT")


def number(raw):
    """A table cell -> float, or None for NA/./empty."""
    if raw is None:
        return None
    text = str(raw).strip()
    if text in ("", "NA", ".", "None", "-"):
        return None
    try:
        value = float(text)
    except ValueError:
        return None
    if value != value or abs(value) == float("inf"):
        return None
    return value


def read_tsv(path):
    if not path or not Path(path).exists():
        return []
    text = Path(path).read_text()
    if not text.strip():
        return []
    delimiter = "\t" if "\t" in text.splitlines()[0] else ","
    return list(csv.DictReader(text.splitlines(), delimiter=delimiter))


def one_row_per_variant(rows):
    """The first row for each variant ID, in order.

    A variant appears once per MHC class in a combined table, and these
    figures are about variants rather than epitopes.
    """
    seen, out = set(), []
    for row in rows:
        key = (row.get("ID") or "").strip()
        if key and key in seen:
            continue
        if key:
            seen.add(key)
        out.append(row)
    return out


# --- data preparation: pure functions, no matplotlib ---------------------

def landscape_points(rows):
    """-> ([point, ...], dropped) for the candidates figure.

    A point needs both axes. Candidates missing either are counted rather
    than dropped silently, so the figure can say how much of the list it
    is actually showing.
    """
    points, dropped = [], 0
    for row in one_row_per_variant(rows):
        ccf = number(row.get("CCF"))
        agretopicity = number(row.get("Agretopicity"))
        if ccf is None or agretopicity is None or agretopicity <= 0:
            dropped += 1
            continue
        expression = number(row.get("Sample TPM"))
        if expression is None:
            expression = number(row.get("RNA Expr"))
        points.append({
            "gene": (row.get("Gene") or "").strip(),
            "ccf": ccf,
            "agretopicity": agretopicity,
            "expression": expression,
            "selection": (row.get("Selection") or "").strip() or "Rejected",
            "rank": number(row.get("Priority Rank")),
        })
    return points, dropped


def flag_counts(rows):
    """-> [(flag, count, is_rejection)], most frequent first."""
    counts = {}
    for row in one_row_per_variant(rows):
        for flag in (row.get("Flags") or "").split(";"):
            flag = flag.strip()
            if flag:
                counts[flag] = counts.get(flag, 0) + 1
    return [(flag, count, flag in REJECT_FLAGS)
            for flag, count in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))]


def loh_entries(rows):
    """-> [(locus, call, label)] for the HLA-LOH figure, class I first."""
    entries = []
    for row in rows:
        locus = (row.get("locus") or "").strip()
        if not locus:
            continue
        call = (row.get("call") or "").strip() or "not_evaluable"
        a1 = (row.get("allele1_pvacseq") or row.get("allele1") or "").strip()
        a2 = (row.get("allele2_pvacseq") or row.get("allele2") or "").strip()
        lost = (row.get("lost_allele_pvacseq") or "").strip()
        label = f"{a1} / {a2}".strip(" /")
        if call == "loss" and lost:
            label = f"{label}  (lost {lost})"
        entries.append((locus, call, label))
    # Class I first, then class II, each alphabetically: the order a
    # reader expects rather than the order SpecHLA happened to write.
    first = ("A", "B", "C")
    return sorted(entries, key=lambda e: (e[0] not in first, e[0]))


# --- rendering -----------------------------------------------------------

def pyplot():
    """matplotlib's pyplot on a non-interactive backend, or None."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    return plt


def save(figure, prefix, name, formats):
    written = []
    for suffix in formats:
        path = f"{prefix}.{name}.{suffix}"
        figure.savefig(path, dpi=150)
        written.append(path)
    return written


def draw_candidates(plt, points, dropped, prefix, formats, top_labels=6):
    figure, axes = plt.subplots(figsize=(7.2, 5.0), constrained_layout=True)
    if not points:
        axes.text(0.5, 0.5, "no candidate had both a CCF and an agretopicity",
                  ha="center", va="center", transform=axes.transAxes)
        axes.set_axis_off()
        return save(figure, prefix, "candidates", formats), figure

    expressions = [p["expression"] for p in points if p["expression"] is not None]
    top_expression = max(expressions) if expressions else None

    def size(point):
        # Area carries expression. A candidate with no expression value is
        # drawn at the floor rather than left out, because missing data is
        # not the same as not expressed.
        if point["expression"] is None or not top_expression:
            return 18.0
        return 18.0 + 180.0 * (point["expression"] / top_expression) ** 0.5

    for selection, colour in SELECTION_COLOURS.items():
        group = [p for p in points if p["selection"] == selection]
        if not group:
            continue
        axes.scatter([p["ccf"] for p in group], [p["agretopicity"] for p in group],
                     s=[size(p) for p in group], c=colour, alpha=0.75,
                     edgecolors="white", linewidths=0.6, label=f"{selection} (n={len(group)})")

    ranked = [p for p in points if p["rank"] is not None]
    for point in sorted(ranked, key=lambda p: p["rank"])[:top_labels]:
        if point["gene"]:
            axes.annotate(point["gene"], (point["ccf"], point["agretopicity"]),
                          textcoords="offset points", xytext=(5, 4), fontsize=8)

    axes.set_yscale("log")
    # Inverted, so a low agretopicity (the good case) is high on the page:
    # the best candidates then sit top right, which is where a reader looks.
    axes.invert_yaxis()
    axes.set_xlim(0, 1.02)
    axes.axhline(1.0, color="#666666", linewidth=0.8, linestyle="--")
    axes.annotate("binds no better than wild type", (0.02, 1.0), fontsize=8, color="#666666",
                  textcoords="offset points", xytext=(0, -11))
    axes.set_xlabel("cancer cell fraction (clonal to the right)")
    axes.set_ylabel("agretopicity: IC50 mutant / wild type\n(log, inverted: more foreign is higher)")
    title = f"Candidate landscape ({len(points)} candidates"
    title += f", {dropped} without both axes)" if dropped else ")"
    axes.set_title(title)
    axes.legend(fontsize="small", loc="upper left")
    axes.grid(alpha=0.25, linewidth=0.5)
    return save(figure, prefix, "candidates", formats), figure


def draw_flags(plt, counts, prefix, formats):
    figure, axes = plt.subplots(figsize=(7.2, max(2.4, 0.42 * len(counts) + 1.4)),
                                constrained_layout=True)
    if not counts:
        axes.text(0.5, 0.5, "no candidate carried a flag", ha="center", va="center",
                  transform=axes.transAxes)
        axes.set_axis_off()
        return save(figure, prefix, "flags", formats), figure

    # Largest at the top, which means plotting in reverse.
    ordered = list(reversed(counts))
    positions = list(range(len(ordered)))
    axes.barh(positions, [c for _, c, _ in ordered],
              color=[REJECT_COLOUR if rejects else CAUTION_COLOUR for _, _, rejects in ordered],
              height=0.65)
    axes.set_yticks(positions)
    axes.set_yticklabels([flag for flag, _, _ in ordered], fontsize=9)
    for position, (_, count, _) in zip(positions, ordered):
        axes.annotate(str(count), (count, position), textcoords="offset points",
                      xytext=(4, 0), va="center", fontsize=8)
    axes.set_xlabel("candidates carrying the flag")
    # Counts, so no half-candidate tick marks.
    axes.xaxis.get_major_locator().set_params(integer=True)
    axes.set_title("Why candidates were flagged (red rejects, amber cautions)")
    axes.grid(axis="x", alpha=0.25, linewidth=0.5)
    return save(figure, prefix, "flags", formats), figure


def draw_hla_loh(plt, entries, prefix, formats):
    figure, axes = plt.subplots(figsize=(7.2, max(2.2, 0.4 * len(entries) + 1.2)),
                                constrained_layout=True)
    if not entries:
        axes.text(0.5, 0.5, "no HLA-LOH result for this patient", ha="center", va="center",
                  transform=axes.transAxes)
        axes.set_axis_off()
        return save(figure, prefix, "hla_loh", formats), figure

    ordered = list(reversed(entries))
    positions = list(range(len(ordered)))
    axes.barh(positions, [1] * len(ordered), height=0.6,
              color=[LOH_COLOURS.get(call, "#9aa0a6") for _, call, _ in ordered])
    for position, (_, call, label) in zip(positions, ordered):
        axes.annotate(f"{call}   {label}", (0.02, position), va="center", fontsize=8, color="white")
    axes.set_yticks(positions)
    axes.set_yticklabels([locus for locus, _, _ in ordered], fontsize=9)
    axes.set_xticks([])
    axes.set_xlim(0, 1)
    lost = sum(1 for _, call, _ in entries if call == "loss")
    axes.set_title(f"HLA loss of heterozygosity: {lost} of {len(entries)} loci lost")
    return save(figure, prefix, "hla_loh", formats), figure


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--prefix", required=True, help="output path prefix")
    parser.add_argument("--prioritized", type=Path, help="ranked table from prioritize_candidates.py")
    parser.add_argument("--hla-loh", type=Path, help="tidy table from parse_hla_loh.py")
    parser.add_argument("--only", help=f"comma-separated subset of: {', '.join(FIGURES)}")
    parser.add_argument("--format", default="png,pdf", help="comma-separated output formats (default png,pdf)")
    parser.add_argument("--top-labels", type=int, default=6, help="how many top-ranked genes to name (default 6)")
    args = parser.parse_args(argv)

    wanted = [f.strip() for f in (args.only or ",".join(FIGURES)).split(",") if f.strip()]
    unknown = [f for f in wanted if f not in FIGURES]
    if unknown:
        parser.error(f"--only: unknown figure {', '.join(unknown)}. Choose from: {', '.join(FIGURES)}")
    formats = [f.strip().lstrip(".") for f in args.format.split(",") if f.strip()]

    plt = pyplot()
    if plt is None:
        print("error: matplotlib is not available, so no figure can be drawn", file=sys.stderr)
        return 1

    candidates = read_tsv(args.prioritized)
    loh = read_tsv(args.hla_loh)
    written = []

    if "candidates" in wanted and args.prioritized:
        points, dropped = landscape_points(candidates)
        paths, figure = draw_candidates(plt, points, dropped, args.prefix, formats, args.top_labels)
        written += paths
        plt.close(figure)
    if "flags" in wanted and args.prioritized:
        paths, figure = draw_flags(plt, flag_counts(candidates), args.prefix, formats)
        written += paths
        plt.close(figure)
    if "hla_loh" in wanted and args.hla_loh:
        paths, figure = draw_hla_loh(plt, loh_entries(loh), args.prefix, formats)
        written += paths
        plt.close(figure)

    for path in written:
        print(path)
    if not written:
        print("note: nothing to draw (no input table given for the selected figures)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
