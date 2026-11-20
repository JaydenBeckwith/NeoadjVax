#!/usr/bin/env python3
"""
How many of a patient's non-synonymous somatic variants, and of its clonal
neoantigens, are seen in tumour RNA at each timepoint.

A DNA variant only becomes antigen if the mutated allele is transcribed,
so the share of mutations that "make it through" to RNA, and whether that
share holds up on treatment (baseline vs week 6), is a per-patient
read-out of how much of the candidate list is real. This is the
single-patient version of the cohort figure (RNA support of clonal
neoantigens, expressed both / week 6 only / baseline only / DNA only).

Inputs are the pipeline's own files:

  --vcf        the VEP + readcount annotated somatic VCF that pVACseq read
               (<id>.somatic.vep.readcounts.vcf.gz): CSQ for the
               consequence, tumour AD/DP/AF for DNA VAF and CCF
  --timepoint  LABEL=VCF, once per RNA sample, in time order. Each VCF has
               that sample's RNA counts (RAD/RDP/RAF) in its tumour
               column: for baseline that is the same readcounts VCF, for
               week 6 the <id>.rna_week6.readcounts.vcf copy
  --class-i / --class-ii
               pVACseq aggregated reports (the .ccf.tsv files are fine):
               a variant is a neoantigen for that class when its best
               peptide has IC50 MT or %ile MT at or under the thresholds

Definitions, all arguments:

  non-synonymous   the canonical transcript's consequence (any transcript
                   when VEP marked none canonical) is missense, nonsense,
                   frameshift, in-frame indel, start/stop lost or
                   protein-altering. Splice-region variants are not
                   counted.
  RNA supported    at least --min-alt-reads reads carry the alt allele in
                   that RNA sample, and RNA VAF >= --min-rna-vaf
  clonal           CCF >= --clonal-ccf, with CCF from DNA VAF, purity,
                   ploidy and the Sequenza segment, exactly as
                   annotate_ccf.py computes it for the pVACseq reports

A variant with no RNA reads at a site is "not supported", not missing:
the gene may simply be silent. The per-variant table keeps the depth so
the two can be told apart.

Outputs (prefix from -o):

  <prefix>.rna_support.variants.tsv  one row per non-synonymous variant
  <prefix>.rna_support.summary.tsv   counts per variant set and category
  <prefix>.rna_support.pdf / .png    the figure (unless --no-plot)

Usage:
    rna_support.py --patient P1 --vcf P1.somatic.vep.readcounts.vcf.gz \\
        --timepoint "Baseline=P1.somatic.vep.readcounts.vcf.gz" \\
        --timepoint "Week 6=P1.rna_week6.readcounts.vcf" \\
        --class-i P1.MHC_I.aggregated.ccf.tsv --class-ii P1.MHC_II.aggregated.ccf.tsv \\
        --purity 0.62 [--ploidy 3.1] [--segments P1_segments.txt] [--sex female] \\
        -o P1
"""
import argparse
import csv
import gzip
import math
import sys
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import annotate_ccf  # noqa: E402  (same directory, on PATH in the pipeline)

NON_SYNONYMOUS = ("missense_variant", "stop_gained", "stop_lost", "start_lost", "frameshift_variant",
                  "inframe_insertion", "inframe_deletion", "protein_altering_variant")

DNA_ONLY = "DNA only"

# Matches the cohort figure, so a single patient reads like one of its columns.
PALETTE_TWO = {"both": "#1F5A96", "later": "#8CCB7E", "first": "#E5AE0C", "none": "#C4B5D9"}
EXTRA_COLOURS = ["#2A9D8F", "#E76F51", "#6D597A", "#B5838D", "#457B9D", "#F4A261"]

AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
       "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
       "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*", "Sec": "U", "Xaa": "X"}


# ---------------------------------------------------------------------------
# VCF reading (shared with patient_report.py)
# ---------------------------------------------------------------------------
def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def parse_float(raw):
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(value) else value


def parse_int(raw):
    value = parse_float(raw)
    return None if value is None else int(round(value))


def short_protein_change(hgvsp):
    """'ENSP0001.1:p.Val600Glu' -> 'V600E'; anything unparseable is returned trimmed."""
    if not hgvsp:
        return ""
    change = hgvsp.split(":")[-1]
    if change.startswith("p."):
        change = change[2:]
    change = change.replace("%3D", "=")
    out, i = [], 0
    while i < len(change):
        token = change[i:i + 3]
        if token in AA3:
            out.append(AA3[token])
            i += 3
        else:
            out.append(change[i])
            i += 1
    return "".join(out)


def pick_csq(entries):
    """One CSQ entry per variant: canonical first, then MANE Select, then the first.

    Among the preferred transcripts a protein-altering consequence wins, so a
    variant that is missense on the canonical transcript is not reported by
    a non-coding transcript that happens to be listed first.
    """
    if not entries:
        return None
    for key in ("CANONICAL", "MANE_SELECT"):
        preferred = [e for e in entries if e.get(key) and e.get(key) not in ("", "-")]
        if preferred:
            break
    else:
        preferred = entries
    coding = [e for e in preferred if is_non_synonymous(e)]
    return (coding or preferred)[0]


def is_non_synonymous(entry):
    return bool(entry) and any(c in NON_SYNONYMOUS for c in (entry.get("Consequence") or "").split("&"))


def protein_change(entry):
    if not entry:
        return ""
    change = short_protein_change(entry.get("HGVSp", ""))
    if change:
        return change
    aas, pos = entry.get("Amino_acids", ""), entry.get("Protein_position", "")
    if aas and pos:
        parts = aas.split("/")
        return f"{parts[0]}{pos}{parts[1] if len(parts) > 1 else ''}"
    return ""


def sample_values(fields, column):
    return dict(zip(fields[8].split(":"), fields[column].split(":"))) if len(fields) > column else {}


def per_allele(values, key, index, number_r):
    """value for ALT index from a Number=A (or Number=R when number_r) field."""
    raw = values.get(key, ".").split(",")
    i = index + 1 if number_r else index
    return raw[i] if i < len(raw) else "."


def read_vcf(path, sample="tumor"):
    """VEP + readcount VCF -> list of variant dicts (one per ALT allele).

    Each has chrom, pos, ref, alt, filter, csq (the picked entry, or None),
    the tumour DNA counts (dna_depth, dna_alt, dna_vaf from DP/AD/AF) and
    RNA counts (rna_depth, rna_alt, rna_vaf from RDP/RAD/RAF), each None
    when absent. An empty file (a stub) yields no variants.
    """
    variants, csq_format, column = [], None, None
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##"):
                if line.startswith("##INFO=<ID=CSQ") and "Format:" in line:
                    csq_format = line.split("Format:")[1].strip().rstrip('">').split("|")
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                if sample not in fields:
                    raise ValueError(f"{path}: no '{sample}' sample column")
                column = fields.index(sample)
                continue
            if len(fields) < 8 or column is None:
                continue
            chrom, pos, _id, ref, alts, _qual, filt, info = fields[:8]
            values = sample_values(fields, column)
            csq_entries = []
            for item in info.split(";"):
                if item.startswith("CSQ=") and csq_format:
                    csq_entries = [dict(zip(csq_format, e.split("|"))) for e in item[4:].split(",")]
            alt_list = alts.split(",")
            for i, alt in enumerate(alt_list):
                entries = csq_entries
                if len(alt_list) > 1:
                    matched = [e for e in csq_entries if e.get("Allele") in (alt, alt[1:] or "-")]
                    entries = matched or csq_entries
                dna_depth = parse_int(values.get("DP"))
                rna_depth = parse_int(values.get("RDP"))
                dna_vaf = parse_float(per_allele(values, "AF", i, False))
                rna_vaf = parse_float(per_allele(values, "RAF", i, False))
                dna_alt = parse_int(per_allele(values, "AD", i, True))
                rna_alt = parse_int(per_allele(values, "RAD", i, True))
                if dna_alt is None and dna_vaf is not None and dna_depth is not None:
                    dna_alt = int(round(dna_vaf * dna_depth))
                if rna_alt is None and rna_vaf is not None and rna_depth is not None:
                    rna_alt = int(round(rna_vaf * rna_depth))
                variants.append({
                    "chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt, "filter": filt,
                    "csq": pick_csq(entries), "csq_all": entries,
                    "dna_depth": dna_depth, "dna_alt": dna_alt, "dna_vaf": dna_vaf,
                    "rna_depth": rna_depth, "rna_alt": rna_alt, "rna_vaf": rna_vaf,
                })
    return variants


def variant_key(v):
    return (v["chrom"], v["pos"], v["ref"], v["alt"])


def variant_id(v):
    return f"{v['chrom']}-{v['pos']}-{v['ref']}-{v['alt']}"


def match_pvacseq_id(pvacseq_id, index):
    """pVACseq 'chrom-start-stop-ref-alt' -> our variant key, or None.

    start is VCF POS - 1 for an SNV, POS for an insertion or deletion
    (see add_rna_timepoint.py), so both are tried.
    """
    try:
        chrom, start, _stop, ref, alt = pvacseq_id.split("-", 4)
        start = int(start)
    except ValueError:
        return None
    for pos in (start + 1, start):
        key = (chrom, pos, ref, alt)
        if key in index:
            return key
    return None


# ---------------------------------------------------------------------------
# Neoantigens from pVACseq aggregated reports
# ---------------------------------------------------------------------------
def read_binders(path, max_ic50, max_percentile):
    """{pVACseq ID: row} for variants whose best peptide passes either threshold."""
    binders = {}
    if not path:
        return binders
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            ic50 = parse_float(row.get("IC50 MT"))
            pct = parse_float(row.get("%ile MT"))
            if (ic50 is not None and ic50 <= max_ic50) or (pct is not None and pct <= max_percentile):
                binders[row.get("ID", "")] = row
    return binders


# ---------------------------------------------------------------------------
# Categories
# ---------------------------------------------------------------------------
def categories(labels):
    """Ordered category names for these timepoint labels, most-supported first."""
    if len(labels) == 1:
        return [f"Expressed ({labels[0]})", DNA_ONLY]
    if len(labels) == 2:
        return ["Expressed both", f"{labels[1]} only", f"{labels[0]} only", DNA_ONLY]
    names = []
    for size in range(len(labels), 0, -1):
        for combo in combinations(labels, size):
            names.append("Expressed all" if size == len(labels) else "Expressed " + " + ".join(combo))
    return names + [DNA_ONLY]


def category_for(supported, labels):
    """supported: the timepoint labels at which the variant is RNA-supported."""
    chosen = [label for label in labels if label in supported]
    if not chosen:
        return DNA_ONLY
    if len(labels) == 1:
        return f"Expressed ({labels[0]})"
    if len(chosen) == len(labels):
        return "Expressed both" if len(labels) == 2 else "Expressed all"
    if len(labels) == 2:
        return f"{chosen[0]} only"
    return "Expressed " + " + ".join(chosen)


def category_colours(labels):
    names = categories(labels)
    if len(labels) == 2:
        return dict(zip(names, [PALETTE_TWO["both"], PALETTE_TWO["later"], PALETTE_TWO["first"], PALETTE_TWO["none"]]))
    colours = {names[0]: PALETTE_TWO["both"], DNA_ONLY: PALETTE_TWO["none"]}
    for i, name in enumerate(n for n in names if n not in colours):
        colours[name] = EXTRA_COLOURS[i % len(EXTRA_COLOURS)]
    return colours


def rna_status(alt, depth, vaf, min_alt, min_vaf):
    if depth is None or depth == 0:
        return "no_coverage"
    if alt is not None and alt >= min_alt and (vaf if vaf is not None else alt / depth) >= min_vaf:
        return "supported"
    return "covered_no_alt"


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------
VARIANT_SETS = [
    ("non_synonymous", "Non-synonymous (all)"),
    ("non_synonymous_clonal", "Non-synonymous (clonal)"),
    ("neoag_class_i_clonal", "Clonal neoAg MHC I"),
    ("neoag_class_ii_clonal", "Clonal neoAg MHC II"),
    ("neoag_class_i", "All neoAg MHC I"),
    ("neoag_class_ii", "All neoAg MHC II"),
]
PLOTTED_SETS = ["non_synonymous", "non_synonymous_clonal", "neoag_class_i_clonal", "neoag_class_ii_clonal"]


def parse_timepoint(spec):
    if "=" not in spec:
        raise argparse.ArgumentTypeError(f"--timepoint must be LABEL=VCF, got '{spec}'")
    label, path = spec.split("=", 1)
    if not label.strip():
        raise argparse.ArgumentTypeError(f"--timepoint has an empty label: '{spec}'")
    return label.strip(), Path(path)


def analyse(args):
    labels = [label for label, _ in args.timepoint]
    if len(set(labels)) != len(labels):
        raise SystemExit("--timepoint labels must be unique")
    segments = annotate_ccf.read_segments(args.segments) if args.segments else None

    dna = [v for v in read_vcf(args.vcf) if v["filter"] in ("PASS", ".")]
    rna = {}
    for label, path in args.timepoint:
        rna[label] = {variant_key(v): v for v in read_vcf(path)}

    index = {variant_key(v): v for v in dna}
    binders = {"I": {}, "II": {}}
    for mhc, path, ic50, pct in (("I", args.class_i, args.class_i_ic50, args.class_i_percentile),
                                 ("II", args.class_ii, args.class_ii_ic50, args.class_ii_percentile)):
        for pid, row in read_binders(path, ic50, pct).items():
            key = match_pvacseq_id(pid, index)
            if key:
                binders[mhc][key] = row

    rows = []
    for v in dna:
        if not is_non_synonymous(v["csq"]):
            continue
        key = variant_key(v)
        ccf_values = annotate_ccf.ccf_columns(v["dna_vaf"], v["chrom"], v["pos"], args.purity,
                                              args.ploidy, args.sex, segments)
        ccf = parse_float(ccf_values[-1])
        clonal = None if ccf is None else ccf >= args.clonal_ccf
        supported, per_tp = set(), []
        for label in labels:
            r = rna[label].get(key)
            depth, alt, vaf = (r["rna_depth"], r["rna_alt"], r["rna_vaf"]) if r else (None, None, None)
            status = rna_status(alt, depth, vaf, args.min_alt_reads, args.min_rna_vaf)
            if status == "supported":
                supported.add(label)
            per_tp.append((depth, alt, vaf, status))
        csq = v["csq"] or {}
        rows.append({
            "key": key, "id": variant_id(v), "gene": csq.get("SYMBOL", ""),
            "consequence": (csq.get("Consequence") or "").split("&")[0],
            "protein_change": protein_change(csq),
            "dna_depth": v["dna_depth"], "dna_vaf": v["dna_vaf"], "ccf": ccf, "clonal": clonal,
            "per_tp": per_tp, "category": category_for(supported, labels),
            "neoag_I": binders["I"].get(key), "neoag_II": binders["II"].get(key),
        })
    return labels, rows


def set_members(rows, name):
    if name == "non_synonymous":
        return rows
    if name == "non_synonymous_clonal":
        return [r for r in rows if r["clonal"]]
    mhc = "I" if "class_i" in name and "class_ii" not in name else "II"
    members = [r for r in rows if r[f"neoag_{mhc}"]]
    return [r for r in members if r["clonal"]] if name.endswith("_clonal") else members


def summarise(labels, rows):
    """[(set name, set label, {category: n}, total)]"""
    summary = []
    for name, label in VARIANT_SETS:
        members = set_members(rows, name)
        counts = {c: 0 for c in categories(labels)}
        for r in members:
            counts[r["category"]] += 1
        summary.append((name, label, counts, len(members)))
    return summary


def fmt(value, digits=3):
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def write_variants(path, labels, rows):
    header = ["ID", "Gene", "Consequence", "Protein Change", "DNA Depth", "DNA VAF", "CCF", "Clonal"]
    for label in labels:
        header += [f"RNA Depth ({label})", f"RNA Alt Reads ({label})", f"RNA VAF ({label})", f"RNA Status ({label})"]
    header += ["RNA Support", "NeoAg MHC I", "MHC I Best Peptide", "MHC I Allele", "MHC I IC50 MT",
               "NeoAg MHC II", "MHC II Best Peptide", "MHC II Allele", "MHC II IC50 MT"]
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for r in rows:
            line = [r["id"], r["gene"], r["consequence"], r["protein_change"], fmt(r["dna_depth"]),
                    fmt(r["dna_vaf"]), fmt(r["ccf"]), fmt(r["clonal"])]
            for depth, alt, vaf, status in r["per_tp"]:
                line += [fmt(depth), fmt(alt), fmt(vaf), status]
            line.append(r["category"])
            for mhc in ("I", "II"):
                b = r[f"neoag_{mhc}"]
                line += ["yes" if b else "no", (b or {}).get("Best Peptide", "") or "",
                         (b or {}).get("Allele", "") or "", (b or {}).get("IC50 MT", "") or ""]
            writer.writerow(line)


def write_summary(path, labels, summary):
    cats = categories(labels)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Set", "Label", "Total"] + cats + [f"% {c}" for c in cats])
        for name, label, counts, total in summary:
            pct = [f"{100 * counts[c] / total:.1f}" if total else "NA" for c in cats]
            writer.writerow([name, label, total] + [counts[c] for c in cats] + pct)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
def plot(prefix, patient, labels, rows, summary, args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    cats = categories(labels)
    colours = category_colours(labels)
    plotted = [s for s in summary if s[0] in PLOTTED_SETS][::-1]   # top-to-bottom in reading order

    fig = plt.figure(figsize=(12, 4.6))
    grid = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 0.9], wspace=0.35)
    ax_n, ax_pct, ax_sc = fig.add_subplot(grid[0]), fig.add_subplot(grid[1]), fig.add_subplot(grid[2])

    names = [s[1] for s in plotted]
    y = range(len(plotted))
    for ax, as_percent in ((ax_n, False), (ax_pct, True)):
        left = [0.0] * len(plotted)
        for cat in cats:
            widths = []
            for _, _, counts, total in plotted:
                value = counts[cat]
                widths.append((100.0 * value / total if total else 0.0) if as_percent else value)
            ax.barh(list(y), widths, left=left, color=colours[cat], edgecolor="white", linewidth=0.6, height=0.68)
            for i, (w, l0) in enumerate(zip(widths, left)):
                if as_percent and w >= 9:
                    ax.text(l0 + w / 2, i, f"{w:.0f}%", ha="center", va="center", fontsize=7.5,
                            color="white" if cat != DNA_ONLY else "#333333")
            left = [a + b for a, b in zip(left, widths)]
        ax.set_yticks(list(y))
        ax.spines[["top", "right"]].set_visible(False)
        if as_percent:
            ax.set_xlim(0, 100)
            ax.set_xlabel("% of variants in set")
            ax.set_yticklabels([])
        else:
            ax.set_yticklabels(names, fontsize=9)
            ax.set_xlabel("Variants (n)")
            xmax = max([s[3] for s in plotted] + [1])
            ax.set_xlim(0, xmax * 1.15)
            for i, s in enumerate(plotted):
                ax.text(s[3] + xmax * 0.01, i, f"n={s[3]}", va="center", fontsize=8)
    ax_n.set_title("Count", fontsize=10, loc="left")
    ax_pct.set_title("Proportion", fontsize=10, loc="left")

    # Right panel: per-variant RNA VAF, baseline vs later when there are
    # two timepoints, otherwise DNA vs RNA VAF.
    def vaf_or_zero(value):
        return 0.0 if value is None else value
    if len(labels) >= 2:
        xs = [vaf_or_zero(r["per_tp"][0][2]) for r in rows]
        ys = [vaf_or_zero(r["per_tp"][-1][2]) for r in rows]
        ax_sc.set_xlabel(f"RNA VAF ({labels[0]})")
        ax_sc.set_ylabel(f"RNA VAF ({labels[-1]})")
    else:
        xs = [vaf_or_zero(r["dna_vaf"]) for r in rows]
        ys = [vaf_or_zero(r["per_tp"][0][2]) for r in rows]
        ax_sc.set_xlabel("DNA VAF")
        ax_sc.set_ylabel(f"RNA VAF ({labels[0]})")
    for cat in cats[::-1]:
        pts = [(x, yv, r) for x, yv, r in zip(xs, ys, rows) if r["category"] == cat]
        if pts:
            ax_sc.scatter([p[0] for p in pts], [p[1] for p in pts], s=[26 if p[2]["clonal"] else 12 for p in pts],
                          color=colours[cat], edgecolor="#333333", linewidth=0.3, alpha=0.85)
    ax_sc.plot([0, 1], [0, 1], ls=":", color="#999999", lw=0.8)
    ax_sc.set_xlim(-0.02, 1.02)
    ax_sc.set_ylim(-0.02, 1.02)
    ax_sc.spines[["top", "right"]].set_visible(False)
    ax_sc.set_title("Per variant (large = clonal)", fontsize=10, loc="left")

    handles = [Patch(facecolor=colours[c], label=c) for c in cats]
    fig.legend(handles=handles, title="RNA support", loc="lower center", ncol=len(cats), frameon=False,
               fontsize=9, title_fontsize=9, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle(f"{patient}: non-synonymous variants and clonal neoantigens supported in tumour RNA",
                 fontsize=12, x=0.01, ha="left", fontweight="bold")
    fig.text(0.01, 0.905,
             f"RNA support: >= {args.min_alt_reads} alt reads"
             + (f" and VAF >= {args.min_rna_vaf:g}" if args.min_rna_vaf > 0 else "")
             + f"; clonal: CCF >= {args.clonal_ccf:g}; neoAg: best peptide IC50 <= {args.class_i_ic50:g} nM or "
               f"%ile <= {args.class_i_percentile:g} (MHC I), IC50 <= {args.class_ii_ic50:g} nM or "
               f"%ile <= {args.class_ii_percentile:g} (MHC II)",
             fontsize=8, color="#555555", ha="left")
    fig.subplots_adjust(left=0.15, right=0.98, top=0.82, bottom=0.25)
    fig.savefig(f"{prefix}.rna_support.pdf")
    fig.savefig(f"{prefix}.rna_support.png", dpi=200)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--patient", required=True)
    parser.add_argument("--vcf", type=Path, required=True, help="VEP + readcount annotated somatic VCF")
    parser.add_argument("--timepoint", type=parse_timepoint, action="append", required=True,
                        help="LABEL=VCF with that RNA sample's counts in the tumour column; repeat, in time order")
    parser.add_argument("--class-i", type=Path, help="pVACseq MHC class I aggregated report")
    parser.add_argument("--class-ii", type=Path, help="pVACseq MHC class II aggregated report")
    parser.add_argument("--purity", type=float, required=True)
    parser.add_argument("--ploidy", type=float)
    parser.add_argument("--segments", type=Path)
    parser.add_argument("--sex", type=str.lower, choices=["male", "female"])
    parser.add_argument("--min-alt-reads", type=int, default=2)
    parser.add_argument("--min-rna-vaf", type=float, default=0.0)
    parser.add_argument("--clonal-ccf", type=float, default=0.8)
    parser.add_argument("--class-i-ic50", type=float, default=500.0)
    parser.add_argument("--class-i-percentile", type=float, default=2.0)
    parser.add_argument("--class-ii-ic50", type=float, default=1000.0)
    parser.add_argument("--class-ii-percentile", type=float, default=10.0)
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument("-o", "--output-prefix", required=True)
    args = parser.parse_args(argv)
    if not 0 < args.purity <= 1:
        parser.error("--purity must be in (0, 1]")
    if args.min_alt_reads < 1:
        parser.error("--min-alt-reads must be at least 1")

    labels, rows = analyse(args)
    summary = summarise(labels, rows)
    write_variants(f"{args.output_prefix}.rna_support.variants.tsv", labels, rows)
    write_summary(f"{args.output_prefix}.rna_support.summary.tsv", labels, summary)
    if not args.no_plot:
        plot(args.output_prefix, args.patient, labels, rows, summary, args)

    for name, label, counts, total in summary:
        detail = ", ".join(f"{c}: {n}" for c, n in counts.items())
        print(f"{label}: {total} ({detail})")


if __name__ == "__main__":
    main()
