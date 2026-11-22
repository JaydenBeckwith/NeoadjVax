#!/usr/bin/env python3
"""
Per-patient melanoma report (LaTeX) from NeoadjVax's own outputs.

Laid out after the NeoDisc report templates (sample summary, HLA, antigen
presentation machinery, variants, expression, candidates), but filled only
from what this pipeline produces and made melanoma-specific:

  * genomic subtype (TCGA SKCM 2015: BRAF hotspot / RAS hotspot / NF1 /
    triple wild-type) and the status of each panel driver, with BRAF
    V600 / non-V600 class called from the protein change;
  * tumour mutational burden, with the 10 mut/Mb TMB-high line, and a UV
    signature read-out (C>T share, C>T at dipyrimidines, CC>TT doublets);
  * purity, ploidy and CCF-based clonality;
  * HLA class I and II alleles, homozygous loci, and gene-level copy
    number / LOH at HLA-A/B/C and B2M from the gene panel;
  * the presentation-machinery and driver panel flags;
  * RNA support of non-synonymous variants and clonal neoantigens across
    timepoints (rna_support.py's figure and table);
  * fusion neoantigens from the optional STAR-Fusion/AGFusion/pVACfuse
    branch, when it ran;
  * GTEx / TCGA expression context and the top MHC I and II candidates.

NeoDisc sections with no NeoadjVax equivalent (virus detection, mass
spectrometry, TAAs, inflammation score, multi-caller Venn diagrams,
genome-wide copy-number plots, allele-specific HLA LOH) are listed in a
"Not assessed" section rather than left out silently.

Every input except --patient, --vcf and --purity is optional; a missing one
turns its section into an explicit "not available" line. This is a
research summary of computational predictions, not a clinical report, and
says so on its first page.

Usage:
    patient_report.py --patient P1 --vcf P1.somatic.vep.readcounts.vcf.gz --purity 0.62 \\
        [--ploidy 3.1] [--sex female] [--segments P1_segments.txt] [--hla 'HLA-A*02:01,...'] \\
        [--fasta genome.fa] [--target-bed targets.bed | --callable-mb 35] \\
        [--gene-panel P1.gene_panel.tsv] \\
        [--class-i P1.MHC_I...tsv] [--class-ii P1.MHC_II...tsv] \\
        [--rna-support-summary P1.rna_support.summary.tsv] \\
        [--pvacfuse-dir pvacfuse/ --fusion-predictions star-fusion.tsv] \\
        [--rna-support-figure P1.rna_support.pdf] [--timepoint-label Baseline ...] \\
        -o P1.report.tex
"""
import argparse
import csv
import datetime
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import annotate_ccf  # noqa: E402
import gene_panel_report  # noqa: E402
import rna_support as rs  # noqa: E402

TMB_HIGH = 10.0

# Codons recurrently mutated in melanoma (TCGA SKCM 2015; Dankort/Yao BRAF
# classes; COSMIC). Used for labelling only: every mutation in a driver gene
# is listed regardless.
HOTSPOTS = {
    "BRAF": {600, 601, 597, 469, 464, 466, 594, 581, 596},
    "NRAS": {12, 13, 61}, "KRAS": {12, 13, 61}, "HRAS": {12, 13, 61},
    "KIT": {559, 560, 576, 642, 816, 822},
    "GNAQ": {183, 209}, "GNA11": {183, 209},
    "RAC1": {29}, "IDH1": {132}, "MAP2K1": {57, 124, 203},
    "CTNNB1": set(range(32, 46)),
}
BRAF_CLASS = {600: "class I", 601: "class II", 597: "class II", 464: "class II", 469: "class II",
              466: "class III", 594: "class III", 581: "class III", 596: "class III"}

CODING_CONSEQUENCES = rs.NON_SYNONYMOUS + ("splice_acceptor_variant", "splice_donor_variant")

PROTEIN_POSITION_RE = re.compile(r"^[A-Z*]?(\d+)")


# ---------------------------------------------------------------------------
# LaTeX helpers
# ---------------------------------------------------------------------------
TEX_ESCAPES = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_",
               "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}


def tex(value):
    text = "NA" if value is None else str(value)
    return "".join(TEX_ESCAPES.get(ch, ch) for ch in text)


def num(value, digits=2, na="NA"):
    if value is None:
        return na
    if isinstance(value, int):
        return str(value)
    return f"{value:.{digits}f}"


def table(header, rows, spec=None, size=r"\small", caption=None, long=False):
    """booktabs table; long=True uses longtable so it can break across pages."""
    spec = spec or "l" * len(header)
    head = " & ".join(rf"\textbf{{{tex(h)}}}" for h in header) + r" \\"
    body = "\n".join(" & ".join(cells) + r" \\" for cells in rows)
    if long:
        cap = rf"\caption{{{caption}}}\\" if caption else ""
        return (f"{{{size}\n\\begin{{longtable}}{{{spec}}}\n{cap}\n\\toprule\n{head}\n\\midrule\n\\endhead\n"
                f"{body}\n\\bottomrule\n\\end{{longtable}}\n}}\n")
    cap = rf"\caption{{{caption}}}" if caption else ""
    return (f"\\begin{{table}}[H]\n\\centering\n{size}\n\\begin{{tabular}}{{{spec}}}\n\\toprule\n{head}\n"
            f"\\midrule\n{body}\n\\bottomrule\n\\end{{tabular}}\n{cap}\n\\end{{table}}\n")


def flag_colour(flag):
    colour = {"LOSS_OF_FUNCTION": "FAIL", "INHIBITORY_GAIN": "FAIL", "ALTERED": "WARN", "OK": "PASS"}.get(flag)
    return rf"\textcolor{{{colour}}}{{\textbf{{{tex(flag)}}}}}" if colour else tex(flag)


def not_available(what):
    return rf"\emph{{{tex(what)}}}" + "\n\n"


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
def read_tsv(path):
    if not path or not Path(path).is_file() or Path(path).stat().st_size == 0:
        return []
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fai(path):
    """{contig: (length, offset, bases per line, bytes per line)}"""
    index = {}
    with open(path) as handle:
        for line in handle:
            name, length, offset, line_bases, line_bytes = line.split("\t")[:5]
            index[name] = (int(length), int(offset), int(line_bases), int(line_bytes))
    return index


class Fasta:
    """Random access to an uncompressed, .fai-indexed FASTA."""

    def __init__(self, path, fai=None):
        self.path = Path(path)
        self.index = read_fai(fai or f"{path}.fai")
        self.handle = open(self.path, "rb")

    def base(self, chrom, pos):
        """1-based base, upper case, or 'N' when outside the contig."""
        entry = self.index.get(chrom)
        if not entry or pos < 1 or pos > entry[0]:
            return "N"
        length, offset, line_bases, line_bytes = entry
        zero = pos - 1
        self.handle.seek(offset + (zero // line_bases) * line_bytes + zero % line_bases)
        return self.handle.read(1).decode("ascii").upper()

    def close(self):
        self.handle.close()


def bed_size_mb(path):
    """Merged size of a BED (0-based) or Picard interval_list (1-based) in Mb."""
    intervals = {}
    interval_list = str(path).endswith(".interval_list")
    with rs.open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "@", "track", "browser")):
                continue
            fields = line.split("\t")
            start, end = int(fields[1]), int(fields[2])
            if interval_list:
                start -= 1
            intervals.setdefault(fields[0], []).append((start, end))
    total = 0
    for spans in intervals.values():
        spans.sort()
        cur_start, cur_end = spans[0]
        for start, end in spans[1:]:
            if start > cur_end:
                total += cur_end - cur_start
                cur_start, cur_end = start, end
            else:
                cur_end = max(cur_end, end)
        total += cur_end - cur_start
    return total / 1e6


def callable_mb(args):
    """(Mb, how it was derived) or (None, reason)."""
    if args.callable_mb:
        return args.callable_mb, "given (--callable-mb)"
    if args.target_bed:
        return bed_size_mb(args.target_bed), f"capture targets ({Path(args.target_bed).name})"
    fai = args.fai or (f"{args.fasta}.fai" if args.fasta else None)
    if fai and Path(fai).is_file():
        index = read_fai(fai)
        contigs = [c for c in (args.contigs.split(",") if args.contigs else index) if c in index]
        if contigs:
            return sum(index[c][0] for c in contigs) / 1e6, f"genome, {len(contigs)} called contigs (WGS)"
    return None, "no target BED, callable size or FASTA index given"


# ---------------------------------------------------------------------------
# Analyses
# ---------------------------------------------------------------------------
def variant_class(v):
    if len(v["ref"]) == len(v["alt"]) == 1:
        return "SNV"
    if len(v["ref"]) == len(v["alt"]):
        return "MNV"
    return "Insertion" if len(v["alt"]) > len(v["ref"]) else "Deletion"


def uv_signature(variants, fasta):
    """C>T share of SNVs, C>T at dipyrimidines, CC>TT doublets (UV hallmarks)."""
    snvs = [v for v in variants if variant_class(v) == "SNV"]
    ct = [v for v in snvs if (v["ref"], v["alt"]) in (("C", "T"), ("G", "A"))]
    doublets = sum(1 for v in variants if (v["ref"], v["alt"]) in (("CC", "TT"), ("GG", "AA")))
    dipyrimidine = None
    if fasta and ct:
        hits = 0
        for v in ct:
            # pyrimidine strand: 5' neighbour of a C is upstream; of a G (C on
            # the minus strand) it is the base downstream on the plus strand
            if v["ref"] == "C":
                hits += fasta.base(v["chrom"], v["pos"] - 1) in "CT"
            else:
                hits += fasta.base(v["chrom"], v["pos"] + 1) in "GA"
        dipyrimidine = hits / len(ct)
    return {"snvs": len(snvs), "ct": len(ct), "ct_fraction": len(ct) / len(snvs) if snvs else None,
            "dipyrimidine_fraction": dipyrimidine, "cc_tt": doublets}


def codon(change):
    match = PROTEIN_POSITION_RE.match(change or "")
    return int(match.group(1)) if match else None


def driver_mutations(variants, genes, ccfs):
    """[(gene, change, consequence, dna_vaf, ccf, rna_vaf, hotspot label)] for panel drivers."""
    found = []
    for v in variants:
        csq = v["csq"] or {}
        gene = csq.get("SYMBOL", "")
        consequences = (csq.get("Consequence") or "").split("&")
        if gene not in genes or not any(c in CODING_CONSEQUENCES for c in consequences):
            continue
        change = rs.protein_change(csq)
        position = codon(change)
        hotspot = ""
        if position in HOTSPOTS.get(gene, ()) and "missense_variant" in consequences:
            hotspot = "hotspot"
            if gene == "BRAF":
                hotspot = f"hotspot, {BRAF_CLASS.get(position, 'class unknown')}"
        found.append((gene, change or "NA", consequences[0], v["dna_vaf"], ccfs.get(rs.variant_key(v)),
                      v["rna_vaf"], hotspot))
    return sorted(found, key=lambda m: (m[0], m[1]))


def melanoma_subtype(drivers):
    """(TCGA SKCM genomic subtype, BRAF status) from the driver mutation list."""
    def hot(gene):
        return [m for m in drivers if m[0] == gene and m[6].startswith("hotspot")]
    braf_hot = hot("BRAF")
    ras_hot = hot("NRAS") + hot("KRAS") + hot("HRAS")
    nf1 = [m for m in drivers if m[0] == "NF1"]
    braf_any = [m for m in drivers if m[0] == "BRAF"]
    if braf_hot:
        braf_status = "; ".join(f"BRAF {m[1]}, {m[6].split(', ', 1)[-1]}" for m in braf_hot)
    elif braf_any:
        braf_status = "BRAF other (" + ", ".join(m[1] for m in braf_any) + ")"
    else:
        braf_status = "BRAF wild-type"
    if braf_hot:
        subtype = "BRAF hotspot"
    elif ras_hot:
        subtype = "RAS hotspot (" + ", ".join(f"{m[0]} {m[1]}" for m in ras_hot) + ")"
    elif nf1:
        subtype = "NF1 mutant"
    else:
        subtype = "Triple wild-type"
    return subtype, braf_status


def split_hla(alleles):
    class_i, class_ii = [], []
    for allele in [a.strip() for a in (alleles or "").replace("|", ",").replace(";", ",").split(",") if a.strip()]:
        gene = allele.replace("HLA-", "").split("*")[0]
        (class_i if gene in ("A", "B", "C") else class_ii).append(allele)
    return class_i, class_ii


def class_i_zygosity(class_i):
    loci = {}
    for allele in class_i:
        loci.setdefault(allele.replace("HLA-", "").split("*")[0], []).append(allele)
    rows = []
    for locus in ("A", "B", "C"):
        listed = loci.get(locus, [])
        if not listed:
            state = "not typed"
        elif len(set(listed)) == 1:
            state = "homozygous" if len(listed) > 1 else "one allele listed (homozygous, or deduplicated)"
        else:
            state = "heterozygous"
        rows.append((f"HLA-{locus}", ", ".join(listed) or "-", state))
    return rows


def binder(row, ic50_max, pct_max):
    ic50, pct = rs.parse_float(row.get("IC50 MT")), rs.parse_float(row.get("%ile MT"))
    return (ic50 is not None and ic50 <= ic50_max) or (pct is not None and pct <= pct_max)


def rank_candidates(rows, ic50_max, pct_max, clonal_ccf, min_rna_vaf_seen):
    """Binders, best first: clonal, seen in baseline RNA, expression OK, then IC50."""
    def key(row):
        ccf = rs.parse_float(row.get("CCF"))
        rna_vaf = rs.parse_float(row.get("RNA VAF"))
        flag = row.get("Expression Flag", "OK")
        ic50 = rs.parse_float(row.get("IC50 MT"))
        return (not (ccf is not None and ccf >= clonal_ccf),
                not (rna_vaf is not None and rna_vaf > min_rna_vaf_seen),
                flag not in ("OK", ""),
                ic50 if ic50 is not None else 1e9)
    return sorted((r for r in rows if binder(r, ic50_max, pct_max)), key=key)


def histogram_rows(values, edges):
    counts = [0] * (len(edges) - 1)
    for value in values:
        for i in range(len(edges) - 1):
            if edges[i] <= value < edges[i + 1] or (i == len(edges) - 2 and value == edges[-1]):
                counts[i] += 1
                break
    peak = max(counts + [1])
    return [(f"{edges[i]:.1f}--{edges[i + 1]:.1f}", counts[i],
             rf"\textcolor{{ceruleanblue}}{{\rule{{{6.0 * counts[i] / peak:.2f}cm}}{{1.6ex}}}}")
            for i in range(len(counts))]


# Fusions recurrently reported in melanoma, and what they mean. Kinase
# fusions are the actionable ones and cluster in Spitz/spitzoid tumours;
# a BRAF fusion is also an acquired route to MAPK-inhibitor resistance.
# Used for labelling only: every fusion with a binder is listed.
FUSION_PARTNERS = {
    "ALK": "kinase fusion (Spitz)", "ROS1": "kinase fusion (Spitz)",
    "NTRK1": "kinase fusion (Spitz)", "NTRK2": "kinase fusion (Spitz)",
    "NTRK3": "kinase fusion (Spitz)", "RET": "kinase fusion (Spitz)",
    "MET": "kinase fusion (Spitz)", "BRAF": "BRAF fusion (MAPK, resistance route)",
    "MAP3K8": "MAP3K8 fusion (MAPK)", "PRKCA": "PRKCA fusion (Spitz)",
    "AXL": "kinase fusion", "FGFR1": "kinase fusion", "PDGFRB": "kinase fusion",
}


def fusion_partners(gene):
    """pVACfuse 'A_B' or 'A-B' gene field -> [partner, ...]."""
    return [p for p in re.split(r"[_\-]", (gene or "").strip()) if p]


def fusion_note(gene):
    for partner in fusion_partners(gene):
        if partner in FUSION_PARTNERS:
            return FUSION_PARTNERS[partner]
    return ""


def read_pvacfuse(directory):
    """pvacfuse/ -> (status, {"I": rows, "II": rows}).

    status is one of:
      not_run        the branch did not run for this patient
      no_proteins    AGFusion produced no coding fusion protein (the module
                     writes no_fusion_proteins.txt and skips pVACfuse)
      ok             pVACfuse wrote its per-class reports
    """
    reports = {"I": [], "II": []}
    if not directory or not Path(directory).is_dir():
        return "not_run", reports
    directory = Path(directory)
    if (directory / "no_fusion_proteins.txt").is_file():
        return "no_proteins", reports
    found = False
    for mhc in ("I", "II"):
        for path in sorted((directory / f"MHC_Class_{mhc}").glob("*.all_epitopes.aggregated.tsv")):
            reports[mhc] = read_tsv(path)
            found = True
            break
    return ("ok" if found else "no_proteins"), reports


def count_fusions(path):
    """Fusions in a STAR-Fusion predictions table, or None."""
    rows = read_tsv(path)
    if not rows:
        return None
    key = next((c for c in rows[0] if c.lstrip("#").strip() in ("FusionName", "fusion_name")), None)
    return len({(r.get(key) or "").strip() for r in rows if (r.get(key) or "").strip()}) if key else len(rows)


def rank_fusions(rows, ic50_max, pct_max, min_reads, min_expr):
    """pVACfuse binders, best first: read support, then expression, then IC50."""
    def key(row):
        reads = rs.parse_float(row.get("Read Support"))
        expr = rs.parse_float(row.get("Expr"))
        ic50 = rs.parse_float(row.get("IC50 MT"))
        return (not (reads is not None and reads >= min_reads),
                not (expr is not None and expr >= min_expr),
                ic50 if ic50 is not None else 1e9)
    return sorted((r for r in rows if binder(r, ic50_max, pct_max)), key=key)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
PREAMBLE = r"""\documentclass[10pt]{article}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage[utf8]{inputenc}
\renewcommand{\familydefault}{\sfdefault}
\usepackage[margin=0.6in]{geometry}
\usepackage{graphicx}
\usepackage[table]{xcolor}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{array}
\usepackage{float}
\usepackage{hyperref}
\definecolor{PASS}{RGB}{0,153,0}
\definecolor{WARN}{RGB}{204,102,0}
\definecolor{FAIL}{RGB}{204,0,0}
\definecolor{ceruleanblue}{rgb}{0.16, 0.32, 0.75}
\definecolor{Gray}{RGB}{224,224,224}
\definecolor{teagreen}{RGB}{226,239,218}
\hypersetup{colorlinks=true, linkcolor=blue, urlcolor=ceruleanblue}
\setlength{\parindent}{0pt}
\setlength{\parskip}{4pt}
"""


def build_report(args):
    # --vcf is absent when the run did not include somatic calling (see
    # --run in docs/RUNNING.md). Everything derived from variants is then
    # reported as not assessed rather than as zero, which would read as a
    # tumour with no mutations.
    somatic_assessed = args.vcf is not None
    variants = [v for v in rs.read_vcf(args.vcf) if v["filter"] in ("PASS", ".")] if somatic_assessed else []
    segments = annotate_ccf.read_segments(args.segments) if args.segments else None
    ccfs = {}
    for v in variants:
        ccfs[rs.variant_key(v)] = rs.parse_float(annotate_ccf.ccf_columns(
            v["dna_vaf"], v["chrom"], v["pos"], args.purity, args.ploidy, args.sex, segments)[-1])
    non_syn = [v for v in variants if rs.is_non_synonymous(v["csq"])]
    fasta = Fasta(args.fasta, args.fai) if args.fasta else None
    uv = uv_signature(variants, fasta) if somatic_assessed else None
    if fasta:
        fasta.close()
    size_mb, size_basis = callable_mb(args)
    tmb = len(variants) / size_mb if (size_mb and somatic_assessed) else None
    tmb_ns = len(non_syn) / size_mb if (size_mb and somatic_assessed) else None

    panel_rows = read_tsv(args.gene_panel)
    panel_by_gene = {r["Gene"]: r for r in panel_rows}
    driver_genes = {g for p, genes in gene_panel_report.DEFAULT_PANEL.items() if p.startswith("driver_") for g in genes}
    driver_genes |= {r["Gene"] for r in panel_rows if r.get("Pathway", "").startswith("driver_")}
    drivers = driver_mutations(variants, driver_genes, ccfs)
    subtype, braf_status = melanoma_subtype(drivers)

    class_i, class_ii = split_hla(args.hla)
    class_i_rows = read_tsv(args.class_i)
    class_ii_rows = read_tsv(args.class_ii)
    rna_summary = read_tsv(args.rna_support_summary)
    fusion_status, fusion_reports = read_pvacfuse(args.pvacfuse_dir)
    timepoints = args.timepoint_label or []

    ns_ccf = [ccfs[rs.variant_key(v)] for v in non_syn if ccfs.get(rs.variant_key(v)) is not None]
    clonal_ns = sum(1 for c in ns_ccf if c >= args.clonal_ccf)

    out = [PREAMBLE, r"\begin{document}", ""]
    out.append(rf"\begin{{center}}{{\LARGE NeoadjVax patient report: \textbf{{{tex(args.patient)}}}}}\\[4pt]"
               rf"{{\large Neoantigen nomination for neoadjuvant melanoma vaccine research}}\\[2pt]"
               rf"NeoadjVax {tex(args.pipeline_version)} \quad run {tex(args.run_date)}\end{{center}}")
    out.append(r"\fcolorbox{WARN}{white}{\parbox{0.96\textwidth}{\small\textbf{Research use only.} "
               r"Every value below is a computational prediction from sequencing data, not a validated or "
               r"clinically accredited result. Binding predictions, clonality and expression calls each carry "
               r"error; treat candidates as hypotheses for experimental testing.}}")

    # 1. Summary -------------------------------------------------------------
    out.append(r"\section{Patient summary}")
    rna_text = ", ".join(timepoints) if timepoints else "none"
    hla_loh_genes = [g for g in ("HLA-A", "HLA-B", "HLA-C", "B2M") if "LOH" in panel_by_gene.get(g, {}).get("CN Call", "")
                     or "LOSS" in panel_by_gene.get(g, {}).get("CN Call", "")]
    summary_rows = [
        ("Patient", tex(args.patient)),
        ("Project", tex(args.project or "NA")),
        ("Diagnosis", tex(args.diagnosis)),
        ("Matched TCGA cohort / GTEx tissue", tex(f"{args.tcga_cancer_type or 'NA'} / {args.gtex_tissue or 'NA'}")),
        ("Tumour RNA timepoints", tex(rna_text)),
        ("Purity (Sequenza)", num(args.purity) + (r" \textcolor{WARN}{(below 0.3: CCF and CN unreliable)}" if args.purity < 0.3 else "")),
        ("Ploidy", num(args.ploidy) if args.ploidy else "NA"),
        ("Sex", tex(args.sex or "NA")),
        ("Genomic subtype (TCGA SKCM)", rf"\textbf{{{tex(subtype)}}}"),
        ("BRAF status", tex(braf_status)),
        ("PASS somatic variants / non-synonymous",
         f"{len(variants)} / {len(non_syn)}" if somatic_assessed else "NA (somatic calling not run)"),
        ("TMB (all PASS, mut/Mb)", (rf"\textbf{{{num(tmb)}}}" + (r" \textcolor{FAIL}{TMB-high}" if tmb >= TMB_HIGH else ""))
         if tmb is not None else "NA"),
        ("Clonal non-synonymous (CCF $\\geq$ " + f"{args.clonal_ccf:g})", f"{clonal_ns} of {len(ns_ccf)}" if ns_ccf else "NA"),
        ("UV signature (C>T of SNVs)",
         f"{100 * uv['ct_fraction']:.0f}\\%" if uv and uv["ct_fraction"] is not None else "NA"),
        ("HLA class I alleles", tex(", ".join(class_i) or "NA")),
        ("Gene-level HLA / B2M copy loss or LOH", tex(", ".join(hla_loh_genes) if hla_loh_genes else
                                                      ("none" if panel_rows else "NA (no gene panel)"))),
        ("Presentation-machinery adverse flags", tex(str(sum(1 for r in panel_rows if r.get("Panel Flag") in
                                                            ("LOSS_OF_FUNCTION", "INHIBITORY_GAIN")
                                                            and not r.get("Pathway", "").startswith("driver_"))))
         if panel_rows else "NA"),
        ("Fusion neoantigens MHC I / MHC II",
         (f"{sum(1 for r in fusion_reports['I'] if binder(r, args.class_i_ic50, args.class_i_percentile))} / "
          f"{sum(1 for r in fusion_reports['II'] if binder(r, args.class_ii_ic50, args.class_ii_percentile))}")
         if fusion_status == "ok" else
         ("none (no coding fusion protein)" if fusion_status == "no_proteins" else "NA (fusion branch not run)")),
        ("Predicted neoantigens MHC I / MHC II (variants)",
         f"{sum(1 for r in class_i_rows if binder(r, args.class_i_ic50, args.class_i_percentile))} / "
         f"{sum(1 for r in class_ii_rows if binder(r, args.class_ii_ic50, args.class_ii_percentile))}"
         if (class_i_rows or class_ii_rows) else "NA"),
    ]
    out.append(table(["Item", "Value"], summary_rows, spec="p{0.38\\textwidth}p{0.55\\textwidth}"))

    out.append(r"\section{Aim of the analysis}")
    out.append("Nominate somatic neoantigen candidates for a personalised neoadjuvant melanoma vaccine, and "
               "record the tumour features that decide whether those candidates can be presented and seen: "
               "driver status, mutational burden, clonality, HLA and presentation-machinery integrity, and "
               "whether each mutation is transcribed before and on treatment.\n")

    # 2. Drivers --------------------------------------------------------------
    out.append(r"\section{Melanoma drivers and genomic subtype}")
    out.append(rf"Genomic subtype: \textbf{{{tex(subtype)}}}. Subtypes follow TCGA SKCM (Cell 2015): a BRAF "
               r"hotspot takes precedence over a RAS hotspot (NRAS/KRAS/HRAS codons 12, 13, 61), then any NF1 "
               r"non-synonymous mutation; otherwise triple wild-type. BRAF class (I: V600; II: K601, L597, "
               r"G464, G469; III: D594, G466, N581, G596) is from the protein position. TERT promoter mutations "
               r"are non-coding and usually absent from an exome call set, so a blank TERT row is not evidence "
               r"of a wild-type promoter." + "\n")
    if drivers:
        rows = []
        for gene, change, consequence, dna_vaf, ccf, rna_vaf, hotspot in drivers:
            p = panel_by_gene.get(gene, {})
            rows.append([rf"\textbf{{{tex(gene)}}}", tex(change), tex(consequence.replace("_variant", "")),
                         tex(hotspot or "-"), num(dna_vaf), num(ccf), num(rna_vaf),
                         tex(p.get("CN Call", "NA")), tex(p.get("GTEx Call", "NA"))])
        out.append(table(["Gene", "Change", "Consequence", "Hotspot", "DNA VAF", "CCF", "RNA VAF", "CN", "Expr vs GTEx"],
                         rows, spec="lllllrrrll", size=r"\footnotesize", caption="Coding mutations in melanoma driver genes"))
    else:
        out.append(not_available("No coding mutation in the panel's driver genes."))
    cn_drivers = [r for r in panel_rows if r.get("Pathway", "").startswith("driver_")
                  and r.get("CN Call") not in (None, "", "NA", "NEUTRAL")]
    if cn_drivers:
        out.append("Driver genes with a copy-number call (gene-level, from Sequenza segments): " +
                   ", ".join(f"{tex(r['Gene'])} ({tex(r['CN Call'])})" for r in cn_drivers) + ".\n")

    # 3. HLA ------------------------------------------------------------------
    out.append(r"\section{HLA typing and HLA loss}")
    if class_i or class_ii:
        out.append(table(["Locus", "Alleles", "Zygosity"], [[tex(a), tex(b), tex(c)] for a, b, c in class_i_zygosity(class_i)],
                         spec="lll", caption="HLA class I (HLA-HD on the normal, or the samplesheet)"))
        out.append("HLA class II: " + (tex(", ".join(class_ii)) if class_ii else "none typed") +
                   ". A homozygous class I locus halves the peptide repertoire that locus can present "
                   "(lower HLA evolutionary divergence).\n")
    else:
        out.append(not_available("No HLA alleles were passed to the report."))
    hla_genes = [g for g in ("HLA-A", "HLA-B", "HLA-C", "B2M") if g in panel_by_gene]
    if hla_genes:
        out.append(table(["Gene", "Tumour CN", "Minor CN", "CN call", "Expr vs GTEx", "Mutations", "Flag"],
                         [[tex(g), tex(panel_by_gene[g].get("Tumor CN")), tex(panel_by_gene[g].get("Minor CN")),
                           tex(panel_by_gene[g].get("CN Call")), tex(panel_by_gene[g].get("GTEx Call")),
                           tex(panel_by_gene[g].get("Mutations")), flag_colour(panel_by_gene[g].get("Panel Flag"))]
                          for g in hla_genes], spec="lrrllll", size=r"\footnotesize",
                         caption="Gene-level copy number and expression at the class I loci and B2M"))
    out.append(r"\textbf{Allele-specific HLA LOH was not assessed.} The rows above are gene-level copy number "
               r"from Sequenza segments; short reads mismap across HLA alleles, so a minor copy number of 0 here "
               r"is a prompt to run an allele-specific caller (LOHHLA, SpecHLA, DASH), not a call. Candidates "
               r"restricted to a lost allele would not be presented." + "\n")

    # 4. APM panel ------------------------------------------------------------
    out.append(r"\section{Antigen presentation machinery and interferon signalling}")
    if panel_rows:
        apm = [r for r in panel_rows if not r.get("Pathway", "").startswith("driver_")]
        pathways = []
        for r in apm:
            if r["Pathway"] not in pathways:
                pathways.append(r["Pathway"])
        rows = []
        for pathway in pathways:
            members = [r for r in apm if r["Pathway"] == pathway]
            adverse = [r["Gene"] for r in members if r.get("Panel Flag") in ("LOSS_OF_FUNCTION", "INHIBITORY_GAIN")]
            altered = [r["Gene"] for r in members if r.get("Panel Flag") == "ALTERED"]
            rows.append([tex(pathway.replace("_", " ")), str(len(members)),
                         rf"\textcolor{{FAIL}}{{{tex(', '.join(adverse))}}}" if adverse else "-",
                         tex(", ".join(altered)) or "-"])
        out.append(table(["Pathway", "Genes", "Adverse (loss of function / inhibitory gain)", "Altered"], rows,
                         spec="lrp{0.38\\textwidth}p{0.25\\textwidth}", size=r"\small",
                         caption="Presentation-machinery panel, by pathway (full table in the appendix)"))
        out.append("Flags are factual, not verdicts: expression is called against the matched GTEx tissue, copy "
                   "number against sample ploidy, and a truncating mutation counts as loss. HLA-E and HLA-G are "
                   "gain-polarity (inhibitory ligands). Low purity pushes expression towards DOWN.\n")
    else:
        out.append(not_available("Gene panel not run (it needs --gtex_db or --tcga_db)."))

    # 5. Variants, TMB, clonality --------------------------------------------
    out.append(r"\section{Somatic variants, mutational burden and clonality}")
    if not somatic_assessed:
        out.append(not_available(
            "Somatic calling was not part of this run, so variant counts, TMB, the UV signature and "
            "clonality were not assessed. Rerun including the somatic analysis to fill this section."))
    else:
        classes = {}
        for v in variants:
            classes[variant_class(v)] = classes.get(variant_class(v), 0) + 1
        consequences = {}
        for v in variants:
            c = ((v["csq"] or {}).get("Consequence") or "no_annotation").split("&")[0]
            consequences[c] = consequences.get(c, 0) + 1
        out.append(table(["Variant class", "n", "%"],
                         [[tex(k), str(n), f"{100 * n / len(variants):.1f}"] for k, n in sorted(classes.items())] or [["none", "0", "-"]],
                         spec="lrr", caption="PASS somatic variants by class"))
        top_consequences = sorted(consequences.items(), key=lambda kv: -kv[1])[:12]
        out.append(table(["Consequence (canonical transcript)", "n"], [[tex(k), str(n)] for k, n in top_consequences] or [["none", "0"]],
                         spec="lr", size=r"\small", caption="Most frequent consequences"))
        tmb_line = (f"TMB is {num(tmb)} mut/Mb over all PASS variants and {num(tmb_ns)} non-synonymous mut/Mb, "
                    f"over {num(size_mb, 1)} Mb ({tex(size_basis)}). ") if tmb is not None else \
            f"TMB not computed: {tex(size_basis)}. "
        out.append(tmb_line + f"The usual TMB-high line is {TMB_HIGH:g} mut/Mb; panel and exome TMB differ, so "
                   "compare within one assay.\n")
        uv_bits = [f"{uv['ct']} of {uv['snvs']} SNVs are C>T"]
        if uv["dipyrimidine_fraction"] is not None:
            uv_bits.append(f"{100 * uv['dipyrimidine_fraction']:.0f}\\% of those at a dipyrimidine")
        uv_bits.append(f"{uv['cc_tt']} CC>TT doublets")
        out.append("UV signature: " + ", ".join(uv_bits) + ". A dominant C>T share at dipyrimidines with CC>TT "
                   "doublets marks a UV-driven (cutaneous) melanoma; acral and mucosal melanomas lack it and carry "
                   "fewer mutations.\n")
        if ns_ccf:
            out.append(rf"Cancer cell fraction of the non-synonymous variants (purity {num(args.purity)}, "
                       rf"ploidy {num(args.ploidy) if args.ploidy else 'NA'}, "
                       rf"{'Sequenza segments' if segments else 'no segments'}): {clonal_ns} of {len(ns_ccf)} clonal "
                       rf"at CCF $\geq$ {args.clonal_ccf:g}. Point estimates without confidence intervals." + "\n")
            out.append(table(["CCF", "n", ""], [[a, str(b), c] for a, b, c in
                                               histogram_rows(ns_ccf, [0, 0.2, 0.4, 0.6, 0.8, 1.0])],
                             spec="lrl", size=r"\small", caption="CCF distribution, non-synonymous variants"))

    # 6. RNA support ----------------------------------------------------------
    out.append(r"\section{RNA support across timepoints}")
    if rna_summary:
        if args.rna_support_figure:
            out.append(r"\begin{figure}[H]\centering" + "\n" +
                       rf"\includegraphics[width=\textwidth]{{{Path(args.rna_support_figure).name}}}" + "\n" +
                       r"\caption{Non-synonymous variants and clonal neoantigens seen in tumour RNA at each "
                       r"timepoint. Left: counts; middle: proportions; right: per-variant RNA VAF.}" + "\n" +
                       r"\end{figure}")
        cats = [c for c in rna_summary[0] if c not in ("Set", "Label", "Total") and not c.startswith("% ")]
        rows = []
        for r in rna_summary:
            rows.append([tex(r["Label"]), tex(r["Total"])] +
                        [f"{tex(r[c])} ({tex(r['% ' + c])}\\%)" if r["Total"] not in ("0", "") else tex(r[c]) for c in cats])
        out.append(table(["Set", "n"] + cats, rows, spec="lr" + "r" * len(cats), size=r"\footnotesize",
                         caption="RNA support by variant set"))
        out.append(rf"Supported means at least {args.min_alt_reads} RNA reads carry the alternate allele. "
                   r"``DNA only'' includes genes that are silent and sites with too little RNA depth to see "
                   r"the allele; the per-variant table keeps the depth to tell them apart. A drop from baseline "
                   r"to week 6 can reflect tumour regression or immune infiltration diluting tumour RNA, not "
                   r"only loss of expression." + "\n")
    else:
        out.append(not_available("No tumour RNA timepoint for this patient, so RNA support was not assessed."))

    # 7. Expression context ---------------------------------------------------
    out.append(r"\section{Expression context (GTEx and TCGA)}")
    annotated = [r for r in class_i_rows + class_ii_rows if "GTEx Specificity" in r or "Expression Flag" in r]
    if annotated:
        def tally(column):
            counts = {}
            for r in annotated:
                if binder(r, max(args.class_i_ic50, args.class_ii_ic50), max(args.class_i_percentile, args.class_ii_percentile)):
                    counts[r.get(column, "NA") or "NA"] = counts.get(r.get(column, "NA") or "NA", 0) + 1
            return counts
        spec_counts, flag_counts = tally("GTEx Specificity"), tally("Expression Flag")
        if spec_counts:
            out.append(table(["GTEx specificity of candidate genes", "n"], [[tex(k), str(v)] for k, v in sorted(spec_counts.items())],
                             spec="lr", size=r"\small"))
        if flag_counts:
            out.append(table(["Expression flag", "n"], [[tex(k), str(v)] for k, v in sorted(flag_counts.items())],
                             spec="lr", size=r"\small"))
        out.append(rf"Normal tissue: GTEx {tex(args.gtex_tissue or 'NA')}; cohort: {tex(args.tcga_cancer_type or 'NA')}. "
                   r"Candidates in broadly expressed genes are more likely tolerised; testis-restricted genes "
                   r"(cancer-testis antigens) are kept visible. The flag uses the patient's baseline RNA when "
                   r"present, otherwise the TCGA cohort median." + "\n")
    else:
        out.append(not_available("Candidate reports carry no expression columns (GTEx/TCGA tables not given)."))

    # 8. Candidates -----------------------------------------------------------
    out.append(r"\section{Top neoantigen candidates}")
    out.append(rf"Ranked binders: clonal (CCF $\geq$ {args.clonal_ccf:g}) first, then seen in baseline RNA, "
               r"then expression flag OK, then best-peptide IC50. Binder thresholds: "
               rf"MHC I IC50 $\leq$ {args.class_i_ic50:g} nM or percentile $\leq$ {args.class_i_percentile:g}; "
               rf"MHC II IC50 $\leq$ {args.class_ii_ic50:g} nM or percentile $\leq$ {args.class_ii_percentile:g}. "
               r"pVACseq's own Tier is shown for comparison." + "\n")
    week6_col = next((c for c in (class_i_rows[0] if class_i_rows else (class_ii_rows[0] if class_ii_rows else {}))
                      if c.startswith("RNA VAF (")), None)
    for label, rows_in, ic50, pct in (("MHC class I", class_i_rows, args.class_i_ic50, args.class_i_percentile),
                                      ("MHC class II", class_ii_rows, args.class_ii_ic50, args.class_ii_percentile)):
        out.append(rf"\subsection{{{label}}}")
        ranked = rank_candidates(rows_in, ic50, pct, args.clonal_ccf, 0.0)[:args.top_n]
        if not ranked:
            out.append(not_available(f"No {label} binders" + ("" if rows_in else " (report not available)") + "."))
            continue
        first = timepoints[0] if timepoints else "RNA"
        header = ["Gene", "Change", "Best peptide", "Allele", "IC50", "%ile", "DNA VAF", "CCF",
                  f"RNA VAF {first}" if timepoints else "RNA VAF"]
        if week6_col:
            header.append(week6_col.replace("(", "").replace(")", ""))
        header += ["Expr", "Tier"]
        rows = []
        for r in ranked:
            cells = [tex(r.get("Gene")), tex(r.get("AA Change")), r"\texttt{" + tex(r.get("Best Peptide")) + "}",
                     tex((r.get("Allele") or "").replace("HLA-", "")), num(rs.parse_float(r.get("IC50 MT")), 0),
                     num(rs.parse_float(r.get("%ile MT")), 2), num(rs.parse_float(r.get("DNA VAF"))),
                     num(rs.parse_float(r.get("CCF"))), num(rs.parse_float(r.get("RNA VAF")))]
            if week6_col:
                cells.append(num(rs.parse_float(r.get(week6_col))))
            cells += [tex(r.get("Expression Flag", "NA")), tex(r.get("Tier", "NA"))]
            rows.append(cells)
        out.append(table(header, rows,
                         spec="ll>{\\raggedright\\arraybackslash}p{0.17\\textwidth}l" + "r" * (len(header) - 6) + "ll",
                         size=r"\scriptsize", long=True,
                         caption=f"Top {len(ranked)} {label} candidates of "
                                 f"{sum(1 for r in rows_in if binder(r, ic50, pct))} binders"))

    # 9. Fusion neoantigens ---------------------------------------------------
    out.append(r"\section{Fusion neoantigens}")
    n_called = count_fusions(args.fusion_predictions)
    if fusion_status == "not_run":
        out.append(not_available("The fusion branch did not run for this patient (--fusion false, or no tumour "
                                 "RNA reads and no STAR-Fusion table in the samplesheet), so fusion neoantigens "
                                 "were not assessed."))
    else:
        out.append("Gene fusions are an uncommon source of neoantigens in cutaneous melanoma, which is driven by "
                   "a high UV point-mutation load rather than by structural rearrangement, so an empty table here "
                   "is the expected result. They matter when they appear: kinase fusions (ALK, ROS1, NTRK1/2/3, "
                   "RET, MET, BRAF) are recurrent in Spitz and spitzoid tumours, which carry few point mutations "
                   "and so few conventional candidates, and a BRAF fusion is also an acquired route to "
                   "MAPK-inhibitor resistance. A fusion junction peptide is wholly tumour-specific, so it has no "
                   "wild-type counterpart to tolerise against.\n")
        if n_called is not None:
            out.append(f"STAR-Fusion called {n_called} fusion{'' if n_called == 1 else 's'} for this patient. ")
        if fusion_status == "no_proteins":
            out.append("None of them gave AGFusion a coding partner pair, so there was no fusion protein for "
                       "pVACfuse to score. This is not a failure: most fusion calls are intergenic, out of frame "
                       "or non-coding on at least one side.\n")
        else:
            counts = {mhc: sum(1 for r in fusion_reports[mhc]
                               if binder(r, getattr(args, f"class_{mhc.lower()}_ic50"),
                                         getattr(args, f"class_{mhc.lower()}_percentile")))
                      for mhc in ("I", "II")}
            out.append(f"pVACfuse scored {len(fusion_reports['I'])} class I and {len(fusion_reports['II'])} class II "
                       f"fusion epitopes, of which {counts['I']} and {counts['II']} pass the binder thresholds "
                       "used elsewhere in this report.\n")
            for mhc, label in (("I", "MHC class I"), ("II", "MHC class II")):
                ranked = rank_fusions(fusion_reports[mhc],
                                      getattr(args, f"class_{mhc.lower()}_ic50"),
                                      getattr(args, f"class_{mhc.lower()}_percentile"),
                                      args.fusion_min_reads, args.fusion_min_expr)[:args.top_n]
                if not ranked:
                    out.append(not_available(f"No {label} fusion binders."))
                    continue
                rows = []
                for r in ranked:
                    gene = r.get("Gene", "")
                    rows.append([tex(gene.replace("_", "--")), r"\texttt{" + tex(r.get("Best Peptide")) + "}",
                                 tex((r.get("Allele") or "").replace("HLA-", "")),
                                 num(rs.parse_float(r.get("IC50 MT")), 0),
                                 num(rs.parse_float(r.get("%ile MT")), 2),
                                 num(rs.parse_float(r.get("Read Support")), 0),
                                 num(rs.parse_float(r.get("Expr")), 2),
                                 tex(r.get("Tier", "NA")), tex(fusion_note(gene) or "-")])
                out.append(table(["Fusion", "Best peptide", "Allele", "IC50", "%ile", "Reads", "Expr", "Tier", "Note"],
                                 rows, spec="l>{\\raggedright\\arraybackslash}p{0.16\\textwidth}lrrrrlp{0.2\\textwidth}",
                                 size=r"\scriptsize", long=True,
                                 caption=f"{label} fusion candidates, best read and expression support first"))
            out.append(rf"Reads is pVACfuse's junction plus spanning read support and Expr its FFPM, both from the "
                       rf"STAR-Fusion table; pVACfuse's own tiering flags read support below {args.fusion_min_reads:g} "
                       rf"and expression below {args.fusion_min_expr:g}. A fusion candidate rests on the junction "
                       r"being called correctly and in frame, which is a harder call than a point mutation: treat a "
                       r"low-read fusion as a lead to confirm (IGV at the breakpoint, or RT-PCR) rather than as a "
                       r"candidate. AGFusion's Ensembl release can differ from the VEP cache used for the SNV and "
                       r"indel candidates, so the two halves of the report may not use identical transcript models." + "\n")

    # 10. Not assessed --------------------------------------------------------
    out.append(r"\section{Not assessed by NeoadjVax}")
    out.append("These sections of the NeoDisc report have no NeoadjVax equivalent yet and are left out rather "
               "than filled with placeholders:\n")
    out.append(r"\begin{itemize}\setlength{\itemsep}{0pt}"
               r"\item Allele-specific HLA LOH (LOHHLA / SpecHLA / DASH)."
               r"\item Genome-wide allele-specific copy-number and HLA copy-number figures."
               r"\item Multi-caller overlap (NeoadjVax calls with Mutect2 only)."
               r"\item Tumour-associated antigens (TAAs) and cancer-testis antigen expression lists."
               r"\item Inflammation / T-cell infiltration score."
               r"\item Virus detection from RNA (not relevant to melanoma)."
               r"\item Splicing, intron-retention and ERV/transposable-element neoantigens."
               r"\item Mass-spectrometry immunopeptidomics."
               r"\item Sequencing QC summaries (FastQC reports are published under \texttt{fastqc/} for FASTQ inputs)."
               r"\end{itemize}")

    # 11. Methods -------------------------------------------------------------
    out.append(r"\section{Methods in brief}")
    out.append(r"\begin{itemize}\setlength{\itemsep}{0pt}"
               r"\item Tumour/normal DNA: BWA-MEM, MarkDuplicates, BQSR; Mutect2 with germline resource, "
               r"orientation model and PASS filter; VEP (canonical and MANE flags)."
               r"\item HLA: HLA-HD on the normal, or alleles from the samplesheet."
               r"\item Read counts: bam-readcount for tumour/normal DNA and each RNA BAM at the somatic sites."
               r"\item Neoantigens: pVACseq with WhatsHap-phased proximal variants and Sequenza purity."
               r"\item CCF: $u = \mathrm{VAF}(p\,\mathrm{CN}_t + \mathrm{CN}_n(1-p))/p$, multiplicity "
               r"$m = \max(1, \mathrm{round}(u))$ capped at the major copy number, CCF $= \min(1, u/m)$."
               r"\item Expression: featureCounts TPM; GTEx and TCGA reference quantiles."
               r"\item Fusions (optional): STAR-Fusion on the baseline tumour RNA, AGFusion, pVACfuse on the same HLA."
               r"\end{itemize}")

    # Appendix: full panel ----------------------------------------------------
    if panel_rows:
        out.append(r"\clearpage\appendix\section{Full gene panel}")
        rows = [[tex(r.get("Pathway", "").replace("_", " ")), tex(r.get("Gene")), tex(r.get("Sample TPM")),
                 tex(r.get("GTEx Call")), tex(r.get("TCGA Call")), tex(r.get("CN Call")),
                 r"\raggedright " + tex(r.get("Mutations")), flag_colour(r.get("Panel Flag"))] for r in panel_rows]
        out.append(table(["Pathway", "Gene", "TPM", "GTEx", "TCGA", "CN", "Mutations", "Flag"], rows,
                         spec="llrlllp{0.22\\textwidth}l", size=r"\scriptsize", long=True))

    out.append(r"\end{document}")
    return "\n".join(out) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--patient", required=True)
    parser.add_argument("--vcf", type=Path,
                        help="VEP + readcount annotated somatic VCF. Omit when the run did not include "
                             "somatic calling: the variant, TMB, UV and clonality sections then say so.")
    parser.add_argument("--purity", type=float, required=True)
    parser.add_argument("--ploidy", type=float)
    parser.add_argument("--sex", type=str.lower, choices=["male", "female"])
    parser.add_argument("--segments", type=Path)
    parser.add_argument("--hla", default="", help="comma-separated alleles as given to pVACseq")
    parser.add_argument("--fasta", type=Path, help="indexed reference FASTA, for the UV dipyrimidine context")
    parser.add_argument("--fai", type=Path)
    parser.add_argument("--target-bed", type=Path, help="capture BED / interval_list, for the TMB denominator")
    parser.add_argument("--callable-mb", type=float, help="TMB denominator in Mb; overrides --target-bed")
    parser.add_argument("--contigs", help="comma-separated called contigs, for a WGS denominator from the .fai")
    parser.add_argument("--gene-panel", type=Path)
    parser.add_argument("--class-i", type=Path, help="MHC I aggregated report (expression- or CCF-annotated)")
    parser.add_argument("--class-ii", type=Path, help="MHC II aggregated report (expression- or CCF-annotated)")
    parser.add_argument("--pvacfuse-dir", type=Path,
                        help="pVACfuse output directory; absent means the fusion branch did not run")
    parser.add_argument("--fusion-predictions", type=Path, help="STAR-Fusion predictions table")
    parser.add_argument("--fusion-min-reads", type=float, default=5.0,
                        help="pVACfuse flags read support below this (default: 5)")
    parser.add_argument("--fusion-min-expr", type=float, default=0.1,
                        help="pVACfuse flags FFPM below this (default: 0.1)")
    parser.add_argument("--rna-support-summary", type=Path)
    parser.add_argument("--rna-support-figure", type=Path)
    parser.add_argument("--timepoint-label", action="append", help="RNA timepoint labels, in order")
    parser.add_argument("--project", default="")
    parser.add_argument("--diagnosis", default="Melanoma")
    parser.add_argument("--tcga-cancer-type", default="")
    parser.add_argument("--gtex-tissue", default="")
    parser.add_argument("--pipeline-version", default="dev")
    parser.add_argument("--run-date", default=datetime.date.today().isoformat())
    parser.add_argument("--clonal-ccf", type=float, default=0.8)
    parser.add_argument("--min-alt-reads", type=int, default=2)
    parser.add_argument("--class-i-ic50", type=float, default=500.0)
    parser.add_argument("--class-i-percentile", type=float, default=2.0)
    parser.add_argument("--class-ii-ic50", type=float, default=1000.0)
    parser.add_argument("--class-ii-percentile", type=float, default=10.0)
    parser.add_argument("--top-n", type=int, default=15)
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if not 0 < args.purity <= 1:
        parser.error("--purity must be in (0, 1]")
    args.output.write_text(build_report(args))
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
