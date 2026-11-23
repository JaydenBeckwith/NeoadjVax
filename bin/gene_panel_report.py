#!/usr/bin/env python3
"""
One-page status table for genes that decide whether a tumour can present
antigen at all, and for the melanoma drivers.

A ranked neoantigen list is worth little if the tumour has lost the
machinery to present it, or if the interferon signalling that upregulates
that machinery is broken: B2M loss, HLA loss, JAK1/JAK2 truncation and
CIITA/NLRC5 silencing are recurrent routes to checkpoint-inhibitor
resistance, and they are invisible in a per-candidate report because the
defect is in the presentation pathway, not in the candidate. The driver
half is context of a different kind: which of the canonical melanoma
drivers this tumour carries, and whether the gene is expressed.

For every panel gene this combines three independent read-outs:

  * expression - this sample's TPM against the matched GTEx normal tissue
    and the matched TCGA cohort, so "low" means low for that tissue rather
    than low in the abstract;
  * copy number - the Sequenza segment covering the gene, with minor copy
    number 0 reported as LOH;
  * mutations - the somatic consequences VEP called in that gene, with
    truncating ones (stop gained, frameshift, splice acceptor/donor, start
    lost) marked, since those are the ones that take a presentation gene
    out.

Each is optional: give what the run produced and the rest reads NA.

Panel Flag is deliberately factual, not a verdict, and which direction
counts as adverse depends on the gene:

    LOSS_OF_FUNCTION  loss-polarity gene: expression DOWN, CN loss/LOH, or
                      a truncating mutation
    INHIBITORY_GAIN   gain-polarity gene: expression UP or amplified
    ALTERED           a non-truncating mutation, or a change in the
                      direction that is not the adverse one for this gene
    OK                nothing of the above
    NA                no data for this gene

Most of the panel is loss polarity: for HLA-A/B/C, B2M, TAP, the
immunoproteasome, the interferon genes, MR1, CD1 and the NKG2D ligands,
losing the gene is what lets a tumour escape. HLA-E and HLA-G are the
other way round - they are inhibitory ligands, so a tumour that
upregulates them is suppressing the response, and gain is the finding.
For the drivers neither direction is "adverse": a BRAF or NRAS hotspot is
activating, not broken, so read the mutation column rather than the flag.

Caveats that matter when reading the table:

  * HLA genes are the most polymorphic in the genome, so short reads
    mismap there: both their expression and their copy number are less
    trustworthy than the rest of the panel, and real HLA allelic loss
    needs a dedicated caller (LOHHLA, SpecHLA), not this. The same applies
    to HLA-E/F/G and to MICA/MICB, which are polymorphic and sit in the
    same region.
  * MICA/MICB surface loss is mostly proteolytic shedding, which is
    post-translational and so invisible to both RNA and DNA: a normal row
    here does not mean the ligand is on the cell surface.
  * HLA-F's immunological direction is unsettled (open conformers,
    KIR3DS1), so it is left at loss polarity by default; do not read its
    flag as strongly as the rest.
  * HLA-DRB3/4/5 are absent from GENCODE and GTEx, so they cannot appear.
  * TERT's recurrent melanoma mutations are promoter, not coding, and are
    usually not in a somatic exome call set at all; an empty TERT mutation
    cell is not evidence of a wild-type promoter.
  * Expression is bulk and gene level: stromal and immune content dilutes
    tumour expression, which pushes genes towards DOWN at low purity.

Usage:
    gene_panel_report.py -o P1.panel.tsv \\
        [--sample-tpm P1.rna_baseline.tpm.tsv] \\
        [--gtex gtex.tsv --gtex-tissue "Skin - Sun Exposed (Lower leg)"] \\
        [--tcga tcga.tsv --tcga-cancer-type TCGA-SKCM] \\
        [--vep-vcf P1.somatic.vep.vcf] [--segments P1_segments.txt --gtf genes.gtf] \\
        [--ploidy 3.1] [--panel my_panel.tsv]
"""
import argparse
import bisect
import csv
import gzip
import math
import re
from pathlib import Path

# Pathway -> genes. Grouped so a defect in one route (say class II
# regulation) is visible as a block rather than as scattered rows.
# Override wholesale with --panel.
DEFAULT_PANEL = {
    "apm_class_i": ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2",
                    "TAPBP", "CALR", "PDIA3"],
    # Class Ib and the lipid/metabolite presenters. HLA-E and HLA-G are
    # inhibitory ligands (NKG2A, ILT2/ILT4), so for those two it is gain,
    # not loss, that signals escape - see GAIN_POLARITY.
    "apm_non_classical": ["HLA-E", "HLA-F", "HLA-G", "MR1",
                          "CD1A", "CD1B", "CD1C", "CD1D", "CD1E"],
    # MHC class I-related, not HLA: NKG2D ligands read by NK and gamma-delta
    # T cells. Losing them is an escape route from NK killing.
    "nkg2d_ligands": ["MICA", "MICB", "ULBP1", "ULBP2", "ULBP3",
                      "RAET1E", "RAET1G", "RAET1L"],
    "apm_class_ii": ["CD74", "HLA-DMA", "HLA-DMB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1",
                     "HLA-DQB1", "HLA-DRA", "HLA-DRB1"],
    "apm_regulation": ["CIITA", "NLRC5", "RFX5", "RFXANK", "RFXAP", "NFYA", "NFYB", "NFYC", "CREB1"],
    "immunoproteasome": ["PSMB5", "PSMB6", "PSMB7", "PSMB8", "PSMB9", "PSMB10",
                         "PSME1", "PSME2", "PSME3", "POMP", "PSMF1"],
    # JAK1/JAK2 and the rest of this route are recurrent checkpoint-inhibitor
    # resistance mechanisms, and they switch the genes above on
    "interferon_signalling": ["IFNG", "IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1", "IRF1", "IRF2"],
    "driver_mapk": ["BRAF", "NRAS", "KRAS", "HRAS", "NF1", "MAP2K1", "MAP2K2", "RAC1",
                    "PPP6C", "KIT", "GNAQ", "GNA11"],
    "driver_cell_cycle_p53": ["TP53", "CDKN2A", "CDK4", "CCND1", "RB1", "MDM2"],
    "driver_pi3k": ["PTEN", "PIK3CA", "AKT1"],
    "driver_other": ["ARID2", "IDH1", "STK19", "DDX3X", "SF3B1", "BAP1", "EIF1AX",
                     "MITF", "TERT", "CTNNB1"],
}

# Genes where gain, not loss, is the adverse finding. HLA-E engages the
# inhibitory receptor NKG2A on NK and CD8 T cells and HLA-G engages
# ILT2/ILT4 and KIR2DL4, so a tumour that upregulates them is actively
# suppressing the response rather than failing to present. Everything else
# in the panel defaults to loss polarity. Override per gene with a third
# column in --panel.
GAIN_POLARITY = {"HLA-E", "HLA-G"}

TRUNCATING = ("stop_gained", "frameshift_variant", "splice_acceptor_variant",
              "splice_donor_variant", "start_lost", "stop_lost", "transcript_ablation")

ATTRIBUTE_RE = re.compile(r'(\S+) "([^"]*)"')


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def parse_float(raw):
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(value) else value


def fmt(value, digits=2):
    return "NA" if value is None else f"{value:.{digits}f}"


def norm_chrom(chrom):
    chrom = chrom.strip().strip('"')
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def read_panel(path):
    """gene<TAB>pathway[<TAB>polarity] -> ({pathway: [genes]}, {gain genes}).

    The third column is optional and is "loss" (the default) or "gain";
    "gain" marks a gene where upregulation rather than loss is the adverse
    finding, as for the inhibitory class Ib ligands.
    """
    panel, gain = {}, set()
    with open_text(path) as handle:
        for line in handle:
            fields = [f.strip() for f in line.rstrip("\n").split("\t")]
            if len(fields) < 2 or not fields[0] or fields[0].lower() in ("gene", "#gene"):
                continue
            panel.setdefault(fields[1], []).append(fields[0])
            if len(fields) > 2 and fields[2].lower() == "gain":
                gain.add(fields[0])
            elif len(fields) > 2 and fields[2] and fields[2].lower() != "loss":
                raise SystemExit(f"{path}: polarity for {fields[0]} must be 'loss' or 'gain', got '{fields[2]}'")
    if not panel:
        raise SystemExit(f"{path}: no gene<TAB>pathway rows found")
    return panel, gain


def read_expression_db(path):
    """build_expression_db.py table -> {gene_key: {group: {quantile: tpm}}}."""
    by_key = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = {"ENSEMBL_ID", "GENE_NAME", "GROUP", "QUANTILE", "TPM"} - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"{path}: missing column(s) {sorted(missing)}: rebuild with bin/build_expression_db.py")
        for row in reader:
            tpm, quantile = parse_float(row["TPM"]), parse_float(row["QUANTILE"])
            if tpm is None or quantile is None:
                continue
            for key in (row["ENSEMBL_ID"], row["GENE_NAME"]):
                if not key or key == "NA":
                    continue
                groups = by_key.setdefault(key, {}).setdefault(row["GROUP"], {})
                if quantile not in groups or tpm > groups[quantile]:
                    groups[quantile] = tpm
    return by_key


def read_sample_tpm(path):
    by_key = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            tpm = parse_float(row.get("TPM"))
            if tpm is None:
                continue
            for key in (row.get("ENSEMBL_ID"), row.get("GENE_NAME")):
                if key and key != "NA" and (key not in by_key or tpm > by_key[key]):
                    by_key[key] = tpm
    return by_key


def read_gene_spans(path):
    """GTF -> {gene_name: (chrom, start, end)} from gene rows."""
    spans = {}
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attributes = dict(ATTRIBUTE_RE.findall(fields[8]))
            name = attributes.get("gene_name")
            if name:
                spans[name] = (norm_chrom(fields[0]), int(fields[3]), int(fields[4]))
    return spans


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


def gene_copy_number(segments, span):
    """Mean total/minor CN over the segments the gene overlaps.

    A gene can straddle a breakpoint, which is itself worth seeing, so the
    overlapping segments are averaged (weighted by overlap) rather than one
    being picked; the minor copy number is the minimum across them, so a
    gene half of which has lost an allele still reports LOH.
    """
    rows = segments.get(span[0], [])
    if not rows:
        return None, None
    start, end = span[1], span[2]
    index = max(0, bisect.bisect_right([r[0] for r in rows], start) - 1)
    total_weight, weighted_cn, minor = 0.0, 0.0, None
    for row in rows[index:]:
        if row[0] > end:
            break
        overlap = min(end, row[1]) - max(start, row[0]) + 1
        if overlap <= 0 or row[2] is None:
            continue
        total_weight += overlap
        weighted_cn += overlap * row[2]
        if row[4] is not None:
            minor = row[4] if minor is None else min(minor, row[4])
    if not total_weight:
        return None, None
    return weighted_cn / total_weight, minor


def read_vep_mutations(path, genes):
    """VEP-annotated VCF -> {gene: [(consequence, protein change)]}.

    Reads SYMBOL, Consequence and HGVSp (falling back to the amino acids
    and protein position) out of the CSQ field, keeping one entry per
    distinct consequence/change so the same variant on many transcripts is
    not listed many times. An empty file (the stub) yields nothing.
    """
    mutations = {}
    csq_format = None
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##"):
                if line.startswith("##INFO=<ID=CSQ") and "Format:" in line:
                    csq_format = line.split("Format:")[1].strip().rstrip('">').split("|")
                continue
            if line.startswith("#") or not line.strip():
                continue
            if not csq_format:
                break
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            for info in fields[7].split(";"):
                if not info.startswith("CSQ="):
                    continue
                for entry in info[4:].split(","):
                    values = entry.split("|")
                    data = dict(zip(csq_format, values))
                    symbol = data.get("SYMBOL")
                    if not symbol or symbol not in genes:
                        continue
                    consequence = (data.get("Consequence") or "").split("&")[0]
                    change = data.get("HGVSp", "")
                    if change:
                        change = change.split(":")[-1]
                    elif data.get("Amino_acids") and data.get("Protein_position"):
                        aas = data["Amino_acids"].split("/")
                        change = f"p.{aas[0]}{data['Protein_position']}" + (aas[1] if len(aas) > 1 else "")
                    pair = (consequence, change or "NA")
                    if pair not in mutations.setdefault(symbol, []):
                        mutations[symbol].append(pair)
    return mutations


def expression_call(sample_tpm, quantiles, fold):
    """(call, basis, low, high) for one gene in one group.

    With two or more quantiles in the table the sample is called against
    that range. With only a median (the GTEx summary file), a fold change
    either side of it stands in - cruder, so the basis is reported rather
    than left implicit.
    """
    if not quantiles:
        return "NA", "no_reference", None, None
    values = sorted(quantiles)
    if len(values) >= 2:
        low, high = quantiles[values[0]], quantiles[values[-1]]
        basis = f"q{values[0]:g}-q{values[-1]:g}"
    else:
        median = quantiles[values[0]]
        low, high, basis = median / fold, median * fold, f"median_x{fold:g}"
    if sample_tpm is None:
        return "NA", basis, low, high
    if sample_tpm < low:
        return "DOWN", basis, low, high
    if sample_tpm > high:
        return "UP", basis, low, high
    return "NORMAL", basis, low, high


def copy_number_call(total, minor, baseline, amp_gain):
    if total is None:
        return "NA"
    calls = []
    if total <= baseline - 1:
        calls.append("LOSS")
    elif total >= baseline + amp_gain:
        calls.append("AMP")
    if minor is not None and minor == 0:
        calls.append("LOH")
    return ",".join(calls) if calls else "NEUTRAL"


def panel_flag(expression, cn_call, mutations, polarity="loss"):
    """The adverse finding for this gene, or ALTERED/OK/NA.

    A gain-polarity gene that is lost or truncated is not flagged
    LOSS_OF_FUNCTION: losing an inhibitory ligand does not help the tumour,
    so it is recorded as ALTERED rather than dressed up as a defect.
    """
    truncating = [m for m in mutations if m[0] in TRUNCATING]
    lost = expression == "DOWN" or "LOSS" in cn_call or "LOH" in cn_call or truncating
    gained = expression == "UP" or "AMP" in cn_call
    if polarity == "gain":
        if gained:
            return "INHIBITORY_GAIN"
    elif lost:
        return "LOSS_OF_FUNCTION"
    if mutations or lost or gained:
        return "ALTERED"
    if expression == "NA" and cn_call == "NA" and not mutations:
        return "NA"
    return "OK"


def parse_timepoint(spec):
    """'Week 6=P1.week6.tpm.tsv' -> ('Week 6', Path). Mirrors agretopicity.py."""
    label, separator, path = str(spec).partition("=")
    if not separator or not label.strip() or not path.strip():
        raise SystemExit(f"--timepoint needs LABEL=TPM_TABLE, got: {spec}")
    return label.strip(), Path(path.strip())


def log2_fold_change(later, baseline):
    """log2((later + 1) / (baseline + 1)). The pseudocount keeps a gene that
    switches on from zero finite, which is exactly the case worth seeing."""
    if later is None or baseline is None:
        return None
    return math.log2((later + 1.0) / (baseline + 1.0))


def columns(timepoint_labels):
    """COLUMNS with a TPM and a log2FC column per extra timepoint.

    The extra columns sit next to Sample TPM rather than at the end, so the
    expression story reads left to right in the table.
    """
    extra = []
    for label in timepoint_labels:
        extra += [f"TPM ({label})", f"Log2FC ({label})"]
    head = ["Pathway", "Gene", "Sample TPM"] + extra
    return head + ["GTEx Range", "GTEx Basis", "GTEx Call", "TCGA Range", "TCGA Call",
                   "Tumor CN", "Minor CN", "CN Call", "Mutations", "Truncating",
                   "Adverse Direction", "Panel Flag"]


def build_rows(panel, gain_polarity, args, sample, gtex, tcga, spans, segments, mutations,
               timepoints=()):
    baseline = args.ploidy if args.ploidy else 2.0
    rows = []
    for pathway in panel:
        for gene in panel[pathway]:
            sample_tpm = sample.get(gene)
            gtex_q = gtex.get(gene, {}).get(args.gtex_tissue, {}) if gtex else {}
            tcga_q = tcga.get(gene, {}).get(args.tcga_cancer_type, {}) if tcga else {}
            call, basis, low, high = expression_call(sample_tpm, gtex_q, args.fold)
            tcga_call, _, tcga_low, tcga_high = expression_call(sample_tpm, tcga_q, args.fold)

            total_cn = minor_cn = None
            if segments and gene in spans:
                total_cn, minor_cn = gene_copy_number(segments, spans[gene])
            cn_call = copy_number_call(total_cn, minor_cn, baseline, args.amp_gain)

            gene_mutations = mutations.get(gene, [])
            formatted = ";".join(f"{consequence}:{change}" for consequence, change in gene_mutations) or "-"
            truncating = "yes" if any(m[0] in TRUNCATING for m in gene_mutations) else "no"

            polarity = "gain" if gene in gain_polarity else "loss"
            timepoint_cells = []
            for _label, table in timepoints:
                later = table.get(gene)
                timepoint_cells += [fmt(later), fmt(log2_fold_change(later, sample_tpm), 2)]
            rows.append([pathway, gene, fmt(sample_tpm)] + timepoint_cells +
                        [f"{fmt(low)}-{fmt(high)}" if low is not None else "NA",
                         basis, call, f"{fmt(tcga_low)}-{fmt(tcga_high)}" if tcga_low is not None else "NA",
                         tcga_call, fmt(total_cn, 2), fmt(minor_cn, 2), cn_call, formatted, truncating,
                         polarity, panel_flag(call, cn_call, gene_mutations, polarity)])
    return rows


COLUMNS = ["Pathway", "Gene", "Sample TPM", "GTEx Range", "GTEx Basis", "GTEx Call",
           "TCGA Range", "TCGA Call", "Tumor CN", "Minor CN", "CN Call", "Mutations",
           "Truncating", "Adverse Direction", "Panel Flag"]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--panel", type=Path,
                        help="gene<TAB>pathway[<TAB>loss|gain] table replacing the built-in panel")
    parser.add_argument("--sample-tpm", type=Path, help="this sample's TPM table from gene_tpm.py")
    parser.add_argument("--timepoint", action="append", default=[], metavar="LABEL=TPM_TABLE",
                        help="a further timepoint's TPM table, e.g. 'Week 6=P1.week6.tpm.tsv'; "
                             "adds a TPM and a log2FC column per timepoint")
    parser.add_argument("--gtex", type=Path)
    parser.add_argument("--gtex-tissue")
    parser.add_argument("--tcga", type=Path)
    parser.add_argument("--tcga-cancer-type")
    parser.add_argument("--vep-vcf", type=Path, help="VEP-annotated somatic VCF, for mutations")
    parser.add_argument("--segments", type=Path, help="Sequenza segments, for copy number")
    parser.add_argument("--gtf", type=Path, help="GTF for gene coordinates (needed with --segments)")
    parser.add_argument("--ploidy", type=float, help="sample ploidy; copy number is called against it (default 2)")
    parser.add_argument("--amp-gain", type=float, default=2.0,
                        help="copies above baseline before a gene is called AMP (default: 2)")
    parser.add_argument("--fold", type=float, default=4.0,
                        help="fold change either side of the median when the table has no range (default: 4)")
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args(argv)

    if args.gtex and not args.gtex_tissue:
        parser.error("--gtex needs --gtex-tissue")
    if args.tcga and not args.tcga_cancer_type:
        parser.error("--tcga needs --tcga-cancer-type")
    if args.segments and not args.gtf:
        parser.error("--segments needs --gtf to locate each gene")

    panel, gain_polarity = read_panel(args.panel) if args.panel else (DEFAULT_PANEL, GAIN_POLARITY)
    genes = {gene for pathway in panel for gene in panel[pathway]}

    gtex = read_expression_db(args.gtex) if args.gtex else {}
    tcga = read_expression_db(args.tcga) if args.tcga else {}
    sample = read_sample_tpm(args.sample_tpm) if args.sample_tpm else {}
    timepoints = []
    for spec in args.timepoint:
        label, path = parse_timepoint(spec)
        if not path.is_file():
            raise SystemExit(f"--timepoint {label}: no such file: {path}")
        timepoints.append((label, read_sample_tpm(path)))
    if timepoints and not args.sample_tpm:
        raise SystemExit("--timepoint needs --sample-tpm: the log2FC is measured against it")
    spans = read_gene_spans(args.gtf) if args.gtf else {}
    segments = read_segments(args.segments) if args.segments else None
    mutations = read_vep_mutations(args.vep_vcf, genes) if args.vep_vcf else {}

    for name, db, group in (("--gtex-tissue", gtex, args.gtex_tissue),
                            ("--tcga-cancer-type", tcga, args.tcga_cancer_type)):
        if db and group not in {g for groups in db.values() for g in groups}:
            available = sorted({g for groups in db.values() for g in groups})
            raise SystemExit(f"{name} '{group}' is not in the table. Available: " + ", ".join(available))

    rows = build_rows(panel, gain_polarity, args, sample, gtex, tcga, spans, segments, mutations,
                      timepoints)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns([label for label, _ in timepoints]))
        writer.writerows(rows)

    flagged = [row for row in rows if row[-1] in ("LOSS_OF_FUNCTION", "INHIBITORY_GAIN")]
    print(f"{len(rows)} panel genes; {len(flagged)} with an adverse finding")
    for row in flagged:
        print(f"  {row[0]}\t{row[1]}\t{row[-1]}\texpression={row[5]}\tCN={row[10]}\tmutations={row[11]}")


if __name__ == "__main__":
    main()
