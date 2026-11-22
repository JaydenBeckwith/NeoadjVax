#!/usr/bin/env python3
"""
Normalise a SpecHLA HLA-LOH result into a tidy per-locus table, with the
alleles renamed to the pVACseq spelling so candidates can be flagged.

Input is what `cal.hla.copy.pl` writes (`merge.hla.copy.txt`), one row per
locus:

    Sample  HLA  Allele1  Allele2  copyratio  KeptHLA  LossHLA  Freq1
    Freq2  Purity  Het_num  LOH  SampleDir  Homozygous  Class

A long collated table with lower-case names (`locus`, `allele1`, ...) is
read as well, so the per-cohort CSVs already produced by the NeoadjLOH
analysis can be fed in directly. Comma-, tab- and semicolon-separated
files are all accepted.

Why this is not just "take the LOH column": SpecHLA emits LOH = Y for any
allelic imbalance, including a 2:1 locus where both alleles are still
present and still presenting. A vaccine target is only invalidated when
its allele is actually gone, so the default call here needs a 0 on one
side of `copyratio` (`--definition copy_loss`). `--definition spechla`
keeps SpecHLA's own call, which is the looser sensitivity analysis.

A locus is only evaluable when the two alleles differ at two-field
resolution (nothing to lose if they are the same), SpecHLA did not mark it
homozygous, and it has at least `--min-het` heterozygous sites. Below that
the call is `not_evaluable` rather than `retained`: absence of evidence,
not evidence of retention. Allelic-imbalance calls are also unreliable at
low purity, so `--min-purity` (default 0.20, as in the cohort analysis)
marks the whole sample `low_purity` instead of trusting it.

Defaults match the NeoadjLOH cohort analysis (min-het 10, min-purity 0.20,
copy_loss), so the pipeline and that analysis call the same patients.

Usage:
    parse_hla_loh.py merge.hla.copy.txt --patient NeoTrio_003 \\
        -o NeoTrio_003.hla_loh.tsv
"""
import argparse
import csv
import re
import sys
from pathlib import Path

# SpecHLA's loci, and which are class I. Class II alleles are named without
# the "HLA-" prefix by pVACseq, class I alleles with it.
CLASS_I = ("A", "B", "C", "E", "F", "G")

OUT_COLUMNS = [
    "patient", "locus", "hla_class", "allele1", "allele2",
    "allele1_pvacseq", "allele2_pvacseq", "copyratio", "cn1", "cn2",
    "freq1", "freq2", "het_num", "purity", "spechla_loh",
    "heterozygous_2field", "evaluable", "call",
    "lost_allele", "lost_allele_pvacseq", "kept_allele_pvacseq",
]

# Header spellings seen in the wild -> our name. Matched case-insensitively.
ALIASES = {
    # A collated table carries both: sample_id is the sequencing id,
    # patient_id ours, and --patient is matched against either.
    "sample": "sample", "sample_id": "sample",
    "patient": "patient_id", "patient_id": "patient_id",
    "hla": "locus", "locus": "locus",
    "allele1": "allele1", "allele2": "allele2",
    "copyratio": "copyratio", "copy_ratio": "copyratio",
    "kepthla": "kept", "kept": "kept",
    "losshla": "lost", "loss": "lost", "lost": "lost",
    "freq1": "freq1", "freq2": "freq2",
    "purity": "purity",
    "het_num": "het_num", "hetnum": "het_num",
    # loh_raw is SpecHLA's own call in the collated analysis tables, where
    # plain `loh` is already a derived TRUE/FALSE. Both are read, and
    # whichever appears first in the header wins.
    "loh_raw": "loh", "loh": "loh",
    "homozygous": "homozygous",
    "class": "class",
}


def two_field(allele):
    """A*02:01:01:01 -> A*02:01; anything unparseable comes back unchanged."""
    if not allele:
        return ""
    match = re.match(r"^([A-Za-z0-9]+\*\d+:\d+)", allele.strip())
    return match.group(1) if match else allele.strip()


def pvacseq_name(allele):
    """Two-field allele in the spelling pVACseq reports, or "" if unparseable.

    Class I gains the HLA- prefix (HLA-A*02:01); class II does not
    (DRB1*15:01), which is how pVACtools names them.
    """
    short = two_field(allele)
    if "*" not in short:
        return ""
    locus = short.split("*", 1)[0].upper()
    if locus.startswith("HLA-"):
        locus, short = locus[4:], short[4:]
    return f"HLA-{short}" if locus in CLASS_I else short


def split_copyratio(raw):
    """'2:0' -> (2.0, 0.0). Anything else -> (None, None)."""
    parts = str(raw or "").strip().split(":")
    if len(parts) != 2:
        return None, None
    try:
        return float(parts[0]), float(parts[1])
    except ValueError:
        return None, None


def number(raw):
    try:
        return float(str(raw).strip())
    except (TypeError, ValueError):
        return None


def read_rows(path):
    """The LOH file -> [{our name: value}], whichever delimiter it uses."""
    text = path.read_text()
    if not text.strip():
        return []
    head = text.splitlines()[0]
    delimiter = "\t"
    for candidate in ("\t", ",", ";"):
        if candidate in head:
            delimiter = candidate
            break
    reader = csv.DictReader(text.splitlines(), delimiter=delimiter)
    if not reader.fieldnames:
        return []
    mapping = {}
    for field in reader.fieldnames:
        key = ALIASES.get((field or "").strip().strip('"').lower())
        # First spelling wins, so a file carrying both `sample_id` and
        # `patient_id` keeps the one named first rather than the last.
        if key and key not in mapping.values():
            mapping[field] = key
    if "locus" not in mapping.values():
        raise SystemExit(f"{path}: no locus column (expected 'HLA' or 'locus'); header was: {', '.join(reader.fieldnames)}")
    rows = []
    for raw in reader:
        row = {our: (raw.get(field) or "").strip() for field, our in mapping.items()}
        if row.get("locus"):
            rows.append(row)
    return rows


def is_yes(raw):
    """SpecHLA writes Y/N; the collated tables write TRUE/FALSE."""
    return str(raw or "").strip().upper() in ("Y", "YES", "TRUE", "T", "1")


def classify(row, min_het, definition, purity_ok):
    """-> (heterozygous_2field, evaluable, call, lost_allele)."""
    a1, a2 = two_field(row.get("allele1")), two_field(row.get("allele2"))
    het = bool(a1) and bool(a2) and a1 != a2
    homozygous_flag = is_yes(row.get("homozygous"))
    lost_raw = (row.get("lost") or "").strip()
    het_num = number(row.get("het_num"))
    evaluable = (
        het
        and not homozygous_flag
        and lost_raw.lower() != "homogeneous"
        and het_num is not None
        and het_num >= min_het
    )
    if not evaluable:
        return het, False, "not_evaluable", ""
    if not purity_ok:
        # Evaluable on its own terms, but the purity is too low for an
        # imbalance call to mean anything, so it is not used either.
        return het, False, "low_purity", ""

    spechla_loh = is_yes(row.get("loh"))
    cn1, cn2 = split_copyratio(row.get("copyratio"))
    zero_side = cn1 == 0 or cn2 == 0
    if definition == "copy_loss":
        is_loss = spechla_loh and zero_side
    else:
        is_loss = spechla_loh
    if not is_loss:
        return het, True, "imbalance" if spechla_loh else "retained", ""

    # Which allele went. SpecHLA's LossHLA is authoritative when it names
    # one; otherwise fall back to the side of copyratio that is at zero.
    lost = lost_raw
    if not lost or lost.lower() in ("na", "homogeneous", "none"):
        lost = row.get("allele1", "") if cn1 == 0 else row.get("allele2", "") if cn2 == 0 else ""
    return het, True, "loss", lost


def tidy(rows, patient, min_het, definition, min_purity):
    out = []
    for row in rows:
        purity = number(row.get("purity"))
        purity_ok = min_purity <= 0 or purity is None or purity >= min_purity
        het, evaluable, call, lost = classify(row, min_het, definition, purity_ok)
        cn1, cn2 = split_copyratio(row.get("copyratio"))
        locus = (row.get("locus") or "").strip().upper()
        if locus.startswith("HLA-"):
            locus = locus[4:]
        a1, a2 = row.get("allele1", ""), row.get("allele2", "")
        kept = row.get("kept", "") if call == "loss" else ""
        if call == "loss" and not kept:
            kept = a2 if two_field(lost) == two_field(a1) else a1
        out.append({
            "patient": patient or row.get("patient_id") or row.get("sample", ""),
            "locus": locus,
            "hla_class": "I" if locus in CLASS_I else "II",
            "allele1": a1,
            "allele2": a2,
            "allele1_pvacseq": pvacseq_name(a1),
            "allele2_pvacseq": pvacseq_name(a2),
            "copyratio": row.get("copyratio", ""),
            "cn1": "NA" if cn1 is None else f"{cn1:g}",
            "cn2": "NA" if cn2 is None else f"{cn2:g}",
            "freq1": row.get("freq1", ""),
            "freq2": row.get("freq2", ""),
            "het_num": row.get("het_num", ""),
            "purity": row.get("purity", ""),
            "spechla_loh": row.get("loh", ""),
            "heterozygous_2field": "TRUE" if het else "FALSE",
            "evaluable": "TRUE" if evaluable else "FALSE",
            "call": call,
            "lost_allele": lost,
            "lost_allele_pvacseq": pvacseq_name(lost),
            "kept_allele_pvacseq": pvacseq_name(kept),
        })
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("loh", type=Path, help="SpecHLA merge.hla.copy.txt, or a collated long table")
    parser.add_argument("--patient", help="patient id to write in the output, and to select on when the file holds several samples")
    parser.add_argument("--min-het", type=float, default=10, help="minimum heterozygous sites for a locus to be evaluable (default 10)")
    parser.add_argument("--min-purity", type=float, default=0.20, help="below this purity no LOH call is trusted; 0 disables (default 0.20)")
    parser.add_argument("--definition", choices=["copy_loss", "spechla"], default="copy_loss",
                        help="copy_loss: one allele at 0 copies (default). spechla: any SpecHLA LOH=Y, including imbalance")
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args(argv)

    rows = read_rows(args.loh)
    if args.patient:
        # Select on the sample column only when it would leave something:
        # a single-sample file often carries the SpecHLA run name rather
        # than our patient_id, and dropping every row there would silently
        # turn into "no LOH anywhere".
        named = [r for r in rows if args.patient in (r.get("sample", ""), r.get("patient_id", ""))]
        if named:
            rows = named
        elif any(r.get("sample") or r.get("patient_id") for r in rows):
            samples = sorted({r.get("sample") or r.get("patient_id") for r in rows if r.get("sample") or r.get("patient_id")})
            if len(samples) > 1:
                raise SystemExit(
                    f"{args.loh}: holds several samples ({', '.join(samples[:5])}"
                    f"{', ...' if len(samples) > 5 else ''}) and none is '{args.patient}'"
                )
            print(f"note: {args.loh} is for sample '{samples[0]}', reporting it as {args.patient}", file=sys.stderr)

    table = tidy(rows, args.patient, args.min_het, args.definition, args.min_purity)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(table)

    lost = [r for r in table if r["call"] == "loss"]
    print(f"{len(table)} loci: {sum(1 for r in table if r['evaluable'] == 'TRUE')} evaluable, {len(lost)} lost", file=sys.stderr)
    for row in lost:
        print(f"  {row['locus']}: lost {row['lost_allele_pvacseq'] or row['lost_allele']} (copyratio {row['copyratio']})", file=sys.stderr)


if __name__ == "__main__":
    main()
