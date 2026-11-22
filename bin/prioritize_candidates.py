#!/usr/bin/env python3
"""
Turn an annotated pVACseq aggregated report into one ranked, auditable
candidate list: a flag per failed criterion, a reason per rejected
candidate, and a rank.

Everything this reads is already in the table by the time it runs - the
binding calls from pVACseq, CCF from annotate_ccf.py, expression from
annotate_expression.py, HLA-LOH from parse_hla_loh.py. What is missing is
the integration: deciding which candidates are actually worth making
peptides for still means reading eight columns per row and holding the
thresholds in your head. This applies them once, writes down what failed,
and orders what is left.

Two derived columns are added that nothing upstream computes:

    Agretopicity   IC50 MT / IC50 WT. Below 1 the mutant binds better
                   than the wild type, so the mutation is what creates
                   the epitope and the T-cell repertoire against it is
                   less likely to have been tolerised. Near or above 1
                   the wild-type peptide is presented by normal cells too
                   (SELF_LIKE).
    DAI            IC50 WT - IC50 MT, the same comparison as a difference
                   rather than a ratio, since the ratio is unstable when
                   both are weak binders.

Flags are either rejections or cautions. A rejection means the candidate
should not go into a vaccine: its allele is gone, it does not bind, or
there is no evidence the gene is expressed. A caution means it is worse
than its neighbours but still a candidate. Nothing is deleted - every row
comes out, with "Flags", "Reject Reason", "Selection" and "Priority Rank",
so a rejection can always be read back and argued with.

The ranking is a sort, not a score: rejected rows last, then fewest
cautions, then clonal before subclonal, then most foreign (agretopicity),
then strongest binder, then highest expressed. No weights to tune, and two
candidates ranked next to each other differ in a way that can be named.

Usage:
    prioritize_candidates.py P1.Combined.aggregated.ccf.tsv \\
        [--hla-loh P1.hla_loh.tsv] -o P1.Combined.prioritized.tsv
"""
import argparse
import csv
import math
import re
import sys
from pathlib import Path

NEW_COLUMNS = ["Agretopicity", "DAI", "HLA LOH", "Flags", "Reject Reason", "Selection", "Priority Rank"]

# Flags that disqualify a candidate, in the order they are reported. The
# first one a row trips becomes its "Reject Reason".
REJECT_FLAGS = ["HLA_LOST", "LOW_BINDING", "NO_EXPRESSION", "NO_RNA_SUPPORT"]
# Flags that rank a candidate down but leave it in.
CAUTION_FLAGS = ["SELF_LIKE", "LOW_EXPRESSION", "NORMAL_TISSUE", "SUBCLONAL",
                 "LOW_DNA_VAF", "LOW_DEPTH", "NO_EXPRESSION_DATA", "NO_RNA_DATA", "HLA_NOT_EVALUABLE"]

EXPLANATIONS = {
    "HLA_LOST": "restricting HLA allele lost in the tumour",
    "LOW_BINDING": "predicted binding too weak",
    "NO_EXPRESSION": "gene not expressed in this tumour",
    "NO_RNA_SUPPORT": "covered by RNA reads but none carry the variant",
    "SELF_LIKE": "mutant binds no better than the wild type",
    "LOW_EXPRESSION": "gene weakly expressed",
    "NORMAL_TISSUE": "gene broadly expressed in normal tissue",
    "SUBCLONAL": "present in a minority of tumour cells",
    "LOW_DNA_VAF": "low variant allele fraction in tumour DNA",
    "LOW_DEPTH": "low DNA coverage at the variant",
    "NO_EXPRESSION_DATA": "no expression evidence either way",
    "NO_RNA_DATA": "no RNA coverage at the variant",
    "HLA_NOT_EVALUABLE": "HLA-LOH could not be called for this allele",
}


def number(raw):
    """A table cell -> float, or None for NA/./empty/inf/nan."""
    if raw is None:
        return None
    text = str(raw).strip().strip('"')
    if text in ("", "NA", ".", "None", "-"):
        return None
    try:
        value = float(text)
    except ValueError:
        return None
    return None if math.isnan(value) or math.isinf(value) else value


def read_loh(path):
    """The tidy LOH table -> {pVACseq allele: call}, class II pairs split.

    Only loci that were actually evaluable contribute a call, so a locus
    SpecHLA could not assess leaves its alleles unknown rather than
    looking retained.
    """
    calls = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            call = (row.get("call") or "").strip()
            lost = (row.get("lost_allele_pvacseq") or "").strip()
            for key in ("allele1_pvacseq", "allele2_pvacseq"):
                allele = (row.get(key) or "").strip()
                if not allele:
                    continue
                if call == "loss":
                    calls[allele] = "lost" if allele == lost else "retained"
                elif call in ("retained", "imbalance"):
                    calls[allele] = "retained"
                else:
                    calls.setdefault(allele, "not_evaluable")
    return calls


def loh_status(allele, calls):
    """The LOH call for a candidate's restricting allele.

    pVACseq names a class II DQ/DP heterodimer as 'DQA1*01:02-DQB1*06:02'.
    Losing either chain loses the molecule, so the pair is lost if either
    half is, and only retained if both halves are.
    """
    if not calls:
        return "NA"
    allele = (allele or "").strip()
    if not allele:
        return "NA"
    parts = [p for p in re.split(r"(?<=\d)-(?=[A-Z])", allele) if p] or [allele]
    statuses = [calls.get(part, "NA") for part in parts]
    if "lost" in statuses:
        return "lost"
    if all(status == "retained" for status in statuses):
        return "retained"
    if any(status == "not_evaluable" for status in statuses):
        return "not_evaluable"
    return "NA"


def agretopicity(row):
    """(ratio, DAI) from the MT/WT IC50 pair, or (None, None)."""
    mt, wt = number(row.get("IC50 MT")), number(row.get("IC50 WT"))
    if mt is None or wt is None or wt <= 0:
        return None, None
    return mt / wt, wt - mt


def flag_row(row, loh_calls, args):
    """-> (flags, hla status, agretopicity, DAI)."""
    flags = []
    ratio, dai = agretopicity(row)

    status = loh_status(row.get("Allele"), loh_calls)
    if status == "lost":
        flags.append("HLA_LOST")
    elif status == "not_evaluable":
        flags.append("HLA_NOT_EVALUABLE")

    # Binding: percentile and IC50 have to both be bad to reject, because
    # either alone is noisy and the two predictors disagree often enough.
    percentile, ic50 = number(row.get("%ile MT")), number(row.get("IC50 MT"))
    bad_percentile = percentile is not None and percentile > args.max_percentile
    bad_ic50 = ic50 is not None and ic50 > args.max_ic50
    if percentile is None and ic50 is None:
        pass
    elif (bad_percentile or percentile is None) and (bad_ic50 or ic50 is None):
        flags.append("LOW_BINDING")

    if ratio is not None and ratio > args.max_agretopicity:
        flags.append("SELF_LIKE")

    expression_flag = (row.get("Expression Flag") or "").strip()
    if expression_flag == "NO_EXPRESSION":
        flags.append("NO_EXPRESSION")
    elif expression_flag == "LOW_EXPRESSION":
        flags.append("LOW_EXPRESSION")
    elif expression_flag in ("NO_DATA", ""):
        flags.append("NO_EXPRESSION_DATA")

    if (row.get("GTEx Specificity") or "").strip() == "broadly_expressed":
        flags.append("NORMAL_TISSUE")

    # RNA: no coverage is a gap in the evidence; coverage with no variant
    # reads is evidence against the variant being transcribed.
    rna_depth, rna_vaf = number(row.get("RNA Depth")), number(row.get("RNA VAF"))
    if rna_depth is None or rna_depth < args.min_rna_depth:
        flags.append("NO_RNA_DATA")
    elif rna_vaf is not None and rna_vaf <= 0:
        flags.append("NO_RNA_SUPPORT")

    ccf = number(row.get("CCF"))
    if ccf is not None and ccf < args.min_ccf:
        flags.append("SUBCLONAL")

    dna_vaf = number(row.get("DNA VAF"))
    if dna_vaf is not None and dna_vaf < args.min_dna_vaf:
        flags.append("LOW_DNA_VAF")

    depth = number(row.get("DNA Depth")) or number(row.get("Tumor DNA Depth"))
    if depth is not None and depth < args.min_dna_depth:
        flags.append("LOW_DEPTH")

    return flags, status, ratio, dai


def sort_key(entry):
    """Rejected last, then fewest cautions, then the best biology first."""
    row, flags = entry["row"], entry["flags"]
    rejected = 1 if any(f in REJECT_FLAGS for f in flags) else 0
    cautions = sum(1 for f in flags if f in CAUTION_FLAGS)
    ccf = number(row.get("CCF"))
    ratio = entry["agretopicity"]
    percentile = number(row.get("%ile MT"))
    tpm = number(row.get("Sample TPM"))
    expr = number(row.get("RNA Expr"))
    return (
        rejected,
        cautions,
        -(ccf if ccf is not None else 0.0),          # clonal first
        ratio if ratio is not None else 1.0,          # most foreign first
        percentile if percentile is not None else 1e9,  # strongest binder first
        -(tpm if tpm is not None else (expr or 0.0)),   # highest expressed first
        str(row.get("ID", "")),                        # stable for ties
    )


def selection(flags):
    """A one-word group, so a wet-lab sheet can be filtered on it."""
    rejects = [f for f in flags if f in REJECT_FLAGS]
    if rejects:
        return "Rejected"
    cautions = [f for f in flags if f in CAUTION_FLAGS]
    if not cautions:
        return "Top ranking"
    return "Caution"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("report", type=Path, help="a CCF/expression-annotated pVACseq aggregated report")
    parser.add_argument("--hla-loh", type=Path, help="tidy LOH table from parse_hla_loh.py")
    parser.add_argument("--max-percentile", type=float, default=2.0, help="above this %%ile MT the binding is weak (default 2.0)")
    parser.add_argument("--max-ic50", type=float, default=500.0, help="above this IC50 MT (nM) the binding is weak (default 500)")
    parser.add_argument("--max-agretopicity", type=float, default=0.5,
                        help="above this IC50 MT/WT ratio the mutant is no better than self (default 0.5)")
    parser.add_argument("--min-ccf", type=float, default=0.5, help="below this CCF a candidate is subclonal (default 0.5)")
    parser.add_argument("--min-dna-vaf", type=float, default=0.05, help="below this tumour DNA VAF (default 0.05)")
    parser.add_argument("--min-dna-depth", type=float, default=20, help="below this tumour DNA depth (default 20)")
    parser.add_argument("--min-rna-depth", type=float, default=5, help="RNA depth needed before absent variant reads mean anything (default 5)")
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args(argv)

    loh_calls = read_loh(args.hla_loh) if args.hla_loh else {}

    with args.report.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        entries = []
        for row in reader:
            flags, status, ratio, dai = flag_row(row, loh_calls, args)
            entries.append({"row": row, "flags": flags, "hla": status, "agretopicity": ratio, "dai": dai})

    entries.sort(key=sort_key)

    fmt = lambda x: "NA" if x is None else f"{x:.3g}"
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fieldnames + NEW_COLUMNS)
        for rank, entry in enumerate(entries, start=1):
            flags = entry["flags"]
            reason = next((EXPLANATIONS[f] for f in REJECT_FLAGS if f in flags), "")
            writer.writerow(
                [entry["row"].get(name, "") for name in fieldnames]
                + [fmt(entry["agretopicity"]), fmt(entry["dai"]), entry["hla"],
                   ";".join(flags), reason, selection(flags), rank]
            )

    groups = {}
    for entry in entries:
        groups[selection(entry["flags"])] = groups.get(selection(entry["flags"]), 0) + 1
    summary = ", ".join(f"{count} {name.lower()}" for name, count in sorted(groups.items()))
    print(f"{len(entries)} candidates: {summary or 'none'}", file=sys.stderr)


if __name__ == "__main__":
    main()
