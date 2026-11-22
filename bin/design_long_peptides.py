#!/usr/bin/env python3
"""
Turn ranked candidates into the long peptides you would actually order.

A pVACseq candidate is a 9- or 15-mer epitope. A synthetic long peptide
vaccine is made of 20-30mers that carry the mutation with natural flanking
sequence, so the peptide is processed and presented by the patient's own
antigen-presenting cells on whichever allele fits, rather than betting the
construct on one predicted epitope. Going from one to the other is the
step that was being done by hand.

Input is the ranked table from prioritize_candidates.py plus the mutant
protein sequences from `pvacseq generate_protein_fasta`, which supplies
the flanking context the candidate tables do not carry.

For each candidate:

  * the mutated residues are located by comparing the MT and WT sequences
    of the same fasta entry, so this does not depend on how a particular
    pVACtools version spells its fasta headers. A substitution gives one
    position; a frameshift or insertion gives a novel stretch running to
    the end of the sequence, and the window is anchored at its start,
    because the junction is the part that is new to the immune system.
  * a window of --max-length is centred on the mutation, pulled inside the
    protein at the termini, and shrunk when the protein is shorter.
  * the ends are trimmed off residues that make synthesis or storage
    unreliable (--bad-nterm, --bad-cterm), as long as the mutation keeps
    --min-flank residues of context on both sides.
  * candidates whose mutations fall inside one window - neighbouring
    mutations in the same protein - are merged into a single peptide,
    listed in "Covers".

Synthesis notes are reported, never acted on: cysteine and methionine
counts (oxidation and disulfide trouble), GRAVY hydrophobicity (above 1.0
tends to be insoluble), net charge, and the longest single-residue run.
They go in a "Warnings" column for whoever places the order.

Rejected candidates are skipped by default (--include-rejected keeps
them), and --top limits how many peptides come out per class.

Usage:
    design_long_peptides.py P1.Combined.prioritized.tsv \\
        --protein-fasta P1.protein.fasta --patient P1 \\
        -o P1.long_peptides.tsv --fasta-out P1.long_peptides.fasta
"""
import argparse
import csv
import re
import sys
from pathlib import Path

# Kyte-Doolittle hydropathy, for GRAVY.
HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
POSITIVE, NEGATIVE = "KR", "DE"

OUT_COLUMNS = [
    "Peptide ID", "Patient", "MHC Class", "Gene", "AA Change", "Long Peptide",
    "Length", "Mutation Position", "Mutation Type", "Covers", "Contains Best Peptide",
    "Best Peptide", "Allele", "HLA LOH", "Agretopicity", "CCF", "Priority Rank",
    "Selection", "Flags", "Cys", "Met", "GRAVY", "Net Charge", "Longest Run", "Warnings",
]


def read_fasta(path):
    """-> [(header, sequence)] in file order."""
    entries, header, chunks = [], None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    entries.append((header, "".join(chunks)))
                header, chunks = line[1:], []
            else:
                chunks.append(line)
    if header is not None:
        entries.append((header, "".join(chunks)))
    return entries


def pair_entries(entries):
    """pVACtools MT/WT fasta -> [{'key', 'gene', 'mt', 'wt'}].

    Headers look like 'MT.<index>.<gene>.<transcript>.<change>' but the
    exact fields move between versions, so only the leading MT/WT and the
    next field (the index that pairs them) are relied on. Anything that
    looks like a gene symbol is kept as a hint for disambiguation.
    """
    pairs = {}
    for header, sequence in entries:
        fields = header.split(".")
        kind = fields[0].upper()
        if kind not in ("MT", "WT"):
            # Not a pVACtools MT/WT fasta: treat every record as a mutant
            # with no wild-type partner to compare against.
            pairs.setdefault(header, {"key": header, "gene": "", "mt": "", "wt": ""})["mt"] = sequence
            continue
        key = fields[1] if len(fields) > 1 else header
        gene = next((f for f in fields[2:] if re.fullmatch(r"[A-Z][A-Z0-9orf\-]{1,14}", f)), "")
        entry = pairs.setdefault(key, {"key": key, "gene": gene, "mt": "", "wt": ""})
        if gene and not entry["gene"]:
            entry["gene"] = gene
        entry["mt" if kind == "MT" else "wt"] = sequence
    return [entry for entry in pairs.values() if entry["mt"]]


def mutated_region(mt, wt):
    """-> (start, end, kind) over mt, 0-based and end-exclusive.

    Substitutions give the changed residues. Anything that shifts the
    frame or changes the length gives the novel stretch from the first
    difference to the end of the mutant sequence, which is what the immune
    system has not seen before.
    """
    if not wt:
        return 0, len(mt), "unknown"
    start = 0
    while start < min(len(mt), len(wt)) and mt[start] == wt[start]:
        start += 1
    if start == min(len(mt), len(wt)) and len(mt) == len(wt):
        return 0, len(mt), "none"          # identical: nothing to centre on
    if len(mt) != len(wt):
        return start, len(mt), "frameshift_or_indel"
    end = len(mt)
    while end > start and mt[end - 1] == wt[end - 1]:
        end -= 1
    return start, end, "substitution" if end - start == 1 else "multi_residue"


def window(sequence, mut_start, mut_end, max_length, anchor_start):
    """A window of up to max_length residues covering the mutation.

    Anchored at the mutation's start for a frameshift (the novel tail can
    be longer than any peptide worth making), otherwise centred on it and
    pushed inside the sequence at either terminus.
    """
    length = min(max_length, len(sequence))
    if anchor_start or (mut_end - mut_start) > length:
        start = mut_start
    else:
        centre = (mut_start + mut_end) / 2
        start = int(round(centre - length / 2))
    start = max(0, min(start, len(sequence) - length))
    return start, start + length


def trim_ends(sequence, start, end, mut_start, mut_end, bad_nterm, bad_cterm, min_flank, min_length):
    """Pull the ends in off awkward residues while keeping the mutation covered."""
    def flank_ok(new_start, new_end):
        if new_end - new_start < min_length:
            return False
        if new_start > mut_start or new_end < mut_end:
            return False
        return (mut_start - new_start) >= min_flank and (new_end - mut_end) >= min_flank

    changed = True
    while changed:
        changed = False
        if bad_nterm and sequence[start] in bad_nterm and flank_ok(start + 1, end):
            start, changed = start + 1, True
        if bad_cterm and sequence[end - 1] in bad_cterm and flank_ok(start, end - 1):
            end, changed = end - 1, True
    return start, end


def peptide_notes(peptide, max_gravy, max_run):
    """-> (metrics dict, [warnings])."""
    cys = peptide.count("C")
    met = peptide.count("M")
    known = [HYDROPATHY[a] for a in peptide if a in HYDROPATHY]
    gravy = sum(known) / len(known) if known else 0.0
    charge = sum(1 for a in peptide if a in POSITIVE) - sum(1 for a in peptide if a in NEGATIVE)
    longest = max((len(m.group(0)) for m in re.finditer(r"(.)\1*", peptide)), default=0)

    warnings = []
    if cys >= 2:
        warnings.append(f"{cys} cysteines (disulfide scrambling)")
    elif cys == 1:
        warnings.append("free cysteine (oxidation)")
    if met >= 2:
        warnings.append(f"{met} methionines (oxidation)")
    if gravy > max_gravy:
        warnings.append(f"GRAVY {gravy:.2f} (poor solubility)")
    if longest >= max_run:
        warnings.append(f"{longest}-residue run")
    if "X" in peptide or "*" in peptide:
        warnings.append("ambiguous or stop residue in sequence")
    return {"Cys": cys, "Met": met, "GRAVY": f"{gravy:.2f}", "Net Charge": charge,
            "Longest Run": longest}, warnings


def match_entry(candidate, pairs):
    """The fasta entry this candidate came from, or None.

    The epitope must be a substring of its own mutant protein, which is a
    stronger and more portable signal than any header convention; the gene
    symbol only breaks ties.
    """
    peptide = (candidate.get("Best Peptide") or "").strip().upper()
    gene = (candidate.get("Gene") or "").strip().upper()
    hits = [e for e in pairs if peptide and peptide in e["mt"].upper()] if peptide else []
    if not hits:
        hits = [e for e in pairs if gene and e["gene"].upper() == gene]
    if not hits:
        return None
    by_gene = [e for e in hits if gene and e["gene"].upper() == gene]
    return (by_gene or hits)[0]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("report", type=Path, help="ranked table from prioritize_candidates.py")
    parser.add_argument("--protein-fasta", type=Path, required=True,
                        help="output of `pvacseq generate_protein_fasta` for this patient")
    parser.add_argument("--patient", default="", help="patient id, used in the peptide names")
    parser.add_argument("--mhc-class", default="", help="label for the MHC Class column (e.g. I, II, Combined)")
    parser.add_argument("--max-length", type=int, default=25, help="longest peptide to design (default 25)")
    parser.add_argument("--min-length", type=int, default=15, help="shortest acceptable peptide after trimming (default 15)")
    parser.add_argument("--min-flank", type=int, default=3, help="residues of context to keep either side of the mutation (default 3)")
    parser.add_argument("--bad-nterm", default="QNCM", help="residues to trim from the N terminus where possible (default QNCM)")
    parser.add_argument("--bad-cterm", default="CMP", help="residues to trim from the C terminus where possible (default CMP)")
    parser.add_argument("--max-gravy", type=float, default=1.0, help="GRAVY above this is warned about (default 1.0)")
    parser.add_argument("--max-run", type=int, default=4, help="single-residue run at or above this is warned about (default 4)")
    parser.add_argument("--top", type=int, help="design at most this many peptides")
    parser.add_argument("--include-rejected", action="store_true", help="design peptides for rejected candidates too")
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument("--fasta-out", type=Path, help="also write the peptides as a fasta (for pVACvector)")
    args = parser.parse_args(argv)
    if args.min_length > args.max_length:
        parser.error("--min-length cannot exceed --max-length")

    pairs = pair_entries(read_fasta(args.protein_fasta)) if args.protein_fasta.exists() else []

    with args.report.open(newline="") as handle:
        candidates = list(csv.DictReader(handle, delimiter="\t"))
    if not args.include_rejected:
        candidates = [c for c in candidates if (c.get("Selection") or "").strip() != "Rejected"]

    designed, unmatched = [], []
    for candidate in candidates:
        if args.top and len(designed) >= args.top:
            break
        entry = match_entry(candidate, pairs)
        if entry is None:
            unmatched.append(candidate)
            continue
        mt = entry["mt"]
        mut_start, mut_end, kind = mutated_region(mt, entry["wt"])
        if kind == "none":
            unmatched.append(candidate)
            continue
        start, end = window(mt, mut_start, mut_end, args.max_length, kind == "frameshift_or_indel")
        start, end = trim_ends(mt, start, end, mut_start, min(mut_end, end), args.bad_nterm,
                               args.bad_cterm, args.min_flank, args.min_length)
        peptide = mt[start:end]

        # Merge a candidate whose mutation already sits inside a designed
        # peptide from the same protein: one peptide, both mutations named.
        merged = next((d for d in designed if d["entry_key"] == entry["key"]
                       and d["start"] <= mut_start and mut_end <= d["end"]), None)
        if merged:
            change = (candidate.get("AA Change") or candidate.get("ID") or "").strip()
            if change and change not in merged["covers"]:
                merged["covers"].append(change)
            continue

        metrics, warnings = peptide_notes(peptide, args.max_gravy, args.max_run)
        best = (candidate.get("Best Peptide") or "").strip().upper()
        contains = "yes" if best and best in peptide.upper() else "no"
        if contains == "no" and best:
            warnings.append("predicted epitope not fully inside the peptide")
        designed.append({
            "entry_key": entry["key"], "start": start, "end": end,
            "covers": [(candidate.get("AA Change") or candidate.get("ID") or "").strip()],
            "row": {
                "Peptide ID": "",
                "Patient": args.patient,
                "MHC Class": args.mhc_class,
                "Gene": candidate.get("Gene", ""),
                "AA Change": candidate.get("AA Change", ""),
                "Long Peptide": peptide,
                "Length": len(peptide),
                "Mutation Position": mut_start - start + 1,
                "Mutation Type": kind,
                "Covers": "",
                "Contains Best Peptide": contains,
                "Best Peptide": candidate.get("Best Peptide", ""),
                "Allele": candidate.get("Allele", ""),
                "HLA LOH": candidate.get("HLA LOH", ""),
                "Agretopicity": candidate.get("Agretopicity", ""),
                "CCF": candidate.get("CCF", ""),
                "Priority Rank": candidate.get("Priority Rank", ""),
                "Selection": candidate.get("Selection", ""),
                "Flags": candidate.get("Flags", ""),
                "Warnings": "; ".join(warnings),
                **metrics,
            },
        })

    for index, item in enumerate(designed, start=1):
        gene = item["row"]["Gene"] or "peptide"
        item["row"]["Peptide ID"] = "_".join(filter(None, [args.patient, f"LP{index:02d}", gene]))
        item["row"]["Covers"] = ";".join(item["covers"])

    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUT_COLUMNS, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for item in designed:
            writer.writerow(item["row"])

    if args.fasta_out:
        with args.fasta_out.open("w", newline="") as handle:
            for item in designed:
                handle.write(f">{item['row']['Peptide ID']}\n{item['row']['Long Peptide']}\n")

    print(f"{len(designed)} long peptides from {len(candidates)} candidates", file=sys.stderr)
    if unmatched:
        print(f"note: {len(unmatched)} candidates had no usable protein sequence in {args.protein_fasta.name}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
