#!/usr/bin/env python3
"""Validated AGFusion -> pVACfuse hand-off. No network downloads or shell evaluation."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import re
import sqlite3
import subprocess


ALGORITHMS = ("MHCflurry", "NetMHCpan", "NetMHCIIpan")
REQUIRED = {
    "starfusion": {"LeftGene", "RightGene", "LeftBreakpoint", "RightBreakpoint",
                   "JunctionReadCount", "SpanningFragCount", "FFPM"},
    "arriba": {"gene1", "gene2", "gene_id1", "gene_id2", "breakpoint1", "breakpoint2",
               "split_reads1", "split_reads2", "discordant_mates"},
}


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def nonnegative(value, name, integer=False):
    try:
        result = int(value) if integer else float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name}: invalid number {value!r}") from exc
    if not math.isfinite(result) or result < 0:
        raise ValueError(f"{name}: expected a finite non-negative number")
    return result


def read_calls(source, caller, min_reads, min_ffpm):
    """Keep caller-native fields; record every exclusion rather than inventing evidence."""
    if caller not in REQUIRED:
        raise ValueError(f"Unsupported caller: {caller}")
    kept, audit = [], []
    with Path(source).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        original = reader.fieldnames
        if not original or len(original) != len(set(original)):
            raise ValueError("Caller TSV must have a non-empty, unique header")
        names = {field.lstrip("#"): field for field in original}
        if not REQUIRED[caller].issubset(names):
            raise ValueError(f"{caller}: missing columns {sorted(REQUIRED[caller] - names.keys())}")
        if caller == "starfusion" and not ({"#FusionName", "#fusion_name"} & set(original)):
            raise ValueError("STAR-Fusion requires #FusionName or #fusion_name")
        for number, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Caller TSV row {number}: wrong number of columns")
            data = {key.lstrip("#"): value for key, value in row.items()}
            bp_keys = ("LeftBreakpoint", "RightBreakpoint") if caller == "starfusion" else ("breakpoint1", "breakpoint2")
            for key in bp_keys:
                pattern = r"[^:\s]+:[1-9][0-9]*:[+-]" if caller == "starfusion" else r"[^:\s]+:[1-9][0-9]*"
                if not re.fullmatch(pattern, data[key]):
                    raise ValueError(f"Row {number}: invalid {key}: {data[key]!r}")
            gene_keys = ("LeftGene", "RightGene") if caller == "starfusion" else ("gene_id1", "gene_id2")
            for key in gene_keys:
                gene = data[key].split("^")[-1] if caller == "starfusion" else data[key]
                if not re.fullmatch(r"ENSG[0-9]+(?:\.[0-9]+)?", gene) or (caller == "starfusion" and "^" not in data[key]):
                    raise ValueError(f"Row {number}: {key} must identify one human Ensembl gene (ambiguous/multiple genes need review)")
            count_keys = ("JunctionReadCount", "SpanningFragCount") if caller == "starfusion" else ("split_reads1", "split_reads2", "discordant_mates")
            reads = sum(nonnegative(data[key], key, integer=True) for key in count_keys)
            ffpm = nonnegative(data["FFPM"], "FFPM") if caller == "starfusion" else None
            reasons = []
            if reads < min_reads:
                reasons.append("low_read_support")
            if ffpm is not None and ffpm < min_ffpm:
                reasons.append("low_ffpm")
            audit.append({"input_row": number, "breakpoint_5prime": data[bp_keys[0]],
                          "breakpoint_3prime": data[bp_keys[1]], "read_support": reads,
                          "ffpm": ffpm if ffpm is not None else "NA",
                          "decision": ";".join(reasons) if reasons else "retained"})
            if not reasons:
                kept.append(row)
    return original, kept, audit


def protein_records(directory):
    records = []
    for fasta in sorted(Path(directory).glob("*/*_protein.fa")):
        exon_base = str(fasta).replace("_protein.fa", ".exons")
        if not any(Path(exon_base + extension).is_file() for extension in (".txt", ".csv")):
            raise ValueError(f"Missing AGFusion exon table for {fasta}")
        header, sequence = None, []
        for line in fasta.read_text(encoding="utf-8").splitlines() + [">"]:
            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(sequence)
                    if sequence.count("*") != 1:
                        raise ValueError(f"{fasta}: expected exactly one --middlestar fusion marker in {header}")
                    records.append((header, sequence))
                header, sequence = line[1:], []
            elif line.strip():
                if header is None:
                    raise ValueError(f"Invalid FASTA: {fasta}")
                sequence.append(line.strip())
    return records


def run_logged(command, log, env=None):
    with Path(log).open("w", encoding="utf-8") as handle:
        result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, env=env, check=False)
    if result.returncode:
        raise RuntimeError(f"{command[0]} failed ({result.returncode}); see {log}")


def check_database(database):
    match = re.fullmatch(r"agfusion\.homo_sapiens\.([0-9]+)\.db", Path(database).name)
    if not match:
        raise ValueError("Keep the AGFusion database name agfusion.homo_sapiens.<release>.db")
    release = int(match[1])
    with sqlite3.connect(Path(database).resolve().as_uri() + "?mode=ro", uri=True) as connection:
        tables = {row[0] for row in connection.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    required = {f"homo_sapiens_{release}", f"homo_sapiens_{release}_transcript"}
    if not required.issubset(tables):
        raise ValueError(f"AGFusion SQLite tables do not match release {release}")
    return release


def annotate(args):
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    fields, calls, audit = read_calls(args.input, args.caller, args.min_reads, args.min_ffpm)
    filtered = output / "filtered_fusions.tsv"
    with filtered.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(calls)
    with (output / "input_audit.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=[
            "input_row", "breakpoint_5prime", "breakpoint_3prime", "read_support", "ffpm", "decision"])
        writer.writeheader()
        writer.writerows(audit)
    release = check_database(args.database)
    command = ["agfusion", "batch", "-f", str(filtered), "-a", args.caller,
               "-db", str(args.database), "-o", str(output), "--middlestar"]
    if args.noncanonical:
        command.append("--noncanonical")
    qc = {"caller": args.caller, "input_fusions": len(audit), "retained_fusions": len(calls),
          "min_read_support": args.min_reads, "min_ffpm": args.min_ffpm,
          "database": Path(args.database).name, "ensembl_release": release,
          "noncanonical": args.noncanonical, "command": command,
          "read_support_definition": "junction + spanning fragments" if args.caller == "starfusion" else "split_reads1 + split_reads2 + discordant_mates",
          "expression_available": args.caller == "starfusion"}
    if calls:
        cache = Path(args.cache).resolve()
        if not cache.is_dir() or not any(cache.rglob("*.db")):
            raise ValueError("PyEnsembl cache is missing its indexed reference; prepare it with pyensembl install")
        env = dict(os.environ, PYENSEMBL_CACHE_DIR=str(cache), MPLBACKEND="Agg")
        run_logged(command, output / "agfusion.log", env)
        # AGFusion catches gene/junction errors and otherwise exits zero. Do not
        # silently turn reference mismatches or skipped fusions into negatives.
        if re.search(r"\bERROR\b", (output / "agfusion.log").read_text(encoding="utf-8")):
            raise RuntimeError("AGFusion logged annotation errors; inspect agfusion.log and reconcile the references")
    proteins = protein_records(output)
    qc["protein_records"] = len(proteins)
    qc["status"] = ("annotated" if proteins else "no_protein_products") if calls else ("no_fusions_after_filter" if audit else "no_fusions")
    write_json(output / "annotation_qc.json", qc)
    return qc


def normalize_alleles(value):
    alleles = list(dict.fromkeys(item.strip() for item in re.split(r"[,|]", value) if item.strip()))
    if not alleles or any(not re.fullmatch(r"[A-Za-z0-9*:_-]+", allele) for allele in alleles):
        raise ValueError("Supply non-empty pVACtools-format HLA alleles")
    return alleles


def prediction_command(args):
    algorithms = [item.strip() for item in args.algorithms.split(",") if item.strip()]
    if not algorithms or len(algorithms) != len(set(algorithms)) or set(algorithms) - set(ALGORITHMS):
        raise ValueError(f"Supported offline predictors: {', '.join(ALGORITHMS)}; no duplicates")
    alleles = normalize_alleles(args.alleles)
    command = ["pvacfuse", "run", str(args.input), args.sample, ",".join(alleles),
               *algorithms, str(args.output), "-t", str(args.threads),
               "-e1", args.lengths_i, "-e2", args.lengths_ii,
               "--read-support", str(args.min_reads), "--expn-val", str(args.min_ffpm),
               "--binding-threshold", str(args.binding_threshold)]
    if args.caller == "starfusion":
        command += ["--starfusion-file", str(Path(args.input) / "filtered_fusions.tsv")]
    if args.iedb:
        command += ["--iedb-install-directory", str(args.iedb)]
    if set(algorithms) & {"NetMHCpan", "NetMHCIIpan"}:
        if not args.iedb:
            raise ValueError("NetMHC predictors require a local IEDB installation; online fallback is disabled")
        for relative in ("mhc_i/src/predict_binding.py", "mhc_ii/mhc_II_binding.py"):
            if not (Path(args.iedb) / relative).is_file():
                raise ValueError(f"Missing IEDB executable: {relative}")
    if "MHCflurry" in algorithms and (not args.models or not (Path(args.models) / "manifest.csv").is_file()):
        raise ValueError("Supply the downloaded MHCflurry affinity model directory containing manifest.csv")
    return command, algorithms, alleles


def validate_prediction_pairs(algorithms, alleles, lengths_i, lengths_ii):
    # Use the pinned pVACtools allele/length registries, not a hard-coded HLA list.
    from pvactools.lib import prediction_class
    supported, ignored = [], []
    for allele in alleles:
        matched = False
        for algorithm in algorithms:
            predictor = getattr(prediction_class, algorithm)()
            if allele in predictor.valid_allele_names():
                lengths = lengths_ii if algorithm == "NetMHCIIpan" else lengths_i
                valid = predictor.valid_lengths_for_allele(allele)
                if not set(map(int, lengths.split(","))).issubset(valid):
                    raise ValueError(f"{algorithm}/{allele}: lengths must be in {valid}")
                supported.append([algorithm, allele])
                matched = True
        if not matched:
            ignored.append(allele)
    if not supported:
        raise ValueError("None of the patient's HLA alleles is supported by the chosen predictors")
    return supported, ignored


def predict(args):
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    command, algorithms, alleles = prediction_command(args)
    if args.models:
        os.environ["MHCFLURRY_MODELS_DIR"] = str(Path(args.models).resolve())
    supported, ignored = validate_prediction_pairs(algorithms, alleles, args.lengths_i, args.lengths_ii)
    qc = {"sample": args.sample, "caller": args.caller, "command": command,
          "supported_algorithm_allele_pairs": supported, "unused_alleles": ignored,
          "read_support_in_pvacfuse_report": args.caller == "starfusion",
          "note": "Arriba read counts were filtered before AGFusion; pVACfuse's AGFusion importer reports Arriba read support and expression as NA. Do not interpret NA as passing expression evidence." if args.caller == "arriba" else ""}
    annotation = json.loads((Path(args.input) / "annotation_qc.json").read_text(encoding="utf-8"))
    if not protein_records(args.input):
        qc.update(status="no_predictable_fusion_proteins", annotation_status=annotation["status"])
    else:
        run_logged(command, output / "pvacfuse.log", os.environ.copy())
        reports = sorted(output.glob("MHC_Class_*/*.all_epitopes.tsv"))
        if not reports:
            # pVACfuse also exits zero when alleles/lengths or fusion sequences
            # cannot be processed. Require a real report for non-empty proteins.
            raise RuntimeError("pVACfuse produced no all_epitopes report; inspect pvacfuse.log")
        qc.update(status="completed", reports=[str(path.relative_to(output)) for path in reports])
    write_json(output / "run_status.json", qc)
    return qc


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    ag = sub.add_parser("annotate")
    ag.add_argument("--database", type=Path, required=True)
    ag.add_argument("--cache", type=Path, required=True)
    ag.add_argument("--noncanonical", action="store_true")
    pv = sub.add_parser("predict")
    pv.add_argument("--sample", required=True)
    pv.add_argument("--alleles", required=True)
    pv.add_argument("--algorithms", default="MHCflurry")
    pv.add_argument("--models", type=Path)
    pv.add_argument("--iedb", type=Path)
    pv.add_argument("--threads", type=int, default=4)
    pv.add_argument("--lengths-i", default="8,9,10,11")
    pv.add_argument("--lengths-ii", default="15")
    pv.add_argument("--binding-threshold", type=int, default=500)
    for child in (ag, pv):
        child.add_argument("--input", type=Path, required=True)
        child.add_argument("--output", type=Path, required=True)
        child.add_argument("--caller", choices=tuple(REQUIRED), required=True)
        child.add_argument("--min-reads", type=int, default=5)
        child.add_argument("--min-ffpm", type=float, default=0.1)
    args = parser.parse_args()
    nonnegative(args.min_reads, "min_reads", integer=True)
    nonnegative(args.min_ffpm, "min_ffpm")
    (annotate if args.action == "annotate" else predict)(args)


if __name__ == "__main__":
    main()
