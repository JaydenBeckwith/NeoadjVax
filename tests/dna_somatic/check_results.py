#!/usr/bin/env python3
"""Fail CI on missing/duplicate tasks, dropped timepoints or overwritten outputs."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path


def check(manifest, trace, outdir):
    manifest = json.loads(Path(manifest).read_text())
    if not manifest.get("synthetic") or not manifest.get("bam_generation_complete"):
        raise ValueError("Smoke test requires successfully generated/validated synthetic BAMs")
    with Path(trace).open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError("Nextflow trace is empty")
    failed = [row for row in rows if row["status"] not in ("COMPLETED", "CACHED") or row["exit"] != "0"]
    if failed:
        raise ValueError(f"Failed/incomplete tasks: {failed}")
    observed = Counter(row["process"] for row in rows)
    patients = len(manifest["patients"])
    samples = len(manifest["rna_samples"])
    expected = {
        "DNA_VARIANT_CALLING:ADD_READ_GROUPS": patients * 2,
        "DNA_VARIANT_CALLING:MARK_DUPLICATES": patients * 2,
        "DNA_VARIANT_CALLING:MUTECT2": patients,
        "DNA_VARIANT_CALLING:FILTER_MUTECT_CALLS": patients,
        "DNA_VARIANT_CALLING:FILTER_STANDARD_CHROMS": patients,
        **{f"RNA_VARIANT_CALLING:{name}": samples for name in (
            "ADD_READ_GROUPS_RNA", "MARK_DUPLICATES_RNA", "SPLIT_NCIGAR_READS",
            "BASE_RECALIBRATOR", "APPLY_BQSR", "HAPLOTYPE_CALLER_RNA",
            "VARIANT_FILTRATION_RNA", "VARIANTS_TO_TABLE")},
        **{f"PVACSEQ_CORE:{name}": samples for name in ("RNA_SUPPORT_FILTER", "VEP_ANNOTATE", "PVACSEQ_RUN")},
    }
    for name, count in expected.items():
        if observed[name] != count:
            raise ValueError(f"{name}: expected {count} tasks, got {observed[name]}")
    forbidden = {"BWA_MEM_ALIGN", "BAM_TO_FASTQ", "STAR_ALIGN", "HLA_TYPING_XHLA"}
    if any(name.split(":")[-1] in forbidden for name in observed):
        raise ValueError("Prealigned/manual-HLA bypass was not respected")
    outdir = Path(outdir)
    for sample in manifest["rna_samples"]:
        patient, tp = sample["patient_id"], sample["timepoint"]
        paths = [
            outdir / patient / "rna_supported_variants" / tp / f"{patient}_{tp}.rna_support_summary.tsv",
            outdir / patient / "dna_variants" / tp / f"{patient}_{tp}.annotated.vcf.gz",
            outdir / patient / "pvacseq" / tp / "pvacseq_output" / "stub.log",
        ]
        for path in paths:
            if not path.is_file():
                raise ValueError(f"Missing per-timepoint output: {path}")
        if "no biological analysis" not in paths[-1].read_text():
            raise ValueError(f"Missing explicit stub marker: {paths[-1]}")
    print(f"PASS: {patients} patients, {samples} RNA timepoints reached pVACseq; wiring only, no biological predictions.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--trace", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    check(args.manifest, args.trace, args.outdir)
