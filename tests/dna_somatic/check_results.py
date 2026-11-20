#!/usr/bin/env python3
"""Fail CI on missing/duplicate tasks, skipped bypasses or missing per-patient outputs."""
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
    fastq = len(manifest["fastq_patients"])
    typed = len(manifest["hla_typed_patients"])
    scatter = len(manifest["scatter"])
    rna = manifest.get("rna_bams", {})
    fusion = manifest.get("fusion_patients", [])
    starfusion = manifest.get("starfusion_patients", [])
    week6 = sum("rna_week6" in timepoints for timepoints in rna.values())
    expected = {
        "PREPROCESS_DNA:TRIM_GALORE": fastq * 2,
        "PREPROCESS_DNA:BWA_MEM_ALIGN": fastq * 2,
        "PREPROCESS_DNA:ADD_READ_GROUPS": (patients - fastq) * 2,
        "PREPROCESS_DNA:MARK_DUPLICATES": patients * 2,
        "PREPROCESS_DNA:APPLY_BQSR": patients * 2,
        "SOMATIC_CALLING:MUTECT2": patients * scatter,
        "SOMATIC_CALLING:MERGE_MUTECT2_SCATTERS": patients,
        "SOMATIC_CALLING:GET_PILEUP_SUMMARIES": patients * 2,
        "SOMATIC_CALLING:CALCULATE_CONTAMINATION": patients,
        "SOMATIC_CALLING:FILTER_MUTECT_CALLS": patients,
        "SOMATIC_CALLING:SELECT_PASS_SOMATIC": patients,
        "HLA_TYPING:EXTRACT_HLA_READS": typed,
        "HLA_TYPING:HLAHD": typed,
        "PROXIMAL_VARIANTS:GERMLINE_HAPLOTYPECALLER": patients,
        "PROXIMAL_VARIANTS:WHATSHAP_PHASE": patients,
        "PROXIMAL_VARIANTS:VEP_ANNOTATE": patients,
        "NEOANTIGENS:VEP_ANNOTATE": patients,
        "NEOANTIGENS:BAM_READCOUNT": patients * 2 + sum(len(t) for t in rna.values()),
        "NEOANTIGENS:VCF_READCOUNT_ANNOTATOR": patients,
        "NEOANTIGENS:RNA_WEEK6_ANNOTATOR": week6,
        "NEOANTIGENS:PVACSEQ": patients,
        "NEOANTIGENS:ANNOTATE_CCF": patients,
        # one featureCounts/TPM task per RNA BAM; annotation runs for every
        # patient, with or without RNA
        "EXPRESSION:FEATURECOUNTS": sum(len(timepoints) for timepoints in rna.values()),
        "EXPRESSION:GENE_TPM": sum(len(timepoints) for timepoints in rna.values()),
        "EXPRESSION:ANNOTATE_EXPRESSION": patients,
        # the panel is per patient and needs neither RNA nor segments
        "EXPRESSION:GENE_PANEL": patients,
        # STAR-Fusion only for a patient with no existing table; AGFusion and
        # pVACfuse for every fusion patient
        "FUSION:STAR_FUSION": len(starfusion),
        "FUSION:AGFUSION": len(fusion),
        "FUSION:PVACFUSE": len(fusion),
        # RNA support only where there is RNA; a report (and PDF) for everyone
        "REPORT:RNA_SUPPORT": sum(1 for timepoints in rna.values() if timepoints),
        "REPORT:PATIENT_REPORT": patients,
        "REPORT:REPORT_PDF": patients,
    }
    for name, count in expected.items():
        if observed[name] != count:
            raise ValueError(f"{name}: expected {count} tasks, got {observed[name]}")
    outdir = Path(outdir)
    for patient in manifest["patients"]:
        paths = [
            outdir / patient / "somatic_variants" / f"{patient}.somatic.pass.vcf.gz",
            outdir / patient / "pvacseq_input" / f"{patient}.somatic.vep.readcounts.vcf.gz.tbi",
            outdir / patient / "pvacseq_input" / f"{patient}.phased_proximal.vcf.gz.tbi",
            outdir / patient / "pvacseq_ccf" / f"{patient}.Combined.aggregated.ccf.tsv",
            outdir / patient / "pvacseq" / "stub.log",
            outdir / patient / "pvacseq_expression" / f"{patient}.Combined.aggregated.ccf.expression.tsv",
            outdir / patient / "gene_panel" / f"{patient}.gene_panel.tsv",
            outdir / patient / "report" / f"{patient}.report.tex",
            outdir / patient / "report" / f"{patient}.report.pdf",
        ]
        if rna.get(patient):
            paths.append(outdir / patient / "report" / f"{patient}.rna_support.summary.tsv")
            paths.append(outdir / patient / "report" / f"{patient}.rna_support.pdf")
        if patient in fusion:
            paths.append(outdir / patient / "fusion" / "agfusion" / "stub.log")
            for mhc in ("I", "II"):
                paths.append(outdir / patient / "pvacfuse" / f"MHC_Class_{mhc}" /
                             f"{patient}.MHC_{mhc}.all_epitopes.aggregated.tsv")
        if patient in starfusion:
            paths.append(outdir / patient / "fusion" / "star_fusion" / "star-fusion.fusion_predictions.tsv")
        for timepoint in rna.get(patient, []):
            paths.append(outdir / patient / "expression" / f"{patient}.{timepoint}.tpm.tsv")
        if patient in manifest["hla_typed_patients"]:
            paths.append(outdir / patient / "hla_typing" / f"{patient}.pvacseq_alleles.txt")
        for path in paths:
            if not path.is_file():
                raise ValueError(f"Missing per-patient output: {path}")
        if "no biological analysis" not in paths[4].read_text():
            raise ValueError(f"Missing explicit stub marker: {paths[4]}")
        # the stub CCF table names the week 6 columns only when it was given that VCF
        has_week6 = "RNA VAF (week 6)" in paths[3].read_text()
        if has_week6 != ("rna_week6" in rna.get(patient, [])):
            raise ValueError(f"{paths[3]}: week 6 RNA columns {'present' if has_week6 else 'missing'} for {patient}")
        # a patient with its own STAR-Fusion table must not have run STAR-Fusion
        ran_starfusion = (outdir / patient / "fusion" / "star_fusion").is_dir()
        if ran_starfusion != (patient in starfusion):
            raise ValueError(f"{patient}: STAR-Fusion {'ran' if ran_starfusion else 'did not run'}, "
                             "which does not match the samplesheet's starfusion_tsv")
        has_support = (outdir / patient / "report" / f"{patient}.rna_support.summary.tsv").exists()
        if has_support != bool(rna.get(patient)):
            raise ValueError(f"{patient}: RNA-support output {'present without' if has_support else 'missing despite'} RNA")
    print(f"PASS: {patients} patients reached pVACseq, CCF, expression annotation, the gene panel, "
          f"{len(fusion)} reached pVACfuse, and every patient got a report; "
          "wiring only, no biological predictions.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--trace", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    check(args.manifest, args.trace, args.outdir)
