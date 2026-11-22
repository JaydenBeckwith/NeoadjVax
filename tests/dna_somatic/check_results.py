#!/usr/bin/env python3
"""Fail CI on missing/duplicate tasks, skipped bypasses or missing per-patient outputs."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path


def check(manifest, trace, outdir, analyses="somatic,fusion,markers"):
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
    # --run selects analyses; a deselected one must produce no tasks and
    # no outputs, which is as much a wiring claim as running does.
    selected = {a.strip() for a in analyses.split(",") if a.strip()}
    somatic, fusion_on, markers = ("somatic" in selected, "fusion" in selected, "markers" in selected)
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
        "HLA_TYPING:PARSE_HLA_LOH": len(manifest.get("hla_loh_patients", [])),
        "EXPRESSION:TAA_REPORT": patients,
        "PRIORITIZE": patients,
        "GENERATE_PROTEIN_FASTA": patients,
        "DESIGN_LONG_PEPTIDES": patients,
        "CLONAL_STRUCTURE": patients,
        "WORKBOOK": patients,
        "PVACVECTOR": patients,
        "PVACBIND": patients,
        "REPORT:RNA_SUPPORT": sum(1 for timepoints in rna.values() if timepoints),
        "REPORT:PATIENT_REPORT": patients,
        "REPORT:REPORT_PDF": patients,
    }
    # Which analysis each task belongs to. A deselected analysis must
    # produce no tasks at all, so its expected count becomes 0 rather than
    # the check being skipped.
    def needs(name):
        if name.startswith(("SOMATIC_CALLING:", "PROXIMAL_VARIANTS:", "NEOANTIGENS:")) or name in (
                "PRIORITIZE", "GENERATE_PROTEIN_FASTA", "DESIGN_LONG_PEPTIDES",
                "CLONAL_STRUCTURE", "PVACVECTOR", "REPORT:RNA_SUPPORT"):
            return somatic
        if name == "EXPRESSION:ANNOTATE_EXPRESSION":
            return somatic and markers      # annotates candidates, so it needs both
        if name.startswith("EXPRESSION:") or name == "PVACBIND":
            return markers
        if name.startswith("FUSION:"):
            return fusion_on
        return True                          # core: alignment, typing, workbook, report

    for name, count in expected.items():
        wanted = count if needs(name) else 0
        if observed[name] != wanted:
            raise ValueError(f"{name}: expected {wanted} tasks, got {observed[name]} (--run {analyses})")
    outdir = Path(outdir)
    for patient in manifest["patients"]:
        # Core outputs: produced whatever --run selected.
        paths = [
            outdir / patient / "report" / f"{patient}.report.tex",
            outdir / patient / "report" / f"{patient}.report.pdf",
            outdir / patient / f"{patient}.xlsx",
        ]
        somatic_paths = [
            outdir / patient / "somatic_variants" / f"{patient}.somatic.pass.vcf.gz",
            outdir / patient / "pvacseq_input" / f"{patient}.somatic.vep.readcounts.vcf.gz.tbi",
            outdir / patient / "pvacseq_input" / f"{patient}.phased_proximal.vcf.gz.tbi",
            outdir / patient / "pvacseq_ccf" / f"{patient}.Combined.aggregated.ccf.tsv",
            outdir / patient / "pvacseq" / "stub.log",
            outdir / patient / "prioritized" / f"{patient}.Combined.prioritized.tsv",
            outdir / patient / "long_peptides" / f"{patient}.protein.fasta",
            outdir / patient / "long_peptides" / f"{patient}.Combined.long_peptides.tsv",
            outdir / patient / "clonal_structure" / f"{patient}.clonal_structure.tsv",
            outdir / patient / "vaccine_vector" / "pvacvector" / "stub.log",
        ]
        marker_paths = [
            outdir / patient / "gene_panel" / f"{patient}.gene_panel.tsv",
            outdir / patient / "taa" / f"{patient}.taa.tsv",
            outdir / patient / "taa" / "pvacbind" / "stub.log",
        ]
        if somatic:
            paths += somatic_paths
        if markers:
            paths += marker_paths
        if somatic and markers:
            paths.append(outdir / patient / "pvacseq_expression" /
                         f"{patient}.Combined.aggregated.ccf.expression.tsv")
        if patient in manifest.get("hla_loh_patients", []):
            paths.append(outdir / patient / "hla_loh" / f"{patient}.hla_loh.tsv")
        if rna.get(patient) and somatic:
            paths.append(outdir / patient / "report" / f"{patient}.rna_support.summary.tsv")
            paths.append(outdir / patient / "report" / f"{patient}.rna_support.pdf")
        if patient in fusion and fusion_on:
            paths.append(outdir / patient / "fusion" / "agfusion" / "stub.log")
            for mhc in ("I", "II"):
                paths.append(outdir / patient / "pvacfuse" / f"MHC_Class_{mhc}" /
                             f"{patient}.MHC_{mhc}.all_epitopes.aggregated.tsv")
        if patient in starfusion and fusion_on:
            paths.append(outdir / patient / "fusion" / "star_fusion" / "star-fusion.fusion_predictions.tsv")
        if markers:
            for timepoint in rna.get(patient, []):
                paths.append(outdir / patient / "expression" / f"{patient}.{timepoint}.tpm.tsv")
        if patient in manifest["hla_typed_patients"]:
            paths.append(outdir / patient / "hla_typing" / f"{patient}.pvacseq_alleles.txt")
        for path in paths:
            if not path.is_file():
                raise ValueError(f"Missing per-patient output: {path}")
        if somatic:
            pvacseq_log = outdir / patient / "pvacseq" / "stub.log"
            if "no biological analysis" not in pvacseq_log.read_text():
                raise ValueError(f"Missing explicit stub marker: {pvacseq_log}")
            # the stub CCF table names the week 6 columns only when it was given that VCF
            ccf = outdir / patient / "pvacseq_ccf" / f"{patient}.Combined.aggregated.ccf.tsv"
            has_week6 = "RNA VAF (week 6)" in ccf.read_text()
            if has_week6 != ("rna_week6" in rna.get(patient, [])):
                raise ValueError(f"{ccf}: week 6 RNA columns {'present' if has_week6 else 'missing'} for {patient}")
        # a patient with its own STAR-Fusion table must not have run STAR-Fusion
        ran_starfusion = (outdir / patient / "fusion" / "star_fusion").is_dir()
        if ran_starfusion != (patient in starfusion and fusion_on):
            raise ValueError(f"{patient}: STAR-Fusion {'ran' if ran_starfusion else 'did not run'}, "
                             "which does not match the samplesheet's starfusion_tsv")
        has_support = (outdir / patient / "report" / f"{patient}.rna_support.summary.tsv").exists()
        if has_support != bool(rna.get(patient) and somatic):
            raise ValueError(f"{patient}: RNA-support output {'present without' if has_support else 'missing despite'} RNA")
    print(f"PASS: --run {analyses}: {patients} patients, "
          f"{'pVACseq, CCF, prioritisation and long peptides; ' if somatic else 'no somatic analysis; '}"
          f"{'expression and the gene panel; ' if markers else 'no expression; '}"
          f"{len(fusion) if fusion_on else 0} reached pVACfuse, and every patient got a report and a workbook; "
          "wiring only, no biological predictions.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--trace", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--analyses", default="somatic,fusion,markers",
                        help="what --run selected for this run; a deselected analysis "
                             "must have produced no tasks and no outputs")
    args = parser.parse_args()
    check(args.manifest, args.trace, args.outdir, args.analyses)
