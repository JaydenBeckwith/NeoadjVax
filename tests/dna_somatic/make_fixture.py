#!/usr/bin/env python3
"""Generate deterministic synthetic paired DNA/RNA reads; no patient data."""
import argparse
import csv
import json
from pathlib import Path
import random
import shutil
import subprocess

PATIENTS = ("P1", "P2")
TIMEPOINTS = ("PRE", "DAY42")
PAIRS = 24


def write_csv(path, fields, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_sam(path, reference, sample, rna=False):
    """SAM SEQ is in alignment orientation, including reverse-strand mates."""
    records = []
    for number in range(PAIRS):
        start = 200 + number * 7
        mate = start + (250 if rna else 150)
        span = mate + 100 - start
        seq = reference[start:start + 50] + reference[start + 150:start + 200] if rna else reference[start:start + 100]
        cigar = "50M100N50M" if rna else "100M"
        name = f"{sample}_pair{number:03d}"
        records += [
            (start, f"{name}\t99\tchr1\t{start + 1}\t60\t{cigar}\t=\t{mate + 1}\t{span}\t{seq}\t{'I' * 100}\tRG:Z:{sample}"),
            (mate, f"{name}\t147\tchr1\t{mate + 1}\t60\t100M\t=\t{start + 1}\t{-span}\t{reference[mate:mate + 100]}\t{'I' * 100}\tRG:Z:{sample}"),
        ]
    header = (f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{len(reference)}\n"
              f"@RG\tID:{sample}\tSM:{sample}\tLB:synthetic\tPL:ILLUMINA\tPU:CI\n")
    path.write_text(header + "\n".join(record for _, record in sorted(records)) + "\n", encoding="ascii")


def generate(root, samtools="samtools", sam_only=False):
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    if any(root.iterdir()):
        raise ValueError(f"Fixture output directory must be empty: {root}")
    if not sam_only and not shutil.which(samtools):
        raise ValueError("Install samtools or pass --samtools /path/to/samtools")
    rng = random.Random(104729)
    reference = "".join(rng.choice("ACGT") for _ in range(10000))
    (root / "genome.fa").write_text(">chr1\n" + "\n".join(reference[i:i + 80] for i in range(0, len(reference), 80)) + "\n")
    (root / "genome.gtf").write_text(
        'chr1\tsynthetic\texon\t201\t500\t.\t+\t.\tgene_id "SYNTHETIC"; transcript_id "SYNTHETIC_T1";\n'
        'chr1\tsynthetic\texon\t601\t1000\t.\t+\t.\tgene_id "SYNTHETIC"; transcript_id "SYNTHETIC_T1";\n')
    vcf_header = '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
    for name in ("known_mills", "known_1000g", "known_dbsnp"):
        (root / f"{name}.vcf").write_text(vcf_header)
    for name in ("vep_cache", "vep_plugins"):
        (root / name).mkdir()
        (root / name / "STUB_ONLY.txt").write_text("Not a reference cache. Used only with nextflow -stub-run.\n")
    dna_rows, rna_rows, bams = [], [], []
    for patient in PATIENTS:
        dna = {"patient_id": patient, "hla_alleles": "HLA-A*02:01|HLA-B*07:02" if patient == "P1" else "HLA-A*01:01|HLA-B*08:01"}
        for sample, rna in [(f"{patient}_tumor", False), (f"{patient}_normal", False)] + [(f"{patient}_{tp}", True) for tp in TIMEPOINTS]:
            sam = root / f"{sample}.sam"
            bam = root / f"{sample}.bam"
            write_sam(sam, reference, sample, rna=rna)
            if not sam_only:
                subprocess.run([samtools, "sort", "-o", str(bam), str(sam)], check=True)
                subprocess.run([samtools, "index", str(bam)], check=True)
                subprocess.run([samtools, "quickcheck", "-v", str(bam)], check=True)
                count = int(subprocess.check_output([samtools, "view", "-c", str(bam)], text=True))
                if count != PAIRS * 2:
                    raise ValueError(f"{bam}: expected {PAIRS * 2} reads, got {count}")
            bams.append({"path": str(bam), "sample": sample, "rna": rna, "reads": PAIRS * 2})
            if rna:
                rna_rows.append({"patient_id": patient, "timepoint": sample[len(patient) + 1:], "rna_bam": str(bam)})
            else:
                dna[sample[len(patient) + 1:] + "_bam"] = str(bam)
        dna_rows.append(dna)
    write_csv(root / "samples_dna.csv", ["patient_id", "tumor_bam", "normal_bam", "hla_alleles"], dna_rows)
    write_csv(root / "samples_rna.csv", ["patient_id", "timepoint", "rna_bam"], rna_rows)
    manifest = {"synthetic": True, "test_type": "wiring_only", "patients": list(PATIENTS),
                "rna_samples": rna_rows, "bams": bams, "bam_generation_complete": not sam_only}
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--sam-only", action="store_true", help="Developer inspection only; CI always builds/validates BAMs")
    args = parser.parse_args()
    result = generate(args.outdir, args.samtools, args.sam_only)
    print(f"Generated {len(result['bams'])} synthetic read sets in {args.outdir}; BAMs complete: {result['bam_generation_complete']}")
