#!/usr/bin/env python3
"""Generate deterministic synthetic tumour/normal DNA inputs; no patient data.

P1 enters as already-aligned BAMs with HLA alleles in the samplesheet, a
Sequenza segments file, an existing STAR-Fusion predictions table, and
baseline and week 6 RNA BAMs (copies of its
tumour BAM; the week 6 index is named <x>.bai rather than <x>.bam.bai); P2 enters as FASTQ pairs with no HLA alleles, so it
goes through alignment and HLA-HD, plus an RNA FASTQ pair so the fusion
branch runs STAR-Fusion for it. A tiny GTF and GTEx/TCGA expression tables
are written too, so the optional expression branch is wired in the stub
run, and placeholder AGFusion/pyensembl/CTAT paths so the fusion branch
passes its launch-time checks. Reference/resource files are empty placeholders: only valid with
nextflow -stub-run.
"""
import argparse
import csv
import gzip
import json
from pathlib import Path
import random
import shutil
import subprocess

PATIENTS = ("P1", "P2")
PAIRS = 24
SCATTER = ("chr1",)


def write_csv(path, fields, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_pairs(sample):
    """(name, start, mate) for each synthetic read pair."""
    return [(f"{sample}_pair{n:03d}", 200 + n * 7, 350 + n * 7) for n in range(PAIRS)]


def write_sam(path, reference, sample):
    """SAM SEQ is in alignment orientation, including reverse-strand mates."""
    records = []
    for name, start, mate in read_pairs(sample):
        span = mate + 100 - start
        records += [
            (start, f"{name}\t99\tchr1\t{start + 1}\t60\t100M\t=\t{mate + 1}\t{span}\t{reference[start:start + 100]}\t{'I' * 100}\tRG:Z:{sample}"),
            (mate, f"{name}\t147\tchr1\t{mate + 1}\t60\t100M\t=\t{start + 1}\t{-span}\t{reference[mate:mate + 100]}\t{'I' * 100}\tRG:Z:{sample}"),
        ]
    header = (f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{len(reference)}\n"
              f"@RG\tID:{sample}\tSM:{sample}\tLB:synthetic\tPL:ILLUMINA\tPU:CI\n")
    path.write_text(header + "\n".join(record for _, record in sorted(records)) + "\n", encoding="ascii")


def write_fastq(r1, r2, reference, sample):
    complement = str.maketrans("ACGT", "TGCA")
    with gzip.open(r1, "wt") as one, gzip.open(r2, "wt") as two:
        for name, start, mate in read_pairs(sample):
            one.write(f"@{name}/1\n{reference[start:start + 100]}\n+\n{'I' * 100}\n")
            two.write(f"@{name}/2\n{reference[mate:mate + 100].translate(complement)[::-1]}\n+\n{'I' * 100}\n")


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
    for name in ("known_dbsnp", "known_mills", "known_1000g", "germline_resource", "panel_of_normals", "contamination_sites"):
        (root / f"{name}.vcf.gz").write_bytes(b"")
        (root / f"{name}.vcf.gz.tbi").write_bytes(b"")
    for name in ("vep_cache", "vep_plugins"):
        (root / name).mkdir()
        (root / name / "STUB_ONLY.txt").write_text("Not a reference cache. Used only with nextflow -stub-run.\n")

    # Fusion branch inputs. P1 brings an existing STAR-Fusion table (so no
    # STAR-Fusion task); P2 brings an RNA FASTQ pair (so STAR-Fusion runs).
    # The reference paths below only have to satisfy main.nf's launch-time
    # checks: every fusion process is stubbed, so nothing reads them.
    sf_header = ("#FusionName\tJunctionReadCount\tSpanningFragCount\tSpliceType\tLeftGene\t"
                 "LeftBreakpoint\tRightGene\tRightBreakpoint\tLargeAnchorSupport\tFFPM\tannots\n")
    starfusion_tsv = root / "P1.star-fusion.fusion_predictions.tsv"
    starfusion_tsv.write_text(sf_header + "EML4--ALK\t12\t4\tONLY_REF_SPLICE\tEML4^ENSG1\tchr2:42000000:+\t"
                                          "ALK^ENSG2\tchr2:29000000:-\tYES\t1.20\t[\"INTERCHROMOSOMAL\"]\n")
    agfusion_db = root / "agfusion.homo_sapiens.111.db"
    agfusion_db.write_bytes(b"")
    pyensembl_cache = root / "pyensembl_cache"
    (pyensembl_cache / "pyensembl" / "GRCh38" / "ensembl111").mkdir(parents=True)
    (pyensembl_cache / "pyensembl" / "GRCh38" / "ensembl111" / "STUB_ONLY.db").write_bytes(b"")
    pyfix = root / "agfusion_pyfix"
    pyfix.mkdir()
    (pyfix / "sitecustomize.py").write_text("# stub shim; the real one patches numpy/datacache\n")
    ctat = root / "ctat_genome_lib"
    (ctat / "ref_genome.fa.star.idx").mkdir(parents=True)
    (ctat / "STUB_ONLY.txt").write_text("Not a CTAT library. Used only with nextflow -stub-run.\n")

    rows, inputs = [], []
    for patient in PATIENTS:
        row = {"patient_id": patient}
        for sample_type in ("tumor", "normal"):
            sample = f"{patient}_{sample_type}"
            if patient == "P1":
                sam, bam = root / f"{sample}.sam", root / f"{sample}.bam"
                write_sam(sam, reference, sample)
                if not sam_only:
                    subprocess.run([samtools, "sort", "-o", str(bam), str(sam)], check=True)
                    subprocess.run([samtools, "index", str(bam)], check=True)
                    subprocess.run([samtools, "quickcheck", "-v", str(bam)], check=True)
                    count = int(subprocess.check_output([samtools, "view", "-c", str(bam)], text=True))
                    if count != PAIRS * 2:
                        raise ValueError(f"{bam}: expected {PAIRS * 2} reads, got {count}")
                row[f"{sample_type}_bam"] = str(bam)
                inputs.append({"path": str(bam), "sample": sample, "kind": "bam", "reads": PAIRS * 2})
            else:
                r1, r2 = root / f"{sample}_R1.fastq.gz", root / f"{sample}_R2.fastq.gz"
                write_fastq(r1, r2, reference, sample)
                row[f"{sample_type}_r1"], row[f"{sample_type}_r2"] = str(r1), str(r2)
                inputs.append({"path": str(r1), "sample": sample, "kind": "fastq", "reads": PAIRS})
        if patient == "P1":
            for timepoint, index_name in (("rna_baseline", "P1_rna_baseline.bam.bai"), ("rna_week6", "P1_rna_week6.bai")):
                rna_bam = root / f"P1_{timepoint}.bam"
                if not sam_only:
                    shutil.copy(root / "P1_tumor.bam", rna_bam)
                    shutil.copy(root / "P1_tumor.bam.bai", root / index_name)
                row[f"{timepoint}_bam"] = str(rna_bam)
                inputs.append({"path": str(rna_bam), "sample": f"P1_{timepoint}", "kind": "rna_bam", "reads": PAIRS * 2})
            segments = root / "P1_segments.txt"
            segments.write_text('"chromosome"\t"start.pos"\t"end.pos"\t"CNt"\t"A"\t"B"\n"chr1"\t1\t10000\t3\t2\t1\n')
            row.update(purity="0.62", ploidy="3.1", sequenza_segments=str(segments), sex="female",
                       hla_alleles="HLA-A*02:01|HLA-B*07:02|DRB1*15:01",
                       starfusion_tsv=str(starfusion_tsv))
        else:
            # P2 has no STAR-Fusion table, so the fusion branch runs STAR-Fusion
            # on this pair. Same synthetic reads as its DNA FASTQ.
            rna_r1, rna_r2 = root / "P2_rna_baseline_R1.fastq.gz", root / "P2_rna_baseline_R2.fastq.gz"
            write_fastq(rna_r1, rna_r2, reference, "P2_rna_baseline")
            inputs.append({"path": str(rna_r1), "sample": "P2_rna_baseline", "kind": "fastq", "reads": PAIRS})
            row.update(purity="0.45", ploidy="", sequenza_segments="", sex="", hla_alleles="",
                       rna_baseline_r1=str(rna_r1), rna_baseline_r2=str(rna_r2))
        rows.append(row)

    # Expression branch inputs: one gene, one tissue, one cohort. Enough for
    # the wiring; the stub processes never read the values.
    (root / "genes.gtf").write_text(
        'chr1\ttest\tgene\t1\t1000\t.\t+\t.\tgene_id "ENSG00000141510.16"; '
        'gene_type "protein_coding"; gene_name "TP53";\n'
        'chr1\ttest\texon\t1\t1000\t.\t+\t.\tgene_id "ENSG00000141510.16"; '
        'gene_type "protein_coding"; gene_name "TP53";\n')
    header = "ENSEMBL_ID\tGENE_NAME\tGROUP\tQUANTILE\tTPM\n"
    (root / "gtex_db.tsv").write_text(
        header + "ENSG00000141510\tTP53\tSkin - Sun Exposed (Lower leg)\t0.5\t12.5000\n"
               + "ENSG00000141510\tTP53\tTestis\t0.5\t30.1000\n")
    (root / "tcga_db.tsv").write_text(
        header + "ENSG00000141510\tTP53\tTCGA-SKCM\t0.5\t15.0000\n")

    fields = ["patient_id", "tumor_r1", "tumor_r2", "normal_r1", "normal_r2", "tumor_bam", "normal_bam",
              "rna_baseline_bam", "rna_week6_bam", "rna_baseline_r1", "rna_baseline_r2", "starfusion_tsv",
              "purity", "ploidy", "sequenza_segments", "sex", "hla_alleles"]
    write_csv(root / "samplesheet.csv", fields, rows)
    manifest = {"synthetic": True, "test_type": "wiring_only", "patients": list(PATIENTS),
                "fastq_patients": ["P2"], "hla_typed_patients": ["P2"],
                "rna_bams": {"P1": ["rna_baseline", "rna_week6"]}, "scatter": list(SCATTER),
                # every fusion patient reaches AGFusion and pVACfuse; only the
                # one without a STAR-Fusion table reaches STAR-Fusion
                "fusion_patients": ["P1", "P2"], "starfusion_patients": ["P2"],
                "inputs": inputs, "bam_generation_complete": not sam_only}
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--sam-only", action="store_true", help="Developer inspection only; CI always builds/validates BAMs")
    args = parser.parse_args()
    result = generate(args.outdir, args.samtools, args.sam_only)
    print(f"Generated {len(result['inputs'])} synthetic read sets in {args.outdir}; BAMs complete: {result['bam_generation_complete']}")
