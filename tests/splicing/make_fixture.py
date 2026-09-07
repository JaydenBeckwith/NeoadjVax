#!/usr/bin/env python3
"""Generate tiny, synthetic splicing inputs. No patient data or real model output."""
import argparse
import gzip
import json
from pathlib import Path


def make_fixture(target):
    target.mkdir(parents=True, exist_ok=True)
    seq = list("A" * 1300)
    for start, text in [(101, "ATG" * 20), (301, "GCT" * 20), (501, "TTT" * 19 + "TAA"),
                        (701, "TTA" + "AAA" * 19), (901, "AGC" * 20), (961, "CCCCCC"), (1101, "CAT" * 20)]:
        seq[start - 1:start - 1 + len(text)] = text
    (target / "genome.fa").write_text(">chr1\n" + "".join(seq) + "\n")
    gtf = []
    for gene, strand, ranges in [("GENEP", "+", [(101, 160), (301, 360), (501, 560)]),
                                  ("GENEN", "-", [(701, 760), (901, 960), (1101, 1160)])]:
        for start, end in ranges:
            attrs = f'gene_id "{gene}"; gene_name "{gene}"; transcript_id "{gene}_TX";'
            for feature, phase in [("exon", "."), ("CDS", "0")]:
                gtf.append(f"chr1\tsynthetic\t{feature}\t{start}\t{end}\t.\t{strand}\t{phase}\t{attrs}")
    (target / "genes.gtf").write_text("\n".join(gtf) + "\n")
    vcf = ['##fileformat=VCFv4.2', '##contig=<ID=chr1,length=1300>',
           '##INFO=<ID=SpliceAI,Number=.,Type=String,Description="Synthetic SpliceAI fixture, not model predictions">',
           '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO']
    for pos, gene, offset in [(290, "GENEP", 5), (291, "GENEP", 7), (292, "GENEP", 9), (1000, "GENEN", -34)]:
        vcf.append(f"chr1\t{pos}\t.\tA\tG\t.\tPASS\tSpliceAI=G|{gene}|0.9|0.1|0|0|{offset}|0|0|0")
    (target / "dna.spliceai.vcf").write_text("\n".join(vcf) + "\n")
    (target / "normal.tsv").write_text("junc_id\nchr1:160-298:+\n")
    with gzip.open(target / "cohort_perind.counts.gz", "wt") as out:
        out.write("chrom P1_day0 P1_day99 P2_resection\n")
        for junc, counts in [
            ("chr1:160:295:clu_1_+", "5/60 0/60 5/60"),
            ("chr1:160:298:clu_1_+", "5/60 5/60 5/60"),
            ("chr1:160:301:clu_1_+", "50/60 55/60 50/60"),
            ("chr1:966:1101:clu_2_-", "5/50 5/50 5/50")]:
            out.write(junc + " " + counts + "\n")
    (target / "dna.csv").write_text("patient_id,spliceai_vcf\n" + "".join(
        f"{patient},{(target / 'dna.spliceai.vcf').as_posix()}\n" for patient in ["P1", "P2"]))
    (target / "rna.csv").write_text("patient_id,timepoint,leafcutter_counts,leafcutter_sample\n" + "".join(
        f"{patient},{tp},{(target / 'cohort_perind.counts.gz').as_posix()},{sample}\n"
        for patient, tp, sample in [("P1", "day0", "P1_day0"), ("P1", "day99", "P1_day99"), ("P2", "resection", "P2_resection")]))
    params = dict(pipelines="splicing", dna_samplesheet=str(target / "dna.csv"),
                  rna_samplesheet=str(target / "rna.csv"), genome_fasta=str(target / "genome.fa"),
                  gtf=str(target / "genes.gtf"), splicing_normal_junctions=str(target / "normal.tsv"),
                  outdir=str(target / "results"))
    (target / "params.json").write_text(json.dumps(params, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output")
    make_fixture(Path(parser.parse_args().output).resolve())
