# Somatic pipeline CI smoke test

This is a **fast wiring test**, not a variant-calling or neoantigen benchmark.

`make_fixture.py` generates two synthetic patients from a deterministic 10 kb reference. P1 enters as coordinate-sorted, indexed tumour/normal BAMs (24 read pairs each, checked with samtools), with HLA alleles, ploidy, a Sequenza segments file and baseline and week 6 RNA BAMs (copies of the tumour BAM; the week 6 index uses the `<x>.bai` name) in the samplesheet. P2 enters as gzipped FASTQ pairs with no HLA alleles, so it goes through trimming, alignment and HLA-HD, and brings an RNA FASTQ pair so the fusion branch runs STAR-Fusion for it; P1 brings an existing STAR-Fusion table instead, so it skips straight to AGFusion. No patient data or external downloads are used.

A one-gene GTF and single-row GTEx/TCGA tables are written too, so the optional expression step is wired in as well, plus empty AGFusion/pyensembl/CTAT placeholders that exist only to satisfy the fusion branch's launch-time checks.

P1 also brings a SpecHLA `merge.hla.copy.txt` shaped like `cal.hla.copy.pl`'s output, with one locus lost (a zero side in `copyratio`, naming an allele P1's `hla_alleles` column also carries) and one retained; P2 brings none, so the run covers both the patient who has an HLA-LOH result and the patient who does not and must still reach the workbook. A two-peptide shared-antigen FASTA is written for pVACbind, and the smoke config turns on both opt-in branches (`--vaccine_vector`, `--taa_peptides`) so their wiring is covered too.

The real `main.nf` runs with `-stub-run`: every process's stub replaces the analysis command. The genome indexes are absent so the reference-building steps run too. The GATK resource VCFs, VEP directories and HLA-HD image name are placeholders and would not work outside a stub run.

`check_results.py` checks the task count for every stage (including that P1 skips alignment and HLA-HD while P2 does not), then that each patient has its PASS VCF, readcount-annotated pVACseq input, phased proximal VCF, pVACseq output, CCF table, expression-annotated report, gene-panel table and patient report (.tex and .pdf), that each patient also has its prioritised candidate table, protein fasta and long peptides, clonal-structure table, shared-antigen table, pVACbind and pVACvector outputs and `.xlsx` workbook, that only P1 has an HLA-LOH table, that P1's RNA timepoints each produced a TPM table and that only P1 (the patient with RNA) has the RNA-support table and figure, that only P1's CCF table carries the week 6 RNA columns, and that each fusion patient reached AGFusion and pVACfuse while only P2 ran STAR-Fusion. A zero exit code alone is not enough.

On Linux, with Python 3, samtools, Java 21 and Nextflow 26.04.6:

```bash
python tests/dna_somatic/make_fixture.py --outdir tests/dna_somatic/generated
export NXF_SYNTAX_PARSER=v1   # modules are not strict-syntax clean yet
nextflow -c tests/dna_somatic/smoke.config run main.nf -stub-run -work-dir test-work/dna_somatic
python tests/dna_somatic/check_results.py \
  --manifest tests/dna_somatic/generated/manifest.json \
  --trace test-results/dna_somatic/trace.tsv \
  --outdir test-results/dna_somatic
```

To check the `--run` selector, run the same fixture again with one analysis dropped and tell the checker what you asked for. A deselected analysis must produce no tasks and no outputs, which is as much a wiring claim as running it:

```bash
nextflow -c tests/dna_somatic/smoke.config run main.nf -stub-run --run markers \
  -work-dir test-work/dna_somatic_markers --outdir test-results/dna_somatic_markers
python tests/dna_somatic/check_results.py --analyses markers \
  --manifest tests/dna_somatic/generated/manifest.json \
  --trace test-results/dna_somatic_markers/trace.tsv \
  --outdir test-results/dna_somatic_markers
```

(The trace path comes from `smoke.config`, so point `-c` at a copy with a different `trace.file` if you want both runs' traces side by side.)

Generation requires an empty destination. Reuse an existing generated fixture with `-resume`; do not rerun the generator over it. Generated data are ignored by Git. Downstream stub BAM/VCF files are placeholders and must never be read as scientific results.
