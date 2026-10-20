# Somatic pipeline CI smoke test

This is a **fast wiring test**, not a variant-calling or neoantigen benchmark.

`make_fixture.py` generates two synthetic patients from a deterministic 10 kb reference. P1 enters as coordinate-sorted, indexed tumour/normal BAMs (24 read pairs each, checked with samtools), with HLA alleles, ploidy and a Sequenza segments file in the samplesheet. P2 enters as gzipped FASTQ pairs with no HLA alleles, so it goes through trimming, alignment and HLA-HD. No patient data or external downloads are used.

The real `main.nf` runs with `-stub-run`: every process's stub replaces the analysis command. The genome indexes are absent so the reference-building steps run too. The GATK resource VCFs, VEP directories and HLA-HD image name are placeholders and would not work outside a stub run.

`check_results.py` checks the task count for every stage (including that P1 skips alignment and HLA-HD while P2 does not), then that each patient has its PASS VCF, readcount-annotated pVACseq input, phased proximal VCF, pVACseq output and CCF table. A zero exit code alone is not enough.

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

Generation requires an empty destination. Reuse an existing generated fixture with `-resume`; do not rerun the generator over it. Generated data are ignored by Git. Downstream stub BAM/VCF files are placeholders and must never be read as scientific results.
