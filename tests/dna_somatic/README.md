# Somatic pipeline CI smoke test

This is a **fast wiring test**, not a variant-calling or neoantigen benchmark.

`make_fixture.py` generates eight valid, coordinate-sorted, indexed BAMs from deterministic synthetic reads: two patients, each with tumour/normal DNA and RNA at PRE and DAY42. Each BAM has 24 read pairs; RNA includes spliced CIGARs. Samtools checks BAM integrity and read counts. No patient data or external genome downloads are used.

The real `main.nf` runs with `--pipelines somatic_neoantigen -stub-run`. Its existing process stubs replace the analysis commands. Prealigned BAMs bypass alignment, and supplied HLA alleles bypass xHLA. The small GTF/known-sites files and VEP directories are structural placeholders, not usable biological references.

`check_results.py` checks successful task counts through DNA calling, RNA calling, RNA-support matching, VEP and pVACseq, then asserts distinct published outputs for **all four RNA samples**. A zero exit code alone is insufficient. CI uploads generated inputs, results, traces and task diagnostics even on failure.

On Linux, with Python 3, samtools, Java 21 and Nextflow 26.04.6:

```bash
python tests/dna_somatic/make_fixture.py --outdir tests/dna_somatic/generated
nextflow -c tests/dna_somatic/smoke.config run main.nf \
  --pipelines somatic_neoantigen -stub-run -work-dir test-work/dna_somatic
python tests/dna_somatic/check_results.py \
  --manifest tests/dna_somatic/generated/manifest.json \
  --trace test-results/dna_somatic/trace.tsv \
  --outdir test-results/dna_somatic
```

Generation requires an empty destination. Reuse an existing generated fixture with `-resume`; do not rerun the generator over it. Generated data are ignored by Git. Input BAMs are real, but downstream stub BAM/VCF files are placeholders and must never be interpreted as scientific results.
