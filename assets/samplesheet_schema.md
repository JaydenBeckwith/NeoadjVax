# Samplesheet schema

The pipeline takes **two** samplesheets, not one, because DNA and RNA don't
share a cardinality: WGS/WES is drawn once per patient at baseline, RNA-seq
is longitudinal (PRE / ED1 / ED2 / CLND, per `[[neoadjuvant-splicing]]` /
NeoTrio conventions). Forcing both into one row-per-timepoint sheet would
mean repeating the DNA columns on every RNA row or leaving them blank on
all but one — `patient_id` is the join key between the two.

## `--dna_samplesheet`

For each of tumor/normal, supply **either** a FASTQ pair **or** an
already-aligned BAM — not both. A row can mix them (e.g. normal already
aligned, tumor still raw FASTQ); each sample_type is handled independently.

Three tiers per patient, cheapest-to-supply first: raw FASTQ (full pipeline
runs), an already-aligned BAM (skips QC/trim/align), or an already-called
VCF (skips everything — the whole patient bypasses `DNA_VARIANT_CALLING`'s
compute). `dna_vcf` is a **patient-level** field (one somatic tumor-vs-normal
VCF, e.g. Mutect2 output), not per-tumor/normal like the other columns.

| column      | required | notes                                                     |
|-------------|----------|-------------------------------------------------------------|
| patient_id  | yes      | join key, e.g. `NeoTrio_003`                                |
| tumor_r1    | if no tumor_bam and no dna_vcf  | tumour WGS/WES FASTQ, read 1        |
| tumor_r2    | if no tumor_bam and no dna_vcf  | tumour WGS/WES FASTQ, read 2        |
| normal_r1   | if no normal_bam and no dna_vcf | matched normal FASTQ, read 1        |
| normal_r2   | if no normal_bam and no dna_vcf | matched normal FASTQ, read 2        |
| tumor_bam   | alt. to tumor_r1/r2, if no dna_vcf  | already-aligned tumour BAM — skips FastQC/Trim Galore/BWA-MEM entirely and goes straight into AddReadGroups+MarkDuplicates |
| normal_bam  | alt. to normal_r1/r2, if no dna_vcf | already-aligned normal BAM — same as above |
| dna_vcf     | optional, overrides all of the above | already-called somatic VCF (tumor vs normal) — skips alignment AND variant calling entirely; still passed through the standard-chroms filter for consistency. **If set and no normal_bam/normal_r1+r2 is also given, `hla_alleles` becomes required** — there's no normal BAM left to type HLA from. |
| hla_alleles | optional | `\|`-separated (a literal comma would be read as another CSV column) — converted internally to the comma-joined format `pvacseq run` expects; if blank, the HLA typing module runs on this patient's normal DNA (requires a normal BAM/FASTQ to exist) |

FASTQ input (original behaviour):

```csv
patient_id,tumor_r1,tumor_r2,normal_r1,normal_r2,tumor_bam,normal_bam,dna_vcf,hla_alleles
NeoTrio_003,/data/NeoTrio_003_tumor_R1.fastq.gz,/data/NeoTrio_003_tumor_R2.fastq.gz,/data/NeoTrio_003_normal_R1.fastq.gz,/data/NeoTrio_003_normal_R2.fastq.gz,,,,
```

Already-aligned BAM input:

```csv
patient_id,tumor_r1,tumor_r2,normal_r1,normal_r2,tumor_bam,normal_bam,dna_vcf,hla_alleles
NeoTrio_004,,,,,/data/NeoTrio_004_tumor.bam,/data/NeoTrio_004_normal.bam,,
```

Already-called VCF input (needs `hla_alleles` since no normal is given here
— use `|` to separate multiple alleles within the one CSV field, since a
bare comma would be parsed as another column):

```csv
patient_id,tumor_r1,tumor_r2,normal_r1,normal_r2,tumor_bam,normal_bam,dna_vcf,hla_alleles
NeoTrio_005,,,,,,,/data/NeoTrio_005_somatic.vcf.gz,HLA-A*02:01|HLA-A*24:02|HLA-B*07:02|HLA-B*15:01|HLA-C*03:04|HLA-C*07:02
```

The BAM path assumes the input isn't yet read-group-tagged or
duplicate-marked (i.e. "aligned" not "GATK-ready") — it still runs
AddReadGroups + MarkDuplicates, just skips trimming and BWA-MEM. If you
also want to skip those, say so and we can add a further flag for
fully-processed BAMs.

`main.nf`'s `validateDnaSamplesheet()` checks all of this at launch (before
any compute starts) and errors out with a specific message if a patient
row is missing everything a given stage needs.

## `--rna_samplesheet`

Same three tiers as DNA, per patient-timepoint: FASTQ, already-aligned BAM,
or already-called VCF.

| column        | required                     | notes                                                              |
|---------------|-------------------------------|---------------------------------------------------------------------|
| patient_id    | yes                            | join key back to the DNA sheet                                      |
| timepoint     | yes                            | e.g. `PRE`, `ED1`, `ED2`, `CLND` — matches the ordering logic in the original RNA PBS script |
| rna_bam       | if no rna_fastq_* and no rna_vcf | pre-aligned RNA-seq BAM, as consumed by the original RNA script — realigned with STAR two-pass unless `--rna_skip_realignment true` |
| rna_fastq_r1  | alt. to rna_bam, if no rna_vcf | use if starting from FASTQ instead of BAM                          |
| rna_fastq_r2  | alt. to rna_bam, if no rna_vcf |                                                                       |
| rna_vcf       | optional, overrides all of the above | already-called+filtered RNA VCF for this patient-timepoint — skips BAM→FASTQ→STAR AND the GATK variant-calling chain entirely |

```csv
patient_id,timepoint,rna_bam,rna_vcf
NeoTrio_003,PRE,/scratch/jo11/NeoTrio_RNA/CAGRF220510736_rnaseq/NeoTrio_003_PRE.bam,
NeoTrio_003,ED1,/scratch/jo11/NeoTrio_RNA/CAGRF220510736_rnaseq/NeoTrio_003_ED1.bam,
NeoTrio_006,PRE,,/data/NeoTrio_006_PRE.filtered.vcf
```

By default, `rna_bam` rows still get realigned through STAR two-pass —
this matches the original PBS script, which always re-derived FASTQs from
the input BAM and realigned rather than trusting the BAM's existing
alignment. Set `--rna_skip_realignment true` to skip BAM→FASTQ→STAR
entirely and feed `rna_bam` straight into AddReadGroups instead — only do
this if those BAMs are already aligned the way this pipeline expects
(STAR two-pass, `--outFilterMultimapNmax 2`); a BAM aligned with different
settings (or a different aligner) can give inconsistent variant calls
downstream. This is a pipeline-wide switch, not per-row — `rna_vcf` is the
per-row way to skip alignment+calling entirely for a specific timepoint.

`main.nf`'s `validateRnaSamplesheet()` checks all of this at launch, same
as the DNA side.

## How they're joined for `PVACSEQ_CORE`

For each RNA row, the pipeline looks up that `patient_id`'s DNA filtered
VCF and checks which somatic loci carry RNA support **at that timepoint**.
So a patient with one baseline DNA draw and 3 RNA timepoints produces 3
RNA-supported neoantigen calls (one per timepoint), each checked against
the same baseline DNA variant set — matching the baseline-vs-week-6 (and
beyond) comparison Ines asked about.

This is a starting proposal, not fixed — flag if any patients have more
than one DNA draw (e.g. a resection sample later in treatment), since the
join logic above assumes exactly one DNA VCF per patient.
