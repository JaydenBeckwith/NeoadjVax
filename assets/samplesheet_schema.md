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

| column      | required | notes                                                     |
|-------------|----------|-------------------------------------------------------------|
| patient_id  | yes      | join key, e.g. `NeoTrio_003`                                |
| tumor_r1    | if no tumor_bam  | tumour WGS/WES FASTQ, read 1                        |
| tumor_r2    | if no tumor_bam  | tumour WGS/WES FASTQ, read 2                        |
| normal_r1   | if no normal_bam | matched normal FASTQ, read 1                        |
| normal_r2   | if no normal_bam | matched normal FASTQ, read 2                        |
| tumor_bam   | if no tumor_r1/r2  | already-aligned tumour BAM — skips FastQC/Trim Galore/BWA-MEM entirely and goes straight into AddReadGroups+MarkDuplicates |
| normal_bam  | if no normal_r1/r2 | already-aligned normal BAM — same as above |
| hla_alleles | optional | pipe- or comma-separated; if blank, the HLA typing module runs on this patient's DNA |

FASTQ input (original behaviour):

```csv
patient_id,tumor_r1,tumor_r2,normal_r1,normal_r2,tumor_bam,normal_bam,hla_alleles
NeoTrio_003,/data/NeoTrio_003_tumor_R1.fastq.gz,/data/NeoTrio_003_tumor_R2.fastq.gz,/data/NeoTrio_003_normal_R1.fastq.gz,/data/NeoTrio_003_normal_R2.fastq.gz,,,
```

Already-aligned BAM input:

```csv
patient_id,tumor_r1,tumor_r2,normal_r1,normal_r2,tumor_bam,normal_bam,hla_alleles
NeoTrio_004,,,,,/data/NeoTrio_004_tumor.bam,/data/NeoTrio_004_normal.bam,
```

The BAM path assumes the input isn't yet read-group-tagged or
duplicate-marked (i.e. "aligned" not "GATK-ready") — it still runs
AddReadGroups + MarkDuplicates, just skips trimming and BWA-MEM. If you
also want to skip those, say so and we can add a second flag for
fully-processed BAMs.

## `--rna_samplesheet`

| column        | required                     | notes                                                              |
|---------------|-------------------------------|---------------------------------------------------------------------|
| patient_id    | yes                            | join key back to the DNA sheet                                      |
| timepoint     | yes                            | e.g. `PRE`, `ED1`, `ED2`, `CLND` — matches the ordering logic in the original RNA PBS script |
| rna_bam       | yes, unless rna_fastq_* given  | pre-aligned RNA-seq BAM, as consumed by the original RNA script — will be realigned with STAR two-pass |
| rna_fastq_r1  | alt. to rna_bam                | use if starting from FASTQ instead of BAM                          |
| rna_fastq_r2  | alt. to rna_bam                |                                                                       |

```csv
patient_id,timepoint,rna_bam
NeoTrio_003,PRE,/scratch/jo11/NeoTrio_RNA/CAGRF220510736_rnaseq/NeoTrio_003_PRE.bam
NeoTrio_003,ED1,/scratch/jo11/NeoTrio_RNA/CAGRF220510736_rnaseq/NeoTrio_003_ED1.bam
```

By default, `rna_bam` rows still get realigned through STAR two-pass —
this matches the original PBS script, which always re-derived FASTQs from
the input BAM and realigned rather than trusting the BAM's existing
alignment. Set `--rna_skip_realignment true` to skip BAM→FASTQ→STAR
entirely and feed `rna_bam` straight into AddReadGroups instead — only do
this if those BAMs are already aligned the way this pipeline expects
(STAR two-pass, `--outFilterMultimapNmax 2`); a BAM aligned with different
settings (or a different aligner) can give inconsistent variant calls
downstream. This is a pipeline-wide switch, not per-row.

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
