# Samplesheet schema

The pipeline takes **two** samplesheets, not one, because DNA and RNA don't
share a cardinality: WGS/WES is drawn once per patient at baseline, RNA-seq
is longitudinal (PRE / ED1 / ED2 / CLND, per `[[neoadjuvant-splicing]]` /
NeoTrio conventions). Forcing both into one row-per-timepoint sheet would
mean repeating the DNA columns on every RNA row or leaving them blank on
all but one — `patient_id` is the join key between the two.

## `--dna_samplesheet`

Splicing accepts one `spliceai_vcf` (already annotated) or `dna_vcf`
(somatic, PASS calls) per patient; `spliceai_vcf` takes precedence within
splicing only. Raw tumour/normal inputs work with
`--pipelines dna_variant_calling,splicing`. No HLA alleles are required for
splicing peptide-context generation. See [SPLICING.md](../docs/SPLICING.md).

The optional ERVcaller DNA branch also uses this sheet. Standalone
`--pipelines erv_dna` requires `patient_id,tumor_bam,normal_bam`, with
coordinate-sorted BAMs and adjacent `.bam.bai` or `.bai` indices. Patient IDs
must be unique and contain only letters, digits, underscores, hyphens and
dots, starting with a letter or digit. Both roles are required.
With `dna_variant_calling,erv_dna`, FASTQs can instead go through the DNA
branch and its resulting BAMs are reused. Rows containing `dna_vcf` must
still provide both indexed BAMs for ERVcaller; the VCF cannot substitute for
reads. See [DNA ERV inputs and references](../docs/ERV_DNA.md).

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

For splicing, supply `rna_bam` (coordinate-sorted, unsplit, with NH tags),
paired FASTQs, or `leafcutter_counts`. BAM rows require `rna_strand`:
`XS` (aligner tags), `RF` (first-strand) or `FR` (second-strand), with
`--splicing_rna_strand` as a cohort fallback. Existing LeafCutter matrices
use `leafcutter_sample` to select the exact sample column; it is required
for a multi-sample matrix. RNA VCFs cannot provide splice-junction evidence.
Every patient/timepoint pair must be unique and match a DNA patient.
Timepoint labels are arbitrary; no PRE/ED1/ED2/CLND ordering is enforced.

Same three tiers as DNA, per patient-timepoint: FASTQ, already-aligned BAM,
or already-called VCF.

| column        | required                     | notes                                                              |
|---------------|-------------------------------|---------------------------------------------------------------------|
| patient_id    | yes                            | join key back to the DNA sheet                                      |
| timepoint     | yes                            | any label, e.g. `PRE`, `day99`, `resection`; not restricted to a predefined sequence |
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

This same `--rna_samplesheet` also feeds `--run_fusion_neoantigens`
(`workflows/fusion_neoantigens.nf`): both `rna_bam` and `rna_fastq_r1`/
`rna_fastq_r2` rows work as input there — `rna_bam` rows go through their
own BAM→FASTQ step (reusing the `BAM_TO_FASTQ` module above) before being
handed to STAR-Fusion/Arriba, each of which then does its own dedicated STAR
alignment pass rather than reusing `RNA_VARIANT_CALLING`'s BAM. `rna_vcf`
rows are not usable for fusion calling (there's no variant file fusion
callers can start from) — a row with only `rna_vcf` set contributes to
`PVACSEQ_CORE` but is skipped by the fusion branch.

Same story again for `--run_erv_neoantigens` (`workflows/erv_neoantigens.nf`):
`rna_bam` (via `BAM_TO_FASTQ`) and `rna_fastq_r1`/`rna_fastq_r2` rows both
work, each realigned through the branch's own dedicated STAR pass (relaxed
multimapping, so ERV/TE loci's characteristically multi-mapping reads
survive) rather than reusing `RNA_VARIANT_CALLING`'s BAM. `rna_vcf`-only
rows are skipped here too, same reasoning as the fusion branch.

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

## `--hla_loh_samplesheet` (optional — only for standalone `--run_hla_loh`)

`--run_hla_loh` normally chains off `--run_dna_variant_calling`'s tumor BAM
and either `--purity_ploidy_csv` or `--run_purity_ploidy`. If you want to
run the HLA-LOH branch entirely on its own — no DNA/RNA branches, no
`dna_samplesheet`/`rna_samplesheet` at all — supply `--hla_loh_samplesheet`
instead and neither of the other two samplesheets is required.

| column     | required | notes                                                  |
|------------|----------|---------------------------------------------------------|
| patient_id | yes      | used as the run tag                                      |
| tumor_bam  | yes      | already-aligned tumor BAM (SpecHLA's LOH module is tumor-only — see workflows/hla_loh.nf); a `.bai` index next to it is assumed |
| purity     | yes      | tumor purity, e.g. from Sequenza's top solution           |
| ploidy     | yes      | tumor ploidy, e.g. from Sequenza's top solution            |

```csv
patient_id,tumor_bam,purity,ploidy
NeoTrio_003,/scratch/jo11/MORPHEUS_DNA/done/NeoTrio_003_tumor.bam,0.62,2.1
```

`main.nf`'s `validateHlaLohSamplesheet()` checks this at launch.

## `--purity_ploidy_csv` (optional — manual purity/ploidy when chaining off DNA_VARIANT_CALLING)

If you're running `--run_hla_loh` chained off `--run_dna_variant_calling`
but haven't (yet) got `--run_purity_ploidy` producing real output (it's
currently a stub — see `docs/ARCHITECTURE.md`), supply purity/ploidy
directly with this flag. Same shape as
[NeoadjLOH](https://github.com/JaydenBeckwith/NeoadjLOH)'s `PURITY_CSV` /
`sequenza_top_solutions_summary.csv`: comma-delimited, header
`sample,purity,ploidy`, where `sample` matches `patient_id` in the
`dna_samplesheet`.

```csv
sample,purity,ploidy
NeoTrio_003,0.62,2.1
```

## `--sequenza_gender_csv` (required for `--run_purity_ploidy`)

Not a pipeline samplesheet — this is your existing
`dna_neotrio_gender_metadata.csv`-shaped file, read exactly the way
`Sequenza_tools`' `sequenza_step2.sh`/`sequenza_step4.sh` do: a numeric
"melpin" id is pulled from the leading digits of each sample's
`patient_id`/`meta.id` and looked up against **column 2** of this CSV;
**column 7** is taken as the gender (lowercased). A header row is assumed
and skipped. A patient whose melpin has no match is skipped with a warning
rather than failing the run — same as the original scripts.

```csv
col1,melpin,col3,col4,col5,col6,gender
...,3,...,...,...,...,Female
```

(Column names above are illustrative — only the *position* of melpin
(2) and gender (7) matters, matching the original `awk -F',' '{...if($2==m)
print tolower($7)...}'` logic exactly.)
