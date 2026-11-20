# Expression context: GTEx and TCGA

An optional sixth step (`workflows/expression.nf`) that answers two
questions pVACseq does not: is this candidate's gene expressed in healthy
tissue, and does this tumour type express it at all?

Both matter for a vaccine target:

- **Normal-tissue expression (GTEx).** A candidate in a gene that healthy
  tissues express is a worse target. The wild-type context of the same
  peptide is presented by normal cells, so T cells against it are more
  likely to have been tolerised, and on-target off-tumour toxicity is a
  real risk. A gene silent outside the tumour (or outside testis, which is
  immune-privileged) is the better target.
- **Cohort expression (TCGA).** A gene the matched cohort expresses highly
  is more plausibly expressed in this patient than one the cohort never
  expresses. This is weak evidence next to the patient's own RNA, but it is
  the only expression evidence available when RNA is missing, so the report
  falls back to it and says so in the `Expression Evidence` column.

It runs after `ANNOTATE_CCF`, so one table carries clonality and
expression together. The two halves are independent: annotation needs only
a reference table and works on a DNA-only patient; quantification needs an
RNA BAM and a GTF.

No reference data ships with this repository. You download GTEx and TCGA
yourself and build the tables once.

## 1. Build the reference tables

`bin/build_expression_db.py` writes a single long-format TSV
(`ENSEMBL_ID  GENE_NAME  GROUP  QUANTILE  TPM`), where `GROUP` is a GTEx
tissue or a TCGA cohort.

### GTEx

The public median-TPM summary (one median per gene per tissue, ~15 MB) is
enough, and is the simplest thing to start from. Download
`GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz` from
the [GTEx portal downloads page](https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression),
then:

```bash
python bin/build_expression_db.py \
    --from-gct GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz \
    -o reference/gtex_v8_median_tpm.tsv
```

Medians only, so `--quantiles` does not apply. For percentiles (e.g. a 95th
percentile, to ask "how high does this gene ever go in normal tissue?"),
use the per-sample GTEx TPM matrix with `--from-matrix` and a
sample-to-tissue table built from the GTEx sample attributes file.

### TCGA

There is no equally tidy official per-cohort median file, so give a
gene-by-sample TPM matrix plus a two-column `sample<TAB>group` table.
A UCSC Xena / Toil recompute matrix works, and its values are
`log2(TPM+1)`, hence `--expm1`. Cohort matrices usually label genes by
Ensembl ID only, so pass the GTF to fill in symbols (pVACseq reports
symbols):

```bash
python bin/build_expression_db.py \
    --from-matrix tcga_gene_tpm.tsv \
    --groups tcga_sample_cohort.tsv \
    --gene-names gencode.v44.annotation.gtf.gz \
    --expm1 --quantiles 0.5 \
    -o reference/tcga_cancertype_tpm.tsv
```

For melanoma the group to report is the SKCM cohort; the exact spelling
(`TCGA-SKCM`, `SKCM`, `Skin Cutaneous Melanoma`) comes from whichever
phenotype table you used, and `--tcga_cancer_type` must match it. Pass a
name that is not in the table and the run stops with the available names
listed, rather than filling the column with `NA`.

The matrix must already be TPM (not counts): nothing is re-normalised.

## 2. Run it

```bash
nextflow run main.nf \
    --samplesheet samplesheets/pilot.csv \
    --rna_gtf /path/to/gencode.v44.annotation.gtf.gz \
    --rna_strandedness 2 \
    --gtex_db reference/gtex_v8_median_tpm.tsv \
    --gtex_tissue 'Skin - Sun Exposed (Lower leg)' \
    --tcga_db reference/tcga_cancertype_tpm.tsv \
    --tcga_cancer_type 'TCGA-SKCM' \
    ...
```

`--rna_strandedness` mirrors featureCounts `-s`: 0 unstranded, 1 forward,
2 reverse. It has no safe default. Most current Illumina kits (TruSeq
Stranded mRNA, NEBNext Ultra II Directional) are reverse-stranded, so 2 is
the usual answer, but set it from the kit: counting a reverse-stranded
library as forward loses most reads, and the pipeline cannot tell that this
happened, it just reports near-zero expression everywhere.

Without `--rna_gtf`, quantification is skipped and expression comes from
the reference tables alone (the run logs a warning). Without either
`--gtex_db` or `--tcga_db`, the whole step is skipped and the CCF reports
are the final output.

## 3. What comes out

`<id>/expression/<id>.<timepoint>.tpm.tsv` — gene-level TPM for every RNA
BAM given, from featureCounts counts and union-exon lengths
(`bin/gene_tpm.py`). `EXPRESSION_QUARTILE` ranks the gene among the genes
expressed in *this* sample (0 = not expressed), so it is a within-sample
rank, not a cross-sample one. Both timepoints are quantified; the week 6
table is published for comparison but the candidate annotation uses
baseline, the sample the candidates were called from.

`<id>/pvacseq_expression/<id>.<class>.aggregated.ccf.expression.tsv` — the
CCF report plus:

| Column | Meaning |
| --- | --- |
| `Sample TPM`, `Sample Expression Quartile` | this patient's baseline RNA (`NA` without RNA) |
| `GTEx <tissue> TPM` | the matched normal tissue |
| `GTEx Median TPM`, `GTEx Max TPM` | across every tissue in the table |
| `GTEx Tissues Expressed` | how many tissues exceed `--gtex_expressed_tpm` |
| `GTEx Specificity` | `not_expressed` / `testis_restricted` / `tissue_restricted` / `broadly_expressed` |
| `TCGA <cohort> TPM` | the matched tumour cohort |
| `Expression Evidence` | `sample_rna`, `tcga_cohort` or `none`: which value the flag used |
| `Expression Flag` | `OK` / `LOW_EXPRESSION` / `NO_EXPRESSION` / `NO_DATA` |

`GTEx Specificity` treats testis separately on purpose: a cancer-testis
antigen such as a CTAG/MAGE family gene is silent in every somatic tissue
but high in testis, and averaging tissues together would hide exactly the
genes worth looking at.

## What this is not

- **The flag is a sort key, not a verdict.** Thresholds are arguments
  (`--expression_min_tpm`, `--gtex_expressed_tpm`), and a flagged candidate
  is annotated, never dropped. `NO_EXPRESSION` means this sample's RNA saw
  the gene at zero TPM, which with shallow RNA is also what a lowly
  expressed gene looks like.
- **Gene level, not allele or transcript level.** A gene can be expressed
  while the mutant allele is not, which is what pVACseq's own RNA VAF
  columns are for, and a gene can be expressed from a transcript that does
  not carry the variant. Neither is corrected for here.
- **Bulk, not tumour-cell-specific.** TPM from a bulk RNA BAM mixes tumour
  and stroma, so low purity drags apparent tumour expression down.
- **No HLA or antigen-presentation panel yet.** NeoDisc also calls each
  antigen-processing gene (HLA-A/B/C, B2M, TAP1/2, the immunoproteasome,
  CIITA) up or down against the same GTEx/TCGA ranges and pairs that with
  copy number, which is how a tumour that cannot present anything gets
  spotted. That is a natural next step here, and it would pair with the
  HLA-LOH work that lives outside this branch.
- **No cross-timepoint expression change.** Week 6 TPM is published but not
  compared to baseline.
