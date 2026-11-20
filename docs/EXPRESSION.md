# Expression context: GTEx, TCGA and the gene panel

An optional sixth step (`workflows/expression.nf`) with two outputs: the
candidate reports gain expression context, and each patient gets a panel
table for the genes that decide whether the tumour can present antigen at
all, plus the melanoma drivers.

The first half answers two questions pVACseq does not: is this candidate's
gene expressed in healthy tissue, and does this tumour type express it at
all?

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

## 4. The gene panel

`<id>/gene_panel/<id>.gene_panel.tsv`, one row per panel gene
(`bin/gene_panel_report.py`).

A ranked candidate list is worth little if the tumour has lost the
machinery to present it. B2M loss, HLA loss, JAK1/JAK2 truncation and
CIITA/NLRC5 silencing are recurrent routes to checkpoint-inhibitor
resistance, and none of them shows up in a per-candidate report, because
the defect is in the pathway rather than in any one candidate. The driver
half is context of a different kind: which canonical melanoma drivers this
tumour carries, and whether those genes are expressed.

For each gene the table combines three independent read-outs, any of which
may be absent:

| Read-out | From | Columns |
| --- | --- | --- |
| Expression | sample TPM vs the matched GTEx tissue and TCGA cohort | `Sample TPM`, `GTEx Range`, `GTEx Basis`, `GTEx Call`, `TCGA Range`, `TCGA Call` |
| Copy number | the Sequenza segment(s) covering the gene (needs `--rna_gtf` for coordinates) | `Tumor CN`, `Minor CN`, `CN Call` |
| Mutations | the VEP-annotated somatic VCF | `Mutations`, `Truncating` |

`GTEx Call` is `DOWN` / `NORMAL` / `UP` against the quantile range in your
table. The GTEx median-TPM summary has only medians, so with that file the
call falls back to a fold change either side of the median
(`--panel_fold`, default 4) and `GTEx Basis` says `median_x4` rather than
`q0.25-q0.75`. Build the table with `--quantiles 0.25 0.5 0.75` from a
per-sample matrix if you want real ranges.

`CN Call` is relative to the sample ploidy when the samplesheet gives it,
else 2: `LOSS` at a copy below baseline, `AMP` at `--panel_amp_gain`
copies above it, and `LOH` whenever the minor copy number is 0. A gene
straddling a breakpoint gets the overlap-weighted mean total copy number
and the lowest minor copy number across its segments, so half a gene
losing an allele still reads `LOH`.

### Which direction is adverse

Not every gene escapes by being lost, so each one carries an
`Adverse Direction` of `loss` or `gain`, and `Panel Flag` reads it:

| Flag | Means |
| --- | --- |
| `LOSS_OF_FUNCTION` | a `loss` gene that is `DOWN`, lost/LOH, or truncated |
| `INHIBITORY_GAIN` | a `gain` gene that is `UP` or amplified |
| `ALTERED` | a mutation, or a change in the non-adverse direction for this gene |
| `OK` | none of the above, with data present |
| `NA` | no data for this gene |

Most of the panel is `loss`: HLA-A/B/C, B2M, TAP, the immunoproteasome,
the regulators, the interferon genes, MR1, CD1 and the NKG2D ligands all
let a tumour escape by going away.

**HLA-E and HLA-G are the exception** and default to `gain`. They are
inhibitory ligands — HLA-E engages NKG2A on NK and CD8 T cells, HLA-G
engages ILT2/ILT4 and KIR2DL4 — so a tumour that *upregulates* them is
actively suppressing the response, and `UP` is the finding. Losing them is
not a defect, so it is recorded as `ALTERED` rather than dressed up as
one. The pairing to look for is a classical-locus `LOSS_OF_FUNCTION`
alongside a non-classical `INHIBITORY_GAIN`.

For the drivers neither direction is adverse: a BRAF V600E or NRAS Q61
hotspot is activating, not broken, so read the `Mutations` column rather
than the flag.

### The panel

94 genes in 11 pathways:

| Pathway | Genes |
| --- | --- |
| `apm_class_i` | HLA-A/B/C, B2M, TAP1/2, TAPBP, CALR, PDIA3 |
| `apm_non_classical` | HLA-E, HLA-F, HLA-G (class Ib), MR1, CD1A-E |
| `nkg2d_ligands` | MICA, MICB, ULBP1-3, RAET1E/G/L |
| `apm_class_ii` | CD74, HLA-DMA/DMB, HLA-DPA1/DPB1, HLA-DQA1/DQB1, HLA-DRA/DRB1 |
| `apm_regulation` | CIITA, NLRC5, RFX5, RFXANK, RFXAP, NFYA/B/C, CREB1 |
| `immunoproteasome` | PSMB5-10, PSME1-3, POMP, PSMF1 |
| `interferon_signalling` | IFNG, IFNGR1/2, JAK1, JAK2, STAT1, IRF1/2 |
| `driver_mapk` | BRAF, NRAS, KRAS, HRAS, NF1, MAP2K1/2, RAC1, PPP6C, KIT, GNAQ, GNA11 |
| `driver_cell_cycle_p53` | TP53, CDKN2A, CDK4, CCND1, RB1, MDM2 |
| `driver_pi3k` | PTEN, PIK3CA, AKT1 |
| `driver_other` | ARID2, IDH1, STK19, DDX3X, SF3B1, BAP1, EIF1AX, MITF, TERT, CTNNB1 |

The NKG2D ligands are MHC class I-*related*, not HLA: they are read by NK
and gamma-delta T cells rather than by CD8 T cells, and they are in the
panel because losing them is a separate escape route that the classical
locus does not show.

Replace the whole panel with `--gene_panel`, a
`gene<TAB>pathway[<TAB>loss|gain]` table; the third column is optional and
defaults to `loss`.

Caveats specific to this table:

- **HLA expression and copy number are the least trustworthy rows in
  it.** These are the most polymorphic genes in the genome, so short reads
  mismap there. Real HLA allelic loss needs a dedicated caller (LOHHLA,
  SpecHLA), not this. It applies to HLA-E/F/G and MICA/MICB too, which are
  polymorphic and sit in the same region.
- **MICA/MICB surface loss is mostly proteolytic shedding.** That is
  post-translational, so neither RNA nor DNA sees it: a normal row here
  does not mean the ligand is on the cell surface.
- **HLA-F's direction is unsettled** (open conformers, KIR3DS1), so it
  sits at `loss` by default; do not read its flag as firmly as the rest.
- **HLA-DRB3/4/5 cannot appear**: they are absent from GENCODE and GTEx.
- **An empty `TERT` mutation cell means nothing.** TERT's recurrent
  melanoma mutations are promoter, not coding, and a somatic exome call
  set usually does not cover them.
- **Low purity pushes genes towards `DOWN`**, because bulk expression is
  diluted by stroma and immune infiltrate.

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
- **The panel is not an HLA-LOH caller.** It sees gene-level copy number
  over the HLA locus, not allele-specific loss; that needs LOHHLA or
  SpecHLA, which live outside this branch.
- **No cross-timepoint expression change.** Week 6 TPM is published but not
  compared to baseline, so a gene the treatment switched off between
  baseline and week 6 is not called out.
