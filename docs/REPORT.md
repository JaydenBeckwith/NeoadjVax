# Per-patient report and RNA-support figure

The last step (`workflows/report.nf`) turns each patient's outputs into
two things: a figure showing how many mutations reach tumour RNA at each
timepoint, and a melanoma report modelled on the NeoDisc LaTeX templates.
It runs by default; turn it off with `--report false`, or keep the `.tex`
and skip PDF compilation with `--report_pdf false`.

Everything here is a computational read-out for research. The report says
so on its first page.

## Where the figures come from

`patient_report.py` writes LaTeX and draws nothing. Every figure is made
by a separate script and handed to the report as a file:

| figure | script | from |
|---|---|---|
| RNA support across timepoints | `bin/rna_support.py` | readcount VCFs, CCF tables |
| Candidate landscape | `bin/report_figures.py` | the prioritised candidate table |
| Flag summary | `bin/report_figures.py` | the prioritised candidate table |
| HLA loss of heterozygosity | `bin/report_figures.py` | the tidy SpecHLA table |

That split is the point: a figure can be regenerated from an existing
results directory and looked at on its own, without pdflatex or the rest
of the report, and the data behind each one is a plain function that is
tested without matplotlib installed at all.

```bash
# redraw one patient's prediction figures from results already on disk
python bin/report_figures.py --prefix /tmp/P1     --prioritized results/P1/prioritized/P1.Combined.aggregated.ccf.expression.prioritized.tsv     --hla-loh results/P1/hla_loh/P1.hla_loh.tsv

# just one of them, PNG only, while iterating on how it looks
python bin/report_figures.py --prefix /tmp/P1 --only candidates --format png     --prioritized results/P1/prioritized/P1.Combined.aggregated.ccf.expression.prioritized.tsv
```

Each is written as `.png` for looking at and `.pdf` for the report, and
each is skipped when its input table is absent. Turn the whole step off
with `--report_figures false`; the report then renders without them.

**Candidate landscape.** Cancer cell fraction against agretopicity
(IC50 mutant / wild type), coloured by what prioritisation decided and
sized by expression. The y axis is inverted, so the candidates worth
making are top right: clonal, and binding far better than their wild
type. It is the one figure that shows *why* the ranking came out as it
did rather than only what it decided.

**Flag summary.** How many candidates carry each flag, rejections in red
and cautions in amber. A patient losing everything to `HLA_LOST` is a
different problem from one losing everything to `NO_EXPRESSION`.

**HLA loss of heterozygosity.** One row per locus, coloured by call, with
the alleles written out. When this table is present the report's HLA
section reports the real calls instead of saying allele-specific LOH was
not assessed.

## 1. RNA support across timepoints (`bin/rna_support.py`)

The single-patient version of the cohort figure of RNA support for clonal
neoantigens. It runs for every patient with at least one RNA BAM; with
baseline and week 6 it gives the four groups of that figure.

For each non-synonymous PASS variant, each RNA sample is checked at the
site:

| Status | Meaning |
|---|---|
| `supported` | at least `--rna_support_min_alt_reads` (default 2) RNA reads carry the alt allele, and RNA VAF ≥ `--rna_support_min_vaf` (default 0) |
| `covered_no_alt` | RNA reads at the site, but not enough carry the alt allele |
| `no_coverage` | no RNA reads at the site (gene silent, or too shallow) |

The variant then falls in one group: **Expressed both**, **Week 6 only**,
**Baseline only** or **DNA only**. With one timepoint the groups are
"Expressed (Baseline)" and "DNA only"; with more than two they generalise
to every combination.

The split is shown for four sets of variants:

- **Non-synonymous (all)**: missense, nonsense, frameshift, in-frame
  indel, start/stop lost or protein-altering on the canonical transcript
  (any transcript when VEP flagged none canonical). Splice-region variants
  are not counted.
- **Non-synonymous (clonal)**: CCF ≥ `--clonal_ccf` (default 0.8). CCF is
  computed from DNA VAF, purity, ploidy and the Sequenza segment with the
  same formula `annotate_ccf.py` applies to the pVACseq reports. It is a
  point estimate with no confidence interval.
- **Clonal neoAg MHC I / MHC II**: clonal variants whose best pVACseq
  peptide has IC50 MT ≤ `--neoag_class_i_ic50` (500 nM) or percentile ≤
  `--neoag_class_i_percentile` (2) for class I, and ≤ 1000 nM or ≤ 10 for
  class II. Counted per variant, not per peptide.

The summary table also has all neoantigens (clonal or not).

Outputs, in `<outdir>/<id>/report/`:

- `<id>.rna_support.pdf` / `.png`: counts (left), proportions (middle),
  and per-variant RNA VAF at baseline against week 6 (right; DNA VAF
  against RNA VAF with one timepoint). Colours match the cohort figure.
- `<id>.rna_support.summary.tsv`: counts and percentages per set and group.
- `<id>.rna_support.variants.tsv`: one row per non-synonymous variant,
  with DNA depth/VAF, CCF, per-timepoint RNA depth, alt reads, VAF and
  status, the group, and the best MHC I and II peptide.

Reading it:

- "DNA only" mixes silent genes and sites too shallow to see the allele.
  The variant table keeps depth so the two can be separated.
- A drop from baseline to week 6 can be tumour regression or immune
  infiltration diluting tumour RNA, not only loss of expression. Read it
  next to purity and, where available, a week 6 purity estimate.
- Nonsense and frameshift alleles are often depleted in RNA by
  nonsense-mediated decay.

## 2. Patient report (`bin/patient_report.py`)

A LaTeX report per patient, compiled with `pdflatex`. The layout follows
the NeoDisc full, clinical and gene-panel templates. Sections:

1. **Patient summary**: purity, ploidy, sex, RNA timepoints, genomic
   subtype, BRAF status, PASS and non-synonymous counts, TMB, clonal
   share, UV signature, HLA class I, HLA/B2M copy loss, presentation
   flags and neoantigen counts.
2. **Aim of the analysis.**
3. **Melanoma drivers and genomic subtype**: TCGA SKCM 2015 subtype
   (BRAF hotspot, then RAS hotspot at NRAS/KRAS/HRAS codons 12, 13 and 61,
   then NF1, otherwise triple wild-type), BRAF class I/II/III from the
   protein position, and every coding mutation in the panel's driver
   genes with DNA VAF, CCF, RNA VAF, copy-number call and expression call.
4. **HLA typing and HLA loss**: class I alleles with zygosity, class II
   alleles, and gene-level copy number and expression at HLA-A/B/C and
   B2M from the gene panel. **Allele-specific HLA LOH is not assessed** and
   the report says so.
5. **Antigen presentation machinery**: the gene panel by pathway, with
   adverse and altered genes named. The full panel is in the appendix.
6. **Somatic variants, mutational burden and clonality**: variant classes,
   consequences, TMB (with the 10 mut/Mb line), UV signature (C>T share,
   C>T at dipyrimidines from the reference FASTA, CC>TT doublets), and a
   CCF histogram.
7. **RNA support across timepoints**: the figure and summary table above.
8. **Expression context**: GTEx specificity and expression flags of the
   candidate genes.
9. **Fusion neoantigens**: only meaningful when the optional fusion
   branch ran (`--fusion true`, `docs/FUSION.md` when that is written).
   Three states, each reading differently: the branch did not run for
   this patient; it ran but AGFusion found no coding fusion protein, so
   pVACfuse had nothing to score; or pVACfuse wrote per-class reports,
   in which case the top binders are tabled with junction/spanning read
   support, FFPM, tier and a note for the melanoma-relevant partners
   (ALK, ROS1, NTRK1/2/3, RET, MET, BRAF, MAP3K8, PRKCA). Fusions are
   uncommon in cutaneous melanoma, so an empty table is the expected
   result; they matter most in Spitz/spitzoid tumours, which have few
   point mutations, and a BRAF fusion is also a MAPK-inhibitor
   resistance route. The state is read from the pVACfuse directory
   itself (`--pvacfuse-dir`): absent means not run, a
   `no_fusion_proteins.txt` inside means no coding partner pair.
10. **Top neoantigen candidates**: top `--report_top_n` (default 15) MHC I
   and MHC II binders, ranked by clonal first, then seen in baseline RNA,
   then expression flag OK, then IC50. pVACseq's own Tier is shown beside
   them.
11. **Not assessed by NeoadjVax**: the NeoDisc sections with no
    equivalent (allele-specific HLA LOH, genome-wide CN figures,
    multi-caller overlap, TAAs, inflammation score, virus detection, mass
    spectrometry, QC summaries).
12. **Methods in brief.**

Every input but the VCF and purity is optional. A patient without RNA, or
a run without `--gtex_db`/`--tcga_db`, still gets a report; the sections
that need the missing input say it was not available.

### TMB denominator

In order of preference: `--tmb_callable_mb`; the size of
`--target_intervals` (merged BED or interval_list); otherwise the total
length of `--scatter_chromosomes` from the FASTA index (WGS). TMB counts
all PASS variants; the report also gives non-synonymous per Mb.

## 3. Containers

- `containers.plots` (default `jupyter/scipy-notebook`) for the figure:
  the pvactools image has no matplotlib. Any Python 3 image with
  matplotlib works.
- `containers.texlive` (default `texlive/texlive:latest-medium`) for the
  PDF. Pin it by digest once pulled.
- The report itself is plain Python and runs in the pvactools image.

Both new images are in `bin/prefetch_containers.sh`.

## 4. What was checked

- Unit tests: `tests/test_rna_support.py` and `tests/test_patient_report.py`.
- The stub run checks wiring: P1 (baseline and week 6 RNA) gets the
  RNA-support outputs, P2 (no RNA) does not, both get a report.
- Neither has been run on real patient data.
