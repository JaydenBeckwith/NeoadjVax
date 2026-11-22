# Shared tumour antigens

A patient with a low mutation burden can simply run out of neoantigens.
Cancer-testis antigens, the melanocyte differentiation antigens and the
broadly overexpressed genes are not patient-specific, so they are not
neoantigens — but they are real targets with characterised epitopes, and
whether one is usable in a given patient is an expression question this
pipeline already answers.

`TAA_REPORT` writes `<patient>/taa/<patient>.taa.tsv`, one row per gene.

## The call

| call | meaning |
|---|---|
| `good_target` | expressed in this tumour, and silent in normal tissue outside testis |
| `target_with_risk` | expressed in this tumour, but normal tissue expresses it too (or there is no GTEx table, so the risk is unknown) |
| `not_expressed` | not expressed in this tumour; nothing to target |
| `no_data` | no RNA for this patient and no cohort table |

Evidence is the patient's own RNA where there is any, and the TCGA cohort
median otherwise — the same fallback, and the same `Evidence` column
convention, as the per-candidate expression annotation.

Testis is treated as immune-privileged, so expression there does not count
against a gene. That is the whole reason the cancer-testis class is
interesting: silent everywhere an immune response could reach, expressed
in the tumour.

The melanocyte differentiation antigens (`MLANA`, `PMEL`, `TYR`, `TYRP1`,
`DCT`, `SOX10`) are in the panel and will come out as `target_with_risk`
in almost every patient, because normal melanocytes express them. They are
listed rather than omitted because they are reliably expressed and
well characterised; the risk column is there to be read, not to be
designed around silently. A response against them can cause vitiligo and
uveitis.

**This call is about expression only.** It says nothing about whether the
patient's alleles can present anything from the gene.

## The panel

[`assets/taa_panel.tsv`](../assets/taa_panel.tsv) ships with the
pipeline: `gene<TAB>class<TAB>note`, with `#` comments. Replace it
wholesale with `--taa_panel my_panel.tsv`. The class only groups the
report; it carries no scoring.

## Epitopes

The other half of this channel is binding predictions on actual peptides,
which depend only on the patient's alleles and not at all on their
variants. Give a peptide FASTA with `--taa_peptides` and pVACbind runs on
it, writing `<patient>/taa/pvacbind/`.

No epitope set ships here. Which published epitopes are worth testing is a
judgement about the literature and about your assay, not a pipeline
default, and shipping a stale list would be worse than shipping none.

## Parameters

| parameter | default | what it does |
|---|---|---|
| `--taa` | `true` | run the expression report |
| `--taa_panel` | `assets/taa_panel.tsv` | the gene list |
| `--taa_peptides` | none | peptide FASTA for pVACbind |
| `--taa_algorithms` | `MHCflurry NetMHCpan BigMHC_IM NetMHCIIpan` | predictors for pVACbind |

Thresholds are shared with the candidate expression annotation
(`--expression_min_tpm`, `--gtex_expressed_tpm`, `--expression_quantile`),
so a gene is not called expressed here and unexpressed there.
