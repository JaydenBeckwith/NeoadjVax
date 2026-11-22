# Prioritisation: one ranked list, with reasons

By the time pVACseq, `ANNOTATE_CCF`, `ANNOTATE_EXPRESSION` and
`PARSE_HLA_LOH` have run, everything needed to judge a candidate is in the
table. What was missing was the integration: deciding what to make still
meant reading eight columns per row with the thresholds held in your head,
and the decision left no record.

`bin/prioritize_candidates.py` applies every criterion once, writes down
what failed, and orders what is left. It runs last on the candidate
tables and produces
`<patient>/prioritized/<patient>.<class>.…prioritized.tsv`.

Nothing is dropped. Every candidate keeps its row, so a rejection can be
read back and argued with.

## Two derived columns

| column | what it is |
|---|---|
| `Agretopicity` | `IC50 MT / IC50 WT`. Below 1 the mutant binds better than the wild type, so the mutation is what creates the epitope and the T-cell repertoire against it is less likely to have been tolerised. Near or above 1 the wild-type peptide is presented by normal cells too. |
| `DAI` | `IC50 WT − IC50 MT`: the same comparison as a difference, which is steadier than the ratio when both peptides are weak binders. |

## Running the two calculations on their own

Agretopicity and cancer cell fraction each have a script that runs
outside the pipeline, against pVACtools output as it comes. Useful for
checking a finished run, for trying a different threshold without
rerunning anything, and because each is then testable on its own.

Point either at a single aggregated report or at a whole `pvacseq/` or
`pvacfuse/` directory; a directory gives one output file per MHC class.
Both take RNA timepoints with `--timepoint LABEL=VCF`, repeated, so the
baseline and week 6 columns land in the same table.

```bash
# agretopicity for every class, from the pVACseq output directory
python bin/agretopicity.py results/P1/pvacseq --outdir P1_agretopicity/

# one report, with both RNA timepoints joined on
python bin/agretopicity.py results/P1/pvacseq/combined/tumor.Combined.all_epitopes.aggregated.tsv     -o P1.agretopicity.tsv     --timepoint "Baseline=results/P1/pvacseq_input/P1.somatic.vep.readcounts.vcf.gz"     --timepoint "Week 6=results/P1/pvacseq_input/P1.rna_week6.readcounts.vcf"

# cancer cell fraction, same shape
python bin/annotate_ccf.py results/P1/pvacseq --purity 0.62 --ploidy 3.1     --sex female --segments P1_segments.txt --outdir P1_ccf/
```

They chain, so running one over the other's output gives a table with
both, which `patient_report.py` takes directly:

```bash
python bin/agretopicity.py results/P1/pvacseq --outdir P1_agr/
python bin/annotate_ccf.py P1_agr/tumor.Combined.all_epitopes.aggregated.agretopicity.tsv     --purity 0.62 -o P1.Combined.annotated.tsv
python bin/patient_report.py --patient P1 --purity 0.62 --callable-mb 35     --class-i P1.Combined.annotated.tsv -o P1.report.tex
```

### What pVACfuse cannot give you

Both scripts read pVACfuse output, and both decline a calculation there,
for reasons that are about the biology rather than a missing column:

- **Agretopicity does not apply.** A fusion junction peptide has no
  wild-type counterpart anywhere in the germline, so there is nothing to
  divide by. The class comes out as `not_applicable`.
- **CCF does not apply.** Fusions are called from RNA and pVACfuse
  reports no DNA VAF, so there is no allele fraction to correct for
  purity and copy number. `CN Source` comes out as `not_applicable`.

Both say so on stderr and still write the table, so the rest of the
columns pass through. Neither fills the field with `NA`, which would read
like a calculation that failed rather than one that does not apply. For
fusions, the columns that carry the evidence are pVACfuse's own
`Read Support` and `Expr`.

## Flags

A **rejection** means the candidate should not go into a vaccine. A
**caution** means it is worse than its neighbours but still a candidate.

### Rejections

| flag | meaning | threshold |
|---|---|---|
| `HLA_LOST` | the restricting allele is deleted in the tumour | see [`docs/HLA_LOH.md`](HLA_LOH.md) |
| `LOW_BINDING` | weak on percentile **and** IC50 | `--priority_max_percentile` 2.0, `--priority_max_ic50` 500 |
| `NO_EXPRESSION` | the gene is not expressed in this tumour | from `Expression Flag` |
| `NO_RNA_SUPPORT` | RNA covers the site but no read carries the variant | `--priority_min_rna_depth` 5 |

Binding needs to fail on both measures because either alone is noisy and
the predictors disagree often enough that one bad number is not a verdict.

### Cautions

| flag | meaning | threshold |
|---|---|---|
| `SELF_LIKE` | the mutant binds no better than the wild type | `--priority_max_agretopicity` 0.5 |
| `SUBCLONAL` | present in a minority of tumour cells | `--priority_min_ccf` 0.5 |
| `LOW_EXPRESSION` | the gene is weakly expressed | from `Expression Flag` |
| `NORMAL_TISSUE` | broadly expressed in healthy tissue, so on-target off-tumour toxicity is a risk | from `GTEx Specificity` |
| `LOW_DNA_VAF`, `LOW_DEPTH` | thin DNA evidence, so every other number is noisy | `--priority_min_dna_vaf` 0.05, `--priority_min_dna_depth` 20 |
| `NO_RNA_DATA`, `NO_EXPRESSION_DATA` | the question was not asked | — |
| `HLA_NOT_EVALUABLE` | LOH could not be called for this allele | — |

Missing data is deliberately never a rejection. "No RNA coverage here" and
"RNA coverage with no variant reads" are different statements, and only
the second is evidence against the candidate.

## Ranking

The order is a sort, not a score, so two candidates ranked next to each
other differ in a way that can be named:

1. rejected candidates last
2. then fewest cautions
3. then highest CCF (clonal before subclonal)
4. then lowest agretopicity (most foreign first)
5. then lowest `%ile MT` (strongest binder)
6. then highest expression
7. then the variant ID, so ties are stable between runs

`Selection` collapses that into one word — `Top ranking`, `Caution`,
`Rejected` — for filtering a spreadsheet, and `Reject Reason` carries the
first rejection in plain words.

## Clonal structure

`bin/clonal_structure.py` groups the variants by CCF and reports whether
each cluster is still transcribed at week 6
(`<patient>/clonal_structure/`). Clusters are cut wherever sorted CCFs are
more than `--clonal_gap` apart, and clusters below `--clonal_min_cluster`
variants are merged into their nearest neighbour; a cluster is labelled
clonal when its mean CCF reaches `--clonal_cluster_ccf` (distinct from the
report's `--clonal_ccf`, which calls a single variant clonal). That is
plainer than a
mixture model on purpose: the input is one point estimate per variant with
no confidence interval, and a fitted model would read as more certain than
that supports. The figure it draws can be checked against the cut by eye.

Persistence comes from the RNA timepoints, because these patients have RNA
at baseline and week 6 but tumour DNA at baseline only. It is weaker
evidence than DNA: losing RNA support can mean the clone is gone, that the
gene stopped being transcribed, or that the biopsy caught less tumour. The
status says which of those the data can distinguish — `not_detected` means
covered and absent, `no_coverage` means the question was never asked.

A tumour evolution tree needs several regions or DNA timepoints. With one
DNA timepoint there is nothing to build one from, so this stops at
clusters.

## The figures

`bin/report_figures.py` draws the candidate landscape (cancer cell
fraction against agretopicity, coloured by `Selection`) and the flag
summary from this same table, and the patient report embeds them. They
can be regenerated on their own from a finished run; see
[`REPORT.md`](REPORT.md).

## The workbook

`bin/write_workbook.py` collects whatever the run produced into one
`<patient>.xlsx`: candidates per MHC class, the designed long peptides,
HLA-LOH, clonal structure, shared antigens, the gene panel, and a legend
sheet that spells out every flag above. The header is frozen and
filterable, rejected rows are red and cautioned rows amber.

It is written with the Python standard library alone, so the pVACtools
container needs no extra packages. Turn it off with `--workbook false`.
