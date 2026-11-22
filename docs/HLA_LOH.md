# HLA loss of heterozygosity

A neoantigen is only a target if the tumour still carries the allele that
presents it. Losing one HLA haplotype is a common immune-escape route in
melanoma, and it is invisible in every other column of a candidate table:
the variant is still there, still clonal, still expressed, and the peptide
still binds — to a molecule the tumour no longer makes.

SpecHLA is run **outside** this pipeline, the same way Sequenza is. The
pipeline reads the result.

## Giving the pipeline a result

Add the patient's `merge.hla.copy.txt` to the samplesheet's `hla_loh`
column. That is the output of `cal.hla.copy.pl`, step 3 of the SpecHLA
run (typing on the HLA reads, then copy number using the Sequenza purity
and ploidy):

```csv
patient_id,...,purity,ploidy,hla_loh
NeoTrio_003,...,0.62,3.1,/scratch/jo11/spechla/loh/NeoTrio_003/merge.hla.copy.txt
```

The column is optional per patient. A patient without one gets `NA` in the
`HLA LOH` column and no flag, which is different from being told the
alleles are retained.

The long collated tables the cohort analysis produces
(`HLA_LOH_per_locus_all_cohorts.csv` and friends) are read too — tab,
comma and semicolon separated files all work, and `--patient` matches
against either a `sample_id` or a `patient_id` column.

## What counts as a loss

SpecHLA reports `LOH = Y` for any allelic imbalance, including a 2:1
locus where both alleles are still present and still presenting. That is
not an escaped tumour; it is a copy-number ripple. By default a loss needs
a zero on one side of `copyratio`:

| `--hla_loh_definition` | a loss is |
|---|---|
| `copy_loss` (default) | SpecHLA says `LOH = Y` **and** one allele is at 0 copies |
| `spechla` | whatever SpecHLA called, imbalance included |

A locus is only **evaluable** when the two alleles differ at two-field
resolution (nothing to lose if they are the same), SpecHLA did not mark it
homozygous, and it has at least `--hla_loh_min_het` heterozygous sites.
Below that the call is `not_evaluable`, not `retained`: absence of
evidence is not evidence of retention, and the candidate gets a caution
rather than a pass.

Allelic imbalance is also unreliable at low tumour content, so a sample
below `--hla_loh_min_purity` is reported as `low_purity` and no call from
it is used.

The defaults — `copy_loss`, 10 heterozygous sites, purity 0.20 — are the
thresholds of the NeoadjLOH cohort analysis, so the pipeline and that
analysis call the same patients.

## What it changes

`bin/parse_hla_loh.py` writes `<patient>/hla_loh/<patient>.hla_loh.tsv`,
one row per locus, with each allele also spelled the way pVACseq spells
it (`HLA-A*02:01` for class I, `DRB1*15:01` for class II).

The prioritisation step then reads that table and, for each candidate,
looks up the allele in its `Allele` column:

- **lost** → flagged `HLA_LOST` and rejected. The epitope cannot be
  presented at all.
- **retained** → no flag.
- **not_evaluable** → flagged `HLA_NOT_EVALUABLE`, a caution.

A class II heterodimer is named as a pair (`DQA1*01:02-DQB1*06:02`).
Losing either chain loses the molecule, so the pair is lost if either half
is, and only retained if both halves are.

Nothing is deleted. A rejected candidate keeps its row, its rank and its
reason, so the call can be read back and argued with.

## Related

- [`docs/PRIORITISATION.md`](PRIORITISATION.md): the other flags, and how the ranking works.
- The presentation-machinery panel in [`docs/EXPRESSION.md`](EXPRESSION.md) answers the
  neighbouring question: B2M, TAP and the interferon route can break
  presentation for every allele at once.
