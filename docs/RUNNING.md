# Choosing what a run does

A full run aligns, calls, types, predicts, prioritises, designs peptides
and writes a report. Most of the time you do not want all of it: you are
rerunning one patient to check a fusion call, or you only want the
presentation-machinery panel across a cohort. `--run` picks the analyses.

```bash
--run all              # everything this run has the reference data for (default)
--run somatic          # SNV/indel neoantigens only
--run fusion           # fusion neoantigens only
--run markers          # expression, gene panel and shared antigens only
--run somatic,fusion   # combine them
```

## What always runs

The core runs whatever you choose, because it is cheap next to the rest
and every read-out depends on it:

- **alignment and duplicate marking** of the tumour and normal DNA
- **HLA typing** (HLA-HD, or the samplesheet's `hla_alleles`)
- **HLA-LOH**, read from the samplesheet's `hla_loh` column
- **Sequenza purity, ploidy and segments**, also read from the samplesheet

None of those three results is computed here, so the last two cost
nothing: they are files you already have.

## The analyses

| `--run` | what it adds |
|---|---|
| `somatic` | Mutect2 and filtering, proximal-variant phasing, VEP, bam-readcount, pVACseq, CCF, prioritisation, long peptides, pVACvector, clonal structure |
| `fusion` | STAR-Fusion (or an existing table), AGFusion, pVACfuse |
| `markers` | featureCounts/TPM, GTEx and TCGA context, the presentation-machinery and driver panel, shared tumour antigens, pVACbind |

`all` means *every analysis this run was given the data for*. The fusion
branch needs a sizeable reference bundle (the AGFusion database, its
pyensembl cache, the CTAT library), so under `all` it is skipped with a
warning when `--agfusion_db` does not exist. Ask for it by name
(`--run fusion`) and a missing reference is an error instead, because then
you asked for something that cannot run.

`--fusion true|false` still works and still wins where it is set, so
existing params files behave as before.

## What you get back

Every selection produces a patient report and a workbook; they contain
what the run actually produced.

- Without `somatic`, the report's variant, TMB, UV-signature and clonality
  sections say **not assessed** rather than reporting zeros. A TMB of 0 is
  a claim about the tumour; "not assessed" is a claim about the run, and
  only the second one is true.
- Without `markers`, candidates carry no expression columns, so
  `NO_EXPRESSION_DATA` appears as a caution rather than expression being
  treated as a pass. See [`PRIORITISATION.md`](PRIORITISATION.md).
- Without `fusion`, the report's fusion section says the branch did not
  run, which is distinct from running and finding no coding fusion.

`--run markers` needs something to compare against: `--gtex_db`,
`--tcga_db`, or `--rna_gtf` with an RNA BAM. Asking for it by name with
none of those is an error rather than an empty run.

## Peptide lengths

| parameter | default |
|---|---|
| `--pvacseq_class_i_lengths` | `8,9,10,11,12,13,14,15` |
| `--pvacseq_class_ii_lengths` | `12,...,25` |

pVACtools runs each predictor only on the lengths it supports and skips
the combinations it cannot, so a wide range costs skipped combinations
rather than a failed run. The long class I lengths come mostly from
MHCflurry and the long class II ones from NetMHCIIpan, so a candidate at
24 residues is supported by fewer predictors than one at 15. The
`Num Passing Peptides` and per-algorithm columns in the pVACseq output
show how many actually scored it.

Long peptides are a separate question: those are designed from the
mutant protein context at `--long_peptide_max_length` (25 by default) and
are not epitope predictions. See [`LONG_PEPTIDES.md`](LONG_PEPTIDES.md).

## Examples

```bash
# the presentation-machinery panel and shared antigens across a cohort
nextflow run main.nf -profile gadi --run markers \
    --samplesheet cohort.csv --gtex_db ... --tcga_db ...

# rerun one patient's fusion branch only, reusing the work directory
nextflow run main.nf -profile gadi --run fusion -resume \
    --samplesheet one_patient.csv --agfusion_db ... --ctat_genome_lib ...

# the usual full run
nextflow run main.nf -profile gadi --run all --samplesheet cohort.csv ...
```

With `-resume`, work already done under a wider selection is reused, so
narrowing `--run` on a rerun skips the steps you dropped rather than
recomputing anything.
