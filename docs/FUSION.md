# Gene-fusion neoantigens

An optional branch (`workflows/fusion.nf`, `--fusion true`) that turns a
patient's tumour RNA into fusion neoantigen candidates:

```
baseline tumour RNA FASTQ ─> STAR-Fusion ─┐
                                          ├─> AGFusion ─> pVACfuse
samplesheet starfusion_tsv ───────────────┘
```

It mirrors the Gadi PBS scripts `run_agfusion.pbs` and `pvacfuse.pbs`. It
is off by default: most patients have no fusion worth vaccinating
against, and the branch needs reference data (a CTAT library, an AGFusion
database and a pyensembl cache) that the rest of the pipeline does not.

## Why bother, in melanoma

Cutaneous melanoma is driven by a high UV point-mutation load rather than
by structural rearrangement, so a run that finds no fusion candidate is
the normal result, not a failure. Two cases make the branch worth having:

- **Spitz and spitzoid tumours** carry recurrent kinase fusions (ALK,
  ROS1, NTRK1/2/3, RET, MET, BRAF, MAP3K8, PRKCA) and few point
  mutations, so conventional candidates are thin exactly where a fusion
  is likely.
- **A BRAF fusion** is a recognised acquired route to MAPK-inhibitor
  resistance, so it is worth knowing about whether or not it makes a good
  epitope.

A fusion junction peptide also has no wild-type counterpart anywhere in
the normal proteome, so there is no central tolerance against it — which
is the strongest theoretical argument for this class of target.

## Inputs

Per samplesheet row, one of:

| Column | Meaning |
|---|---|
| `rna_baseline_r1`, `rna_baseline_r2` | tumour RNA FASTQ pair; STAR-Fusion runs on it |
| `starfusion_tsv` | the **full** predictions table of an existing STAR-Fusion run |

Giving both for one patient is an error; giving neither skips that
patient. Note these are separate from `rna_baseline_bam`, which feeds
bam-readcount and featureCounts: STAR-Fusion needs reads, not an
alignment to this pipeline's reference.

Use the full `star-fusion.fusion_predictions.tsv`, not the abridged one:
pVACfuse reads junction and spanning read counts and FFPM out of it.

## Reference data

| Param | What it is |
|---|---|
| `ctat_genome_lib` | the `ctat_genome_lib_build_dir` of a CTAT plug-n-play bundle. Only needed when some patient brings FASTQ. Its STAR index must be readable by the image's STAR (2.7.11b), so use a bundle released for STAR-Fusion 1.13 or later. Checked at launch for `ref_genome.fa.star.idx` |
| `agfusion_db` | `agfusion.homo_sapiens.<release>.db`. The filename must carry the Ensembl release, because the pyensembl cache has to match it |
| `agfusion_pyensembl_cache` | `PYENSEMBL_CACHE_DIR`. Must already hold an indexed cache for the same release at `<cache>/pyensembl/GRCh38/ensembl<release>` |
| `agfusion_pythonpath` | directory holding the `sitecustomize.py` numpy/datacache shim the AGFusion 1.4.3 image needs. Set it to `''` if yours does not |

Build the pyensembl cache **on a login node**, since compute nodes have
no network and AGFusion will otherwise try to download it mid-job:

```bash
PYENSEMBL_CACHE_DIR=/scratch/jo11/agfusion/pyensembl \
  pyensembl install --release 111 --species homo_sapiens
```

`main.nf` checks all of this before submitting a single task, so a
missing cache fails in seconds rather than after STAR-Fusion has run.

**Keep the releases aligned.** Ensembl 111 is GENCODE 45 and matches the
VEP 111 cache the SNV/indel half uses. The pipeline warns rather than
errors when `agfusion_db`'s release differs from `vep_cache_version`, but
the two halves of a report then rest on different transcript models.

## Steps

1. **STAR-Fusion** (`modules/local/fusion/star_fusion.nf`) aligns with
   its own STAR against the CTAT library, not this pipeline's
   `genome_fasta`. Both prediction tables are published; the full one
   goes downstream.
2. **AGFusion** (`modules/local/fusion/agfusion.nf`) runs
   `agfusion batch -a starfusion --middlestar` plus `--noncanonical`
   unless `--agfusion_noncanonical false`. `--middlestar` marks the
   junction in each fusion protein, which is how pVACfuse knows where the
   novel sequence starts. A predictions table with no calls is a valid
   negative result and is written out as `agfusion/no_fusions.txt`.
3. **pVACfuse** (`modules/local/pvactools/pvacfuse.nf`) uses the same
   pVACtools image, IEDB install and HLA alleles as pVACseq, with
   `--starfusion-file` so each candidate carries its read support and
   FFPM. pVACfuse's own tiering flags read support below 5 and expression
   below 0.1 FFPM. Its defaults here are `pvacfuse.pbs`'s, not the
   `pvacseq_*` ones: longer peptide lengths both classes, a different
   algorithm set, and `--problematic-amino-acids C`.

   If AGFusion produced no `*_protein.fa` — no fusion had a coding
   partner pair — pVACfuse is skipped and
   `pvacfuse/no_fusion_proteins.txt` is written instead. That is a
   normal outcome: most fusion calls are intergenic, out of frame, or
   non-coding on at least one side.

## Outputs

```
<outdir>/<id>/fusion/star_fusion/star-fusion.fusion_predictions.tsv
<outdir>/<id>/fusion/star_fusion/star-fusion.fusion_predictions.abridged.tsv
<outdir>/<id>/fusion/agfusion/
<outdir>/<id>/pvacfuse/MHC_Class_I/<id>.MHC_I.all_epitopes.aggregated.tsv
<outdir>/<id>/pvacfuse/MHC_Class_II/<id>.MHC_II.all_epitopes.aggregated.tsv
<outdir>/<id>/pvacfuse/combined/
```

The patient report reads this directory and gives it a section of its
own; see [`REPORT.md`](REPORT.md).

## Reading the candidates

- **Read support is the first filter, not the IC50.** A fusion candidate
  rests on the junction having been called correctly and in frame, which
  is a harder call than a point substitution. Confirm a low-read fusion
  at the breakpoint (IGV, or RT-PCR) before treating it as a target.
- A fusion found only in RNA says nothing about how many tumour cells
  carry it. There is no CCF equivalent here: clonality would need the
  genomic breakpoint from DNA, which this pipeline does not call.
- Read-through transcripts between neighbouring genes are a common
  artefact class and are biologically uninteresting; STAR-Fusion
  annotates them, so check the `annots` column of the predictions table.
- Expression (FFPM) comes from the fusion's own reads, so it is
  independent of the gene-level TPM in the expression branch.

## What this is not

- No DNA-level structural variant calling: fusions come from RNA only.
- No fusion-specific clonality, and no week 6 comparison — the branch
  runs on the baseline RNA.
- No splicing, intron-retention or ERV/transposable-element neoantigens.
- The stub test covers wiring only: the fusion processes never run a real
  STAR-Fusion, AGFusion or pVACfuse in CI.

## Tests

`tests/dna_somatic` gives P1 an existing `starfusion_tsv` and P2 an RNA
FASTQ pair, so one patient exercises STAR-Fusion and both reach AGFusion
and pVACfuse. `check_results.py` asserts those task counts, the published
outputs, and that a patient with its own table did **not** re-run
STAR-Fusion.
