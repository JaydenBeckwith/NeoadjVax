# NeoadjVax

Neoadjuvant melanoma pipelines for determining neoantigen candidates for
vaccine development.

Orchestrated with Nextflow DSL2. See [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md)
for the full pipeline design, a diagram of how the branches fit together,
and a list of open decisions that need input before the newer branches
(splicing/fusion/ERV neoantigens) are ready to run for real. See
[`docs/RUNNING_ON_GADI.md`](docs/RUNNING_ON_GADI.md) for how to actually
launch this on Gadi — it's not a `qsub script.sh` like the original
scripts; Nextflow submits every step as its own PBS job for you, but
compute nodes having no external network means a one-time container
pre-caching step is required first.

## Quick start

```bash
# once, from a Gadi LOGIN node (see docs/RUNNING_ON_GADI.md)
bash bin/prefetch_containers.sh

# then either qsub bin/submit_nextflow_gadi.pbs, or run this directly
# in a screen/tmux session on the login node:
nextflow run main.nf -profile gadi \
    --dna_samplesheet dna_samples.csv \
    --rna_samplesheet rna_samples.csv \
    --vep_cache /path/to/vep_cache \
    --vep_plugins /path/to/vep_plugins \
    --outdir results/
```

See `assets/samplesheet_schema.md` for the samplesheet formats, and
`nextflow.config` for every tunable default (reference paths, container
images, HLA/pVACseq settings, per-branch on/off switches).

## Branches

| `--run_*` flag | Status |
|---|---|
| `run_dna_variant_calling` | Ported from the original DNA script, default on |
| `run_rna_variant_calling` | Ported from the original RNA PBS script, default on |
| `run_pvacseq_core` | New RNA-support matching + pVACseq, default on |
| `run_splicing_neoantigens` | New, default off — peptide step is a stub |
| `run_fusion_neoantigens` | New, default off — annotation step is a stub |
| `run_erv_neoantigens` | New, default off — peptide step is a stub |
| `run_purity_ploidy` | Fully ported from `Sequenza_tools`, default off — WGS/200bp bins only (WES errors at launch, see `docs/ARCHITECTURE.md`) |
| `run_hla_loh` | Ported from [`NeoadjLOH`](https://github.com/JaydenBeckwith/NeoadjLOH), default off — independently runnable, see `assets/samplesheet_schema.md` |

`run_purity_ploidy` and `run_hla_loh` don't need the DNA/RNA branches above
turned on — `run_hla_loh` can run entirely standalone via
`--hla_loh_samplesheet`, or chained off `--run_dna_variant_calling` with
purity/ploidy from `--purity_ploidy_csv`.
