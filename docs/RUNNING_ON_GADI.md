# Running NeoadjVax on NCI Gadi

This isn't like Jayden's original PBS scripts, where `qsub script.sh`
submitted the whole job directly. With Nextflow, `qsub` mostly disappears
from view: `-profile gadi` sets `process.executor = 'pbspro'`
(`conf/gadi.config`), so **Nextflow submits every single pipeline step as
its own PBS job automatically**: no per-tool PBS script to hand-write,
no manual `-W depend=afterok:...` chaining like the original scripts did.
You only ever launch one thing yourself: the Nextflow *driver* process,
which then does all the `qsub`ing for you as the pipeline runs.

The one real complication is specific to Gadi: **compute nodes have no
external network access.** Nextflow's Singularity integration normally
pulls a container the first time a process needs it: fine on most
clusters, but on Gadi that first pull would happen inside a PBS job on a
compute node with no way to reach Docker Hub / depot.galaxyproject.org,
and it would just fail. Login nodes DO have network. So containers need
pre-caching once, from the login node, before any real run.

## 1. One-time setup: pre-cache containers

From a Gadi **login** node (not inside a job):

```bash
module load singularity
PBS_PROJECT=jo11 bash bin/prefetch_containers.sh
```

This pulls every container the pipeline references (everything in
`nextflow.config`'s `containers {}` block except the four `.sif` files
that already exist on Jayden's `/scratch/jo11`) into
`/scratch/jo11/$USER/nxf_singularity_cache`: the same path
`conf/gadi.config`'s `singularity.cacheDir` points Nextflow at, using the
exact filenames Nextflow itself would generate (confirmed against
Nextflow's `SingularityCache` source), so the real run finds them already
there instead of trying to fetch them again on a compute node.

Re-run this any time a container reference in `nextflow.config` changes
(new tool, version bump).

## 2. Launching the pipeline: two options

**Option A: dedicated driver PBS job (NCI's documented recommendation).**
Edit `bin/submit_nextflow_gadi.pbs` with your actual samplesheet paths,
then:

```bash
qsub bin/submit_nextflow_gadi.pbs
```

This submits a small, cheap job (1 cpu, a few GB, long walltime) that runs
only as the Nextflow orchestrator: it in turn `qsub`s every real pipeline
step as its own job. Log out any time; the driver job keeps running and
keeps submitting/monitoring steps until the pipeline finishes or the
walltime runs out.

**Option B: a login-node `screen`/`tmux` session**, if you'd rather watch
it live or iterate quickly without waiting in the job queue for the driver
itself:

```bash
screen -S neoadjvax
module load nextflow singularity
nextflow run main.nf -profile gadi \
    --dna_samplesheet dna_samples.csv \
    --rna_samplesheet rna_samples.csv \
    --vep_cache /path/to/vep_cache \
    --vep_plugins /path/to/vep_plugins \
    --outdir results/
# Ctrl+A, D to detach: log back in and `screen -r neoadjvax` to reattach
```

Either way, the actual variant calling / alignment / pVACseq work always
runs as real PBS jobs on compute nodes, submitted automatically: these
two options only differ in where the lightweight Nextflow *orchestrator*
itself lives.

## 3. Monitoring

- `qstat -u $USER`: see the driver job (if using Option A) and every
  currently-running pipeline-step job Nextflow has submitted.
- Nextflow's own `work/` directory has one subfolder per task; `.command.log`,
  `.command.err`, `.command.sh` inside each are the same kind of per-step
  logs the original scripts wrote to `pbs_logs/`.
- `results/pipeline_info/execution_report.html` (generated automatically:
  see `nextflow.config`'s `report {}` block) gives a run-level summary
  after completion: per-process walltime, memory, CPU usage.
- `-resume` (already in `bin/submit_nextflow_gadi.pbs`) picks up from
  wherever the last run left off: this replaces the original scripts'
  hand-rolled `checkpoint_run` file-existence checks entirely.

## Troubleshooting

For the optional ERVcaller DNA branch, build and transfer a custom SIF before
submission and pass its local path with `--ervcaller_container`. The normal
container prefetch helper does not provide this image. See
[ERVcaller on DNA](ERV_DNA.md) for the build recipe and required references.

**A step fails immediately with a network/connection error, and its log
mentions pulling or fetching an image**: a container wasn't pre-cached, or
`nextflow.config` was changed since the last `prefetch_containers.sh` run.
Re-run step 1 from a login node, then `-resume`.

**Two runs seem to be fighting over the same work directory**: use a
different `--outdir` (and ideally a separate `-w`/work directory) per
concurrent run: same as running two of the original PBS scripts against
the same `VARIANT_DIR` would have collided.

---

Sources consulted for the Gadi-specific guidance above:
- [nf-core/configs: nci_gadi.md](https://github.com/nf-core/configs/blob/master/docs/nci_gadi.md)
- [NCI Opus Confluence: Nextflow on Gadi](https://opus.nci.org.au/spaces/DAE/pages/138903678/Nextflow)
- [Nextflow Singularity container docs (Seqera)](https://docs.seqera.io/nextflow/container/singularity)
- [nextflow-io/nextflow discussion #5110: Singularity image naming/caching](https://github.com/nextflow-io/nextflow/discussions/5110)
