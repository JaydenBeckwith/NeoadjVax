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

**Option B: a tmux session**, if you'd rather watch it live or iterate
quickly without waiting in the queue for the driver itself.
`bin/run_somatic_gadi.sh` runs a preflight check, then starts the driver
in a detached tmux session. The preflight checks the Nextflow version,
the samplesheets, the VEP cache and plugins, the reference files and the
cached containers.

Run it inside an NCI **persistent session**, not on a plain login node:
NCI kills long-running processes on the login nodes, and persistent
sessions exist for exactly this kind of workflow-manager process. Check
the current syntax with `persistent-sessions --help`.

```bash
# from a Gadi login node
persistent-sessions start -p jo11 neoadjvax
ssh neoadjvax.${USER}.jo11.ps.gadi.nci.org.au

cd /scratch/jo11/$USER/NeoadjVax          # wherever the repo is cloned
cp assets/params_somatic_gadi.example.yaml /scratch/jo11/$USER/somatic.yaml
# edit somatic.yaml: samplesheets, outdir, vep_cache, vep_plugins

bash bin/run_somatic_gadi.sh --check /scratch/jo11/$USER/somatic.yaml   # preflight only
bash bin/run_somatic_gadi.sh /scratch/jo11/$USER/somatic.yaml           # launch in tmux

tmux attach -t neoadjvax_somatic          # Ctrl+B then D to detach
```

The tmux session keeps going after you log out. Its output is also teed
to `logs/somatic_<timestamp>.log`. Rerunning the script after a failure
resumes from the last completed step: it always passes `-resume` with the
same work directory (`/scratch/jo11/$USER/neoadjvax_work`).

The same params file works for Option A:
`qsub -v PARAMS_FILE=/scratch/jo11/$USER/somatic.yaml bin/submit_nextflow_gadi.pbs`.

### Mutational track checklist

`--pipelines somatic_neoantigen` covers DNA calling, RNA calling and
pVACseq. Before its first real run:

1. `bash bin/prefetch_containers.sh somatic` on a login node. This pulls
   only BWA, Trim Galore, FastQC, VEP 111, pVACtools 7.1.3 and xHLA.
2. Use a VEP cache that matches the container: `homo_sapiens/111_GRCh38`.
   The plugin directory needs pVACtools' `Wildtype.pm` and
   `Frameshift.pm`.
3. Use Nextflow ≥ 24.04. `conf/gadi.config` uses `resourceLimits`, which
   caps retried steps to one normal-queue node (48 cpus, 190 GB, 48 h)
   instead of letting qsub reject them.
4. pVACseq reads the tumour and normal VCF columns `tumor` and `normal`.
   These come from Mutect2's read-group sample names. If you supply your
   own `dna_vcf`, set `pvacseq_tumor_sample_name` and
   `pvacseq_normal_sample_name` to match it.

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
