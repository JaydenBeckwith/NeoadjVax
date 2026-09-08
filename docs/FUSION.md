# Gene fusion neoantigen calling

The fusion pipeline (alias gene_fusion / fusion_neoantigen) runs STAR-Fusion
and/or Arriba on tumour RNA, annotates the resulting calls with AGFusion, and
predicts MHC-binding peptides spanning each fusion junction with pVACfuse.
Each caller keeps its own dedicated STAR alignment pass: STAR-Fusion needs
its own CTAT genome resource lib, Arriba needs relaxed multimapping plus
chimeric-alignment flags neither this pipeline's core RNA branch nor the
other caller's STAR pass sets. Results are never merged or forced into
cross-caller consensus: each caller's calls, AGFusion annotation and
pVACfuse prediction are kept in their own per-caller output subdirectory.

## Inputs and execution

`--rna_samplesheet` is the same sheet used by the core RNA/pVACseq branch
(`assets/samplesheet_schema.md`), with two additional optional columns this
branch reads: `starfusion_tsv` and `arriba_tsv`. Set one of these to a
caller's own native output TSV to bypass that caller's STAR pass and calling
step entirely for that patient/timepoint: any caller left without a TSV
still runs its own dedicated STAR alignment and calling from `rna_bam` or
`rna_fastq_r1`/`rna_fastq_r2` on that row. `rna_vcf`-only rows contribute
nothing to this branch: there is no variant file a fusion caller can start
from, so they participate only in `PVACSEQ_CORE`.

~~~csv
patient_id,timepoint,rna_bam,rna_fastq_r1,rna_fastq_r2,rna_vcf,hla_alleles,starfusion_tsv,arriba_tsv
P001,PRE,,/data/P001_PRE_R1.fastq.gz,/data/P001_PRE_R2.fastq.gz,,,,
P002,PRE,,,,,,/data/P002.star-fusion.fusion_predictions.abridged.tsv,/data/P002.arriba_fusions.tsv
~~~

HLA alleles come from one of three sources, in this order of precedence:
the RNA sheet's own `hla_alleles` column (`|`-separated, same convention as
the DNA sheet); `--hla_alleles_manual`, which is restricted to single-patient
runs (a cohort with more than one patient must use the RNA-sheet column
instead, since one global value cannot be per-patient); or, when neither is
set for a patient, reuse of that patient's typed HLA from
`--run_pvacseq_core` (the patient must also appear in `--dna_samplesheet`
for this to work). A patient with none of the three at launch fails
validation before any compute starts.

~~~bash
nextflow run main.nf -profile gadi \
    --pipelines gene_fusion \
    --rna_samplesheet rna_samples.csv \
    --fusion_callers starfusion,arriba \
    --ctat_resource_lib /data/reference/ctat_genome_lib \
    --arriba_blacklist /data/reference/arriba/blacklist_hg38_GRCh38_GENCODE46.tsv.gz \
    --agfusion_container /scratch/jo11/singularity_images/agfusion_1.5.0.sif \
    --agfusion_database /data/reference/agfusion.homo_sapiens.112.db \
    --agfusion_pyensembl_cache /data/reference/pyensembl_cache \
    --pvacfuse_mhcflurry_models /data/reference/mhcflurry_models \
    --hla_alleles_manual "HLA-A*02:01,HLA-A*24:02,HLA-B*07:02,HLA-B*15:01,HLA-C*03:04,HLA-C*07:02" \
    --outdir results
~~~

`--fusion_callers` (default `starfusion,arriba`) selects which caller(s) run;
either alone or both. Arriba's blacklist is not an optional flag here: a
missing `--arriba_blacklist` fails launch validation, since running Arriba
candidate discovery without one is not a safe default (Arriba's own docs
describe this as materially more false positives from read-through
transcription and homology, not a lower-confidence-but-fine mode). The
blacklist, known-fusions and protein-domains files are genome- and
annotation-build-specific and are not bundled with the arriba binary or this
repository: fetch them for your build with Arriba's own
`download_references.sh`, e.g. `download_references.sh hg38+GENCODE46` to
match this pipeline's GENCODE v46 GTF.

To combine fusion calling with the core somatic neoantigen pipeline in one
run:

~~~bash
nextflow run main.nf -profile gadi \
    --pipelines somatic_neoantigen,gene_fusion \
    --dna_samplesheet dna_samples.csv \
    --rna_samplesheet rna_samples.csv \
    --fusion_callers starfusion,arriba \
    --ctat_resource_lib /data/reference/ctat_genome_lib \
    --arriba_blacklist /data/reference/arriba/blacklist_hg38_GRCh38_GENCODE46.tsv.gz \
    --agfusion_container /scratch/jo11/singularity_images/agfusion_1.5.0.sif \
    --agfusion_database /data/reference/agfusion.homo_sapiens.112.db \
    --agfusion_pyensembl_cache /data/reference/pyensembl_cache \
    --pvacfuse_mhcflurry_models /data/reference/mhcflurry_models \
    --outdir results
~~~

Here no `--hla_alleles_manual`/RNA-sheet `hla_alleles` is needed: every RNA
patient that also has a DNA samplesheet row gets its HLA typing reused from
`PVACSEQ_CORE`'s own typing step automatically.

## Software and references

AGFusion (PyPI `agfusion==1.5.0`) needs a local container: build the
supplied recipe, same pattern as the ERVcaller/Telescope custom images.

~~~bash
docker build -t neoadjvax-agfusion:1.5.0 -f containers/agfusion/Dockerfile .
docker save -o agfusion_1.5.0.tar neoadjvax-agfusion:1.5.0
singularity build agfusion_1.5.0.sif docker-archive://agfusion_1.5.0.tar
~~~

pVACfuse runs from the existing pVACtools image, pinned separately from the
core pVACseq branch's own `pvactools` container tag
(`containers.pvacfuse = docker://griffithlab/pvactools:7.1.3`) so upgrading
one branch cannot silently change the other's results.

AGFusion's database is the part most likely to trip up a first run.
AGFusion ships pre-built, S3-hosted databases via `agfusion download`, but
those pre-built databases only cover Ensembl releases up to **92**
([`murphycj/AGFusionDB`](https://github.com/murphycj/AGFusionDB)'s own
listing). This pipeline's reference GTF is GENCODE v46, which corresponds to
**Ensembl release 112** (per GENCODE's own release-history table): too new
for a pre-built download. Build the database directly instead:

~~~bash
agfusion build -d . -s homo_sapiens -r 112 --pfam Pfam-A.clans.tsv
pyensembl install --species homo_sapiens --release 112
~~~

`agfusion build` needs network access to a public Ensembl MySQL server and a
Pfam-to-clan mapping file (`Pfam-A.clans.tsv`, from the Pfam/InterPro FTP
site); `pyensembl install` populates a matching local cache so PyEnsembl's
lookups agree with the AGFusion database's release. Point
`--agfusion_database` at the resulting `agfusion.homo_sapiens.112.db` file
(launch validation enforces this exact `agfusion.homo_sapiens.<release>.db`
naming so a mismatched release can't be pointed at by accident) and
`--agfusion_pyensembl_cache` at the PyEnsembl cache directory. Do this once
per Ensembl release, not per run: cache the built database and PyEnsembl
cache directory alongside the other Gadi reference files.

AGFusion's CLI, as called by `bin/fusion_tools.py annotate` for each caller's
TSV:

~~~bash
agfusion batch -f <caller_tsv> -a <starfusion|arriba> \
    -db agfusion.homo_sapiens.112.db -o agfusion_out --middlestar
~~~

`--middlestar` marks the fusion junction within the predicted protein
sequence; pVACfuse needs it to locate where a fusion-derived peptide differs
from wild type. `--noncanonical` (`--agfusion_noncanonical true`) additionally
explores non-canonical transcript combinations, which can substantially
increase the number of predicted fusion transcripts.

pVACfuse's CLI, as called by `bin/fusion_tools.py predict`:

~~~bash
pvacfuse run agfusion_out <sample> <alleles> <algorithm(s)> pvacfuse_output \
    -e1 8,9,10,11 -e2 15 \
    --binding-threshold 500 \
    --n-threads 4 \
    [--iedb-install-directory /path/to/iedb] \
    [--starfusion-file star-fusion.fusion_predictions.abridged.tsv]
~~~

`--starfusion-file` is only meaningful for STAR-Fusion-derived AGFusion
output: it lets pVACfuse forward STAR-Fusion's own junction/spanning read
counts and FFPM into the aggregated report. Arriba-derived AGFusion output
has no equivalent upstream field, so pVACfuse reports read support and FFPM
as `NA` for Arriba calls in the aggregated report; this is expected, not a
bug in this pipeline's wiring. `--pvacfuse_algorithms` defaults to
`MHCflurry` (offline, needs only `--pvacfuse_mhcflurry_models` pointed at a
model directory containing `manifest.csv`); adding `NetMHCpan` and/or
`NetMHCIIpan` requires a locally installed, licensed IEDB (`mhc_i/` and
`mhc_ii/`) at `--pvacfuse_iedb_install_dir`.

## Outputs and interpretation

Each `results/PATIENT/fusion/TIMEPOINT/CALLER/` directory contains:

- `agfusion_out/`: AGFusion's own per-fusion transcript/CDS/protein/exon
  annotation output, plus `annotation_qc.json` recording the caller, the
  min-read-support/min-FFPM thresholds applied, and whether any fusions
  survived filtering.
- `results/PATIENT/pvacfuse/TIMEPOINT/CALLER/pvacfuse_output/`: pVACfuse's
  standard output tree (filtered/unfiltered aggregate reports, per-length
  binding predictions), plus `run_status.json` recording the sample,
  caller, algorithms and allele set used, and whether any AGFusion protein
  records were available to predict from.

An empty AGFusion annotation (no fusions passed the read-support/FFPM
filters) is a valid negative result, not a failure: `bin/fusion_tools.py`
writes `annotation_qc.json` with `status: no_fusions` and does not invoke
pVACfuse for that patient/timepoint/caller, rather than launching pVACfuse
against nothing. The reverse (a non-empty AGFusion result reporting zero
Predictable fusion protein records) is written to `run_status.json` as
`status: no_predictable_fusion_proteins` rather than a silently-empty
pVACfuse report. A zero-exit AGFusion run that logged an internal error is
treated as a failure, not a successful negative result: see
`tests/test_fusion_tools.py` for the exact conditions this covers.

Filter thresholds (`--fusion_min_read_support`, default 5;
`--fusion_min_ffpm`, default 0.1, STAR-Fusion only, since Arriba reports no
FFPM) are applied identically to both callers' inputs before annotation, not
inside AGFusion or pVACfuse themselves, so the same evidence bar applies
across callers even though their native read-support columns differ.

## Verification

~~~bash
python3 -B -m unittest discover -s tests -v
nextflow -c tests/fusion/smoke.config run tests/fusion/smoke.nf -stub-run
~~~

Run from the repository root. The smoke fixture uses one patient/timepoint
with a precomputed `starfusion_tsv` (so STAR-Fusion calling is bypassed) and
placeholder FASTQs for Arriba (so Arriba's own dedicated STAR pass, AGFusion
annotation and pVACfuse prediction are all structurally exercised): it tests
that every process wires together and every declared output exists, not
biological accuracy. Never run it without `-stub-run`, and never point
`--agfusion_database`/`--agfusion_pyensembl_cache`/`--pvacfuse_mhcflurry_models`
at the fixture's placeholder files outside of `-stub-run`.

Before a cohort run, build the AGFusion database once for your Ensembl
release, build and cache the AGFusion/pVACfuse containers, fetch a real
Arriba blacklist for your genome+annotation build, and check one real
sample through the branch end to end.

## Sources

- [pVACtools documentation: pvacfuse](https://pvactools.readthedocs.io/en/latest/pvacfuse.html)
- [AGFusion (murphycj/AGFusion)](https://github.com/murphycj/AGFusion)
- [AGFusionDB pre-built database listing](https://github.com/murphycj/AGFusionDB)
- [GENCODE release history](https://www.gencodegenes.org/human/releases.html)
- [Arriba (suhrig/arriba)](https://github.com/suhrig/arriba)
- [STAR-Fusion (STAR-Fusion/STAR-Fusion)](https://github.com/STAR-Fusion/STAR-Fusion)
