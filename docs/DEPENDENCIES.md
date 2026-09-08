# External reference data and dependencies

This is a single index of every external file NeoadjVax needs that isn't
part of this repository: where to get or build it, and which `--param`
passes it into the pipeline. It doesn't repeat each branch's full setup
walkthrough: those live in [`docs/FUSION.md`](FUSION.md),
[`docs/ERV_DNA.md`](ERV_DNA.md), [`docs/SPLICING.md`](SPLICING.md) and
[`docs/ARCHITECTURE.md`](ARCHITECTURE.md), linked from each section below.
For container images specifically (as opposed to reference data), see
[`docs/RUNNING_ON_GADI.md`](RUNNING_ON_GADI.md)'s prefetch step: most are
pulled automatically, a handful need a one-off custom build, noted in
context below and summarised again in the final section.

`nextflow.config`'s `params {}` block is the single source of truth for
every default value mentioned here; treat this doc as a map onto that
file, not a replacement for reading it.

## 1. Core reference genome and annotation

Used by every branch (`--genome_fasta`/`--gtf` directly; `--star_index_dir`
by any RNA-aligning branch).

| Param | What it is | Where to get it | Default (Jayden's Gadi scratch) |
|---|---|---|---|
| `--genome_fasta` | GRCh38 primary-assembly FASTA, uncompressed | [GENCODE Human Release 46](https://www.gencodegenes.org/human/release_46.html): `GRCh38.primary_assembly.genome.fa.gz` | `/scratch/jo11/GRCh38.fa/GRCh38.primary_assembly.genome.fa` |
| `--gtf` | Matching GENCODE annotation | Same release page: `gencode.v46.chr_patch_hapl_scaff.annotation.gtf.gz` | `/scratch/jo11/GRCh38.fa/gencode.v46.chr_patch_hapl_scaff.annotation.gtf` |
| `--star_index_dir` | Pre-built STAR genome index | Build once with `STAR --runMode genomeGenerate` against the FASTA+GTF above (`--sjdbOverhang` = read length minus one), or leave the branch to build it inline from `--genome_fasta`/`--gtf` the first time (`modules/local/rna/star_index.nf`: links from here instead of rebuilding if this path already has a populated `SA` file) | `/scratch/jo11/STAR_index_GRCh38_gencode46` |

**Every reference-consuming param in this pipeline expects this exact
build** (GRCh38, GENCODE v46 gene models): mixing in a different genome
build or annotation version anywhere downstream (known-sites VCFs, VEP
cache, splicing/fusion/ERV reference files) silently produces wrong or
missing calls rather than an error. Confirm contig naming matches too
(`chr1` vs `1`): GENCODE's primary-assembly FASTA already uses `chr`-prefixed
names, so anything else pulled from elsewhere may need renaming first.

## 2. GATK known-sites bundle (RNA branch: BQSR + HaplotypeCaller)

| Param | What it is | Where to get it | Default |
|---|---|---|---|
| `--known_mills` | Mills & 1000G gold-standard indels, hg38 | `gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` | `/scratch/jo11/NeoTrio_RNA/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` |
| `--known_1000g` | 1000 Genomes phase 1 high-confidence SNPs, hg38 | Same bucket, `hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz` (confirm the exact filename against the bucket listing at [GATK's Resource Bundle article](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle) before relying on it: this is the one file of the three not present in GATK's own WGS-germline workflow inputs, so it's worth a direct check) | `/scratch/jo11/NeoTrio_RNA/1000G_phase1.snps.high_confidence.hg38.vcf.gz` |
| `--known_dbsnp` | dbSNP 138, hg38 | Same bucket, `hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz` | `/scratch/jo11/NeoTrio_RNA/dbsnp_138.hg38.vcf.gz` (a renamed copy of the file above) |

The Broad's public GCS bucket for this bundle moved at some point from
`gs://genomics-public-data/...` to `gs://gcp-public-data--broad-references/...`;
users have reported access issues on the older path as recently as this
year. If a `gsutil`/`gcloud storage` pull 403s, check GATK's own
[Resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)
page for the current location rather than assuming the path above is
still right; the GATK FTP mirror that used to also host this bundle was
retired in 2020. Each VCF needs a `.tbi` next to it (`TABIX_INDEX`,
`modules/local/rna/index_known_sites.nf`, builds one automatically via
`gatk IndexFeatureFile` if it's missing, so a bare `.vcf.gz` without an
index still works, just costs one extra step per run).

## 3. VEP cache and plugins (core pVACseq branch)

| Param | What it is | Where to get it | Default |
|---|---|---|---|
| `--vep_cache` | Local offline VEP cache (GRCh38, matching Ensembl release) | `vep_install` (bundled with the `ensembl-vep` container/conda package) with `--AUTO cf --SPECIES homo_sapiens --ASSEMBLY GRCh38 --CACHEDIR <dir>`, or see the [VEP cache download instructions](http://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache) | `null`, must be supplied |
| `--vep_plugins` | A directory containing `Wildtype.pm` and `Frameshift.pm` | Either `pvacseq install_vep_plugin <dir>` (pVACtools ships them directly, see [pVACtools' VEP annotation guide](https://pvactools.readthedocs.io/en/latest/pvacseq/input_file_prep/vep.html)), or clone [`Ensembl/VEP_plugins`](https://github.com/Ensembl/VEP_plugins) and use its copies | `null`, must be supplied |

Both plugins are required, not optional: `VEP_ANNOTATE`
(`modules/local/dna/vep_annotate.nf`) always passes
`--plugin Wildtype,<fasta> --plugin Frameshift,<fasta>`, and pVACseq needs
their annotation fields for correct wild-type/mutant peptide pairing and
frameshift neoantigen calling. `pvacseq install_vep_plugin` is the easier
path since it also guarantees the exact plugin version pVACtools was
tested against, rather than whatever's current on `Ensembl/VEP_plugins`'
main branch.

## 4. HLA typing (xHLA)

No separate reference download: xHLA types class I + class II (DRB1) from
one indexed normal-DNA BAM inside its own container
(`params.containers.xhla`, prefetched like any other container, see
[`docs/RUNNING_ON_GADI.md`](RUNNING_ON_GADI.md)). It needs the BAM aligned
against a primary-assembly reference **without ALT contigs** (GENCODE's
`GRCh38.primary_assembly.genome.fa` from section 1 already qualifies): a
GATK-bundle-style reference with ALT contigs can trap reads that should
land in the HLA region, breaking typing. See
`modules/local/hla_typing/xhla.nf`'s header comment for the one detail
here (the in-container entrypoint) that's sourced from a third-party
wrapper rather than xHLA's own docs, worth a one-off sanity check before
a real run.

Skip it entirely per-patient by supplying `hla_alleles` directly in the
DNA samplesheet (`|`-separated, see
[`assets/samplesheet_schema.md`](../assets/samplesheet_schema.md)).

## 5. Splicing branch (`--pipelines splicing`)

See [`docs/SPLICING.md`](SPLICING.md) for the full input contract; the two
external reference dependencies are:

| Param | What it is | Where to get it |
|---|---|---|
| `--spliceai_annotation` | SpliceAI gene annotation table | `grch37`/`grch38` (SpliceAI's own bundled tables, GENCODE v24) for fresh SpliceAI annotation, or a custom table matching your GTF version if you need v46-consistent gene models |
| `--splicing_normal_junctions` | Curated normal/GTEx splice-junction panel | Not bundled: build your own from a reference-matched normal RNA-seq cohort (e.g. GTEx), one-column TSV, header `junc_id` (see SPLICING.md for the exact coordinate convention) |

`--splicing_tools_container` (splice2neo + supporting Python) is pulled
automatically like any other container; no manual build needed.

## 6. Gene fusion branch (`--pipelines gene_fusion`)

Full walkthrough: [`docs/FUSION.md`](FUSION.md). Summary of every external
dependency and its param:

| Param | What it is | Where to get it |
|---|---|---|
| `--ctat_resource_lib` | STAR-Fusion's CTAT genome resource lib | STAR-Fusion's own [CTAT genome lib downloads](https://github.com/STAR-Fusion/STAR-Fusion/wiki) (pick the GRCh38 build matching your GENCODE release) |
| `--arriba_blacklist` (mandatory once fusion calling is on), `--arriba_known_fusions`, `--arriba_protein_domains` | Arriba's genome+annotation-build-specific reference files | Arriba's own bundled `download_references.sh hg38+GENCODE46` (matches this pipeline's GENCODE v46 GTF) |
| `--agfusion_database` | AGFusion's transcript/CDS/protein annotation DB | Build locally: `agfusion build -d . -s homo_sapiens -r 112 --pfam Pfam-A.clans.tsv` (pre-built downloads only cover Ensembl release <=92; GENCODE v46 is Ensembl 112). `Pfam-A.clans.tsv` comes from the [Pfam/InterPro FTP site](https://www.ebi.ac.uk/interpro/download/Pfam/) |
| `--agfusion_pyensembl_cache` | Matching PyEnsembl cache | `pyensembl install --species homo_sapiens --release 112`, same Ensembl release as the AGFusion database above |
| `--agfusion_container` | AGFusion, containerized (no public image) | Build `containers/agfusion/Dockerfile` in this repo, convert to a `.sif` |
| `--pvacfuse_mhcflurry_models` | MHCflurry class-I affinity models | `mhcflurry-downloads fetch models_class1_presentation` (ships a `manifest.csv` in the resulting directory, which is what this param should point at) |
| `--pvacfuse_iedb_install_dir` | Local IEDB install (only if using NetMHCpan/NetMHCIIpan instead of the default MHCflurry) | IEDB's own licensed offline install packages (`mhc_i`/`mhc_ii`), not redistributable, fetch directly from IEDB |

## 7. ERV branches

**RNA quantification** (`--pipelines erv`, `--run_erv_neoantigens`):

| Param | What it is | Where to get it |
|---|---|---|
| `--erv_annotation_gtf` | Telescope-format ERV/TE annotation GTF | One of Telescope's own official builds, [`mlbendall/telescope_annotation_db`](https://github.com/mlbendall/telescope_annotation_db/tree/master/builds); `retro.hg38.v1` is the recommended default. **Not** a generic RepeatMasker dump |

**DNA non-reference insertions** (`--pipelines erv_dna`): see
[`docs/ERV_DNA.md`](ERV_DNA.md) for the full picture.

| Param | What it is | Where to get it |
|---|---|---|
| `--ervcaller_te_fasta` | TE/ERV reference library FASTA | Pick a library that matches the question you're asking (an ERV-focused library differs from a mixed Alu/LINE1/SVA/HERV library): not bundled, source per-study |
| `--ervcaller_container` | ERVcaller, containerized (no public image) | Build `containers/ervcaller/Dockerfile` in this repo (pins upstream commit `802dbf4`), convert to a `.sif` |
| `--ervcaller_reference_dir` | Optional: a pre-indexed genome+TE reference directory, to skip rebuilding BWA/FASTA indices every run | Copy from a completed `ERVCALLER_INDEX` task's work directory once you have one; see ERV_DNA.md for the exact expected file layout |
| `--ervcaller_script` | Path to `ERVcaller_v1.4.pl` inside the container | Already baked into the custom image above; only change this if you rebuild with a different upstream commit |

## 8. Tumor purity/ploidy (Sequenza, `--run_purity_ploidy`)

Fully ported from [`Sequenza_tools`](https://github.com/JaydenBeckwith/Sequenza_tools);
see `docs/ARCHITECTURE.md`'s purity/ploidy section for the branch design.

| Param | What it is | Where to get it | Default |
|---|---|---|---|
| `--sequenza_gc_file` | Genome-wide GC-content reference (`sequenza-utils gc_wiggle`) | Generate once from `--genome_fasta` with `sequenza-utils gc_wiggle -w 50 --fasta <fasta> -o hg38.gc50Base.txt.gz` (matches this branch's 50-base window) | `/scratch/jo11/neoadjuvant/hg38.gc50Base.txt.gz` |
| `--sequenza_gender_csv` | Per-patient gender/melpin lookup, exact column-position format from the original `sequenza_step2.sh`/`step4.sh` | Your own cohort metadata file, not a downloadable reference: see `assets/samplesheet_schema.md` for the exact column-position convention (melpin at column 2, gender spelled `male`/`female` at column 7) | `null`, must be supplied |
| `--sequenza_conda_sh`, `--sequenza_conda_bin_fallback`, `--sequenza_conda_env` | Path to the `r_sequenza` conda environment `SEQUENZA_MERGE_BINS`/`SEQUENZA_FIT` run in (outside the Singularity container: see those modules' comments for why) | Wherever the `r_sequenza` conda env is actually built on your system: the checked-in defaults point at a personal Gadi home directory (`/home/562/jb1592/...`) and need editing before a real run | see `nextflow.config` |

The `sequenza` container itself is an existing local `.sif` on Jayden's
Gadi scratch (`params.containers.sequenza`), not something
`prefetch_containers.sh` pulls.

## 9. HLA loss-of-heterozygosity (SpecHLA, `--run_hla_loh`)

Ported from [`NeoadjLOH`](https://github.com/JaydenBeckwith/NeoadjLOH). No
separate reference-file download beyond the genome FASTA (section 1): the
three in-container command paths below and the SpecHLA container itself
are carried over from `NeoadjLOH`'s own `config.sh`, which by its own
admission was never fully inspected inside the `.sif` (see
`docs/ARCHITECTURE.md`'s "unverified detail" note).

| Param | What it is | Default |
|---|---|---|
| `--spechla_extract_cmd`, `--spechla_typing_cmd`, `--spechla_loh_cmd` | In-container script paths for SpecHLA's extract/type/LOH stages | `/opt/SpecHLA/script/ExtractHLAread.sh`, `bash /opt/SpecHLA/script/whole/SpecHLA.sh`, `perl /opt/SpecHLA/script/cal.hla.copy.pl` |
| `--spechla_ref_build` | Reference build tag SpecHLA expects | `hg38` |

Before a real run, confirm these three paths actually resolve inside the
`.sif` (`params.containers.spechla`, already a local file on Jayden's
scratch): `singularity exec <spechla.sif> bash -lc 'find / -xdev -iname "ExtractHLAread.sh" -o -iname "SpecHLA.sh" -o -iname "cal.hla.copy.pl"'`.

## 10. Containers, at a glance

Most containers this pipeline references are public images/`.sif` URIs in
`nextflow.config`'s `containers {}` block, pulled automatically by
[`bin/prefetch_containers.sh`](../bin/prefetch_containers.sh) (run once
from a Gadi login node, see [`docs/RUNNING_ON_GADI.md`](RUNNING_ON_GADI.md)).
Three need a one-off manual build first, since no suitable public image
exists:

| Container | Recipe | Used by |
|---|---|---|
| `--ervcaller_container` | `containers/ervcaller/Dockerfile` | `--pipelines erv_dna` |
| `--agfusion_container` | `containers/agfusion/Dockerfile` | `--pipelines gene_fusion` |
| `params.splicing_tools_container` | `containers/splicing-tools/Dockerfile` | `--pipelines splicing` |

Four more (`star`, `picard`, `gatk`, `samtools`) and two (`sequenza`,
`spechla`) are already local `.sif` files on Jayden's Gadi scratch and
need no pulling or building at all.
