# External reference data and containers

Every file the pipeline needs that is not in this repository, where to get
it, and the `--param` that passes it in. `nextflow.config` holds the
defaults; this page maps onto it.

Everything must be **GRCh38 with `chr`-prefixed contig names**. Mixing
builds (or `1` vs `chr1`) does not error; it silently loses calls.

## 1. Genome

| Param | What | Where | Default |
|---|---|---|---|
| `--genome_fasta` | GRCh38 primary-assembly FASTA, uncompressed | [GENCODE release 46](https://www.gencodegenes.org/human/release_46.html): `GRCh38.primary_assembly.genome.fa.gz` | `/scratch/jo11/GRCh38.fa/GRCh38.primary_assembly.genome.fa` |

The primary assembly has no ALT or HLA contigs, which is what the
HLA-region extraction for HLA-HD assumes (`--hlahd_regions`).

`.fai`, `.dict` and the BWA index (`.amb .ann .bwt .pac .sa`) are reused if
they sit beside the FASTA (the `.dict` named `<fasta without .fa>.dict`).
Missing ones are built at the start of each run; BWA indexing alone takes
a few CPU-hours, so build it once on Gadi:

```bash
bwa index GRCh38.primary_assembly.genome.fa
samtools faidx GRCh38.primary_assembly.genome.fa
gatk CreateSequenceDictionary -R GRCh38.primary_assembly.genome.fa
```

## 2. GATK resource bundle (hg38)

All from the Broad public bucket `gs://gcp-public-data--broad-references/hg38/v0/`
(GATK [Resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)
page has the current location). Each needs its `.tbi` beside it; the
pipeline checks at launch and stops if one is missing
(`gatk IndexFeatureFile -I <vcf>` creates it).

| Param | Used by | File |
|---|---|---|
| `--known_dbsnp` | BQSR | `Homo_sapiens_assembly38.dbsnp138.vcf.gz` (default `/scratch/jo11/NeoTrio_RNA/dbsnp_138.hg38.vcf.gz`) |
| `--known_mills` | BQSR | `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` (default under `/scratch/jo11/NeoTrio_RNA/`) |
| `--known_1000g` | BQSR | `1000G_phase1.snps.high_confidence.hg38.vcf.gz` (default under `/scratch/jo11/NeoTrio_RNA/`) |
| `--germline_resource` | Mutect2 (required) | `somatic-hg38/af-only-gnomad.hg38.vcf.gz` |
| `--panel_of_normals` | Mutect2 (recommended) | `somatic-hg38/1000g_pon.hg38.vcf.gz`, or a PoN built from your own normals with `CreateSomaticPanelOfNormals` |
| `--contamination_sites` | GetPileupSummaries (recommended) | `somatic-hg38/small_exac_common_3.hg38.vcf.gz` |
| `--target_intervals` | Mutect2, GetPileupSummaries (WES only) | Your capture kit's BED for GRCh38 |

Without `--panel_of_normals` or `--contamination_sites` the run still
goes ahead, with a warning at launch.

## 3. VEP cache and plugins

| Param | What | Where |
|---|---|---|
| `--vep_cache` | Offline VEP cache, `homo_sapiens/111_GRCh38` (must match the VEP 111 container; change `--vep_cache_version` only together with the container) | `vep_install --AUTO c --SPECIES homo_sapiens --ASSEMBLY GRCh38 --CACHE_VERSION 111 --CACHEDIR <dir>`, or the [Ensembl cache download](https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html) |
| `--vep_plugins` | Directory holding `Wildtype.pm` and `Frameshift.pm` | `pvacseq install_vep_plugin <dir>` (inside the pVACtools container), which gives the versions pVACtools was tested with |

## 4. HLA-HD

HLA-HD is free for academic use but has to be downloaded after
registering, so there is no public image. Build `containers/hlahd/` from
your copy (see [its README](../containers/hlahd/README.md)), convert to a
`.sif`, and set `--hlahd_container` (or `hlahd_container` in the params
file). The HLA dictionary ships inside the HLA-HD download; its IMGT
version is whatever that release bundled.

Patients with `hla_alleles` in the samplesheet do not need HLA-HD; if
every patient has them, the container need not exist.

## 5. pVACtools

`griffithlab/pvactools:7.1.3` includes MHCflurry (with models), and
standalone IEDB MHC class I and II at `/opt/iedb`, which NetMHCpan and
NetMHCIIpan use (`--pvacseq_iedb_install_dir`). Nothing else to download.

## 6. Containers

| `containers.` key | Image | Notes |
|---|---|---|
| `fastqc`, `trimgalore`, `vep` | depot.galaxyproject.org Biocontainers | as before |
| `bwa`, `samtools`, `bcftools`, `whatshap` | Seqera community images (bwa+samtools, samtools 1.24, bcftools 1.23.1, WhatsHap 2.8) | the same images nf-core/modules pins; `bwa` holds samtools too because alignment pipes into `samtools sort` |
| `picard`, `gatk` | local `.sif` on `/scratch/jo11/singularity_images/` (GATK 4.5.0.0) | no pull needed |
| `bam_readcount` | `docker://mgibio/bam_readcount_helper-cwl:1.1.1` | the image pVACtools' docs point to; helper arguments confirmed from the script inside this tag |
| `vatools` | `docker://griffithlab/vatools:5.1.1` | `vcf-readcount-annotator` |
| `pvactools` | `docker://griffithlab/pvactools:7.1.3` | |

`bin/prefetch_containers.sh` pulls every public one into the shared
Singularity cache from a login node.
