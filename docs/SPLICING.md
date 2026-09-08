# Splicing peptide candidates

`--pipelines splicing` implements somatic DNA VCF → SpliceAI → splice2neo
junction/transcript annotation, combined with RNA RegTools/LeafCutter evidence,
then GENCODE/normal filtering and splice2neo CDS/ORF peptide-context prediction.
One baseline DNA sample fans out to **every RNA timepoint** for that patient.
No fixed labels, baseline RNA or minimum number of timepoints are required.

## Input contracts

DNA, one unique row per patient:

~~~csv
patient_id,spliceai_vcf
P1,/data/P1.somatic.spliceai.vcf.gz
~~~

Or supply `patient_id,dna_vcf` for fresh SpliceAI annotation. With
`--pipelines dna_variant_calling,splicing`, use the usual tumour/normal BAM or
FASTQ columns to call somatic variants first. `spliceai_vcf` takes precedence
for splicing only. VCFs must contain tumour-versus-normal somatic calls:
only FILTER=PASS records are considered (not FILTER=.). No HLA alleles are
needed for peptide-context generation. WES may miss intronic variants
outside capture regions.

RNA, one unique row per patient/timepoint, using existing LeafCutter:

~~~csv
patient_id,timepoint,leafcutter_counts,leafcutter_sample
P1,day0,/data/cohort_perind.counts.gz,P1_day0.bam
P1,day99,/data/cohort_perind.counts.gz,P1_day99.bam
~~~

`leafcutter_sample` selects one exact matrix header column and is required
for cohort matrices. Omit only for a single-sample file. Input must be the
**regtools-based, strand-aware** `*_perind.counts.gz`: header `chrom SAMPLE...`,
IDs `chr:start:end:clu_N_+/-`, and `numerator/denominator` counts. Older
strandless output and `perind_numers` are rejected. These coordinates already
describe 1-based **joined exon bases**; do not apply a BED offset again.
Zero counts in another cohort sample never provide support for this sample.

Alternatively supply `rna_bam` or paired `rna_fastq_r1,rna_fastq_r2`.
BAMs must be coordinate-sorted, with intact spliced `N` CIGARs and `NH` tags.
Only primary `NH=1` alignments are used. **Do not supply GATK
SplitNCigarReads/ApplyBQSR BAMs.** BAM indexing occurs within the task.
RNA VCFs are not junction evidence.

For BAMs, specify `rna_strand` in the row or `--splicing_rna_strand`:
`XS` uses aligner tags, `RF` means first-strand, `FR` second-strand.
Missing XS tags cannot support strand-specific junctions in XS mode.
FASTQs receive a dedicated STAR two-pass alignment with intron-motif XS
tags and unique mappings, even when RNA variant calling is also enabled.
Use a matching `--star_index_dir`, or set it to null in a params file to
build from FASTA/GTF. `splicing_star_overhang` defaults to 149; set to
read length minus one when building an index.

## References and normal filtering

Use the **same genome build and contig naming** for DNA, RNA, FASTA, GTF and
normal junctions. No liftover or automatic chromosome renaming occurs.
FASTA must be uncompressed; GTF must contain gene_id, gene_name, exons and
CDS. Reference preparation reports excluded GTF contigs absent from the
FASTA (e.g. patch/haplotype annotations with a primary-assembly genome).

Fresh annotation requires `--spliceai_annotation`: a matching custom
SpliceAI table, or explicitly `grch37`/`grch38`. The bundled annotation is
GENCODE v24, **not** the project's v46 GTF default. Harmonize annotation
versions for your study; they affect sensitivity and canonical definitions.
Unmapped SpliceAI genes fail rather than being silently discarded.

`--splicing_normal_junctions` is a curated, reference-matched normal/GTEx
panel: **one-column TSV**, header `junc_id`, values
`chromosome:joinedExonBase1-joinedExonBase2:strand`.
Select tissues and detection thresholds for your study and record the
panel's provenance. Raw GTEx files with unspecified coordinate conventions
are not accepted. Convert a panel according to its documented convention;
splice2neo's `bed_to_junc` supports compatible BED representations.

Normal filtering is required by default. For exploration only,
`--splicing_allow_missing_normal true` allows an absent panel and records
`normal_filter_applied=false` in results. A supplied empty panel is an
error, not an implicit opt-out. Absence from a panel does not prove
tumour specificity.

## Containers and running

R processes use the upstream
`ghcr.io/tron-bioinformatics/splice2neo:v0.6.14` image. Python adapters use
3.11.11. These suffice for precomputed inputs.

Fresh SpliceAI or BAM/FASTQ-based LeafCutter processing additionally needs
`--splicing_tools_container`. Build off Gadi:

~~~bash
docker build -t neoadjvax-splicing-tools:1.0 -f containers/splicing-tools/Dockerfile .
~~~

The recipe pins SpliceAI 1.3.1, TensorFlow 2.15.1, RegTools 1.0.0,
samtools/bcftools 1.20 and LeafCutter's clusterer commit
`2c9907ef66adf0bfb3092f0ceb6886ee5c046fbb`. Review
[SpliceAI's code/model licensing](https://github.com/Illumina/SpliceAI)
before use or redistribution. The custom image build and production
container run still require validation. Transfer a built SIF to Gadi and
pass its absolute path. `bin/prefetch_containers.sh` includes the two public
images, not this custom image; compute nodes must not need network access.

With precomputed inputs:

~~~bash
nextflow run main.nf -profile gadi --pipelines splicing \
  --dna_samplesheet dna.csv --rna_samplesheet rna.csv \
  --genome_fasta /refs/GRCh38.fa --gtf /refs/gencode.annotation.gtf \
  --splicing_normal_junctions /refs/normal_junctions.tsv --outdir results
~~~

For fresh annotation add `--spliceai_annotation /refs/spliceai.annotation.tsv`
and `--splicing_tools_container /images/neoadjvax-splicing-tools.sif`.
Use `dna_variant_calling,splicing` for raw DNA inputs. Standalone `splicing`
does not automatically enable variant calling.

## Outputs and scope

Per patient, `splicing/dna/` contains annotations when newly run, parsed
effects/QC and predicted junction/transcript hypotheses.
Per `patient/splicing/timepoint/`:

- `junction_evidence.tsv`: predictions, canonical/normal flags and RNA counts.
- `peptide_candidates.tsv`: supported noncanonical hypotheses, CDS/ORF annotations
  and `peptide_exported`; hypotheses lacking usable peptides remain visible.
- `peptide_contexts.fasta`: candidate IDs map back to the table.
- `summary.tsv`, `versions.txt`; LeafCutter evidence is under `leafcutter/`.

Defaults: per-effect SpliceAI score ≥0.5; ≥3 junction reads in the selected
sample; LeafCutter `-m 3` and ratio ≥0.001; peptide context
flanks of 13 wild-type amino acids and minimum exported length 8.
A high score for one effect does not qualify a different low-scoring effect.
Empty/no-support cases yield header-only tables and an empty FASTA.
Malformed inputs fail; SpliceAI QC records unannotated and low-score variants.
The pinned LeafCutter clusterer applies `-m` to individual junctions **and**
cluster totals, despite its option name. `leafcutter_min_reads=3` preserves
junctions at our RNA threshold; setting it to 50 is stringent per sample.
Its initial pooling also requires at least three reads, independently of `-m`.

This implements the supplied figure's **SpliceAI/LeafCutter intersection**
path, not EasyQuant targeted re-quantification. MMSplice, SplAdder, intron
retention, RNA-only peptide discovery, HLA binding/ranking, normal-proteome
peptide screening and experimental validation are not included. LeafCutter
per-sample clustering is not differential-splicing analysis; isolated and
low-support junctions may be missed. Peptides use annotated coding frames,
without phased nearby variants or full isoform reconstruction. Predictions
do not establish causality, expressed isoform identity or immunogenicity.

## Tests

~~~bash
python -B -m unittest discover -s tests -p 'test_*.py' -v
python tests/splicing/make_fixture.py /tmp/neoadjvax-splicing
python bin/prepare_splicing_inputs.py spliceai \
  --input /tmp/neoadjvax-splicing/dna.spliceai.vcf \
  --output /tmp/neoadjvax-splicing/effects.tsv --qc /tmp/neoadjvax-splicing/qc.json
python bin/prepare_splicing_inputs.py leafcutter \
  --input /tmp/neoadjvax-splicing/cohort_perind.counts.gz --sample P1_day0 \
  --output /tmp/neoadjvax-splicing/selected.gz --evidence /tmp/neoadjvax-splicing/evidence.tsv
Rscript tests/splicing/test_splice2neo.R /tmp/neoadjvax-splicing
# Linux/Gadi channel/staging test; stubs do not test tool execution:
nextflow run main.nf -c tests/splicing/stub.config \
  -params-file /tmp/neoadjvax-splicing/params.json -stub-run
~~~

The smoke test must produce P1/day0, P1/day99 and P2/resection outputs.
The fixture's SpliceAI annotations are synthetic, not model predictions.
The actual upstream clusterer can also be tested with
`python tests/splicing/test_leafcutter.py /path/to/leafcutter_cluster_regtools.py`
using the pinned commit above. It verifies both strands, coordinate offsets,
read counts and the upstream `-m` threshold behaviour.

### Verification status

Checked locally on 8 September 2026: 31 Python regression tests; real
splice2neo 0.6.14 translation on both strands, canonical/normal exclusions
and empty/no-CDS cases; and the pinned upstream LeafCutter clusterer on
synthetic BED12 junctions, including its per-junction read cutoff.

Full Nextflow orchestration (including the stub run), the fresh SpliceAI,
STAR and RegTools binaries, the custom container build and real-data
validation have not been executed here. Run the Linux/Gadi smoke test and
a reference-matched pilot before production use. The README SVG was also
checked for label overlap/clipping in light and dark modes; its canvas is transparent.

## Methods

See the [splice2neo workflow/API](https://tron-bioinformatics.github.io/splice2neo/articles/splice2neo_workflow.html),
[LeafCutter workflow](https://davidaknowles.github.io/leafcutter/articles/Usage.html),
[RegTools strand options](https://regtools.readthedocs.io/en/latest/commands/junctions-extract/)
and [Lang et al., 2024](https://doi.org/10.1093/bioadv/vbae080).
