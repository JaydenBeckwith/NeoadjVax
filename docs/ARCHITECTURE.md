# Architecture

One path, per patient: tumour and matched-normal DNA to pVACseq
neoantigen candidates with clonality. Entry point `main.nf`; one
sub-workflow per stage in `workflows/`; processes in `modules/local/`;
helper scripts in `bin/` (on `PATH` inside every task).

## 1. PREPROCESS_DNA (`workflows/preprocess_dna.nf`)

FASTQ pairs go through FastQC (raw and trimmed), Trim Galore and BWA-MEM,
with the read group set at alignment. BAM inputs skip those and get
`AddOrReplaceReadGroups` instead. Both then go through Picard
MarkDuplicates and GATK BQSR (dbSNP, Mills, 1000G known sites).

The read-group sample name is always the literal `tumor` or `normal`.
Mutect2, bam-readcount, WhatsHap and pVACseq all identify samples by
that name, so fixing it here means no renaming later.

Published: `<id>/bam/<id>.{tumor,normal}.recal.bam(.bai)`, `<id>/qc/`, `<id>/fastqc/`.

## 2. SOMATIC_CALLING (`workflows/somatic_calling.nf`)

GATK's tumour/matched-normal best practice:

- Mutect2 with `--germline-resource` (gnomAD AF-only) and, if given,
  `--panel-of-normals`, one task per contig in `--scatter_chromosomes`
  (chr1-22, X, Y by default). With `--target_intervals` each task is the
  intersection of its contig and the capture targets.
- The scatters are merged (`MergeVcfs`, `MergeMutectStats`) and the
  F1R2 counts become a read-orientation model (`LearnReadOrientationModel`),
  which removes FFPE and OxoG artefacts.
- With `--contamination_sites`, `GetPileupSummaries` on both samples and
  `CalculateContamination` give a contamination estimate and tumour
  segmentation.
- `FilterMutectCalls` uses all of the above.
- PASS records are left-aligned and split to one ALT per record
  (`bcftools norm -m -any`): pVACtools asks for decomposed records before
  bam-readcount, and the main pVACseq VCF and the phasing VCF need
  identical records.

Contigs outside `--scatter_chromosomes` are never called, which replaces
the old standard-chromosome filter.

Published: `<id>/somatic_variants/` (filtered VCF, PASS VCF, contamination
and filtering stats).

## 3. HLA_TYPING (`workflows/hla_typing.nf`)

HLA-HD on the normal BAM, so tumour HLA loss cannot bias typing. Only
reads aligned to the MHC (`--hlahd_regions`, chr6:28.5-33.5 Mb) plus
unplaced unmapped reads are extracted to FASTQ first; HLA-HD on a whole
WGS FASTQ gives the same answer far more slowly.

`bin/parse_hlahd_result.py` turns `<id>_final.result.txt` into pVACseq
names at two-field resolution:

- class I `HLA-A*02:01`, `HLA-B…`, `HLA-C…`
- DR beta chains `DRB1*15:01` (and DRB3/4/5 when typed)
- DQ and DP as alpha-beta pairs, `DQA1*01:02-DQB1*06:02`. Phase between
  the alpha and beta loci is unknown from short reads, so every cis and
  trans combination is listed.

Null alleles are dropped. A patient with no usable A, B or C call fails
with a message rather than going ahead with a partial class I set.
A samplesheet `hla_alleles` value skips HLA-HD for that patient.

Published: `<id>/hla_typing/` (raw HLA-HD result, allele table, pVACseq string).

## 4. PROXIMAL_VARIANTS (`workflows/proximal_variants.nf`)

pVACseq builds each peptide as if the somatic variant were the only
change nearby unless it is given a phased VCF of proximal variants
(`--phased-proximal-variants-vcf`). This follows pVACtools'
"Creating a phased VCF of proximal variants" recipe:

1. A BED of `--phasing_window_bp` (200 bp) either side of every PASS somatic site.
2. HaplotypeCaller on the normal, only inside those windows (minutes,
   not a whole-genome germline run). Calls with QUAL below
   `--germline_min_qual` are dropped.
3. The tumour column of the somatic VCF and the germline calls (renamed
   to `tumor` to match the tumour BAM) are combined and sorted.
4. WhatsHap phases them against the tumour BAM with `--tag HP`.
   pVACtools' recipe uses GATK 3.6 ReadBackedPhasing, which no current
   GATK ships; WhatsHap's HP tags are the same field pVACseq reads
   (`pvactools/lib/proximal_variant.py`) to decide what is in phase.
5. VEP with the same options as the main VCF, then bgzip and tabix.

pVACseq only applies proximal missense variants; candidates near other
proximal variant types still deserve a look in IGV.

Published: `<id>/pvacseq_input/<id>.phased_proximal.vcf.gz(.tbi)`.

## 5. NEOANTIGENS (`workflows/neoantigens.nf`)

1. VEP on the PASS somatic VCF with pVACseq's required options and the
   Wildtype and Frameshift plugins.
2. bam-readcount on the tumour and normal BAMs, via
   `bam_readcount_helper.py` from `mgibio/bam_readcount_helper-cwl:1.1.1`,
   as pVACtools' "Adding coverage data to your VCF" describes (min base
   quality 20 by default).
3. `vcf-readcount-annotator` (VAtools) writes the counts into AD/DP/AF of
   the `tumor` and `normal` columns, which pVACseq uses for its DNA
   VAF/depth filters and its clonality tiering.
4. pVACseq with the HLA alleles, the phased proximal VCF and
   `--tumor-purity` from the samplesheet. Default predictors are
   MHCflurry, MHCflurryEL, NetMHCpan, NetMHCpanEL, NetMHCIIpan and
   NetMHCIIpanEL, using the standalone IEDB install inside the pVACtools
   image (`--pvacseq_iedb_install_dir /opt/iedb`).
5. `bin/annotate_ccf.py` adds copy-number and CCF columns to each
   aggregated report. pVACseq's own clonality call is
   `purity × 0.5` as the clonal VAF, i.e. it assumes every variant sits
   on a diploid segment; this step corrects for copy number using the
   Sequenza segment at the variant, else the sample ploidy:

   ```
   u   = VAF × (purity × CN_tumour + CN_normal × (1 − purity)) / purity
   m   = max(1, round(u)), capped at the major copy number when known
   CCF = min(1, u / m)
   ```

   It is a point estimate without a confidence interval.

Published: `<id>/pvacseq/` (the full pVACseq output),
`<id>/pvacseq_input/<id>.somatic.vep.readcounts.vcf.gz`,
`<id>/pvacseq_ccf/<id>.{MHC_I,MHC_II,Combined}.aggregated.ccf.tsv`.

## Not covered

- No RNA: pVACseq gets no expression (TX/GX) or RNA readcounts, so its
  expression filters do not apply. Adding an RNA FASTQ column for
  quantification and RNA bam-readcount is the natural next step.
- No HLA-LOH and no in-pipeline Sequenza; purity, ploidy and segments
  come from existing Sequenza runs.
- One FASTQ pair per sample.
- Tests check wiring (stubbed runs, unit tests on the two helper
  scripts), never biological accuracy. No real-patient run yet.
