# Architecture

Which of the steps below run is chosen with `--run`
(see [`RUNNING.md`](RUNNING.md)). Steps 1 and 3 are core and always run,
along with the HLA-LOH result and the Sequenza values read from the
samplesheet. Steps 2, 4, 5 and 8 to 10 are the `somatic` analysis, step 6
is `markers`, and step 7 is `fusion`. Step 11 runs either way and reports
what the run produced.

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
4. WhatsHap phases them against the tumour BAM with `--tag HP --indels
   --ignore-read-groups`, the settings of the Gadi PBS phasing script,
   using `/scratch/jo11/singularity_images/whatshap.sif`. pVACtools'
   recipe uses GATK 3.6 ReadBackedPhasing, which no current GATK ships
   (and hangs on Gadi). pVACseq 7.1.3 reads only the HP tag
   (`pvactools/lib/proximal_variant.py`); with WhatsHap's default PS tags
   it would treat every heterozygous neighbour as out of phase and use
   only homozygous ones.
5. VEP with the same options as the main VCF, then bgzip and tabix.

pVACseq only applies proximal missense variants; candidates near other
proximal variant types still deserve a look in IGV.

Published: `<id>/pvacseq_input/<id>.phased_proximal.vcf.gz(.tbi)`.

## 5. NEOANTIGENS (`workflows/neoantigens.nf`)

1. VEP on the PASS somatic VCF with pVACseq's required options and the
   Wildtype and Frameshift plugins.
2. bam-readcount on the tumour and normal DNA BAMs and on any RNA BAMs in
   the samplesheet (`rna_baseline_bam`, `rna_week6_bam`), via
   `bam_readcount_helper.py` from `mgibio/bam_readcount_helper-cwl:1.1.1`,
   as pVACtools' "Adding coverage data to your VCF" describes (min base
   quality 20 by default). RNA is counted at the ALT alleles of the
   `tumor` genotype column.
3. `vcf-readcount-annotator` (VAtools) writes the DNA counts into AD/DP/AF
   of the `tumor` and `normal` columns, which pVACseq uses for its DNA
   VAF/depth filters and its clonality tiering, and the baseline RNA
   counts into RAD/RDP/RAF of the `tumor` column. pVACseq reads one RNA
   sample, so only baseline RNA reaches it. With baseline RNA present,
   pVACseq's defaults apply: `filtered.tsv` keeps only candidates with
   tumour RNA depth > 10 and RNA VAF > 0.25 (`--trna-cov`, `--trna-vaf`;
   change them through `--pvacseq_extra_args`), and the aggregated report
   tiers a candidate with RNA VAF 0 as NoExpr. `all_epitopes.tsv` and the
   aggregated report keep every candidate. Week 6 RNA counts go into a
   separate copy of the VCF (not passed to pVACseq).
4. pVACseq with the HLA alleles, the phased proximal VCF and
   `--tumor-purity` from the samplesheet. Default predictors are
   MHCflurry, MHCflurryEL, NetMHCpan, NetMHCpanEL, NetMHCIIpan and
   NetMHCIIpanEL, using the standalone IEDB install inside the pVACtools
   image (`--pvacseq_iedb_install_dir /opt/iedb`).
5. When the patient has week 6 RNA, `bin/add_rna_timepoint.py` adds
   `RNA Depth (week 6)` and `RNA VAF (week 6)` to each aggregated report,
   matching on pVACseq's variant ID. pVACseq's own `RNA Depth` and
   `RNA VAF` columns are the baseline.
6. `bin/annotate_ccf.py` adds copy-number and CCF columns to each
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

## 6. EXPRESSION (`workflows/expression.nf`, optional)

Runs when `--gtex_db` or `--tcga_db` is set; see
[`EXPRESSION.md`](EXPRESSION.md) for the reference tables and the output
columns.

1. featureCounts over `--rna_gtf` for every RNA BAM in the samplesheet,
   then `bin/gene_tpm.py` for TPM from counts and union-exon length
   (`RPK = counts / kb`, `TPM = RPK / (sum RPK / 1e6)`) plus a
   within-sample expression quartile. Skipped without `--rna_gtf`.
   `--rna_strandedness` mirrors featureCounts `-s` and has no safe
   default.
2. `bin/annotate_expression.py` adds GTEx normal-tissue columns (matched
   tissue, median and max across tissues, how many tissues express the
   gene, and a `not_expressed` / `testis_restricted` / `tissue_restricted`
   / `broadly_expressed` call), the matched TCGA cohort's TPM, this
   sample's baseline TPM, and an `Expression Flag`. The patient's own RNA
   is the evidence for the flag when present, else the TCGA cohort, which
   is what makes the step useful on a DNA-only patient.

3. `bin/gene_panel_report.py` writes a per-patient status table for 94
   genes in 11 pathways: class I and II presentation, the non-classical
   class Ib genes (HLA-E, HLA-F, HLA-G) with MR1 and CD1A-E, the
   MHC class I-related NKG2D ligands (MICA, MICB, ULBP1-3, RAET1E/G/L),
   their regulators (CIITA, NLRC5, RFX, NFY), the immunoproteasome, the
   interferon route (IFNG, IFNGR1/2, JAK1, JAK2, STAT1, IRF1/2) and the canonical melanoma
   drivers (BRAF, NRAS, KRAS, HRAS, NF1, MAP2K1/2, RAC1, PPP6C, KIT,
   GNAQ/11, TP53, CDKN2A, CDK4, CCND1, RB1, MDM2, PTEN, PIK3CA, AKT1,
   ARID2, IDH1, STK19, DDX3X, SF3B1, BAP1, EIF1AX, MITF, TERT, CTNNB1).
   Each gene gets an expression call against the GTEx/TCGA range, copy
   number from the Sequenza segments covering it (gene coordinates come
   from `--rna_gtf`, with minor copy number 0 reported as LOH), and the
   somatic consequences VEP called in it, with truncating ones marked.
   Each gene carries an `Adverse Direction`: most are `loss`, but HLA-E
   and HLA-G are inhibitory ligands (NKG2A, ILT2/ILT4), so for them gain
   is the finding and `Panel Flag` reports `INHIBITORY_GAIN`. Otherwise
   the flag is `LOSS_OF_FUNCTION` / `ALTERED` / `OK` / `NA`. This is
   the step that catches a tumour that cannot present antigen at all
   (B2M or HLA loss, JAK1/JAK2 truncation), which no per-candidate report
   would show. Replace the panel with `--gene_panel`. It is not an
   HLA-LOH caller: copy number here is gene level, not allele specific.

Reference tables are built once per cohort by
`bin/build_expression_db.py` from public GTEx/TCGA downloads; no reference
data ships in this repository.

Published: `<id>/expression/<id>.<timepoint>.tpm.tsv`,
`<id>/pvacseq_expression/<id>.{MHC_I,MHC_II,Combined}.aggregated.ccf.expression.tsv`,
`<id>/gene_panel/<id>.gene_panel.tsv`.

## 7. FUSION (`workflows/fusion.nf`, optional)

Runs with `--fusion true` for every patient that has a baseline tumour
RNA FASTQ pair or an existing STAR-Fusion table: STAR-Fusion (or that
table) -> AGFusion -> pVACfuse, on the same HLA alleles, pVACtools image
and IEDB install as pVACseq. See [`FUSION.md`](FUSION.md) for the
reference data it needs, why it is off by default, and how to read the
candidates.

Published: `<id>/fusion/star_fusion/`, `<id>/fusion/agfusion/`,
`<id>/pvacfuse/`.

## 8. PRIORITIZE (`modules/local/pvactools/prioritize.nf`, on by default)

See [`PRIORITISATION.md`](PRIORITISATION.md).

Runs last on the candidate tables, on the expression-annotated ones where
step 6 ran and the CCF ones otherwise, joined with the patient's HLA-LOH
table where the samplesheet gave one. `bin/prioritize_candidates.py`
applies every criterion once and records what failed: `HLA_LOST`,
`LOW_BINDING`, `NO_EXPRESSION` and `NO_RNA_SUPPORT` reject a candidate;
`SELF_LIKE`, `SUBCLONAL`, `NORMAL_TISSUE`, the thin-evidence flags and the
missing-data flags are cautions. It also adds `Agretopicity` (IC50 MT/WT)
and `DAI`, which nothing upstream computes. Nothing is dropped: the order
is a documented sort and every row keeps its rank, flags and reason.

Published: `<id>/prioritized/`.

## 9. LONG_PEPTIDES (`modules/local/pvactools/long_peptides.nf`, on by default)

See [`LONG_PEPTIDES.md`](LONG_PEPTIDES.md).

`pvacseq generate_protein_fasta` regenerates the mutant protein context
from the same VCF and phasing pVACseq used, since no candidate table
carries flanking sequence. `bin/design_long_peptides.py` then places a
window of up to `--long_peptide_max_length` on each ranked candidate,
locating the mutation by comparing each MT/WT pair rather than by parsing
fasta headers, anchoring frameshifts at the junction, trimming awkward
termini while the mutation keeps its flanks, and merging neighbouring
mutations that fall inside one window. With `--vaccine_vector true`,
pVACvector orders the result into a string-of-beads construct.

Published: `<id>/long_peptides/`, `<id>/vaccine_vector/`.

## 10. CLONAL_STRUCTURE (`modules/local/report/clonal_structure.nf`, on by default)

`bin/clonal_structure.py` cuts the sorted CCFs into clusters wherever
consecutive values are more than `--clonal_gap` apart, merges clusters
below `--clonal_min_cluster` variants into their nearest neighbour, and
reports for each cluster how many of its variants are still detected in
the week 6 RNA. `not_detected` (covered, no variant reads) is kept
distinct from `no_coverage`. No tree: that needs several regions or DNA
timepoints.

Published: `<id>/clonal_structure/`.

## 11. REPORT (`workflows/report.nf`, on by default)

See [`REPORT.md`](REPORT.md).

1. `bin/report_figures.py`, for every patient: the candidate landscape,
   the flag summary and the HLA-LOH overview. The report draws nothing
   itself, so each figure can be regenerated from an existing results
   directory on its own. Turn it off with `--report_figures false`.
2. `bin/rna_support.py`, for every patient with at least one RNA BAM:
   each non-synonymous PASS variant is called RNA-supported at a
   timepoint when at least `--rna_support_min_alt_reads` RNA reads carry
   the alt allele, and grouped as expressed both / week 6 only / baseline
   only / DNA only. The same split is shown for clonal variants (CCF at or
   above `--clonal_ccf`) and for clonal MHC I and II neoantigens (best
   peptide under the `--neoag_class_*` IC50 or percentile thresholds).
   Writes a per-variant table, a summary table and the figure.
3. `bin/patient_report.py`, for every patient: a LaTeX report laid out
   after the NeoDisc templates but filled from this pipeline: TCGA SKCM
   genomic subtype and driver status, TMB and UV signature, purity, ploidy
   and clonality, HLA alleles and gene-level HLA/B2M copy number, the gene
   panel, the RNA-support figure, GTEx/TCGA context and the top
   candidates. NeoDisc sections with no equivalent are listed as not
   assessed.
4. `pdflatex` (twice) unless `--report_pdf false`.

The report also reads the pVACfuse directory when the fusion branch ran,
and distinguishes "branch did not run" from "ran, but no fusion had a
coding partner pair".

Published: `<id>/report/<id>.report.{tex,pdf}`,
`<id>/report/<id>.rna_support.{summary.tsv,variants.tsv,pdf,png}`.

## Not covered

- No RNA alignment: RNA comes in as aligned BAMs. Gene-level TPM is
  computed in step 6 but is not written back into the VCF, so pVACseq
  itself still gets no gene or transcript expression (GX/TX) and its own
  expression filters and Allele Expr tier check do not apply; expression
  arrives as extra columns afterwards instead.
- No cross-timepoint expression change: week 6 TPM is published but not
  compared to baseline.
- No in-pipeline Sequenza and no in-pipeline SpecHLA; purity, ploidy,
  segments and the HLA-LOH result come from existing runs through the
  samplesheet.
- No tumour evolution tree: one DNA timepoint gives CCF clusters, not a
  phylogeny.
- No immunopeptidomics: nothing here checks by mass spectrometry that a
  predicted peptide is actually presented.
- No consensus variant calling: Mutect2 with a panel of normals and the
  read-orientation model is the only caller.
- One FASTQ pair per sample.
- Tests check wiring (stubbed runs, unit tests on the helper scripts),
  never biological accuracy. No real-patient run yet.
