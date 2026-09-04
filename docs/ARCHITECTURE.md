# NeoadjVax architecture

Neoadjuvant melanoma vaccine neoantigen pipeline: DNA + RNA variant calling
feeding pVACseq (the part with an existing design — see
`[[neoantigen-score]]`), plus three new discovery branches — splicing,
gene fusion, and ERV neoantigens — scaffolded at roughly equal depth.

Orchestration: Nextflow DSL2. Execution: PBS Pro + Singularity on NCI Gadi
(`-profile gadi`), reusing Jayden's existing `.sif` images and PBS project
settings where they already exist, matching how `spliceai-variant-pipeline`
and the original RNA script already run.

**Containers: Singularity, not Docker**, under `-profile gadi` or
`-profile singularity` — Gadi's compute nodes have no Docker (no root on
HPC). Every image in `nextflow.config`'s `containers` block, including the
`docker://griffithlab/pvactools` reference, is pulled and converted to a
`.sif` by Singularity itself; Docker is never invoked. `-profile docker`
exists only as an optional local-dev convenience (mirrors the original DNA
script's Docker calls) for iterating off Gadi — not how this runs in
production. See the block's own comment in `nextflow.config` for detail.

**Already-aligned input**: both branches accept pre-aligned BAMs and skip
the alignment step — DNA per-row via `tumor_bam`/`normal_bam` columns
(skips FastQC/Trim Galore/BWA-MEM, still runs AddReadGroups+MarkDuplicates),
RNA via `--rna_skip_realignment true` (skips BAM→FASTQ→STAR entirely,
pipeline-wide — only safe if those BAMs already match this pipeline's STAR
settings). See `assets/samplesheet_schema.md`.

## Pipeline graph

```mermaid
flowchart TB
    subgraph DNA["DNA branch — dna_variant_calling.nf"]
        d1["Trim Galore + FastQC"] --> d2["BWA-MEM align"]
        d2 --> d3["AddReadGroups → MarkDuplicates"]
        d3 --> d4["Mutect2 (tumor vs normal)"]
        d4 --> d5["FilterMutectCalls"]
        d5 --> d6["Filter to standard chroms"]
    end

    subgraph RNA["RNA branch — rna_variant_calling.nf (per timepoint)"]
        r1["BAM → FASTQ (or FASTQ direct)"] --> r2["STAR two-pass align"]
        r2 --> r3["AddReadGroups → MarkDuplicates"]
        r3 --> r4["SplitNCigarReads"]
        r4 --> r5["BQSR"]
        r5 --> r6["HaplotypeCaller"]
        r6 --> r7["VariantFiltration"]
    end

    d6 --> match["RNA_SUPPORT_FILTER<br/>(DNA loci ∩ RNA support, per timepoint)"]
    r7 --> match
    match --> vep["VEP annotate<br/>(Wildtype + Frameshift plugins)"]
    normalbam["Normal DNA reads"] --> hla["OptiType HLA typing<br/>(class I; class II TODO)"]
    hla --> pvacseq["pVACseq"]
    vep --> pvacseq
    pvacseq --> core["Core neoantigen calls<br/>(per patient × timepoint)"]

    subgraph SPLICE["splicing_neoantigens.nf — stub past ingestion"]
        s1["Existing SpliceAI pipeline output<br/>(reused, not re-run)"] --> s2["Filter by delta score"]
        s2 --> s3["Splice→peptide — TODO"]
    end

    subgraph FUSION["fusion_neoantigens.nf"]
        f1["STAR-Fusion"] --> f2["AGFusion annotate — TODO"]
        f2 --> f3["pVACfuse"]
    end

    subgraph ERV["erv_neoantigens.nf"]
        e1["Telescope quantification"] --> e2["ERV→peptide — TODO,<br/>design decision open"]
    end

    style s3 fill:#4a2a2a,stroke:#c66
    style e2 fill:#4a2a2a,stroke:#c66
    style f2 fill:#4a2a2a,stroke:#c66
```

Red boxes = explicit TODO stubs (process exits 1 with a comment explaining
the open decision) — everything else is ported from your working scripts
or is new code built to the same standard.

## Where your original scripts went

| Original | Where it landed |
|---|---|
| `RNA_variant_pipeline.sh` (PBS) — STAR index, per-sample BAM→FASTQ→STAR→GATK RNA variant calling | `modules/local/rna/*.nf` + `workflows/rna_variant_calling.nf`. Same steps, same tool flags (e.g. `--dont-use-soft-clipped-bases`, `FS > 30.0 \|\| QD < 2.0`, `sjdbOverhang 149`) — checkpoint-by-file-existence replaced by Nextflow's own `-resume` cache, which does the same job without the manual `checkpoint_run` wrapper. |
| DNA pVACseq script (Python/Docker) — trim/align/Mutect2/VEP/pVACseq | Split across `modules/local/dna/*.nf` (through `FILTER_STANDARD_CHROMS`) and `workflows/pvacseq_core.nf` (VEP + pVACseq) — see below for why it's split. Docker calls became Singularity containers per `nextflow.config`'s `containers` block; swap back to `-profile docker` for local iteration if that's still useful. |

**One structural change from the original DNA script**: VEP annotation and
`pvacseq run` no longer happen immediately after Mutect2. They now happen
*after* `RNA_SUPPORT_FILTER`, on the smaller RNA-supported subset —
because `[[neoantigen-score]]`'s design is "RNA-supported VCF subset then
run through pVACseq," not "every raw DNA call through pVACseq." The
original script (DNA-only) had no RNA step to wait for, so this wasn't a
choice it needed to make.

## New code that didn't exist in either script

- **`bin/match_rna_support.py`** + `modules/local/matching/rna_support_filter.nf`
  — the RNA-support check itself. Matches DNA/RNA calls by exact
  `(CHROM, POS, REF, ALT)`, optionally thresholding on RNA ALT-allele read
  depth (`--min_rna_alt_reads`, default 1). This is the one part of the
  core backbone that's genuinely new rather than ported — worth a close
  look before trusting it on real data.
- **HLA typing** (`modules/local/hla_typing/optitype.nf`) — the DNA script
  took `--hla` as a hand-supplied list; nothing typed it. OptiType now runs
  on the **normal/germline** DNA BAM specifically, not tumor, to avoid
  bias from tumor HLA-LOH — which `[[loh-analysis]]` is already tracking
  separately, so typing off tumor reads here would fold that effect into
  the neoantigen calls unintentionally. Class II typing
  (`arcashla.nf`) is a TODO stub — `[[neoantigen-score]]` calls for both
  classes but neither script did class II at all.
- **Two-samplesheet input design** (`--dna_samplesheet` / `--rna_samplesheet`,
  joined on `patient_id`) — because DNA is one baseline draw per patient and
  RNA is longitudinal (PRE/ED1/ED2/CLND). See `assets/samplesheet_schema.md`.

## Open decisions — need your input, not guessed at

1. **Splicing → peptide.** No standard tool converts a novel splice
   junction into a candidate peptide. Options: NeoSplice (purpose-built,
   narrower toolchain) vs. custom ORF translation built on the
   IsoformSwitchAnalyzeR output `[[neoadjuvant-splicing]]` already
   produces. See `modules/local/splicing/splice_to_peptide.nf`.
2. **ERV → peptide.** Same shape of problem, less prior art. Options:
   treat expressed loci as novel ORFs (reuses pVACtools infrastructure,
   hand-rolled and unvalidated scoring) vs. a dedicated published
   TE-neoantigen scorer (more defensible, more work). See
   `modules/local/erv/erv_to_peptide.nf`.
3. **Fusion branch's FASTQ source.** `fusion_neoantigens.nf` currently only
   reads `rna_fastq_r1/r2` from the samplesheet; most rows will only have
   `rna_bam` per the schema. Needs either its own BAM→FASTQ step or to
   share `RNA_VARIANT_CALLING`'s output — not decided yet.
4. **ERV branch's alignment.** Telescope wants multi-mapping reads that
   `RNA_VARIANT_CALLING`'s `--outFilterMultimapNmax 2` (carried over from
   your original script) discards — so it can't just reuse that BAM as-is.
   Needs its own STAR pass with relaxed multimapping, not yet built.
5. **RNA_variant_pipeline.sh tail was truncated** in what got pasted in —
   it references `${bn}.braf_allelic_counts.tsv` in its cleanup step and a
   trailing "PIPELINE SUMMARY" / `COLLECTOREOF` heredoc that never actually
   runs anything. The BRAF V600 allelic-count extraction itself isn't
   reflected anywhere in this scaffold — flagged rather than guessed at;
   send the missing piece if you still want it in the pipeline.
6. **DNA-RNA matching depth threshold.** `--min_rna_alt_reads` defaults to
   1 (any RNA support counts). Whether that's the right bar — vs. a higher
   depth threshold, or a VAF-based check instead of a raw read count — is
   worth deciding once there's real data to look at.

## What's real vs. stubbed, at a glance

**Fully ported / working shape**, ready to run once reference paths, VEP
cache/plugins, and containers are confirmed on Gadi:
DNA variant calling, RNA variant calling, RNA-support matching, VEP
annotation, pVACseq core.

**New but complete modules**, not yet run against real data:
OptiType HLA typing (class I), SpliceAI-output filtering, STAR-Fusion +
pVACfuse, Telescope ERV quantification.

**Explicit TODO stubs** (exit 1, comment explains the decision needed):
splice→peptide, ERV→peptide, AGFusion annotation, arcasHLA class II.

## Next steps

1. Confirm reference/container paths in `nextflow.config` actually resolve
   on Gadi (`params.genome_fasta`, VEP cache/plugin dirs, CTAT resource lib,
   ERV annotation GTF — none of these have real paths filled in yet).
2. Fill in `--dna_samplesheet` / `--rna_samplesheet` for one real patient
   and do a dry run of just the core backbone
   (`--run_splicing_neoantigens false --run_fusion_neoantigens false
   --run_erv_neoantigens false`) before touching the new branches.
3. Pick a direction on the three open peptide-generation decisions above —
   those block real progress on splicing/fusion/ERV more than any missing
   code does.
