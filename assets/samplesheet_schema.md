# Samplesheet

`--samplesheet` is a CSV with one row per patient. `main.nf` checks every
row at launch, before any job is submitted, and stops with a message
naming the patient and column if something is missing or malformed.

| column | required | notes |
|---|---|---|
| `patient_id` | yes | Unique; letters, digits, `.`, `_`, `-`; starts with a letter or digit. Used in every output path. |
| `tumor_r1`, `tumor_r2` | one of these or `tumor_bam` | Tumour DNA FASTQ pair (gzipped is fine). |
| `tumor_bam` | one of these or the FASTQ pair | Already-aligned tumour BAM (GRCh38, same reference as `--genome_fasta`). Skips trimming and BWA-MEM but still gets read groups reset, duplicate marking and BQSR. |
| `normal_r1`, `normal_r2` / `normal_bam` | as for tumour | Matched normal (blood) DNA. |
| `rna_baseline_bam` | no | Tumour RNA BAM at baseline, already aligned (e.g. STAR) to the same GRCh38 reference with `chr` names, with an index (`<x>.bam.bai` or `<x>.bai`) beside it. Counted with bam-readcount at the somatic sites; pVACseq gets these counts as its tumour RNA depth/VAF. |
| `rna_week6_bam` | no | Tumour RNA BAM at week 6, same requirements. Counted the same way and added to the CCF tables as `RNA Depth (week 6)` and `RNA VAF (week 6)`; pVACseq does not filter on it. |
| `purity` | yes | Sequenza cellularity, a number in (0, 1]. Passed to `pvacseq run --tumor-purity` and used for CCF. |
| `ploidy` | no | Sequenza ploidy. Used as the tumour copy number for CCF when no segments file is given, or for variants outside every segment. |
| `sequenza_segments` | no | Path to the Sequenza `<sample>_segments.txt`. Gives local copy number (CNt, A, B) per variant for CCF. |
| `sex` | no | `male` or `female`. Only affects CCF on chrX/chrY (normal copy number 1 in males). |
| `hla_alleles` | no | Skip HLA-HD for this patient and use these alleles. Separate with `\|` or `;` (a comma would split the CSV column). Use pVACseq names, e.g. `HLA-A*02:01\|HLA-B*07:02\|DRB1*15:01\|DQA1*01:02-DQB1*06:02`. |

Give each sample either a FASTQ pair or a BAM, not both. A row can mix
them (e.g. tumour FASTQ, normal BAM). One FASTQ pair per sample is
supported; merge lanes beforehand if a sample was sequenced on several.

Purity and ploidy are kept out of the cached task identity on purpose:
correcting a purity value and rerunning with `-resume` reruns only
pVACseq and the CCF step, not alignment or calling.

RNA BAMs are optional and independent: a patient can have baseline only,
week 6 only, both or neither. Only read counts are taken from them; there
is no expression quantification, so pVACseq's gene/transcript expression
columns stay NA.

```csv
patient_id,tumor_r1,tumor_r2,normal_r1,normal_r2,tumor_bam,normal_bam,rna_baseline_bam,rna_week6_bam,purity,ploidy,sequenza_segments,sex,hla_alleles
NeoTrio_003,/data/003_T_R1.fq.gz,/data/003_T_R2.fq.gz,/data/003_N_R1.fq.gz,/data/003_N_R2.fq.gz,,,/data/rna/003_baseline.bam,/data/rna/003_wk6.bam,0.62,3.1,/data/sequenza/003/003_segments.txt,female,
NeoTrio_004,,,,,/data/004_T.bam,/data/004_N.bam,,,0.41,2.0,,male,HLA-A*01:01|HLA-A*03:01|HLA-B*08:01|HLA-B*07:02|HLA-C*07:01|HLA-C*07:02
```
