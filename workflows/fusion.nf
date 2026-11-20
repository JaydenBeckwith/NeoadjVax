/*
 * Optional gene-fusion neoantigens (--fusion true), per patient:
 *
 *   baseline tumour RNA FASTQ ─> STAR-Fusion (CTAT genome lib)
 *     or the samplesheet's starfusion_tsv (an existing STAR-Fusion run)
 *   ─> AGFusion (Ensembl DB, --middlestar --noncanonical)
 *   ─> pVACfuse (same HLA, image and IEDB as pVACseq; --starfusion-file)
 *
 * Mirrors the Gadi PBS scripts run_agfusion.pbs and pvacfuse.pbs. Patients
 * with neither a FASTQ pair nor a STAR-Fusion table are skipped. See
 * docs/FUSION.md.
 */

include { STAR_FUSION } from '../modules/local/fusion/star_fusion'
include { AGFUSION    } from '../modules/local/fusion/agfusion'
include { PVACFUSE    } from '../modules/local/pvactools/pvacfuse'

workflow FUSION {

    take:
    fastq_ch          // [meta, r1, r2]           patients that need STAR-Fusion
    calls_ch          // [meta, predictions.tsv]  patients with an existing STAR-Fusion table
    hla_ch            // [meta, allele string]
    ctat_genome_lib   // CTAT ctat_genome_lib_build_dir, or [] when nobody needs STAR-Fusion
    agfusion_db

    main:
    STAR_FUSION(fastq_ch, ctat_genome_lib)
    predictions_ch = STAR_FUSION.out.predictions.mix(calls_ch)

    AGFUSION(predictions_ch, agfusion_db)

    PVACFUSE(AGFUSION.out.annotated
        .join(predictions_ch)
        .join(hla_ch))

    emit:
    predictions = predictions_ch       // [meta, full STAR-Fusion table]
    agfusion    = AGFUSION.out.annotated
    pvacfuse    = PVACFUSE.out.results // [meta, pvacfuse/]
}
