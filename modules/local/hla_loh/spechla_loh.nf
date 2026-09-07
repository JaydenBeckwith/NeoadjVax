// Ported from: NeoadjLOH/specHLA/run_sample.pbs, Step 3 ("HLA LOH")
// https://github.com/JaydenBeckwith/NeoadjLOH
//
// Combines the typing result + freq list from SPECHLA_TYPING with tumor
// purity/ploidy to compute per-locus allelic copy ratio and call LOH.
// Purity/ploidy come from Sequenza (see workflows/purity_ploidy.nf: a
// stub pending Jayden's Sequenza_tools code) or, in the meantime, from a
// manually-supplied CSV / standalone samplesheet: see
// assets/samplesheet_schema.md and docs/ARCHITECTURE.md.
//
// -C 5 is a fixed flag in the original run_sample.pbs (not templated
// there), kept hardcoded here for the same reason.

process SPECHLA_LOH {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.spechla
    publishDir "${params.outdir}/${meta.id}/hla_loh/loh", mode: 'copy'

    input:
    tuple val(meta), path(typing_out), path(hla_result), path(freq_list), val(purity), val(ploidy)

    output:
    tuple val(meta), path("loh_out/merge.hla.copy.txt"), emit: loh_calls
    tuple val(meta), path("loh_out"), emit: loh_dir

    script:
    """
    mkdir -p loh_out
    ${params.spechla_loh_cmd} -S ${meta.id} -C 5 -purity ${purity} -ploidy ${ploidy} \\
        -F ${freq_list} -T ${hla_result} -O loh_out

    if [[ ! -s loh_out/merge.hla.copy.txt ]]; then
        echo "ERROR: expected loh_out/merge.hla.copy.txt was not produced for ${meta.id}" >&2
        exit 1
    fi
    echo "[INFO] SUCCESS: loh_out/merge.hla.copy.txt (purity=${purity}, ploidy=${ploidy})"
    """
}
