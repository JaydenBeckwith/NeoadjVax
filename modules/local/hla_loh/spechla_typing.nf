// Ported from: NeoadjLOH/specHLA/run_sample.pbs, Step 2 ("SpecHLA typing")
// https://github.com/JaydenBeckwith/NeoadjLOH
//
// Calls 4-field HLA alleles (A/B/C/DPA1/DPB1/DQA1/DQB1/DRB1) from the
// extracted read pair. Also selects the typing result file and builds the
// freq-list Step 3 (cal.hla.copy.pl) needs — both moved into THIS process
// rather than staying separate shell steps like the original PBS script,
// which is the one structural change from the original: run_sample.pbs
// could freely `find` across the whole job's $OUTROOT because every step
// shared one working directory, but each Nextflow process gets its own
// isolated work dir, so the result-selection and freq-list logic has to
// run right after SpecHLA.sh produces its output, in the same task, while
// still landing in this same directory. The selection logic itself, and
// its gotcha comments, are unchanged from the original.
//
// The whole "typing_out" directory is emitted (not just the two files
// picked out below) because SPECHLA_LOH's -F freq list has to reference
// the underlying *_freq.txt files by path, and re-declaring each one
// individually as a process output isn't practical when the exact set of
// freq files SpecHLA writes isn't fixed ahead of time.

process SPECHLA_TYPING {
    tag "${meta.id}"
    label 'process_high'
    container params.containers.spechla
    publishDir "${params.outdir}/${meta.id}/hla_loh/typing", mode: 'copy'

    input:
    tuple val(meta), path(fq1), path(fq2)

    output:
    tuple val(meta), path("typing_out"), path("hla_result_selected.txt"), path("freq.list"), emit: typing

    script:
    """
    mkdir -p typing_out
    ${params.spechla_typing_cmd} -n ${meta.id} -1 ${fq1} -2 ${fq2} -o typing_out

    # Prefer the exact plain filename cal.hla.copy.pl is documented to expect
    # (SpecHLA's own usage example is "-T hla.result.txt"). Do NOT just glob
    # "*hla*result*.txt" and take the alphabetically-first hit: SpecHLA also
    # writes a more verbose "hla.result.details.txt", which sorts before the
    # plain file ("." < "t" in "details" vs "txt") and silently gets fed to
    # cal.hla.copy.pl instead — it then fails to parse it (visible as repeated
    # "Use of uninitialized value \$hla1/\$hla2" warnings) and leaves
    # Allele1/Allele2/KeptHLA/LossHLA blank in the output while copyratio
    # still looks plausible, so the failure is easy to miss.
    HLA_RESULT=\$(find typing_out -iname "hla.result.txt" | sort | head -n1)
    if [[ -z "\$HLA_RESULT" ]]; then
        HLA_RESULT=\$(find typing_out -iname "*hla*result*.txt" ! -iname "*details*" | sort | head -n1)
    fi
    if [[ -z "\$HLA_RESULT" ]]; then
        HLA_RESULT=\$(find typing_out -iname "*hla*result*.txt" | sort | head -n1)
        if [[ -n "\$HLA_RESULT" ]]; then
            echo "WARNING: only found '\$HLA_RESULT' (no plain hla.result.txt) for ${meta.id} — this previously caused cal.hla.copy.pl to leave KeptHLA/LossHLA blank; check the typing output naming" >&2
        fi
    fi
    if [[ -z "\$HLA_RESULT" ]]; then
        echo "ERROR: no hla result file found under typing_out for ${meta.id} — inspect the typing output and fix the 'find' pattern in this module if the naming differs" >&2
        exit 1
    fi
    echo "[INFO] Typing result: \$HLA_RESULT"
    cp "\$HLA_RESULT" hla_result_selected.txt

    # freq list entries are written relative to this process's work dir
    # ("typing_out/<file>"), not absolute paths, so the same list still
    # resolves correctly once "typing_out" is re-staged into SPECHLA_LOH's
    # own (different) work dir.
    find typing_out -iname "*_freq.txt" | sort > freq.list
    if [[ ! -s freq.list ]]; then
        echo "ERROR: no *_freq.txt files found under typing_out for ${meta.id} for the -F freq list" >&2
        exit 1
    fi
    echo "[INFO] Freq list: \$(wc -l < freq.list) files"
    """
}
