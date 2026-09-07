// Ported from: DNA pVACseq script -> bwa_index_if_missing()
// Builds a BWA index once, reused across all tumor/normal samples.
// NOTE: on Gadi, prefer pre-building this out-of-band once and pointing
// params.genome_fasta at an already-indexed FASTA: this process exists so
// the pipeline is self-contained for smaller/local runs too.

process BWA_INDEX {
    tag "${fasta.simpleName}"
    label 'process_medium'
    container params.containers.bwa

    input:
    path fasta

    output:
    path "${fasta}.*", emit: index
    path fasta,         emit: fasta

    script:
    """
    if [ -f "${fasta}.bwt" ]; then
        echo "[INFO] BWA index already present next to ${fasta}"
    else
        bwa index ${fasta}
    fi
    """
}
