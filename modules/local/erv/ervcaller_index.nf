// Isolated, cached reference preparation. Never write beside the source FASTAs.
process ERVCALLER_INDEX {
    tag 'genome-and-TE-library'
    label 'process_high'
    container params.ervcaller_container

    input:
    path 'source-genome.fa'
    path 'source-te.fa'

    output:
    path 'ervcaller_refs', emit: references

    script:
    """
    set -euo pipefail
    mkdir ervcaller_refs
    cp -L source-genome.fa ervcaller_refs/genome.fa
    cp -L source-te.fa ervcaller_refs/te.fa
    bwa index ervcaller_refs/genome.fa
    bwa index ervcaller_refs/te.fa
    samtools faidx ervcaller_refs/genome.fa
    samtools faidx ervcaller_refs/te.fa
    sha256sum ervcaller_refs/genome.fa ervcaller_refs/te.fa > ervcaller_refs/references.sha256
    """

    stub:
    """
    mkdir ervcaller_refs
    for ref in genome te; do
        for suffix in fa fa.amb fa.ann fa.bwt fa.pac fa.sa fa.fai; do
            touch ervcaller_refs/\${ref}.\${suffix}
        done
    done
    """
}
