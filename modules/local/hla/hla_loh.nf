// HLA loss of heterozygosity, read in rather than computed.
//
// SpecHLA (typing + cal.hla.copy.pl, with Sequenza purity and ploidy) is
// run outside this pipeline, so all that is needed here is to normalise
// its merge.hla.copy.txt into a per-locus table with the alleles spelled
// the way pVACseq spells them (bin/parse_hla_loh.py). The prioritisation
// step then flags any candidate whose restricting allele is gone.
//
// This matters more than any ranking change: a candidate predicted on an
// allele the tumour has deleted cannot be presented at all, and nothing
// else in the candidate table would reveal it.

process PARSE_HLA_LOH {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.pvactools
    publishDir "${params.outdir}/${meta.id}/hla_loh", mode: 'copy'

    input:
    tuple val(meta), path(loh)

    output:
    tuple val(meta), path("${meta.id}.hla_loh.tsv"), emit: tsv

    script:
    def definition  = params.hla_loh_definition  ? "--definition ${params.hla_loh_definition}" : ''
    def min_het     = params.hla_loh_min_het != null ? "--min-het ${params.hla_loh_min_het}" : ''
    def min_purity  = params.hla_loh_min_purity != null ? "--min-purity ${params.hla_loh_min_purity}" : ''
    """
    parse_hla_loh.py ${loh} \\
        --patient ${meta.id} \\
        ${min_het} \\
        ${min_purity} \\
        ${definition} \\
        -o ${meta.id}.hla_loh.tsv
    """

    stub:
    """
    printf 'patient\\tlocus\\thla_class\\tcall\\tlost_allele_pvacseq\\n' > ${meta.id}.hla_loh.tsv
    printf '${meta.id}\\tA\\tI\\tretained\\t\\n' >> ${meta.id}.hla_loh.tsv
    """
}
