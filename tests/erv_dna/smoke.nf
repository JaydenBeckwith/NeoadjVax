// Run from repository root, using tests/erv_dna/smoke.config.
nextflow.enable.dsl = 2
include { ERV_DNA; validateErvDnaSamplesheet } from '../../workflows/erv_dna'

params.dna_samplesheet = 'tests/erv_dna/samples.csv'
params.run_dna_variant_calling = false

workflow {
    validateErvDnaSamplesheet(params.dna_samplesheet, false)
    rows = Channel.fromPath(params.dna_samplesheet).splitCsv(header: true)
    ERV_DNA(rows, Channel.empty(), Channel.empty(), false)
}
