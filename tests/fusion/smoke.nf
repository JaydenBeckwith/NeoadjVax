// Run from repository root, using tests/fusion/smoke.config.
nextflow.enable.dsl = 2
include { FUSION_NEOANTIGENS; validateFusionInputs } from '../../workflows/fusion_neoantigens'

params.rna_samplesheet = 'tests/fusion/samples.csv'

workflow {
    def settings = validateFusionInputs(params.rna_samplesheet, '', false, params)
    rows_ch = Channel.fromPath(params.rna_samplesheet).splitCsv(header: true)
    FUSION_NEOANTIGENS(rows_ch, Channel.empty(), settings)
}
