// pVACfuse on AGFusion's fusion proteins, as pvacfuse.pbs runs it: same
// pvactools image and IEDB install as pVACseq, --starfusion-file with the
// full STAR-Fusion table so each candidate carries its junction plus
// spanning read support and FFPM (pVACfuse tiers on both: Read Support < 5
// and Expr < 0.1 are flagged).
//
// AGFusion writes no *_protein.fa when none of a patient's fusions has a
// coding partner pair; pVACfuse would then stop before writing its class I
// and II reports and fail building the combined one, so that case is
// written out as a note and skipped.

process PVACFUSE {
    tag "${meta.id}"
    label 'process_pvacfuse'
    container params.containers.pvactools
    publishDir "${params.outdir}/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(agfusion_dir), path(predictions), val(hla_alleles)

    output:
    tuple val(meta), path("pvacfuse"), emit: results

    script:
    def iedb = params.pvacseq_iedb_install_dir ? "--iedb-install-directory ${params.pvacseq_iedb_install_dir}" : ''
    def problematic = params.pvacfuse_problematic_amino_acids ? "--problematic-amino-acids ${params.pvacfuse_problematic_amino_acids}" : ''
    """
    mkdir -p pvacfuse tmp
    if ! ls ${agfusion_dir}/*/*_protein.fa > /dev/null 2>&1; then
        echo "AGFusion produced no fusion protein sequences for ${meta.id}: nothing for pVACfuse to predict" > pvacfuse/no_fusion_proteins.txt
        exit 0
    fi

    # keep pVACtools' temp files off the node's /tmp
    export TMPDIR="\$PWD/tmp"
    pvacfuse run \\
        ${agfusion_dir} \\
        ${meta.id} \\
        '${hla_alleles}' \\
        ${params.pvacfuse_algorithms} \\
        pvacfuse \\
        -e1 ${params.pvacfuse_class_i_lengths} \\
        -e2 ${params.pvacfuse_class_ii_lengths} \\
        ${problematic} \\
        --starfusion-file ${predictions} \\
        ${iedb} \\
        -t ${task.cpus} \\
        ${params.pvacfuse_extra_args}
    rm -rf tmp
    """

    stub:
    def header = 'ID\\tGene\\tBest Peptide\\tBest Transcript\\tAllele\\tProb Pos\\tNum Included Peptides\\tNum Passing Peptides\\tIC50 MT\\t%ile MT\\tExpr\\tRead Support\\tTier\\tEvaluation'
    """
    mkdir -p pvacfuse/MHC_Class_I pvacfuse/MHC_Class_II pvacfuse/combined
    printf 'Stub run: no biological analysis performed.\\n' > pvacfuse/stub.log
    for c in I II; do
        printf '${header}\\n' > pvacfuse/MHC_Class_\$c/${meta.id}.MHC_\$c.all_epitopes.aggregated.tsv
    done
    printf '${header}\\n' > pvacfuse/combined/${meta.id}.Combined.all_epitopes.aggregated.tsv
    """
}
