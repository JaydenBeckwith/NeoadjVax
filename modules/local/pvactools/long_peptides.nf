// The vaccine deliverable: synthetic long peptides, not epitopes.
//
// pVACseq ranks 8-11mers (class I) and 12-18mers (class II). A peptide
// vaccine is made of 20-30mers that carry the mutation with its natural
// flanking sequence, so the patient's own antigen-presenting cells
// process and present it on whichever allele fits rather than the
// construct betting on one predicted epitope. Turning the candidate
// table into that peptide list was the step still being done by hand.
//
// Two processes, because the flanking sequence is not in any table:
// pVACtools regenerates the mutant protein context, then
// bin/design_long_peptides.py places, trims and checks the window.

process GENERATE_PROTEIN_FASTA {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.pvactools
    publishDir "${params.outdir}/${meta.id}/long_peptides", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi), path(phased_vcf), path(phased_tbi)

    output:
    tuple val(meta), path("${meta.id}.protein.fasta"), emit: fasta

    script:
    // The flanking length is per side, so max_length is more context than
    // any designed peptide needs - deliberately, since trimming can only
    // move the window inside sequence that exists.
    def phased = phased_vcf ? "--phased-proximal-variants-vcf ${phased_vcf}" : ''
    """
    pvacseq generate_protein_fasta \\
        ${vcf} \\
        ${params.long_peptide_max_length} \\
        ${meta.id}.protein.fasta \\
        --sample-name tumor \\
        --pass-only \\
        ${phased}
    """

    stub:
    """
    printf '>MT.1.GENE1.ENST001.missense\\n' > ${meta.id}.protein.fasta
    printf 'AAAAAGGGGGSSSSSTTTTTVVVVVWLLLLKKKKKDDDDD\\n' >> ${meta.id}.protein.fasta
    printf '>WT.1.GENE1.ENST001.missense\\n' >> ${meta.id}.protein.fasta
    printf 'AAAAAGGGGGSSSSSTTTTTVVVVVLLLLLKKKKKDDDDD\\n' >> ${meta.id}.protein.fasta
    """
}

process DESIGN_LONG_PEPTIDES {
    tag "${meta.id}"
    label 'process_low'
    container params.containers.pvactools
    publishDir "${params.outdir}/${meta.id}/long_peptides", mode: 'copy'

    input:
    tuple val(meta), path(reports), path(protein_fasta)

    output:
    tuple val(meta), path("*.long_peptides.tsv"), emit: tsv
    tuple val(meta), path("*.long_peptides.fasta"), emit: fasta

    script:
    def top = params.long_peptide_top ? "--top ${params.long_peptide_top}" : ''
    """
    for report in ${reports}; do
        # tumor.MHC_Class_I.<...>.prioritized.tsv -> the class, for the label
        mhc_class=\$(basename "\$report" | sed -E 's/.*\\.(MHC_Class_[IV]+|MHC_[IV]+|Combined)\\..*/\\1/; t; s/.*/unknown/')
        design_long_peptides.py "\$report" \\
            --protein-fasta ${protein_fasta} \\
            --patient ${meta.id} \\
            --mhc-class "\$mhc_class" \\
            --max-length ${params.long_peptide_max_length} \\
            --min-length ${params.long_peptide_min_length} \\
            --min-flank ${params.long_peptide_min_flank} \\
            --bad-nterm '${params.long_peptide_bad_nterm}' \\
            --bad-cterm '${params.long_peptide_bad_cterm}' \\
            ${top} \\
            -o "${meta.id}.\${mhc_class}.long_peptides.tsv" \\
            --fasta-out "${meta.id}.\${mhc_class}.long_peptides.fasta"
    done
    """

    stub:
    """
    printf 'Peptide ID\\tPatient\\tGene\\tLong Peptide\\tLength\\tWarnings\\n' > ${meta.id}.Combined.long_peptides.tsv
    printf '${meta.id}_LP01_GENE1\\t${meta.id}\\tGENE1\\tSSSSSTTTTTVVVVVWLLLLKKKKK\\t25\\t\\n' >> ${meta.id}.Combined.long_peptides.tsv
    printf '>${meta.id}_LP01_GENE1\\nSSSSSTTTTTVVVVVWLLLLKKKKK\\n' > ${meta.id}.Combined.long_peptides.fasta
    """
}

// Optional: order the designed peptides into one string-of-beads construct
// (for a DNA/RNA vaccine or a single long synthetic insert), with spacers
// chosen so the junctions between peptides do not themselves create
// neoepitopes. Off by default: it is a search over peptide orderings and
// runs the binding predictors many times over.
process PVACVECTOR {
    tag "${meta.id}"
    label 'process_pvacseq'
    container params.containers.pvactools
    publishDir "${params.outdir}/${meta.id}/vaccine_vector", mode: 'copy'

    input:
    tuple val(meta), path(peptide_fasta), val(hla_alleles)

    output:
    tuple val(meta), path("pvacvector"), optional: true, emit: results

    script:
    def iedb = params.pvacseq_iedb_install_dir ? "--iedb-install-directory ${params.pvacseq_iedb_install_dir}" : ''
    """
    pvacvector run \\
        ${peptide_fasta} \\
        ${meta.id} \\
        '${hla_alleles}' \\
        ${params.pvacvector_algorithms} \\
        pvacvector \\
        -e1 ${params.pvacseq_class_i_lengths} \\
        ${iedb} \\
        -t ${task.cpus} \\
        ${params.pvacvector_extra_args}
    """

    stub:
    """
    mkdir -p pvacvector
    printf 'Stub run: no biological analysis performed.\\n' > pvacvector/stub.log
    """
}
