// Ported from: RNA_variant_pipeline.sh, Stage 1 (star_index_build PBS job)
// sjdbOverhang 149 was hardcoded for the original read length — exposed as
// a param since future cohorts may use a different read length.

process STAR_INDEX {
    tag "${fasta.simpleName}"
    label 'process_high'
    container params.containers.star

    input:
    path fasta
    path gtf

    output:
    path "star_index", emit: index

    script:
    def overhang = params.star_sjdb_overhang ?: 149
    """
    mkdir -p star_index
    if [ -f "${params.star_index_dir}/SA" ] && [ -s "${params.star_index_dir}/SA" ]; then
        echo "[INFO] STAR index already exists at ${params.star_index_dir} — linking instead of rebuilding"
        ln -s ${params.star_index_dir}/* star_index/
    else
        STAR \\
            --runMode genomeGenerate \\
            --genomeDir star_index \\
            --genomeFastaFiles ${fasta} \\
            --sjdbGTFfile ${gtf} \\
            --sjdbOverhang ${overhang} \\
            --runThreadN ${task.cpus}
    fi
    """
}
