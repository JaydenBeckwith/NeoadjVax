process SPLICING_STAR_INDEX {
    label 'process_high'
    container params.containers.star
    input:
    path fasta
    path gtf
    output:
    path 'star_index', emit: index
    script:
    """
    mkdir star_index
    STAR --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles '${fasta}' \\
        --sjdbGTFfile '${gtf}' --sjdbOverhang ${params.splicing_star_overhang} --runThreadN ${task.cpus}
    """
    stub:
    """
    mkdir star_index
    touch star_index/SA
    """
}

process SPLICING_STAR_ALIGN {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_high'
    container params.containers.star
    publishDir "${params.outdir}/${meta.id}/splicing/${meta.timepoint}/star", mode: 'copy', pattern: '*.out'
    input:
    tuple val(meta), path(r1), path(r2)
    path index
    output:
    tuple val(meta), path('rna.Aligned.sortedByCoord.out.bam'), emit: bam
    path 'rna.Log.final.out', emit: log
    script:
    def decompress = r1.name.endsWith('.gz') ? '--readFilesCommand zcat' : ''
    """
    STAR --runThreadN ${task.cpus} --genomeDir '${index}' --readFilesIn '${r1}' '${r2}' ${decompress} \\
        --twopassMode Basic --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes NH HI AS NM MD --outFilterMultimapNmax 1 \\
        --limitBAMsortRAM 40000000000 --outFileNamePrefix rna.
    """
    stub:
    """
    touch rna.Aligned.sortedByCoord.out.bam rna.Log.final.out
    """
}

process LEAFCUTTER {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_medium'
    container params.splicing_tools_container
    publishDir "${params.outdir}/${meta.id}/splicing/${meta.timepoint}/leafcutter", mode: 'copy'
    input:
    tuple val(meta), path(bam)
    output:
    tuple val(meta), path('leafcutter/sample_perind.counts.gz'), emit: counts
    path 'sample.junc', emit: junctions
    script:
    """
    samtools quickcheck -v '${bam}'
    samtools view -@ ${task.cpus} -b -F 2308 -e '[NH] == 1' -o unique.bam '${bam}'
    samtools index -@ ${task.cpus} unique.bam
    regtools junctions extract -a 8 -m 50 -M 500000 -s '${meta.strand}' -o sample.junc unique.bam
    mkdir leafcutter
    if [ -s sample.junc ]; then
        printf '%s\\n' sample.junc > junctions.txt
        python /opt/leafcutter/leafcutter_cluster_regtools.py -j junctions.txt -o sample -r leafcutter \\
            -m ${params.leafcutter_min_reads} -p ${params.leafcutter_min_ratio} -l 500000
    else
        printf 'chrom sample\\n' | gzip -c > leafcutter/sample_perind.counts.gz
    fi
    test -s leafcutter/sample_perind.counts.gz
    """
    stub:
    """
    mkdir leafcutter
    printf 'chrom sample\\n' | gzip -c > leafcutter/sample_perind.counts.gz
    touch sample.junc
    """
}

process LEAFCUTTER_SELECT {
    tag "${meta.id}:${meta.timepoint}"
    label 'process_low'
    container params.containers.splicing_python
    publishDir "${params.outdir}/${meta.id}/splicing/${meta.timepoint}/leafcutter", mode: 'copy'
    input:
    tuple val(meta), path(counts), val(sample)
    output:
    tuple val(meta), path('selected_perind.counts.gz'), path('rna_junctions.tsv'), emit: selected
    script:
    """
    python '${projectDir}/bin/prepare_splicing_inputs.py' leafcutter --input '${counts}' --sample '${sample}' \\
        --output selected_perind.counts.gz --evidence rna_junctions.tsv --min-reads ${params.splicing_min_rna_reads}
    """
    stub:
    """
    touch selected_perind.counts.gz rna_junctions.tsv
    """
}
