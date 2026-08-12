process STAR_INDEX {
    publishDir "${params.outdir}/star_index", mode: 'copy'

    container 'mordziarz/pipeline_rnaseq:latest'

    input:
    path fasta
    path gtf

    output:
    path "star_index_dir", emit: index

    script:
    """
    mkdir -p star_index_dir
    STAR --runMode genomeGenerate \
        --runThreadN ${params.threads_star_index} \
        --genomeDir star_index_dir \
        --genomeFastaFiles ${fasta} \
        --genomeSAindexNbases 12 \
        --sjdbGTFfile ${gtf}
    """
}
