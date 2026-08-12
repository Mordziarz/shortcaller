process BWA_INDEX {
    tag "$fasta"
    container 'mordziarz/pipeline_rnaseq:latest'
    publishDir "${params.outdir}/bwa_index", mode: 'copy'

    input:
    path fasta

    output:
    path "bwa_idx"

    script:
    """
    mkdir bwa_idx
    cp ${fasta} bwa_idx/
    bwa index bwa_idx/${fasta.name}
    """
}
