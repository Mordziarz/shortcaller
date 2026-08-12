process BWA_MEM {
    tag "$sample"
    container 'mordziarz/pipeline_rnaseq:latest'
    publishDir "${params.outdir}/bwa_mem", mode: 'copy'

    input:
    tuple val(sample), val(group), path(reads)
    path fasta
    path bwa_idx

    output:
    tuple val(sample), val(group), path("${sample}.sam"), emit: sam

    script:
    """
    bwa mem -T 16 ${bwa_idx}/${fasta.name} ${reads[0]} ${reads[1]} > ${sample}.sam
    """
}
