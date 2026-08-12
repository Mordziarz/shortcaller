process STRINGTIE_RUN {
    tag "$sample"
    publishDir "${params.outdir}/stringtie", mode: 'copy'
    container 'mordziarz/pipeline_rnaseq:latest'

    input:
    tuple val(sample), val(group), path(bam), path(gtf), val(strand_flag)

    output:
    tuple val(sample), path("${sample}.gtf"), emit: gtf
    tuple val(sample), path("${sample}.gtf"), emit: gtf_path

    script:
    """
    stringtie ${bam} -p ${task.cpus} --${strand_flag} -G ${gtf} -o ${sample}.gtf -l ${sample}
    """
}

process STRINGTIE_MERGE {
    publishDir "${params.outdir}/stringtie", mode: 'copy'
    container 'mordziarz/pipeline_rnaseq:latest'

    input:
    path gtf_files
    path gtf
    val strand_flag

    output:
    path "stringtie_merged.gtf", emit: merged_gtf

    script:
    """
    ls *.gtf > mergelist.txt
    stringtie --merge -p ${task.cpus} --${strand_flag} -G ${gtf} -o stringtie_merged.gtf -l CL mergelist.txt
    """
}
