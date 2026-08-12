process FASTP {
    tag "$sample"
    publishDir "${params.outdir}/fastp", mode: 'copy'

    container 'mordziarz/pipeline_rnaseq:latest'

    input:
    tuple val(sample), val(group), path(reads)

    output:
    tuple val(sample), val(group), path("${sample}_1.trim.fastq.gz"), path("${sample}_2.trim.fastq.gz"), emit: reads
    path "${sample}.fastp.json", emit: json
    path "${sample}.fastp.html", emit: html

    script:
    """
    fastp -w ${task.cpus} \\
        -i ${reads[0]} -I ${reads[1]} \\
        -o ${sample}_1.trim.fastq.gz -O ${sample}_2.trim.fastq.gz \\
        --cut_front --cut_tail \\
        --cut_mean_quality ${params.trim_avgqual} \\
        --length_required ${params.trim_minlen} \\
        --json ${sample}.fastp.json \\
        --html ${sample}.fastp.html
    """
}
