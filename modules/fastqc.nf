process FASTQC {
    tag "$sample"
    publishDir "${params.outdir}/fastqc_results", mode: 'copy'

    container 'mordziarz/pipeline_rnaseq:latest'

    input:
    tuple val(sample), val(group), path(reads)

    output:
    path "*.{html,zip}", emit: html_zip

    script:
    """
    fastqc -o . -t ${params.threads_fastqc} ${reads[0]} ${reads[1]}
    """
}
