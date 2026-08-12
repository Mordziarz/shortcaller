process MULTIQC {
    tag "multiqc"
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    container 'mordziarz/pipeline_rnaseq:latest'

    input:
    path fastqc_dir

    output:
    path "multiqc_report.html", emit: html
    path "multiqc_data",        emit: data

    script:
    """
    multiqc ${fastqc_dir}
    """
}
