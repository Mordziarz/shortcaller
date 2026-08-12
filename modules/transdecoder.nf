process TRANSDECODER {
    tag "TransDecoder ORFs"
    container 'mordziarz/pipeline_rnaseq:latest'
    publishDir "${params.outdir}/functional_annotation/transdecoder", mode: 'copy'

    input:
    path transcripts_fa

    output:
    path "target_transcripts.fa.transdecoder_dir*"

    script:
    """
    perl /opt/conda/opt/transdecoder/util/TransDecoder.LongOrfs -t ${transcripts_fa}
    """
}
