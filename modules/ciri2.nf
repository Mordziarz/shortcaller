process CIRI2 {
    tag "$sample"
    publishDir "${params.outdir}/ciri2", mode: 'copy'
    container 'mordziarz/pipeline_rnaseq:latest'

    input:
    tuple val(sample), val(group), path(sam), path(fasta)
    path index
    path merged_gtf

    output:
    path "${sample}_circRNA_results.txt", emit: circrna

    script:
    """
    perl /opt/conda/bin/CIRI2.pl \
        -I ${sam} \
        -O ${sample}_circRNA_results.txt \
        -F ${fasta} \
        -A ${merged_gtf} \
        -T 4
    """
}
