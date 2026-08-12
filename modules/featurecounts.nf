process FEATURECOUNTS {
    publishDir "${params.outdir}/featurecounts", mode: 'copy'
    container 'mordziarz/pipeline_rnaseq:latest'

    input:
    path bams        
    path merged_gtf  

    output:
    path "Feature_counts.txt", emit: counts
    path "Feature_counts.txt.summary", emit: summary

    script:
    """
    featureCounts -T ${params.threads_featurecounts} \
        -p \
        -a ${merged_gtf} \
        -o Feature_counts.txt \
        *.bam
    """
}
