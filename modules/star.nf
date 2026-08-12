process STAR {
    tag "$sample"
    publishDir "${params.outdir}/star", mode: 'copy'

    container 'mordziarz/pipeline_rnaseq:latest'

    input:
    tuple val(sample), val(group), path(reads)
    path genome_dir

    output:
    tuple val(sample), val(group), path("${sample}_Aligned.sortedByCoord.out.bam"), emit: bam
    path "${sample}_Log.final.out"                                                 , emit: log_final
    path "${sample}_Log.out"                                                       , emit: log_out
    path "${sample}_Log.progress.out"                                              , emit: log_progress

    script:
    """
    STAR --runThreadN ${params.threads_star} \
        --genomeDir ${genome_dir} \
        --readFilesIn ${reads[0]} ${reads[1]} \
        --readFilesCommand gunzip -c \
        --outFileNamePrefix ${sample}_ \
        --outSAMmapqUnique 50 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs RemoveNoncanonical \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.04 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
    """
}
