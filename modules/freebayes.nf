process FREEBAYES {
    tag "$meta.id"
    container 'mordziarz/pipeline_rnaseq:latest'
    publishDir "${params.outdir}/freebayes", mode: 'copy'

    input:
    tuple val(meta), path(bam)
    path fasta

    output:
    tuple val(meta), path("*_variants.vcf"), emit: vcf

    script:
    """
    if [ ! -f "${fasta}.fai" ]; then
        samtools faidx ${fasta}
    fi

    samtools index ${bam}
    
    freebayes -f ${fasta} \
        --min-mapping-quality 20 \
        --min-base-quality 20 \
        --min-coverage 10 \
        --min-alternate-count 2 \
        ${bam} > ${meta.id}_variants.vcf
    """
}
