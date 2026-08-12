process GFFREAD {
    tag "Extracting fasta via gffread"
    container 'mordziarz/pipeline_rnaseq:latest'
    publishDir "${params.outdir}/functional_annotation", mode: 'copy'

    input:
    path target_ids
    path gtf
    path fasta

    output:
    path "target_transcripts.fa", emit: fasta

    script:
    """
    awk '{print "gene_id \\"" \$1 "\\""}' ${target_ids} > patterns.txt

    grep -F -f patterns.txt ${gtf} > filtered_annotation.gtf || true

    if [ ! -s filtered_annotation.gtf ]; then
        echo "Warning: No matching entries in GTF. Using the full file."
        cp ${gtf} filtered_annotation.gtf
    fi

    gffread -w target_transcripts.fa -g ${fasta} filtered_annotation.gtf

    echo "Number of extracted transcript sequences:"
    grep -c ">" target_transcripts.fa || echo "0"
    """
}
