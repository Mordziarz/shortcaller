process CDS_UTR_ANNOTATION {
    tag "Running TransDecoder on StringTie merged GTF"
    publishDir "${params.outdir}/transdecoder_cds_utr_predict", mode: 'copy'
    container 'mordziarz/pipeline_rnaseq:latest'

    input:
    path gtf
    path fasta

    output:
    path "genome.transdecoder.gff3", emit: gff3
    path "transcriptome.fa.transdecoder.cds", emit: cds
    path "transcriptome.fa.transdecoder.pep", emit: pep
    path "*", emit: all

    script:
    """
    gffread -w transcriptome.fa -g ${fasta} ${gtf}

    if [ -f "/opt/conda/opt/transdecoder/util/gtf_to_alignment_gff3.pl" ]; then
        CONV_SCRIPT="/opt/conda/opt/transdecoder/util/gtf_to_alignment_gff3.pl"
        MAP_SCRIPT="/opt/conda/opt/transdecoder/util/cdna_alignment_orf_to_genome_orf.pl"
    else
        CONV_SCRIPT="gtf_to_alignment_gff3.pl"
        MAP_SCRIPT="cdna_alignment_orf_to_genome_orf.pl"
    fi

    \$CONV_SCRIPT ${gtf} > transcripts.gff3

    /opt/conda/opt/transdecoder/util/TransDecoder.LongOrfs -t transcriptome.fa
    /opt/conda/opt/transdecoder/util/TransDecoder.Predict -t transcriptome.fa --single_best_only

    \$MAP_SCRIPT \
        transcriptome.fa.transdecoder.gff3 \
        transcripts.gff3 \
        transcriptome.fa > genome.transdecoder.gff3
    """
}
