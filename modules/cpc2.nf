process CPC2 {
    tag "CPC2 coding potential"
    container 'mordziarz/pipeline_rnaseq:latest'
    publishDir "${params.outdir}/functional_annotation/cpc2", mode: 'copy'

    input:
    path transcripts_fa

    output:
    path "cpc2_results.txt"

    script:
    """
    # Uruchomienie poprawnej, skompilowanej wersji CPC2 ze ścieżki /opt/cpc2
    /opt/cpc2/bin/CPC2.py -i ${transcripts_fa} -o cpc2_results
    """
}
