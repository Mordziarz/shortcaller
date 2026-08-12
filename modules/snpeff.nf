process SNPEFF {
    tag "Annotating variants with SnpEff"
    container 'mordziarz/pipeline_rnaseq:latest'
    publishDir "${params.outdir}/snpeff_annotation", mode: 'copy'

    input:
    path vcf
    path gff
    path fasta

    output:
    path "annotated_variants.vcf.gz", emit: vcf
    path "annotated_variants.vcf.gz.tbi", emit: index
    path "snpEff_summary.html", emit: html

    script:
    """
    mkdir -p data/my_organism
    
    cp ${fasta} data/my_organism/sequences.fa
    cp ${gff} data/my_organism/genes.gff
    
    echo "my_organism.genome : my_organism" >> snpEff.config
    
    snpEff build -gff3 -v -noCheckCds -noCheckProtein my_organism
    
    snpEff -v my_organism ${vcf} > temp_annotated.vcf

    bgzip -c temp_annotated.vcf > annotated_variants.vcf.gz
    tabix -p vcf annotated_variants.vcf.gz
    rm temp_annotated.vcf
    """
}
