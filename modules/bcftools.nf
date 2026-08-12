process MERGE_VCFS {
    publishDir "${params.outdir}/variants_merged", mode: 'copy'
    container 'mordziarz/pipeline_rnaseq:latest'

    input:
    path vcfs 

    output:
    path "merged_variants.vcf.gz", emit: merged_vcf
    path "merged_variants.vcf.gz.tbi", emit: merged_index

    script:
    """
    for file in *.vcf; do
        SAMPLENAME=\$(basename \$file _variants.vcf)
        echo "\$SAMPLENAME" > \${SAMPLENAME}_name.txt
        bcftools view -O z -o \${file}.gz \$file
        bcftools reheader -s \${SAMPLENAME}_name.txt -o \${SAMPLENAME}_ready.vcf.gz \${file}.gz
        bcftools index -t \${SAMPLENAME}_ready.vcf.gz
    done

    bcftools merge -m all -0 -O z -o temp_merged.vcf.gz *_ready.vcf.gz
    bcftools norm -m -any -O z -o merged_variants.vcf.gz temp_merged.vcf.gz
    bcftools index -t merged_variants.vcf.gz
    """
}

process VCF_TO_MATRIX {
    publishDir "${params.outdir}/variants_matrix", mode: 'copy'
    container "staphb/bcftools:1.19"

    input:
    path annotated_vcf
    path annotated_index

    output:
    path "variants_matrix.tsv", emit: matrix

    script:
    """
    echo -n -e "CHROM\\tPOS\\tREF\\tALT\\tANN" > header.txt
    
    bcftools query -l ${annotated_vcf} | while read sample; do
        echo -n -e "\\t\${sample}_GT\\t\${sample}_DP\\t\${sample}_RO\\t\${sample}_AO" >> header.txt
    done
    echo "" >> header.txt

    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/ANN[\\t%GT\\t%DP\\t%RO\\t%AO]\\n' ${annotated_vcf} > data.txt

    cat header.txt data.txt > variants_matrix.tsv
    """
}
