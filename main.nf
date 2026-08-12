nextflow.enable.dsl=2

include { FASTP } from './modules/fastp'
include { FASTQC as FASTQC_RAW  } from './modules/fastqc'
include { FASTQC as FASTQC_TRIM } from './modules/fastqc'
include { STAR_INDEX  } from './modules/star_index'
include { BWA_INDEX } from './modules/bwa_index'
include { BWA_MEM }   from './modules/bwa_mem' 
include { STAR         } from './modules/star'
include { CIRI2 } from './modules/ciri2' 
include { STRINGTIE_RUN; STRINGTIE_MERGE } from './modules/stringtie'
include { FEATURECOUNTS } from './modules/featurecounts'
include { PREPARE_RMATS; RMATS; RMATS2SASHIMIPLOT } from './modules/rmats'
include { MULTIQC     } from './modules/multiqc'
include { GFFREAD      } from './modules/gffread'
include { CPC2         } from './modules/cpc2'
include { TRANSDECODER } from './modules/transdecoder'
include { FREEBAYES } from './modules/freebayes'
include { MERGE_VCFS; VCF_TO_MATRIX } from './modules/bcftools'
include { CDS_UTR_ANNOTATION         } from './modules/cds_utr_annotation'
include { SNPEFF } from './modules/snpeff'

workflow {
    if (!params.sample_table) {
        error "Missing required parameter: --sample_table."
    }
    if (!params.fasta) {
        error "Missing required parameter: --fasta."
    }
    if (!params.gtf) {
        error "Missing required parameter: --gtf."
    }
    
    input_ch = Channel
        .fromPath(params.sample_table, checkIfExists: true)
        | splitCsv(header: true)
        | map { row -> 
            if (!row.sample || !row.fastq_1 || !row.fastq_2 || !row.group) {
                error "Error in CSV file! Check if headers are: sample, fastq_1, fastq_2, group"
            }
            return [ row.sample, row.group, [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ] ]
        }

    fasta_ch = Channel.fromPath(params.fasta, checkIfExists: true).first()
    gtf_ch   = Channel.fromPath(params.gtf, checkIfExists: true).first()
    
    STAR_INDEX(fasta_ch, gtf_ch)

    FASTQC_RAW(input_ch)

    FASTP(input_ch)

    fastq_trim_ch = FASTP.out.reads.map { sample, group, r1, r2 ->
        return [ sample, group, [r1, r2] ]
    }

    FASTQC_TRIM(fastq_trim_ch)

    STAR(fastq_trim_ch, STAR_INDEX.out.index.collect())

    stringtie_ready_ch = STAR.out.bam.combine(gtf_ch).combine(Channel.value(params.strand))
    
    STRINGTIE_RUN(stringtie_ready_ch)
    
    STRINGTIE_MERGE(STRINGTIE_RUN.out.gtf.map { sample, gtf -> gtf }.collect(), gtf_ch, params.strand)
    
    BWA_INDEX(fasta_ch)

    BWA_MEM(
        fastq_trim_ch, 
        fasta_ch, 
        BWA_INDEX.out.collect()
    )
    
    ciri_ready_ch = BWA_MEM.out.sam.combine(fasta_ch)
    
    CIRI2(
        ciri_ready_ch, 
        BWA_INDEX.out.collect(), 
        STRINGTIE_MERGE.out.merged_gtf
    )
    
    CDS_UTR_ANNOTATION(STRINGTIE_MERGE.out.merged_gtf, fasta_ch)
    
    bam_files_ch = STAR.out.bam.map { sample, group, bam -> bam }.collect()
    
    FEATURECOUNTS(bam_files_ch, STRINGTIE_MERGE.out.merged_gtf)

    bam_for_rmats_ch = STAR.out.bam.map { sample, group, bam -> bam }.collect()
    
    PREPARE_RMATS(file(params.sample_table), bam_for_rmats_ch)

    RMATS(
        PREPARE_RMATS.out.b1,
        PREPARE_RMATS.out.b2,
        STRINGTIE_MERGE.out.merged_gtf,
        params.trim_minlen,
        params.rmats_lib_type
    )


    multiqc_input_ch = Channel.value(file("${params.outdir}/fastqc_results"))
        .mix(STAR.out.log_final.collect())
        .collect()

    MULTIQC(multiqc_input_ch)

    PLOT_SPLICING(RMATS.out.results, file(params.sample_table))
    
    PLOT_EXPRESSION(FEATURECOUNTS.out.counts, file(params.sample_table), PLOT_SPLICING.out[0])
    
    RUN_DESEQ2(FEATURECOUNTS.out.counts, file(params.sample_table))

    EXTRACT_IDS(
        RUN_DESEQ2.out[0],
        PLOT_SPLICING.out[0]
    )
    
    GFFREAD(EXTRACT_IDS.out, STRINGTIE_MERGE.out.merged_gtf, fasta_ch)
    
    CPC2(GFFREAD.out)
    
    TRANSDECODER(GFFREAD.out)
    
    SPLICING_CONSEQUENCE(
        STRINGTIE_MERGE.out.merged_gtf,
        TRANSDECODER.out,                    
        PLOT_SPLICING.out[0],                
        CPC2.out,                        
        PLOT_SPLICING.out[1]                 
    )
    
    DEGS_DELS_HYBRID(
        RUN_DESEQ2.out[0],
        RUN_DESEQ2.out[1],
        CPC2.out,
        STRINGTIE_MERGE.out.merged_gtf,
        file(params.sample_table)
    )
    
    sashimi_ch = PLOT_SPLICING.out.sashimi_files
        .flatten()
        .map { file -> 
            def type = file.name.replaceAll('_sashimi.txt', '')
            return tuple(type, file) 
        }
        
     RMATS2SASHIMIPLOT(
       PREPARE_RMATS.out.b1, 
       PREPARE_RMATS.out.b2, 
       PREPARE_RMATS.out.g1_name,
       PREPARE_RMATS.out.g2_name,
       sashimi_ch
    )  
    
     freebayes_input_ch = STAR.out.bam.map { sample, group, bam -> 
        return [ [id: sample, group: group], bam ] 
    }

     FREEBAYES(
        freebayes_input_ch,
        fasta_ch
    )
    
    vcf_files_ch = FREEBAYES.out.vcf.map { meta, vcf -> vcf }.collect()

    MERGE_VCFS(vcf_files_ch)
    
    SNPEFF(
        MERGE_VCFS.out.merged_vcf, 
        CDS_UTR_ANNOTATION.out.gff3, 
        fasta_ch
    )
    
    VCF_TO_MATRIX(
        SNPEFF.out.vcf, 
        SNPEFF.out.index
    )
    
    SNV_POSTPROCESSING(
        VCF_TO_MATRIX.out.matrix,
        file(params.sample_table)
    )
    
}

process PLOT_SPLICING {
    container 'mordziarz/pipeline_rnaseq:latest'
    publishDir "${params.outdir}/splicing_plots_and_tables", mode: 'copy'

    input:
    path rmats_dir
    path samples_csv

    output:
    path "Splicing_significant_results.csv", optional: true  
    path "All_splicing_events.csv", optional: true  
    path "*.png", optional: true                           
    path "*_sashimi.txt", emit: sashimi_files, optional: true

    script:
    """
    Rscript ${baseDir}/bin/splicing_postprocessing.R . ${samples_csv}
    """
}

process PLOT_EXPRESSION {
    container 'mordziarz/pipeline_rnaseq:latest' 
    publishDir "${params.outdir}/kmeans_results", mode: 'copy'

    input:
    path feature_counts_txt
    path samples_csv
    path splicing_results_dir

    output:
    path "*.png"
    path "*.csv", optional: true

    script:
    """
    Rscript ${baseDir}/bin/expression_postprocessing.R ${feature_counts_txt} ${samples_csv} ${splicing_results_dir}
    """
}

process RUN_DESEQ2 {
    container 'mordziarz/pipeline_rnaseq:latest'
    publishDir "${params.outdir}/differential_expression", mode: 'copy'

    input:
    path feature_counts_txt
    path samples_csv

    output:
    path "Genes_sig.csv"     // Indeks [0]
    path "Expression.csv"    // Indeks [1]
    path "*.png", optional: true

    script:
    """
    Rscript ${baseDir}/bin/deseq2_analysis.R ${feature_counts_txt} ${samples_csv} .
    """
}

process EXTRACT_IDS {
    tag "Extracting significant IDs"
    container 'mordziarz/pipeline_rnaseq:latest'
    publishDir "${params.outdir}/functional_annotation", mode: 'copy'

    input:
    path genes_sig_file
    path splicing_csv

    output:
    path "target_ids.txt"

    script:
    """
    if [ -f "${genes_sig_file}" ]; then
        awk -F';' 'NR>1 {print \$1}' "${genes_sig_file}" | tr -d '"' > degs_ids.txt
    else
        touch degs_ids.txt
    fi

    if [ -f "${splicing_csv}" ]; then
        awk -F';' 'NR>1 {print \$1}' "${splicing_csv}" | tr -d '"' > splic_ids.txt
    else
        touch splic_ids.txt
    fi

    cat degs_ids.txt splic_ids.txt | grep -v '^\\s*\$' | grep -v 'GeneID' | sort -u > target_ids.txt

    echo "Wyciągnięto unikalnych ID:"
    wc -l target_ids.txt
    """
}

process SPLICING_CONSEQUENCE {
    container 'mordziarz/pipeline_rnaseq:latest'
    publishDir "${params.outdir}/consequence_analysis", mode: 'copy'

    input:
    path gtf_file
    path longest_orfs_gff
    path sig_splicing
    path cpc2_file
    path all_splicing

    output:
    path "consequence_matrix.csv"
    path "splicing_consequence_*.png", optional: true

    script:
    """
    mkdir -p output_consequences

    Rscript ${baseDir}/bin/splicing_consequences.R \
        "${gtf_file}" \
        "${longest_orfs_gff}/longest_orfs.gff3" \
        "${sig_splicing}" \
        "${cpc2_file}" \
        "${all_splicing}" \
        "output_consequences"

    if [ -d "output_consequences" ]; then
        cp -r output_consequences/* .
    fi
    """
}

process DEGS_DELS_HYBRID {
    container 'mordziarz/pipeline_rnaseq:latest'
    publishDir "${params.outdir}/DEGs_DELs_Hybrid", mode: 'copy'

    input:
    path genes_sig_file
    path expression_file
    path cpc2_file
    path gtf_file
    path samples_csv_path

    output:
    path "*", emit: results

    script:
    """
    Rscript ${baseDir}/bin/degs_dels_hybrid.R \
        "${genes_sig_file}" \
        "${expression_file}" \
        "${cpc2_file}" \
        "${gtf_file}" \
        "${samples_csv_path}" \
        "DEGs_DELs_Hybrid"
    """
}

process SNV_POSTPROCESSING {
    container 'mordziarz/pipeline_rnaseq:latest'
    publishDir "${params.outdir}/snv_plots_and_tables", mode: 'copy'

    input:
    path variants_matrix
    path samples_csv

    output:
    path "All_snv.csv"
    path "SNV_sig.csv", optional: true
    path "SNV_significant_results.csv", optional: true
    path "*.png", optional: true

    script:
    """
    Rscript ${baseDir}/bin/snv_postprocessing.R ${variants_matrix} ${samples_csv}
    """
}
