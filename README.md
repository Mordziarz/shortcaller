# shortcaller
A Nextflow pipeline for short-read RNA-seq analysis, supporting differential expression and alternative splicing.

## Current Features
* **Quality Control & Trimming:** FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ ) and Fastp (https://doi.org/10.1093/bioinformatics/bty560).
* **Alignment:** STAR (https://doi.org/10.1093/bioinformatics/bts635).
* **Transcript Assembly & Quantification:** StringTie (https://doi.org/10.1038/nbt.3122) and FeatureCounts (https://doi.org/10.1093/bioinformatics/btt656).
* **Differential Expression Analysis:** DESeq2 (https://doi.org/10.1186/s13059-014-0550-8).
* **Alternative Splicing Analysis:** rMATS (https://doi.org/10.1073/pnas.1419161111) , sashimi plots (https://github.com/Xinglab/rmats2sashimiplot).
* **Downstream Analysis:** CPC2 (coding potential; https://doi.org/10.1093/nar/gkx428) and TransDecoder (NMD prediciton; https://github.com/TransDecoder/TransDecoder).

## Upcoming Features
* Single Nucleotide Variant (SNV) calling
* Alternative Polyadenylation (APA)
* Circular RNA (circRNA) detection

## Prerequisites & Environment
The pipeline relies on a pre-built Docker container containing all necessary dependencies:
**[mordziarz/pipeline_rnaseq on Docker Hub](https://hub.docker.com/r/mordziarz/pipeline_rnaseq)**

**Recommended Setup:**
1. Install **Nextflow** via Conda.
2. Ensure **Docker** is installed and running on your system.

## Input Data Format
Create a simple CSV file defining your samples, paths to paired-end reads, and experimental groups:

Note: quotation marks are not allowed in the table

```csv
sample,fastq_1,fastq_2,group
A_1,A1_1.test.fq.gz,A1_2.test.fq.gz,A
A_2,A2_1.test.fq.gz,A2_2.test.fq.gz,A
A_3,A3_1.test.fq.gz,A3_2.test.fq.gz,A
A_4,A4_1.test.fq.gz,A4_2.test.fq.gz,A
B_5,B1_1.test.fq.gz,B1_2.test.fq.gz,B
B_6,B2_1.test.fq.gz,B2_2.test.fq.gz,B
B_7,B3_1.test.fq.gz,B3_2.test.fq.gz,B
B_8,B4_1.test.fq.gz,B4_2.test.fq.gz,B
```
## Basic usage
```
nextflow run main.nf --sample_table samples_table.csv --fasta reference_genome.fasta --gtf reference_gtf.gtf -profile standard
```

## Parameters that can be defined

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `sample_table` | `null` | Path to the samples CSV metadata table (Required). |
| `outdir` | `"results"` | Output directory for pipeline results. |
| `fasta` | `null` | Path to the reference genome FASTA file (Required). |
| `gtf` | `null` | Path to the reference annotation GTF file (Required). |
| `strand` | `"rf"` | Library strand-specificity (`rf` for reverse-stranded, `fr` for forward, etc.). |
| `trim_minlen` | `140` | Minimum length requirement for reads to pass filtering. |
| `trim_avgqual` | `20` | Minimum average quality score (Phred) required for reads. |
| `rmats_lib_type` | `"fr-firststrand"` | Library type parameter for rMATS. |
| `threads_fastqc` | `2` | CPU threads allocated to FastQC. |
| `threads_fastp` | `8` | CPU threads allocated to Fastp trimming. |
| `threads_star_index` | `8` | CPU threads allocated for building the STAR index. |
| `threads_star` | `8` | CPU threads allocated to STAR alignment. |
| `threads_stringtie` | `4` | CPU threads allocated to StringTie. |
| `threads_featurecounts` | `4` | CPU threads allocated to FeatureCounts. |
| `threads_rmats` | `8` | CPU threads allocated to rMATS. |

## Differentially Expressed Genes Analysis

## Differentially Alternative Splicing Analysis

## Differentially Expressed Genes Analysis (DEGs, DELs, Hybrid)

## Kmeans Genes Analysis
