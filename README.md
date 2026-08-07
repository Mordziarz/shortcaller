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

The `results/DEGs_DELs_Hybrid` directory contains visualization plots and two primary tables:

1. `Genes_sig_with_coding_potential.csv` – A table containing DESeq2 statistical significance results along with functional classifications distinguishing whether a gene is classified as a DEL, DEG, or Hybrid (i.e., genes where not all transcripts possess coding potential).

**Example:**

| GeneID | baseMean | log2FoldChange | lfcSE | pvalue | padj | LAND_2 | LAND_3 | LAND_4 | WATER_1 | WATER_2 | WATER_3 | WATER_4 | Expression | Total_Isoforms | Protein_Coding_Isoforms | Non_Coding_Isoforms | Transcripts |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| CL.10056 | 400,59 | -3,34 | 0,40 | 4,98e-18 | 8,57e-16 | 1129,10 | 670,60 | 695,81 | 72,85 | 72,94 | 103,56 | 59,27 | Downregulated | 4 | 3 | 1 | Hybrid |
| CL.1007 | 1756,87 | -1,21 | 0,24 | 4,37e-08 | 1,93e-06 | 3012,97 | 2994,86 | 1933,54 | 1010,95 | 1132,17 | 1097,20 | 1116,41 | Downregulated | 6 | 5 | 1 | Hybrid |
| CL.10171 | 277,20 | 1,19 | 0,59 | 0,0027 | 0,0295 | 171,82 | 174,50 | 47,69 | 270,58 | 445,90 | 464,63 | 365,24 | Upregulated | 3 | 1 | 2 | Hybrid |
| CL.10202 | 2265,56 | -2,09 | 0,28 | 1,04e-14 | 1,23e-12 | 4718,89 | 3551,25 | 3939,70 | 1306,81 | 643,16 | 1069,21 | 629,88 | Downregulated | 4 | 2 | 2 | Hybrid |
| CL.10348 | 22,63 | 7,46 | 3,24 | 5,99e-05 | 0,0013 | 0 | 0 | 0 | 37,17 | 45,58 | 34,99 | 40,66 | Upregulated | 1 | 0 | 1 | DELs |
| CL.10553 | 2085,90 | 1,43 | 0,25 | 1,37e-09 | 7,85e-08 | 935,80 | 1134,72 | 984,11 | 2047,18 | 3016,90 | 2625,44 | 3857,12 | Upregulated | 2 | 2 | 0 | DEGs |
| CL.10555 | 50,50 | 4,36 | 1,39 | 5,23e-05 | 0,0012 | 0 | 7,31 | 0 | 83,26 | 130,95 | 69,97 | 62,02 | Upregulated | 1 | 0 | 1 | DELs |
| CL.1059 | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - |

2. `per_gene_pairwise_categories_correlation.csv` – A table showing pairwise expression profile correlations across categories (DEGs vs. DELs, DEGs vs. Hybrid, and DELs vs. Hybrid). It includes gene IDs, their respective categories, correlation values, raw p-values, and FDR-adjusted p-values.

**Example:**
| GeneID1 | GeneID2 | Correlation | P_Value | Category1 | Category2 | FDR |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| CL.10553 | CL.10348 | 0,9038 | 0,0052 | DEGs | DELs | 0,0286 |
| CL.10553 | CL.10555 | 0,7607 | 0,0471 | DEGs | DELs | 0,0822 |
| CL.10553 | CL.10795 | -0,6383 | 0,1229 | DEGs | DELs | 0,1534 |
| CL.10553 | CL.10801 | -0,6585 | 0,1078 | DEGs | DELs | 0,1396 |
| CL.10553 | CL.10802 | -0,9140 | 0,0040 | DEGs | DELs | 0,0257 |
| CL.10553 | CL.10895 | -0,6554 | 0,1100 | DEGs | DELs | 0,1417 |
| CL.10553 | CL.10965 | 0,7039 | 0,0775 | DEGs | DELs | 0,1117 |
| CL.10553 | CL.11114 | -0,8336 | 0,0198 | DEGs | DELs | 0,0516 |
| CL.10553 | CL.11337 | 0,8577 | 0,0136 | DEGs | DELs | 0,0430 |
| CL.10553 | CL.11702 | 0,5394 | 0,2115 | DEGs | DELs | 0,2344 |
| CL.10553 | CL.12460 | 0,9900 | 1,91e-05 | DEGs | DELs | 0,0063 |
| CL.10553 | CL.12757 | -0,7110 | 0,0733 | DEGs | DELs | 0,1078 |
| CL.10553 | CL.12834 | 0,8304 | 0,0207 | DEGs | DELs | 0,0527 |
| CL.10553 | CL.12940 | 0,2234 | 0,6301 | DEGs | DELs | 0,6347 |
| CL.10553 | CL.12984 | -0,7387 | 0,0579 | DEGs | DELs | 0,0930 |
| CL.10553 | CL.12989 | 0,6493 | 0,1145 | DEGs | DELs | 0,1458 |

### Visualizations

* **Volcano plot** highlighting DEGs, DELs, Hybrid, as well as up- and down-regulated genes:
  ![Volcano_DEGs](plots/Volcano_DEGs_DELs_Hybrid.png)

* **MA plot** highlighting DEGs, DELs, Hybrid, as well as up- and down-regulated genes:
  ![MA_DEGs](plots/MA_plot_DEGs_DELs_Hybrid.png)

* **Example heatmap** showing DEGs (outputs for DELs and Hybrid are also generated):
  ![Heatmap_DEGs](plots/Heatmap_DEGs.png)

* **Example heatmap** displaying the top 10 most significant DEGs (outputs for DELs and Hybrid are also generated):
  ![Heatmap_DEGs](plots/Heatmap_Top10_DELs.png)

* **Circos plot** integrating the results from this analysis step:
  1. The 3 outermost heatmap tracks corresponding to DEGs, DELs, and Hybrid.
  2. A dot plot track indicating whether a given molecule is up- or down-regulated.
  3. Internal links connecting molecules with correlated expression profiles (default: FDR < 0.05 and r > 0.99).
  ![circos_degs_dels_hybrid](plots/circos_DEGs_DELs_Hybrid.png)

## Kmeans Genes Analysis

![expression_cluster1](plots/expression_cluster_g1.png)
![expression_cluster2](plots/expression_cluster_g2.png)
![expression_cluster3](plots/expression_cluster.png)
![expression_cluster4](plots/expression_cluster_type.png)