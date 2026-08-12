#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
featurecounts_file <- args[1]
samples_csv_path   <- args[2]
rmats_dir          <- args[3]

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)
library(purrr)

cat("--- STARTING EXPRESSION ANALYSIS AND SPLICING INTEGRATION ---\n")

samples_df <- read.csv(samples_csv_path)

unique_groups <- unique(samples_df$group)
if(length(unique_groups) < 2) {
  stop("The samples.csv file must contain at least two different groups!")
}

g1_name <- unique_groups[1]
g2_name <- unique_groups[2]

samples_g1 <- samples_df %>% filter(group == g1_name) %>% pull(sample)
samples_g2 <- samples_df %>% filter(group == g2_name) %>% pull(sample)

counts_data <- read.table(featurecounts_file, header=TRUE, skip=1)

colnames(counts_data) <- gsub("_Aligned.sortedByCoord.out.bam", "", colnames(counts_data))
colnames(counts_data) <- gsub("_Aligned.out.bam", "", colnames(counts_data))

counts_to_tpm <- function(counts, lengths) {
  rpk <- counts / (lengths / 1000)
  scaling_factor <- sum(rpk, na.rm = TRUE) / 1e6
  if(scaling_factor == 0) return(rep(0, length(counts)))
  tpm <- rpk / scaling_factor
  return(tpm)
}

tpm_data <- counts_data

all_sample_cols <- c(samples_g1, samples_g2)
for (s in all_sample_cols) {
  if (s %in% colnames(tpm_data)) {
    tpm_col_name <- paste0("TPM_", s)
    tpm_data[[tpm_col_name]] <- counts_to_tpm(tpm_data[[s]], tpm_data$Length)
  }
}

tpm_cols_g1 <- paste0("TPM_", samples_g1)
tpm_cols_g2 <- paste0("TPM_", samples_g2)

tpm_data$Mean_TPM_G1 <- rowMeans(tpm_data[, tpm_cols_g1, drop=FALSE], na.rm = TRUE)
tpm_data$Mean_TPM_G2 <- rowMeans(tpm_data[, tpm_cols_g2, drop=FALSE], na.rm = TRUE)

tpm_summary <- tpm_data %>%
  dplyr::select(Geneid, Length, Mean_TPM_G1, Mean_TPM_G2)

g1_for_km <- tpm_summary %>% filter(Mean_TPM_G1 > 0) %>% mutate(log_TPM = log2(Mean_TPM_G1 + 1))
g2_for_km <- tpm_summary %>% filter(Mean_TPM_G2 > 0) %>% mutate(log_TPM = log2(Mean_TPM_G2 + 1))

run_kmeans_stratification <- function(df) {
  set.seed(123)
  km_result <- kmeans(df$log_TPM, centers = 3, nstart = 25)
  df$Cluster <- km_result$cluster
  
  cluster_order <- df %>% 
    group_by(Cluster) %>% 
    summarise(m = mean(log_TPM)) %>% 
    arrange(m) %>% 
    pull(Cluster)
  
  df %>% mutate(
    Expression_Class = case_when(
      Cluster == cluster_order[1] ~ "Low",
      Cluster == cluster_order[2] ~ "Medium",
      Cluster == cluster_order[3] ~ "High"
    )
  )
}

g1_classified <- run_kmeans_stratification(g1_for_km) %>% dplyr::select(Geneid, Class_G1 = Expression_Class)
g2_classified <- run_kmeans_stratification(g2_for_km) %>% dplyr::select(Geneid, Class_G2 = Expression_Class)

final_classification <- tpm_summary %>%
  left_join(g1_classified, by = "Geneid") %>%
  left_join(g2_classified, by = "Geneid") %>%
  mutate(
    Class_G1 = replace_na(Class_G1, "No Expression"),
    Class_G2 = replace_na(Class_G2, "No Expression")
  )

g1_res <- run_kmeans_stratification(g1_for_km)
g2_res <- run_kmeans_stratification(g2_for_km)

pca_like_data <- g1_res %>%
  dplyr::select(Geneid, log_TPM_G1 = log_TPM, Class_G1 = Expression_Class) %>%
  inner_join(
    g2_res %>% dplyr::select(Geneid, log_TPM_G2 = log_TPM, Class_G2 = Expression_Class),
    by = "Geneid"
  )

pca_like_data$Class_G1 <- factor(pca_like_data$Class_G1, levels = c("Low", "Medium", "High"))
pca_like_data$Class_G2 <- factor(pca_like_data$Class_G2, levels = c("Low", "Medium", "High"))

p1 <- ggplot(pca_like_data, aes(x = log_TPM_G1, y = log_TPM_G2, color = Class_G1)) +
  geom_point(size = 0.3, alpha = 0.2) +
  scale_color_manual(values = c("Low" = "plum3", "Medium" = "aquamarine2", "High" = "magenta3")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  labs(x = paste(g1_name, "expression (mean TPM)"), y = paste(g2_name, "expression (mean TPM)"), color = "Class", title = "") +
  theme_bw() + theme(legend.position = "none")

p2 <- ggplot(pca_like_data, aes(x = log_TPM_G1, y = log_TPM_G2, color = Class_G2)) +
  geom_point(size = 0.3, alpha = 0.2) +
  scale_color_manual(values = c("Low" = "plum3", "Medium" = "aquamarine2", "High" = "magenta3")) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  labs(x = paste(g1_name, "expression (mean TPM)"), y = paste(g2_name, "expression (mean TPM)"), color = "Class", title = "") +
  theme_bw() + theme(legend.position = "none")

combined_plot <- p1 + p2 + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom") & 
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

ggsave("expression_cluster_g1.png", plot = p1, width = 4, height = 4, dpi = 600, units = "in")
ggsave("expression_cluster_g2.png", plot = p2, width = 4, height = 4, dpi = 600, units = "in")
ggsave("expression_cluster.png", plot = combined_plot, width = 6, height = 6, dpi = 600, units = "in")

cat("--- INTEGRATION WITH SIGNIFICANT SPLICING RESULTS ---\n")
sig_splicing_files <- list.files(".", pattern = "^Splicing_significant_results.*\\.csv$", full.names = TRUE)

if(length(sig_splicing_files) > 0) {
  sig_file_to_read <- sig_splicing_files[order(file.mtime(sig_splicing_files), decreasing = TRUE)][1]
  cat("Found significant splicing table:", sig_file_to_read, "\n")
  
  wszystkie_splicingi <- read.csv2(sig_file_to_read, stringsAsFactors = FALSE)
  
  if("GeneID" %in% colnames(wszystkie_splicingi) && nrow(wszystkie_splicingi) > 0) {
    
    if(!"Type" %in% colnames(wszystkie_splicingi) && "Event" %in% colnames(wszystkie_splicingi)) {
      wszystkie_splicingi$Type <- wszystkie_splicingi$Event
    }
    
    splicing_expressed <- wszystkie_splicingi %>%
      left_join(final_classification, by = c("GeneID" = "Geneid"))
    
    plot_data_final <- splicing_expressed %>%
      rename(Class = Class_G1) %>%
      filter(Class != "No Expression" & !is.na(Class)) %>%
      distinct(Type, GeneID, .keep_all = TRUE) %>%
      dplyr::select(Type, Class)
    
    if(nrow(plot_data_final) > 0) {
      plot_data_final$Class <- factor(plot_data_final$Class, levels = c("Low", "Medium", "High"))
      
      p_type <- ggplot(plot_data_final, aes(x = Class, fill = Class)) +
        geom_bar(width = 0.7) +
        facet_wrap(~Type, scales = "free_y") +
        scale_y_continuous(labels = scales::label_number(accuracy = 1)) +
        scale_fill_manual(values = c("Low" = "plum3", "Medium" = "aquamarine2", "High" = "magenta3")) +
        labs(title = "", x = "K-means Expression Class", y = "Number of Splicing Events", fill = "Expression Class") +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 9),
          axis.text.y = element_text(size = 9),
          strip.text = element_text(face = "bold", size = 10),
          panel.spacing = unit(1.2, "lines"),
          legend.position = "bottom",
          panel.grid.minor = element_blank()
        )
      
      ggsave("expression_cluster_type.png", plot = p_type, width = 8, height = 4, dpi = 600, units = "in")
    }
  } else {
    cat("Significant splicing table is empty. Skipping integrative plot generation.\n")
  }
} else {
  cat("Missing significant splicing file (Splicing_significant_results_*.csv). Skipping integrative plot generation.\n")
}

write.csv2(final_classification, "expression_classification.csv", row.names = FALSE)
cat("--- EXPRESSION ANALYSIS COMPLETED SUCCESSFULLY ---\n")
