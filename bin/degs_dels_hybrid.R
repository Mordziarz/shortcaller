#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
genes_sig_file     <- args[1]
expression_file    <- args[2]
cpc2_file          <- args[3]
gtf_file           <- args[4]
samples_csv_path   <- args[5]

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(corrplot)
  library(rtracklayer)
  library(gridBase)
})

options(scipen = 999)

cat("--- Loading samples.csv to determine groups ---\n")
samples_df <- read.csv(samples_csv_path)
unique_groups <- unique(samples_df$group)

if(length(unique_groups) < 2) {
  stop("The samples.csv file must contain at least two different groups for comparison!")
}

ref_group <- unique_groups[1]
test_group <- unique_groups[2]

cat("Reference group (ref):", ref_group, "\n")
cat("Test group (test):", test_group, "\n")

samples_ref <- samples_df %>% filter(group == ref_group) %>% pull(sample)
samples_test <- samples_df %>% filter(group == test_group) %>% pull(sample)
sample_names <- c(samples_ref, samples_test)

cat("--- Loading GTF file ---\n")
gtf <- import(gtf_file)
gtf_df <- as.data.frame(gtf)
gtf_df <- gtf_df[gtf_df$type %in% "transcript",]
gtf_df_genes <- gtf_df[,c("transcript_id","gene_id")]
gtf_df_genes <- gtf_df_genes[!duplicated(gtf_df_genes$transcript_id),]

cat("--- Loading input files ---\n")
genes_sig <- read.csv2(genes_sig_file, stringsAsFactors = FALSE)
expression_df <- read.csv2(expression_file, stringsAsFactors = FALSE)

genes_sig$Transcripts <- NULL

gene_hybrid_status <- data.frame()

if(file.exists(cpc2_file)) {
  cpc2 <- read.csv(cpc2_file, header=TRUE, sep="\t")
  cpc2 <- merge(cpc2, gtf_df_genes, by.x="X.ID", by.y="transcript_id")
  
  gene_hybrid_status <- cpc2 %>%
    mutate(is_nc = tolower(label) %in% c("noncoding", "non-coding")) %>%
    group_by(gene_id) %>%
    summarise(
      Total_Isoforms = n(),
      Protein_Coding_Isoforms = sum(!is_nc),
      Non_Coding_Isoforms = sum(is_nc),
      Transcripts = case_when(
        Protein_Coding_Isoforms > 0 & Non_Coding_Isoforms > 0 ~ "Hybrid",
        Protein_Coding_Isoforms > 0 & Non_Coding_Isoforms == 0 ~ "DEGs",
        Protein_Coding_Isoforms == 0 & Non_Coding_Isoforms > 0 ~ "DELs",
        TRUE ~ "Other"
      ),
      .groups = 'drop'
    )
}

if(nrow(gene_hybrid_status) > 0) {
  genes_sig <- merge(genes_sig, gene_hybrid_status, by.x="GeneID", by.y="gene_id", all.x = TRUE)
  genes_sig$Transcripts[is.na(genes_sig$Transcripts)] <- "DEGs"
} else {
  genes_sig$Transcripts <- "DEGs"
}

write.csv2(genes_sig, file.path("Genes_sig_with_coding_potential.csv"))

Final_ekspresja <- expression_df

genes_sig <- genes_sig[,c("GeneID","Total_Isoforms","Protein_Coding_Isoforms","Non_Coding_Isoforms","Expression","Transcripts")]

Final_ekspresja <- merge(Final_ekspresja, genes_sig, by.x="GeneID", by.y="GeneID", all.x=T)

Final_ekspresja$Expression[is.na(Final_ekspresja$Expression)] <- "No significant"
Final_ekspresja$Transcripts[is.na(Final_ekspresja$Transcripts)] <- "No significant"

cat("--- Generating Volcano plot with Hybrid / DELs / DEGs split ---\n")
options(ggrepel.max.overlaps = Inf)

Final_ekspresja$Symbol <- NA

cat("--- Generating Volcano plot with Hybrid Categories ---\n")

options(ggrepel.max.overlaps = Inf)
volcano_p <- ggplot(Final_ekspresja, aes(x = log2FoldChange, y = -log10(padj), color = Transcripts, shape = Expression, fill = Transcripts)) +
  geom_point(alpha = 0.7, size = 2) + 
  geom_text_repel(show.legend = FALSE, aes(label = Symbol), force = 2, nudge_y = 1, nudge_x = -0.25, color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", linewidth = 0.5, color = "red") +
  geom_hline(yintercept = 1.30103, linetype = "dotted", linewidth = 0.4, color = "red") +
  xlab(expression("log2FoldChange")) + ylab("-log10(padj)") +
  scale_color_manual(values = c("DELs" = "blue3", "DEGs" = "green3", "Hybrid" = "purple3", "No significant" = "grey25")) +
  scale_shape_manual(values = c("Upregulated" = 21, "Downregulated" = 24, "No significant" = 22)) + 
  scale_fill_manual(values = c("DELs" = "blue", "DEGs" = "green", "Hybrid" = "purple", "No significant" = "grey50")) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_blank(), 
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16), 
    legend.position = "right"
  ) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 4)))

ggsave(file.path("Volcano_DEGs_DELs_Hybrid.png"), plot = volcano_p, width = 6, height = 6, dpi = 300)

cat("--- Generating MA plot with Hybrid / DELs / DEGs split ---\n")

MA_plot <- Final_ekspresja

match_samples <- intersect(sample_names, colnames(MA_plot))

if(length(match_samples) > 0) {
  MA_plot$Means <- log2(rowMeans(sapply(MA_plot[, match_samples, drop = FALSE], as.numeric), na.rm = TRUE) + 1)
} else if(all(MA_plot$GeneID %in% expression_df$GeneID)) {
  merged_ma <- merge(MA_plot, expression_df, by = "GeneID", suffixes = c("", "_expr"))
  match_samples_expr <- intersect(sample_names, colnames(merged_ma))
  if(length(match_samples_expr) > 0) {
    MA_plot$Means <- log2(rowMeans(sapply(merged_ma[, match_samples_expr, drop = FALSE], as.numeric), na.rm = TRUE) + 1)
  } else { 
    MA_plot$Means <- 0 
  }
} else { 
  MA_plot$Means <- 0 
}

if(!"Symbol" %in% colnames(MA_plot)) {
  MA_plot$Symbol <- NA
}

maplot_p <- ggplot(MA_plot, aes(x = Means, y = log2FoldChange, color = Transcripts, shape = Expression, fill = Transcripts)) +
  geom_point(alpha = 0.7, size = 2) + 
  xlab(expression("log2(Mean of normalized counts)")) + 
  ylab("log2FoldChange") +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted", linewidth = 0.4, color = "red") +
  scale_color_manual(values = c("DELs" = "blue3", "DEGs" = "green3", "Hybrid" = "purple3", "No significant" = "grey25")) +
  scale_shape_manual(values = c("Upregulated" = 21, "Downregulated" = 24, "No significant" = 22)) + 
  scale_fill_manual(values = c("DELs" = "blue", "DEGs" = "green", "Hybrid" = "purple", "No significant" = "grey50")) +
  geom_text_repel(show.legend = FALSE, aes(label = Symbol), force = 2, nudge_y = 1, nudge_x = -0.25, color = "black") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_blank(), 
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16), 
    legend.position = "right"
  ) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 4)))

ggsave(file.path("MA_plot_DEGs_DELs_Hybrid.png"), plot = maplot_p, width = 10, height = 6, dpi = 300)

cat("--- Generating Per-Gene Pairwise Correlation Analysis (Optimized Vectorized Version) ---\n")

if(all(c("GeneID") %in% colnames(Final_ekspresja)) && all(sample_names %in% colnames(expression_df))) {
  
  gene_lists <- list(
    DEGs = Final_ekspresja %>% filter(Transcripts == "DEGs") %>% pull(GeneID),
    DELs = Final_ekspresja %>% filter(Transcripts == "DELs") %>% pull(GeneID),
    Hybrid = Final_ekspresja %>% filter(Transcripts == "Hybrid") %>% pull(GeneID)
  )
  
  expr_mat <- as.matrix(expression_df[, intersect(sample_names, colnames(expression_df)), drop = FALSE])
  rownames(expr_mat) <- expression_df$GeneID
  
  summary_per_gene_corr <- data.frame()
  category_pairs <- combn(names(gene_lists), 2, simplify = FALSE)
  
  for(pair in category_pairs) {
    cat1 <- pair[1]
    cat2 <- pair[2]
    
    genes1 <- intersect(gene_lists[[cat1]], rownames(expr_mat))
    genes2 <- intersect(gene_lists[[cat2]], rownames(expr_mat))
    
    if(length(genes1) > 0 && length(genes2) > 0) {
      sub1 <- expr_mat[genes1, , drop = FALSE]
      sub2 <- expr_mat[genes2, , drop = FALSE]
      
      cor_matrix <- cor(t(sub1), t(sub2), method = "pearson", use = "complete.obs")
      
      n <- ncol(sub1)
      if(n > 2) {
        t_stat <- cor_matrix * sqrt((n - 2) / (1 - cor_matrix^2))
        pval_matrix <- 2 * pt(-abs(t_stat), df = n - 2)
        
        df_cor <- as.data.frame(as.table(cor_matrix))
        colnames(df_cor) <- c("GeneID1", "GeneID2", "Correlation")
        
        df_pval <- as.data.frame(as.table(pval_matrix))
        colnames(df_pval) <- c("GeneID1", "GeneID2", "P_Value")
        
        pair_df <- merge(df_cor, df_pval, by = c("GeneID1", "GeneID2"))
        pair_df$Category1 <- cat1
        pair_df$Category2 <- cat2
        
        pair_df <- pair_df %>% 
          filter(!is.na(Correlation)) %>%
          mutate(Correlation = round(Correlation, 4),
                 FDR = p.adjust(P_Value, method = "fdr"))
        
        summary_per_gene_corr <- rbind(summary_per_gene_corr, pair_df)
      }
    }
  }
  
  if(nrow(summary_per_gene_corr) > 0) {
    output_per_gene_csv <- file.path("per_gene_pairwise_categories_correlation.csv")
    write.csv2(summary_per_gene_corr, output_per_gene_csv, row.names = FALSE)
    cat("Saved fast per-gene correlation table to:", output_per_gene_csv, "\n")
  } else {
    message("Note: No valid per-gene correlations could be calculated.")
  }
  
} else {
  message("Note: Missing required columns for per-gene correlation analysis.")
}

cat("--- Generating Separate Heatmaps and 3-Sector Circos Plots for DEGs, DELs, Hybrid ---\n")

categories_list <- c("DEGs", "DELs", "Hybrid")
cat_colors_map <- c("DEGs" = "green3", "DELs" = "blue3", "Hybrid" = "purple3")

for(cat_name in categories_list) {
  sub_df <- Final_ekspresja %>% filter(Transcripts == cat_name & Expression != "No significant")
  if(nrow(sub_df) > 1) {
    mat_data <- as.matrix(expression_df[expression_df$GeneID %in% sub_df$GeneID, intersect(sample_names, colnames(expression_df)), drop = FALSE])
    rownames(mat_data) <- sub_df$GeneID[match(rownames(mat_data), sub_df$GeneID)]
    
    if(nrow(mat_data) > 1) {
      mat_scaled <- t(scale(t(log2(mat_data + 1))))
      cond_colors <- c("steelblue", "tomato")
      names(cond_colors) <- c(ref_group, test_group)
      
      ha <- HeatmapAnnotation(Condition = rep(c(ref_group, test_group), c(length(samples_ref), length(samples_test))),
                              col = list(Condition = cond_colors), simple_anno_size = unit(0.4, "cm"))
      
      png(file.path(paste0("Heatmap_", cat_name, ".png")), width = 9, height = 8, units = "in", res = 300)
      ht <- Heatmap(mat_scaled, name = "Z-score", show_row_names = FALSE,
                    column_title = paste0("Heatmap for ", cat_name), top_annotation = ha, heatmap_width = unit(6, "inch"))
      draw(ht, heatmap_legend_side = "right", merge_legends = TRUE)
      dev.off()
    }
  }
}

rownames(Final_ekspresja) <- Final_ekspresja$GeneID

for(cat_name in categories_list) {
  top_sub <- Final_ekspresja %>% 
    filter(Transcripts == cat_name) %>% 
    arrange(padj) %>% 
    head(10)
  
  if(nrow(top_sub) > 0) {
    sample_cols <- intersect(sample_names, colnames(top_sub))
    mat_data <- as.matrix(top_sub[, sample_cols, drop = FALSE])
    
    row_labels <- ifelse(is.na(top_sub$Symbol) | top_sub$Symbol == "", top_sub$GeneID, top_sub$Symbol)
    rownames(mat_data) <- row_labels
    
    if(nrow(mat_data) > 1) {
      mat_log <- log2(mat_data + 1)
      mat_scaled <- t(scale(t(mat_log)))
      mat_scaled[is.na(mat_scaled)] <- 0
      
      cond_colors <- c("steelblue", "tomato")
      names(cond_colors) <- c(ref_group, test_group)
      
      ha_top <- HeatmapAnnotation(Condition = rep(c(ref_group, test_group), c(length(samples_ref), length(samples_test))),
                                  col = list(Condition = cond_colors), simple_anno_size = unit(0.4, "cm"))
      
      png(file.path(paste0("Heatmap_Top10_", cat_name, ".png")), width = 9, height = 8, units = "in", res = 300)
      ht_top <- Heatmap(mat_scaled, name = "Z-score", 
                        show_row_names = TRUE,
                        row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                        column_title = paste0("Top ", nrow(mat_scaled), " ", cat_name, " Heatmap"), 
                        top_annotation = ha_top)
      draw(ht_top, heatmap_legend_side = "right", merge_legends = TRUE)
      dev.off()
    }
  }
}

cat("--- Generating Circos Plot for DEGs, DELs, Hybrid with exact legend style from reference script ---\n")

categories_list <- c("DEGs", "DELs", "Hybrid")
filtered_df <- Final_ekspresja %>% filter(Transcripts %in% categories_list & Expression != "No significant")

types_subset <- intersect(categories_list, unique(filtered_df$Transcripts))
filtered_df$Transcripts <- factor(filtered_df$Transcripts, levels = types_subset)

All <- filtered_df[order(filtered_df$Transcripts), ]
All$Numer <- 1:nrow(All)
row.names(All) <- All$Numer

All_gene_circos <- All[,c("GeneID","Numer")]

summary_per_gene_corr_sig <- summary_per_gene_corr[summary_per_gene_corr$FDR < 0.05,]
summary_per_gene_corr_sig <- summary_per_gene_corr_sig[abs(summary_per_gene_corr_sig$Correlation) > 0.99,]

summary_per_gene_corr_sig <- merge(summary_per_gene_corr_sig,All_gene_circos,by.x="GeneID1",by.y="GeneID")
names(summary_per_gene_corr_sig)[names(summary_per_gene_corr_sig) == "Numer"] <- "Numer_1"

summary_per_gene_corr_sig <- merge(summary_per_gene_corr_sig,All_gene_circos,by.x="GeneID2",by.y="GeneID")
names(summary_per_gene_corr_sig)[names(summary_per_gene_corr_sig) == "Numer"] <- "Numer_2"

target_samples <- c(samples_ref, samples_test)
sample_cols <- intersect(target_samples, colnames(All))

DE_all_FPKM <- as.matrix(All[, sample_cols])
rownames(DE_all_FPKM) <- All$Numer

DE_all_FPKM_null_log = log2(DE_all_FPKM + 1)
DE_all_FPKM_null_log = t(scale(t(DE_all_FPKM_null_log), scale = FALSE))
DE_all_FPKM_null_log[is.na(DE_all_FPKM_null_log)] <- 0

split_sizes <- sapply(types_subset, function(t) sum(All$Transcripts == t))
split <- factor(rep(types_subset, split_sizes), levels = types_subset)

idx_list <- split(1:nrow(All), All$Transcripts)

matrix_parts <- list()
for(t_name in types_subset) {
  if(length(idx_list[[t_name]]) > 0) {
    matrix_parts[[t_name]] <- DE_all_FPKM_null_log[idx_list[[t_name]], , drop = FALSE]
  }
}

offset_vals <- c(0, 200, 300)
names(offset_vals) <- types_subset

for(t_name in types_subset) {
  if(t_name != types_subset[1] && length(matrix_parts[[t_name]]) > 0) {
    off <- offset_vals[t_name]
    matrix_parts[[t_name]] <- ifelse(matrix_parts[[t_name]] >= 0, 
                                     matrix_parts[[t_name]] + off, 
                                     matrix_parts[[t_name]] - off)
  }
}

sample_cols_rev <- rev(sample_cols)

Heatmap_all <- do.call(rbind, matrix_parts)
Heatmap_all <- Heatmap_all[, sample_cols_rev, drop = FALSE]

col_fun_main = colorRamp2(
  c(-302, -300, -202, -200, -2, 0, 2, 200, 202, 300, 302),
  c("green", "black", "blue", "black", "green", "black", "red", "black", "red", "black", "yellow")
)

png(file.path("circos_DEGs_DELs_Hybrid.png"), width = 14, height = 16, units = 'in', res = 600)
plot.new()
circle_size = unit(1, "snpc") 
pushViewport(viewport(x = 0.5, y = 1, width = circle_size, height = circle_size, just = c("center", "top")))
par(omi = gridOMI(), new = TRUE)

circos.par(gap.after = rep(8, length(types_subset)), track.margin = c(0.02, 0.01))

circos.heatmap(Heatmap_all, split = split, col = col_fun_main, dend.side = NULL, 
               track.height = 0.4)

n_ref <- length(samples_ref)
n_test <- length(samples_test)
total_samples <- n_ref + n_test

for(i in seq_along(types_subset)) {
  circos.track(track.index = 1, ylim = c(0, total_samples), panel.fun = function(x, y) {
    if(CELL_META$sector.numeric.index == i) {
      circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), 0,
                  CELL_META$cell.xlim[2] + convert_x(6.5, "mm"), n_ref,
                  col = "white", border = "black")
      circos.text(CELL_META$cell.xlim[2] + convert_x(3.5, "mm"), n_ref / 2,
                  ref_group, cex = 0.9, facing = "clockwise")
      
      circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), n_ref,
                  CELL_META$cell.xlim[2] + convert_x(6.5, "mm"), total_samples,
                  col = "white", border = "black")
      circos.text(CELL_META$cell.xlim[2] + convert_x(3.5, "mm"), n_ref + (n_test / 2),
                  test_group, cex = 0.9, facing = "clockwise")
    }
  }, bg.border = NA)
}

row_mean <- as.numeric(All$log2FoldChange)
names(row_mean) <- 1:nrow(All)

circos.track(ylim = range(row_mean, na.rm = TRUE), track.height = 0.15, panel.fun = function(x, y) {
  current_subset <- CELL_META$subset
  if (length(current_subset) > 0) {
    y_sub = row_mean[current_subset]
    if (!is.null(CELL_META$row_order)) {
      y_sub = y_sub[CELL_META$row_order]
    }
    circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "grey")
    circos.points(pch = 16, seq_along(y_sub) - 0.5, y_sub, col = ifelse(y_sub > 0, "red", "blue"))
  }
}, cell.padding = c(0.02, 0, 0.02, 0))

has_pos_degs_hybrid <- FALSE
has_pos_dels_hybrid <- FALSE
has_pos_degs_dels  <- FALSE
has_neg_degs_hybrid <- FALSE
has_neg_dels_hybrid <- FALSE
has_neg_degs_dels  <- FALSE

if (nrow(summary_per_gene_corr_sig) > 0) {
  for(i in seq_len(nrow(summary_per_gene_corr_sig))) {
    c1 <- as.character(summary_per_gene_corr_sig$Category1[i])
    c2 <- as.character(summary_per_gene_corr_sig$Category2[i])
    r  <- summary_per_gene_corr_sig$Correlation[i]
    pair_str <- paste(sort(c(c1, c2)), collapse = "-")
    
    if(r > 0 && pair_str == "DEGs-Hybrid") {
      circos.heatmap.link(summary_per_gene_corr_sig$Numer_1[i], summary_per_gene_corr_sig$Numer_2[i], col = "cyan4")
      has_pos_degs_hybrid <- TRUE
    } else if(r > 0 && pair_str == "DELs-Hybrid") {
      circos.heatmap.link(summary_per_gene_corr_sig$Numer_1[i], summary_per_gene_corr_sig$Numer_2[i], col = "blue")
      has_pos_dels_hybrid <- TRUE
    } else if(r > 0 && pair_str == "DEGs-DELs") {
      circos.heatmap.link(summary_per_gene_corr_sig$Numer_1[i], summary_per_gene_corr_sig$Numer_2[i], col = "springgreen3")
      has_pos_degs_dels <- TRUE
    } else if(r < 0 && pair_str == "DEGs-Hybrid") {
      circos.heatmap.link(summary_per_gene_corr_sig$Numer_1[i], summary_per_gene_corr_sig$Numer_2[i], col = "purple")
      has_neg_degs_hybrid <- TRUE
    } else if(r < 0 && pair_str == "DELs-Hybrid") {
      circos.heatmap.link(summary_per_gene_corr_sig$Numer_1[i], summary_per_gene_corr_sig$Numer_2[i], col = "orange3")
      has_neg_dels_hybrid <- TRUE
    } else if(r < 0 && pair_str == "DEGs-DELs") {
      circos.heatmap.link(summary_per_gene_corr_sig$Numer_1[i], summary_per_gene_corr_sig$Numer_2[i], col = "brown")
      has_neg_degs_dels <- TRUE
    }
  }
}

col_fun_degs   <- colorRamp2(c(-2, 0, 2), c("green", "black", "red"))
col_fun_dels   <- colorRamp2(c(-2, 0, 2), c("blue", "black", "red"))
col_fun_hybrid <- colorRamp2(c(-2, 0, 2), c("green", "black", "yellow"))

active_legends <- list()

if ("DEGs" %in% types_subset && sum(All$Transcripts == "DEGs") > 0) {
  active_legends$degs <- ComplexHeatmap::Legend(title = "DEGs Z-score", col_fun = col_fun_degs, direction = "horizontal", labels_gp = gpar(fontsize = 10), title_gp = gpar(fontsize = 12), grid_width = unit(2, "cm"), grid_height = unit(4, "mm"), background = "white")
}
if ("DELs" %in% types_subset && sum(All$Transcripts == "DELs") > 0) {
  active_legends$dels <- ComplexHeatmap::Legend(title = "DELs Z-score", col_fun = col_fun_dels, direction = "horizontal", labels_gp = gpar(fontsize = 10), title_gp = gpar(fontsize = 12), grid_width = unit(2, "cm"), grid_height = unit(4, "mm"), background = "white")
}
if ("Hybrid" %in% types_subset && sum(All$Transcripts == "Hybrid") > 0) {
  active_legends$hybrid <- ComplexHeatmap::Legend(title = "Hybrid Z-score", col_fun = col_fun_hybrid, direction = "horizontal", labels_gp = gpar(fontsize = 10), title_gp = gpar(fontsize = 12), grid_width = unit(2, "cm"), grid_height = unit(4, "mm"), background = "white")
}

if (nrow(All) > 0) {
  active_legends$kola <- Legend(at = c("Downregulated", "Upregulated"), type = "points", legend_gp = gpar(col = c("blue", "red")), title = "log2FC", labels_gp = gpar(fontsize = 12), title_gp = gpar(fontsize = 16), grid_width = unit(0.5, "cm"), grid_height = unit(5, "mm"), background = "white")
}

pos_link_labels <- c("DEGs - Hybrid", "DELs - Hybrid", "DEGs - DELs")[c(has_pos_degs_hybrid, has_pos_dels_hybrid, has_pos_degs_dels)]
pos_link_colors <- c("cyan4", "blue", "springgreen3")[c(has_pos_degs_hybrid, has_pos_dels_hybrid, has_pos_degs_dels)]
if (length(pos_link_labels) > 0) {
  active_legends$links_pos <- Legend(at = pos_link_labels, type = "lines", legend_gp = gpar(col = pos_link_colors, lwd = 2), title = "Pos. Links", labels_gp = gpar(fontsize = 12), title_gp = gpar(fontsize = 16), grid_width = unit(0.5, "cm"), grid_height = unit(5, "mm"), background = "white")
}

neg_link_labels <- c("DEGs - Hybrid", "DELs - Hybrid", "DEGs - DELs")[c(has_neg_degs_hybrid, has_neg_dels_hybrid, has_neg_degs_dels)]
neg_link_colors <- c("purple", "orange3", "brown")[c(has_neg_degs_hybrid, has_neg_dels_hybrid, has_neg_degs_dels)]
if (length(neg_link_labels) > 0) {
  active_legends$links_neg <- Legend(at = neg_link_labels, type = "lines", legend_gp = gpar(col = neg_link_colors, lwd = 2), title = "Neg. Links", labels_gp = gpar(fontsize = 12), title_gp = gpar(fontsize = 16), grid_width = unit(0.5, "cm"), grid_height = unit(5, "mm"), background = "white")
}

upViewport()

if (length(active_legends) > 0) {
  lgd_list = do.call(packLegend, c(active_legends, list(direction = "horizontal", gap = unit(5, "mm"))))
  draw(lgd_list, y = unit(1, "npc") - circle_size, just = "top")
}

circos.clear()
dev.off()

cat("--- Circos plot generated successfully with reference legend style ---\n")
