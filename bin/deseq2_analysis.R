#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
featurecounts_file <- args[1]
samples_csv_path   <- args[2]

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(gridBase)
  library(tibble)
})

cat("--- STARTING PAIRWISE DIFFERENTIAL EXPRESSION ANALYSIS (DESeq2) ---\n")

'%!in%' <- function(x,y)!('%in%'(x,y))
options(scipen = 999)

samples_df <- read.csv(samples_csv_path)
unique_groups <- unique(samples_df$group)

if(length(unique_groups) < 2) {
  stop("The samples.csv file must contain at least two different groups for a 1 vs 1 comparison!")
}

ref_group <- unique_groups[1]
test_group <- unique_groups[2]

cat("Reference group (ref):", ref_group, "\n")
cat("Test group (test):", test_group, "\n")

samples_ref <- samples_df %>% filter(group == ref_group) %>% pull(sample)
samples_test <- samples_df %>% filter(group == test_group) %>% pull(sample)
target_samples <- c(samples_ref, samples_test)

count_data <- read.delim(featurecounts_file, header = TRUE, skip = 1, stringsAsFactors = FALSE)

colnames(count_data) <- gsub("_Aligned.sortedByCoord.out.bam", "", colnames(count_data))
colnames(count_data) <- gsub("_Aligned.out.bam", "", colnames(count_data))

id_col_name <- intersect(c("Geneid", "gene_id"), colnames(count_data))[1]
if(is.na(id_col_name)) {
  stop("Could not find 'Geneid' or 'gene_id' column in the count matrix!")
}

available_samples <- intersect(target_samples, colnames(count_data))

if(length(available_samples) == 0) {
  stop("No matching sample names between samples.csv and featureCounts column headers!")
}

count_subset <- count_data[, c(id_col_name, available_samples)]
rownames(count_subset) <- count_subset[[id_col_name]]
count_subset[[id_col_name]] <- NULL

sample_names <- colnames(count_subset)
condition <- sapply(sample_names, function(s) {
  samples_df$group[samples_df$sample == s][1]
})

colData <- data.frame(condition = as.factor(condition))
rownames(colData) <- sample_names

dds <- DESeqDataSetFromMatrix(countData = as.matrix(count_subset),
                              colData = colData,
                              design = ~ condition)

dds$condition <- relevel(dds$condition, ref = ref_group)
dds <- DESeq(dds)

dds_counts <- counts(dds, normalized = TRUE)

coef_name <- paste0("condition_", test_group, "_vs_", ref_group)

resLFC <- tryCatch({
  lfcShrink(dds, coef = coef_name, type = "apeglm")
}, error = function(e) {
  res <- results(dds, contrast = c("condition", test_group, ref_group))
  lfcShrink(dds, contrast = c("condition", test_group, ref_group), res = res, type = "normal")
})

Ekspresja <- as.data.frame(resLFC)
Ekspresja <- Ekspresja[!is.na(Ekspresja$log2FoldChange), ]
Ekspresja <- merge(Ekspresja, dds_counts, by = 0)
rownames(Ekspresja) <- Ekspresja$Row.names
Ekspresja$Row.names <- NULL

genes_sig <- Ekspresja
genes_sig$padj <- as.numeric(genes_sig$padj)
genes_sig$log2FoldChange <- as.numeric(genes_sig$log2FoldChange)
genes_sig <- genes_sig[!is.na(genes_sig$padj), ]

genes_sig$Expression <- "No significant"
genes_sig$Expression <- ifelse(genes_sig$log2FoldChange < -1 & genes_sig$padj < 0.05, "Downregulated", genes_sig$Expression)
genes_sig$Expression <- ifelse(genes_sig$log2FoldChange > 1 & genes_sig$padj < 0.05, "Upregulated", genes_sig$Expression)
genes_sig$Transcripts <- "No significant"
genes_sig$Transcripts <- ifelse(abs(genes_sig$log2FoldChange) > 1 & genes_sig$padj < 0.05, "DEGs", genes_sig$Transcripts)
genes_sig$order <- ifelse(genes_sig$Transcripts %in% "No significant", 0, 1)
genes_sig$Symbol <- NA

degs_df <- genes_sig[genes_sig$Transcripts %in% "DEGs", ]

degs_df$GeneID <- rownames(degs_df)
degs_df <- degs_df[, c("GeneID", setdiff(names(degs_df), "GeneID"))]

Ekspresja$GeneID <- rownames(Ekspresja)
Ekspresja <- Ekspresja[, c("GeneID", setdiff(names(Ekspresja), "GeneID"))]

Ekspresja_save <- Ekspresja
Ekspresja_save$order <- NULL
Ekspresja_save$Symbol <- NULL

degs_df_save <- degs_df
degs_df_save$order <- NULL
degs_df_save$Symbol <- NULL

prefix_file <- paste0(test_group, "_vs_", ref_group)
write.csv2(Ekspresja_save, paste0("Expression", ".csv"), row.names = F)
write.csv2(degs_df_save, paste0("Genes_sig", ".csv"), row.names = F)

cat("--- DIFFERENTIAL EXPRESSION ANALYSIS COMPLETED SUCCESSFULLY ---\n")

cat("--- Generating PCA plot ---\n")
vsd <- tryCatch({
  vst(dds, blind = FALSE)
}, error = function(e) {
  message("Note: Too few genes for default VST. Using varianceStabilizingTransformation() instead.")
  varianceStabilizingTransformation(dds, blind = FALSE)
})

pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 4, show.legend = FALSE) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 14) +
  ggtitle("PCA Plot of Samples") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("PCA_plot.png", plot = pca_p, width = 8, height = 6, dpi = 300)


cat("--- Generating Volcano plot ---\n")
Final_ekspresja <- genes_sig
Final_ekspresja$project <- prefix_file

options(ggrepel.max.overlaps = Inf)
volcano_p <- ggplot(Final_ekspresja, aes(x = log2FoldChange, y = -log10(padj), color = Expression, shape = Expression, fill = Expression)) +
  geom_point(alpha = 0.6, size = 2) + 
  geom_text_repel(show.legend = FALSE, aes(label = Symbol), force = 2, nudge_y = 1, nudge_x = -0.25, color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", linewidth = 0.5, color = "red") +
  geom_hline(yintercept = 1.30103, linetype = "dotted", linewidth = 0.4, color = "red") +
  xlab(expression("log2FoldChange")) + ylab("-log10(padj)") +
  scale_color_manual(values = c("Upregulated" = "firebrick", "Downregulated" = "dodgerblue", "No significant" = "grey25")) +
  scale_shape_manual(values = c("Upregulated" = 21, "Downregulated" = 21, "No significant" = 21)) +
  scale_fill_manual(values = c("Upregulated" = "firebrick1", "Downregulated" = "dodgerblue1", "No significant" = "grey50")) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "right"
  ) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 4)))

ggsave("Volcano_DEGs.png", plot = volcano_p, width = 10, height = 6, dpi = 300)

cat("--- Generating MA plot ---\n")
MA_plot <- Final_ekspresja

MA_plot$Means <- log2(rowMeans(sapply(MA_plot[, sample_names, drop = FALSE], as.numeric), na.rm = TRUE) + 1)

maplot_p <- ggplot(MA_plot, aes(x = Means, y = log2FoldChange, color = Expression, shape = Expression, fill = Expression)) +
  geom_point(alpha = 0.6, size = 2) + 
  xlab(expression("log2(Mean of normalized counts)")) + ylab("log2FoldChange") +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted", linewidth = 0.4, color = "red") +
  scale_color_manual(values = c("Upregulated" = "firebrick", "Downregulated" = "dodgerblue", "No significant" = "grey25")) +
  scale_shape_manual(values = c("Upregulated" = 21, "Downregulated" = 21, "No significant" = 21)) +
  scale_fill_manual(values = c("Upregulated" = "firebrick1", "Downregulated" = "dodgerblue1", "No significant" = "grey50")) +
  geom_text_repel(show.legend = FALSE, aes(label = Symbol), force = 2, nudge_y = 1, nudge_x = -0.25, color = "black") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "right"
  ) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 4)))

ggsave("MA_plot_DEGs.png", plot = maplot_p, width = 10, height = 6, dpi = 300)


cat("--- Generating standard Heatmap ---\n")
sig_matrix_data <- as.matrix(count_subset[rownames(degs_df), sample_names, drop = FALSE])

if(nrow(sig_matrix_data) > 1) {
  sig_matrix_scaled <- t(scale(t(log2(sig_matrix_data + 1))))
  
  cond_colors <- c("steelblue", "tomato")
  names(cond_colors) <- c(ref_group, test_group)
  
  ha <- HeatmapAnnotation(
    Condition = colData$condition,
    col = list(Condition = cond_colors),
    show_legend = c(Condition = TRUE),
    annotation_legend_param = list(
      Condition = list(title = "Condition", direction = "vertical")
    ),
    simple_anno_size = unit(0.4, "cm")
  )
  
  png("Heatmap_DEGs.png", width = 9, height = 8, units = "in", res = 300)
  
  ht <- Heatmap(sig_matrix_scaled, 
                name = "Z-score", 
                show_row_names = FALSE,
                column_title = paste0("Differentially Expressed Genes Heatmap"),
                top_annotation = ha,
                heatmap_legend_param = list(direction = "vertical"),
                heatmap_width = unit(6, "inch"))
  
  draw(ht, heatmap_legend_side = "right", merge_legends = TRUE, padding = unit(c(10, 15, 10, 10), "mm"))
  
  dev.off()
} else {
  message("Note: Not enough significant genes to generate a standard heatmap.")
}

cat("--- Generating Top 10 DEGs Heatmap ---\n")

if(nrow(degs_df) > 0) {
  top10_degs <- degs_df[order(degs_df$padj), ]
  top10_n <- min(10, nrow(top10_degs))
  top10_degs <- top10_degs[1:top10_n, ]
  
  top10_matrix_data <- as.matrix(count_subset[top10_degs$GeneID, sample_names, drop = FALSE])
  
  if(nrow(top10_matrix_data) > 1) {
    top10_matrix_scaled <- t(scale(t(log2(top10_matrix_data + 1))))
    
    rownames(top10_matrix_scaled) <- top10_degs$GeneID
    
    cond_colors <- c("steelblue", "tomato")
    names(cond_colors) <- c(ref_group, test_group)
    
    ha_top10 <- HeatmapAnnotation(
      Condition = colData$condition,
      col = list(Condition = cond_colors),
      show_legend = c(Condition = TRUE),
      annotation_legend_param = list(
        Condition = list(title = "Condition", direction = "vertical")
      ),
      simple_anno_size = unit(0.4, "cm")
    )
    
    png("Heatmap_Top10_DEGs.png", width = 9, height = 8, units = "in", res = 300)
    
    ht_top10 <- Heatmap(top10_matrix_scaled, 
                        name = "Z-score", 
                        show_row_names = TRUE,
                        row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                        column_title = paste0("Top ", top10_n, " DEGs Heatmap"),
                        top_annotation = ha_top10,
                        heatmap_legend_param = list(direction = "vertical"),
                        heatmap_width = unit(6, "inch"))
    
    draw(ht_top10, heatmap_legend_side = "right", merge_legends = TRUE, padding = unit(c(10, 15, 10, 10), "mm"))
    
    dev.off()
  } else {
    message("Note: Not enough genes to generate Top 10 heatmap.")
  }
} else {
  message("Note: No significant DEGs found to generate Top 10 heatmap.")
}

cat("--- Generating automated Circos plot for all DEGs ---\n")

n_ref <- length(samples_ref)
n_test <- length(samples_test)

if (nrow(degs_df) == 0) {
  message("Not significant differentially expressed genes found (degs_df is empty). Aborting Circos plot generation.")
} else {
  
  All <- degs_df[order(degs_df$log2FoldChange), ]
  All$Numer <- 1:nrow(All)
  row.names(All) <- All$Numer
  
  target_samples <- c(samples_ref, samples_test) 
  circos_mat <- log2(as.matrix(All[, rev(target_samples), drop = FALSE]) + 1)
  
  circos_mat <- t(scale(t(circos_mat), scale = FALSE))
  circos_mat[is.na(circos_mat)] <- 100
  rownames(circos_mat) <- All$GeneID
  
  col_fun_circos = colorRamp2(c(-2, 0, 2), c("blue", "black", "red"))
  col_fun_circos2 = colorRamp2(c(-2, 0, 2, 100), c("blue", "black", "red", "white"))
  
  png("circoss_all_degs.png", width = 12, height = 14, units = 'in', res = 600)
  
  plot.new()
  circle_size = unit(1, "snpc") 
  pushViewport(viewport(x = 0.5, y = 1, width = circle_size, height = circle_size, just = c("center", "top")))
  par(omi = gridOMI(), new = TRUE)
  
  circos.par(track.margin = c(0.02, 0.01))
  
  circos.heatmap(circos_mat, 
                 col = col_fun_circos2, 
                 dend.side = NULL, 
                 rownames.side = NULL, 
                 track.height = 0.3, 
                 cluster = FALSE)
  
  row_mean_vals <- as.numeric(All$log2FoldChange)
  names(row_mean_vals) <- All$Numer
  
  circos.track(ylim = range(row_mean_vals, na.rm = TRUE), track.height = 0.15, panel.fun = function(x, y) {
    current_subset <- CELL_META$subset
    if (length(current_subset) > 0) {
      y_sub = row_mean_vals[current_subset]
      if (!is.null(CELL_META$row_order)) {
        y_sub = y_sub[CELL_META$row_order]
      }
      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "grey")
      circos.points(pch = 16, x = seq_along(y_sub) - 0.5, y = y_sub, col = ifelse(y_sub > 0, "red", "blue"))
    }
  }, cell.padding = c(0.02, 0, 0.02, 0))
  
  mid_height <- n_test 
  max_height <- n_test + n_ref
  
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), 0, 
                CELL_META$cell.xlim[2] + convert_x(6.5, "mm"), n_ref, col = "white", border = "black")
    circos.text(CELL_META$cell.xlim[2] + convert_x(3.5, "mm"), n_ref / 2, ref_group, cex = 1, facing = "clockwise")
    
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), n_ref, 
                CELL_META$cell.xlim[2] + convert_x(6.5, "mm"), max_height, col = "white", border = "black")
    circos.text(CELL_META$cell.xlim[2] + convert_x(3.5, "mm"), n_ref + (n_test / 2), test_group, cex = 1, facing = "clockwise")
  }, bg.border = NA)
  
  lgd_ekspresja = ComplexHeatmap::Legend(title = "Z-score", col_fun = col_fun_circos, direction = "horizontal", labels_gp = gpar(fontsize = 10), title_gp = gpar(fontsize = 12), grid_width = unit(2, "cm"), grid_height = unit(4, "mm"))
  lgd_kola = Legend(at = c("Downregulated", "Upregulated"), type = "points", legend_gp = gpar(col = c("blue", "red")), title = "log2FoldChange", labels_gp = gpar(fontsize = 10), title_gp = gpar(fontsize = 12), grid_width = unit(0.4, "cm"), grid_height = unit(4, "mm"), background = "white")
  
  upViewport()
  lgd_list = packLegend(lgd_ekspresja, lgd_kola, direction = "horizontal")
  draw(lgd_list, y = unit(1, "npc") - circle_size, just = "top")
  
  circos.clear()
  dev.off()
  
  cat("--- Circos plot generated successfully ---\n")
}

cat("--- Generating automated Circos plot for Top 50 DEGs ---\n")

if (nrow(degs_df) == 0) {
  message("No significant differentially expressed genes found (degs_df is empty). Aborting Circos plot generation.")
} else {
  
  top50_degs <- degs_df[order(degs_df$padj), ]
  n_top50 <- min(50, nrow(top50_degs))
  top50_degs <- top50_degs[1:n_top50, ]
  
  All <- top50_degs[order(top50_degs$log2FoldChange), ]
  All$Numer <- 1:nrow(All)
  row.names(All) <- All$Numer
  
  target_samples <- c(samples_ref, samples_test)
  circos_mat <- log2(as.matrix(All[, rev(target_samples), drop = FALSE]) + 1)
  
  circos_mat <- t(scale(t(circos_mat), scale = FALSE))
  circos_mat[is.na(circos_mat)] <- 100
  rownames(circos_mat) <- All$GeneID
  
  col_fun_circos = colorRamp2(c(-2, 0, 2), c("blue", "black", "red"))
  col_fun_circos2 = colorRamp2(c(-2, 0, 2, 100), c("blue", "black", "red", "white"))
  
  png("circoss_top50_degs.png", width = 12, height = 14, units = 'in', res = 600)
  
  plot.new()
  circle_size = unit(1, "snpc") 
  pushViewport(viewport(x = 0.5, y = 1, width = circle_size, height = circle_size, just = c("center", "top")))
  par(omi = gridOMI(), new = TRUE)
  
  circos.par(track.margin = c(0.02, 0.01))
  
  circos.heatmap(circos_mat, 
                 col = col_fun_circos2, 
                 dend.side = NULL, 
                 rownames.side = "outside", 
                 track.height = 0.3, 
                 rownames.cex = 0.8,
                 cluster = FALSE)
  
  row_mean_vals <- as.numeric(All$log2FoldChange)
  names(row_mean_vals) <- All$Numer
  
  circos.track(ylim = range(row_mean_vals, na.rm = TRUE), track.height = 0.15, panel.fun = function(x, y) {
    current_subset <- CELL_META$subset
    if (length(current_subset) > 0) {
      y_sub = row_mean_vals[current_subset]
      if (!is.null(CELL_META$row_order)) {
        y_sub = y_sub[CELL_META$row_order]
      }
      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "grey")
      circos.points(pch = 16, x = seq_along(y_sub) - 0.5, y = y_sub, col = ifelse(y_sub > 0, "red", "blue"))
    }
  }, cell.padding = c(0.02, 0, 0.02, 0))
  
  mid_height <- n_test 
  max_height <- n_test + n_ref
  
  circos.track(track.index = 2, panel.fun = function(x, y) {
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), 0, 
                CELL_META$cell.xlim[2] + convert_x(6.5, "mm"), n_ref, col = "white", border = "black")
    circos.text(CELL_META$cell.xlim[2] + convert_x(3.5, "mm"), n_ref / 2, ref_group, cex = 1, facing = "clockwise")
    
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), n_ref, 
                CELL_META$cell.xlim[2] + convert_x(6.5, "mm"), max_height, col = "white", border = "black")
    circos.text(CELL_META$cell.xlim[2] + convert_x(3.5, "mm"), n_ref + (n_test / 2), test_group, cex = 1, facing = "clockwise")
  }, bg.border = NA)
  
  lgd_ekspresja = ComplexHeatmap::Legend(title = "Z-score", col_fun = col_fun_circos, direction = "horizontal", labels_gp = gpar(fontsize = 10), title_gp = gpar(fontsize = 12), grid_width = unit(2, "cm"), grid_height = unit(4, "mm"))
  lgd_kola = Legend(at = c("Downregulated", "Upregulated"), type = "points", legend_gp = gpar(col = c("blue", "red")), title = "log2FoldChange", labels_gp = gpar(fontsize = 10), title_gp = gpar(fontsize = 12), grid_width = unit(0.4, "cm"), grid_height = unit(4, "mm"), background = "white")
  
  upViewport()
  lgd_list = packLegend(lgd_ekspresja, lgd_kola, direction = "horizontal")
  draw(lgd_list, y = unit(1, "npc") - circle_size, just = "top")
  
  circos.clear()
  dev.off()
  
  cat("--- Circos plot generated successfully ---\n")
}

cat("--- ALL DESeq2 ANALYSIS AND PLOTS COMPLETED SUCCESSFULLY ---\n")
