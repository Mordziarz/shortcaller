#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(purrr)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(scales)
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
  library(gridBase)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript snv_postprocessing.R <variants_matrix.tsv> <samples_table.csv>")
}

matrix_file <- args[1]
samples_file <- args[2]

samples_info <- fread(samples_file)
if (!all(c("sample", "group") %in% names(samples_info))) {
  stop("CSV must contain 'sample' and 'group' columns.")
}

group_map <- setNames(samples_info$group, samples_info$sample)
unique_groups <- unique(samples_info$group)

if (length(unique_groups) != 2) {
  stop("Exactly two groups are required for comparison.")
}

group_name_1 <- unique_groups[1]
group_name_2 <- unique_groups[2]
message(sprintf("Group 1: %s | Group 2: %s", group_name_1, group_name_2))

dt <- fread(matrix_file, sep = "\t", header = TRUE)

dt <- dt[nchar(REF) == 1 & nchar(ALT) == 1]
dt[, ANN_clean := gsub('["\']', '', ANN)]
dt[, ANN_clean := sub(",.*", "", ANN_clean)]

dt[, c("Ann_Allele", "Effect", "Impact", "Gene_Name", "Gene_ID", 
       "Feature_Type", "transcript_id", "BioType", "Rank", "HGVS_c", "HGVS_p") := 
     tstrsplit(ANN_clean, "\\|", keep = 1:11, fill = "")]

dt[Effect == "intergenic_region", c("Gene_Name", "Gene_ID", "transcript_id") := "intergenic"]

dt[Effect != "intergenic_region", Gene_Name := sub(".*(CL\\.[0-9]+).*", "\\1", Gene_ID)]
dt[Effect != "intergenic_region", Gene_ID   := sub(".*(CL\\.[0-9]+).*", "\\1", Gene_ID)]
dt[Effect != "intergenic_region", transcript_id := sub("\\.p[0-9]+$", "", transcript_id)]

dt[, c("ANN", "ANN_clean") := NULL]

id_cols <- c("CHROM", "POS", "REF", "ALT")
annotation_cols <- c("Ann_Allele", "Effect", "Impact", "Gene_Name", "Gene_ID", 
                     "Feature_Type", "transcript_id", "BioType", "Rank", "HGVS_c", "HGVS_p")

dt_annotations <- unique(dt[, .SD, .SDcols = c(id_cols, annotation_cols)])
sample_cols <- setdiff(names(dt), c(id_cols, annotation_cols))

dt_long <- melt(
  dt, 
  id.vars = c(id_cols, annotation_cols), 
  measure.vars = sample_cols,
  variable.name = "Sample_Metric", 
  value.name = "Value"
)

dt_long[, c("Sample", "Metric") := tstrsplit(Sample_Metric, "_(?=[^_]+$)", perl = TRUE)]
dt_long[, Sample_Metric := NULL]

dt_wide_samples <- dcast(dt_long, CHROM + POS + REF + ALT + Sample + Ann_Allele + Effect + Impact + Gene_Name + Gene_ID + Feature_Type + transcript_id + BioType + Rank + HGVS_c + HGVS_p ~ Metric, value.var = "Value")


dt_wide_samples[, `:=` (
  DP = suppressWarnings(as.numeric(ifelse(DP == ".", 0, DP))),
  RO = suppressWarnings(as.numeric(ifelse(RO == ".", 0, RO))),
  AO = suppressWarnings(as.numeric(ifelse(AO == ".", 0, AO))),
  Group = group_map[Sample]
)]

dt_wide_samples[, `:=` (
  RF = ifelse(DP == 0, 0, RO / DP),
  AF = ifelse(DP == 0, 0, AO / DP),
  Successes = AO,
  Failures = DP - AO
)]

dt_model_prep <- dt_wide_samples[DP > 0]

#############################################
#############################################
dt_model_prep[, has_AF := AF > 0]

variant_group_stats <- dt_model_prep[, .(
    prop_AF_pos = sum(has_AF) / .N
), by = .(CHROM, POS, REF, ALT, Group)]


stats_wide <- dcast(variant_group_stats, CHROM + POS + REF + ALT ~ Group, value.var = "prop_AF_pos", fill = 0)

stats_wide[, keep := (get(group_name_1) >= 0.5) | (get(group_name_2) >= 0.5)]

dt_model_prep <- dt_model_prep[
    stats_wide[keep == TRUE], 
    on = id_cols, 
    nomatch = NULL
]

dt_model_prep[, has_AF := NULL]
#####################################################################
#####################################################################

variant_list_chisq <- split(dt_model_prep[Group %in% c(group_name_1, group_name_2)], by = id_cols)

valid_variants_chisq <- Filter(function(sub_df) {
  uniqueN(sub_df$Group) == 2 && all(table(sub_df$Group) >= 1)
}, variant_list_chisq)

message(sprintf("Number of variants qualified for the Chi-square test: %d", length(valid_variants_chisq)))

p_values <- map_dbl(valid_variants_chisq, function(sub_df) {
  agg <- sub_df[, .(Successes = sum(Successes), Failures = sum(Failures)), by = Group]
  
  if (nrow(agg) < 2) return(1.0)
  
  mat <- matrix(
    c(agg[Group == group_name_1, Successes], agg[Group == group_name_1, Failures],
      agg[Group == group_name_2, Successes], agg[Group == group_name_2, Failures]),
    nrow = 2, byrow = TRUE
  )
  
  if (any(rowSums(mat) == 0) || any(colSums(mat) == 0)) return(1.0)
  
  tryCatch({ chisq.test(mat)$p.value }, error = function(e) { 1.0 })
})

#results_chisq <- data.table(choosen_variant = names(p_values), p_value = p_values)
#results_chisq[, c("CHROM", "POS", "REF", "ALT") := tstrsplit(choosen_variant, "\\.")]
#results_chisq[, POS := as.integer(POS)]
#results_chisq[, choosen_variant := NULL]
#results_chisq[, p_adj := p.adjust(p_value, method = "BH")]

chisq_list_results <- lapply(names(p_values), function(variant_key) {
  sub_df <- valid_variants_chisq[[variant_key]]
  data.table(
    CHROM = sub_df$CHROM[1],
    POS = as.integer(sub_df$POS[1]),
    REF = sub_df$REF[1],
    ALT = sub_df$ALT[1],
    p_value = p_values[[variant_key]]
  )
})

results_chisq <- rbindlist(chisq_list_results)
results_chisq[, p_adj := p.adjust(p_value, method = "BH")]

dt_wide_all_metrics <- dcast(
  dt_wide_samples, 
  CHROM + POS + REF + ALT ~ Sample, 
  value.var = c("GT", "DP", "RO", "AO", "RF", "AF")
)

results_with_annot <- merge(results_chisq, dt_annotations, by = id_cols, all.x = TRUE)
final_output <- merge(results_with_annot, dt_wide_all_metrics, by = id_cols, all.x = TRUE)


samples_g1 <- names(group_map)[group_map == group_name_1]
samples_g2 <- names(group_map)[group_map == group_name_2]

af_g1_cols <- intersect(paste0("AF_", samples_g1), names(final_output))
af_g2_cols <- intersect(paste0("AF_", samples_g2), names(final_output))
rf_g1_cols <- intersect(paste0("RF_", samples_g1), names(final_output))
rf_g2_cols <- intersect(paste0("RF_", samples_g2), names(final_output))

final_output[, mean_AF_G1 := rowMeans(.SD, na.rm = TRUE), .SDcols = af_g1_cols]
final_output[, mean_AF_G2 := rowMeans(.SD, na.rm = TRUE), .SDcols = af_g2_cols]
final_output[, mean_RF_G1 := rowMeans(.SD, na.rm = TRUE), .SDcols = rf_g1_cols]
final_output[, mean_RF_G2 := rowMeans(.SD, na.rm = TRUE), .SDcols = rf_g2_cols]

final_output[, delta_AF := mean_AF_G2 - mean_AF_G1]
final_output[, delta_RF := mean_RF_G2 - mean_RF_G1]

final_output[, c("mean_AF_G1", "mean_AF_G2", "mean_RF_G1", "mean_RF_G2") := NULL]


final_output[, Simple_Effect := fcase(
  grepl("missense_variant", Effect), "Missense",
  grepl("5_prime_UTR_variant", Effect), "5' UTR",
  grepl("3_prime_UTR_variant", Effect), "3' UTR",
  grepl("synonymous_variant", Effect), "Synonymous",
  default = "Other"
)]

setorder(final_output, p_adj)
write.csv2(final_output, "All_snv.csv", row.names = FALSE)


final_output$delta_AF <- as.numeric(final_output$delta_AF)
final_output$p_adj <- as.numeric(final_output$p_adj)

final_output$AltAllelleFracDifference <- "Not significant"
final_output$AltAllelleFracDifference[final_output$delta_AF > 0.1 & final_output$p_adj < 0.05] <- "Higher"
final_output$AltAllelleFracDifference[final_output$delta_AF < -0.1 & final_output$p_adj < 0.05] <- "Lower"

final_output$Plot_Effect <- ifelse(
  final_output$AltAllelleFracDifference == "Not significant", 
  "Not significant", 
  final_output$Simple_Effect
)

effect_colors <- c(
  "Missense" = "blue4", 
  "5' UTR" = "cyan4", 
  "3' UTR" = "yellow4", 
  "Synonymous" = "palevioletred1", 
  "Other" = "magenta3",
  "Not significant" = "grey25"
)

effect_fills <- c(
  "Missense" = "blue1", 
  "5' UTR" = "cyan1", 
  "3' UTR" = "yellow1", 
  "Synonymous" = "palevioletred2", 
  "Other" = "magenta1",
  "Not significant" = "grey50"
)

# --- Barplot 1 ---
barplot_all_data <- final_output %>% count(Simple_Effect)

p_bar_all <- ggplot(barplot_all_data, aes(x = Simple_Effect, y = n, fill = Simple_Effect)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = effect_colors) +
  theme_bw(base_size = 14) +
  labs(title = "Total Detected Single Nucleotide Variant", x = "Effect Type", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none")

ggsave("single_nucleotide_variant_total_barplot.png", plot = p_bar_all, width = 8, height = 6, dpi = 300)

# --- Barplot 2 ---
final_output$Significance <- ifelse(final_output$AltAllelleFracDifference == "Not significant", "Not Significant", "Significant")
barplot_data <- final_output %>% count(Simple_Effect, Significance)

p_bar <- ggplot(barplot_data, aes(x = Simple_Effect, y = n, fill = Significance)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("Not Significant" = "grey70", "Significant" = "firebrick")) +
  theme_minimal(base_size = 14) +
  labs(title = "Detected Single Nucleotide Variant", x = "Effect Type", y = "Count", fill = "Status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), panel.grid.minor = element_blank())

ggsave("single_nucleotide_variant_barplot.png", plot = p_bar, width = 8, height = 6, dpi = 300)

# --- Volcano Plot ---
final_output$Name <- NA 

final_output$Plot_Effect <- factor(
  final_output$Plot_Effect, 
  levels = c("Missense", "5' UTR", "3' UTR", "Synonymous", "Other", "Not significant")
)

volcano_snv <- ggplot(final_output, aes(x = delta_AF, y = -log10(p_adj), color = Plot_Effect, shape = AltAllelleFracDifference, fill = Plot_Effect)) +
  geom_point(alpha = 0.5) + 
  geom_label_repel(show.legend = FALSE, aes(label = ifelse(!is.na(Name) & Name != "", Name, ""))) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dotted", linewidth = 0.5) +
  geom_hline(yintercept = 1.30103, linetype = "dotted", linewidth = 0.4) +
  xlab(expression("Alternative Allele Difference")) + 
  ylab("-log10(padj)") +
  scale_color_manual("Effect",values = effect_colors) +
  scale_fill_manual("Effect",values = effect_fills) +
  scale_shape_manual("Alternative Allele Difference",values = c("Higher" = 21, "Lower" = 22, "Not significant" = 24)) +
  scale_x_continuous(breaks = pretty_breaks(n = 5), limits = c(-1, 1)) +
  theme_bw() +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20, angle = 90),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text.align = 0,
    legend.title.align = 0,
    axis.title = element_text(size = rel(4.0))
  ) +
  guides(
    fill = guide_legend(reverse = FALSE),
    color = guide_legend(override.aes = list(shape = 15, size = 4)),
    shape = guide_legend(override.aes = list(size = 4))
  )

ggsave("volcano_single_nucleotide_variant.png", plot = volcano_snv, width = 7, height = 6, units = "in", dpi = 300)

SNV_sig <- final_output[final_output$p_adj < 0.05 & abs(final_output$delta_AF) > 0.1, ]

message("Starting Circos plot generation...")

all_possible_types <- c("Missense", "5' UTR", "3' UTR", "Synonymous", "Other")
counts_per_type <- sapply(all_possible_types, function(t) sum(SNV_sig$Simple_Effect == t))
valid_types <- all_possible_types[counts_per_type >= 2]

output_filename <- "SNV_significant_results.csv"

if (nrow(SNV_sig) == 0) {
  message("Not significant SNV events found (SNV_sig is empty). Aborting further script execution.")
  write.csv2(SNV_sig, output_filename, row.names = FALSE)
} else {
  
  g1_cols <- af_g1_cols
  g2_cols <- af_g2_cols
  
  write.csv2(SNV_sig, output_filename, row.names = FALSE)
  
  if (length(valid_types) < length(all_possible_types)) {
    message("Not every SNV type has at least 2 significant events. Skipping the circlize circular plot, but saving the CSV file.")
  } else {
    
    links_map <- list(
      list(t1="Missense", t2="5' UTR", col="blue4"),
      list(t1="Missense", t2="3' UTR", col="red"),
      list(t1="Missense", t2="Synonymous", col="palevioletred1"),
      list(t1="Missense", t2="Other", col="yellow"),
      list(t1="5' UTR", t2="3' UTR", col="magenta4"),
      list(t1="5' UTR", t2="Synonymous", col="cyan1"),
      list(t1="5' UTR", t2="Other", col="darkgreen"),
      list(t1="3' UTR", t2="Synonymous", col="yellow4"),
      list(t1="3' UTR", t2="Other", col="orange"),
      list(t1="Synonymous", t2="Other", col="black")
    )
    
    types_subset <- intersect(all_possible_types, unique(SNV_sig$Simple_Effect))
    
    SNV_sig$Simple_Effect <- factor(SNV_sig$Simple_Effect, levels = types_subset)
    All <- SNV_sig[order(SNV_sig$Simple_Effect), ]
    All$Numer <- 1:nrow(All)
    row.names(All) <- All$Numer
    
    all_sample_cols <- c(g1_cols, g2_cols)
    
    Heatmap_circlize <- as.matrix(All[, ..all_sample_cols])
    rownames(Heatmap_circlize) <- All$Numer
    
    if ("5' UTR" %in% types_subset) {
      Heatmap_circlize[All$Simple_Effect == "5' UTR", ] <- ifelse(Heatmap_circlize[All$Simple_Effect == "5' UTR", ] >= 0, Heatmap_circlize[All$Simple_Effect == "5' UTR", ] + 200, Heatmap_circlize[All$Simple_Effect == "5' UTR", ] - 200)
    }
    if ("3' UTR" %in% types_subset) {
      Heatmap_circlize[All$Simple_Effect == "3' UTR", ] <- ifelse(Heatmap_circlize[All$Simple_Effect == "3' UTR", ] >= 0, Heatmap_circlize[All$Simple_Effect == "3' UTR", ] + 300, Heatmap_circlize[All$Simple_Effect == "3' UTR", ] - 300)
    }
    if ("Synonymous" %in% types_subset) {
      Heatmap_circlize[All$Simple_Effect == "Synonymous", ] <- ifelse(Heatmap_circlize[All$Simple_Effect == "Synonymous", ] >= 0, Heatmap_circlize[All$Simple_Effect == "Synonymous", ] + 400, Heatmap_circlize[All$Simple_Effect == "Synonymous", ] - 400)
    }
    if ("Other" %in% types_subset) {
      Heatmap_circlize[All$Simple_Effect == "Other", ] <- ifelse(Heatmap_circlize[All$Simple_Effect == "Other", ] >= 0, Heatmap_circlize[All$Simple_Effect == "Other", ] + 500, Heatmap_circlize[All$Simple_Effect == "Other", ] - 500)
    }
    
    row_mean <- as.numeric(All$delta_AF)
    names(row_mean) <- All$Numer
    
    get_links <- function(t1, t2) {
      df1 <- All[All$Simple_Effect == t1, ]
      df2 <- All[All$Simple_Effect == t2, ]
      common_genes <- intersect(df1$Gene_ID, df2$Gene_ID)
      common_genes <- common_genes[!is.na(common_genes) & common_genes != "intergenic" & common_genes != ""]
      if (length(common_genes) == 0) return(data.frame())
      
      link_df <- data.frame()
      for (g in common_genes) {
        idx1 <- df1$Numer[df1$Gene_ID == g]
        idx2 <- df2$Numer[df2$Gene_ID == g]
        expand <- expand.grid(Numerkolo.x = idx1, Numerkolo.y = idx2)
        link_df <- rbind(link_df, expand)
      }
      return(link_df)
    }
    
    links_table_list <- lapply(seq_along(links_map), function(i) {
      lm <- links_map[[i]]
      if (lm$t1 %in% types_subset && lm$t2 %in% types_subset) {
        df_link <- get_links(lm$t1, lm$t2)
        if (nrow(df_link) > 0) {
          df_link$Link_Type <- paste("Intersected", lm$t1, "&", lm$t2)
          df_link$Color <- lm$col
          return(df_link)
        }
      }
      return(NULL)
    })
    Links_All <- do.call(rbind, links_table_list[!sapply(links_table_list, is.null)])
    
    split_sizes <- sapply(types_subset, function(t) sum(All$Simple_Effect == t))
    split <- factor(rep(types_subset, split_sizes), levels = types_subset)
    
    col_fun6 = colorRamp2(
      c(0, 0.2, 0.4, 0.6, 0.8, 1,
        200, 200.2, 200.4, 200.6, 200.8, 201,
        300, 300.2, 300.4, 300.6, 300.8, 301,
        400, 400.2, 400.4, 400.6, 400.8, 401,
        500, 500.2, 500.4, 500.6, 500.8, 501),
      c("#e5f5e0","#c7e9c0","#a1d99b","#74c476","#31a354","#006d2c", 
        "#eff3ff","#c6dbef","#9ecae1","#6baed6","#3182bd","#08519c", 
        "#f1eef6","#d4b9da","#c994c7","#df65b0","#dd1c77","#980043", 
        "#f2f0f7","#dadaeb","#cbc9e2","#9e9ac8","#756bb1","#54278f", 
        "#ffffb2","#fed976","#feb24c","#fd8d3c","#f03b20","#bd0026")  
    )
    
    png("SNV_significant.png", width = 14.7, height = 17, units = 'in', res = 600)
    plot.new()
    circle_size = unit(1, "snpc") 
    pushViewport(viewport(x = 0.5, y = 1, width = circle_size, height = circle_size, just = c("center", "top")))
    par(omi = gridOMI(), new = TRUE)
    circos.par(gap.after = rep(8, length(types_subset)), track.margin = c(0.02, 0.01))
    
    col_order_names <- c(g2_cols, g1_cols)
    circos.heatmap(Heatmap_circlize[, col_order_names], split = split, col = col_fun6, dend.side = NULL, rownames.side = NULL, track.height = 0.3, cluster = TRUE)
    
    circos.track(ylim = range(row_mean, na.rm = TRUE), track.height = 0.15, panel.fun = function(x, y) {
      current_subset <- CELL_META$subset
      if (length(current_subset) > 0) {
        y_sub = row_mean[current_subset]
        if (!is.null(CELL_META$row_order)) {
          y_sub = y_sub[CELL_META$row_order]
        }
        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "grey")
        circos.points(pch = 16, x = seq_along(y_sub) - 0.5, y = y_sub, col = ifelse(y_sub > 0, "red", "blue"))
      }
    }, cell.padding = c(0.02, 0, 0.02, 0))
    
    n1 <- length(g1_cols)
    n2 <- length(g2_cols)
    max_height <- n2 + n1
    
    for (i in seq_along(types_subset)) {
      circos.track(track.index = 1, panel.fun = function(x, y) {
        if(CELL_META$sector.numeric.index == i) { 
          circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), 0, 
                      CELL_META$cell.xlim[2] + convert_x(6.5, "mm"), n1, col = "white", border = "black")
          circos.text(CELL_META$cell.xlim[2] + convert_x(3.5, "mm"), n1 / 2, group_name_1, cex = 1, facing = "clockwise")
          
          circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), n1, 
                      CELL_META$cell.xlim[2] + convert_x(6.5, "mm"), max_height, col = "white", border = "black")
          circos.text(CELL_META$cell.xlim[2] + convert_x(3.5, "mm"), n1 + (n2 / 2), group_name_2, cex = 1, facing = "clockwise")
        }
      }, bg.border = NA)
    }
    
    if (!is.null(Links_All) && nrow(Links_All) > 0) {
      for (lm in links_map) {
        if (lm$t1 %in% types_subset && lm$t2 %in% types_subset) {
          df_link <- Links_All[Links_All$Link_Type == paste("Intersected", lm$t1, "&", lm$t2), ]
          if (nrow(df_link) > 0) {
            for (i in seq_len(nrow(df_link))) {
              circos.heatmap.link(df_link$Numerkolo.x[i], df_link$Numerkolo.y[i], col = lm$col)
            }
          }
        }
      }
    }
    
    col_fun1 = colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1), c("#ffffb2","#fed976","#feb24c","#fd8d3c","#f03b20","#bd0026")) 
    col_fun2 = colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1), c("#e5f5e0","#c7e9c0","#a1d99b","#74c476","#31a354","#006d2c")) 
    col_fun3 = colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1), c("#eff3ff","#c6dbef","#9ecae1","#6baed6","#3182bd","#08519c")) 
    col_fun4 = colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1), c("#f1eef6","#d4b9da","#c994c7","#df65b0","#dd1c77","#980043")) 
    col_fun5 = colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1), c("#f2f0f7","#dadaeb","#cbc9e2","#9e9ac8","#756bb1","#54278f")) 
    
    active_legends <- list()
    
    if ("Other" %in% All$Simple_Effect && sum(All$Simple_Effect == "Other") > 0) {
      active_legends$Other <- ComplexHeatmap::Legend(title = "AF Other", col_fun = col_fun1, labels = c(0,0.2,0.4,0.6,0.8,1), at = c(0,0.2,0.4,0.6,0.8,1), direction = "horizontal", labels_gp = gpar(fontsize=12), title_gp = gpar(fontsize=16), grid_width = unit(3, "cm"), grid_height = unit(5,"mm"))
    }
    if ("Missense" %in% All$Simple_Effect && sum(All$Simple_Effect == "Missense") > 0) {
      active_legends$Missense <- ComplexHeatmap::Legend(title = "AF Missense", col_fun = col_fun2, labels = c(0,0.2,0.4,0.6,0.8,1), at = c(0,0.2,0.4,0.6,0.8,1), direction = "horizontal", labels_gp = gpar(fontsize=12), title_gp = gpar(fontsize=16), grid_width = unit(3, "cm"), grid_height = unit(5,"mm"))
    }
    if ("5' UTR" %in% All$Simple_Effect && sum(All$Simple_Effect == "5' UTR") > 0) {
      active_legends$UTR5 <- ComplexHeatmap::Legend(title = "AF 5' UTR", col_fun = col_fun3, labels = c(0,0.2,0.4,0.6,0.8,1), at = c(0,0.2,0.4,0.6,0.8,1), direction = "horizontal", labels_gp = gpar(fontsize=12), title_gp = gpar(fontsize=16), grid_width = unit(3, "cm"), grid_height = unit(5,"mm"))
    }
    if ("3' UTR" %in% All$Simple_Effect && sum(All$Simple_Effect == "3' UTR") > 0) {
      active_legends$UTR3 <- ComplexHeatmap::Legend(title = "AF 3' UTR", col_fun = col_fun4, labels = c(0,0.2,0.4,0.6,0.8,1), at = c(0,0.2,0.4,0.6,0.8,1), direction = "horizontal", labels_gp = gpar(fontsize=12), title_gp = gpar(fontsize=16), grid_width = unit(3, "cm"), grid_height = unit(5,"mm"))
    }
    if ("Synonymous" %in% All$Simple_Effect && sum(All$Simple_Effect == "Synonymous") > 0) {
      active_legends$Synonymous <- ComplexHeatmap::Legend(title = "AF Synonymous", col_fun = col_fun5, labels = c(0,0.2,0.4,0.6,0.8,1), at = c(0,0.2,0.4,0.6,0.8,1), direction = "horizontal", labels_gp = gpar(fontsize=12), title_gp = gpar(fontsize=16), grid_width = unit(3, "cm"), grid_height = unit(5,"mm"))
    }
    
    if (nrow(All) > 0) {
      active_legends$kola <- Legend(at = c("Lower", "Higher"), type = "points", legend_gp = gpar(col = c("blue","red")), title = "AF Difference", labels_gp = gpar(fontsize=12), title_gp = gpar(fontsize=16), grid_width = unit(0.5, "cm"), grid_height = unit(5,"mm"), background = "white")
    }
    
    all_link_labels <- c("Intersected Missense & 5' UTR", "Intersected Missense & 3' UTR", "Intersected Missense & Synonymous", "Intersected Missense & Other",
                         "Intersected 5' UTR & 3' UTR", "Intersected 5' UTR & Synonymous", "Intersected 5' UTR & Other",
                         "Intersected 3' UTR & Synonymous", "Intersected 3' UTR & Other", "Intersected Synonymous & Other")
    all_link_colors <- c("blue4", "red", "palevioletred1", "yellow", "magenta4", "cyan1", "darkgreen", "yellow4", "orange", "black")
    
    if (!is.null(Links_All) && nrow(Links_All) > 0) {
      present_links <- intersect(all_link_labels, unique(Links_All$Link_Type))
      if (length(present_links) > 0) {
        matched_indices <- match(present_links, all_link_labels)
        active_legends$link <- Legend(
          at = all_link_labels[matched_indices], 
          type = "lines", 
          legend_gp = gpar(col = all_link_colors[matched_indices], lwd = 2),
          title = "Gene Intersection", 
          labels_gp = gpar(fontsize=12), title_gp = gpar(fontsize=16), 
          grid_width = unit(0.5, "cm"), grid_height = unit(5,"mm"), background = "white"
        )
      }
    }
    
    upViewport()
    if (length(active_legends) > 0) {
      lgd_list = do.call(packLegend, c(active_legends, list(direction = "horizontal", gap = unit(7, "mm"))))
      draw(lgd_list, y = unit(1, "npc") - circle_size, just = "top")
    }
    circos.clear()
    dev.off()
  }
}
