#!/usr/bin/env Rscript

library(ggplot2)
library(ggrepel)
library(scales)
library(stringr)
library(circlize)
library(ComplexHeatmap)
library(grid)
library(dplyr)
library(gridBase)

args <- commandArgs(trailingOnly = TRUE)
r_dir <- args[1]      
samples_csv <- args[2]

samples_info <- read.table(samples_csv, sep = ",", header = TRUE, stringsAsFactors = FALSE) 

unique_groups <- unique(samples_info$group)
group_name_1 <- unique_groups[1] 
group_name_2 <- unique_groups[2] 

event_types <- c("SE", "MXE", "A5SS", "A3SS", "RI")
data_list <- list()

for (et in event_types) {
  file_jcec <- file.path(r_dir, paste0(et, ".MATS.JCEC.txt"))
  file_jc   <- file.path(r_dir, paste0(et, ".MATS.JC.txt"))
  
  if (file.exists(file_jcec) && file.exists(file_jc) && 
      file.info(file_jcec)$size > 0 && file.info(file_jc)$size > 0) {
    
    df_jcec <- tryCatch(read.table(file_jcec, header = TRUE, stringsAsFactors = FALSE), error = function(e) NULL)
    df_jc   <- tryCatch(read.table(file_jc, header = TRUE, stringsAsFactors = FALSE), error = function(e) NULL)
    
    if (!is.null(df_jcec) && !is.null(df_jc) && nrow(df_jcec) > 0 && nrow(df_jc) > 0) {
      df_tot  <- rbind(df_jcec, df_jc)
      df_tot  <- df_tot[!duplicated(df_tot$ID), ]
      df_tot$Type <- et
      data_list[[et]] <- df_tot
    } else {
      message(paste("Nextflow info: Files for type", et, "are empty. Skipping."))
    }
  } else {
    message(paste("Nextflow info: No files found for type", et, ". Skipping."))
  }
}

event_types <- names(data_list)

for (et in event_types) {
  df <- as.data.frame(data_list[[et]])
  
  if (!is.null(df) && nrow(df) > 0) {
    
    if ("FDR" %in% colnames(df) & "IncLevelDifference" %in% colnames(df)) {
      df$FDR <- as.numeric(df$FDR)
      df$IncLevelDifference <- as.numeric(df$IncLevelDifference)
      df$IncLevelDifference <- -(df$IncLevelDifference)
      df$Type <- NULL
      
      df_sig <- df[!is.na(df$FDR) & df$FDR < 0.05 & !is.na(df$IncLevelDifference) & abs(df$IncLevelDifference) > 0.1, ]
      
      df_sig$geneSymbol <- df_sig$GeneID
      
      write.table(df_sig, paste0(et, "_sashimi.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    }
  }
}

event_types_present <- names(data_list)
processed_list <- list()

for (et in event_types_present) {
  df <- data_list[[et]]
  
  if (is.null(df) || nrow(df) == 0) {
    next
  }
  
  if (et == "RI") {
    df$index <- paste0(df$GeneID, "_", df$riExonStart_0base, "_", df$riExonEnd, "_", df$upstreamES, "_", df$upstreamEE, "_", df$downstreamES, "_", df$downstreamEE)
  } else if (et == "SE") {
  df$index <- paste0(df$GeneID, "_", df$exonStart_0base, "_", df$exonEnd, "_", df$upstreamES, "_", df$upstreamEE, "_", df$downstreamES, "_", df$downstreamEE)
  } else if (et %in% c("A3SS", "A5SS")) {
    df$index <- paste0(df$GeneID, "_", df$longExonStart_0base, "_", df$longExonEnd, "_", df$shortES, "_", df$shortEE, "_", df$flankingES, "_", df$flankingEE)
  } else if (et == "MXE") {
    df$index <- paste0(df$GeneID, "_", df$X1stExonStart_0base, "_", df$X1stExonEnd, "_", df$X2ndExonStart_0base, "_", df$X2ndExonEnd, "_", df$upstreamES, "_", df$upstreamEE, "_", df$downstreamES, "_", df$downstreamEE)
  }
  
  sub_df <- data.frame(
    GeneID             = df$GeneID,
    IncLevel1          = df$IncLevel1,
    IncLevel2          = df$IncLevel2,
    FDR                = df$FDR,
    IncLevelDifference = df$IncLevelDifference,
    ID                 = df$ID,
    index              = df$index,
    Type               = df$Type,
    stringsAsFactors   = FALSE
  )
  
  processed_list[[et]] <- sub_df
}

Wszystkie <- do.call(rbind, processed_list)

colnames(Wszystkie) <- c("GeneID", "IncLevel1", "IncLevel2", "FDR", "IncLevelDifference", "ID", "index", "Type")

Wszystkie$IncLevelDifference <- as.numeric(Wszystkie$IncLevelDifference)
Wszystkie$FDR <- as.numeric(Wszystkie$FDR)
Wszystkie$IncLevelDifference <- -(Wszystkie$IncLevelDifference)
Wszystkie$OriginalType <- Wszystkie$Type

Wszystkie$`Inclusion level` <- "Not significant"

idx_higher <- which(Wszystkie$IncLevelDifference > 0.1 & Wszystkie$FDR < 0.05)
if (length(idx_higher) > 0) {
  Wszystkie$`Inclusion level`[idx_higher] <- "Higher"
}

idx_lower <- which(Wszystkie$IncLevelDifference < -0.1 & Wszystkie$FDR < 0.05)
if (length(idx_lower) > 0) {
  Wszystkie$`Inclusion level`[idx_lower] <- "Lower"
}

Wszystkie$order <- ifelse(Wszystkie$`Inclusion level` %in% "Not significant", 0, 1)
Wszystkie$Name <- NA
Wszystkie$Type <- ifelse(abs(Wszystkie$IncLevelDifference) < 0.1 | Wszystkie$FDR > 0.05, "Not significant", Wszystkie$Type)
Wszystkie$FDR2log <- -log10(Wszystkie$FDR)

barplot_all_data <- Wszystkie %>%
  count(OriginalType)

p_bar_all <- ggplot(barplot_all_data, aes(x = OriginalType, y = n, fill = OriginalType)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("RI" = "red3", "SE" = "green3", "A3SS" = "blue3", "MXE" = "magenta3", "A5SS" = "gold3")) +
  theme_bw(base_size = 14) +
  labs(
    title = "Total Detected Alternative Splicing Events",
    x = "Event Type",
    y = "Number of Events"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none"
  )

ggsave("splicing_events_total_barplot.png", plot = p_bar_all, width = 8, height = 6, units = "in", dpi = 300)

Wszystkie$Significance <- ifelse(Wszystkie$`Inclusion level` == "Not significant", "Not Significant", "Significant")

barplot_data <- Wszystkie %>%
  count(OriginalType, Significance)

p_bar <- ggplot(barplot_data, aes(x = OriginalType, y = n, fill = Significance)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("Not Significant" = "grey70", "Significant" = "firebrick")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Detected Alternative Splicing Events",
    x = "Event Type",
    y = "Number of Events",
    fill = "Status"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("splicing_events_barplot.png", plot = p_bar, width = 8, height = 6, dpi = 300)

volcano_splicing <- ggplot(Wszystkie, aes(x = IncLevelDifference, y = FDR2log, color = Type, shape = `Inclusion level`, fill = Type)) +
  geom_point(alpha = 0.5) + 
  geom_label_repel(show.legend = FALSE, aes(label = ifelse(!is.na(Name) & Name != "", Name, ""))) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dotted", linewidth = 0.5) +
  geom_hline(yintercept = 1.30103, linetype = "dotted", linewidth = 0.4) +
  xlab(expression("InclusionLevelDifference")) + 
  ylab("-log10(FDR)") +
  scale_color_manual(values = c("RI" = "red4", "Not significant" = "grey25", "SE" = "green4", "A3SS" = "blue4", "MXE" = "magenta4", "A5SS" = "gold4")) +
  scale_shape_manual(values = c("Higher" = 21, "Lower" = 22, "Not significant" = 24)) +
  scale_fill_manual(values = c("RI" = "red1", "Not significant" = "grey50", "SE" = "green1", "A3SS" = "blue1", "MXE" = "magenta1", "A5SS" = "gold1")) +
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
    color = guide_legend(override.aes = list(shape = 15, size = 4))
  )

ggsave("volcano_splicing.png", plot = volcano_splicing, width = 6, height = 6, units = "in", dpi = 300)

Wszystkie$order <- NULL
Wszystkie$FDR2log <- NULL
Wszystkie$Name <- NULL

write.csv2(Wszystkie, "All_splicing_events.csv", row.names = FALSE)

Splic_sig <- Wszystkie[Wszystkie$FDR < 0.05 & abs(Wszystkie$IncLevelDifference) > 0.1, ]

all_possible_types <- c("RI", "SE", "A5SS", "A3SS", "MXE")
counts_per_type <- sapply(all_possible_types, function(t) sum(Splic_sig$Type == t))
valid_types <- all_possible_types[counts_per_type >= 2]

if (nrow(Splic_sig) == 0) {
  message("Not significant splicing events found (Splic_sig is empty). Aborting further script execution.")
} else {

  n_per_group_1 <- sum(samples_info$group == group_name_1)
  n_per_group_2 <- sum(samples_info$group == group_name_2)

  g1_cols <- paste0(group_name_1, "_", 1:n_per_group_1)
  g2_cols <- paste0(group_name_2, "_", 1:n_per_group_2)
  
  Splic_sig[, g1_cols] <- str_split_fixed(Splic_sig$IncLevel1, ",", n_per_group_1)
  Splic_sig[, g2_cols] <- str_split_fixed(Splic_sig$IncLevel2, ",", n_per_group_2)
  
  Splic_sig$order <- NULL
  Splic_sig$FDR2log <- NULL
    
  output_filename <- paste0("Splicing_significant_results.csv")
  Splic_sig$Name <- NULL
  write.csv2(Splic_sig, output_filename, row.names = FALSE)

  if (length(valid_types) < length(all_possible_types)) {
    message("Not every splicing type has at least 2 significant events. Skipping the circlize circular plot, but saving the CSV file with all significant results.")
    
  } else {
  
  Splic_sig$order <- NULL
  Splic_sig$FDR2log <- NULL
  output_filename <- paste0("Splicing_significant_results.csv")
  Splic_sig$Name <- NULL
  write.csv2(Splic_sig, output_filename, row.names = FALSE)
  
  links_map <- list(
  list(t1="RI", t2="SE", col="blue4"),
  list(t1="RI", t2="A5SS", col="red"),
  list(t1="RI", t2="A3SS", col="palevioletred1"),
  list(t1="RI", t2="MXE", col="yellow"),
  list(t1="SE", t2="A5SS", col="magenta4"),
  list(t1="SE", t2="A3SS", col="cyan1"),
  list(t1="SE", t2="MXE", col="darkgreen"),
  list(t1="A5SS", t2="A3SS", col="yellow4"),
  list(t1="A5SS", t2="MXE", col="orange"),
  list(t1="A3SS", t2="MXE", col="black")
  )

  all_possible_types <- c("RI", "SE", "A5SS", "A3SS", "MXE")
  types_subset <- intersect(all_possible_types, unique(Splic_sig$Type))
  
  Splic_sig$Type <- factor(Splic_sig$Type, levels = types_subset)
  All <- Splic_sig[order(Splic_sig$Type), ]
  All$Numer <- 1:nrow(All)
  row.names(All) <- All$Numer
  
  all_sample_cols <- c(g1_cols, g2_cols)
  
  Heatmap_circlize <- as.matrix(sapply(All[, all_sample_cols], as.numeric))
  rownames(Heatmap_circlize) <- All$Numer
  
  if ("SE" %in% types_subset) {
    Heatmap_circlize[All$Type == "SE",   ] <- ifelse(Heatmap_circlize[All$Type == "SE",   ] >= 0, Heatmap_circlize[All$Type == "SE",   ] + 200, Heatmap_circlize[All$Type == "SE",   ] - 200)
  }
  if ("A5SS" %in% types_subset) {
    Heatmap_circlize[All$Type == "A5SS", ] <- ifelse(Heatmap_circlize[All$Type == "A5SS", ] >= 0, Heatmap_circlize[All$Type == "A5SS", ] + 300, Heatmap_circlize[All$Type == "A5SS", ] - 300)
  }
  if ("A3SS" %in% types_subset) {
    Heatmap_circlize[All$Type == "A3SS", ] <- ifelse(Heatmap_circlize[All$Type == "A3SS", ] >= 0, Heatmap_circlize[All$Type == "A3SS", ] + 400, Heatmap_circlize[All$Type == "A3SS", ] - 400)
  }
  if ("MXE" %in% types_subset) {
    Heatmap_circlize[All$Type == "MXE",  ] <- ifelse(Heatmap_circlize[All$Type == "MXE",  ] >= 0, Heatmap_circlize[All$Type == "MXE",  ] + 500, Heatmap_circlize[All$Type == "MXE",  ] - 500)
  }
  
  row_mean <- as.numeric(All$IncLevelDifference)
  names(row_mean) <- All$Numer
  
  get_links <- function(t1, t2) {
  df1 <- All[All$Type == t1, ]
  df2 <- All[All$Type == t2, ]
  common_genes <- intersect(df1$GeneID, df2$GeneID)
  
  if (length(common_genes) == 0) return(data.frame())
  
  link_df <- data.frame()
  for (g in common_genes) {
    idx1 <- df1$Numer[df1$GeneID == g]
    idx2 <- df2$Numer[df2$GeneID == g]
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
  
  split_sizes <- sapply(types_subset, function(t) sum(All$Type == t))
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
  
  png("circlize_splicing_significant.png", width = 14, height = 16, units = 'in', res = 600)
  
  plot.new()
  circle_size = unit(1, "snpc") 
  pushViewport(viewport(x = 0.5, y = 1, width = circle_size, height = circle_size,
                        just = c("center", "top")))
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
  
  mid_height <- n2 
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
  
  if ("MXE" %in% All$Type && sum(All$Type == "MXE") > 0) {
    active_legends$MXE <- ComplexHeatmap::Legend(title = "Inc Level MXE", col_fun = col_fun1, labels = c(0,0.2,0.4,0.6,0.8,1), at = c(0,0.2,0.4,0.6,0.8,1), direction = "horizontal", labels_gp = gpar(fontsize=12), title_gp = gpar(fontsize=16), grid_width = unit(3, "cm"), grid_height = unit(5,"mm"))
  }
  if ("RI" %in% All$Type && sum(All$Type == "RI") > 0) {
    active_legends$RI <- ComplexHeatmap::Legend(title = "Inc Level RI", col_fun = col_fun2, labels = c(0,0.2,0.4,0.6,0.8,1), at = c(0,0.2,0.4,0.6,0.8,1), direction = "horizontal", labels_gp = gpar(fontsize=12), title_gp = gpar(fontsize=16), grid_width = unit(3, "cm"), grid_height = unit(5,"mm"))
  }
  if ("SE" %in% All$Type && sum(All$Type == "SE") > 0) {
    active_legends$SE <- ComplexHeatmap::Legend(title = "Inc Level SE", col_fun = col_fun3, labels = c(0,0.2,0.4,0.6,0.8,1), at = c(0,0.2,0.4,0.6,0.8,1), direction = "horizontal", labels_gp = gpar(fontsize=12), title_gp = gpar(fontsize=16), grid_width = unit(3, "cm"), grid_height = unit(5,"mm"))
  }
  if ("A5SS" %in% All$Type && sum(All$Type == "A5SS") > 0) {
    active_legends$A5SS <- ComplexHeatmap::Legend(title = "Inc Level A5SS", col_fun = col_fun4, labels = c(0,0.2,0.4,0.6,0.8,1), at = c(0,0.2,0.4,0.6,0.8,1), direction = "horizontal", labels_gp = gpar(fontsize=12), title_gp = gpar(fontsize=16), grid_width = unit(3, "cm"), grid_height = unit(5,"mm"))
  }
  if ("A3SS" %in% All$Type && sum(All$Type == "A3SS") > 0) {
    active_legends$A3SS <- ComplexHeatmap::Legend(title = "Inc Level A3SS", col_fun = col_fun5, labels = c(0,0.2,0.4,0.6,0.8,1), at = c(0,0.2,0.4,0.6,0.8,1), direction = "horizontal", labels_gp = gpar(fontsize=12), title_gp = gpar(fontsize=16), grid_width = unit(3, "cm"), grid_height = unit(5,"mm"))
  }
  
  if (nrow(All) > 0) {
    active_legends$kola <- Legend(at = c("Lower", "Higher"), type = "points", legend_gp = gpar(col = c("blue","red")), title = "IncLevelDiff", labels_gp = gpar(fontsize=12), title_gp = gpar(fontsize=16), grid_width = unit(0.5, "cm"), grid_height = unit(5,"mm"), background = "white")
  }
  
  all_link_labels <- c("Intersected RI & SE", "Intersected RI & A5SS", "Intersected RI & A3SS", "Intersected RI & MXE",
                       "Intersected SE & A5SS", "Intersected SE & A3SS", "Intersected SE & MXE",
                       "Intersected A5SS & A3SS", "Intersected A5SS & MXE", "Intersected A3SS & MXE")
  
  all_link_colors <- c("blue4", "red", "palevioletred1", "yellow", "magenta4", "cyan1", "darkgreen", "yellow4", "orange", "black")
  
  if (!is.null(Links_All) && nrow(Links_All) > 0) {
    present_links <- intersect(all_link_labels, unique(Links_All$Link_Type))
    if (length(present_links) > 0) {
      matched_indices <- match(present_links, all_link_labels)
      active_legends$link <- Legend(
        at = all_link_labels[matched_indices], 
        type = "lines", 
        legend_gp = gpar(col = all_link_colors[matched_indices], lwd = 2),
        title = "Correlation links", 
        labels_gp = gpar(fontsize=12), 
        title_gp = gpar(fontsize=16), 
        grid_width = unit(0.5, "cm"), 
        grid_height = unit(5,"mm"), 
        background = "white"
      )
    }
  }
  
  upViewport()
  
  if (length(active_legends) > 0) {
    lgd_list = do.call(packLegend, c(active_legends, list(direction = "horizontal", gap = unit(5, "mm"))))
    draw(lgd_list, y = unit(1, "npc") - circle_size, just = "top")
  }
  
  circos.clear()
  dev.off()
}}

if (!file.exists("Splicing_significant_results.csv")) {
  header <- "GeneID;IncLevel1;IncLevel2;FDR;IncLevelDifference;ID\n"
  writeLines(header, "Splicing_significant_results.csv")
}
