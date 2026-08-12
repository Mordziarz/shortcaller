#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
gtf_file         <- args[1]
longest_orfs_gff <- args[2]
sig_splicing     <- args[3]
cpc2_file        <- args[4]
all_splicing     <- args[5]  
out_dir          <- args[6]

suppressPackageStartupMessages({
  library(rtracklayer)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
})

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

gtf <- import(gtf_file)
gtf_df <- as.data.frame(gtf)

orf_gff_df <- data.frame(
  seqnames = character(), source = character(), type = character(), 
  start = integer(), end = integer(), score = character(), 
  strand = character(), phase = character(), attributes = character(),
  stringsAsFactors = FALSE
)

if(length(longest_orfs_gff) > 0 && nzchar(longest_orfs_gff) && file.exists(longest_orfs_gff) && file.info(longest_orfs_gff)$size > 0) {
  gff_lines <- readLines(longest_orfs_gff, warn = FALSE)
  gff_lines <- gff_lines[!grepl("^#", gff_lines) & nchar(trimws(gff_lines)) > 0]
  
  if(length(gff_lines) > 0) {
    temp_df <- tryCatch(
      read.delim(text = gff_lines, header = FALSE, stringsAsFactors = FALSE, quote = "", fill = TRUE),
      error = function(e) data.frame()
    )
    if(nrow(temp_df) > 0) {
      if(ncol(temp_df) >= 9) {
        colnames(temp_df)[1:9] <- c("seqnames", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
      } else {
        for(i in (ncol(temp_df)+1):9) temp_df[[i]] <- ""
        colnames(temp_df)[1:9] <- c("seqnames", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
      }
      orf_gff_df <- temp_df
    }
  }
}

load_splicing_file <- function(filepath) {
  df <- data.frame()
  if(!file.exists(filepath) || file.info(filepath)$size == 0) return(df)
  
  for(sep_char in c(",", ";", "\t")) {
    temp <- tryCatch(read.csv(filepath, sep = sep_char, stringsAsFactors = FALSE), error = function(e) data.frame())
    if(nrow(temp) > 0 && ncol(temp) > 1) {
      df <- temp
      break
    }
  }
  if(nrow(df) == 0) {
    df <- tryCatch(read.csv2(filepath, stringsAsFactors = FALSE), error = function(e) data.frame())
  }
  
  if(nrow(df) > 0) {
    gene_col_candidates <- c("GeneID", "gene_id", "GeneID.", "geneSymbol", "gene_symbol", "ID")
    found_gene_col <- intersect(gene_col_candidates, colnames(df))[1]
    if(!is.na(found_gene_col)) {
      df$GeneID <- df[[found_gene_col]]
    } else {
      df$GeneID <- as.character(df[,1])
    }
    
    type_col_candidates <- c("OriginalType", "Type", "type", "event_type", "Event", "Splicing_Type")
    found_type_col <- intersect(type_col_candidates, colnames(df))[1]
    if(!is.na(found_type_col)) {
      df$Type <- df[[found_type_col]]
    } else {
      df$Type <- "Alternative_Splicing"
    }
  }
  return(df)
}

sig_df <- load_splicing_file(sig_splicing)
interest_genes <- if("GeneID" %in% colnames(sig_df)) unique(sig_df$GeneID) else c()

gtf_filtered <- if(length(interest_genes) > 0) gtf_df[gtf_df$gene_id %in% interest_genes, ] else gtf_df
if(nrow(gtf_filtered) == 0) {
  gtf_filtered <- gtf_df 
}

cds_data <- orf_gff_df[orf_gff_df$type == "CDS", ]
if(nrow(cds_data) > 0) {
  cds_data$Parent <- sub(".*Parent=([^;]+).*", "\\1", cds_data$attributes)
  cds_data$transcript_id <- gsub("\\.p[0-9]+$", "", cds_data$Parent)
  cds_data$width <- cds_data$end - cds_data$start + 1
  orf_lengths <- aggregate(width ~ transcript_id, data = cds_data, FUN = max)
  colnames(orf_lengths)[2] <- "orf_length"
} else {
  orf_lengths <- data.frame(transcript_id = character(), orf_length = integer())
}

exons <- gtf_filtered[gtf_filtered$type == "exon", ]
all_tx_ids <- unique(gtf_filtered$transcript_id)

get_pos_safe <- function(tx_id) {
  tx_exons <- exons[exons$transcript_id == tx_id, ]
  if (nrow(tx_exons) < 2) return(0)
  
  if (all(tx_exons$strand == "-")) {
    tx_exons <- tx_exons[order(tx_exons$start, decreasing = TRUE), ]
  } else {
    tx_exons <- tx_exons[order(tx_exons$start), ]
  }
  return(sum(tx_exons$width[1:(nrow(tx_exons) - 1)]))
}

all_junction_pos <- sapply(all_tx_ids, get_pos_safe)

df_junctions <- data.frame(
  transcript_id = all_tx_ids,
  last_junction_pos = all_junction_pos,
  stringsAsFactors = FALSE
)

nmd_analysis <- merge(df_junctions, orf_lengths, by = "transcript_id", all.x = TRUE)
nmd_analysis$orf_length[is.na(nmd_analysis$orf_length)] <- 0

nmd_analysis$is_NMD <- ifelse(
  nmd_analysis$last_junction_pos == 0,
  FALSE,
  nmd_analysis$orf_length < (nmd_analysis$last_junction_pos - 50)
)

if(file.exists(cpc2_file)) {
  cpc2 <- read.csv(cpc2_file, header=TRUE, sep="\t")
  id_col <- intersect(c("X.ID", "transcript_id", "ID"), colnames(cpc2))[1]
  cpc2_coding <- cpc2[cpc2$label == "coding", id_col]
  nmd_analysis$is_coding <- nmd_analysis$transcript_id %in% cpc2_coding
} else {
  nmd_analysis$is_coding <- TRUE
}

tx_to_gene <- unique(gtf_filtered[, c("transcript_id", "gene_id")])
nmd_analysis_full <- merge(nmd_analysis, tx_to_gene, by = "transcript_id")

gene_isoform_stats <- nmd_analysis_full %>%
  group_by(gene_id) %>%
  summarise(
    Total_Isoforms = n(),
    Protein_Coding = sum(is_coding),
    Non_Coding = sum(!is_coding),
    NMD_Isoforms = sum(is_NMD),
    .groups = 'drop'
  )

if(nrow(sig_df) > 0 && "GeneID" %in% colnames(sig_df) && "Type" %in% colnames(sig_df)) {
  sig_events_per_gene <- sig_df %>%
    group_by(GeneID, Type) %>%
    summarise(Count = n(), .groups = "drop") %>%
    pivot_wider(names_from = Type, values_from = Count, values_fill = 0, names_prefix = "Sig_")
} else {
  sig_events_per_gene <- data.frame(GeneID = character())
}

all_df <- load_splicing_file(all_splicing)

if(nrow(all_df) > 0 && "GeneID" %in% colnames(all_df) && "Type" %in% colnames(all_df)) {
  tot_events_per_gene <- all_df %>%
    group_by(GeneID, Type) %>%
    summarise(Total_Count = n(), .groups = "drop") %>%
    pivot_wider(
      names_from = Type, 
      values_from = Total_Count, 
      names_prefix = "All_",
      values_fill = 0
    )
  
  tot_events_per_gene <- tot_events_per_gene %>%
    mutate(total_splicing_events = rowSums(select(., where(is.numeric)), na.rm = TRUE))
} else {
  tot_events_per_gene <- data.frame(GeneID = character())
}

master_matrix_final <- gene_isoform_stats

if(nrow(sig_events_per_gene) > 0 && "GeneID" %in% colnames(sig_events_per_gene)) {
  master_matrix_final <- master_matrix_final %>%
    left_join(sig_events_per_gene, by = c("gene_id" = "GeneID"))
}
if(nrow(tot_events_per_gene) > 0 && "GeneID" %in% colnames(tot_events_per_gene)) {
  master_matrix_final <- master_matrix_final %>%
    left_join(tot_events_per_gene, by = c("gene_id" = "GeneID"))
}

master_matrix_final[is.na(master_matrix_final)] <- 0

sig_cols_check <- grep("^Sig_", colnames(master_matrix_final), value = TRUE)
for(sig_col in sig_cols_check) {
  base_name <- gsub("^Sig_", "", sig_col)
  all_col <- paste0("All_", base_name)
  nosig_col <- paste0("NoSig_", base_name)
  
  if(all_col %in% colnames(master_matrix_final)) {
    master_matrix_final[[nosig_col]] <- pmax(0, master_matrix_final[[all_col]] - master_matrix_final[[sig_col]])
  } else {
    master_matrix_final[[nosig_col]] <- 0
  }
}

output_csv <- file.path(out_dir, "consequence_matrix.csv")
write.csv2(master_matrix_final, output_csv, row.names = FALSE)
cat("Saved complete consequence matrix to:", output_csv, "\n")

all_genes <- unique(master_matrix_final$gene_id)

if(length(all_genes) > 0) {
  for(g in all_genes) {
    df_geny <- master_matrix_final %>%
      filter(gene_id == g)
    
    if(nrow(df_geny) == 0) next
    
    rownames(df_geny) <- df_geny$gene_id
    
    sig_cols <- grep("^Sig_", colnames(df_geny), value = TRUE)
    nosig_cols <- grep("^NoSig_", colnames(df_geny), value = TRUE)
    
    if(length(sig_cols) > 0) {
      mat_sig <- t(as.matrix(df_geny[, sig_cols, drop = FALSE]))
      rownames(mat_sig) <- gsub("Sig_", "", rownames(mat_sig))
      rownames(mat_sig) <- gsub("_", " ", rownames(mat_sig))
      
      max_sig_val <- max(mat_sig, na.rm = TRUE)
      if(max_sig_val == 0) max_sig_val <- 1
      col_sig <- colorRamp2(c(0, max_sig_val), c("#f7fbff", "#084594"))
      
      if(length(nosig_cols) > 0) {
        mat_nosig <- t(as.matrix(df_geny[, nosig_cols, drop = FALSE]))
        rownames(mat_nosig) <- gsub("NoSig_", "", rownames(mat_nosig))
        rownames(mat_nosig) <- gsub("_", " ", rownames(mat_nosig))
      } else {
        mat_nosig <- matrix(0, nrow = nrow(mat_sig), ncol = ncol(mat_sig), dimnames = dimnames(mat_sig))
      }
      
      max_nosig_val <- max(mat_nosig, na.rm = TRUE)
      if(max_nosig_val == 0) max_nosig_val <- 1
      col_nosig <- colorRamp2(c(0, max_nosig_val), c("#fffff5", "#990000"))
      
      df_barplot <- as.matrix(df_geny[, intersect(c("Total_Isoforms", "Protein_Coding", "NMD_Isoforms"), colnames(df_geny)), drop = FALSE])
      
      barplot_colors <- c("Total_Isoforms" = "#41ab5d", "Protein_Coding" = "#ef3b2c", "NMD_Isoforms" = "#fec44f")
      barplot_colors <- barplot_colors[intersect(names(barplot_colors), colnames(df_barplot))]
      if(length(barplot_colors) == 0) barplot_colors <- c("Total_Isoforms" = "#41ab5d")
      
      dol_barplot <- HeatmapAnnotation(
        "Isoform Profile" = anno_barplot(
          df_barplot, 
          beside = TRUE, 
          gp = gpar(fill = barplot_colors), 
          bar_width = 0.8,
          height = unit(4, "cm")
        ),
        annotation_name_side = "left",
        annotation_name_rot = 90,
        which = "column"
      )
      
      breaks_sig <- unique(round(seq(0, max_sig_val, length.out = min(3, max_sig_val + 1))))
      breaks_nosig <- unique(round(seq(0, max_nosig_val, length.out = min(3, max_nosig_val + 1))))
      
      lgd_barplot = Legend(labels = names(barplot_colors), title = "Isoform Composition", legend_gp = gpar(fill = barplot_colors), ncol = 1)
      lgd_sig = Legend(title = "Significant AS", col_fun = col_sig, at = breaks_sig, labels = breaks_sig)
      lgd_nosig = Legend(title = "Not significant AS", col_fun = col_nosig, at = breaks_nosig, labels = breaks_nosig)
      
      all_legends = packLegend(lgd_sig, lgd_nosig, lgd_barplot, direction = "vertical", gap = unit(4, "mm"))
      
      ht_list <- 
        Heatmap(mat_sig, 
                name = "Significant AS", 
                col = col_sig,
                rect_gp = gpar(col = "black", lwd = 1), 
                cluster_rows = FALSE, cluster_columns = FALSE,
                column_title = "",
                row_names_gp = gpar(fontsize = 10),
                show_column_names = FALSE, 
                heatmap_legend_param = list(show=FALSE), show_heatmap_legend = FALSE) %v%
        
        Heatmap(mat_nosig, 
                name = "Not significant AS", 
                col = col_nosig,
                rect_gp = gpar(col = "black", lwd = 1), 
                cluster_rows = FALSE, cluster_columns = FALSE,
                column_title = "",
                row_names_gp = gpar(fontsize = 10),
                column_names_gp = gpar(fontface = "bold.italic", fontsize = 13), 
                bottom_annotation = dol_barplot,
                heatmap_legend_param = list(show=FALSE), show_heatmap_legend = FALSE)
      
      png_path <- file.path(out_dir, paste0("splicing_consequence_", g, ".png"))
      
      png(png_path, width = 6, height = 6, res = 600, units = "in")
      
      draw(ht_list, 
           column_title = paste0("Gene: ", g),
           column_title_gp = gpar(fontsize = 14, fontface = "bold"),
           heatmap_legend_side = "right", 
           annotation_legend_side = "right",
           annotation_legend_list = list(all_legends),
           padding = unit(c(2, 2, 2, 15), "mm"))
      dev.off()
    }
  }
  cat("Generated separate plots for all", length(all_genes), "genes.\n")
}