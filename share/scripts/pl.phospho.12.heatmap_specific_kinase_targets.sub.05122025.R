###############################################################
## SIMPLIFIED KINASE TARGET HEATMAPS
## Top 10 targets per kinase - LogFC display
## Now supports kinases with single targets (e.g., MARK2)
###############################################################

library(pheatmap)
library(dplyr)
library(stringr)
library(PhosR)

setwd("D:/Research/CXCR7_platelet_analysis/SubProjects/CXCR7_validation/phosphoproteomics/scripts/R")

# Load PhosphoSitePlus database
data("PhosphoSitePlus")

###############################################################
# 1. SETUP - DEFINE YOUR KINASES ----
###############################################################

kinases_of_interest <- c(
  "LYN", "SRC", "FYN",           # Src family kinases
  "PRKD2", "PRKCA", "PRKACA",    # PKC/PKA family
  "PRKG1",                       # cGMP pathway (PKG)
  "MARK2",                       # AMPK-related (cytoskeleton)
  "MTOR", "MAPK14"               # mTOR/p38 MAPK
)

###############################################################
# 2. PREPARE DATA STRUCTURES ----
###############################################################

# Your normalized intensity data
dataset_df <- as.data.frame(norm_intensity_filter)

# Extract gene and site for matching
name <- sapply(strsplit(rownames(dataset_df), ";"), `[`, 2)
site <- sapply(strsplit(rownames(dataset_df), ";"), `[`, 3)
dataset_df$namesite <- paste0(name, ";", site, ";")

# Use CXCR7 vs DMSO comparison
top.all.input <- top.all

###############################################################
# 3. FUNCTION: PLOT ONE KINASE ----
###############################################################

plot_kinase_targets <- function(kinase_name, top_n = 10) {
  
  cat("Processing:", kinase_name, "\n")
  
  ## ----------------------------------------------------------
  ## A) Get known targets from PhosphoSitePlus
  ## ----------------------------------------------------------
  kinase_targets <- PhosphoSite.human[[kinase_name]]
  
  if (is.null(kinase_targets) || length(kinase_targets) == 0) {
    message("  No targets found for ", kinase_name)
    return(NULL)
  }
  
  ## ----------------------------------------------------------
  ## B) Filter dataset to detected targets
  ## ----------------------------------------------------------
  targets_data <- dataset_df[dataset_df$namesite %in% kinase_targets, ]
  
  if (nrow(targets_data) == 0) {
    message("  No detected targets for ", kinase_name)
    return(NULL)
  }
  
  ## ----------------------------------------------------------
  ## C) Get logFC and p-values
  ## ----------------------------------------------------------
  targets_stats <- top.all.input[rownames(top.all.input) %in% rownames(targets_data), ]
  
  # Calculate max absolute logFC for ranking
  targets_stats$max_absFC <- pmax(
    abs(targets_stats$logFC.10),
    abs(targets_stats$logFC.600),
    abs(targets_stats$logFC.1800),
    na.rm = TRUE
  )
  
  # Select top N targets (or all available if less than N)
  n_available <- min(top_n, nrow(targets_stats))
  top_targets <- targets_stats %>%
    arrange(desc(max_absFC)) %>%
    head(n_available)
  
  cat("  Found", nrow(top_targets), "targets\n")
  
  ## ----------------------------------------------------------
  ## D) Build logFC matrix
  ## ----------------------------------------------------------
  logfc_matrix <- top_targets[, c("logFC.10", "logFC.600", "logFC.1800")]
  logfc_matrix$logFC.00 <- 0  # Add baseline
  logfc_matrix <- logfc_matrix[, c(4, 1, 2, 3)]
  colnames(logfc_matrix) <- c("0s", "10s", "600s", "1800s")
  
  # Clean row names
  rownames(logfc_matrix) <- sapply(strsplit(rownames(logfc_matrix), ";"), 
                                   function(x) paste0(x[2], " (", x[3], ")"))
  
  ## ----------------------------------------------------------
  ## E) Build significance matrix
  ## ----------------------------------------------------------
  sig_matrix <- matrix("", nrow = nrow(logfc_matrix), ncol = ncol(logfc_matrix),
                       dimnames = dimnames(logfc_matrix))
  
  for (i in 1:nrow(sig_matrix)) {
    if (!is.na(top_targets$adj.P.Val.10[i]) && top_targets$adj.P.Val.10[i] < 0.05) 
      sig_matrix[i, "10s"] <- "*"
    if (!is.na(top_targets$adj.P.Val.600[i]) && top_targets$adj.P.Val.600[i] < 0.05) 
      sig_matrix[i, "600s"] <- "*"
    if (!is.na(top_targets$adj.P.Val.1800[i]) && top_targets$adj.P.Val.1800[i] < 0.05) 
      sig_matrix[i, "1800s"] <- "*"
  }
  
  ## ----------------------------------------------------------
  ## F) Create heatmap (works with 1 or more targets)
  ## ----------------------------------------------------------
  colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
  breaks <- seq(-2, 2, length.out = 101)
  
  # Determine if we can cluster (need at least 2 rows)
  can_cluster <- nrow(logfc_matrix) >= 2
  
  # Save as PDF
  out_dir <- "D:/Research/CXCR7_platelet_analysis/figures/Supplementary_Figures/Kinase_Targets"
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  pdf(file.path(out_dir, paste0(kinase_name, "_top10_logFC.pdf")),
      width = 6, height = max(4, nrow(logfc_matrix) * 0.3))
  
  pheatmap(
    logfc_matrix,
    color = colors,
    breaks = breaks,
    cluster_rows = can_cluster,  # Only cluster if >= 2 rows
    cluster_cols = FALSE,
    scale = "none",
    display_numbers = sig_matrix,
    number_color = "black",
    fontsize = 12,
    fontsize_row = 9,
    fontsize_col = 11,
    fontsize_number = 12,
    cellwidth = 30,
    cellheight = 14,
    border_color = "grey70",
    main = paste0(kinase_name, " - Top ", nrow(logfc_matrix), " Target", 
                  ifelse(nrow(logfc_matrix) > 1, "s", ""), " (logFC)"),
    na_col = "grey90"
  )
  
  dev.off()
  
  # Also save as TIFF
  tiff(file.path(out_dir, paste0(kinase_name, "_top10_logFC.tiff")),
       width = 6, height = max(4, nrow(logfc_matrix) * 0.3), 
       units = "in", res = 300, compression = "lzw")
  
  pheatmap(
    logfc_matrix,
    color = colors,
    breaks = breaks,
    cluster_rows = can_cluster,  # Only cluster if >= 2 rows
    cluster_cols = FALSE,
    scale = "none",
    display_numbers = sig_matrix,
    number_color = "black",
    fontsize = 12,
    fontsize_row = 9,
    fontsize_col = 11,
    fontsize_number = 12,
    cellwidth = 30,
    cellheight = 14,
    border_color = "grey70",
    main = paste0(kinase_name, " - Top ", nrow(logfc_matrix), " Target", 
                  ifelse(nrow(logfc_matrix) > 1, "s", ""), " (logFC)"),
    na_col = "grey90"
  )
  
  dev.off()
  
  cat("  Saved:", file.path(out_dir, paste0(kinase_name, "_top10_logFC.pdf")), "\n")
  
  return(logfc_matrix)
}

###############################################################
# 4. RUN FOR ALL KINASES ----
###############################################################

cat("\n========================================\n")
cat("GENERATING KINASE TARGET HEATMAPS\n")
cat("========================================\n\n")

results <- list()

for (kinase in kinases_of_interest) {
  results[[kinase]] <- plot_kinase_targets(kinase, top_n = 10)
  cat("\n")
}

###############################################################
# 5. CREATE COMBINED MULTI-PANEL FIGURE (PNG) ----
###############################################################

cat("Creating combined multi-panel figure...\n")

# Filter out NULL results
valid_results <- results[!sapply(results, is.null)]

if (length(valid_results) > 0) {
  
  library(gridExtra)
  library(grid)
  
  # Calculate layout
  n_kinases <- length(valid_results)
  n_cols <- 3
  n_rows <- ceiling(n_kinases / n_cols)
  
  # Create list of plots
  plot_list <- lapply(names(valid_results), function(kinase_name) {
    
    logfc_mat <- valid_results[[kinase_name]]
    
    # Get original row names to retrieve p-values
    gene_names <- sapply(strsplit(rownames(logfc_mat), " \\("), `[`, 1)
    site_names <- gsub(".*\\(|\\)", "", rownames(logfc_mat))
    
    # Build search patterns
    search_patterns <- paste0(gene_names, ";", site_names)
    
    # Find matching rows in top.all.input
    matching_rows <- sapply(search_patterns, function(pattern) {
      matches <- grep(pattern, rownames(top.all.input), fixed = TRUE)
      if (length(matches) > 0) matches[1] else NA
    })
    
    # Create significance matrix
    sig_mat <- matrix("", nrow = nrow(logfc_mat), ncol = ncol(logfc_mat),
                      dimnames = dimnames(logfc_mat))
    
    for (i in 1:nrow(sig_mat)) {
      if (!is.na(matching_rows[i])) {
        row_data <- top.all.input[matching_rows[i], ]
        if (!is.na(row_data$adj.P.Val.10) && row_data$adj.P.Val.10 < 0.05) 
          sig_mat[i, "10s"] <- "*"
        if (!is.na(row_data$adj.P.Val.600) && row_data$adj.P.Val.600 < 0.05) 
          sig_mat[i, "600s"] <- "*"
        if (!is.na(row_data$adj.P.Val.1800) && row_data$adj.P.Val.1800 < 0.05) 
          sig_mat[i, "1800s"] <- "*"
      }
    }
    
    colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
    breaks <- seq(-2, 2, length.out = 101)
    
    pheatmap(
      logfc_mat,
      color = colors,
      breaks = breaks,
      cluster_rows = (nrow(logfc_mat) >= 2),  # Only cluster if >= 2 rows
      cluster_cols = FALSE,
      scale = "none",
      display_numbers = sig_mat,
      number_color = "black",
      fontsize = 8,
      fontsize_row = 6,
      fontsize_col = 7,
      fontsize_number = 7,
      cellwidth = 15,
      cellheight = 6,
      border_color = "grey80",
      main = paste0(kinase_name),
      legend = FALSE,
      silent = TRUE
    )
  })
  
  # Save as high-res PNG (compressed layout with minimal spacing)
  png("D:/Research/CXCR7_platelet_analysis/figures/Supplementary_Figures/Kinase_Targets/All_Kinases_Top10_Combined.png",
      width = 12, height = 3.5 * n_rows, units = "in", res = 600)
  
  grid.arrange(grobs = lapply(plot_list, function(x) x[[4]]), 
               ncol = n_cols, nrow = n_rows,
               padding = unit(0.5, "line"))
  
  dev.off()
  
  # Also save TIFF version for publication
  tiff("D:/Research/CXCR7_platelet_analysis/figures/Supplementary_Figures/Kinase_Targets/All_Kinases_Top10_Combined.tiff",
       width = 12, height = 3.5 * n_rows, units = "in", res = 600, compression = "lzw")
  
  grid.arrange(grobs = lapply(plot_list, function(x) x[[4]]), 
               ncol = n_cols, nrow = n_rows,
               padding = unit(0.5, "line"))
  
  dev.off()
  
  # Also save as PDF (individual pages)
  pdf("D:/Research/CXCR7_platelet_analysis/figures/Supplementary_Figures/Kinase_Targets/All_Kinases_Top10_Combined.pdf",
      width = 6, height = 8)
  
  for (kinase_name in names(valid_results)) {
    logfc_mat <- valid_results[[kinase_name]]
    
    # Get significance
    gene_names <- sapply(strsplit(rownames(logfc_mat), " \\("), `[`, 1)
    site_names <- gsub(".*\\(|\\)", "", rownames(logfc_mat))
    search_patterns <- paste0(gene_names, ";", site_names)
    
    matching_rows <- sapply(search_patterns, function(pattern) {
      matches <- grep(pattern, rownames(top.all.input), fixed = TRUE)
      if (length(matches) > 0) matches[1] else NA
    })
    
    sig_mat <- matrix("", nrow = nrow(logfc_mat), ncol = ncol(logfc_mat),
                      dimnames = dimnames(logfc_mat))
    
    for (i in 1:nrow(sig_mat)) {
      if (!is.na(matching_rows[i])) {
        row_data <- top.all.input[matching_rows[i], ]
        if (!is.na(row_data$adj.P.Val.10) && row_data$adj.P.Val.10 < 0.05) 
          sig_mat[i, "10s"] <- "*"
        if (!is.na(row_data$adj.P.Val.600) && row_data$adj.P.Val.600 < 0.05) 
          sig_mat[i, "600s"] <- "*"
        if (!is.na(row_data$adj.P.Val.1800) && row_data$adj.P.Val.1800 < 0.05) 
          sig_mat[i, "1800s"] <- "*"
      }
    }
    
    colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
    breaks <- seq(-2, 2, length.out = 101)
    
    print(pheatmap(
      logfc_mat,
      color = colors,
      breaks = breaks,
      cluster_rows = (nrow(logfc_mat) >= 2),  # Only cluster if >= 2 rows
      cluster_cols = FALSE,
      scale = "none",
      display_numbers = sig_mat,
      number_color = "black",
      fontsize = 12,
      fontsize_row = 9,
      fontsize_col = 11,
      fontsize_number = 11,
      cellwidth = 30,
      cellheight = 14,
      border_color = "grey70",
      main = paste0(kinase_name, " - Top Target", 
                    ifelse(nrow(logfc_mat) > 1, "s", "")),
      silent = FALSE
    ))
  }
  
  dev.off()
  
  cat("\n✓ Combined figures saved!\n")
  cat("  PNG (grid):   All_Kinases_Top10_Combined.png\n")
  cat("  TIFF (grid):  All_Kinases_Top10_Combined.tiff\n")
  cat("  PDF (pages):  All_Kinases_Top10_Combined.pdf\n")
}

###############################################################
# 6. CREATE SUMMARY TABLE ----
###############################################################

cat("\nCreating summary table...\n")

if (length(valid_results) > 0) {
  
  summary_data <- lapply(names(valid_results), function(k) {
    
    logfc_mat <- valid_results[[k]]
    
    # Extract gene and site names
    gene_names <- sapply(strsplit(rownames(logfc_mat), " \\("), `[`, 1)
    site_names <- gsub(".*\\(|\\)", "", rownames(logfc_mat))
    
    data.frame(
      Kinase = k,
      Gene = gene_names,
      Site = site_names,
      Target = rownames(logfc_mat),
      logFC_10s = round(logfc_mat[, "10s"], 2),
      logFC_600s = round(logfc_mat[, "600s"], 2),
      logFC_1800s = round(logfc_mat[, "1800s"], 2),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
  
  write.csv(summary_data,
            "D:/Research/CXCR7_platelet_analysis/figures/Supplementary_Figures/Kinase_Targets/Kinase_Targets_Summary.csv",
            row.names = FALSE)
  
  cat("  Summary table saved: Kinase_Targets_Summary.csv\n")
}

cat("\n========================================\n")
cat("✓ ANALYSIS COMPLETE!\n")
cat("========================================\n")
cat("\nOutput files:\n")
cat("  Individual PDFs/TIFFs: Kinase_Targets/\n")
cat("  Combined PNG (grid):   All_Kinases_Top10_Combined.png\n")
cat("  Combined TIFF (grid):  All_Kinases_Top10_Combined.tiff\n")
cat("  Combined PDF (pages):  All_Kinases_Top10_Combined.pdf\n")
cat("  Summary table:         Kinase_Targets_Summary.csv\n\n")