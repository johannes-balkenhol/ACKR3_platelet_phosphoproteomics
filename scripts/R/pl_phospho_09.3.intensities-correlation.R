###############################################################
## Publication-Ready Correlation Heatmap
## CXCR7 Normalized Intensities — Initial vs Validation
## 
## Run AFTER pl_phospho_09_2_intensities-comparative-analysis.R
## Requires: val_norm_clean, init_norm_clean (harmonized)
###############################################################

library(pheatmap)
library(dplyr)

###############################################################
# 1) SELECT CXCR7 SAMPLES AT 10s, 600s, 1800s
###############################################################

# After harmonization, column names look like:
#   VAL:  "CXCR7.10s_01", "CXCR7.600s_03", "DMSO.10s_05", ...
#   INIT: "CXCR7.10s_01", "CXCR7.600s_07", ...

# Select only CXCR7 at shared timepoints
select_cxcr7 <- function(mat, times = c("10s", "600s", "1800s")) {
  cols <- grep("^CXCR7\\.", colnames(mat), value = TRUE)
  # Filter to desired timepoints
  keep <- cols[sapply(cols, function(c) {
    tp <- sub("^CXCR7\\.(\\d+s)_.*", "\\1", c)
    tp %in% times
  })]
  mat[, keep, drop = FALSE]
}

val_cxcr7  <- select_cxcr7(val_norm_clean)
init_cxcr7 <- select_cxcr7(init_norm_clean)

cat("VAL CXCR7 samples: ", ncol(val_cxcr7),  "→", paste(colnames(val_cxcr7), collapse=", "), "\n")
cat("INIT CXCR7 samples:", ncol(init_cxcr7), "→", paste(colnames(init_cxcr7), collapse=", "), "\n")

###############################################################
# 2) INTERSECT PHOSPHOSITES & COMBINE
###############################################################

shared_ids <- intersect(rownames(val_cxcr7), rownames(init_cxcr7))
cat("Shared phosphosites:", length(shared_ids), "\n")

val_sub  <- val_cxcr7[shared_ids, ]
init_sub <- init_cxcr7[shared_ids, ]

# Prefix columns to distinguish datasets
colnames(val_sub)  <- paste0("VAL_",  colnames(val_sub))
colnames(init_sub) <- paste0("INIT_", colnames(init_sub))

combined <- cbind(val_sub, init_sub)

###############################################################
# 3) COMPUTE SAMPLE-LEVEL CORRELATION
###############################################################

cor_mat <- cor(combined, method = "pearson", use = "pairwise.complete.obs")

###############################################################
# 4) PUBLICATION-READY LABELS
###############################################################

# Make clean labels: "VAL_CXCR7.600s_03" → "Val CXCR7 600s #03"
make_pretty <- function(x) {
  dataset <- ifelse(grepl("^VAL_", x), "Val", "Init")
  time    <- sub(".*CXCR7\\.(\\d+s)_.*", "\\1", x)
  rep     <- sub(".*_(\\d+)$", "\\1", x)
  paste0(dataset, " ", time, " #", rep)
}

pretty_labs <- sapply(colnames(combined), make_pretty)
rownames(cor_mat) <- pretty_labs
colnames(cor_mat) <- pretty_labs

###############################################################
# 5) ANNOTATION
###############################################################

# Extract metadata for annotations
annotation_df <- data.frame(
  Dataset   = factor(ifelse(grepl("^Val", pretty_labs), "Validation", "Initial"),
                     levels = c("Validation", "Initial")),
  Timepoint = factor(sub("^(Val|Init) (\\d+s) .*", "\\2", pretty_labs),
                     levels = c("10s", "600s", "1800s")),
  row.names = pretty_labs
)

# Color scheme
ann_colors <- list(
  Dataset   = c("Validation" = "#2C7BB6",    # blue
                "Initial"    = "#D7191C"),    # red
  Timepoint = c("10s"   = "#ABDDA4",         # light green
                "600s"  = "#FDAE61",          # orange
                "1800s" = "#D53E4F")          # dark red
)

###############################################################
# 6) COLOR PALETTE (diverging, centered)
###############################################################

n_colors <- 200

# Determine symmetric breaks around median
cor_range <- range(cor_mat, na.rm = TRUE)
cat("Correlation range:", round(cor_range, 3), "\n")

color_palette <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE",
                                    "#D1E5F0", "#F7F7F7",
                                    "#FDDBC7", "#F4A582",
                                    "#D6604D", "#B2182B"))(n_colors)

# For intensities, correlations are typically all positive (0.7–1.0)
# Use a sequential-ish palette that zooms into that range
breaks <- seq(floor(cor_range[1] * 20) / 20,   # round down to nearest 0.05
              1,
              length.out = n_colors + 1)

###############################################################
# 7) PLOT — TIFF (300 dpi)
###############################################################

tiff(filename = "analysis/PCA/correlation_heatmap_cxcr7_intensities.tiff",
     width = 10 * 300, height = 9 * 300,
     res = 300, compression = "lzw")

pheatmap(
  cor_mat,
  
  # Clustering
  clustering_method        = "complete",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  
  # Colors
  color  = color_palette,
  breaks = breaks,
  
  # Annotations
  annotation_col    = annotation_df,
  annotation_row    = annotation_df,
  annotation_colors = ann_colors,
  
  # Display correlation values
  display_numbers  = TRUE,
  number_format    = "%.2f",
  number_color     = "black",
  fontsize_number  = 6,
  
  # Fonts
  fontsize     = 10,
  fontsize_row = 8,
  fontsize_col = 8,
  
  # Appearance
  border_color   = "grey90",
  treeheight_row = 35,
  treeheight_col = 35,
  
  # Title
  main = "CXCR7 Normalized Intensity Correlation\nValidation vs Initial Experiment"
)

dev.off()

###############################################################
# 8) PLOT — PDF (vector, for journal submission)
###############################################################

png(filename = "analysis/PCA/correlation_heatmap_cxcr7_intensities.png",
    width = 10 * 300, height = 9 * 300,
    res = 300)

pheatmap(
  cor_mat,
  clustering_method        = "complete",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  color  = color_palette,
  breaks = breaks,
  annotation_col    = annotation_df,
  annotation_row    = annotation_df,
  annotation_colors = ann_colors,
  display_numbers  = TRUE,
  number_format    = "%.2f",
  number_color     = "black",
  fontsize_number  = 6,
  fontsize     = 10,
  fontsize_row = 8,
  fontsize_col = 8,
  border_color   = "grey90",
  treeheight_row = 35,
  treeheight_col = 35,
  main = "CXCR7 Normalized Intensity Correlation\nValidation vs Initial Experiment"
)

dev.off()

cat("\n✓ Saved: correlation_heatmap_cxcr7_intensities.tiff\n")
cat("✓ Saved: correlation_heatmap_cxcr7_intensities.pdf\n")












###############################################################
# 10) SAME ANALYSIS ON RAW INTENSITIES
###############################################################

val_cxcr7_raw  <- select_cxcr7(val_raw_clean)
init_cxcr7_raw <- select_cxcr7(init_raw_clean)

shared_ids_raw <- intersect(rownames(val_cxcr7_raw), rownames(init_cxcr7_raw))
cat("Shared phosphosites (raw):", length(shared_ids_raw), "\n")

val_sub_raw  <- val_cxcr7_raw[shared_ids_raw, ]
init_sub_raw <- init_cxcr7_raw[shared_ids_raw, ]

# Z-score per experiment
val_z_raw  <- t(scale(t(val_sub_raw)))
init_z_raw <- t(scale(t(init_sub_raw)))

good_rows_raw <- complete.cases(val_z_raw) & complete.cases(init_z_raw)
val_z_raw  <- val_z_raw[good_rows_raw, ]
init_z_raw <- init_z_raw[good_rows_raw, ]
cat("Phosphosites after z-score filter (raw):", sum(good_rows_raw), "\n")

colnames(val_z_raw)  <- paste0("VAL_",  colnames(val_z_raw))
colnames(init_z_raw) <- paste0("INIT_", colnames(init_z_raw))

combined_z_raw <- cbind(val_z_raw, init_z_raw)
cor_mat_raw <- cor(combined_z_raw, method = "pearson", use = "pairwise.complete.obs")

cat("Correlation range (raw z-scored):", round(range(cor_mat_raw), 3), "\n")

# Labels & annotations (reuse make_pretty and annotation setup)
pretty_labs_raw <- sapply(colnames(combined_z_raw), make_pretty)
rownames(cor_mat_raw) <- pretty_labs_raw
colnames(cor_mat_raw) <- pretty_labs_raw

annotation_df_raw <- data.frame(
  Dataset   = factor(ifelse(grepl("^Val", pretty_labs_raw), "Validation", "Initial"),
                     levels = c("Validation", "Initial")),
  Timepoint = factor(sub("^(Val|Init) (\\d+s) .*", "\\2", pretty_labs_raw),
                     levels = c("10s", "600s", "1800s")),
  row.names = pretty_labs_raw
)

# Use same color scale as normalized for direct comparison
cor_range_raw <- range(cor_mat_raw[cor_mat_raw < 1], na.rm = TRUE)
max_abs_raw <- max(abs(cor_range_raw))
breaks_raw <- seq(-max_abs_raw, 1, length.out = n_colors + 1)

png(filename = "analysis/PCA/correlation_heatmap_cxcr7_zscore_combined_RAW.png",
    width = 10 * 300, height = 9 * 300, res = 300)

pheatmap(
  cor_mat_raw,
  clustering_method        = "complete",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  color  = color_palette,
  breaks = breaks_raw,
  annotation_col    = annotation_df_raw,
  annotation_row    = annotation_df_raw,
  annotation_colors = ann_colors,
  display_numbers   = FALSE,
  fontsize          = 10,
  fontsize_row      = 7,
  fontsize_col      = 7,
  border_color      = NA,
  treeheight_row    = 35,
  treeheight_col    = 35,
  main = "CXCR7 RAW Intensity Correlation (Z-scored per Experiment)\nValidation vs Initial"
)

dev.off()

cat("✓ correlation_heatmap_cxcr7_zscore_combined_RAW.png\n")








###############################################################
# CORRELATION OF AVERAGED TIMEPOINT PROFILES
# (median per timepoint, then correlate)
###############################################################

###############################################################
# 1) Helper: average replicates per condition.timepoint
###############################################################

average_by_timepoint <- function(mat) {
  # colnames like "CXCR7.10s_01" → group = "CXCR7.10s"
  groups <- sub("_\\d+$", "", colnames(mat))
  unique_groups <- unique(groups)
  
  avg_mat <- sapply(unique_groups, function(g) {
    cols <- which(groups == g)
    rowMedians <- apply(mat[, cols, drop = FALSE], 1, median, na.rm = TRUE)
    rowMedians
  })
  
  colnames(avg_mat) <- unique_groups
  avg_mat
}

###############################################################
# 2) Shared phosphosites (same as before)
###############################################################

shared_ids <- intersect(rownames(val_cxcr7), rownames(init_cxcr7))

# Also get raw CXCR7 subsets
val_cxcr7_raw  <- select_cxcr7(val_raw_clean)
init_cxcr7_raw <- select_cxcr7(init_raw_clean)
shared_ids_raw <- intersect(rownames(val_cxcr7_raw), rownames(init_cxcr7_raw))

###############################################################
# 3) NORMALIZED: average, z-score, combine, correlate
###############################################################

val_avg_norm  <- average_by_timepoint(val_cxcr7[shared_ids, ])
init_avg_norm <- average_by_timepoint(init_cxcr7[shared_ids, ])

# Z-score per experiment
val_avg_norm_z  <- t(scale(t(val_avg_norm)))
init_avg_norm_z <- t(scale(t(init_avg_norm)))

good_norm <- complete.cases(val_avg_norm_z) & complete.cases(init_avg_norm_z)

colnames(val_avg_norm_z)  <- paste0("Val_",  colnames(val_avg_norm_z))
colnames(init_avg_norm_z) <- paste0("Init_", colnames(init_avg_norm_z))

combined_avg_norm <- cbind(val_avg_norm_z[good_norm, ], init_avg_norm_z[good_norm, ])
cor_avg_norm <- cor(combined_avg_norm, method = "pearson", use = "pairwise.complete.obs")

###############################################################
# 4) RAW: average, z-score, combine, correlate
###############################################################

val_avg_raw  <- average_by_timepoint(val_cxcr7_raw[shared_ids_raw, ])
init_avg_raw <- average_by_timepoint(init_cxcr7_raw[shared_ids_raw, ])

val_avg_raw_z  <- t(scale(t(val_avg_raw)))
init_avg_raw_z <- t(scale(t(init_avg_raw)))

good_raw <- complete.cases(val_avg_raw_z) & complete.cases(init_avg_raw_z)

colnames(val_avg_raw_z)  <- paste0("Val_",  colnames(val_avg_raw_z))
colnames(init_avg_raw_z) <- paste0("Init_", colnames(init_avg_raw_z))

combined_avg_raw <- cbind(val_avg_raw_z[good_raw, ], init_avg_raw_z[good_raw, ])
cor_avg_raw <- cor(combined_avg_raw, method = "pearson", use = "pairwise.complete.obs")

###############################################################
# 5) Pretty labels
###############################################################

make_pretty_avg <- function(x) {
  dataset <- ifelse(grepl("^Val_", x), "Val", "Init")
  time    <- sub(".*CXCR7\\.(\\d+s)$", "\\1", sub("^(Val|Init)_", "", x))
  paste0(dataset, " CXCR7 ", time)
}

rownames(cor_avg_norm) <- sapply(colnames(combined_avg_norm), make_pretty_avg)
colnames(cor_avg_norm) <- rownames(cor_avg_norm)

rownames(cor_avg_raw) <- sapply(colnames(combined_avg_raw), make_pretty_avg)
colnames(cor_avg_raw) <- rownames(cor_avg_raw)

###############################################################
# 6) Annotations
###############################################################

make_annotation <- function(labs) {
  data.frame(
    Dataset   = factor(ifelse(grepl("^Val", labs), "Validation", "Initial"),
                       levels = c("Validation", "Initial")),
    Timepoint = factor(sub("^(Val|Init) CXCR7 ", "", labs),
                       levels = c("10s", "600s", "1800s")),
    row.names = labs
  )
}

ann_norm <- make_annotation(rownames(cor_avg_norm))
ann_raw  <- make_annotation(rownames(cor_avg_raw))

ann_colors <- list(
  Dataset   = c("Validation" = "#2C7BB6", "Initial" = "#D7191C"),
  Timepoint = c("10s" = "#ABDDA4", "600s" = "#FDAE61", "1800s" = "#D53E4F")
)

###############################################################
# 7) PLOT — NORMALIZED (averaged)
###############################################################

n_colors <- 200
palette_div <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE",
                                  "#D1E5F0", "#F7F7F7",
                                  "#FDDBC7", "#F4A582",
                                  "#D6604D", "#B2182B"))(n_colors)

png(filename = "analysis/PCA/correlation_heatmap_cxcr7_NORM_averaged.png",
    width = 6 * 300, height = 5.5 * 300, res = 300)

pheatmap(
  cor_avg_norm,
  clustering_method = "complete",
  color  = palette_div,
  breaks = seq(-1, 1, length.out = n_colors + 1),
  annotation_col    = ann_norm,
  annotation_row    = ann_norm,
  annotation_colors = ann_colors,
  display_numbers   = TRUE,
  number_format     = "%.2f",
  number_color      = "black",
  fontsize_number   = 10,
  fontsize          = 12,
  fontsize_row      = 11,
  fontsize_col      = 11,
  border_color      = "grey80",
  cellwidth         = 60,
  cellheight        = 60,
  treeheight_row    = 30,
  treeheight_col    = 30,
  main = "RUV-Normalized Intensities (Median per Timepoint)\nCXCR7 — Validation vs Initial"
)

dev.off()

###############################################################
# 8) PLOT — RAW (averaged)
###############################################################

png(filename = "analysis/PCA/correlation_heatmap_cxcr7_RAW_averaged.png",
    width = 6 * 300, height = 5.5 * 300, res = 300)

pheatmap(
  cor_avg_raw,
  clustering_method = "complete",
  color  = palette_div,
  breaks = seq(-1, 1, length.out = n_colors + 1),
  annotation_col    = ann_raw,
  annotation_row    = ann_raw,
  annotation_colors = ann_colors,
  display_numbers   = TRUE,
  number_format     = "%.2f",
  number_color      = "black",
  fontsize_number   = 10,
  fontsize          = 12,
  fontsize_row      = 11,
  fontsize_col      = 11,
  border_color      = "grey80",
  cellwidth         = 60,
  cellheight        = 60,
  treeheight_row    = 30,
  treeheight_col    = 30,
  main = "RAW Intensities (Median per Timepoint)\nCXCR7 — Validation vs Initial"
)

dev.off()

cat("\n✓ correlation_heatmap_cxcr7_NORM_averaged.png\n")
cat("✓ correlation_heatmap_cxcr7_RAW_averaged.png\n")

# Print the matrices
cat("\n=== NORMALIZED (averaged, z-scored) ===\n")
print(round(cor_avg_norm, 3))

cat("\n=== RAW (averaged, z-scored) ===\n")
print(round(cor_avg_raw, 3))






###############################################################
# BETTER APPROACH: Per-timepoint scatter correlation
# For each timepoint, correlate Val median vs Init median
# across all shared phosphosites
###############################################################

# Reuse val_avg_norm / init_avg_norm (before z-scoring!)
# These are: rows = phosphosites, cols = timepoints (CXCR7.10s, CXCR7.600s, CXCR7.1800s)

shared_ids_both <- intersect(rownames(val_avg_norm), rownames(init_avg_norm))

# Also for raw
shared_ids_raw_both <- intersect(rownames(val_avg_raw), rownames(init_avg_raw))

###############################################################
# 1) Per-timepoint Pearson across phosphosites
###############################################################

timepoints <- c("CXCR7.10s", "CXCR7.600s", "CXCR7.1800s")

cat("\n=== NORMALIZED: Per-timepoint correlation (across phosphosites) ===\n")
cor_norm_tp <- sapply(timepoints, function(tp) {
  cor(val_avg_norm[shared_ids_both, tp],
      init_avg_norm[shared_ids_both, tp],
      use = "pairwise.complete.obs")
})
names(cor_norm_tp) <- c("10s", "600s", "1800s")
print(round(cor_norm_tp, 4))

cat("\n=== RAW: Per-timepoint correlation (across phosphosites) ===\n")
cor_raw_tp <- sapply(timepoints, function(tp) {
  cor(val_avg_raw[shared_ids_raw_both, tp],
      init_avg_raw[shared_ids_raw_both, tp],
      use = "pairwise.complete.obs")
})
names(cor_raw_tp) <- c("10s", "600s", "1800s")
print(round(cor_raw_tp, 4))

###############################################################
# 2) Publication figure: scatter plots per timepoint
###############################################################

library(ggplot2)
library(patchwork)

make_scatter <- function(val_mat, init_mat, shared, tp, tp_label, data_type) {
  df <- data.frame(
    Val  = val_mat[shared, tp],
    Init = init_mat[shared, tp]
  )
  r <- cor(df$Val, df$Init, use = "pairwise.complete.obs")
  
  ggplot(df, aes(Val, Init)) +
    geom_point(alpha = 0.15, size = 0.8, color = "#2C7BB6") +
    geom_smooth(method = "lm", color = "#D7191C", se = FALSE, linewidth = 0.8) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
             label = paste0("r = ", round(r, 3)),
             size = 4.5, fontface = "bold") +
    theme_bw(base_size = 12) +
    labs(title = paste0(tp_label, " (", data_type, ")"),
         x = "Validation (median intensity)",
         y = "Initial (median intensity)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11))
}

# Normalized scatter plots
p1n <- make_scatter(val_avg_norm, init_avg_norm, shared_ids_both, "CXCR7.10s",   "10s",   "Normalized")
p2n <- make_scatter(val_avg_norm, init_avg_norm, shared_ids_both, "CXCR7.600s",  "600s",  "Normalized")
p3n <- make_scatter(val_avg_norm, init_avg_norm, shared_ids_both, "CXCR7.1800s", "1800s", "Normalized")

# Raw scatter plots
p1r <- make_scatter(val_avg_raw, init_avg_raw, shared_ids_raw_both, "CXCR7.10s",   "10s",   "Raw")
p2r <- make_scatter(val_avg_raw, init_avg_raw, shared_ids_raw_both, "CXCR7.600s",  "600s",  "Raw")
p3r <- make_scatter(val_avg_raw, init_avg_raw, shared_ids_raw_both, "CXCR7.1800s", "1800s", "Raw")

# Combine: top row = Normalized, bottom row = Raw
combined_plot <- (p1n | p2n | p3n) / (p1r | p2r | p3r) +
  plot_annotation(
    title = "CXCR7 Intensity Correlation: Validation vs Initial",
    subtitle = "Top: RUV-Normalized | Bottom: Raw",
    theme = theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11)
    )
  )

ggsave("analysis/PCA/scatter_correlation_val_vs_init_CXCR7.png",
       combined_plot, width = 14, height = 9, dpi = 300)

cat("\n✓ scatter_correlation_val_vs_init_CXCR7.png\n")








