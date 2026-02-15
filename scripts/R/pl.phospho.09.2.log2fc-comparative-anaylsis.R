###############################################################
## COMBINED PIPELINE FOR NEW + OLD PHOSPHOPROTEOMICS DATASET
## Harmonization • Collapse by Site • Build Unified Input
###############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
})

###############################################################
# 1) HARMONIZATION FUNCTIONS
###############################################################

## ------------------------------------------------------------
## A) Harmonize NEW validation datasets
## ------------------------------------------------------------
harmonize_new <- function(df) {
  df %>%
    dplyr::rename(
      uniprot_id = uniprot,
      name       = symbol,
      Average    = AveExpr,
      PValue     = adj.P.Val,
      PSite      = psite
    ) %>%
    dplyr::select(uniprot_id, name, Average, logFC, PValue, PSite) %>%
    as.data.frame()
}


## ------------------------------------------------------------
## B) Harmonize OLD initial datasets  
##    structure from df$id: "UNIPROT;GENE;PSITE;PEPTIDE;INDEX"
## ------------------------------------------------------------
harmonize_old <- function(df) {
  
  parts <- strsplit(df$id, ";")
  
  df$uniprot_id <- sapply(parts, `[`, 1)
  df$name       <- sapply(parts, `[`, 2)
  df$PSite      <- sapply(parts, `[`, 3)
  df$Peptide    <- sapply(parts, `[`, 4)
  
  df %>%
    dplyr::transmute(
      uniprot_id,
      name,
      Average = AveExpr,
      logFC,
      PValue = P.Value,
      PSite
    ) %>%
    as.data.frame()
}


## ------------------------------------------------------------
## C) Collapse AFTER intersection (not used yet here)
##    - group by uniprot_id
##    - select phosphosite with highest |logFC|
## ------------------------------------------------------------
collapse_uniprot <- function(df) {
  df %>%
    mutate(abs_logFC = abs(logFC)) %>%
    group_by(uniprot_id) %>%
    slice_max(abs_logFC, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(uniprot_id, name, PSite, Average, logFC, PValue)
}



###############################################################
# 2) LOAD & HARMONIZE ALL DATA
###############################################################

## ----------------------------  
## NEW validation datasets
## ----------------------------

dfs_new <- list(
  top.10, top.600, top.1800,
  top.10.dmso.vs.0s, top.600.dmso.vs.0s, top.1800.dmso.vs.0s,
  top.10.cxcr7.vs.0s, top.600.cxcr7.vs.0s, top.1800.cxcr7.vs.0s
)

names(dfs_new) <- c(
  "val_10.dmso.vs.cxcr7", "val_600.dmso.vs.cxcr7", "val_1800.dmso.vs.cxcr7",
  "val_10.dmso.vs.0s",    "val_600.dmso.vs.0s",    "val_1800.dmso.vs.0s",
  "val_10.cxcr7.vs.0s",   "val_600.cxcr7.vs.0s",   "val_1800.cxcr7.vs.0s"
)

dfs_new_raw <- lapply(dfs_new, harmonize_new)


## ----------------------------  
## OLD initial datasets
## ----------------------------

old_path  <- "SubProjects/CXCR7_initial/phosphoproteomics/data/processed_data"
old_files <- list.files(old_path, pattern = "cxcr7", full.names = TRUE)

old_tables_raw        <- lapply(old_files, read.delim)
old_tables_harmonized <- lapply(old_tables_raw, harmonize_old)

names(old_tables_harmonized) <- c(
  "init_10.cxcr7.vs.0s","init_30.cxcr7.vs.0s","init_60.cxcr7.vs.0s",
  "init_300.cxcr7.vs.0s","init_600.cxcr7.vs.0s","init_900.cxcr7.vs.0s",
  "init_1800.cxcr7.vs.0s"
)



###############################################################
# 3) INTERSECT PHOSPHOSITES ACROSS ALL DATASETS
###############################################################

## ----------------------------  
## A) Add true phosphosite key (uniprot_id + PSite)  
## ----------------------------
add_key <- function(df) {
  df$phospho_id <- paste(df$uniprot_id, df$PSite, sep = "@")
  df
}

new_keyed  <- lapply(dfs_new_raw, add_key)
init_keyed <- lapply(old_tables_harmonized, add_key)


## ----------------------------  
## B) Compute intersection across ALL datasets  
## ----------------------------
common_phosphosites <- Reduce(
  intersect,
  lapply(c(new_keyed, init_keyed), function(df) df$phospho_id)
)

cat("Number of common phosphosites:", length(common_phosphosites), "\n")


## ----------------------------  
## C) Filter each dataset to intersected sites  
## ----------------------------
filter_common <- function(df) {
  df[df$phospho_id %in% common_phosphosites, ]
}

dfs_new_intersect  <- lapply(new_keyed,  filter_common)
init_intersect     <- lapply(init_keyed, filter_common)



###############################################################
# 4) BUILD FINAL all_inputs (INTERSECTED, NOT COLLAPSED)
###############################################################

## Combine validation + initial intersected tables
all_inputs <- c(dfs_new_intersect, init_intersect)

## Sort each dataset by Uniprot + PSite for consistency
all_inputs <- lapply(all_inputs, function(df) {
  df[order(df$uniprot_id, df$PSite), ]
})

## Report
cat("\nFinal all_inputs created:\n")
cat("Validation datasets:", length(dfs_new_intersect), "\n")
cat("Initial datasets:   ", length(init_intersect), "\n")
cat("Total tables:       ", length(all_inputs), "\n\n")

## Optional: show dimensions of each table
for (nm in names(all_inputs)) {
  cat(sprintf("%s → %d rows, %d columns\n",
              nm, nrow(all_inputs[[nm]]), ncol(all_inputs[[nm]])))
}



###############################################################
# 5) COLLAPSE AFTER INTERSECTION (UNIPROT LEVEL)
###############################################################
###############################################################
# 1) SELECT PSite per Uniprot using ONLY val_dmso_vs_cxcr7 triplet
###############################################################

val_triplet <- c("val_10.cxcr7.vs.0s",
                 "val_600.cxcr7.vs.0s",
                 "val_1800.cxcr7.vs.0s")

# merge the three validation datasets
merged_triplet <- dplyr::bind_rows(lapply(val_triplet, function(nm){
  df <- all_inputs[[nm]]
  df$source <- nm
  df
}))

# Choose PSite with strongest regulation
selected_sites <- merged_triplet %>%
  mutate(abs_logFC = abs(logFC)) %>%
  group_by(uniprot_id) %>%
  slice_max(abs_logFC, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(uniprot_id,
            PSite,
            phospho_id = paste(uniprot_id, PSite, sep="@"))

selected_phosphosites <- selected_sites$phospho_id
length(selected_phosphosites)  # should be 907


###############################################################
## 1) Use selected_phosphosites from VAL_dmso_vs_cxcr7 triplet
###############################################################

# selected_phosphosites was produced earlier:
# length(selected_phosphosites)  # should be 907

###############################################################
## 2) Filter all 16 datasets by selected phosphosites
###############################################################

filter_to_selected <- function(df, selected_phosphosites) {
  df$phospho_id <- paste(df$uniprot_id, df$PSite, sep = "@")
  df[df$phospho_id %in% selected_phosphosites, ]
}

all_inputs_filtered <- lapply(all_inputs, function(df)
  filter_to_selected(df, selected_phosphosites)
)

###############################################################
## 3) Collapse EVERY dataset by Uniprot (max |logFC|)
##    → ensures perfect alignment (907 rows)
###############################################################

collapse_uniprot <- function(df) {
  df %>%
    mutate(abs_logFC = abs(logFC)) %>%
    group_by(uniprot_id) %>%
    slice_max(abs_logFC, with_ties = FALSE) %>%
    ungroup() %>%
    select(uniprot_id, name, PSite, Average, logFC, PValue)
}

all_inputs_collapsed <- lapply(all_inputs_filtered, collapse_uniprot)

###############################################################
## 4) Sort consistently
###############################################################

all_inputs_collapsed <- lapply(all_inputs_collapsed, function(df) {
  df[order(df$uniprot_id), ]
})

###############################################################
## 5) Report result
###############################################################

cat("\n=== FINAL COLLAPSED DATASETS ===\n")
for (nm in names(all_inputs_collapsed)) {
  cat(sprintf("%s → %d rows, %d cols\n",
              nm, nrow(all_inputs_collapsed[[nm]]), 
              ncol(all_inputs_collapsed[[nm]])))
}





###############################################################
## VERIFY TRUE OVERLAP AFTER COLLAPSING
###############################################################

# reconstruct the uniprot@PSite identifier
add_key <- function(df) {
  paste(df$uniprot_id, df$PSite, sep="@")
}

collapsed_keys <- lapply(all_inputs_collapsed, add_key)

# compare overlap across all 16 datasets
pairwise_overlap <- combn(names(collapsed_keys), 2, function(idx) {
  n <- length(intersect(collapsed_keys[[idx[1]]], collapsed_keys[[idx[2]]]))
  sprintf("%s <-> %s :  %d sites", idx[1], idx[2], n)
})

cat("\n=== Overlap of collapsed datasets (should all be 907) ===\n")
cat(paste(pairwise_overlap, collapse="\n"))







###############################################################
# BUILD UNIFIED LOGFC MATRIX (907 phosphosites × 16 datasets)
###############################################################

# extract logFC columns only
logfc_matrix <- do.call(cbind, lapply(all_inputs_collapsed, function(df) df$logFC))

# name columns by dataset names
colnames(logfc_matrix) <- names(all_inputs_collapsed)

# add rownames = uniprot_id@PSite
rownames(logfc_matrix) <- paste(
  all_inputs_collapsed[[1]]$uniprot_id,
  all_inputs_collapsed[[1]]$PSite,
  sep = "@"
)

dim(logfc_matrix)   # should be 907 × 16
head(logfc_matrix)


###############################################################
# UMAP
###############################################################


library(uwot)
set.seed(1)
umap_res <- umap(t(logfc_matrix), n_neighbors = 10, min_dist = 0.3)

# uwot returns a matrix directly, not $layout
plot(umap_res[,1], umap_res[,2],
     col = "steelblue", pch = 19, cex = 2,
     main = "UMAP of 16 phosphoproteomics profiles",
     xlab = "UMAP-1", ylab = "UMAP-2")

text(umap_res[,1], umap_res[,2],
     labels = colnames(logfc_matrix),
     pos = 3, cex = 0.7)


###############################################################
# heatmap correlation
###############################################################


library(pheatmap)

cor_mat <- cor(logfc_matrix, method = "pearson")

pheatmap(cor_mat,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "Correlation of 16 logFC profiles",
         border_color = NA,
         color = colorRampPalette(c("navy", "white", "firebrick"))(200))



###############################################################
## Build top psites
###############################################################

###############################################################
## Compute phosphosite statistics across all contrasts
###############################################################
stats_df <- data.frame(
  phospho_id = rownames(logfc_matrix),
  mean_abs_logFC = rowMeans(abs(logfc_matrix)),
  max_abs_logFC = apply(abs(logfc_matrix), 1, max),
  sd_logFC = apply(logfc_matrix, 1, sd)
)

###############################################################
## Select the top N strongest phosphosites
###############################################################
top_n <- 100   # change to 30, 50, 100, etc.

top_sites <- stats_df[order(-stats_df$max_abs_logFC), ][1:top_n, ]
top_ids <- top_sites$phospho_id
length(top_ids)


###############################################################
## Subset logFC matrix to top N dynamic phosphosites
###############################################################
logfc_top <- logfc_matrix[top_ids, ]
dim(logfc_top)
# e.g., 200 × 16



###############################################################
## Correlation heatmap using only top sites
###############################################################
library(pheatmap)

cor_mat_top <- cor(logfc_top, method = "pearson")


dev.new()

pheatmap(
  cor_mat_top,
  clustering_method = "complete",
  color = colorRampPalette(c("blue","white","red"))(200),
  main = paste("Correlation of 16 contrasts using top", top_n, "phosphosites"),
  fontsize_row = 10,
  fontsize_col = 10
)



###############################################################
## UMAP using top phosphosites only
###############################################################
library(uwot)

dev.new()
set.seed(1)
umap_res_top <- umap(t(logfc_top), n_neighbors = 10, min_dist = 0.3)

plot(
  umap_res_top[,1],
  umap_res_top[,2],
  pch = 19, cex = 1.5,
  col = "darkred",
  main = paste("UMAP of 16 contrasts (top", top_n, "sites)"),
  xlab = "UMAP-1", ylab = "UMAP-2"
)

text(
  umap_res_top[,1],
  umap_res_top[,2],
  labels = colnames(logfc_top),
  pos = 3, cex = 0.7
)




###############################################################
## Publication-Ready logFC Correlation Heatmap
## CXCR7 vs 0s — All timepoints, Val + Init
## 907 shared phosphosites
###############################################################

library(pheatmap)

###############################################################
# 1) SELECT ONLY CXCR7 vs 0s CONTRASTS
###############################################################

cxcr7_vs_0s <- c(
  "val_10.cxcr7.vs.0s",
  "val_600.cxcr7.vs.0s",
  "val_1800.cxcr7.vs.0s",
  "init_10.cxcr7.vs.0s",
  "init_30.cxcr7.vs.0s",
  "init_60.cxcr7.vs.0s",
  "init_300.cxcr7.vs.0s",
  "init_600.cxcr7.vs.0s",
  "init_900.cxcr7.vs.0s",
  "init_1800.cxcr7.vs.0s"
)

logfc_cxcr7 <- logfc_matrix[, cxcr7_vs_0s]
cor_cxcr7 <- cor(logfc_cxcr7, method = "pearson")

###############################################################
# 2) PRETTY LABELS
###############################################################

pretty_labels <- c(
  "val_10.cxcr7.vs.0s"     = "Val 10s",
  "val_600.cxcr7.vs.0s"    = "Val 600s",
  "val_1800.cxcr7.vs.0s"   = "Val 1800s",
  "init_10.cxcr7.vs.0s"    = "Init 10s",
  "init_30.cxcr7.vs.0s"    = "Init 30s",
  "init_60.cxcr7.vs.0s"    = "Init 60s",
  "init_300.cxcr7.vs.0s"   = "Init 300s",
  "init_600.cxcr7.vs.0s"   = "Init 600s",
  "init_900.cxcr7.vs.0s"   = "Init 900s",
  "init_1800.cxcr7.vs.0s"  = "Init 1800s"
)

rownames(cor_cxcr7) <- pretty_labels[rownames(cor_cxcr7)]
colnames(cor_cxcr7) <- pretty_labels[colnames(cor_cxcr7)]

###############################################################
# 3) ANNOTATIONS
###############################################################

annotation_df <- data.frame(
  Dataset = factor(
    ifelse(grepl("^Val", pretty_labels[cxcr7_vs_0s]), "Validation", "Initial"),
    levels = c("Validation", "Initial")
  ),
  row.names = pretty_labels[cxcr7_vs_0s]
)

ann_colors <- list(
  Dataset = c("Validation" = "#2C7BB6", "Initial" = "#D7191C")
)

###############################################################
# 4) COLOR PALETTE
###############################################################

n_colors <- 200
palette_div <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE",
                                  "#D1E5F0", "#F7F7F7",
                                  "#FDDBC7", "#F4A582",
                                  "#D6604D", "#B2182B"))(n_colors)

breaks <- seq(-1, 1, length.out = n_colors + 1)

###############################################################
# 5) PLOT
###############################################################

png(filename = "analysis/PCA/correlation_heatmap_logFC_cxcr7_vs_0s.png",
    width = 8 * 300, height = 7 * 300, res = 300)

pheatmap(
  cor_cxcr7,
  clustering_method        = "complete",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  color  = palette_div,
  breaks = breaks,
  annotation_col    = annotation_df,
  annotation_row    = annotation_df,
  annotation_colors = ann_colors,
  display_numbers   = TRUE,
  number_format     = "%.2f",
  number_color      = "black",
  fontsize_number   = 8,
  fontsize          = 12,
  fontsize_row      = 11,
  fontsize_col      = 11,
  border_color      = "grey80",
  cellwidth         = 50,
  cellheight        = 50,
  treeheight_row    = 30,
  treeheight_col    = 30,
  main = "CXCR7 vs Control — logFC Correlation\n907 Shared Phosphosites"
)

dev.off()

cat("\n✓ correlation_heatmap_logFC_cxcr7_vs_0s.png\n")
print(round(cor_cxcr7, 3))








###############################################################
## logFC Correlation — Filtered by significance & effect size
###############################################################

library(pheatmap)

###############################################################
# 1) BUILD P-VALUE MATRIX (same structure as logfc_matrix)
###############################################################

pval_matrix <- do.call(cbind, lapply(all_inputs_collapsed, function(df) df$PValue))
colnames(pval_matrix) <- names(all_inputs_collapsed)
rownames(pval_matrix) <- rownames(logfc_matrix)

###############################################################
# 2) DEFINE FILTERS
###############################################################

cxcr7_vs_0s <- c(
  "val_10.cxcr7.vs.0s", "val_600.cxcr7.vs.0s", "val_1800.cxcr7.vs.0s",
  "init_10.cxcr7.vs.0s", "init_30.cxcr7.vs.0s", "init_60.cxcr7.vs.0s",
  "init_300.cxcr7.vs.0s", "init_600.cxcr7.vs.0s", "init_900.cxcr7.vs.0s",
  "init_1800.cxcr7.vs.0s"
)

logfc_cxcr7 <- logfc_matrix[, cxcr7_vs_0s]
pval_cxcr7  <- pval_matrix[, cxcr7_vs_0s]

# Filter A: p < 0.05 in at least one contrast
sig_any <- apply(pval_cxcr7, 1, function(r) any(r < 0.05, na.rm = TRUE))

# Filter B: |logFC| > 0.5 in at least one contrast
fc_any <- apply(abs(logfc_cxcr7), 1, function(r) any(r > 0.5, na.rm = TRUE))

# Filter C: BOTH sig AND strong FC
both <- sig_any & fc_any

cat("All phosphosites:        ", nrow(logfc_cxcr7), "\n")
cat("p < 0.05 (any contrast): ", sum(sig_any), "\n")
cat("|logFC| > 0.5 (any):     ", sum(fc_any), "\n")
cat("Both p<0.05 & |FC|>0.5:  ", sum(both), "\n")

###############################################################
# 3) HELPER: make publication heatmap
###############################################################

pretty_labels <- c(
  "val_10.cxcr7.vs.0s"    = "Val 10s",
  "val_600.cxcr7.vs.0s"   = "Val 600s",
  "val_1800.cxcr7.vs.0s"  = "Val 1800s",
  "init_10.cxcr7.vs.0s"   = "Init 10s",
  "init_30.cxcr7.vs.0s"   = "Init 30s",
  "init_60.cxcr7.vs.0s"   = "Init 60s",
  "init_300.cxcr7.vs.0s"  = "Init 300s",
  "init_600.cxcr7.vs.0s"  = "Init 600s",
  "init_900.cxcr7.vs.0s"  = "Init 900s",
  "init_1800.cxcr7.vs.0s" = "Init 1800s"
)

annotation_df <- data.frame(
  Dataset = factor(
    ifelse(grepl("^Val", pretty_labels[cxcr7_vs_0s]), "Validation", "Initial"),
    levels = c("Validation", "Initial")
  ),
  row.names = pretty_labels[cxcr7_vs_0s]
)

ann_colors <- list(
  Dataset = c("Validation" = "#2C7BB6", "Initial" = "#D7191C")
)

n_colors <- 200
palette_div <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE",
                                  "#D1E5F0", "#F7F7F7",
                                  "#FDDBC7", "#F4A582",
                                  "#D6604D", "#B2182B"))(n_colors)

plot_cor_heatmap <- function(logfc_sub, filter_name, n_sites, filename) {
  
  cor_sub <- cor(logfc_sub, method = "pearson")
  rownames(cor_sub) <- pretty_labels[colnames(logfc_sub)]
  colnames(cor_sub) <- pretty_labels[colnames(logfc_sub)]
  
  png(filename = filename, width = 8 * 300, height = 7 * 300, res = 300)
  
  pheatmap(
    cor_sub,
    clustering_method        = "complete",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    color  = palette_div,
    breaks = seq(-1, 1, length.out = n_colors + 1),
    annotation_col    = annotation_df,
    annotation_row    = annotation_df,
    annotation_colors = ann_colors,
    display_numbers   = TRUE,
    number_format     = "%.2f",
    number_color      = "black",
    fontsize_number   = 8,
    fontsize          = 12,
    fontsize_row      = 11,
    fontsize_col      = 11,
    border_color      = "grey80",
    cellwidth         = 50,
    cellheight        = 50,
    treeheight_row    = 30,
    treeheight_col    = 30,
    main = paste0("CXCR7 vs Control — logFC Correlation\n",
                  filter_name, " (n = ", n_sites, " phosphosites)")
  )
  
  dev.off()
  
  cat("\n===", filter_name, "===\n")
  print(round(cor_sub, 3))
}

###############################################################
# 4) GENERATE ALL HEATMAPS
###############################################################

# A) p < 0.05
plot_cor_heatmap(
  logfc_cxcr7[sig_any, ],
  "p < 0.05", sum(sig_any),
  "analysis/PCA/correlation_logFC_cxcr7_pval05.png"
)

# B) |logFC| > 0.5
plot_cor_heatmap(
  logfc_cxcr7[fc_any, ],
  "|logFC| > 0.5", sum(fc_any),
  "analysis/PCA/correlation_logFC_cxcr7_fc05.png"
)

# C) Both
plot_cor_heatmap(
  logfc_cxcr7[both, ],
  "p < 0.05 & |logFC| > 0.5", sum(both),
  "analysis/PCA/correlation_logFC_cxcr7_sig_and_fc.png"
)

cat("\n✓ All three heatmaps saved\n")






###############################################################
# Timepoint-specific cross-experiment correlation
# Raw vs Normalized — matching CXCR7 timepoints
###############################################################

timepoints <- c("10s", "600s", "1800s")

# Helper: get median per timepoint
get_median_tp <- function(mat, tp) {
  cols <- grep(paste0("CXCR7\\.", tp, "_"), colnames(mat), value = TRUE)
  if (length(cols) == 0) return(NULL)
  rowMedians <- apply(mat[, cols, drop = FALSE], 1, median, na.rm = TRUE)
  rowMedians
}

###############################################################
# NORMALIZED
###############################################################
shared_norm <- intersect(rownames(val_norm_clean), rownames(init_norm_clean))

cat("\n=== NORMALIZED (n =", length(shared_norm), "phosphosites) ===\n")
for (tp in timepoints) {
  v <- get_median_tp(val_norm_clean[shared_norm, ], tp)
  i <- get_median_tp(init_norm_clean[shared_norm, ], tp)
  r <- cor(v, i, use = "pairwise.complete.obs")
  cat(sprintf("  Val %s ↔ Init %s:  r = %.3f\n", tp, tp, r))
}

###############################################################
# RAW
###############################################################
shared_raw <- intersect(rownames(val_raw_clean), rownames(init_raw_clean))

cat("\n=== RAW (n =", length(shared_raw), "phosphosites) ===\n")
for (tp in timepoints) {
  v <- get_median_tp(val_raw_clean[shared_raw, ], tp)
  i <- get_median_tp(init_raw_clean[shared_raw, ], tp)
  r <- cor(v, i, use = "pairwise.complete.obs")
  cat(sprintf("  Val %s ↔ Init %s:  r = %.3f\n", tp, tp, r))
}














###############################################################
# REPRODUCE FIGURE S5 C+D
# Overlap of significant phosphosites between experiments
###############################################################

library(ggplot2)
library(dplyr)
library(tidyr)

###############################################################
# 1) GET SIGNIFICANT PHOSPHOSITES PER TIMEPOINT
#    Criteria: p < 0.05 AND |logFC| > 0.5
###############################################################

# Build combined data with logFC and pvalue
get_sig_sites <- function(dataset_name) {
  df <- all_inputs_collapsed[[dataset_name]]
  df$phospho_id <- paste(df$uniprot_id, df$PSite, sep = "@")
  sig <- df %>% filter(PValue < 0.05, abs(logFC) > 0.5)
  list(
    up   = sig$phospho_id[sig$logFC > 0],
    down = sig$phospho_id[sig$logFC < 0],
    all  = sig$phospho_id
  )
}

# Matched timepoints
timepoints <- c("10", "600", "1800")

init_names <- paste0("init_", timepoints, ".cxcr7.vs.0s")
val_names  <- paste0("val_", timepoints, ".cxcr7.vs.0s")

# Get sig sites per experiment and timepoint
init_sig <- lapply(setNames(init_names, timepoints), get_sig_sites)
val_sig  <- lapply(setNames(val_names, timepoints), get_sig_sites)

###############################################################
# 2) COUNT: total and overlapping
###############################################################

cat("\n=== INITIAL EXPERIMENT (Panel C perspective) ===\n")
cat("Timepoint | Init UP | Init DOWN | Shared UP | Shared DOWN\n")
cat("----------|---------|-----------|-----------|------------\n")

for (tp in timepoints) {
  i_up   <- init_sig[[tp]]$up
  i_down <- init_sig[[tp]]$down
  v_up   <- val_sig[[tp]]$up
  v_down <- val_sig[[tp]]$down
  
  shared_up   <- length(intersect(i_up, v_up))
  shared_down <- length(intersect(i_down, v_down))
  
  cat(sprintf("%7ss   | %7d | %9d | %9d | %11d\n",
              tp, length(i_up), length(i_down), shared_up, shared_down))
}

cat("\n=== VALIDATION EXPERIMENT (Panel D perspective) ===\n")
cat("Timepoint | Val UP | Val DOWN | Shared UP | Shared DOWN\n")
cat("----------|--------|----------|-----------|------------\n")

for (tp in timepoints) {
  i_up   <- init_sig[[tp]]$up
  i_down <- init_sig[[tp]]$down
  v_up   <- val_sig[[tp]]$up
  v_down <- val_sig[[tp]]$down
  
  shared_up   <- length(intersect(i_up, v_up))
  shared_down <- length(intersect(i_down, v_down))
  
  cat(sprintf("%7ss   | %6d | %8d | %9d | %11d\n",
              tp, length(v_up), length(v_down), shared_up, shared_down))
}

###############################################################
# 3) BUILD PLOT DATA
###############################################################

build_bar_data <- function(experiment_sig, other_sig, panel_label) {
  rows <- list()
  for (tp in timepoints) {
    e_up   <- experiment_sig[[tp]]$up
    e_down <- experiment_sig[[tp]]$down
    o_up   <- other_sig[[tp]]$up
    o_down <- other_sig[[tp]]$down
    
    rows[[length(rows) + 1]] <- data.frame(
      timepoint   = paste0(tp, "s"),
      direction   = "up",
      total       = length(e_up),
      shared      = length(intersect(e_up, o_up)),
      panel       = panel_label
    )
    rows[[length(rows) + 1]] <- data.frame(
      timepoint   = paste0(tp, "s"),
      direction   = "down",
      total       = -length(e_down),      # negative for downregulated
      shared      = -length(intersect(e_down, o_down)),
      panel       = panel_label
    )
  }
  do.call(rbind, rows)
}

df_C <- build_bar_data(init_sig, val_sig, "C) Initial Experiment")
df_D <- build_bar_data(val_sig, init_sig, "D) Validation Experiment")

df_all <- rbind(df_C, df_D)
df_all$timepoint <- factor(df_all$timepoint, levels = c("10s", "600s", "1800s"))

###############################################################
# 4) PLOT — matching the style of Figure S5 C+D
###############################################################

plot_panel <- function(df, title) {
  ggplot(df, aes(x = timepoint)) +
    # Light bars = total significant
    geom_col(aes(y = total), fill = ifelse(df$direction == "up", "#FABCAC", "#ABC8E8"),
             width = 0.65) +
    # Solid bars = shared with other experiment
    geom_col(aes(y = shared), fill = ifelse(df$direction == "up", "#D32F2F", "#1565C0"),
             width = 0.65) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    scale_y_continuous(limits = c(-450, 450),
                       breaks = seq(-400, 400, 200)) +
    labs(x = "Time", y = "", title = title) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0, face = "bold", size = 14),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text = element_text(size = 11),
      plot.margin = margin(10, 15, 10, 10)
    )
}

p_C <- plot_panel(df_all[df_all$panel == "C) Initial Experiment", ],
                  "C")
p_D <- plot_panel(df_all[df_all$panel == "D) Validation Experiment", ],
                  "D")

library(patchwork)
combined <- p_C | p_D

ggsave("analysis/PCA/FigS5_CD_logFC_overlap.png",
       combined, width = 12, height = 5.5, dpi = 300)

cat("\n✓ FigS5_CD_logFC_overlap.png saved\n")









###############################################################
# PRINT NUMBERS — COLLAPSED (907 proteins)
###############################################################

cat("\n=== COLLAPSED DATA (907 proteins, p<0.05 & |logFC|>0.5) ===\n\n")

cat("--- Panel C: Initial Experiment ---\n")
cat("Timepoint | Init UP | Init DOWN | Shared UP | Shared DOWN\n")
for (tp in timepoints) {
  i_up   <- init_sig[[tp]]$up
  i_down <- init_sig[[tp]]$down
  v_up   <- val_sig[[tp]]$up
  v_down <- val_sig[[tp]]$down
  
  cat(sprintf("  %5ss  |  %5d  |  %7d  |  %7d  |  %9d\n",
              tp, length(i_up), length(i_down),
              length(intersect(i_up, v_up)),
              length(intersect(i_down, v_down))))
}

cat("\n--- Panel D: Validation Experiment ---\n")
cat("Timepoint | Val UP | Val DOWN | Shared UP | Shared DOWN\n")
for (tp in timepoints) {
  i_up   <- init_sig[[tp]]$up
  i_down <- init_sig[[tp]]$down
  v_up   <- val_sig[[tp]]$up
  v_down <- val_sig[[tp]]$down
  
  cat(sprintf("  %5ss  |  %5d  |  %7d  |  %7d  |  %9d\n",
              tp, length(v_up), length(v_down),
              length(intersect(i_up, v_up)),
              length(intersect(i_down, v_down))))
}

###############################################################
# NOW TRY WITH UNCOLLAPSED DATA (1896 phosphosites)
###############################################################

get_sig_sites_uncollapsed <- function(dataset_name) {
  df <- all_inputs[[dataset_name]]  # NOT collapsed!
  df$phospho_id <- paste(df$uniprot_id, df$PSite, sep = "@")
  sig <- df %>% filter(PValue < 0.05, abs(logFC) > 0.5)
  list(
    up   = sig$phospho_id[sig$logFC > 0],
    down = sig$phospho_id[sig$logFC < 0],
    all  = sig$phospho_id
  )
}

init_sig_unc <- lapply(setNames(init_names, timepoints), get_sig_sites_uncollapsed)
val_sig_unc  <- lapply(setNames(val_names, timepoints), get_sig_sites_uncollapsed)

cat("\n\n=== UNCOLLAPSED DATA (1896 phosphosites, p<0.05 & |logFC|>0.5) ===\n\n")

cat("--- Panel C: Initial Experiment ---\n")
cat("Timepoint | Init UP | Init DOWN | Shared UP | Shared DOWN\n")
for (tp in timepoints) {
  i_up   <- init_sig_unc[[tp]]$up
  i_down <- init_sig_unc[[tp]]$down
  v_up   <- val_sig_unc[[tp]]$up
  v_down <- val_sig_unc[[tp]]$down
  
  cat(sprintf("  %5ss  |  %5d  |  %7d  |  %7d  |  %9d\n",
              tp, length(i_up), length(i_down),
              length(intersect(i_up, v_up)),
              length(intersect(i_down, v_down))))
}

cat("\n--- Panel D: Validation Experiment ---\n")
cat("Timepoint | Val UP | Val DOWN | Shared UP | Shared DOWN\n")
for (tp in timepoints) {
  i_up   <- init_sig_unc[[tp]]$up
  i_down <- init_sig_unc[[tp]]$down
  v_up   <- val_sig_unc[[tp]]$up
  v_down <- val_sig_unc[[tp]]$down
  
  cat(sprintf("  %5ss  |  %5d  |  %7d  |  %7d  |  %9d\n",
              tp, length(v_up), length(v_down),
              length(intersect(i_up, v_up)),
              length(intersect(i_down, v_down))))
}

###############################################################
# ALSO TRY: p<0.05 only (no FC cutoff)
###############################################################

get_sig_ponly <- function(dataset_name, collapsed = FALSE) {
  df <- if (collapsed) all_inputs_collapsed[[dataset_name]] else all_inputs[[dataset_name]]
  df$phospho_id <- paste(df$uniprot_id, df$PSite, sep = "@")
  sig <- df %>% filter(PValue < 0.05)
  list(
    up   = sig$phospho_id[sig$logFC > 0],
    down = sig$phospho_id[sig$logFC < 0],
    all  = sig$phospho_id
  )
}

init_sig_ponly <- lapply(setNames(init_names, timepoints), get_sig_ponly, collapsed = FALSE)
val_sig_ponly  <- lapply(setNames(val_names, timepoints), get_sig_ponly, collapsed = FALSE)

cat("\n\n=== UNCOLLAPSED, p<0.05 ONLY (no FC cutoff) ===\n\n")

cat("--- Panel C: Initial Experiment ---\n")
cat("Timepoint | Init UP | Init DOWN | Shared UP | Shared DOWN\n")
for (tp in timepoints) {
  i_up   <- init_sig_ponly[[tp]]$up
  i_down <- init_sig_ponly[[tp]]$down
  v_up   <- val_sig_ponly[[tp]]$up
  v_down <- val_sig_ponly[[tp]]$down
  
  cat(sprintf("  %5ss  |  %5d  |  %7d  |  %7d  |  %9d\n",
              tp, length(i_up), length(i_down),
              length(intersect(i_up, v_up)),
              length(intersect(i_down, v_down))))
}

cat("\n--- Panel D: Validation Experiment ---\n")
cat("Timepoint | Val UP | Val DOWN | Shared UP | Shared DOWN\n")
for (tp in timepoints) {
  i_up   <- init_sig_ponly[[tp]]$up
  i_down <- init_sig_ponly[[tp]]$down
  v_up   <- val_sig_ponly[[tp]]$up
  v_down <- val_sig_ponly[[tp]]$down
  
  cat(sprintf("  %5ss  |  %5d  |  %7d  |  %7d  |  %9d\n",
              tp, length(v_up), length(v_down),
              length(intersect(i_up, v_up)),
              length(intersect(i_down, v_down))))
}







###############################################################
# REPLOT with UNCOLLAPSED, p<0.05 only
###############################################################

build_bar_data_v2 <- function(experiment_sig, other_sig, panel_label) {
  rows <- list()
  for (tp in timepoints) {
    e_up   <- experiment_sig[[tp]]$up
    e_down <- experiment_sig[[tp]]$down
    o_up   <- other_sig[[tp]]$up
    o_down <- other_sig[[tp]]$down
    
    rows[[length(rows) + 1]] <- data.frame(
      timepoint = paste0(tp, "s"), direction = "up",
      total = length(e_up),
      shared = length(intersect(e_up, o_up)),
      panel = panel_label
    )
    rows[[length(rows) + 1]] <- data.frame(
      timepoint = paste0(tp, "s"), direction = "down",
      total = -length(e_down),
      shared = -length(intersect(e_down, o_down)),
      panel = panel_label
    )
  }
  do.call(rbind, rows)
}

df_C <- build_bar_data_v2(init_sig_ponly, val_sig_ponly, "C")
df_D <- build_bar_data_v2(val_sig_ponly, init_sig_ponly, "D")

df_all <- rbind(df_C, df_D)
df_all$timepoint <- factor(df_all$timepoint, levels = c("10s", "600s", "1800s"))

plot_panel_v2 <- function(df, title) {
  ggplot(df, aes(x = timepoint)) +
    geom_col(aes(y = total),
             fill = ifelse(df$direction == "up", "#FABCAC", "#ABC8E8"),
             width = 0.65) +
    geom_col(aes(y = shared),
             fill = ifelse(df$direction == "up", "#D32F2F", "#1565C0"),
             width = 0.65) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    scale_y_continuous(limits = c(-450, 450),
                       breaks = seq(-400, 400, 200)) +
    labs(x = "Time", y = "", title = title) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0, face = "bold", size = 14),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text = element_text(size = 11),
      plot.margin = margin(10, 15, 10, 10)
    )
}

p_C <- plot_panel_v2(df_all[df_all$panel == "C", ], "C")
p_D <- plot_panel_v2(df_all[df_all$panel == "D", ], "D")

combined <- p_C | p_D

ggsave("analysis/PCA/FigS5_CD_logFC_overlap_v2.png",
       combined, width = 12, height = 5.5, dpi = 300)

cat("\n✓ FigS5_CD_logFC_overlap_v2.png saved\n")




###############################################################
# IDENTIFY the actual shared phosphosites
###############################################################

# Using adj.P.Val from original data
for (tp in c("10", "600", "1800")) {
  
  init_df <- old_tables_raw[[match(paste0("init_", tp, ".cxcr7.vs.0s"),
                                   names(old_tables_harmonized))]]
  val_df  <- dfs_new[[paste0("val_", tp, ".cxcr7.vs.0s")]]
  
  i_sig <- get_sig_adjp(init_df)
  v_sig <- get_sig_adjp(val_df)
  
  shared_up   <- intersect(i_sig$up, v_sig$up)
  shared_down <- intersect(i_sig$down, v_sig$down)
  
  # Also check: discordant (sig in both but opposite direction)
  discordant_1 <- intersect(i_sig$up, v_sig$down)   # Init UP, Val DOWN
  discordant_2 <- intersect(i_sig$down, v_sig$up)    # Init DOWN, Val UP
  
  cat(sprintf("\n=== %ss ===\n", tp))
  cat(sprintf("  Concordant UP:   %d\n", length(shared_up)))
  cat(sprintf("  Concordant DOWN: %d\n", length(shared_down)))
  cat(sprintf("  Discordant:      %d (Init UP/Val DOWN: %d, Init DOWN/Val UP: %d)\n",
              length(discordant_1) + length(discordant_2),
              length(discordant_1), length(discordant_2)))
  
  if (length(shared_up) > 0) {
    cat("\n  Shared UP phosphosites:\n")
    for (s in sort(shared_up)) cat("    ", s, "\n")
  }
  if (length(shared_down) > 0) {
    cat("\n  Shared DOWN phosphosites:\n")
    for (s in sort(shared_down)) cat("    ", s, "\n")
  }
}








###############################################################
# MAP shared phosphosites to gene symbols
###############################################################

# Build lookup from all_inputs (has both uniprot_id and name)
lookup <- unique(rbind(
  all_inputs[[1]][, c("uniprot_id", "name", "PSite")],
  old_tables_harmonized[[1]][, c("uniprot_id", "name", "PSite")]
))

# Function to convert uniprot@PSite to GeneName_PSite
id_to_gene <- function(phospho_ids) {
  parts <- strsplit(phospho_ids, "@")
  uniprot <- sapply(parts, `[`, 1)
  psite   <- sapply(parts, `[`, 2)
  
  gene <- lookup$name[match(uniprot, lookup$uniprot_id)]
  paste0(gene, "_", psite)
}

# Reprint with gene symbols
for (tp in c("10", "600", "1800")) {
  
  init_df <- old_tables_raw[[match(paste0("init_", tp, ".cxcr7.vs.0s"),
                                   names(old_tables_harmonized))]]
  val_df  <- dfs_new[[paste0("val_", tp, ".cxcr7.vs.0s")]]
  
  i_sig <- get_sig_adjp(init_df)
  v_sig <- get_sig_adjp(val_df)
  
  shared_up   <- sort(intersect(i_sig$up, v_sig$up))
  shared_down <- sort(intersect(i_sig$down, v_sig$down))
  discordant  <- c(intersect(i_sig$up, v_sig$down), intersect(i_sig$down, v_sig$up))
  
  cat(sprintf("\n=== %ss ===\n", tp))
  cat(sprintf("  Concordant UP: %d | Concordant DOWN: %d | Discordant: %d\n",
              length(shared_up), length(shared_down), length(discordant)))
  
  if (length(shared_up) > 0) {
    cat("\n  Shared UP:\n")
    cat("    ", paste(id_to_gene(shared_up), collapse = "\n     "), "\n")
  }
  if (length(shared_down) > 0) {
    cat("\n  Shared DOWN:\n")
    cat("    ", paste(id_to_gene(shared_down), collapse = "\n     "), "\n")
  }
}











###############################################################
# COMPREHENSIVE DIFFERENTIAL EXPRESSION COUNTS
# Init vs Val — CXCR7 vs 0s — with overlap
# Using adj.P.Val from original data
###############################################################

timepoints <- c("10", "600", "1800")

###############################################################
# 1) adj.P.Val < 0.05 only
###############################################################

cat("\n========================================================\n")
cat("  adj.P.Val < 0.05 (no FC cutoff)\n")
cat("========================================================\n\n")

cat("--- Per experiment ---\n")
cat(sprintf("%-12s | %6s | %6s | %6s | %6s\n", 
            "Dataset", "UP", "DOWN", "Total", "% of all"))

for (tp in timepoints) {
  # Init
  nm <- paste0("init_", tp, ".cxcr7.vs.0s")
  df <- old_tables_raw[[match(nm, names(old_tables_harmonized))]]
  i_up   <- sum(df$adj.P.Val < 0.05 & df$logFC > 0, na.rm = TRUE)
  i_down <- sum(df$adj.P.Val < 0.05 & df$logFC < 0, na.rm = TRUE)
  cat(sprintf("Init %5ss   | %6d | %6d | %6d | %5.1f%%\n",
              tp, i_up, i_down, i_up + i_down, 100*(i_up+i_down)/nrow(df)))
}
cat("\n")
for (tp in timepoints) {
  nm <- paste0("val_", tp, ".cxcr7.vs.0s")
  df <- dfs_new[[nm]]
  v_up   <- sum(df$adj.P.Val < 0.05 & df$logFC > 0, na.rm = TRUE)
  v_down <- sum(df$adj.P.Val < 0.05 & df$logFC < 0, na.rm = TRUE)
  cat(sprintf("Val  %5ss   | %6d | %6d | %6d | %5.1f%%\n",
              tp, v_up, v_down, v_up + v_down, 100*(v_up+v_down)/nrow(df)))
}

cat("\n--- Overlap (adj.P.Val < 0.05) ---\n")
cat(sprintf("%-8s | %9s | %11s | %11s | %9s | %s\n",
            "TP", "Conc. UP", "Conc. DOWN", "Total Conc.", "Discord.", "Concordance%"))

for (tp in timepoints) {
  init_df <- old_tables_raw[[match(paste0("init_", tp, ".cxcr7.vs.0s"),
                                   names(old_tables_harmonized))]]
  val_df  <- dfs_new[[paste0("val_", tp, ".cxcr7.vs.0s")]]
  
  i_sig <- get_sig_adjp(init_df)
  v_sig <- get_sig_adjp(val_df)
  
  conc_up   <- length(intersect(i_sig$up, v_sig$up))
  conc_down <- length(intersect(i_sig$down, v_sig$down))
  disc_1    <- length(intersect(i_sig$up, v_sig$down))
  disc_2    <- length(intersect(i_sig$down, v_sig$up))
  total     <- conc_up + conc_down + disc_1 + disc_2
  conc_pct  <- if (total > 0) 100 * (conc_up + conc_down) / total else 0
  
  cat(sprintf("%5ss    | %9d | %11d | %11d | %9d | %5.1f%%\n",
              tp, conc_up, conc_down, conc_up + conc_down,
              disc_1 + disc_2, conc_pct))
}


###############################################################
# 2) adj.P.Val < 0.05 AND |logFC| > 0.5
###############################################################

cat("\n\n========================================================\n")
cat("  adj.P.Val < 0.05 AND |logFC| > 0.5\n")
cat("========================================================\n\n")

cat("--- Per experiment ---\n")
cat(sprintf("%-12s | %6s | %6s | %6s | %6s\n",
            "Dataset", "UP", "DOWN", "Total", "% of all"))

for (tp in timepoints) {
  nm <- paste0("init_", tp, ".cxcr7.vs.0s")
  df <- old_tables_raw[[match(nm, names(old_tables_harmonized))]]
  i_up   <- sum(df$adj.P.Val < 0.05 & df$logFC > 0.5, na.rm = TRUE)
  i_down <- sum(df$adj.P.Val < 0.05 & df$logFC < -0.5, na.rm = TRUE)
  cat(sprintf("Init %5ss   | %6d | %6d | %6d | %5.1f%%\n",
              tp, i_up, i_down, i_up + i_down, 100*(i_up+i_down)/nrow(df)))
}
cat("\n")
for (tp in timepoints) {
  nm <- paste0("val_", tp, ".cxcr7.vs.0s")
  df <- dfs_new[[nm]]
  v_up   <- sum(df$adj.P.Val < 0.05 & df$logFC > 0.5, na.rm = TRUE)
  v_down <- sum(df$adj.P.Val < 0.05 & df$logFC < -0.5, na.rm = TRUE)
  cat(sprintf("Val  %5ss   | %6d | %6d | %6d | %5.1f%%\n",
              tp, v_up, v_down, v_up + v_down, 100*(v_up+v_down)/nrow(df)))
}

cat("\n--- Overlap (adj.P.Val < 0.05 & |FC| > 0.5) ---\n")
cat(sprintf("%-8s | %9s | %11s | %11s | %9s | %s\n",
            "TP", "Conc. UP", "Conc. DOWN", "Total Conc.", "Discord.", "Concordance%"))

for (tp in timepoints) {
  init_df <- old_tables_raw[[match(paste0("init_", tp, ".cxcr7.vs.0s"),
                                   names(old_tables_harmonized))]]
  val_df  <- dfs_new[[paste0("val_", tp, ".cxcr7.vs.0s")]]
  
  i_sig <- get_sig_adjp(init_df, fc_cutoff = 0.5)
  v_sig <- get_sig_adjp(val_df, fc_cutoff = 0.5)
  
  conc_up   <- length(intersect(i_sig$up, v_sig$up))
  conc_down <- length(intersect(i_sig$down, v_sig$down))
  disc_1    <- length(intersect(i_sig$up, v_sig$down))
  disc_2    <- length(intersect(i_sig$down, v_sig$up))
  total     <- conc_up + conc_down + disc_1 + disc_2
  conc_pct  <- if (total > 0) 100 * (conc_up + conc_down) / total else 0
  
  cat(sprintf("%5ss    | %9d | %11d | %11d | %9d | %5.1f%%\n",
              tp, conc_up, conc_down, conc_up + conc_down,
              disc_1 + disc_2, conc_pct))
}


###############################################################
# 3) LIST SHARED SITES with gene names (both cutoffs)
###############################################################

cat("\n\n========================================================\n")
cat("  Shared phosphosites with gene names\n")
cat("  adj.P.Val < 0.05 & |logFC| > 0.5\n")
cat("========================================================\n")

for (tp in timepoints) {
  init_df <- old_tables_raw[[match(paste0("init_", tp, ".cxcr7.vs.0s"),
                                   names(old_tables_harmonized))]]
  val_df  <- dfs_new[[paste0("val_", tp, ".cxcr7.vs.0s")]]
  
  i_sig <- get_sig_adjp(init_df, fc_cutoff = 0.5)
  v_sig <- get_sig_adjp(val_df, fc_cutoff = 0.5)
  
  conc_up   <- sort(intersect(i_sig$up, v_sig$up))
  conc_down <- sort(intersect(i_sig$down, v_sig$down))
  disc      <- c(intersect(i_sig$up, v_sig$down), intersect(i_sig$down, v_sig$up))
  
  cat(sprintf("\n=== %ss (adj.P.Val<0.05 & |FC|>0.5) ===\n", tp))
  cat(sprintf("  Concordant UP: %d | DOWN: %d | Discordant: %d\n",
              length(conc_up), length(conc_down), length(disc)))
  
  if (length(conc_up) > 0)
    cat("  UP:  ", paste(id_to_gene(conc_up), collapse = ", "), "\n")
  if (length(conc_down) > 0)
    cat("  DOWN:", paste(id_to_gene(conc_down), collapse = ", "), "\n")
  if (length(disc) > 0)
    cat("  DISC:", paste(id_to_gene(disc), collapse = ", "), "\n")
}


