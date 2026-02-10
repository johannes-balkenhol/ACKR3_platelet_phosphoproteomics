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


library(umap)
dev.new()
set.seed(1)
umap_res <- umap(t(logfc_matrix))   # transpose = samples are columns

plot(umap_res$layout,
     col = "steelblue", pch = 19, cex = 2,
     main = "UMAP of 16 phosphoproteomics profiles")

text(umap_res$layout, labels = colnames(logfc_matrix),
     pos = 3, cex = 0.8)


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
library(umap)


dev.new()
umap_res_top <- umap(t(logfc_top))

plot(
  umap_res_top$layout[,1],
  umap_res_top$layout[,2],
  pch = 19, cex = 1.5,
  col = "darkred",
  main = paste("UMAP of 16 contrasts (top", top_n, "sites)"),
  xlab = "UMAP-1", ylab = "UMAP-2"
)

text(
  umap_res_top$layout[,1],
  umap_res_top$layout[,2],
  labels = colnames(logfc_top),
  pos = 3, cex = 0.7
)


###############################################################
## Pathway enrichment











