### ----------------------------------------------------------
### QC intensities
### ----------------------------------------------------------

library(tidyverse)
library(umap)
library(ggplot2)
library(patchwork)

###############################################################
# 1) LOAD RAW + NORMALIZED INTENSITIES FOR INIT + VAL
###############################################################


val_raw  <- read.delim("D:/Research/CXCR7_platelet_analysis/SubProjects/CXCR7_validation/phosphoproteomics/data/raw_data/raw_intensities.txt",
                       check.names = FALSE)

val_norm <- read.delim("D:/Research/CXCR7_platelet_analysis/SubProjects/CXCR7_validation/phosphoproteomics/data/processed_data/norm_intensity.txt",
                       check.names = FALSE)

init_raw <- read.delim("D:/Research/CXCR7_platelet_analysis/SubProjects/CXCR7_initial/phosphoproteomics/data/raw_data/raw_intensities.txt",
                       check.names = FALSE)

init_norm <- read.delim("D:/Research/CXCR7_platelet_analysis/SubProjects/CXCR7_initial/phosphoproteomics/data/processed_data/norm_intensity.txt",
                        check.names = FALSE)


###############################################################
# 2) INSPECT STRUCTURE
###############################################################

cat("\n--- VALIDATION RAW ---\n")
print(dim(val_raw))
print(head(colnames(val_raw)))

cat("\n--- INITIAL RAW ---\n")
print(dim(init_raw))
print(head(colnames(init_raw)))

cat("\n--- VALIDATION NORMALIZED ---\n")
print(dim(val_norm))

cat("\n--- INITIAL NORMALIZED ---\n")
print(dim(init_norm))

cat("\n===== VALIDATION RAW =====\n")
cat("dim: "); print(dim(val_raw))
cat("colnames:\n"); print(head(colnames(val_raw), 20))
cat("rownames:\n"); print(head(rownames(val_raw), 20))

cat("\n===== INITIAL RAW =====\n")
cat("dim: "); print(dim(init_raw))
cat("colnames:\n"); print(head(colnames(init_raw), 20))
cat("rownames:\n"); print(head(rownames(init_raw), 20))

cat("\n===== VALIDATION NORMALIZED =====\n")
cat("dim: "); print(dim(val_norm))
cat("colnames:\n"); print(head(colnames(val_norm), 20))
cat("rownames:\n"); print(head(rownames(val_norm), 20))

cat("\n===== INITIAL NORMALIZED =====\n")
cat("dim: "); print(dim(init_norm))
cat("colnames:\n"); print(head(colnames(init_norm), 20))
cat("rownames:\n"); print(head(rownames(init_norm), 20))


###############################################################
# 2) harmonize
###############################################################

extract_uid_gene_psite <- function(df) {
  
  raw_ids <- rownames(df)
  
  clean_ids <- sapply(raw_ids, function(x) {
    parts <- strsplit(x, ";", fixed = TRUE)[[1]]
    
    # Must contain UniProt ; Gene ; Site
    if (length(parts) < 3) return(NA)
    
    uni  <- parts[1]
    gene <- parts[2]
    site <- parts[3]              # S123, T55, Y999, S17|T21, …
    
    # Keep ONLY if contains S/T/Y
    if (!grepl("[STY]", site)) return(NA)
    
    paste(uni, gene, site, sep = ";")
  })
  
  df$clean_id <- clean_ids
  df
}


collapse_by_clean_id <- function(df) {
  
  df2 <- df %>%
    dplyr::filter(!is.na(clean_id)) %>%   # remove invalid rows
    dplyr::group_by(clean_id) %>%
    dplyr::summarise(
      dplyr::across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    as.data.frame()
  
  rownames(df2) <- df2$clean_id
  df2$clean_id <- NULL
  
  df2
}


###############################################################
# HARMONIZE ALL
###############################################################

# 1) extract IDs
val_raw_step1  <- extract_uid_gene_psite(val_raw)
init_raw_step1 <- extract_uid_gene_psite(init_raw)
val_norm_step1 <- extract_uid_gene_psite(val_norm)
init_norm_step1<- extract_uid_gene_psite(init_norm)

# 2) collapse
val_raw_clean  <- collapse_by_clean_id(val_raw_step1)
init_raw_clean <- collapse_by_clean_id(init_raw_step1)
val_norm_clean <- collapse_by_clean_id(val_norm_step1)
init_norm_clean<- collapse_by_clean_id(init_norm_step1)


cat("\nVAL RAW CLEAN:");    print(dim(val_raw_clean)); 
sum(duplicated(rownames(val_raw_clean)))

cat("\nINIT RAW CLEAN:");   print(dim(init_raw_clean));
sum(duplicated(rownames(init_raw_clean)))

cat("\nVAL NORM CLEAN:");   print(dim(val_norm_clean));
sum(duplicated(rownames(val_norm_clean)))

cat("\nINIT NORM CLEAN:");  print(dim(init_norm_clean));
sum(duplicated(rownames(init_norm_clean)))



dim(val_raw_clean)
dim(init_raw_clean)
dim(val_norm_clean)
dim(init_norm_clean)



clean_colnames_unified <- function(df) {
  cn <- colnames(df)
  
  new <- sapply(cn, function(x) {
    
    #### INIT format: X0030_7 ####
    if (grepl("^X\\d{4}_\\d+$", x)) {
      time_raw <- sub("^X(\\d{4}).*", "\\1", x)
      time_sec <- as.numeric(time_raw)      # 0010 → 10, 0600 → 600
      rep      <- sub(".*_(\\d+)$", "\\1", x)
      rep      <- sprintf("%02d", as.numeric(rep))
      
      return(paste0("CXCR7.", time_sec, "s_", rep))
    }
    
    
    #### VALIDATION format: x10sek_CXCR7_03 ####
    if (grepl("^x\\d+sek_", x)) {
      time_sec <- sub("^x(\\d+)sek_.*", "\\1", x)
      
      treat <- sub("^x\\d+sek_([^_]+).*", "\\1", x)
      treat <- toupper(treat)
      treat[treat %in% c("CTRL","CONTROL")] <- "DMSO"
      
      rep <- sub(".*_(\\d+)$", "\\1", x)
      rep <- sprintf("%02d", as.numeric(rep))
      
      return(paste0(treat, ".", time_sec, "s_", rep))
    }
    
    #### fallback, keep original ####
    return(x)
  })
  
  colnames(df) <- new
  df
}


val_raw_clean  <- clean_colnames_unified(val_raw_clean)
val_norm_clean <- clean_colnames_unified(val_norm_clean)
init_raw_clean <- clean_colnames_unified(init_raw_clean)
init_norm_clean<- clean_colnames_unified(init_norm_clean)


head(colnames(val_raw_clean), 20)
head(colnames(init_norm_clean), 20)


############################################
## Reproducibility
###########################################
###############################################################################
# 0) PACKAGES
###############################################################################
library(dplyr)
library(pheatmap)

###############################################################################
# 1) CORRELATION FUNCTION (SAFE VERSION)
###############################################################################

compute_cor_matrix_safe <- function(mat) {
  
  # remove rows with all NA
  mat <- mat[rowSums(is.na(mat)) < ncol(mat), , drop = FALSE]
  
  # replace columns with zero variance / full NA with tiny noise
  bad_cols <- sapply(as.data.frame(mat), function(x) sd(x, na.rm=TRUE) == 0 | all(is.na(x)))
  if (any(bad_cols)) {
    message("Replacing ", sum(bad_cols), " zero-variance columns with small noise.")
    mat[, bad_cols] <- matrix(rnorm(sum(bad_cols) * nrow(mat), sd=1e-6),
                              nrow = nrow(mat))
  }
  
  # compute correlation
  cor(t(mat), method = "pearson", use = "pairwise.complete.obs")
}

###############################################################################
# 2) RUN CORRELATION MATRICES
###############################################################################

cor_val_raw  <- compute_cor_matrix_safe(val_raw_clean)
cor_init_raw <- compute_cor_matrix_safe(init_raw_clean)

cor_val_norm  <- compute_cor_matrix_safe(val_norm_clean)
cor_init_norm <- compute_cor_matrix_safe(init_norm_clean)




###############################################
# CROSS-DATASET CORRELATION: RAW
###############################################

# intersect phosphosites
common_ids_raw <- intersect(rownames(val_raw_clean), rownames(init_raw_clean))

val_raw_aligned  <- val_raw_clean[common_ids_raw, ]
init_raw_aligned <- init_raw_clean[common_ids_raw, ]

cross_raw_cor <- cor(
  val_raw_aligned,
  init_raw_aligned,
  use = "pairwise.complete.obs",
  method = "pearson"
)

###############################################
# CROSS-DATASET CORRELATION: NORMALIZED
###############################################

common_ids_norm <- intersect(rownames(val_norm_clean), rownames(init_norm_clean))

val_norm_aligned  <- val_norm_clean[common_ids_norm, ]
init_norm_aligned <- init_norm_clean[common_ids_norm, ]

cross_norm_cor <- cor(
  val_norm_aligned,
  init_norm_aligned,
  use = "pairwise.complete.obs",
  method = "pearson"
)

dev.new()
pheatmap(
  cross_raw_cor,
  main = "Initial vs Validation — RAW intensities",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  color = colorRampPalette(c("blue","white","red"))(200)
)

dev.new()
pheatmap(
  cross_norm_cor,
  main = "Initial vs Validation — RUV Normalized intensities",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  color = colorRampPalette(c("blue","white","red"))(200)
)



###############################################################################
# umap

###############################################################
# 0) FAST correlation helper
###############################################################
compute_cor <- function(mat) {
  cor(t(mat), method = "pearson", use = "pairwise.complete.obs")
}

###############################################################
# 1) Compute correlations WITHOUT plotting huge heatmaps
###############################################################

cor_val_raw  <- compute_cor(val_raw_clean)
cor_init_raw <- compute_cor(init_raw_clean)

cor_val_norm  <- compute_cor(val_norm_clean)
cor_init_norm <- compute_cor(init_norm_clean)

# summary metrics:
cor_stats <- function(cmat) {
  vals <- cmat[upper.tri(cmat)]
  c(
    mean = mean(vals, na.rm=TRUE),
    median = median(vals, na.rm=TRUE),
    sd = sd(vals, na.rm=TRUE),
    q25 = quantile(vals, 0.25, na.rm=TRUE),
    q75 = quantile(vals, 0.75, na.rm=TRUE)
  )
}

cat("\n=== REPRODUCIBILITY SUMMARY ===\n")
print(list(
  VAL_RAW  = cor_stats(cor_val_raw),
  INIT_RAW = cor_stats(cor_init_raw),
  VAL_NORM = cor_stats(cor_val_norm),
  INIT_NORM = cor_stats(cor_init_norm)
))

###############################################################
# 2) INTERSECT proteins & combine VAL + INIT for UMAP
###############################################################

# RAW
common_raw <- intersect(rownames(val_raw_clean), rownames(init_raw_clean))
raw_combined <- cbind(
  VAL  = val_raw_clean[common_raw, ],
  INIT = init_raw_clean[common_raw, ]
)

# NORM
common_norm <- intersect(rownames(val_norm_clean), rownames(init_norm_clean))
norm_combined <- cbind(
  VAL  = val_norm_clean[common_norm, ],
  INIT = init_norm_clean[common_norm, ]
)

###############################################################
# 3) UMAP
###############################################################
library(uwot)
set.seed(123)

# prepare
umap_raw_in  <- t(raw_combined)
umap_norm_in <- t(norm_combined)

# run
umap_raw  <- umap(umap_raw_in,  n_neighbors = 15, min_dist = 0.3)
umap_norm <- umap(umap_norm_in, n_neighbors = 15, min_dist = 0.3)

# convert to df
df_raw <- data.frame(
  UMAP1 = umap_raw[,1],
  UMAP2 = umap_raw[,2],
  sample = rownames(umap_raw_in),
  dataset = ifelse(grepl("^VAL", rownames(umap_raw_in)), "VAL", "INIT")
)

df_norm <- data.frame(
  UMAP1 = umap_norm[,1],
  UMAP2 = umap_norm[,2],
  sample = rownames(umap_norm_in),
  dataset = ifelse(grepl("^VAL", rownames(umap_norm_in)), "VAL", "INIT")
)

###############################################################
# 4) PLOT
###############################################################
library(ggplot2)

p_raw <- ggplot(df_raw, aes(UMAP1, UMAP2, color = dataset)) +
  geom_point(size=3) +
  theme_minimal(base_size=14) +
  ggtitle("UMAP – RAW")

p_norm <- ggplot(df_norm, aes(UMAP1, UMAP2, color = dataset)) +
  geom_point(size=3) +
  theme_minimal(base_size=14) +
  ggtitle("UMAP – NORM")

p_raw
p_norm





###############################################################
# 1) intersect shared samples between VAL and INIT
###############################################################

shared_samples <- intersect(
  colnames(val_norm_clean),
  colnames(init_norm_clean)
)

length(shared_samples)
print(shared_samples)


###############################################################
# 2) intersect shared phosphosites
###############################################################

shared_psites <- intersect(
  rownames(val_norm_clean),
  rownames(init_norm_clean)
)



val_sub  <- val_norm_clean[shared_psites, shared_samples, drop=FALSE]
init_sub <- init_norm_clean[shared_psites, shared_samples, drop=FALSE]

###############################################################
# 3) compute replicate correlations
###############################################################

cor_val  <- cor(val_sub,  use="pairwise.complete.obs", method="pearson")
cor_init <- cor(init_sub, use="pairwise.complete.obs", method="pearson")

library(pheatmap)

pheatmap(cor_val,  main="VAL reproducibility", fontsize=8)
pheatmap(cor_init, main="INIT reproducibility", fontsize=8)


###############################################################
# 4) UMAP for RAW and NORM
###############################################################

# RAW
common_raw <- intersect(rownames(val_raw_clean), rownames(init_raw_clean))
raw_mat <- cbind(val_raw_clean[common_raw, ], init_raw_clean[common_raw, ])

# NORM
common_norm <- intersect(rownames(val_norm_clean), rownames(init_norm_clean))
norm_mat <- cbind(val_norm_clean[common_norm, ], init_norm_clean[common_norm, ])







###############################################################
# UMAP on VAL + INIT (normalized data)
###############################################################
###############################################################
# UMAP on VAL + INIT (normalized data)
###############################################################

library(uwot)
library(ggplot2)
library(dplyr)
library(stringr)

# ---------- 1) restrict to shared phosphosites ----------
shared_ids <- intersect(rownames(val_norm_clean), rownames(init_norm_clean))

val_use  <- val_norm_clean[shared_ids, ]
init_use <- init_norm_clean[shared_ids, ]

# ---------- 2) combine ----------
combined_norm <- cbind(val_use, init_use)

# ---------- 3) feature filter (remove all-NA or zero-variance) ----------
make_umap_safe <- function(mat, min_feat = 20) {
  mat <- mat[rowSums(is.na(mat)) == 0, ]
  mat <- mat[apply(mat, 1, sd) > 0, ]
  if(nrow(mat) < min_feat) stop("Not enough rows for UMAP.")
  mat
}

umap_input <- make_umap_safe(combined_norm)

# ---------- 4) UMAP ----------
set.seed(123)
umap_res <- umap(
  t(umap_input),
  n_neighbors = 15,
  min_dist = 0.3,
  n_components = 2,
  scale = FALSE,
  pca = 50
)

# ---------- 5) build metadata ----------
parse_sample <- function(x) {
  # e.g. "CXCR7.600s_05"
  parts <- str_split(x, "\\.", simplify = TRUE)
  cond  <- parts[1]
  rest  <- parts[2]
  
  tparts <- str_split(rest, "_", simplify = TRUE)
  time  <- tparts[1]
  rep   <- tparts[2]
  
  data.frame(sample = x, condition = cond, time = time, rep = rep)
}

meta <- bind_rows(lapply(colnames(combined_norm), parse_sample))

# ---------- 6) dataframe for ggplot ----------
umap_df <- data.frame(
  UMAP1 = umap_res[,1],
  UMAP2 = umap_res[,2],
  sample = colnames(combined_norm)
) %>% left_join(meta, by = "sample")

# ---------- 7) Plot ----------
dev.new()
ggplot(umap_df, aes(UMAP1, UMAP2, label = sample,
                    color = condition, shape = time)) +
  geom_point(size = 3) +
  geom_text(size = 2, vjust = -0.5) +
  theme_bw() +
  ggtitle("UMAP – VAL + INIT (normalized)") +
  theme(plot.title = element_text(hjust = 0.5))




###############################################################
# UMAP on VAL + INIT (raw data) — NA-safe version
###############################################################

library(uwot)
library(ggplot2)
library(dplyr)
library(stringr)

# ---------- 1) restrict to shared phosphosites ----------
shared_ids_raw <- intersect(rownames(val_raw_clean), rownames(init_raw_clean))

valR_use  <- val_raw_clean[shared_ids_raw, ]
initR_use <- init_raw_clean[shared_ids_raw, ]

# ---------- 2) combine ----------
combined_raw <- cbind(valR_use, initR_use)

# ---------- 3) feature filter (remove all-NA or zero-variance) ----------
make_umap_safe_raw <- function(mat, min_feat = 20) {
  
  # remove any row that contains *any* NA (UMAP requirement)
  mat <- mat[rowSums(is.na(mat)) == 0, ]
  
  # remove zero-variance rows
  mat <- mat[apply(mat, 1, sd, na.rm = TRUE) > 0, ]
  
  if(nrow(mat) < min_feat)
    stop("Not enough rows for UMAP (raw).")
  
  mat
}

raw_umap_input <- make_umap_safe_raw(combined_raw)

# ---------- 4) log2 + scale ----------
raw_umap_input <- log2(raw_umap_input + 1)
raw_umap_input <- t(scale(t(raw_umap_input)))

# ---------- 5) UMAP ----------
set.seed(123)
raw_umap_res <- umap(
  t(raw_umap_input),
  n_neighbors = 15,
  min_dist = 0.3,
  n_components = 2,
  scale = FALSE,
  pca = 50
)

# ---------- 6) metadata parser ----------
parse_sample <- function(x) {
  parts <- str_split(x, "\\.", simplify = TRUE)
  cond  <- parts[1]
  rest  <- parts[2]
  
  tparts <- str_split(rest, "_", simplify = TRUE)
  time  <- tparts[1]
  rep   <- tparts[2]
  
  data.frame(sample = x, condition = cond, time = time, rep = rep)
}

meta_raw <- bind_rows(lapply(colnames(combined_raw), parse_sample))

# ---------- 7) UMAP dataframe ----------
raw_umap_df <- data.frame(
  UMAP1 = raw_umap_res[,1],
  UMAP2 = raw_umap_res[,2],
  sample = colnames(combined_raw)
) %>%
  left_join(meta_raw, by = "sample")

dev.new()
# ---------- 8) Plot ----------
ggplot(raw_umap_df, aes(UMAP1, UMAP2, label = sample,
                        color = condition, shape = time)) +
  geom_point(size = 3) +
  geom_text(size = 2, vjust = -0.5) +
  theme_bw() +
  ggtitle("UMAP – RAW VAL + INIT (no NA)") +
  theme(plot.title = element_text(hjust = 0.5))





###############################################################
# REPRODUCIBILITY ANALYSIS – RAW vs NORM, VAL vs INIT
###############################################################

library(dplyr)
library(ggplot2)
library(stringr)


###############################################################
# 1) Function: compute within-condition/time replicate correlations
###############################################################
# sample name example: "CXCR7.600s_05"
# group key:          "CXCR7.600s"

replicate_cor <- function(mat) {
  
  samples <- colnames(mat)
  
  # remove replicate number → “CXCR7.600s”
  group_id <- sub("_[0-9]+$", "", samples)
  
  cors <- c()
  
  for (g in unique(group_id)) {
    
    grp_samples <- samples[group_id == g]
    
    if (length(grp_samples) > 1) {
      
      submat <- mat[, grp_samples, drop = FALSE]
      cmat   <- suppressWarnings(cor(submat, use = "pairwise.complete.obs"))
      
      cors <- c(cors, cmat[upper.tri(cmat)])
    }
  }
  
  return(cors)
}


###############################################################
# 2) Compute reproducibility inside each dataset
###############################################################

cor_raw_val   <- replicate_cor(val_raw_clean)
cor_raw_init  <- replicate_cor(init_raw_clean)
cor_norm_val  <- replicate_cor(val_norm_clean)
cor_norm_init <- replicate_cor(init_norm_clean)


###############################################################
# 3) VAL–INIT reproducibility (pairing identical replicates)
###############################################################

# shared sample names
shared_cols <- intersect(colnames(val_norm_clean), colnames(init_norm_clean))
cat("Matched VAL–INIT replicates:", length(shared_cols), "\n")

# intersect phosphosite IDs
shared_ids <- intersect(rownames(val_norm_clean), rownames(init_norm_clean))
cat("Shared phosphosites between VAL and INIT:", length(shared_ids), "\n")

# subset rows + columns
val_s  <- val_norm_clean[shared_ids, shared_cols, drop = FALSE]
init_s <- init_norm_clean[shared_ids, shared_cols, drop = FALSE]

# sanity check
stopifnot(
  all(colnames(val_s) == colnames(init_s)),
  all(rownames(val_s) == rownames(init_s))
)

# compute one correlation per sample
pairwise_cor <- sapply(shared_cols, function(s) {
  cor(val_s[, s], init_s[, s], use = "pairwise.complete.obs")
})

pairwise_cor



###############################################################
# 4) Plot replicate reproducibility across datasets
###############################################################

df <- data.frame(
  cor = c(cor_raw_val, cor_raw_init, cor_norm_val, cor_norm_init),
  dataset = factor(c(
    rep("RAW_VAL",  length(cor_raw_val)),
    rep("RAW_INIT", length(cor_raw_init)),
    rep("NORM_VAL", length(cor_norm_val)),
    rep("NORM_INIT", length(cor_norm_init))
  ), levels = c("RAW_VAL", "RAW_INIT", "NORM_VAL", "NORM_INIT"))
)

ggplot(df, aes(dataset, cor, fill = dataset)) +
  geom_boxplot(outlier.alpha = 0.2) +
  theme_bw(base_size = 14) +
  ggtitle("Reproducibility – RAW vs NORMALIZED (VAL + INIT)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 35, hjust = 1)
  ) +
  ylab("Replicate Pearson correlation") +
  xlab("Dataset")


###############################################################
# 5) VAL–INIT similarity density plot
###############################################################

df_valinit <- data.frame(cor = pairwise_cor)

ggplot(df_valinit, aes(x = cor)) +
  geom_density(fill = "steelblue", alpha = 0.35, linewidth = 0.8) +
  geom_vline(xintercept = median(df_valinit$cor), color = "red", lwd = 1) +
  theme_bw(base_size = 14) +
  ggtitle("VAL vs INIT reproducibility (normalized – paired replicates)") +
  xlab("Pearson correlation") +
  theme(plot.title = element_text(hjust = 0.5))

###############################################################
###############################################################
###############################################################
# Investigate tie traces
###############################################################


###############################################################
# KINETIC REPRODUCIBILITY – shared timepoints VAL vs INIT
###############################################################

library(dplyr)
library(stringr)
library(dtw)
library(ggplot2)

###############################################################
# 0) Restrict to shared phosphosites
###############################################################

shared_ids <- intersect(rownames(val_norm_clean), rownames(init_norm_clean))

val_norm_filt  <- val_norm_clean[shared_ids, ]
init_norm_filt <- init_norm_clean[shared_ids, ]

cat("Shared phosphosites:", length(shared_ids), "\n")


###############################################################
# 1) Extract all timepoints
###############################################################

get_time <- function(x) sub(".*\\.(.*)_.*", "\\1", x)

val_times  <- unique(get_time(colnames(val_norm_filt)))
init_times <- unique(get_time(colnames(init_norm_filt)))

cat("VAL times:",  val_times,  "\n")
cat("INIT times:", init_times, "\n")

# Only times present in both experiments
shared_times <- intersect(val_times, init_times)
cat("Shared times:", shared_times, "\n")


###############################################################
# 2) Function: average replicates per timepoint
###############################################################

mean_by_time <- function(mat, times) {
  res <- sapply(times, function(t) {
    cols <- grep(paste0("\\.", t, "_"), colnames(mat), value = TRUE)
    if (length(cols) == 0)
      stop(paste("No columns found for time", t))
    rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE)
  })
  colnames(res) <- times
  return(res)
}


###############################################################
# 3) Build time-series matrices (VAL & INIT)
###############################################################

val_ts  <- mean_by_time(val_norm_filt,  shared_times)
init_ts <- mean_by_time(init_norm_filt, shared_times)

cat("VAL TS dim: ",  dim(val_ts),  "\n")
cat("INIT TS dim:",  dim(init_ts), "\n")

# safety check
stopifnot(identical(rownames(val_ts), rownames(init_ts)))
stopifnot(identical(colnames(val_ts), colnames(init_ts)))


###############################################################
# 4) Pearson correlation kinetics per phosphosite
###############################################################

pearson_kin <- sapply(rownames(val_ts), function(pid) {
  cor(val_ts[pid, ], init_ts[pid, ], use = "pairwise.complete.obs")
})


###############################################################
# 5) DTW distance per phosphosite
###############################################################

dtw_dist <- sapply(rownames(val_ts), function(pid) {
  v <- val_ts[pid, ]
  i <- init_ts[pid, ]
  dtw(v, i, keep = FALSE)$distance
})


###############################################################
# 6) Histograms
###############################################################
dev.new()
par(mfrow = c(1, 2))
hist(pearson_kin, breaks = 40, main = "Pearson kinetics (VAL vs INIT)",
     xlab = "Pearson correlation")
hist(dtw_dist, breaks = 40, main = "DTW distance (VAL vs INIT)",
     xlab = "DTW distance")


###############################################################
# 7) Scatter plot: DTW vs Pearson
###############################################################
dev.new()
df_compare <- data.frame(
  pearson = pearson_kin,
  dtw = dtw_dist
)

ggplot(df_compare, aes(dtw, pearson)) +
  geom_point(alpha = 0.3) +
  theme_bw() +
  ggtitle("Kinetic reproducibility: VAL vs INIT") +
  xlab("DTW distance") +
  ylab("Pearson correlation")




###############################################################
###############################################################
###############################################################

###############################################################
# OPTION 1 — Interpolate VAL → INIT time structure
# Produces VAL and INIT time-series matrices with identical grids
###############################################################

library(dplyr)
library(stringr)
library(splines2)

###############################################################
# 1) Helper: extract "10s" from "CXCR7.10s_03"
###############################################################

get_time <- function(x) {
  sub(".*\\.(.*)_.*", "\\1", x)
}

###############################################################
# 2) Extract available timepoints (character times → numeric seconds)
###############################################################

val_times_raw  <- get_time(colnames(val_norm_clean))
init_times_raw <- get_time(colnames(init_norm_clean))

val_times_unique  <- sort(unique(val_times_raw))
init_times_unique <- sort(unique(init_times_raw))

cat("VAL times: ",  paste(val_times_unique, collapse=", "), "\n")
cat("INIT times:", paste(init_times_unique, collapse=", "), "\n")

# convert "600s" → 600 numeric
to_sec <- function(x) as.numeric(sub("s$", "", x))

val_t_sec  <- to_sec(val_times_unique)
init_t_sec <- to_sec(init_times_unique)

###############################################################
# 3) Function: average replicates at each timepoint
###############################################################

mean_by_time <- function(mat, times_char) {
  out <- sapply(times_char, function(t) {
    cols <- grep(paste0("\\.", t, "_"), colnames(mat), value = TRUE)
    if (length(cols) == 0)
      stop(paste("No columns for", t))
    rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE)
  })
  colnames(out) <- times_char
  return(out)
}

###############################################################
# 4) Build the raw kinetic matrices (before interpolation)
###############################################################

val_kin_raw  <- mean_by_time(val_norm_clean,  val_times_unique)
init_kin_raw <- mean_by_time(init_norm_clean, init_times_unique)

# convert column names to times in seconds
colnames(val_kin_raw)  <- val_t_sec
colnames(init_kin_raw) <- init_t_sec

###############################################################
# 5) Interpolate VAL onto INIT time grid (spline)
###############################################################

target_grid <- init_t_sec
cat("Interpolation target grid (INIT):", target_grid, "\n")

val_int <- matrix(NA, nrow=nrow(val_kin_raw), ncol=length(target_grid))
rownames(val_int) <- rownames(val_kin_raw)
colnames(val_int) <- target_grid

init_int <- init_kin_raw  # already on the grid

for (pid in rownames(val_kin_raw)) {
  
  x_val <- as.numeric(colnames(val_kin_raw))
  y_val <- as.numeric(val_kin_raw[pid, ])
  
  # spline interpolation
  val_int[pid, ] <- spline(x = x_val, y = y_val, xout = target_grid)$y
}

###############################################################
# 6) Output summary
###############################################################

cat("\nInterpolated VAL dim:  ", dim(val_int), "\n")
cat("INIT kinetic dim:      ", dim(init_int), "\n")

cat("\nExample site check:\n")
print(rbind(
  VAL_before = val_kin_raw[1, ],
  VAL_after  = val_int[1, ],
  INIT       = init_int[1, ]
))



###############################################################
###############################################################
###############################################################

###############################################################
# OPTION 1 — Interpolate VAL → INIT time structure
# With cross-experiment z-score normalization (recommended)
###############################################################

library(dplyr)
library(stringr)
library(splines2)

###############################################################
# 0) Z-score normalize each experiment independently
# (fixes scale mismatch created by separate RUV3 runs)
###############################################################

val_z  <- t(scale(t(val_norm_clean)))
init_z <- t(scale(t(init_norm_clean)))

###############################################################
# 1) Helper: extract "10s" from "CXCR7.10s_03"
###############################################################

get_time <- function(x) {
  sub(".*\\.(.*)_.*", "\\1", x)
}

###############################################################
# 2) Extract timepoints
###############################################################

val_times_raw  <- get_time(colnames(val_z))
init_times_raw <- get_time(colnames(init_z))

val_times_unique  <- sort(unique(val_times_raw))
init_times_unique <- sort(unique(init_times_raw))

cat("VAL times: ",  paste(val_times_unique, collapse=", "), "\n")
cat("INIT times:", paste(init_times_unique, collapse=", "), "\n")

to_sec <- function(x) as.numeric(sub("s$", "", x))

val_t_sec  <- to_sec(val_times_unique)
init_t_sec <- to_sec(init_times_unique)

###############################################################
# 3) Function: average replicates at each timepoint
###############################################################

mean_by_time <- function(mat, times_char) {
  out <- sapply(times_char, function(t) {
    cols <- grep(paste0("\\.", t, "_"), colnames(mat), value = TRUE)
    if (length(cols) == 0)
      stop(paste("No columns for", t))
    rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE)
  })
  colnames(out) <- times_char
  return(out)
}

###############################################################
# 4) Build raw kinetic matrices after z-scoring
###############################################################

val_kin_raw  <- mean_by_time(val_z,  val_times_unique)
init_kin_raw <- mean_by_time(init_z, init_times_unique)

colnames(val_kin_raw)  <- val_t_sec
colnames(init_kin_raw) <- init_t_sec

###############################################################
# 5) Interpolate VAL → INIT time grid
###############################################################

target_grid <- init_t_sec
cat("Interpolation target grid (INIT):", target_grid, "\n")

val_int <- matrix(NA, nrow=nrow(val_kin_raw), ncol=length(target_grid))
rownames(val_int) <- rownames(val_kin_raw)
colnames(val_int) <- target_grid

init_int <- init_kin_raw  # already on INIT grid

for (pid in rownames(val_kin_raw)) {
  
  x_val <- as.numeric(colnames(val_kin_raw))
  y_val <- as.numeric(val_kin_raw[pid, ])
  
  val_int[pid, ] <- spline(x = x_val, y = y_val, xout = target_grid)$y
}

###############################################################
# 6) Output summary
###############################################################

cat("\nInterpolated VAL dim:  ", dim(val_int), "\n")
cat("INIT kinetic dim:      ", dim(init_int), "\n")

cat("\nExample site check:\n")
print(rbind(
  VAL_before = val_kin_raw[1, ],
  VAL_after  = val_int[1, ],
  INIT       = init_int[1, ]
))


###############################################################
# OPTION A — KINETIC REPRODUCIBILITY (Pearson + DTW)
###############################################################

library(dtw)
library(ggplot2)

# val_int : interpolated VAL kinetics (rows = sites, cols = timepoints)
# init_int: INIT kinetics on same grid

# Ensure same phosphosites (rows)
shared_pids <- intersect(rownames(val_int), rownames(init_int))
val_k  <- val_int[shared_pids, , drop=FALSE]
init_k <- init_int[shared_pids, , drop=FALSE]

cat("Shared phosphosites =", length(shared_pids), "\n")

###############################################################
# 1) Pearson correlation for each phosphosite
###############################################################

pearson_kin <- sapply(shared_pids, function(pid) {
  cor(val_k[pid, ], init_k[pid, ], use = "pairwise.complete.obs")
})

###############################################################
# 2) DTW distance for each phosphosite
###############################################################

dtw_dist <- sapply(shared_pids, function(pid) {
  v <- val_k[pid, ]
  i <- init_k[pid, ]
  dtw(v, i, keep = FALSE)$distance
})

###############################################################
# 3) Composite reproducibility score (optional)
# Higher score = more reproducible
###############################################################

pearson_scaled <- (pearson_kin - min(pearson_kin)) / (max(pearson_kin) - min(pearson_kin))
dtw_scaled     <- 1 - (dtw_dist - min(dtw_dist)) / (max(dtw_dist) - min(dtw_dist))

repro_score <- 0.5 * pearson_scaled + 0.5 * dtw_scaled

###############################################################
# 4) Summary plots
###############################################################

par(mfrow=c(1,2))
hist(pearson_kin, breaks=40, col="steelblue", main="Pearson kinetics", xlab="cor")
hist(dtw_dist, breaks=40, col="tomato", main="DTW distance", xlab="distance")

###############################################################
# 5) Scatter Pearson vs DTW
###############################################################

df_compare <- data.frame(
  PID = shared_pids,
  pearson = pearson_kin,
  dtw = dtw_dist,
  score = repro_score
)

ggplot(df_compare, aes(dtw, pearson)) +
  geom_point(alpha=0.4) +
  theme_bw() +
  labs(title="Kinetic reproducibility (VAL vs INIT)",
       x="DTW distance",
       y="Pearson correlation")

###############################################################
# 6) Rank reproducible phosphosites
###############################################################

df_ranked <- df_compare[order(df_compare$score, decreasing = TRUE), ]
head(df_ranked, 20)


###############################################################
# UMAP of kinetic features (VAL + INIT) — fixed rownames
###############################################################

library(uwot)
library(ggplot2)
library(dplyr)

# val_k and init_k are the kinetic matrices (rows = phosphosites, cols = timepoints)

###############################################################
# 1) Z-score each kinetic vector (per phosphosite)
###############################################################

z_val  <- t(scale(t(val_k)))     # VAL kinetics normalized
z_init <- t(scale(t(init_k)))    # INIT kinetics normalized

# Make rownames unique
rownames(z_val)  <- paste0(rownames(z_val), "_VAL")
rownames(z_init) <- paste0(rownames(z_init), "_INIT")

###############################################################
# 2) Combine matrices
###############################################################

all_mat <- rbind(z_val, z_init)

meta <- data.frame(
  PID = c(rownames(val_k), rownames(init_k)),
  set = c(rep("VAL",  nrow(val_k)),
          rep("INIT", nrow(init_k))),
  row.names = rownames(all_mat)
)

###############################################################
# 3) Run UMAP
###############################################################

set.seed(123)
umap_res <- umap(
  all_mat,
  n_neighbors = 20,
  min_dist = 0.2,
  metric = "euclidean",
  scale = FALSE,
  n_components = 2
)

umap_df <- data.frame(
  UMAP1 = umap_res[,1],
  UMAP2 = umap_res[,2],
  PID   = meta$PID,
  set   = meta$set
)

###############################################################
# 4) Plot UMAP
###############################################################

ggplot(umap_df, aes(UMAP1, UMAP2, color = set)) +
  geom_point(alpha = 0.4) +
  theme_bw() +
  ggtitle("UMAP of Kinetic Features (VAL vs INIT)") +
  theme(plot.title = element_text(hjust = 0.5))



# Z-score per phosphosite
val_z  <- t(scale(t(val_int)))
init_z <- t(scale(t(init_int)))

# Remove rows with NA (some sites are constant)
val_z  <- val_z[complete.cases(val_z), ]
init_z <- init_z[complete.cases(init_z), ]

dim(val_z)
dim(init_z)


# Clustering
val_dist <- dist(val_z, method="euclidean")
val_clust <- hclust(val_dist, method="ward.D2")

dev.new()
# Heatmap
heatmap(val_z[order(val_clust$order), ],
        Colv=NA, scale="none",
        main="VAL – Z-scored kinetics (clustered by phosphosite)",
        labRow=FALSE)



init_dist <- dist(init_z, method="euclidean")
init_clust <- hclust(init_dist, method="ward.D2")
dev.new()
heatmap(init_z[order(init_clust$order), ],
        Colv=NA, scale="none",
        main="INIT – Z-scored kinetics (clustered by phosphosite)",
        labRow=FALSE)



k <- 5
val_groups  <- cutree(val_clust,  k=k)
init_groups <- cutree(init_clust, k=k)

table(val_groups)
table(init_groups)


library(ggplot2)
library(reshape2)

plot_cluster_profiles <- function(zmat, groups, title) {
  df <- as.data.frame(zmat)
  df$cluster <- factor(groups)
  df_m <- melt(df, id.vars="cluster")
  colnames(df_m) <- c("cluster","time","value")
  
  ggplot(df_m, aes(x = time, y = value, group=cluster, color=cluster)) +
    stat_summary(fun="mean", geom="line", size=1.2) +
    theme_bw() + ggtitle(title)
}

dev.new()
plot_cluster_profiles(val_z,  val_groups,  "VAL – cluster mean kinetics")
dev.new()
plot_cluster_profiles(init_z, init_groups, "INIT – cluster mean kinetics")



# rank VAL
df_ranked_VAL <- df_compare_VAL[order(df_compare_VAL$score, decreasing=TRUE), ]

# rank INIT
df_ranked_INIT <- df_compare_INIT[order(df_compare_INIT$score, decreasing=TRUE), ]





library(dtw)
library(dplyr)

###############################################################
# 1) k-Means clustering on INIT (gold standard experiment)
###############################################################

set.seed(1)
k <- 5
km_init <- kmeans(init_int, centers = k, nstart = 50)
init_clusters <- km_init$cluster

###############################################################
# 2) Assign VAL sites to closest INIT cluster using DTW
###############################################################

assign_val_to_cluster <- function(val_vec) {
  dists <- sapply(1:k, function(cl) {
    init_centroid <- km_init$centers[cl, ]
    dtw(val_vec, init_centroid)$distance
  })
  which.min(dists)
}

val_clusters <- sapply(1:nrow(val_int), function(i) {
  assign_val_to_cluster(val_int[i, ])
})

names(val_clusters) <- rownames(val_int)

###############################################################
# 3) Identify reproducible phosphosites
###############################################################

shared_ids <- intersect(names(val_clusters), names(init_clusters))

df_cluster_compare <- data.frame(
  PID = shared_ids,
  init_cluster = init_clusters[shared_ids],
  val_cluster  = val_clusters[shared_ids],
  reproducible = init_clusters[shared_ids] == val_clusters[shared_ids]
)

table(df_cluster_compare$reproducible)

###############################################################
# 4) Extract reproducible modules
###############################################################

repro_sites <- df_cluster_compare %>% filter(reproducible)

cat("Reproducible phosphosites:", nrow(repro_sites), "\n")



