### ----------------------------------------------------------
### QC intensities
### ----------------------------------------------------------

library(tidyverse)
library(umap)
library(ggplot2)
library(patchwork)
library(ggplot2)
library(patchwork)
library(reshape2)
library(matrixStats)
library(scales)
library(RColorBrewer)
library(factoextra)

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


raw  <- init_raw

norm <- init_norm



###############################################################
# ALIGN ROWS BETWEEN RAW AND NORMALIZED MATRICES
###############################################################

common_rows <- intersect(rownames(raw), rownames(norm))

raw_aligned  <- raw[common_rows, , drop = FALSE]
norm_aligned <- norm[common_rows, , drop = FALSE]

cat("Number of common phosphosites:", length(common_rows), "\n")

raw <- raw_aligned
norm <- norm_aligned


##############################################################
# 2) CHECK BASIC CONSISTENCY
###############################################################

cat("Raw dimensions: "); print(dim(raw))
cat("Norm dimensions: "); print(dim(norm))

# Check row and column identity
print(all(rownames(raw) == rownames(norm)))
print(all(colnames(raw) == colnames(norm)))

###############################################################
# 3) PCA FUNCTION
###############################################################

run_pca <- function(mat, title){
  mat2 <- mat[complete.cases(mat), ]     # remove NA rows for PCA
  pca <- prcomp(t(mat2), scale. = TRUE)
  
  df <- data.frame(pca$x[,1:2], Condition = colnames(mat2))
  
  ggplot(df, aes(PC1, PC2, label = Condition)) +
    geom_point(size = 3) +
    geom_text(hjust = 1.2, vjust = 1.2) +
    ggtitle(title) +
    theme_bw() +
    theme(plot.title = element_text(size = 16, face = "bold"))
}

###############################################################
# 4) PCA BEFORE AND AFTER NORMALIZATION
###############################################################

p1 <- run_pca(raw,  "PCA – Raw Intensities")
p2 <- run_pca(norm, "PCA – Normalized (RUV) Intensities")

dev.new()
p1
dev.new()
p2


###############################################################
# 5) DENSITY DISTRIBUTIONS (RAW VS NORM)
###############################################################

plot_density <- function(mat, title){
  df <- as.data.frame(mat)
  df$Site <- rownames(mat)
  
  df_melt <- reshape2::melt(df, id.vars = "Site",
                            variable.name = "Sample",
                            value.name = "Intensity")
  
  ggplot(df_melt, aes(Intensity, color = Sample, group = Sample)) +
    geom_density(alpha = 0.2) +
    theme_bw() +
    ggtitle(title) +
    theme(legend.position = "none")
}

d1 <- plot_density(raw_aligned, "Density – Raw (aligned)")
d2 <- plot_density(norm_aligned, "Density – Norm (aligned)")

d1 | d2



###############################################################
# 6) MEAN–VARIANCE TREND
###############################################################

plot_mean_var <- function(mat, title){
  means <- rowMeans(mat, na.rm = TRUE)
  vars  <- rowVars(as.matrix(mat), na.rm = TRUE)
  
  df <- data.frame(mean = means, var = vars)
  
  ggplot(df, aes(mean, var)) +
    geom_point(alpha = 0.3) +
    scale_y_log10() +
    scale_x_log10() +
    ggtitle(title) +
    theme_bw()
}

mv1 <- plot_mean_var(raw, "Mean–Variance (Raw)")
mv2 <- plot_mean_var(norm, "Mean–Variance (Normalized)")

mv1 | mv2


###############################################################
# 7) CV DISTRIBUTION
###############################################################

plot_cv <- function(mat, title){
  cv <- rowSds(as.matrix(mat), na.rm = TRUE) / 
    rowMeans(as.matrix(mat), na.rm = TRUE)
  
  ggplot(data.frame(CV = cv), aes(CV)) +
    geom_histogram(bins = 60, fill = "steelblue", color = "white") +
    ggtitle(title) +
    theme_bw()
}

cv1 <- plot_cv(raw, "CV – Raw")
cv2 <- plot_cv(norm, "CV – Normalized")

cv1 | cv2



plot_cv <- function(mat, title){
  m <- rowMeans(mat, na.rm = TRUE)
  s <- rowSds(as.matrix(mat), na.rm = TRUE)
  
  cv <- s / m
  
  # Remove impossible values
  cv <- cv[is.finite(cv) & cv >= 0]
  
  ggplot(data.frame(CV = cv), aes(CV)) +
    geom_histogram(bins = 60, fill = "steelblue", color = "white") +
    ggtitle(title) +
    theme_bw()
}

cv1 <- plot_cv(raw_aligned, "CV – Raw")
cv2 <- plot_cv(norm_aligned, "CV – Normalized")

cv1 | cv2




plot_cv_log <- function(mat, title){
  m <- rowMeans(mat, na.rm = TRUE)
  s <- rowSds(as.matrix(mat), na.rm = TRUE)
  cv <- s / m
  cv <- cv[is.finite(cv) & cv >= 0]
  
  ggplot(data.frame(CV = cv), aes(CV)) +
    geom_histogram(bins = 60, fill = "steelblue", color = "white") +
    scale_x_log10() +
    ggtitle(title) +
    theme_bw()
}

plot_cv_log(raw_aligned, "CV – Raw (log scale)") |
  plot_cv_log(norm_aligned, "CV – Norm (log scale)")



###############################################################
# log2fc
###############################################################

# ALIGN RAW AND NORMALIZED MATRICES
common_rows <- intersect(rownames(raw), rownames(norm))

raw_aligned  <- raw[common_rows, , drop = FALSE]
norm_aligned <- norm[common_rows, , drop = FALSE]


cond_list <- list(
  T0      = grep("0000",     colnames(raw_aligned), value = TRUE),
  T10     = grep("0010",    colnames(raw_aligned), value = TRUE),
  T600    = grep("0600",   colnames(raw_aligned), value = TRUE),
  T1800   = grep("1800",  colnames(raw_aligned), value = TRUE)
)


compute_fc <- function(mat, groupA, groupB){
  A <- mat[, groupA, drop = FALSE]
  B <- mat[, groupB, drop = FALSE]
  
  rowMeans(B, na.rm = TRUE) - rowMeans(A, na.rm = TRUE)
}


fc_raw_10    <- compute_fc(raw_aligned,  cond_list$T0, cond_list$T10)
fc_norm_10   <- compute_fc(norm_aligned, cond_list$T0, cond_list$T10)

fc_raw_600   <- compute_fc(raw_aligned,  cond_list$T0, cond_list$T600)
fc_norm_600  <- compute_fc(norm_aligned, cond_list$T0, cond_list$T600)

fc_raw_1800  <- compute_fc(raw_aligned,  cond_list$T0, cond_list$T1800)
fc_norm_1800 <- compute_fc(norm_aligned, cond_list$T0, cond_list$T1800)




plot_fc <- function(raw_fc, norm_fc, title){
  df <- data.frame(raw = raw_fc, norm = norm_fc)
  
  ggplot(df, aes(raw, norm)) +
    geom_point(alpha = 0.2) +
    geom_smooth(method = "lm", color = "red") +
    theme_bw() +
    xlab("Raw log2FC") +
    ylab("Normalized log2FC") +
    ggtitle(title)
}

p10   <- plot_fc(fc_raw_10,   fc_norm_10,   "FC: 0s → 10s")
p600  <- plot_fc(fc_raw_600,  fc_norm_600,  "FC: 0s → 600s")
p1800 <- plot_fc(fc_raw_1800, fc_norm_1800, "FC: 0s → 1800s")

p10 | p600 | p1800



cat("Corr (0→10s):   ", cor(fc_raw_10,   fc_norm_10,   use="pairwise"), "\n")
cat("Corr (0→600s):  ", cor(fc_raw_600,  fc_norm_600,  use="pairwise"), "\n")
cat("Corr (0→1800s): ", cor(fc_raw_1800, fc_norm_1800, use="pairwise"), "\n")



###############################################################
# get numbers for differential regulated psites
###############################################################


# 1. Load initial + validation datasets


val_raw  <- read.delim(
  "D:/Research/CXCR7_platelet_analysis/SubProjects/CXCR7_validation/phosphoproteomics/data/raw_data/raw_intensities.txt",
  row.names = 1, check.names = FALSE)

val_norm <- read.delim(
  "D:/Research/CXCR7_platelet_analysis/SubProjects/CXCR7_validation/phosphoproteomics/data/processed_data/norm_intensity.txt",
  row.names = 1, check.names = FALSE)

init_raw <- read.delim(
  "D:/Research/CXCR7_platelet_analysis/SubProjects/CXCR7_initial/phosphoproteomics/data/raw_data/raw_intensities.txt",
  row.names = 1, check.names = FALSE)

init_norm <- read.delim(
  "D:/Research/CXCR7_platelet_analysis/SubProjects/CXCR7_initial/phosphoproteomics/data/processed_data/norm_intensity.txt",
  row.names = 1, check.names = FALSE)



# 2. Align rows between raw + normalized for both datasets


align_matrices <- function(raw_mat, norm_mat){
  common <- intersect(rownames(raw_mat), rownames(norm_mat))
  raw_aligned  <- raw_mat[common, , drop=FALSE]
  norm_aligned <- norm_mat[common, , drop=FALSE]
  list(raw = raw_aligned, norm = norm_aligned, rows = common)
}

init_aligned <- align_matrices(init_raw, init_norm)
val_aligned  <- align_matrices(val_raw,  val_norm)

cat("Initial dataset common phosphosites:", length(init_aligned$rows), "\n")
cat("Validation dataset common phosphosites:", length(val_aligned$rows), "\n")


# Extract aligned matrices
raw_init  <- init_aligned$raw
norm_init <- init_aligned$norm

raw_val   <- val_aligned$raw
norm_val  <- val_aligned$norm



# 3. Build condition list automatically (T0, T10, T600, T1800)


build_cond_list <- function(colnames_vec){
  list(
    T0    = grep("0000", colnames_vec, value = TRUE),
    T10   = grep("0010", colnames_vec, value = TRUE),
    T600  = grep("0600", colnames_vec, value = TRUE),
    T1800 = grep("1800", colnames_vec, value = TRUE)
  )
}

cond_init <- build_cond_list(colnames(norm_init))
cond_val  <- build_cond_list(colnames(norm_val))



# 4. Compute log2 fold changes (same function for all datasets)


compute_fc <- function(mat, groupA, groupB){
  A <- mat[, groupA, drop=FALSE]
  B <- mat[, groupB, drop=FALSE]
  rowMeans(B, na.rm=TRUE) - rowMeans(A, na.rm=TRUE)
}

# Initial experiment FC
fc_init <- list(
  T10   = compute_fc(norm_init, cond_init$T0, cond_init$T10),
  T600  = compute_fc(norm_init, cond_init$T0, cond_init$T600),
  T1800 = compute_fc(norm_init, cond_init$T0, cond_init$T1800)
)

# Validation experiment FC
fc_val <- list(
  T10   = compute_fc(norm_val, cond_val$T0, cond_val$T10),
  T600  = compute_fc(norm_val, cond_val$T0, cond_val$T600),
  T1800 = compute_fc(norm_val, cond_val$T0, cond_val$T1800)
)

cat("Computed log2FC for initial + validation datasets.\n")





###############################################################
# RUVphospho
###############################################################
library("EDASeq")
library("ggplot2")
library("ggpubr")
library("limma")
library("matrixStats")
library("SummarizedExperiment")
library("tidyr")
library("PhosR")


# Load scaled + imputed intensities
scaled_file <- "D:/Research/CXCR7_platelet_analysis/SubProjects/CXCR7_validation/phosphoproteomics/data/processed_data/intensities_imputed_scaled.txt"

scaled_mat <- as.matrix(read.delim(scaled_file, row.names = 1, check.names = FALSE))

scaled_log2 <- log2(scaled_mat)


ctl_file <- "D:/Research/CXCR7_platelet_analysis/SubProjects/CXCR7_validation/phosphoproteomics/data/processed_data/empirical_topall.txt"
empirical_topall <- scan(ctl_file, what = character())


ctl2 <- which(rownames(scaled_log2) %in% empirical_topall)
length(ctl2)

grps <- factor(gsub("(.*?)(_\\d+)$", "\\1", colnames(scaled_log2)))  # adjust pattern if needed
design <- model.matrix(~ grps - 1)




ruv_out <- RUVphospho(scaled_log2,
                      M = design,
                      k = 16,
                      ctl = ctl2)

ruv_out <- RUVphospho(scaled_log2,
                      M = design,
                      k = 16,
                      ctl = ctl2)

ruv_full <- PhosR:::RUV(
  mat = scaled_log2,
  M   = design,
  ctl = ctl2,
  k   = 16,
  m   = 1.6,
  s   = 0.6,
  keepImpute = FALSE
)



W <- ruv_full$W
dim(W)





library(ruv)

ruv_out <- RUVIII(
  Y   = scaled_log2_t,   # <— FIX: transpose used here
  M   = design,
  ctl = ctl2,            # careful: ctl2 indexes now refer to columns (features)
  k   = 16
)




eigen_vals <- eigen(t(W) %*% W)$values

plot(eigen_vals, type="b", pch=19,
     xlab="RUV Factor Index (1–16)",
     ylab="Eigenvalue",
     main="Variance Explained by RUV Unwanted Factors")


pca_raw <- prcomp(t(scaled_log2), scale.=TRUE)

cor_mat <- cor(W, pca_raw$x[, 1:5])    # correlation with top 5 PCs

library(corrplot)
corrplot(cor_mat, method="color", tl.col="black",
         main="Correlation of RUV Factors (W) with Raw PCA")


library(pheatmap)

annotation <- data.frame(Group = grps)
rownames(annotation) <- colnames(scaled_log2)

pheatmap(W,
         annotation_row = annotation,
         clustering_method = "ward.D2",
         scale = "row",
         main = "Heatmap of RUV Unwanted Factors (W Matrix)",
         fontsize_row = 6)



############################################
## Publication-Ready Multi-Panel Figure Code
############################################

library(ggplot2)
library(reshape2)
library(matrixStats)
library(patchwork)

### -------- 1) Density plots (Raw vs Normalized) -------- ###

plot_density <- function(mat, title) {
  df <- as.data.frame(mat)
  df$Site <- rownames(mat)
  df_melt <- melt(df, id.vars = "Site",
                  variable.name = "Sample",
                  value.name = "Intensity")
  
  ggplot(df_melt, aes(Intensity, color = Sample, group = Sample)) +
    geom_density(alpha = 0.25, linewidth = 0.7) +
    theme_bw(base_size = 12) +
    ggtitle(title) +
    theme(legend.position = "none")
}

p_density_raw  <- plot_density(raw_aligned,  "Density – Raw")
p_density_norm <- plot_density(norm_aligned, "Density – Normalized")


### -------- 2) Mean–Variance plots -------- ###

plot_mean_var <- function(mat, title) {
  means <- rowMeans(mat, na.rm = TRUE)
  vars  <- rowVars(as.matrix(mat), na.rm = TRUE)
  
  df <- data.frame(mean = means, var = vars)
  
  ggplot(df, aes(mean, var)) +
    geom_point(alpha = 0.3, size = 0.7) +
    scale_x_log10() + scale_y_log10() +
    theme_bw(base_size = 12) +
    ggtitle(title)
}

p_mv_raw  <- plot_mean_var(raw_aligned,  "Mean–Variance – Raw")
p_mv_norm <- plot_mean_var(norm_aligned, "Mean–Variance – Normalized")


### -------- 3) CV distributions (log-scale) -------- ###

plot_cv_log <- function(mat, title){
  m <- rowMeans(mat, na.rm = TRUE)
  s <- rowSds(as.matrix(mat), na.rm = TRUE)
  cv <- s / m
  cv <- cv[is.finite(cv) & cv >= 0]
  
  ggplot(data.frame(CV = cv), aes(CV)) +
    geom_histogram(bins = 60, fill = "steelblue", color = "white") +
    scale_x_log10() +
    theme_bw(base_size = 12) +
    ggtitle(title)
}

p_cv_raw  <- plot_cv_log(raw_aligned,  "CV – Raw (log scale)")
p_cv_norm <- plot_cv_log(norm_aligned, "CV – Normalized (log scale)")


### -------- 4) Fold-change comparison -------- ###

plot_fc <- function(raw_fc, norm_fc, title){
  df <- data.frame(raw = raw_fc, norm = norm_fc)
  
  ggplot(df, aes(raw, norm)) +
    geom_point(alpha = 0.25, size = 0.7) +
    geom_smooth(method = "lm", color = "red", linewidth = 0.7) +
    theme_bw(base_size = 12) +
    xlab("Raw log2FC") +
    ylab("Normalized log2FC") +
    ggtitle(title)
}

p_fc_10   <- plot_fc(fc_raw_10,   fc_norm_10,   "FC: 0 → 10 s")
p_fc_600  <- plot_fc(fc_raw_600,  fc_norm_600,  "FC: 0 → 600 s")
p_fc_1800 <- plot_fc(fc_raw_1800, fc_norm_1800, "FC: 0 → 1800 s")


### -------- 5) Combine into multi-panel publication figure -------- ###

final_plot <- 
  (p_density_raw + p_density_norm) /
  (p_mv_raw + p_mv_norm) /
  (p_cv_raw + p_cv_norm) /
  (p_fc_10 + p_fc_600 + p_fc_1800)

final_plot


