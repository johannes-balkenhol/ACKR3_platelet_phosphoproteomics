############################################
#### script for annotation of phosphoproteom data 
## Kinase Substrate enrichment analysis
## with detected-kinase filtering
############################################

library(KSEAapp)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)


############################################################
## Load Data
############################################################
setwd("D:/Research/CXCR7_platelet_analysis/SubProjects/CXCR7_validation/phosphoproteomics/scripts/R")

KSData <- read.csv("../../analysis/KSEA/PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv")

############################################################
## Detect kinases present in your proteome
############################################################

extract_genes <- function(x) sapply(strsplit(x, ";"), `[`, 2)
detected_genes <- extract_genes(rownames(norm_intensity))

detected_kinases <- intersect(unique(KSData$GENE), detected_genes)
cat("Detected kinases:", length(detected_kinases), "\n")

## Filter KSData BEFORE running KSEA
KSData_filtered <- KSData %>% filter(GENE %in% detected_kinases)

cat("Kinases remaining after KSData filtering:",
    length(unique(KSData_filtered$GENE)), "\n")

############################################################
## Input lists
############################################################

dfs_input = list(
  top.filter.10, top.filter.600, top.filter.1800,
  top.filter.10.dmso.vs.0s, top.filter.600.dmso.vs.0s, top.filter.1800.dmso.vs.0s,
  top.filter.10.cxcr7.vs.0s, top.filter.600.cxcr7.vs.0s, top.filter.1800.cxcr7.vs.0s
)

names_input = c(
  "10","600","1800",
  "10.dmso.vs.0s","600.dmso.vs.0s","1800.dmso.vs.0s",
  "10.cxcr7.vs.0s","600.cxcr7.vs.0s","1800.cxcr7.vs.0s"
)

############################################################
## RUN KSEA USING FILTERED KSData
############################################################

for (i in 1:length(dfs_input)) {
  
  my_df1 <- dfs_input[[i]]
  
  my_df1$Peptide <- sapply(strsplit(rownames(my_df1), ";"), "[[", 4)
  my_df1$FC <- 2^my_df1$logFC
  
  ## required KSEA format
  my_df1 <- my_df1[c(1,2,7,6,5,8)]
  colnames(my_df1) <- c("Protein", "Gene", "Peptide", "Residue.Both", "p", "FC")
  
  PX <- my_df1
  PX$Residue.Both <- str_replace_all(PX$Residue.Both, fixed("|"), ";")
  rownames(PX) <- NULL
  
  condition <- names_input[[i]]
  message("\n--- Running: ", condition, " ---")
  
  ############################################################
  ## KS table — using filtered KSData
  ############################################################
  
  KSData.dataset <- KSEA.KS_table(KSData_filtered, PX, 
                                  NetworKIN = FALSE, NetworKIN.cutoff = 5)
  
  write.table(KSData.dataset,
              file = paste0("../../analysis/KSEA/", condition, "_Kinase_Substrate_Links.csv"),
              sep = ",", quote = FALSE, row.names = FALSE)
  
  ############################################################
  ## KSEA scores — using filtered KSData
  ############################################################
  
  Scores <- KSEA.Scores(KSData_filtered, PX, 
                        NetworKIN = FALSE, NetworKIN.cutoff = 5)
  
  assign(paste0("Scores", condition), Scores)
  
  write.table(Scores,
              file = paste0("../../analysis/KSEA/", condition, "_KSEA_Kinase_Scores-validation.csv"),
              sep = ",", quote = FALSE, row.names = FALSE)
  
  ############################################################
  ## Barplot (default KSEA barplot)
  ############################################################
  
  tiff(filename = paste0("../../analysis/KSEA/", condition, "_barplot.tiff"),
       width = 5 * 300, height = 5 * 300, res = 300, compression = "lzw")
  
  KSEA.Barplot(KSData_filtered, PX, 
               NetworKIN = FALSE, NetworKIN.cutoff = 5, 
               m.cutoff = 5, p.cutoff = 0.05, export = FALSE)
  
  dev.off()
}

############################################################
## Heatmap (example with 3 samples)
############################################################

KSEA.Heatmap(list(Scores10, Scores600, Scores1800),
             sample.labels = c("10", "600", "1800"),
             m.cutoff = 5, stats = "p.value",
             p.cutoff = 0.05, sample.cluster = FALSE)





############################################################
## Collect all scores into a single matrix
############################################################

all_scores <- list()

for (nm in names_input) {
  all_scores[[nm]] <- get(paste0("Scores", nm))
}



############################################################
## Collect all scores into a single matrix
############################################################

# List of all kinases across all conditions
all_kinases <- sort(unique(unlist(
  lapply(all_scores, function(df) df$Kinase.Gene)
)))

heatmap_mat <- matrix(NA,
                      nrow = length(all_kinases),
                      ncol = length(all_scores),
                      dimnames = list(all_kinases, names(all_scores)))

signif_mat <- matrix("",
                     nrow = length(all_kinases),
                     ncol = length(all_scores),
                     dimnames = list(all_kinases, names(all_scores)))

for (nm in names(all_scores)) {
  
  df <- all_scores[[nm]]
  
  # fill z-scores
  heatmap_mat[df$Kinase.Gene, nm] <- df$z.score
  
  # fill stars
  signif_idx <- df$Kinase.Gene[df$p.value < 0.05]
  signif_mat[signif_idx, nm] <- "*"
}

############################################################
##  heatmap all
############################################################


library(pheatmap)

# Color scale centered at 0 (blue to red)
breaks <- seq(-4, 4, length.out = 100)
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(99)

tiff("../../analysis/KSEA/KSEA_Heatmap_DETECTED-validation.tiff",
     width = 300*10,
     height = 300*14,
     res = 300,
     compression = "lzw")

pheatmap(heatmap_mat,
         cluster_rows = TRUE,
         cluster_cols = FALSE,   # no sample clustering
         color = colors,
         breaks = breaks,
         na_col = "grey90",
         display_numbers = signif_mat,
         number_color = "black",
         fontsize = 10,
         fontsize_row = 8,
         main = "KSEA Z-score Heatmap (Detected Kinases Only)\n* = p < 0.05")

dev.off()


#KSEA.Heatmap(list(Scores10, Scores600, Scores1800),
#             sample.labels = c("10", "600", "1800"),
#             m.cutoff = 5, stats="p.value",p.cutoff=0.05, sample.cluster=F)



############################################################
##  heatmap all cxcr7 vs dmso
############################################################
cols_raw <- c("10", "600", "1800")

heatmap_raw <- heatmap_mat[, cols_raw, drop = FALSE]
signif_raw  <- signif_mat[, cols_raw, drop = FALSE]

keep_rows <- rowSums(!is.na(heatmap_raw)) > 0
heatmap_raw <- heatmap_raw[keep_rows, ]
signif_raw  <- signif_raw[keep_rows, ]


colors <- colorRampPalette(c("navy", "white", "firebrick3"))(99)
breaks <- seq(-4, 4, length.out = 100)

tiff("../../analysis/KSEA/KSEA_Heatmap_raw_10_600_1800.tiff",
     width = 300*10,
     height = 300*14,
     res = 300,
     compression = "lzw")

pheatmap(
  heatmap_raw,
  cluster_rows = TRUE,
  cluster_cols = TRUE,     # allow clustering of 10/600/1800
  color = colors,
  breaks = breaks,
  na_col = "grey90",
  display_numbers = signif_raw,
  number_color = "black",
  fontsize = 10,
  fontsize_row = 8,
  main = "KSEA Z-score Heatmap (10 / 600 / 1800)\n* = p < 0.05"
)

dev.off()



############################################################
##  heatmap all cxcr7 vs dmso top
############################################################
# extract the 1800 column
z1800 <- heatmap_mat[, "1800"]

# keep only non-NA values
z1800_clean <- z1800[!is.na(z1800)]

# get names of top 20 kinases
top20 <- names(sort(abs(z1800_clean), decreasing = TRUE))[1:30]

cols_raw <- c("10", "600", "1800")

heatmap_top20 <- heatmap_mat[top20, cols_raw, drop = FALSE]
signif_top20  <- signif_mat[top20, cols_raw, drop = FALSE]

colors <- colorRampPalette(c("navy", "white", "firebrick3"))(99)
breaks <- seq(-4, 4, length.out = 100)

tiff("../../analysis/KSEA/KSEA_Heatmap_top20_1800.tiff",
     width = 100*10,
     height = 300*12,
     res = 300,
     compression = "lzw")

pheatmap(
  heatmap_top20,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colors,
  breaks = breaks,
  na_col = "grey90",
  display_numbers = signif_top20,
  number_color = "black",
  fontsize = 10,
  fontsize_row = 9,
  main = "Top 20 Kinases at 1800 (|Z|)\nHeatmap for 10 / 600 / 1800"
)

dev.off()






############################################################
## Harmonize old tables
############################################################


harmonize_old_table <- function(df) {
  
  # The annotation is in the 8th column
  annot <- df[[8]]
  
  # Split annotation
  split_annot <- strsplit(annot, ";")
  
  df$uniprot_id <- sapply(split_annot, `[`, 1)
  df$name       <- sapply(split_annot, `[`, 2)
  df$PSite      <- sapply(split_annot, `[`, 3)
  df$Peptide    <- sapply(split_annot, `[`, 4)
  
  # Harmonize column names to match new tables
  df_harmonized <- df %>%
    transmute(
      uniprot_id = uniprot_id,
      name       = name,
      Average    = AveExpr,
      logFC      = logFC,
      PValue     = P.Value,
      PSite      = PSite
    )
  
  return(df_harmonized)
}




############################################################
## Load OLD CXCR7 stimulation-only tables
############################################################

old_files <- list.files(
  "SubProjects/CXCR7_initial/phosphoproteomics/data/processed_data",
  pattern = "*.cxcr7.vs.0s.txt",
  full.names = TRUE
)

old_tables <- lapply(old_files, read.delim)
names(old_tables) <- gsub(".txt$", "", basename(old_files))


old_tables_harmonized <- lapply(old_tables, harmonize_old_table)

###########################################################
## Merge NEW + OLD tables
###########################################################

dfs_input <- c(
  dfs_input,            # your existing 9 tables
  old_tables_harmonized # the harmonized old CXCR7 tables
)

names_input <- c(
  names_input,
  names(old_tables_harmonized)
)










############################################
#### Load Libraries
############################################

library(KSEAapp)
library(dplyr)
library(tidyr)
library(stringr)
library(pheatmap)


############################################
#### Load KSData
############################################

KSData <- read.csv("../../analysis/KSEA/PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv")


############################################
#### DETECTED KINASES FROM PROTEOME
############################################

extract_genes <- function(x) sapply(strsplit(x, ";"), `[`, 2)
detected_genes <- extract_genes(rownames(norm_intensity))

detected_kinases <- intersect(unique(KSData$GENE), detected_genes)
cat("Detected kinases:", length(detected_kinases), "\n")

KSData_filtered <- KSData %>% filter(GENE %in% detected_kinases)
cat("Kinases remaining after KSData filtering:",
    length(unique(KSData_filtered$GENE)), "\n")




############################################
#### Load and HARMONIZE OLD CXCR7 TABLES
############################################

old_path <- "../../../../../SubProjects/CXCR7_initial/phosphoproteomics/data/processed_data"

old_files <- list.files(
  old_path,
  pattern = "cxcr7",
  full.names = TRUE,
  ignore.case = TRUE
)

old_files[1:3]

old_tables <- lapply(old_files, read.delim, stringsAsFactors = FALSE)
names(old_tables) <- gsub(".txt$", "", basename(old_files))

old_tables[1:3]



harmonize_old <- function(df) {
  annot <- df[[8]]  # Annotation column
  parts <- strsplit(annot, ";")
  
  df$uniprot_id <- sapply(parts, `[`, 1)
  df$name       <- sapply(parts, `[`, 2)
  df$PSite      <- sapply(parts, `[`, 3)
  df$Peptide    <- sapply(parts, `[`, 4)
  
  out <- df %>%
    transmute(
      uniprot_id,
      name,
      Average = AveExpr,
      logFC,
      PValue = P.Value,
      PSite
    )
  
  return(out)
}

colnames(old_tables)
old_tables_harmonized <- lapply(old_tables, harmonize_old)
colnames(old_tables_harmonized)


############################################
#### NEW TABLES (already in correct format)
############################################

dfs_new <- list(
  top.filter.10, top.filter.600, top.filter.1800,
  top.filter.10.dmso.vs.0s, top.filter.600.dmso.vs.0s, top.filter.1800.dmso.vs.0s,
  top.filter.10.cxcr7.vs.0s, top.filter.600.cxcr7.vs.0s, top.filter.1800.cxcr7.vs.0s
)

names_new <- c(
  "10","600","1800",
  "10.dmso.vs.0s","600.dmso.vs.0s","1800.dmso.vs.0s",
  "10.cxcr7.vs.0s","600.cxcr7.vs.0s","1800.cxcr7.vs.0s"
)


############################################
#### MERGE NEW + OLD INTO ONE LIST + Harmonize names
############################################

dfs_input <- c(dfs_new, old_tables_harmonized)
names_input <- c(names_new, names(old_tables_harmonized))


names_input_clean <- c(
  "cxcr7VSdmso.10",
  "cxcr7VSdmso.600",
  "cxcr7VSdmso.1800",
  
  "dmso.10.validation",
  "dmso.600.validation",
  "dmso.1800.validation",
  
  "cxcr7.10.validation",
  "cxcr7.600.validation",
  "cxcr7.1800.validation",
  
  "cxcr7.10.initial",
  "cxcr7.1800.initial",
  "cxcr7.30.initial",
  "cxcr7.300.initial",
  "cxcr7.60.initial",
  "cxcr7.600.initial",
  "cxcr7.900.initial"
)


names_input_clean

names(dfs_input) <- names_input_clean


############################################
#### RUN KSEA FOR ALL TABLES (with clean names)
############################################

Scores_list <- list()

for (i in seq_along(dfs_input)) {
  
  df <- dfs_input[[i]]
  condition <- names(dfs_input)[i]   # <-- use the CLEAN NAME
  message("\n--- Running: ", condition, " ---")
  
  ### Build PX table for KSEA
  df$Peptide <- sapply(strsplit(df$PSite, ";"), `[`, 4)
  df$FC <- 2^df$logFC
  
  PX <- df %>%
    transmute(
      Protein = uniprot_id,
      Gene = name,
      Peptide,
      Residue.Both = str_replace_all(PSite, "\\|", ";"),
      p = PValue,
      FC
    )
  
  rownames(PX) <- NULL
  
  ### K-S Table
  KS_table <- KSEA.KS_table(KSData_filtered, PX,
                            NetworKIN = FALSE, NetworKIN.cutoff = 5)
  
  write.csv(KS_table,
            paste0("../../analysis/KSEA/", condition, "_KS_links.csv"),
            row.names = FALSE)
  
  ### KSEA Scores
  Scores <- KSEA.Scores(KSData_filtered, PX,
                        NetworKIN = FALSE, NetworKIN.cutoff = 5)
  
  Scores_list[[condition]] <- Scores
  
  write.csv(Scores,
            paste0("../../analysis/KSEA/", condition, "_KSEA_scores.csv"),
            row.names = FALSE)
  
  ### Barplot
  tiff(filename = paste0("../../analysis/KSEA/", condition, "_barplot.tiff"),
       width = 300*5, height = 300*5, res = 300)
  
  KSEA.Barplot(KSData_filtered, PX,
               NetworKIN = FALSE,
               NetworKIN.cutoff = 5,
               m.cutoff = 5,
               p.cutoff = 0.05,
               export = FALSE)
  
  dev.off()
}



############################################
#### COMBINE ALL SCORES INTO HEATMAP
############################################

# collect all kinases across all conditions
all_kinases <- sort(unique(unlist(
  lapply(Scores_list, function(df) df$Kinase.Gene)
)))

# initialize matrices
heatmap_mat <- matrix(
  NA,
  nrow = length(all_kinases),
  ncol = length(Scores_list),
  dimnames = list(all_kinases, names(Scores_list))
)

signif_mat <- matrix(
  "",
  nrow = length(all_kinases),
  ncol = length(Scores_list),
  dimnames = list(all_kinases, names(Scores_list))
)

# fill matrices
for (nm in names(Scores_list)) {
  df <- Scores_list[[nm]]
  
  # fill z-scores
  heatmap_mat[df$Kinase.Gene, nm] <- df$z.score
  
  # significance stars
  signif_idx <- df$Kinase.Gene[df$p.value < 0.05]
  signif_mat[signif_idx, nm] <- "*"
}




############################################
#### filter NA
############################################
heatmap_mat_clean <- heatmap_mat[complete.cases(heatmap_mat), ]
signif_mat_clean  <- signif_mat[rownames(heatmap_mat_clean), ]


# remove any row containing ANY NA
keep_rows <- complete.cases(heatmap_mat)

heatmap_mat_clean <- heatmap_mat[keep_rows, ]
signif_mat_clean  <- signif_mat[keep_rows, ]

# verify
print(which(!is.finite(heatmap_mat_clean), arr.ind = TRUE))



############################################
#### HEATMAP all
############################################


colors <- colorRampPalette(c("navy", "white", "firebrick3"))(99)
breaks <- seq(-4, 4, length.out = 100)

tiff("../../analysis/KSEA/KSEA_Heatmap_DETECTED.tiff",
     width = 300*14, height = 300*18, res = 300)

pheatmap(
  heatmap_mat_clean,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colors,
  breaks = breaks,
  na_col = "grey90",
  display_numbers = signif_mat_clean,
  number_color = "black",
  fontsize = 10,
  fontsize_row = 7,
  main = "KSEA Z-score Heatmap (Detected Kinases Only)\n* = p < 0.05"
)

dev.off()

############################################
#### HEATMAP dmsoVScxcr7
############################################

colors <- colorRampPalette(c("navy", "white", "firebrick3"))(99)
breaks <- seq(-4, 4, length.out = 100)


cols_A <- c("cxcr7VSdmso.10", "cxcr7VSdmso.600", "cxcr7VSdmso.1800")

heatmap_A <- heatmap_mat_clean[, cols_A, drop = FALSE]
signif_A  <- signif_mat_clean[, cols_A, drop = FALSE]

tiff("../../analysis/KSEA/KSEA_Heatmap_A_VSdmso.tiff",
     width = 300*10, height = 300*14, res = 300)

pheatmap(
  heatmap_A,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colors,
  breaks = breaks,
  na_col = "grey90",
  display_numbers = signif_A,
  number_color = "black",
  fontsize = 10,
  fontsize_row = 7,
  main = "KSEA Z-score (CXCR7 vs DMSO)\n* = p < 0.05"
)

dev.off()


############################################
#### HEATMAP initial vs validation
############################################


cols_B <- c(
  "cxcr7.10.validation", "cxcr7.600.validation", "cxcr7.1800.validation",
  "cxcr7.10.initial", "cxcr7.600.initial", "cxcr7.1800.initial"
)

heatmap_B <- heatmap_mat_clean[, cols_B, drop = FALSE]
signif_B  <- signif_mat_clean[, cols_B, drop = FALSE]

tiff("../../analysis/KSEA/KSEA_Heatmap_B_validation_initial.tiff",
     width = 300*12, height = 300*18, res = 300)

pheatmap(
  heatmap_B,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colors,
  breaks = breaks,
  na_col = "grey90",
  display_numbers = signif_B,
  number_color = "black",
  fontsize = 10,
  fontsize_row = 7,
  main = "KSEA Z-score (Validation + Initial)\n* = p < 0.05"
)

dev.off()




###############################################################
## PUBLICATION-READY KSEA HEATMAP
## CXCR7 vs DMSO across timepoints
###############################################################

library(pheatmap)
library(RColorBrewer)

## ------------------------------------------------------------
## 1. Define color scheme and breaks
## ------------------------------------------------------------

# Blue-White-Red color palette (publication standard)
colors <- colorRampPalette(c("navy", "white", "firebrick3"))(99)
breaks <- seq(-4, 4, length.out = 100)

## ------------------------------------------------------------
## 2. Select columns for CXCR7 vs DMSO comparison
## ------------------------------------------------------------

cols_A <- c("cxcr7VSdmso.10", "cxcr7VSdmso.600", "cxcr7VSdmso.1800")

# Extract heatmap data and significance matrix
heatmap_A <- heatmap_mat_clean[, cols_A, drop = FALSE]
signif_A  <- signif_mat_clean[, cols_A, drop = FALSE]

# Rename columns for publication (cleaner labels)
colnames(heatmap_A) <- c("10s", "600s", "1800s")
colnames(signif_A) <- c("10s", "600s", "1800s")

## ------------------------------------------------------------
## 3. Create publication-ready heatmap
## ------------------------------------------------------------

# Display in R
pheatmap(
  heatmap_A,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colors,
  breaks = breaks,
  na_col = "grey90",
  display_numbers = signif_A,
  number_color = "black",
  fontsize = 14,
  fontsize_row = 10,
  fontsize_col = 12,
  fontsize_number = 10,
  cellwidth = 30,
  cellheight = 14,
  main = "KSEA Z-score (CXCR7 vs DMSO)\n* = p < 0.05",
  border_color = "grey80",
  legend_breaks = c(-4, -2, 0, 2, 4),
  legend_labels = c("-4", "-2", "0", "2", "4")
)

## ------------------------------------------------------------
## 4. Save as high-resolution TIFF
## ------------------------------------------------------------

tiff("../../analysis/KSEA/KSEA_Heatmap_CXCR7vsDMSO.tiff",
     width = 8, height = 12, units = "in", res = 300, compression = "lzw")

pheatmap(
  heatmap_A,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colors,
  breaks = breaks,
  na_col = "grey90",
  display_numbers = signif_A,
  number_color = "black",
  fontsize = 14,
  fontsize_row = 10,
  fontsize_col = 12,
  fontsize_number = 10,
  cellwidth = 30,
  cellheight = 14,
  main = "KSEA Z-score (CXCR7 vs DMSO)\n* = p < 0.05",
  border_color = "grey80",
  legend_breaks = c(-4, -2, 0, 2, 4),
  legend_labels = c("-4", "-2", "0", "2", "4")
)

dev.off()

## ------------------------------------------------------------
## 5. Save as PDF (vector format)
## ------------------------------------------------------------

pdf("../../analysis/KSEA/KSEA_Heatmap_CXCR7vsDMSO.pdf",
    width = 8, height = 12)

pheatmap(
  heatmap_A,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colors,
  breaks = breaks,
  na_col = "grey90",
  display_numbers = signif_A,
  number_color = "black",
  fontsize = 14,
  fontsize_row = 10,
  fontsize_col = 12,
  fontsize_number = 10,
  cellwidth = 30,
  cellheight = 14,
  main = "KSEA Z-score (CXCR7 vs DMSO)\n* = p < 0.05",
  border_color = "grey80",
  legend_breaks = c(-4, -2, 0, 2, 4),
  legend_labels = c("-4", "-2", "0", "2", "4")
)

dev.off()

## ------------------------------------------------------------
## 6. Save as PNG (for presentations)
## ------------------------------------------------------------

png("../../analysis/KSEA/KSEA_Heatmap_CXCR7vsDMSO.png",
    width = 8, height = 12, units = "in", res = 300)

pheatmap(
  heatmap_A,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colors,
  breaks = breaks,
  na_col = "grey90",
  display_numbers = signif_A,
  number_color = "black",
  fontsize = 14,
  fontsize_row = 10,
  fontsize_col = 12,
  fontsize_number = 10,
  cellwidth = 30,
  cellheight = 14,
  main = "KSEA Z-score (CXCR7 vs DMSO)\n* = p < 0.05",
  border_color = "grey80",
  legend_breaks = c(-4, -2, 0, 2, 4),
  legend_labels = c("-4", "-2", "0", "2", "4")
)

dev.off()

cat("\n✓ KSEA heatmap saved successfully!\n")
cat("  Formats: TIFF (publication), PDF (vector), PNG (presentation)\n")
cat("  Location: D:/Research/CXCR7_platelet_analysis/figures/Supplementary_Figures/\n")
