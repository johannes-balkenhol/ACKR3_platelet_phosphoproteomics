############################################
#### Clean KSEA workflow (detected kinases only)
####  → KSEA tables
####  → filtered results
####  → barplots
############################################

library(KSEAapp)
library(dplyr)
library(stringr)

############################################################
## 1) Load data
############################################################

KSData <- read.csv("../../analysis/KSEA/PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv")

# Extract gene symbols from norm_intensity rownames
extract_genes <- function(x) sapply(strsplit(x, ";"), `[`, 2)
detected_genes <- extract_genes(rownames(norm_intensity))

# Intersect: detected kinases present in KSData
detected_kinases <- intersect(unique(KSData$GENE), detected_genes)
cat("Detected kinases:", length(detected_kinases), "\n")

############################################################
## 2) Helper: filter to detected kinases
############################################################

filter_ksea_to_detected <- function(scores) {
  scores %>% filter(Kinase.Gene %in% detected_kinases)
}

############################################################
## 3) Simple barplot for detected kinases only  (FIXED)
############################################################

plot_ksea_bar <- function(df, title, outfile) {
  
  df <- df %>%
    filter(!is.na(z.score)) %>%
    arrange(desc(z.score))
  
  if (nrow(df) == 0) {
    message("No kinases to plot for: ", title)
    return(NULL)
  }
  
  # FIX 1: ensure character vector (prevents barplot error)
  df$Kinase.Gene <- as.character(df$Kinase.Gene)
  
  # FIX 2: protect length consistency
  n <- min(length(df$z.score), length(df$Kinase.Gene))
  zvals  <- df$z.score[1:n]
  labels <- df$Kinase.Gene[1:n]
  
  # Barplot
  tiff(outfile, width = 5*300, height = 5*300, res = 300, compression = "lzw")
  
  barplot(
    zvals,
    names.arg = labels,
    las = 2,
    cex.names = 0.6,
    main = paste("KSEA (Detected) -", title)
  )
  
  abline(h = c(-2, 2), col = "red", lty = 2)
  dev.off()
}

############################################################
## 4) Input lists
############################################################

dfs_input <- list(
  top.filter.10,
  top.filter.600,
  top.filter.1800,
  top.filter.10.dmso.vs.0s,
  top.filter.600.dmso.vs.0s,
  top.filter.1800.dmso.vs.0s,
  top.filter.10.cxcr7.vs.0s,
  top.filter.600.cxcr7.vs.0s,
  top.filter.1800.cxcr7.vs.0s
)

names_input <- c(
  "10","600","1800",
  "10.dmso.vs.0s","600.dmso.vs.0s","1800.dmso.vs.0s",
  "10.cxcr7.vs.0s","600.cxcr7.vs.0s","1800.cxcr7.vs.0s"
)

############################################################
## 5) MAIN KSEA LOOP
############################################################


scores_filtered_list <- list()

for (i in seq_along(dfs_input)) {
  
  condition <- names_input[i]
  message("\n---- Running:", condition, "----")
  
  df <- dfs_input[[i]]
  
  ## Build PX table
  df$Peptide <- sapply(strsplit(rownames(df), ";"), `[`, 4)
  df$FC      <- 2^df$logFC
  
  PX <- df %>%
    transmute(
      Protein = uniprot_id,
      Gene = name,
      Peptide = Peptide,
      Residue.Both = str_replace_all(PSite, "\\|", ";"),
      p = PValue,
      FC = FC
    )
  
  ## KS table (filtered KSData)
  KS_dataset <- KSEA.KS_table(KSData_filtered, PX, NetworKIN = FALSE)
  write.csv(KS_dataset,
            paste0("../../analysis/KSEA/", condition, "_KS_links.csv"),
            row.names = FALSE)
  
  ## KSEA scores (filtered KSData)
  Scores <- KSEA.Scores(KSData_filtered, PX, NetworKIN = FALSE)
  write.csv(Scores,
            paste0("../../analysis/KSEA/", condition, "_KSEA_ALL.csv"),
            row.names = FALSE)
  
  ## NO need to filter again — already filtered upstream!
  scores_filtered_list[[condition]] <- Scores
  
  ## Barplot (Scores is already filtered!)
  plot_ksea_bar(
    Scores,
    condition,
    paste0("../../analysis/KSEA/", condition, "_KSEA_barplot_DETECTED.tiff")
  )
  
  message("Kinases in filtered KSData:", 
          length(unique(KSData_filtered$GENE)),
          "  → KSEA output:", nrow(Scores))
}


message("\n===== DONE: KSEA on detected kinases =====\n")



############################################################
## PUBLICATION HEATMAP (Detected kinases only)
## - No column clustering
## - Rows clustered
## - Using absolute Enrichment score |Enrichment|
## - Significance (p < 0.05) marked with *
## - Clean kinase names
## - Output: hetamp_KSEA_FILTERED.tiff
############################################################

library(pheatmap)
library(dplyr)
library(stringr)

# Collect all detected kinases across conditions
all_kinases <- sort(unique(unlist(
  lapply(scores_filtered_list, function(s) s$Kinase.Gene)
)))

# Initialize matrices
heatmap_mat <- matrix(
  NA,
  nrow = length(all_kinases),
  ncol = length(scores_filtered_list),
  dimnames = list(all_kinases, names(scores_filtered_list))
)

signif_mat <- matrix(
  "",
  nrow = length(all_kinases),
  ncol = length(scores_filtered_list),
  dimnames = list(all_kinases, names(scores_filtered_list))
)

# Fill matrices
for (nm in names(scores_filtered_list)) {
  
  df <- scores_filtered_list[[nm]]
  
  # Ensure clean kinase names (remove line breaks)
  df$Kinase.Gene <- gsub("\n|\r", "", df$Kinase.Gene)
  
  # Use Abs(Enrichment) as requested
  heatmap_mat[df$Kinase.Gene, nm] <- abs(df$Enrichment)
  
  # Mark significance
  signif_mat[df$Kinase.Gene[df$p.value < 0.05], nm] <- "*"
}

# OPTIONAL: remove all-NA rows
row_keep <- rowSums(!is.na(heatmap_mat)) > 0
heatmap_mat <- heatmap_mat[row_keep, ]
signif_mat <- signif_mat[row_keep, ]

# Custom color palette (blue → white → red)
heat_colors <- colorRampPalette(c("#313695", "#74ADD1", "#FFFFBF", "#F46D43", "#A50026"))(100)

tiff("../../analysis/KSEA/hetamp_KSEA_FILTERED.tiff",
     width = 9 * 300,
     height = 13 * 300,
     res = 300,
     compression = "lzw")

pheatmap(
  mat               = heatmap_mat,
  color             = heat_colors,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,        # <-- NO COLUMN CLUSTERING
  display_numbers   = signif_mat,   # <-- SIGNIFICANCE STARS
  number_color      = "black",
  fontsize          = 10,
  fontsize_number   = 8,
  border_color      = NA,
  na_col            = "grey90",
  main              = "Kinase Enrichment Heatmap (|Enrichment|, detected kinases only)"
)

dev.off()

message("===== PUBLICATION HEATMAP COMPLETE =====")
