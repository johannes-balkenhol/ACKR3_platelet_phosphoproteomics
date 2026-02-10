# Johannes Balkenhol, MIT License

#### Script Overview ####
# This script reconstructs an interaction network from IntAct data, filters interactions for specific proteins,
# and prepares the network for export to Cytoscape. It also includes functions for collapsing datasets
# with different strategies, with a preferred method set as the default.

#### Inputs ####
# 1. Interaction data files (e.g., intact.txt, intact_negative.txt).
# 2. Phosphoproteomics and proteomics data frames (e.g., top.all.init.cxcr7.vs.0s).
# 3. Custom protein lists.

#### Outputs ####
# 1. Processed interaction data (e.g., interactions.intact.all.refined.platelets.collapsed.omni.txt).
# 2. Node and edge files for Cytoscape export (e.g., node4.txt, edge.txt).

#### Required Packages ####
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("stringr")) install.packages("stringr")
if (!require("dplyr")) install.packages("dplyr")
if (!require("readr")) install.packages("readr")
if (!require("httr")) install.packages("httr")
if (!require("jsonlite")) install.packages("jsonlite")
if (!require("VennDiagram")) install.packages("VennDiagram")
if (!require("grid")) install.packages("grid")
if (!require("tidyr")) install.packages("tidyr")

# Load necessary libraries
library(stringr)
library(dplyr)
library(readr)
library(httr)
library(jsonlite)
library(VennDiagram)
library(grid)
library(tidyr)

#### Set Working Directory ####
# Modify the path as necessary
setwd("D:/Research/CXCR7/CXCR7_platelet_analysis/CXCR7_validation/phosphoproteomics/scripts/R")

#### Load Local Data ####
# Load prepared data from previous scripts (e.g., pl.phospho.10.0.network.prepare.tables.R)
# This step assumes that the following data frames are already loaded into the environment:
# - top.all.init.cxcr7.vs.0s
# - top.all.val.cxcr7.vs.0s
# - top.all.val.dmso.vs.0s
# - top.all.val.dmso.vs.cxcr7
# - top.all.prot.val
# - top.all.prot.init

#### Customize Nodes for Network ####
node_proteom1 <- unique(rownames(top.all.prot.val))
node_proteom2 <- unique(rownames(top.all.prot.init))
node_psite <- unique(top.all.init.cxcr7.vs.0s$uniprot)
node_psite2 <- unique(top.all.val.cxcr7.vs.0s$uniprot)

uniprot_id <- unique(c(node_proteom1, node_proteom2, node_psite, node_psite2))

# Test if specific proteins are present in the combined list
proteins_to_test <- c("P25106", "P61073", "P49407", "P32121", "P62993", "P42345", "Q03135", "P51636", "Q7Z460", "O75122", "Q15942", "Q05209")
gene_symbols <- c("CXCR7", "CXCR4", "ARRB1", "ARRB2", "GRB2", "MTOR", "CAV1", "CAV2", "CLASP1", "CLASP2", "ZYX", "PTN12")
proteins_df <- data.frame(Protein_ID = proteins_to_test, Gene_Symbol = gene_symbols)
proteins_df$In_List <- proteins_to_test %in% uniprot_id
print(proteins_df)

#### Download and Process IntAct Data ####

# Function to download datasets from the IntAct FTP server
download_intact_data <- function(url, destfile, metadata_file) {
  cat("Downloading data from IntAct...\n")
  httr::GET(url, httr::write_disk(destfile, overwrite = TRUE))
  cat("Download completed: ", destfile, "\n")
  
  # Create metadata
  metadata <- list(
    url = url,
    destfile = destfile,
    download_date = Sys.time()
  )
  
  # Save metadata to JSON file
  write_json(metadata, metadata_file, pretty = TRUE)
  cat("Metadata saved to: ", metadata_file, "\n")
}

# Define the IntAct FTP URLs and the destination file paths
intact_urls <- list(
  intact = "https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.txt",
  intact_negative = "https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact_negative.txt"
)

destfiles <- list(
  intact = "../../../phosphoproteomics/data/raw_data/intact.txt",
  intact_negative = "../../../phosphoproteomics/data/raw_data/intact_negative.txt"
)

metadata_files <- list(
  intact = "../../../phosphoproteomics/data/raw_data/intact_metadata.json",
  intact_negative = "../../../phosphoproteomics/data/raw_data/intact_negative_metadata.json"
)

# Download the datasets if not already downloaded
for (dataset in names(intact_urls)) {
  if (!file.exists(destfiles[[dataset]])) {
    download_intact_data(intact_urls[[dataset]], destfiles[[dataset]], metadata_files[[dataset]])
  } else {
    cat("File already exists: ", destfiles[[dataset]], "\n")
  }
}

#### Load and Process Data ####
# Load IntAct data in chunks, applying filters and selecting columns
process_chunk_1 <- function(data, pos) {
  data %>%
    filter(grepl("taxid:9606", `Host organism(s)`)) %>%  # Filter for human interactions
    select(1, 2, 5:13, 15, 21, 22)                       # Select relevant columns
}

# Read and process the 'intact.txt' dataset in chunks
interactions.intact <- read_delim_chunked(
  destfiles$intact, 
  delim = "\t", 
  callback = DataFrameCallback$new(process_chunk_1), 
  chunk_size = 100000
)

# Additional processing for other data chunks can be added here

#### Filter and Refine Interactions ####
interactions.all <- interactions.intact

# Rename columns for clarity
colnames(interactions.all) <- c("A", "B", "aliasA", "aliasB", "detectionMethod", "firstAuthor", "publicationID", "taxonA", "taxonB", "type", "sourceDatabases", "confidenceScore", "typeA", "typeB")

# Filter for relevant interactions
detection_methods <- c("psi-mi:MI:0407(direct interaction)", "psi-mi:MI:0915(physical association)")
interactions.filtered <- interactions.all %>%
  filter(type %in% detection_methods) %>%
  filter(grepl("uniprotkb", A) & grepl("uniprotkb", B)) %>%
  mutate(A = gsub("uniprotkb:", "", A),
         B = gsub("uniprotkb:", "", B)) %>%
  select(A, B, aliasA, aliasB, confidenceScore, type)

# Save refined interaction data
write.table(interactions.filtered, "../../data/processed_data/interactions.intact.all.refined.txt", sep="\t", row.names=FALSE, quote = FALSE, eol = "\n")





##### export to cytoscpae
## network sif file
interactions.intact.all.refined.platelets.b <- interactions.intact.all.refined.platelets.collapsed.omni
sif <- data.frame(interactions.intact.all.refined.platelets.b$uniprot_id1, "interacts_with", interactions.intact.all.refined.platelets.b$uniprot_id2)
colnames(sif) <- c("uniprot_id1", "type", "uniprot_id2")
write.table(sif, "../../analysis/intact/sif.txt", sep="\t", row.names=FALSE, quote = FALSE, eol = "\n")



### edge info
edge <- data.frame(paste0(interactions.intact.all.refined.platelets.b$uniprot_id1, " (", "interacts_with", ") ", interactions.intact.all.refined.platelets.b$uniprot_id2),
interactions.intact.all.refined.platelets.b$gene_name1, interactions.intact.all.refined.platelets.b$gene_name2,
interactions.intact.all.refined.platelets.b$curation_effort, interactions.intact.all.refined.platelets.b$type,
interactions.intact.all.refined.platelets.b$score)
colnames(edge) <- c("interaction_id", "gene_name1", "gene_name2", "curation_effort", "type", "score")
write.table(edge, "../../analysis/intact/edge.txt", sep="\t", row.names=FALSE, quote = FALSE)



## node information file
# Function to collapse the dataset by selecting the site with the maximum logFC
collapse_dataset <- function(dataset, logFC_columns, adj_p_val_columns, fixed_columns, columns_to_remove) {
  # Remove max_logFC if it already exists to avoid conflicts
  if ("max_logFC" %in% colnames(dataset)) {
    dataset <- dataset %>%
      dplyr::select(-max_logFC)
  }
  
  # Calculate max_logFC
  dataset <- dataset %>%
    mutate(max_logFC = do.call(pmax, c(.[logFC_columns], list(na.rm = TRUE))))
  
  # Collapse the dataset
  collapsed_df <- dataset %>%
    group_by(uniprot, symbol) %>%
    filter(max_logFC == max(max_logFC, na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::select(-max_logFC)
  
  # Reorder columns and remove specified columns
  collapsed_df <- collapsed_df %>%
    dplyr::select(-one_of(columns_to_remove))

  other_columns <- setdiff(colnames(collapsed_df), fixed_columns)
  sorted_columns <- sort(other_columns)
  final_column_order <- c(fixed_columns, sorted_columns)
  collapsed_df <- collapsed_df[, final_column_order]
    
  return(collapsed_df)
}

# Function to extend the collapsed dataset with missing uniprot_ids
# Function to extend the collapsed dataset with missing uniprot_ids
extend_with_missing_uniprot_ids <- function(collapsed_df, interactions_df, logFC_columns, adj_p_val_columns, p_value_columns, expr_columns) {
    part1 <- interactions_df %>%
        dplyr::select(uniprot_id1, gene_name1)
    colnames(part1) <- c("uniprot", "symbol")
    part2 <- interactions_df %>%
        dplyr::select(uniprot_id2, gene_name2)
    colnames(part2) <- c("uniprot", "symbol")   
    unique_uniprot_symbols <- bind_rows(part1, part2) %>%
        distinct()  
    
    missing_uniprot_ids <- setdiff(unique_uniprot_symbols$uniprot, collapsed_df$uniprot)
    new_rows <- anti_join(unique_uniprot_symbols, collapsed_df, by = "uniprot")
    collapsed_df_extended <- bind_rows(collapsed_df, new_rows)

    # Fill missing values
    collapsed_df_extended <- collapsed_df_extended %>%
      mutate(across(all_of(logFC_columns), ~ ifelse(is.na(.), 0, .))) %>%
      mutate(across(all_of(adj_p_val_columns), ~ ifelse(is.na(.), 1, .))) %>%
      mutate(across(all_of(p_value_columns), ~ ifelse(is.na(.), 1, .))) %>%
      mutate(across(all_of(expr_columns), ~ ifelse(is.na(.), 0, .)))

      return(collapsed_df_extended)
}




#### Collapsing and Extending Datasets ####

# Example usage of collapsing function for DMSO vs CXCR7
collapsed_df_dmso_vs_cxcr7 <- collapse_dataset(
  dataset = top.all.val.dmso.vs.cxcr7,
  logFC_columns = c("logFC_10_dmso.vs.cxcr7", "logFC_600_dmso.vs.cxcr7", "logFC_1800_dmso.vs.cxcr7"),
  adj_p_val_columns = c("adj.P.Val_10_dmso.vs.cxcr7", "adj.P.Val_600_dmso.vs.cxcr7", "adj.P.Val_1800_dmso.vs.cxcr7"),
  fixed_columns = c("uniprot", "symbol", "psite", "id", 
                    "logFC_10_dmso.vs.cxcr7", "logFC_600_dmso.vs.cxcr7", "logFC_1800_dmso.vs.cxcr7",
                    "adj.P.Val_10_dmso.vs.cxcr7", "adj.P.Val_600_dmso.vs.cxcr7", "adj.P.Val_1800_dmso.vs.cxcr7"),
  columns_to_remove = c("t_10_dmso.vs.cxcr7", "t_1800_dmso.vs.cxcr7", "t_600_dmso.vs.cxcr7", 
                        "B_10_dmso.vs.cxcr7", "B_1800_dmso.vs.cxcr7", "B_600_dmso.vs.cxcr7")
)

# Example usage of collapsing function for DMSO vs 0s
collapsed_df_dmso_vs_0s <- collapse_dataset(
  dataset = top.all.val.dmso.vs.0s,
  logFC_columns = c("logFC_10_dmso.vs.0s", "logFC_600_dmso.vs.0s", "logFC_1800_dmso.vs.0s"),
  adj_p_val_columns = c("adj.P.Val_10_dmso.vs.0s", "adj.P.Val_600_dmso.vs.0s", "adj.P.Val_1800_dmso.vs.0s"),
  fixed_columns = c("uniprot", "symbol", "psite", "id", 
                    "logFC_10_dmso.vs.0s", "logFC_600_dmso.vs.0s", "logFC_1800_dmso.vs.0s",
                    "adj.P.Val_10_dmso.vs.0s", "adj.P.Val_600_dmso.vs.0s", "adj.P.Val_1800_dmso.vs.0s"),
  columns_to_remove = c("t_10_dmso.vs.0s", "t_1800_dmso.vs.0s", "t_600_dmso.vs.0s", 
                        "B_10_dmso.vs.0s", "B_1800_dmso.vs.0s", "B_600_dmso.vs.0s")
)

# Extend collapsed data frames with missing uniprot_ids
final_collapsed_df_dmso_vs_cxcr7 <- extend_with_missing_uniprot_ids(
  collapsed_df_dmso_vs_cxcr7, 
  interactions.intact.all.refined.platelets.b, 
  logFC_columns = c("logFC_10_dmso.vs.cxcr7", "logFC_600_dmso.vs.cxcr7", "logFC_1800_dmso.vs.cxcr7"),
  adj_p_val_columns = c("adj.P.Val_10_dmso.vs.cxcr7", "adj.P.Val_600_dmso.vs.cxcr7", "adj.P.Val_1800_dmso.vs.cxcr7"),
  p_value_columns = c("P.Value_10_dmso.vs.cxcr7", "P.Value_1800_dmso.vs.cxcr7", "P.Value_600_dmso.vs.cxcr7"),
  expr_columns = c("AveExpr_10_dmso.vs.cxcr7", "AveExpr_1800_dmso.vs.cxcr7", "AveExpr_600_dmso.vs.cxcr7")
)

final_collapsed_df_dmso_vs_0s <- extend_with_missing_uniprot_ids(
  collapsed_df_dmso_vs_0s, 
  interactions.intact.all.refined.platelets.b, 
  logFC_columns = c("logFC_10_dmso.vs.0s", "logFC_600_dmso.vs.0s", "logFC_1800_dmso.vs.0s"),
  adj_p_val_columns = c("adj.P.Val_10_dmso.vs.0s", "adj.P.Val_600_dmso.vs.0s", "adj.P.Val_1800_dmso.vs.0s"),
  p_value_columns = c("P.Value_10_dmso.vs.0s", "P.Value_1800_dmso.vs.0s", "P.Value_600_dmso.vs.0s"),
  expr_columns = c("AveExpr_10_dmso.vs.0s", "AveExpr_1800_dmso.vs.0s", "AveExpr_600_dmso.vs.0s")
)

# Write the final data frames to files
write.table(final_collapsed_df_dmso_vs_cxcr7, "../../analysis/intact/node_dmso_vs_cxcr7.txt", sep="\t", row.names=FALSE, quote = FALSE)
write.table(final_collapsed_df_dmso_vs_0s, "../../analysis/intact/node_dmso_vs_0s.txt", sep="\t", row.names=FALSE, quote = FALSE)
