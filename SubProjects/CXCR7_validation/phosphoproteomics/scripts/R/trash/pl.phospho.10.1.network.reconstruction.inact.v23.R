################################################################
BiocManager::install("PSICQUIC")
if (!requireNamespace("stringr", quietly = TRUE)) {
  install.packages("stringr")



##### intact network reconstruction 
#library("PSICQUIC")
#library(config)       # For configuration management
#library(SparkR)
library(stringr)
library(readr)        # For efficient reading of large files
library(dplyr)        # For data manipulation
library(httr)         # For downloading files from the internet
library(jsonlite)    # For creating and handling JSON metadata files
library(VennDiagram)
library(grid)
library(tidyr)



##### set working directory
# e.g. 
setwd("C:/Users/Ozge/Nextcloud/Phosphoproteom/phosphoproteom validation/scripts") ##ÖO folder structure should fit.
setwd("D:/Research/CXCR7_platelet_analysis/CXCR7_validation/phosphoproteomics/scripts/R")
setwd("D:/Research/CXCR7/CXCR7_platelet_analysis/CXCR7_validation/phosphoproteomics/scripts/R")


##### load local data
# use pl.phospho.10.0.network.prepare.tables.R
# load proteom data
#top.all.init.cxcr7.vs.0s
#top.all.val.cxcr7.vs.0s
#top.all.val.dmso.vs.0s
#top.all.val.dmso.vs.cxcr7
#top.all.prot.val
#top.all.prot.init


#### customize the nodes (uniprot_id, proteins) for your network
node_proteom1 <- unique(c(rownames(top.all.prot.val)))
node_proteom2 <- unique(c(rownames(top.all.prot.init)))
node_psite <- unique(top.all.init.cxcr7.vs.0s$uniprot)
node_psite2 <- unique(top.all.val.cxcr7.vs.0s$uniprot)

intersection_proteom <- intersect(node_proteom1, node_proteom2)
intersection_psite <- intersect(node_psite, node_psite2)

uniprot_id  <-unique(c(node_proteom1, node_proteom2, node_psite, node_psite2))

#P25106 cxcr7
#P61073 cxcr4
#P49407 ARRB1
#P32121 ARRB2
#"P49407" ARRB
#"P62993" GRB2
#"P42345" MTOR
#"Q03135" CAV1
#"P51636" CAV2
#"Q7Z460" CLASP1
#"O75122" CLASP2
#"Q15942" ZYX
#"Q05209" PTN12

## test if proteins of interest are within the diffrent lists (node_proteom1, node_proteom2, node_psite, node_psite2)
proteins_to_test <- c("P25106", "P61073", "P49407", "P32121", "P49407", "P62993", "P42345", "Q03135", "P51636", "Q7Z460", "O75122", "Q15942", "Q05209")
gene_symbols <- c("CXCR7", "CXCR4", "ARRB1", "ARRB2", "ARRB", "GRB2", "MTOR", "CAV1", "CAV2", "CLASP1", "CLASP2", "ZYX", "PTN12")
proteins_df <- data.frame(Protein_ID = proteins_to_test, Gene_Symbol = gene_symbols)
proteins_in_list <- proteins_to_test %in% uniprot_id
proteins_df$In_List <- proteins_in_list
print(proteins_df)



## download from intact (note: dataset is big, dont use when u have it already)
# https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/

## memory limit 
# https://www.rdocumentation.org/packages/utils/versions/3.4.1/topics/memory.size
# memory.size(max = FALSE)
# memory.limit(size = NA)
#  gc()
# read with sparkR
# install.packages("SparkR")

# https://docs.databricks.com/spark/latest/sparkr/overview.html
# install.packages("config")
# https://cran.r-project.org/web/packages/config/vignettes/introduction.html


# Function to download datasets from the IntAct FTP server and log metadata
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
  
  # Write metadata to JSON file
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


# Download the datasets if not already downloaded and save metadata
for (dataset in names(intact_urls)) {
  if (!file.exists(destfiles[[dataset]])) {
    download_intact_data(intact_urls[[dataset]], destfiles[[dataset]], metadata_files[[dataset]])
  } else {
    cat("File already exists: ", destfiles[[dataset]], "\n")
  }
}



## Load data into R
# Function for custom chunk processing
process_chunk_1 <- function(data, pos) {
  data %>%
    filter(grepl("taxid:9606", `Host organism(s)`)) %>%  # Filter for human interactions
    select(1, 2, 5:13, 15, 21, 22)                       # Select relevant columns
}

process_chunk_2 <- function(data, pos) {
  data[, c(1, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 21, 22)]  # Select specific columns
}

process_chunk_3 <- function(data, pos) {
  subset(data, grepl("9606", data$`Host organism(s)`))   # Filter using grepl
}

# Read and process the 'intact.txt' dataset in chunks
cat("Reading and processing the IntAct dataset in chunks...\n")
interactions.intact <- read_delim_chunked(
  destfiles$intact, 
  delim = "\t", 
  callback = DataFrameCallback$new(process_chunk_1), 
  chunk_size = 100000
)
cat("Data processing for intact.txt completed.\n")

# Read and process the 'intact_negative.txt' dataset in chunks using a different filter
cat("Reading and processing the IntAct Negative dataset in chunks...\n")
interactions.intact_negative <- read_delim_chunked(
  destfiles$intact_negative, 
  delim = "\t", 
  callback = DataFrameCallback$new(process_chunk_2), 
  chunk_size = 1000000
)
cat("Data processing for intact_negative.txt completed.\n")


# Display the first few rows of the processed data
print(head(interactions.intact))
print(head(interactions.intact_negative))


##### parse and pre-filter the imported interactions
interactions.all <- interactions.intact

colnames(interactions.all)
colnames(interactions.all) <- c("A","B","aliasA","aliasB","detectionMethod","firstAuthor","publicationID","taxonA","taxonB","type","sourceDatabases","confidenceScore","typeA","typeB")

#table(interactions.all$provider)
table(interactions.all$detectionMethod)
table(interactions.all$type)

detectionmethod <- c("psi-mi:MI:0407(direct interaction)", "psi-mi:MI:0915(physical association)", 
"psi-mi:MI:0407(direct interaction)")
interactions.all.a <- interactions.all
interactions.all.a$type <- gsub("\"","",interactions.all.a$type)
interactions.all.a <- interactions.all.a[which(interactions.all.a$type %in% detectionmethod),]
interactions.all.a <- interactions.all.a[grep("uniprotkb", interactions.all.a$A),]
interactions.all.a <- interactions.all.a[grep("uniprotkb", interactions.all.a$B),]
interactions.all.a$A <- gsub("uniprotkb:", "", interactions.all.a$A)
interactions.all.a$B <- gsub("uniprotkb:", "", interactions.all.a$B)

interactions.all.a$gene_name1 <- str_extract(interactions.all.a$aliasA, "\\|uniprotkb:\\S+(gene name)")
interactions.all.a$gene_name1 <- gsub("\\|uniprotkb:(.+)\\(gene name","\\1",interactions.all.a$gene_name1)

interactions.all.a$gene_name2 <- str_extract(interactions.all.a$aliasB, "\\|uniprotkb:\\S+(gene name)")
interactions.all.a$gene_name2 <- gsub("\\|uniprotkb:(.+)\\(gene name","\\1",interactions.all.a$gene_name2)

interactions.all.a$score <- gsub(".*intact-miscore:(\\d+)","\\1",interactions.all.a$confidenceScore,,perl=TRUE)

## check how tables look like
colnames(interactions.all.a)
# [1] "A"               "B"               "aliasA"          "aliasB"         
# [5] "detectionMethod" "firstAuthor"     "publicationID"   "taxonA"         
# [9] "taxonB"          "type"            "sourceDatabases" "confidenceScore"
#[13] "typeA"           "typeB"           "gene_name1"      "gene_name2"     
#[17] "score" 

interactions.all.b <- interactions.all.a[,c(1,2,15,16,17,10,5,6,7,11)]
#interactions.all.b <- interactions.all.a[,c(interactions.all.a$A,interactions.all.a$B,interactions.all.a$gene_name1,
#interactions.all.a$gene_name2,interactions.all.a$score,interactions.all.a$type,interactions.all.a$detectionMethod,
#interactions.all.a$firstAuthor,interactions.all.a$publicationID,interactions.all.a$sourceDatabase)]

colnames(interactions.all.b)[1] = "uniprot_id1"
colnames(interactions.all.b)[2] = "uniprot_id2"

colnames(interactions.all.b)
# [1] "uniprot_id1"     "uniprot_id2"     "gene_name1"      "gene_name2"     
# [5] "score"           "type"            "detectionMethod" "firstAuthor"    
# [9] "publicationID"   "sourceDatabases"


## save interactions as file
write.table(interactions.all.b, "../../data/processed_data/interactions.intact.all.refined.txt", sep="\t", row.names=FALSE, , quote = FALSE, , eol = "\n")
## if interactions.all.b.txt already exists just go on to read it drectly
interactions.all.b <- read.delim("../../data/processed_data//interactions.intact.all.refined.txt", sep="\t", header=TRUE, skip = 0)


## filter interactions
interactions.all.c <- interactions.all.b[interactions.all.b$score > 0.3,]
dim(interactions.all.c)

##### now filter for the PL proteins
## test whether proteins of interest are present in data, and if not add manually
#P25106 cxcr7
#P61073 cxcr4
#P49407 ARRB1
#P32121 ARRB2

interactions.all.c[interactions.all.c$uniprot_id1 == "P25106",]
UniprotID <- c(uniprot_id,"P25106","P61073")


## filter the interactions for the nodes of interest 
interactions.all.d <- unique(interactions.all.c[which(interactions.all.c$uniprot_id1 %in% UniprotID),])
interactions.all.d <- unique(interactions.all.d[which(interactions.all.d$uniprot_id2 %in% UniprotID),])
dim(interactions.all.d)


## save interactions as file
write.table(interactions.all.d, "../../data/processed_data/interactions.intact.all.refined.platelets.txt", sep="\t", row.names=FALSE, , quote = FALSE, , eol = "\n")
## if interactions.all.b.txt already exists just go on to read it drectly
interactions.intact.all.refined.platelets <- read.delim("../../data/processed_data/interactions.intact.all.refined.platelets.txt", sep="\t", header=TRUE, skip = 0)


# Create the reverse interactions
reverse_interactions <- interactions.intact.all.refined.platelets.b %>%
  mutate(temp_uniprot = uniprot_id1,
         temp_gene = gene_name1) %>%
  mutate(uniprot_id1 = uniprot_id2,
         gene_name1 = gene_name2,
         uniprot_id2 = temp_uniprot,
         gene_name2 = temp_gene) %>%
  dplyr::select(-temp_uniprot, -temp_gene)

# Combine the original and reverse interactions
combined_interactions <- bind_rows(interactions.intact.all.refined.platelets.b, reverse_interactions)

# Remove duplicate rows to ensure the table is undirected and unique
unique_interactions <- combined_interactions %>%
  distinct()

# View the final unique interactions
View(unique_interactions)


# Collapse by the specified columns
collapsed_interactions <- unique_interactions %>%
  group_by(uniprot_id1, uniprot_id2, gene_name1, gene_name2) %>%
  summarise(
    score = max(score),
    interaction_type = paste(unique(type), collapse = ";"),
    detectionMethod = paste(unique(detectionMethod), collapse = ";"),
    firstAuthor = paste(unique(firstAuthor), collapse = ";"),
    publicationID = paste(unique(publicationID), collapse = ";"),
    sourceDatabases = paste(unique(sourceDatabases), collapse = ";"),
    .groups = 'drop'
  )

# View the final collapsed interactions
collapsed_interactions

interactions.intact.all.refined.platelets.collapsed <- collapsed_interactions[collapsed_interactions$score > 0.4, ]

write.table(interactions.intact.all.refined.platelets.collapsed, "../../data/processed_data/interactions.intact.all.refined.platelets.collpased.txt", sep="\t", row.names=FALSE, , quote = FALSE, , eol = "\n")



# Perform the left join: uniprot_id1 = source and uniprot_id2 = target
final_df <- interactions.intact.all.refined.platelets.collapsed %>%
  left_join(filtered_interactions_omni2 %>%
              dplyr::select(source, target, type, curation_effort),
            by = c("uniprot_id1" = "source", "uniprot_id2" = "target"))

# fill NAs
final_df <- final_df %>%
  mutate(
    curation_effort = replace_na(curation_effort, 0),
    type = replace_na(type, "unknown")
  )

# View the final dataframe
head(final_df)
# View the final dataframe
colnames(final_df)
final_df


interactions.intact.all.refined.platelets.collapsed.omni <- final_df

write.table(interactions.intact.all.refined.platelets.collapsed.omni, "../../data/processed_data/interactions.intact.all.refined.platelets.collpased.omni.txt", sep="\t", row.names=FALSE, quote = FALSE, eol = "\n")

interactions.intact.all.refined.platelets.collapsed.omni <- read.table(
  "../../data/processed_data/interactions.intact.all.refined.platelets.collpased.omni.txt", 
  sep = "\t",              
  header = TRUE,           
  stringsAsFactors = FALSE,
  fill = TRUE
)



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
#####  get top tables data and abundance data
#### collapse sites
### 1. (prefered) collpase top.all.val.dmso.vs.cxcr7 from psite to uniprot by the site that has the top adj.P.val in "adj.P.Val_10_dmso.vs.cxcr7", "adj.P.Val_600_dmso.vs.cxcr7", "adj.P.Val_1800_dmso.vs.cxcr7"
# Function to collapse the dataset by selecting the site with the maximum logFC
collapse_dataset <- function(dataset, logFC_columns, adj_p_val_columns, fixed_columns, columns_to_remove) {
  dataset <- dataset %>%
    mutate(max_logFC = do.call(pmax, c(.[logFC_columns], na.rm = TRUE)))
  
  collapsed_df <- dataset %>%
    group_by(uniprot, symbol) %>%
    filter(max_logFC == max(max_logFC, na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::select(-max_logFC)
  
  other_columns <- setdiff(colnames(collapsed_df), fixed_columns)
  sorted_columns <- sort(other_columns)
  final_column_order <- c(fixed_columns, sorted_columns)
  collapsed_df <- collapsed_df[, final_column_order]
  
  collapsed_df <- collapsed_df %>%
    dplyr::select(-one_of(columns_to_remove))
  
  return(collapsed_df)
}

# Example usage of collapsing function
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


### 2. collpase by uniprot, choose row with max(abs(logFC))

library(dplyr)

# Define the logFC columns and corresponding adj.P.Val and P.Value columns
logFC_columns <- c("logFC_10_dmso.vs.cxcr7", "logFC_600_dmso.vs.cxcr7", "logFC_1800_dmso.vs.cxcr7")
adj_pval_columns <- c("adj.P.Val_10_dmso.vs.cxcr7", "adj.P.Val_600_dmso.vs.cxcr7", "adj.P.Val_1800_dmso.vs.cxcr7")
pval_columns <- c("P.Value_10_dmso.vs.cxcr7", "P.Value_600_dmso.vs.cxcr7", "P.Value_1800_dmso.vs.cxcr7")

# Create a list to store the collapsed data frames
collapsed_list <- list()

# Iterate over the logFC columns and collapse the data frame for each
for (i in seq_along(logFC_columns)) {
  logFC_column <- logFC_columns[i]
  adj_pval_column <- adj_pval_columns[i]
  pval_column <- pval_columns[i]
  
  # Collapse the data frame by taking the max absolute logFC for each uniprot
  collapsed_df <- top.all.val.dmso.vs.cxcr7 %>%
    group_by(uniprot) %>%
    filter(abs(.data[[logFC_column]]) == max(abs(.data[[logFC_column]]))) %>%
    ungroup() %>%
    dplyr::select(uniprot, symbol, id, 
           psite, 
           !!sym(logFC_column), 
           !!sym(adj_pval_column), 
           !!sym(pval_column))
  
  # Rename the psite column specifically for this time point
  psite_column_name <- paste0("psite_", sub("logFC_", "", logFC_column))
  collapsed_df <- collapsed_df %>%
    rename(!!psite_column_name := psite)
  
  # Store the collapsed data frame in the list
  collapsed_list[[i]] <- collapsed_df
}

# Merge all collapsed data frames by "uniprot", "symbol", and "id"
top.all.logfc.collapse <- Reduce(function(x, y) merge(x, y, by = c("uniprot", "symbol"), all = TRUE), collapsed_list)

# View the final collapsed data frame
top.all.logfc.collapse


# Rename the columns
top.all.logfc.collapse <- top.all.logfc.collapse %>%
  rename(id_10 = id.x,
         id_600 = id.y,
         id_1800 = id)

# Check the updated column names
colnames(top.all.logfc.collapse)

write.table(top.all.logfc.collapse, "../../analysis/intact/node.txt", sep="\t", row.names=FALSE, quote = FALSE)




### 3. with phoscollpase with seem to be faily
top.all.logfc <- top.all.val.dmso.vs.cxcr7[,c(3,6,9)]
library(PhosR)
top.all.logfc.collapse.b <- phosCollapse(top.all.logfc, id=sapply(strsplit(rownames(norm_intensity), ";"), "[[", 1), 
                        stat=apply(abs(top.all.logfc), 1, max), by = "max")

top.all.pval <- top.all[,c(5,8,11)]
library(PhosR)
top.all.pval.collapse <- phosCollapse(top.all.pval, id=sapply(strsplit(rownames(norm_intensity), ";"), "[[", 1), 
                        stat=apply(abs(top.all.pval), 1, min), by = "min")



### 4. or collpase self made
top.all.logfc.collapse <- ddply(top.all.val.dmso.vs.cxcr7, .(uniprot), summarise,
              symbol = max(symbol),
              min.logFC.10 = min(logFC.10),
              max.logFC.10 = max(logFC.10),
              Abs.logFC.10 = max(abs(logFC.10)),
              min.logFC.600 = min(logFC.600),
              max.logFC.600 = max(logFC.600),
              Abs.logFC.600 = max(abs(logFC.600)),
              min.logFC.1800 = min(logFC.1800),
              max.logFC.1800 = max(logFC.1800),
              Abs.logFC.1800 = max(abs(logFC.1800))
             )

#top.collapse[abs(top.collapse[, "min.logFC.10"]) < abs(top.collapse[, "max.logFC.10"]),]

## 10 sec
top.all.logfc.collapse$logFC.10 <- 0
top.all.logfc.collapse$logFC.10[abs(top.all.logfc.collapse[, "min.logFC.10"]) < abs(top.all.logfc.collapse[, "max.logFC.10"])] = top.all.logfc.collapse$max.logFC.10[abs(top.all.logfc.collapse[, "min.logFC.10"]) < abs(top.all.logfc.collapse[, "max.logFC.10"])]

top.all.logfc.collapse$logFC.10[abs(top.all.logfc.collapse[, "min.logFC.10"]) > abs(top.all.logfc.collapse[, "max.logFC.10"])] = top.all.logfc.collapse$min.logFC.10[abs(top.all.logfc.collapse[, "min.logFC.10"]) > abs(top.all.logfc.collapse[, "max.logFC.10"])]

top.all.logfc.collapse$logFC.10[abs(top.all.logfc.collapse[, "min.logFC.10"]) == abs(top.all.logfc.collapse[, "max.logFC.10"])] = top.all.logfc.collapse$min.logFC.10[abs(top.all.logfc.collapse[, "min.logFC.10"]) == abs(top.all.logfc.collapse[, "max.logFC.10"])]
## 600 sec
top.all.logfc.collapse$logFC.600 <- 0
top.all.logfc.collapse$logFC.600[abs(top.all.logfc.collapse[, "min.logFC.600"]) < abs(top.all.logfc.collapse[, "max.logFC.600"])] = top.all.logfc.collapse$max.logFC.600[abs(top.all.logfc.collapse[, "min.logFC.600"]) < abs(top.all.logfc.collapse[, "max.logFC.600"])]

top.all.logfc.collapse$logFC.600[abs(top.all.logfc.collapse[, "min.logFC.600"]) > abs(top.all.logfc.collapse[, "max.logFC.600"])] = top.all.logfc.collapse$min.logFC.600[abs(top.all.logfc.collapse[, "min.logFC.600"]) > abs(top.all.logfc.collapse[, "max.logFC.600"])]

top.all.logfc.collapse$logFC.600[abs(top.all.logfc.collapse[, "min.logFC.600"]) == abs(top.all.logfc.collapse[, "max.logFC.600"])] = top.all.logfc.collapse$min.logFC.600[abs(top.all.logfc.collapse[, "min.logFC.600"]) == abs(top.all.logfc.collapse[, "max.logFC.600"])]
## 1800 sec
top.all.logfc.collapse$logFC.1800 <- 0
top.all.logfc.collapse$logFC.1800[abs(top.all.logfc.collapse[, "min.logFC.1800"]) < abs(top.all.logfc.collapse[, "max.logFC.1800"])] = top.all.logfc.collapse$max.logFC.1800[abs(top.all.logfc.collapse[, "min.logFC.1800"]) < abs(top.all.logfc.collapse[, "max.logFC.1800"])]

top.all.logfc.collapse$logFC.1800[abs(top.all.logfc.collapse[, "min.logFC.1800"]) > abs(top.all.logfc.collapse[, "max.logFC.1800"])] = top.all.logfc.collapse$min.logFC.1800[abs(top.all.logfc.collapse[, "min.logFC.1800"]) > abs(top.all.logfc.collapse[, "max.logFC.1800"])]

top.all.logfc.collapse$logFC.1800[abs(top.all.logfc.collapse[, "min.logFC.1800"]) == abs(top.all.logfc.collapse[, "max.logFC.1800"])] = top.all.logfc.collapse$min.logFC.1800[abs(top.all.logfc.collapse[, "min.logFC.1800"]) == abs(top.all.logfc.collapse[, "max.logFC.1800"])]

top.all.logfc.collapse.b <- as.data.frame(top.all.logfc.collapse[,c("logFC.10", "logFC.600", "logFC.1800")])
rownames(top.all.logfc.collapse.b) <- top.all.logfc.collapse$uniprot
top.all.logfc.collapse.b$uniprot <- top.all.logfc.collapse$uniprot
top.all.logfc.collapse.b$symbol <- top.all.logfc.collapse$symbol


interactors <- unique(c(interactions.all.d$uniprot_id1,interactions.all.d$uniprot_id2))

node <- as.data.frame(unique(top.all.logfc.collapse.b[which(rownames(top.all.logfc.collapse.b) %in% interactors),]))
node$uniprot <- rownames(node)
node1 <- merge(node, top.all[,c(1,3,5)], by = c("uniprot", "logFC.10"))
node2 <- merge(node1, top.all[,c(1,6,8)], by = c("uniprot", "logFC.600"))
node3 <- merge(node2, top.all[,c(1,9,11)], by = c("uniprot", "logFC.1800"))
node3 <- node3[,c(1,5,4,3,2,6,7,8)]
rownames(node3) <- node3$uniprot

## to get gene_name on top of uniprot sybols
library(dplyr)
df1 <- as.data.frame(interactions.all.d[,c(1,3)])
df2 <- as.data.frame(interactions.all.d[,c(2,4)])
colnames(df1) <- c("uniprot","gene_name")
colnames(df2) <- c("uniprot","gene_name")
int_names <- unique(rbind(df1,df2))
rownames(int_names) <- int_names$uniprot

node4 <- merge(node3,int_names,by="row.names")
node4 <- node4[,-10]




## export
write.table(node4, "../../analysis/cytsocape/node4.txt", sep="\t", row.names=FALSE, quote = FALSE)

