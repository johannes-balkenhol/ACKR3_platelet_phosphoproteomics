
#### packages
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
install.packages("RCy3")
BiocManager::install("KEGGREST")


library(igraph)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)  # Assuming you're working with human data
library(ReactomePA)
library(KEGGREST)
library(fgsea)
library(MCL)
#library(RCy3)



#### source data
#interactions.platelet.combined
#top.all.val.dmso.vs.cxcr7


#### convert to graph
# Extract the edge list
#edges <- interactions.platelet.combined[, c("source", "target")]
edges <- interactions.intact.all.refined.platelets[, c("uniprot_id1", "uniprot_id2")]
interactions.intact.all.refined.platelets.b <- interactions.intact.all.refined.platelets[interactions.intact.all.refined.platelets$score > 0.6, ]
original_pairs <- interactions.intact.all.refined.platelets.b[, c("uniprot_id1", "uniprot_id2")]
reversed_pairs <- original_pairs[, c("uniprot_id2", "uniprot_id1")]
colnames(reversed_pairs) <- c("uniprot_id1", "uniprot_id2")
combined_pairs <- unique(rbind(original_pairs, reversed_pairs))
dim(combined_pairs) 
# 2677


#### test if key proteins are represented in the network 
network_nodes <- unique(c(interactions.platelet.combined$source, interactions.platelet.combined$target))
network_nodes <- unique(c(filtered_interactions_omni2$source, filtered_interactions_omni2$target))
network_nodes <- unique(c(final_merged$source, final_merged$target))
network_nodes <- unique(c(interactions.intact.all.refined.platelets.b$uniprot_id1, interactions.intact.all.refined.platelets.b$uniprot_id2))
length(network_nodes) 
# 1032

proteins_to_test <- c("P25106", "P61073", "P49407", "P32121", "P49407", "P62993", "P42345", "Q03135", "P51636", "Q7Z460", "O75122", "Q15942", "Q05209")
gene_symbols <- c("CXCR7", "CXCR4", "ARRB1", "ARRB2", "ARRB", "GRB2", "MTOR", "CAV1", "CAV2", "CLASP1", "CLASP2", "ZYX", "PTN12")
proteins_df <- data.frame(Protein_ID = proteins_to_test, Gene_Symbol = gene_symbols)
proteins_in_list <- proteins_to_test %in% network_nodes
proteins_df$In_List <- proteins_in_list
print(proteins_df)



# Save the edge list to a file that MCL can read
write.table(edges, "network_edges.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

graph <- graph_from_data_frame(d = edges, directed = FALSE)



#### cluster algorithms

## 1 Infomap
clusters <- cluster_infomap(graph)
membership_vector <- membership(clusters)


## 2 Walktrap
clusters <- cluster_walktrap(graph)
membership_vector <- membership(clusters)


## 3 Louvain
clusters_louvain <- cluster_louvain(graph)
membership_vector <- membership(clusters_louvain)


## 4 (MCODE)
# Assume edges is your edge list data frame
createNetworkFromDataFrames(edges = edges, title = "My Network", collection = "My Collection")
# Apply MCODE clustering
mcode_clusters <- cyrestPOST('mcode/v1/apply')
# Retrieve the cluster results
mcode_results <- cyrestGET(paste0('mcode/v1/getResults/', mcode_clusters$SUID))
# Extract node membership
mcode_cluster_df <- do.call(rbind, lapply(mcode_results$clusters, function(cluster) {
  data.frame(node = cluster$nodes, cluster = cluster$score)
}))
# Print the result
print(mcode_cluster_df)


#ä 5 MCL 
adj_matrix <- as_adjacency_matrix(graph, sparse = FALSE)
mcl_result <- mcl(adj_matrix, addLoops = TRUE, expansion = 3, inflation = 2)
clusters <- mcl_result$Cluster
# Get the cluster assignment for each node
cluster_membership <- clusters
# Add the cluster membership to your original graph
V(graph)$cluster <- cluster_membership
# View the first few nodes and their cluster assignments
head(data.frame(Node = V(graph)$name, Cluster = V(graph)$cluster))
# Create the data frame with Node and Cluster information
membership_df <- data.frame(Node = V(graph)$name, Cluster = V(graph)$cluster)
# Convert to a membership vector
membership_vector <- setNames(membership_df$Cluster, membership_df$Node)
# View the membership vector
membership_vector
# Plot the graph with nodes colored by cluster
plot(graph, vertex.color = V(graph)$cluster, vertex.label = NA, main = "MCL Clustering of PPI Network")

# Create the data frame
cluster_df <- data.frame(
  Uniprot_ID = names(membership_vector),
  Cluster_Number = as.integer(membership_vector)
)
head(cluster_df)
write.csv(cluster_df, "cluster_results.csv", row.names = FALSE)

#### investigate cluster size
cluster_size_df <- cluster_df %>%
  group_by(Cluster_Number) %>%
  summarise(Protein_Count = n())
cluster_size_df <- cluster_size_df %>%
  arrange(desc(Protein_Count))
print(cluster_size_df)



## 6 Cytoscape enrichment by hand: (visual feedback helps great by choosen parameters)
# MCL
# source script: pl.phospho-10.1.network.reconstruction.intact.v2.R
# source tables (pl.network.cluster.nodes, interactions.intact.all.refined.platelets.b, collapsed_df)
pl.network.cluster.nodes <- read.delim("../../analysis/intact/pl.network.cluster.mcl.1.7.v3.csv", sep=",", header=TRUE, skip = 0)

cluster_df <- pl.network.cluster.nodes %>%
  dplyr::select(Uniprot_ID = shared.name, Cluster_Number = X__mclCluster, Symbol = symbol)

filtered_cluster_df <- cluster_df %>%
  filter(Cluster_Number == 1)

# View the filtered data
filtered_cluster_df

#### investigate cluster size
cluster_size_df <- cluster_df %>%
  group_by(Cluster_Number) %>%
  summarise(Protein_Count = n())
cluster_size_df <- cluster_size_df %>%
  arrange(desc(Protein_Count))
print(cluster_size_df)


cluster_df <- na.omit(cluster_df)
cluster_size_df <- na.omit(cluster_size_df)



gene_sets <- split(cluster_df$Uniprot_ID, cluster_df$Cluster_Number)
uniprot_ids <- gene_sets[[1]]
num_uniprot_ids <- length(uniprot_ids)
num_uniprot_ids
gene_symbols <- bitr(uniprot_ids, fromType = "UNIPROT", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
gene_symbols




#### GO enrichment of clusters v2

background_ids <- unique(cluster_df$Uniprot_ID)

# Set a cluster size cutoff
cluster_size_cutoff <- 5
filtered_clusters <- cluster_size_df %>% 
  filter(Protein_Count >= cluster_size_cutoff)

# Prepare the data for the filtered clusters
filtered_proteins <- cluster_df %>% 
  filter(Cluster_Number %in% filtered_clusters$Cluster_Number)

# Perform GO Enrichment Analysis for Each Cluster
go_results <- list()

# Loop through each filtered cluster
for (cluster_num in filtered_clusters$Cluster_Number) {
  # Get the Uniprot IDs for the current cluster
  protein_ids <- filtered_proteins %>% 
    filter(Cluster_Number == cluster_num) %>%
    pull(Uniprot_ID)
  
  # Perform GO enrichment analysis using the "Biological Process" ontology
  ego <- enrichGO(
    gene         = protein_ids,
    OrgDb        = org.Hs.eg.db,
    keyType      = "UNIPROT",
    ont          = "BP",            
    pAdjustMethod = "BH",            
    universe     = background_ids,   
    qvalueCutoff = 0.05,             
    minGSSize    = 10,               
    maxGSSize    = 100,              
    readable     = TRUE              
  )
  
  # Get the top 3 GO terms
  top_terms <- ego@result %>% 
    arrange(p.adjust) %>% 
    head(3)
  
  # Store the results in the list
  go_results[[as.character(cluster_num)]] <- top_terms
}

# Sort clusters by highest p.adjust and print the results
sorted_clusters <- names(go_results)[order(sapply(go_results, function(x) min(x$p.adjust, na.rm = TRUE)))]

for (cluster_num in sorted_clusters) {
  cat("\nCluster:", cluster_num, "\n")
  print(go_results[[cluster_num]][, c("ID", "Description", "p.adjust")])
}




## KEGG pathway enrichment for each cluster

background_ids <- unique(top.all.val.dmso.vs.cxcr7$uniprot)
# Translate Uniprot IDs to Entrez Gene IDs
bg_entrez_mapping <- bitr(background_ids, 
                    fromType = "UNIPROT", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)

# Initialize an empty list to store results
kegg_results <- list()
search_kegg_organism('ece', by='kegg_code')
hsa <- search_kegg_organism('Homo sapiens', by='scientific_name')
dim(hsa)

# Loop through each of the top 10 clusters
for (cluster_num in top_10_clusters$Cluster_Number) {
  # Get the Uniprot IDs for the current cluster
  protein_ids <- top_10_proteins %>% 
    filter(Cluster_Number == cluster_num) %>%
    pull(Uniprot_ID)
  
  # Convert Uniprot IDs to Entrez Gene IDs
  entrez_mapping <- bitr(protein_ids, 
                         fromType = "UNIPROT", 
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)
  
  # Extract Entrez Gene IDs
  entrez_ids <- entrez_mapping$ENTREZID
  
  # Perform KEGG pathway enrichment analysis
  ekegg <- enrichKEGG(
    gene         = entrez_ids,      # Use Entrez Gene IDs
    organism     = "hsa",           # "hsa" for Homo sapiens
    keyType      = "kegg",          # Key type for KEGG (Entrez ID expected)
    pAdjustMethod = "BH",
    universe     = bg_entrez_mapping$ENTREZID,  # Convert background IDs to Entrez as well
    qvalueCutoff = 0.05,
    minGSSize    = 5,              # Minimum gene set size
    maxGSSize    = 100              # Maximum gene set size
  )
  
  # Get the top 3 KEGG pathways
  top_kegg <- ekegg@result %>% 
    arrange(p.adjust) %>% 
    head(3)
  
  # Store the results in the list
  kegg_results[[as.character(cluster_num)]] <- top_kegg
}

# View and Report the Results
# Print the results
for (cluster_num in names(kegg_results)) {
  cat("\nCluster:", cluster_num, "\n")
  print(kegg_results[[cluster_num]][, c("ID", "Description", "p.adjust")])
}



## KEGG pathway enrichment for each cluster v2
# Retrieve the KEGG pathway data
kegg_pathways <- keggLink("pathway", "hsa")
kegg_pathway_names <- keggList("pathway", "hsa")

# Convert the kegg_pathways list to a data frame
kegg_df <- data.frame(
    EntrezID = sub("hsa:", "", names(kegg_pathways)),
    PathwayID = sub("path:", "", kegg_pathways),
    stringsAsFactors = FALSE
)

# Merge with pathway names
kegg_df <- merge(kegg_df, kegg_pathway_names, by.x = "PathwayID", by.y = "pathway_id", all.x = TRUE)

# Check the first few rows
head(kegg_df)

# Function to perform custom enrichment
perform_kegg_enrichment <- function(entrez_ids, background_ids, kegg_df) {
  pathway_counts <- kegg_df %>%
    group_by(PathwayID, Description) %>%
    summarise(
      GeneCount = sum(EntrezID %in% entrez_ids),
      BgCount = sum(EntrezID %in% background_ids)
    ) %>%
    filter(GeneCount > 0) %>%
    mutate(
      Pvalue = phyper(GeneCount - 1, BgCount, length(background_ids) - BgCount, length(entrez_ids), lower.tail = FALSE),
      FDR = p.adjust(Pvalue, method = "BH")
    ) %>%
    arrange(FDR)
  
  return(pathway_counts)
}

# Perform KEGG enrichment for each cluster
kegg_results <- list()

for (cluster_num in top_10_clusters$Cluster_Number) {
  # Get the Uniprot IDs for the current cluster
  protein_ids <- top_10_proteins %>% 
    filter(Cluster_Number == cluster_num) %>%
    pull(Uniprot_ID)
  
  # Convert Uniprot IDs to Entrez Gene IDs
  entrez_mapping <- bitr(protein_ids, 
                         fromType = "UNIPROT", 
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)
  
  # Extract Entrez Gene IDs
  entrez_ids <- entrez_mapping$ENTREZID
  
  # Perform custom KEGG enrichment analysis
  enriched_kegg <- perform_kegg_enrichment(entrez_ids, bg_entrez_mapping$ENTREZID, kegg_df)
  
  # Store the top 3 pathways in the results list
  kegg_results[[as.character(cluster_num)]] <- head(enriched_kegg, 3)
}

# Print the results
for (cluster_num in names(kegg_results)) {
  cat("\nCluster:", cluster_num, "\n")
  print(kegg_results[[cluster_num]][, c("PathwayID", "Description", "FDR")])
}




## Reactome enrichment v2
background_ids <- unique(cluster_df$Uniprot_ID)
# Translate Uniprot IDs to Entrez Gene IDs
bg_entrez_mapping <- bitr(background_ids, 
                    fromType = "UNIPROT", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)

# Set a cluster size cutoff
cluster_size_cutoff <- 7

# Filter clusters based on the size cutoff
filtered_clusters <- cluster_size_df %>% 
  filter(Protein_Count >= cluster_size_cutoff)

# Prepare the data for the filtered clusters
filtered_proteins <- cluster_df %>% 
  filter(Cluster_Number %in% filtered_clusters$Cluster_Number)

# Reactome Enrichment Analysis for Each Cluster
reactome_results <- list()

# Loop through each filtered cluster
for (cluster_num in filtered_clusters$Cluster_Number) {
  # Get the Uniprot IDs for the current cluster
  protein_ids <- filtered_proteins %>% 
    filter(Cluster_Number == cluster_num) %>%
    pull(Uniprot_ID)

  # Translate Uniprot IDs to Entrez Gene IDs
  entrez_mapping <- bitr(protein_ids, 
                         fromType = "UNIPROT", 
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)

  # Perform Reactome pathway enrichment analysis
  ereactome <- enrichPathway(
    gene         = entrez_mapping$ENTREZID,
    organism     = "human",           
    pAdjustMethod = "BH",              
    universe     = bg_entrez_mapping$ENTREZID,    
    qvalueCutoff = 0.05,               
    minGSSize    = 5,                 
    maxGSSize    = 200,               
    readable     = TRUE               
  )
  
  # Get the top 3 Reactome pathways
  top_terms <- ereactome@result %>% 
    arrange(p.adjust) %>% 
    head(3)
  
  # Store the results in the list
  reactome_results[[as.character(cluster_num)]] <- top_terms
}

# Sort clusters by the highest p.adjust and print the results
sorted_clusters <- names(reactome_results)[order(sapply(reactome_results, function(x) min(x$p.adjust, na.rm = TRUE)))]

for (cluster_num in sorted_clusters) {
  cat("\nCluster:", cluster_num, "\n")
  print(reactome_results[[cluster_num]][, c("ID", "Description", "p.adjust")])
}




##### network clustering 
#top.all.val.dmso.vs.cxcr7

# Collapse the dataset by `uniprot`, keeping the row with the maximum absolute the column `logFC_10_dmso.vs.cxcr7`

# Define the variable for the column name
input_data = "1800_dmso.vs.cxcr7"
logFC_column <- "logFC_1800_dmso.vs.cxcr7"
adj_p_val_column <- "adj.P.Val_1800_dmso.vs.cxcr7"


## Collapse the DataFrame by Uniprot ID
# 1 collpase here
#collapsed_df <- top.all.val.dmso.vs.cxcr7 %>%
#  group_by(uniprot) %>%
#  filter(abs(.data[[logFC_column]]) == max(abs(.data[[logFC_column]]))) %>%
#  ungroup()

#head(collapsed_df)


# 2 collpase as in pl.phospho.10.1.network.reconstruction.intact.v2.R (collpase by   unirpot choosing a psite with highest logFC over all timepoints)
#collapsed_df 
head(collapsed_df)


## Prepare Gene Sets Based on Clusters
filtered_cluster_numbers <- cluster_size_df %>%
  filter(Protein_Count > 7) %>%
  pull(Cluster_Number)
filtered_cluster_numbers
filtered_cluster_df <- cluster_df %>%
  filter(Cluster_Number %in% filtered_cluster_numbers) %>%
  dplyr::select(Uniprot_ID, Cluster_Number, Symbol)
filtered_cluster_df


gene_sets <- split(filtered_cluster_df$Uniprot_ID, filtered_cluster_df$Cluster_Number)


# View the first few gene sets
str(gene_sets)

## Prepare Ranked List for GSEA
# Create a ranked list using the specified logFC column
ranked_list <- collapsed_df %>%
  arrange(desc(.data[[logFC_column]])) %>%
  dplyr::select(uniprot, .data[[logFC_column]])

# Convert to a named vector where the names are Uniprot IDs and the values are logFC
ranked_list_vector <- ranked_list[[logFC_column]]
names(ranked_list_vector) <- ranked_list$uniprot

# Check the ranked list to ensure correct ordering
print(head(ranked_list_vector))

background = collapsed_df$uniprot

# Run GSEA with your custom gene sets
fgsea_results <- fgsea(pathways = gene_sets, 
                       stats = ranked_list_vector, 
                       minSize = 6,   # Minimum gene set size
                       maxSize = 80,  # Maximum gene set size
                       nperm = 1000000)

# View the GSEA results
head(fgsea_results)

# Sort the fgsea_results by padj (adjusted p-value)
#fgsea_results <- fgsea_results %>%
#  arrange(padj)

# View the top results
fgsea_results %>% 
  arrange(padj) %>% 
  head(100)



# Extract the leading edge genes for the top pathway
top_leading_edge <- fgsea_results$leadingEdge[[1]]

# Print the leading edge genes
print(top_leading_edge)

# Translate Uniprot IDs to Gene Symbols
leading_edge_symbols <- bitr(
  top_leading_edge,
  fromType = "UNIPROT",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
)

# Print the resulting gene symbols
print(leading_edge_symbols$SYMBOL)

# Merge leading_edge_symbols with top.all.val.dmso.vs.cxcr7 or collapsed_df based on Uniprot IDs
input = top.all.val.dmso.vs.cxcr7
#input = collapsed_df

merged_df <- merge(leading_edge_symbols, input, 
                   by.x = "UNIPROT", by.y = "uniprot")

# Select the relevant columns using dplyr::select
result_df <- dplyr::select(merged_df, UNIPROT, SYMBOL, !!logFC_column, !!adj_p_val_column)



# View the resulting dataframe
# Dynamically name the result_df based on input_data
output_name <- paste0("result_df_", gsub("\\.", "_", input_data))
assign(output_name, result_df)

output_name <- paste0("fgsea_results_", gsub("\\.", "_", input_data))
assign(output_name, fgsea_results)

output_name <- paste0("ranked_list_", gsub("\\.", "_", input_data))
assign(output_name, ranked_list)

output_name <- paste0("ranked_list_df_", gsub("\\.", "_", input_data))
ranked_list = as.data.frame(ranked_list)
assign(output_name, ranked_list)

output_name <- paste0("ranked_list_vector", gsub("\\.", "_", input_data))
assign(output_name, ranked_list_vector)


# Plot enrichment for the top pathway
dev.new()
library(ggplot2)

fgsea_results = fgsea_results_1800_dmso_vs_cxcr7

# Sort the results by adjusted p-value (padj)
fgsea_results <- fgsea_results[order(fgsea_results$pathway), ]

# Extract the top pathway
top_pathway <- fgsea_results$pathway[[1]]

# Retrieve the gene set for the top pathway
gene_set <- gene_sets[[top_pathway]]

# Generate the enrichment plot for the top pathway
enrichment_plot <- plotEnrichment(
  pathway = gene_set, 
  stats = ranked_list_vector10_dmso_vs_cxcr7
) + labs(title = paste("Enrichment Plot for Pathway:", top_pathway))

# Print the plot
print(enrichment_plot)



## Enrichment of pathway info add to nodes: add NS values to the nodes table and save to import into cytoscape
## chosse the leadingedge proteins and their NES score and adjusted pvalue and add to nodes tables

pl.network.cluster.nodes <- read.delim("../../analysis/intact/pl.network.cluster.mcl.1.5.v2.csv", sep=",", header=TRUE, skip = 0)

