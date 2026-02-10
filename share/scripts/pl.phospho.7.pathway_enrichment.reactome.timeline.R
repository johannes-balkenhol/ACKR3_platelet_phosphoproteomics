#############################
## Refined Phosphoproteomics Enrichment and Visualization Script
#############################

## --- INSTALL & LOAD PACKAGES ---
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("org.Hs.eg.db", "reactome.db", "clusterProfiler", "ReactomePA"))
install.packages(c("ggplot2", "cowplot", "pheatmap", "dplyr", "tidyr", "RColorBrewer", "sjmisc"))

suppressPackageStartupMessages({
  library(annotate)       # provides getSYMBOL
  library(basicPlotteR)   # if needed for plotting helpers
  library(calibrate)
  library(clusterProfiler)
  library(cowplot)
  library(directPA)
  library(dplyr)
  library(enrichplot)
  library(ggplot2)
  library(limma)
  library(org.Hs.eg.db)
  library(pheatmap)
  library(PhosR)
  library(plyr)
  library(RColorBrewer)
  library(reactome.db)
  library(ReactomePA)
  library(remotes)
  library(rlist)
  library(sjmisc)
  library(stringr)
  library(tidyr)
})

#############################
## 1. Prepare the Pathway Annotation
#############################

# Get Reactome pathways and convert IDs to human-readable names
pathways <- as.list(reactomePATHID2EXTID)
path_names <- as.list(reactomePATHID2NAME)
name_id <- match(names(pathways), names(path_names))
names(pathways) <- unlist(path_names)[name_id]

# Restrict to Homo sapiens pathways
pathways <- pathways[grepl("Homo sapiens", names(pathways), ignore.case = TRUE)]

# Convert each pathway’s Entrez IDs to gene symbols
# (getSYMBOL from the annotate package is now available)
pathways <- lapply(pathways, function(path) {
  gene_names <- unname(getSYMBOL(path, data = "org.Hs.eg"))
  toupper(unique(gene_names))
})

#############################
## 2. Define Input Datasets
#############################

# (A) The full phosphosite tables (for later annotation & heatmap)
# These are your original filtered data (with multiple entries per protein)
phospho_inputs <- list(top.filter.10, top.filter.600, top.filter.1800, 
                       top.filter.10.dmso.vs.0s, top.filter.600.dmso.vs.0s, top.filter.1800.dmso.vs.0s,
                       top.filter.10.cxcr7.vs.0s, top.filter.600.cxcr7.vs.0s, top.filter.1800.cxcr7.vs.0s)
names_input <- c("10", "600", "1800", 
                 "10.dmso.vs.0s", "600.dmso.vs.0s", "1800.dmso.vs.0s",
                 "10.cxcr7.vs.0s", "600.cxcr7.vs.0s", "1800.cxcr7.vs.0s")
for (i in 1:length(phospho_inputs)) {
  phospho_inputs[[i]] <- phospho_inputs[[i]][order(rownames(phospho_inputs[[i]])), ]
}
names(phospho_inputs) <- names_input

# (B) The collapsed tables (one row per phosphoprotein/psite; maximum |logFC|)
# These collapsed datasets will be used for pathway enrichment.
collapsed_inputs <- list(top.collapse.10, top.collapse.600, top.collapse.1800, 
                         top.collapse.10.dmso.vs.0s, top.collapse.600.dmso.vs.0s, top.collapse.1800.dmso.vs.0s,
                         top.collapse.10.cxcr7.vs.0s, top.collapse.600.cxcr7.vs.0s, top.collapse.1800.cxcr7.vs.0s)
for (i in 1:length(collapsed_inputs)) {
  collapsed_inputs[[i]] <- collapsed_inputs[[i]][order(rownames(collapsed_inputs[[i]])), ]
}
names(collapsed_inputs) <- names_input

# Build a matrix of logFC values from the collapsed table (using second element of the rowname, e.g., gene symbol)
Tc.gene <- as.matrix(cbind(collapsed_inputs[[1]]$logFC, collapsed_inputs[[2]]$logFC, collapsed_inputs[[3]]$logFC,
                           collapsed_inputs[[4]]$logFC, collapsed_inputs[[5]]$logFC, collapsed_inputs[[6]]$logFC,
                           collapsed_inputs[[7]]$logFC, collapsed_inputs[[8]]$logFC, collapsed_inputs[[9]]$logFC))
rownames(Tc.gene) <- paste(sapply(strsplit(rownames(collapsed_inputs[[1]]), ";"), "[[", 2))
colnames(Tc.gene) <- names_input


#############################
## 3. Pathway Enrichment using Collapsed Data
#############################

# Function to run enrichment analysis for one dataset and one regulation type.
# Here we use enrichPathway (ReactomePA) as an example.
run_enrichment <- function(data, regulation = c("up", "down"), pathway_list) {
  regulation <- match.arg(regulation)
  
  # Select sites based on logFC sign
  if(regulation == "up") {
    selected_sites <- rownames(data)[data$logFC > 0]
  } else {
    selected_sites <- rownames(data)[data$logFC < 0]
  }
  
  if(length(selected_sites) < 1) return(NULL)
  
  # Run enrichment (adjust pvalue cutoff and method as desired)
  enrich_res <- enrichPathway(gene         = selected_sites,
                              organism     = "human",
                              pvalueCutoff = 0.05,
                              readable     = TRUE)
  if(is.null(enrich_res) || nrow(as.data.frame(enrich_res)) < 1) return(NULL)
  
  res_df <- as.data.frame(enrich_res)
  
  # For each pathway, get pathway size from the pathway list and calculate ratio = (Count/pathway_size)
  res_df$pw_size <- sapply(res_df$ID, function(pid) {
    pw <- pathway_list[[pid]]
    if(is.null(pw)) NA else length(pw)
  })
  res_df$ratio <- res_df$Count / res_df$pw_size
  res_df$regulation <- regulation
  return(res_df)
}

# Run enrichment for each timepoint (using collapsed_inputs) for up and down separately.
enrich_results <- list()
for(time in names_input) {
  dat <- collapsed_inputs[[time]]
  up_enrich   <- run_enrichment(dat, regulation = "up", pathway_list = pathways)
  down_enrich <- run_enrichment(dat, regulation = "down", pathway_list = pathways)
  
  # Select top 5 pathways (by adjusted p-value) if available
  if(!is.null(up_enrich)) {
    up_top <- up_enrich %>% arrange(p.adjust) %>% head(5)
    up_top$time <- time
  } else { up_top <- NULL }
  
  if(!is.null(down_enrich)) {
    down_top <- down_enrich %>% arrange(p.adjust) %>% head(5)
    down_top$time <- time
  } else { down_top <- NULL }
  
  enrich_results[[time]] <- bind_rows(up_top, down_top)
}
# Combine all enrichment results
enrich_summary <- bind_rows(enrich_results)

#############################
## 4. Bar Plot: Enrichment (Percentage of Sites in Pathway) with p-value
#############################

# For plotting, assign positive ratio for up-regulated and negative for down-regulated.
enrich_summary <- enrich_summary %>% 
  mutate(plot_ratio = ifelse(regulation == "up", ratio, -ratio))

# Bar plot: x-axis = time; y-axis = % (number of sites/ pathway size); annotate with p-value.
bar_plot <- ggplot(enrich_summary, aes(x = time, y = plot_ratio, fill = regulation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5)) +
  geom_text(aes(label = paste0("p=", round(p.adjust, 4))),
            vjust = ifelse(enrich_summary$plot_ratio > 0, -0.5, 1.5),
            position = position_dodge(width = 0.5),
            size = 3) +
  labs(title = "Pathway Enrichment: % of Phosphosites per Pathway",
       y = "Enrichment Ratio (Phosphosites / Pathway Size)",
       x = "Time (sec)") +
  theme_cowplot() +
  scale_fill_manual(values = c("up" = "red", "down" = "blue"))

# Save bar plot
tiff(filename = "../analysis/Reactome_enrichment/enrichment_barplot.tiff",
     width = 8 * 300, height = 6 * 300, res = 300, compression = "lzw")
print(bar_plot)
dev.off()

#############################
## 5. Retrieve Phosphosite Details from Filtered Data for Top Enriched Pathways
#############################

# From the enrichment summary, get the unique set of enriched pathway IDs.
top_pathways <- unique(enrich_summary$ID)

# For each timepoint (using the full phosphosite table) extract phosphosites in the top pathways.
# Also, annotate each site with its pathway membership. If a site is present in multiple pathways,
# concatenate the pathway names.
phospho_details <- list()

for(time in names_input) {
  dat <- phospho_inputs[[time]]
  
  # Assume rownames are in the format "protein;gene_symbol" (or similar).
  gene_symbols <- sapply(strsplit(rownames(dat), ";"), function(x) {
    if(length(x) >= 2) x[2] else x[1]
  })
  
  # Create an empty data.frame to store annotated sites
  annotated_sites <- data.frame()
  
  for(pid in top_pathways) {
    pathway_genes <- pathways[[pid]]
    if(is.null(pathway_genes)) next
    
    # Select sites where the gene symbol is in the pathway
    sel_idx <- which(gene_symbols %in% pathway_genes)
    if(length(sel_idx) > 0) {
      temp <- dat[sel_idx, , drop = FALSE]
      # Add a column for pathway membership (if already exists, concatenate)
      temp$Pathway <- pid
      # If a phosphosite belongs to more than one pathway, we will combine them later.
      annotated_sites <- rbind(annotated_sites, temp)
    }
  }
  
  # Aggregate pathway membership for each phosphosite (using rownames as unique identifier)
  if(nrow(annotated_sites) > 0) {
    annotated_sites$siteID <- rownames(annotated_sites)
    aggregated <- annotated_sites %>% 
      group_by(siteID) %>% 
      summarise(logFC = mean(logFC),   # or choose appropriate summary if needed
                Pathway = paste(unique(Pathway), collapse = ";"))
    rownames(aggregated) <- aggregated$siteID
    phospho_details[[time]] <- as.data.frame(aggregated)
  }
}

#############################
## 6. Build a Combined Matrix and Generate Heatmap
#############################

# Merge the logFC values across timepoints for the annotated phosphosites.
# Note: Not all timepoints may have the same sites; use a full merge.
phospho_matrix_list <- lapply(names(phospho_details), function(t) {
  df <- phospho_details[[t]][, "logFC", drop = FALSE]
  colnames(df) <- t
  df
})
phospho_matrix <- Reduce(function(x, y) merge(x, y, by = "row.names", all = TRUE),
                           phospho_matrix_list)
rownames(phospho_matrix) <- phospho_matrix$Row.names
phospho_matrix <- as.matrix(phospho_matrix[,-1])

# Generate an overall heatmap with clustering.
tiff(filename = "../analysis/Reactome_enrichment/phosphosite_heatmap.tiff",
     width = 8 * 300, height = 10 * 300, res = 300, compression = "lzw")
pheatmap(phospho_matrix, 
         main = "Log2FC of Phosphosites in Top Enriched Pathways\n(Pathway membership annotated in row names)",
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         fontsize_row = 6)
dev.off()

#############################
## 7. Optional: Cluster Dendrogram and Subgroup Heatmaps
#############################

# If you want to extract subgroups from the dendrogram and plot separate heatmaps:
# dend <- as.dendrogram(hclust(dist(phospho_matrix)))
# groups <- cutree(dend, k = 3)  # adjust number of clusters as needed
# for(grp in unique(groups)) {
#   sub_matrix <- phospho_matrix[names(groups)[groups == grp], ]
#   if(nrow(sub_matrix) >= 5) {  # only plot if the subgroup has enough sites
#     tiff(filename = paste0("../analysis/Reactome_enrichment/phosphosite_heatmap_cluster_", grp, ".tiff"),
#          width = 6 * 300, height = 6 * 300, res = 300, compression = "lzw")
#     pheatmap(sub_matrix,
#              main = paste("Cluster", grp),
#              clustering_distance_rows = "euclidean",
#              clustering_method = "complete",
#              fontsize_row = 8)
#     dev.off()
#   }
# }
