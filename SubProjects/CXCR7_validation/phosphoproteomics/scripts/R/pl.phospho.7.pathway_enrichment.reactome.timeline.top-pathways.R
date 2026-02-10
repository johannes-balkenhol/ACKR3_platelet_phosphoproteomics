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
## 1. Prepare Collapsed Input & Pathway Annotation
#############################

## Define input2 (collapsed tables)
input <- list(top.collapse.10, top.collapse.600, top.collapse.1800, 
              top.collapse.10.dmso.vs.0s, top.collapse.600.dmso.vs.0s, top.collapse.1800.dmso.vs.0s,
              top.collapse.10.cxcr7.vs.0s, top.collapse.600.cxcr7.vs.0s, top.collapse.1800.cxcr7.vs.0s)

names_input <- c("10", "600", "1800", 
                 "10.dmso.vs.0s", "600.dmso.vs.0s", "1800.dmso.vs.0s",
                 "10.cxcr7.vs.0s", "600.cxcr7.vs.0s", "1800.cxcr7.vs.0s")

for (i in seq_along(input)) {
  input[[i]] <- input[[i]][order(rownames(input[[i]])), ]
}
names(input) <- names_input

## Build a logFC matrix from the collapsed tables.
Tc.gene <- as.matrix(cbind(input[[1]]$logFC, input[[2]]$logFC, input[[3]]$logFC,
                           input[[4]]$logFC, input[[5]]$logFC, input[[6]]$logFC,
                           input[[7]]$logFC, input[[8]]$logFC, input[[9]]$logFC))
rownames(Tc.gene) <- paste(sapply(strsplit(rownames(input[[1]]), ";"), "[[", 2))
colnames(Tc.gene) <- names_input

## Prepare Reactome pathways (only Homo sapiens)
pathways <- as.list(reactomePATHID2EXTID)
path_names <- as.list(reactomePATHID2NAME)
name_id <- match(names(pathways), names(path_names))
names(pathways) <- unlist(path_names)[name_id]
pathways <- pathways[grepl("Homo sapiens", names(pathways), ignore.case = TRUE)]
pathways <- lapply(pathways, function(path) {
  gene_name <- unname(getSYMBOL(path, data = "org.Hs.eg"))
  toupper(unique(gene_name))
})

#############################
## 2. Rank-Based Enrichment per Timepoint: Up vs. Down
#############################

# We'll process only the three main timepoints (first three elements)
up_results <- list()
down_results <- list()


# Loop over timepoints (i=1:3) for UP-regulated (alter = "greater")
for (i in 1:3) {
  res <- pathwayRankBasedEnrichment(Tc.gene[, i], annotation = pathways, alter = "greater")
  
  # Merge enrichment result with pathway sizes
  pathways2 <- as.data.frame(t(data.frame(lapply(pathways, length))))
  path3 <- as.data.frame(res)
  path4 <- merge(path3, pathways2, by = "row.names", all = FALSE)
  colnames(path4) <- c("pathway", "pvalue", "number.substrates", "substrates", "pw.size")
  
  # Clean and convert numeric columns
  path4 <- na.omit(path4)
  path4$number.substrates <- as.numeric(path4$number.substrates)
  path4$pw.size <- as.numeric(path4$pw.size)
  path4$pvalue <- as.numeric(path4$pvalue)
  path4$ratio <- path4$number.substrates / path4$pw.size
  
  # Clean pathway names
  path4$pathway <- gsub("_", " ", gsub("REACTOME_|Homo.sapiens..", "", path4$pathway))
  
  # Filter significant pathways (pvalue < 0.05) and sort by ratio (descending)
  #path4 <- path4[path4$pvalue < 0.05, ]
  path4 <- path4[order(path4$ratio, decreasing = TRUE), ]
  
  up_results[[ names_input[i] ]] <- path4
}

# Loop over timepoints (i=1:3) for DOWN-regulated (alter = "less")
for (i in 1:3) {
  res <- pathwayRankBasedEnrichment(Tc.gene[, i], annotation = pathways, alter = "less")
  
  pathways2 <- as.data.frame(t(data.frame(lapply(pathways, length))))
  path3 <- as.data.frame(res)
  path4 <- merge(path3, pathways2, by = "row.names", all = FALSE)
  colnames(path4) <- c("pathway", "pvalue", "number.substrates", "substrates", "pw.size")
  
  path4 <- na.omit(path4)
  path4$number.substrates <- as.numeric(path4$number.substrates)
  path4$pw.size <- as.numeric(path4$pw.size)
  path4$pvalue <- as.numeric(path4$pvalue)
  path4$ratio <- path4$number.substrates / path4$pw.size
  
  path4$pathway <- gsub("_", " ", gsub("REACTOME_|Homo.sapiens..", "", path4$pathway))
  
  #path4 <- path4[path4$pvalue < 0.05, ]
  path4 <- path4[order(path4$ratio, decreasing = TRUE), ]
  
  down_results[[ names_input[i] ]] <- path4
}

#############################
## 3. Save Enrichment Results
#############################

# Save UP-regulated enrichment tables
for (time in names(up_results)) {
  write.table(up_results[[time]],
              file = paste0("D:/Research/CXCR7_platelet_analysis/CXCR7_validation/phosphoproteomics/analysis/Reactome_enrichment/greater_", time, "_pathways.txt"),
              sep = "\t", row.names = FALSE)
}

# Save DOWN-regulated enrichment tables
for (time in names(down_results)) {
  write.table(down_results[[time]],
              file = paste0("D:/Research/CXCR7_platelet_analysis/CXCR7_validation/phosphoproteomics/analysis/Reactome_enrichment/less_", time, "_pathways.txt"),
              sep = "\t", row.names = FALSE)
}




#############################
## 4. Extract Top 5 Pathways (by ratio) per timepoint for each condition and Visualize
#############################

# For up-regulated:
top_up <- lapply(names(up_results)[1:3], function(time) {
  df <- up_results[[time]]
  if(nrow(df) > 0) {
    df_top <- df %>% arrange((pvalue)) %>% head(5)
    df_top$time <- time
    df_top$regulation <- "up"
    return(df_top)
  } else {
    return(NULL)
  }
})
top_up <- bind_rows(top_up)

# For down-regulated:
top_down <- lapply(names(down_results)[1:3], function(time) {
  df <- down_results[[time]]
  if(nrow(df) > 0) {
    df_top <- df %>% arrange((pvalue)) %>% head(5)
    df_top$time <- time
    df_top$regulation <- "down"
    # For plotting, convert ratio to negative values:
    df_top$ratio <- -df_top$ratio
    return(df_top)
  } else {
    return(NULL)
  }
})
top_down <- bind_rows(top_down)


# --- Combine top pathway data for each condition ---
# For plotting, we want a common grouping by pathway. If the same pathway is found at multiple timepoints, we'll show them together.


## First, extract the top pathways (from your previous step):
up_summary <- top_up %>% 
  dplyr::select(pathway, ratio, pvalue, time, regulation)
down_summary <- top_down %>% 
  dplyr::select(pathway, ratio, pvalue, time, regulation)



## Build a long-format summary for all up/down-regulated results from timepoints 10, 600, 1800.
up_summary_long_all <- lapply(c("10", "600", "1800"), function(t) {
  if (!is.null(up_results[[t]])) {
    df <- up_results[[t]] %>% 
      dplyr::select(pathway, pvalue, ratio, substrates) %>% 
      mutate(time = t)
    return(df)
  } else {
    return(NULL)
  }
}) %>% bind_rows()

down_summary_long_all <- lapply(c("10", "600", "1800"), function(t) {
  if (!is.null(down_results[[t]])) {
    df <- down_results[[t]] %>% 
      dplyr::select(pathway, pvalue, ratio, substrates) %>% 
      mutate(time = t)
    return(df)
  } else {
    return(NULL)
  }
}) %>% bind_rows()



## Filter long-format up/down data to include only pathways found in the top summary.
up_summary_long <- up_summary_long_all %>% 
  filter(pathway %in% up_summary$pathway) %>% 
  mutate(pvalue = as.numeric(pvalue),
         ratio = as.numeric(ratio),
         neg_log10_p = -log10(pvalue))

down_summary_long <- down_summary_long_all %>% 
  filter(pathway %in% down_summary$pathway) %>% 
  mutate(pvalue = as.numeric(pvalue),
         ratio = as.numeric(ratio),
         neg_log10_p = -log10(pvalue))



# Transform to wide format:
up_summary_long <- up_summary_long %>% 
  mutate(ratio = as.numeric(ratio),
         neg_log10_p = as.numeric(neg_log10_p))

down_summary_long <- down_summary_long %>% 
  mutate(ratio = as.numeric(ratio),
         neg_log10_p = as.numeric(neg_log10_p))


up_wide <- dplyr::select(up_summary_long, time, pathway, neg_log10_p) %>% 
  pivot_wider(names_from = pathway, 
              values_from = neg_log10_p,
              values_fill = list(neg_log10_p = 0))

down_wide <- dplyr::select(down_summary_long, time, pathway, neg_log10_p) %>% 
  pivot_wider(names_from = pathway, 
              values_from = neg_log10_p,
              values_fill = list(neg_log10_p = 0))

  

# Convert to matrix: Remove the "time" column and then assign rownames from it
up_mat <- as.matrix(up_wide[,-1])
rownames(up_mat) <- up_wide$time

# Define a color vector for each timepoint (adjust colors as needed)
time_colors <- c("#E6E6FA", "#BA55D3", "#4B0082")

dev.new()
# Create grouped barplot
barplot(up_mat, 
         beside = TRUE, 
         col = time_colors, 
         border = "white", 
         legend.text = rownames(up_mat), 
         xlab = "Pathway",
         main = "Up-regulated Pathway Enrichment (-log10(pvalue))",
         font.axis = 6,
         font.lab = 6,
         las = 2,         # makes axis labels vertical
         cex.names = 0.8) # reduces the font size of axis names

# Convert to matrix: Remove the "time" column and then assign rownames from it
down_mat <- as.matrix(down_wide[,-1])
rownames(down_mat) <- down_wide$time

# Define a color vector for each timepoint (adjust colors as needed)
time_colors <- c("#E6E6FA", "#BA55D3", "#4B0082")

dev.new()
# Create grouped barplot
# Create grouped barplot
barplot(down_mat, 
         beside = TRUE, 
         col = time_colors, 
         border = "white", 
         legend.text = rownames(up_mat), 
         xlab = "Pathway",
         main = "Up-regulated Pathway Enrichment (-log10(pvalue))",
         font.axis = 6,
         font.lab = 6,
         las = 2,         # makes axis labels vertical
         cex.names = 0.8) # reduces the font size of axis names
