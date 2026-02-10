#############################
## Simplified Phosphoproteomic Analysis
## 1. Extract Clathrin pathway phosphosites
## 2. Hallmark-only GSEA (ROS, cAMP, PKC relevant sets)
#############################

## --- INSTALL & LOAD PACKAGES ---
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("annotate", "org.Hs.eg.db", "reactome.db", 
                       "clusterProfiler", "ReactomePA", "msigdbr", "fgsea"))

install.packages(c("ggplot2", "cowplot", "pheatmap", "dplyr", "tidyr", 
                   "RColorBrewer", "stringr"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(annotate)
  library(reactome.db)
  library(org.Hs.eg.db)
  library(RColorBrewer)
  library(pheatmap)
  library(ggplot2)
  library(enrichplot)
  library(cowplot)
  library(msigdbr)
  library(fgsea)
})

#############################
## 1. Prepare Input Data
#############################

## Your collapsed input dataframes
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

## Build logFC matrix
Tc.gene <- as.matrix(cbind(input[[1]]$logFC, input[[2]]$logFC, input[[3]]$logFC,
                           input[[4]]$logFC, input[[5]]$logFC, input[[6]]$logFC,
                           input[[7]]$logFC, input[[8]]$logFC, input[[9]]$logFC))

rownames(Tc.gene) <- paste(sapply(strsplit(rownames(input[[1]]), ";"), "[[", 2))
colnames(Tc.gene) <- names_input

## Create reference dataframe with phosphosite info
phospho.info <- data.frame(
  phosphosite = rownames(input[[1]]),
  gene = rownames(Tc.gene),
  logFC.10 = input[[1]]$logFC,
  logFC.600 = input[[2]]$logFC,
  logFC.1800 = input[[3]]$logFC,
  stringsAsFactors = FALSE
)

#############################
## 2. Load Reactome Pathways
#############################

pathways <- as.list(reactomePATHID2EXTID)
path_names <- as.list(reactomePATHID2NAME)
name_id <- match(names(pathways), names(path_names))
names(pathways) <- unlist(path_names)[name_id]
pathways <- pathways[grepl("Homo sapiens", names(pathways), ignore.case = TRUE)]
reactome_list <- lapply(pathways, function(path) {
  gene_name <- unname(getSYMBOL(path, data = "org.Hs.eg"))
  toupper(unique(gene_name))
})

message("Total Reactome pathways loaded: ", length(reactome_list))

#############################
## 3. Extract Clathrin Pathway Phosphosites
#############################

# Find Clathrin pathway in Reactome
clathrin_pathway_name <- names(reactome_list)[
  grepl("clathrin|endocytosis", names(reactome_list), ignore.case = TRUE)
][1]  # Take first match

if (!is.na(clathrin_pathway_name)) {
  
  message("\n=== CLATHRIN PATHWAY EXTRACTION ===")
  message("Pathway found: ", clathrin_pathway_name)
  
  # Get genes in Clathrin pathway
  clathrin_genes <- reactome_list[[clathrin_pathway_name]]
  message("Total genes in pathway: ", length(clathrin_genes))
  
  # Extract phosphosites from these genes
  clathrin_phospho <- phospho.info %>% 
    filter(toupper(gene) %in% toupper(clathrin_genes)) %>% 
    arrange(desc(abs(logFC.10)))
  
  message("Phosphosites found in Clathrin pathway: ", nrow(clathrin_phospho))
  
  # Save Clathrin phosphosites
  write.csv(clathrin_phospho, 
            file = "clathrin_pathway_phosphosites.csv", 
            row.names = FALSE)
  
  message("\nTop Clathrin phosphosites:")
  print(head(clathrin_phospho, 10))
  
  # Visualize Clathrin phosphosites across timepoints
  if (nrow(clathrin_phospho) > 0) {
    
    clathrin_mat <- cbind(
      clathrin_phospho$logFC.10,
      clathrin_phospho$logFC.600,
      clathrin_phospho$logFC.1800
    )
    rownames(clathrin_mat) <- clathrin_phospho$phosphosite
    colnames(clathrin_mat) <- c("10s", "600s", "1800s")
    
    # Z-score
    clathrin_mat_z <- scale(clathrin_mat)
    
    pdf("heatmap_clathrin_pathway_phosphosites.pdf", width = 8, height = 10)
    pheatmap(
      clathrin_mat_z,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      main = paste0("Clathrin Pathway Phosphosites (", nrow(clathrin_phospho), " sites)"),
      breaks = seq(-2, 2, by = 0.1),
      color = colorRampPalette(c("blue", "white", "red"))(40),
      display_numbers = round(clathrin_mat, 2),
      fontsize_number = 8
    )
    dev.off()
    
    message("\nSaved: heatmap_clathrin_pathway_phosphosites.pdf")
  }
  
} else {
  message("WARNING: Clathrin pathway not found in Reactome database")
  message("Available pathways with 'clathrin' or 'endocytosis':")
  print(names(reactome_list)[grepl("clathrin|endocytosis", names(reactome_list), ignore.case = TRUE)])
}

#############################
## 4. Load Hallmark Gene Sets
#############################

hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)

message("\n=== HALLMARK GENE SETS ===")
message("Total Hallmark sets: ", length(hallmark_list))

# Show Hallmark sets and which are relevant to ROS, cAMP, PKC
hallmark_names <- names(hallmark_list)

ros_hallmarks <- hallmark_names[grepl("OXIDATIVE|REACTIVE_OXYGEN|HYPOXIA|STRESS|APOPTOSIS|MITOCHONDRIA", 
                                      hallmark_names, ignore.case = TRUE)]

camp_hallmarks <- hallmark_names[grepl("CAMP|SIGNALING|RESPONSE", 
                                       hallmark_names, ignore.case = TRUE)]

message("\nHallmark sets related to ROS/Stress:")
print(ros_hallmarks)

message("\nHallmark sets related to signaling:")
print(camp_hallmarks)

#############################
## 5. Hallmark-only GSEA (3 main timepoints)
#############################

gsea_hallmark_results <- list()

timepoints_main <- c("10", "600", "1800")

for (tp in timepoints_main) {
  
  message("\n--- Processing timepoint: ", tp, "s ---")
  
  # Get logFC for this timepoint
  logfc_vec <- Tc.gene[, tp]
  names(logfc_vec) <- rownames(Tc.gene)
  
  # Sort by logFC (descending) for GSEA
  logfc_ranked <- sort(logfc_vec, decreasing = TRUE)
  
  # Perform GSEA with Hallmark gene sets only
  gsea_res <- fgsea(
    pathways = hallmark_list,
    stats = logfc_ranked,
    minSize = 3,
    maxSize = 5000,
    eps = 1e-10
  )
  
  # Add timepoint
  gsea_res$timepoint <- tp
  gsea_hallmark_results[[tp]] <- as.data.frame(gsea_res)
  
  message("  Significant pathways (padj < 0.05): ", 
          sum(gsea_res$padj < 0.05, na.rm = TRUE))
}

# Combine all results
gsea_hallmark_all <- bind_rows(gsea_hallmark_results)

# Save full results
write.csv(gsea_hallmark_all, 
          file = "hallmark_GSEA_enrichment_all.csv", 
          row.names = FALSE)

message("\n=== TOP SIGNIFICANT HALLMARK PATHWAYS (padj < 0.05) ===")
gsea_top <- gsea_hallmark_all %>% 
  filter(padj < 0.05) %>% 
  arrange(padj) %>% 
  head(15)
print(gsea_top[, c("pathway", "ES", "NES", "padj", "timepoint")])

#############################
## 6. Classify Hallmark sets by category
#############################

gsea_hallmark_all <- gsea_hallmark_all %>% 
  mutate(
    category = case_when(
      grepl("OXIDATIVE|REACTIVE_OXYGEN|MITOCHONDRIAL|HYPOXIA", pathway, ignore.case = TRUE) ~ "ROS/Stress",
      grepl("APOPTOSIS|INFLAMMATORY|IMMUNE|TNF", pathway, ignore.case = TRUE) ~ "Apoptosis/Immune",
      grepl("ANGIOGENESIS|METASTASIS|INVASION", pathway, ignore.case = TRUE) ~ "EMT/Metastasis",
      grepl("PROLIFERATION|CELL_CYCLE", pathway, ignore.case = TRUE) ~ "Proliferation",
      grepl("EPITHELIAL|MESENCHYMAL", pathway, ignore.case = TRUE) ~ "EMT",
      grepl("MYC|E2F|RB", pathway, ignore.case = TRUE) ~ "Oncogenic",
      grepl("WONT|HEDGEHOG|NOTCH|TGF", pathway, ignore.case = TRUE) ~ "Signaling Pathways",
      grepl("DNA|UV|RADIATION|P53", pathway, ignore.case = TRUE) ~ "DNA Damage",
      TRUE ~ "Other"
    )
  )

# Save with categories
write.csv(gsea_hallmark_all, 
          file = "hallmark_GSEA_enrichment_categorized.csv", 
          row.names = FALSE)

#############################
## 7. Visualization 1: Bubble plot by category and timepoint
#############################

gsea_sig <- gsea_hallmark_all %>% filter(padj < 0.05)

if (nrow(gsea_sig) > 0) {
  
  p1 <- ggplot(gsea_sig, aes(x = NES, y = -log10(padj), 
                             color = category, size = abs(ES))) +
    geom_point(alpha = 0.6) +
    facet_wrap(~timepoint, nrow = 1) +
    scale_color_brewer(palette = "Set2") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    labs(title = "Hallmark Pathway Enrichment Across Timepoints",
         x = "Normalized Enrichment Score (NES)",
         y = "-log10(padj)",
         color = "Category",
         size = "Effect Size") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 10),
      legend.position = "bottom"
    )
  
  pdf("hallmark_bubble_plot.pdf", width = 12, height = 6)
  print(p1)
  dev.off()
  
  message("\nSaved: hallmark_bubble_plot.pdf")
}

#############################
## 8. Visualization 2: Heatmap of Hallmark pathways
#############################

gsea_hm <- gsea_hallmark_all %>% 
  filter(padj < 0.1) %>% 
  arrange(desc(abs(NES)))

if (nrow(gsea_hm) > 0) {
  
  # Pivot to wide format
  gsea_hm_wide <- gsea_hm %>% 
    dplyr::select(pathway, NES, timepoint) %>% 
    pivot_wider(names_from = timepoint, values_from = NES, values_fill = 0) %>% 
    column_to_rownames("pathway") %>% 
    as.matrix()
  
  # Reorder columns
  gsea_hm_wide <- gsea_hm_wide[, c("10", "600", "1800")]
  
  pdf("hallmark_heatmap_NES.pdf", width = 6, height = 10)
  pheatmap(
    gsea_hm_wide,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    main = "Hallmark Pathway Enrichment (NES)",
    breaks = seq(-3, 3, by = 0.1),
    color = colorRampPalette(c("blue", "white", "red"))(60),
    display_numbers = round(gsea_hm_wide, 2),
    fontsize_number = 8
  )
  dev.off()
  
  message("Saved: hallmark_heatmap_NES.pdf")
}

#############################
## 9. Summary Statistics
#############################

summary_stats <- data.frame(
  Analysis = c(
    "Total phosphosites",
    "Clathrin pathway phosphosites",
    "Total Hallmark sets tested",
    "Significant Hallmark pathways (padj<0.05) - 10s",
    "Significant Hallmark pathways (padj<0.05) - 600s",
    "Significant Hallmark pathways (padj<0.05) - 1800s",
    "ROS/Stress-related pathways enriched",
    "Apoptosis/Immune pathways enriched"
  ),
  Count = c(
    nrow(phospho.info),
    if (exists("clathrin_phospho")) nrow(clathrin_phospho) else 0,
    length(hallmark_list),
    sum(gsea_hallmark_all$padj < 0.05 & gsea_hallmark_all$timepoint == "10", na.rm = TRUE),
    sum(gsea_hallmark_all$padj < 0.05 & gsea_hallmark_all$timepoint == "600", na.rm = TRUE),
    sum(gsea_hallmark_all$padj < 0.05 & gsea_hallmark_all$timepoint == "1800", na.rm = TRUE),
    sum(gsea_sig$category == "ROS/Stress", na.rm = TRUE),
    sum(gsea_sig$category == "Apoptosis/Immune", na.rm = TRUE)
  )
)

print("\n=== SUMMARY ===")
print(summary_stats)
write.csv(summary_stats, file = "analysis_summary.csv", row.names = FALSE)

#############################
## 10. Detailed pathway tables
#############################

# ROS/Stress pathways
ros_results <- gsea_hallmark_all %>% 
  filter(category == "ROS/Stress") %>% 
  arrange(timepoint, padj)

write.csv(ros_results, 
          file = "hallmark_ROS_stress_pathways.csv", 
          row.names = FALSE)

message("\n=== ROS/STRESS PATHWAYS ===")
print(ros_results)

# Apoptosis/Immune pathways
apoptosis_results <- gsea_hallmark_all %>% 
  filter(category == "Apoptosis/Immune") %>% 
  arrange(timepoint, padj)

write.csv(apoptosis_results, 
          file = "hallmark_apoptosis_immune_pathways.csv", 
          row.names = FALSE)

# All categories summary
category_summary <- gsea_sig %>% 
  group_by(timepoint, category) %>% 
  summarise(
    count = n(),
    mean_NES = mean(NES),
    mean_padj = mean(padj),
    .groups = "drop"
  ) %>% 
  arrange(timepoint)

write.csv(category_summary, 
          file = "hallmark_category_summary.csv", 
          row.names = FALSE)

message("\n=== PATHWAY CATEGORIES BY TIMEPOINT ===")
print(category_summary)

#############################
## FINAL OUTPUT
#############################

message("\n\n=== ANALYSIS COMPLETE ===")
message("\nOutput files created:")
message("  Clathrin Pathway:")
message("  ✓ clathrin_pathway_phosphosites.csv")
message("  ✓ heatmap_clathrin_pathway_phosphosites.pdf")
message("\n  Hallmark GSEA Results:")
message("  ✓ hallmark_GSEA_enrichment_all.csv (all results)")
message("  ✓ hallmark_GSEA_enrichment_categorized.csv (with categories)")
message("  ✓ hallmark_ROS_stress_pathways.csv")
message("  ✓ hallmark_apoptosis_immune_pathways.csv")
message("  ✓ hallmark_category_summary.csv")
message("\n  Visualizations:")
message("  ✓ hallmark_bubble_plot.pdf")
message("  ✓ hallmark_heatmap_NES.pdf")
message("\n  Summary:")
message("  ✓ analysis_summary.csv")