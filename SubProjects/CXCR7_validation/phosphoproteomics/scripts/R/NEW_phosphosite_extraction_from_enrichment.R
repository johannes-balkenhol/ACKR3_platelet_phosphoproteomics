###############################################################
## NEW SCRIPT: EXTRACT & ANALYZE PHOSPHOSITES FROM SELECTED PATHWAYS
## For use AFTER pl_phospho_7_pathway_enrichment_reactome_timeline_selected-pathways_comparative_04122025.R
## 
## This script:
## 1. Uses enrichment results (up_results, down_results)
## 2. Uses selected pathways defined in main script
## 3. Extracts phosphosites from those pathways
## 4. Creates visualizations: Panels A-F (pathway p-values + phosphosite heatmaps)
## 5. Generates publication-ready tables
###############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(stringr)
  library(ComplexHeatmap)
})

###############################################################
## PREREQUISITE: Must have run the main enrichment script first!
## Required objects in R environment:
## - up_results, down_results (from enrichment)
## - selected_pathways (pathway list)
## - dfs_new_intersect (validation datasets - phosphosite level)
## - all_inputs_collapsed (gene level)
###############################################################

cat("\n", "="*70, "\n")
cat("PHOSPHOSITE EXTRACTION FROM SELECTED PATHWAYS\n")
cat("="*70, "\n\n")

# Verify prerequisites
required_objects <- c("up_results", "down_results", "selected_pathways", "dfs_new_intersect")
missing <- required_objects[!sapply(required_objects, exists)]

if (length(missing) > 0) {
  cat("⚠️ ERROR: Missing required objects:\n")
  print(missing)
  cat("\nPlease run main enrichment script first:\n")
  cat("source('pl_phospho_7_pathway_enrichment_reactome_timeline_selected-pathways_comparative_04122025.R')\n")
  stop("Missing prerequisites")
}

cat("✓ All prerequisites found\n")
cat("✓ Selected pathways:", length(selected_pathways), "\n")
cat("✓ Enrichment results loaded\n\n")

###############################################################
## 1. CREATE OUTPUT DIRECTORY
###############################################################

outdir <- file.path(".", "phosphosite_selected_pathways_analysis")
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

cat("Output directory:", outdir, "\n\n")

###############################################################
## 2. EXTRACT PATHWAY GENES FROM ENRICHMENT RESULTS
###############################################################

cat("Step 1: Extracting genes for each selected pathway...\n")

pathway_genes_list <- list()

for (pathway in selected_pathways) {
  
  genes_list <- character()
  
  # From UP-regulated results
  for (contrast in names(up_results)) {
    if (!is.null(up_results[[contrast]])) {
      res <- up_results[[contrast]]
      # Handle different column names (Gene vs gene)
      gene_col <- if ("Gene" %in% colnames(res)) "Gene" else if ("gene" %in% colnames(res)) "gene" else NULL
      
      if (!is.null(gene_col)) {
        genes <- res[res$Pathway == pathway, gene_col]
        genes_list <- c(genes_list, genes)
      }
    }
  }
  
  # From DOWN-regulated results
  for (contrast in names(down_results)) {
    if (!is.null(down_results[[contrast]])) {
      res <- down_results[[contrast]]
      gene_col <- if ("Gene" %in% colnames(res)) "Gene" else if ("gene" %in% colnames(res)) "gene" else NULL
      
      if (!is.null(gene_col)) {
        genes <- res[res$Pathway == pathway, gene_col]
        genes_list <- c(genes_list, genes)
      }
    }
  }
  
  pathway_genes_list[[pathway]] <- unique(genes_list)
}

cat("✓ Extracted genes for", length(pathway_genes_list), "pathways\n")

# Count genes per pathway
genes_per_path <- sapply(pathway_genes_list, length)
cat("  Range:", min(genes_per_path), "-", max(genes_per_path), "genes per pathway\n\n")

###############################################################
## 3. EXTRACT PHOSPHOSITES FROM SELECTED PATHWAYS
###############################################################

cat("Step 2: Extracting phosphosites from selected pathways...\n")

# Define the 9 contrast names (validation only - 3 comparisons × 3 timepoints)
contrast_names <- c(
  "val_10.dmso.vs.cxcr7",   "val_600.dmso.vs.cxcr7",   "val_1800.dmso.vs.cxcr7",
  "val_10.dmso.vs.0s",      "val_600.dmso.vs.0s",      "val_1800.dmso.vs.0s",
  "val_10.cxcr7.vs.0s",     "val_600.cxcr7.vs.0s",     "val_1800.cxcr7.vs.0s"
)

# Map to actual names in dfs_new_intersect
actual_names <- names(dfs_new_intersect)[1:9]

phosphosite_all <- data.frame()

for (i in seq_along(actual_names)) {
  
  contrast <- actual_names[i]
  contrast_label <- contrast_names[i]
  
  df <- dfs_new_intersect[[contrast]]
  
  if (is.null(df) || nrow(df) == 0) {
    cat("⚠️ Skipping:", contrast, "\n")
    next
  }
  
  cat("Processing:", contrast_label, "...")
  
  # For each pathway, extract phosphosites
  for (pathway in selected_pathways) {
    
    genes_in_pathway <- pathway_genes_list[[pathway]]
    
    if (length(genes_in_pathway) == 0) next
    
    # Filter to genes in this pathway
    phos_subset <- df %>%
      filter(toupper(name) %in% toupper(genes_in_pathway)) %>%
      mutate(
        condition = contrast_label,
        pathway = pathway,
        psite_id = paste(name, PSite, sep = "_")
      ) %>%
      select(psite_id, name, PSite, logFC, PValue, condition, pathway)
    
    phosphosite_all <- rbind(phosphosite_all, phos_subset)
  }
  
  cat(" ✓\n")
}

cat("\n✓ EXTRACTION COMPLETE\n")
cat("  Total phosphosites:", nrow(phosphosite_all), "\n")
cat("  Unique pathways:", n_distinct(phosphosite_all$pathway), "\n")
cat("  Unique phosphosites:", n_distinct(phosphosite_all$psite_id), "\n\n")

###############################################################
## 4. SAVE MAIN PHOSPHOSITE DATA
###############################################################

write.csv(phosphosite_all, 
          file.path(outdir, "phosphosite_selected_pathways_all_conditions.csv"),
          row.names = FALSE)

cat("✓ Saved: phosphosite_selected_pathways_all_conditions.csv\n\n")

###############################################################
## 5. CREATE PATHWAY SUMMARY TABLE
###############################################################

cat("Step 3: Creating summary statistics...\n")

pathway_summary <- phosphosite_all %>%
  group_by(pathway, condition) %>%
  summarise(
    n_phosphosites = n(),
    mean_logFC = mean(logFC),
    sd_logFC = sd(logFC),
    median_logFC = median(logFC),
    mean_pval = mean(PValue),
    n_sig_pval05 = sum(PValue < 0.05),
    n_strong_logfc = sum(abs(logFC) > 0.5),
    .groups = "drop"
  ) %>%
  arrange(pathway, condition)

write.csv(pathway_summary,
          file.path(outdir, "pathway_summary_by_condition.csv"),
          row.names = FALSE)

cat("✓ Saved: pathway_summary_by_condition.csv\n\n")

###############################################################
## 6. EXTRACT TOP PHOSPHOSITES PER PATHWAY
###############################################################

cat("Step 4: Identifying top phosphosites per pathway...\n")

top_phosphosites <- phosphosite_all %>%
  group_by(pathway) %>%
  arrange(desc(abs(logFC))) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  select(pathway, name, PSite, logFC, PValue, condition)

write.csv(top_phosphosites,
          file.path(outdir, "top_5_phosphosites_per_pathway.csv"),
          row.names = FALSE)

cat("✓ Saved: top_5_phosphosites_per_pathway.csv\n\n")

###############################################################
## 7. CLATHRIN-SPECIFIC DETAILED ANALYSIS
###############################################################

cat("Step 5: CLATHRIN-MEDIATED ENDOCYTOSIS - DETAILED ANALYSIS\n")
cat("-" * 60, "\n\n")

clathrin_pathways <- grep("clathrin", selected_pathways, ignore.case = TRUE, value = TRUE)

if (length(clathrin_pathways) > 0) {
  
  clathrin_all <- phosphosite_all %>%
    filter(pathway %in% clathrin_pathways) %>%
    arrange(PValue)
  
  cat("Clathrin pathway(s):", length(clathrin_pathways), "\n")
  cat("Total phosphosites in Clathrin:", nrow(clathrin_all), "\n")
  cat("Unique genes:", n_distinct(clathrin_all$name), "\n")
  
  # Save clathrin-specific
  write.csv(clathrin_all,
            file.path(outdir, "clathrin_phosphosites_detailed.csv"),
            row.names = FALSE)
  
  cat("\n✓ Saved: clathrin_phosphosites_detailed.csv\n")
  
  # Check for PACS1
  pacs1_present <- "PACS1" %in% clathrin_all$name
  
  if (pacs1_present) {
    cat("\n✅ PACS1 FOUND IN CLATHRIN PATHWAY!\n")
    pacs1_data <- clathrin_all %>%
      filter(name == "PACS1") %>%
      select(name, PSite, logFC, PValue, condition)
    print(pacs1_data)
    
    # Extract logFC pattern
    pacs1_summary <- pacs1_data %>%
      pivot_wider(names_from = condition, values_from = logFC, values_fill = NA)
    cat("\nPACS1 logFC across conditions:\n")
    print(pacs1_summary)
  } else {
    cat("\n⚠️ PACS1 NOT FOUND in Clathrin pathway\n")
    cat("This suggests SELECTIVE pathway regulation!\n")
    
    # List genes that ARE in clathrin
    clathrin_genes <- unique(clathrin_all$name)
    cat("\nGenes in Clathrin pathway:\n")
    print(head(clathrin_genes, 15))
  }
  
  # Clathrin summary by condition
  cat("\n\nClathrin pathway summary by condition:\n")
  clathrin_summary <- clathrin_all %>%
    group_by(condition) %>%
    summarise(
      n_phosphosites = n(),
      mean_logFC = mean(logFC),
      median_logFC = median(logFC),
      mean_pval = mean(PValue),
      n_sig = sum(PValue < 0.05)
    )
  print(clathrin_summary)
  
} else {
  cat("⚠️ No Clathrin pathway found in selected_pathways\n")
}

cat("\n\n")

###############################################################
## 8. CREATE HEATMAPS FOR EACH PATHWAY
###############################################################

cat("Step 6: Creating phosphosite heatmaps...\n")

pdf(file.path(outdir, "heatmaps_all_pathways.pdf"), width = 14, height = 9)

for (pathway in selected_pathways) {
  
  phos_data <- phosphosite_all %>%
    filter(pathway == pathway)
  
  if (nrow(phos_data) < 3) next
  
  # Pivot to wide format (phosphosites × conditions)
  hm_data <- phos_data %>%
    select(psite_id, condition, logFC) %>%
    pivot_wider(
      id_cols = psite_id,
      names_from = condition,
      values_from = logFC,
      values_fill = 0
    ) %>%
    column_to_rownames("psite_id") %>%
    as.matrix()
  
  # Keep top 25 by variance
  if (nrow(hm_data) > 25) {
    row_vars <- apply(hm_data, 1, var)
    hm_data <- hm_data[order(row_vars, decreasing = TRUE)[1:25], ]
  }
  
  # Create heatmap
  pheatmap(
    hm_data,
    main = pathway,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    display_numbers = round(hm_data, 2),
    fontsize_number = 7,
    fontsize_row = 7,
    fontsize = 10,
    breaks = seq(-2, 2, by = 0.1),
    color = colorRampPalette(c("blue", "white", "red"))(40),
    margins = c(12, 8)
  )
}

dev.off()

cat("✓ Saved: heatmaps_all_pathways.pdf\n\n")

###############################################################
## 9. CREATE COMPARISON PLOTS (3-WAY: CXCR7 vs 0s | DMSO vs 0s | CXCR7 vs DMSO)
###############################################################

cat("Step 7: Creating 3-way comparison plots...\n")

if (length(clathrin_pathways) > 0) {
  
  clath_data <- phosphosite_all %>%
    filter(pathway == clathrin_pathways[1]) %>%
    mutate(
      comparison = case_when(
        grepl("dmso.vs.cxcr7", condition, ignore.case = TRUE) ~ "DMSO vs CXCR7",
        grepl("dmso.vs.0", condition, ignore.case = TRUE) ~ "DMSO vs 0s",
        grepl("cxcr7.vs.0", condition, ignore.case = TRUE) ~ "CXCR7 vs 0s"
      ),
      timepoint = case_when(
        grepl("_10", condition) ~ "10s",
        grepl("_600", condition) ~ "600s",
        grepl("_1800", condition) ~ "1800s"
      )
    ) %>%
    drop_na(comparison, timepoint) %>%
    group_by(name, comparison, timepoint) %>%
    summarise(mean_logFC = mean(logFC), .groups = "drop") %>%
    filter(n() > 0)
  
  if (nrow(clath_data) > 0) {
    
    # Get top genes for visualization
    top_genes <- clath_data %>%
      group_by(name) %>%
      summarise(max_abs_logfc = max(abs(mean_logFC))) %>%
      arrange(desc(max_abs_logfc)) %>%
      slice_head(n = 8) %>%
      pull(name)
    
    clath_data_top <- clath_data %>%
      filter(name %in% top_genes)
    
    p <- ggplot(clath_data_top, aes(x = timepoint, y = mean_logFC, 
                                     color = comparison, group = comparison)) +
      geom_point(size = 3, alpha = 0.7) +
      geom_line(alpha = 0.5, size = 1) +
      facet_wrap(~name, scales = "free_y", nrow = 2) +
      theme_minimal() +
      labs(title = "Clathrin Pathway - 3-Way Comparison",
           x = "Timepoint",
           y = "Mean logFC",
           color = "Comparison") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        strip.text = element_text(size = 9),
        plot.title = element_text(face = "bold")
      )
    
    ggsave(file.path(outdir, "comparison_clathrin_3way.pdf"), 
           p, width = 12, height = 8)
    
    cat("✓ Saved: comparison_clathrin_3way.pdf\n")
  }
}

cat("\n")

###############################################################
## 10. CREATE PATHWAY P-VALUE PANEL (PANEL A EQUIVALENT)
###############################################################

cat("Step 8: Creating pathway significance panel...\n")

# Get p-values from enrichment results for selected pathways
pval_data <- data.frame()

for (pathway in selected_pathways) {
  
  for (contrast in names(up_results)) {
    
    if (is.null(up_results[[contrast]])) next
    
    up_df <- up_results[[contrast]]
    down_df <- down_results[[contrast]]
    
    # Get UP result
    up_row <- up_df[up_df$Pathway == pathway, ]
    pval_up <- if (nrow(up_row) > 0) up_row$pval[1] else NA
    log10p_up <- if (!is.na(pval_up)) -log10(pval_up) else 0
    
    # Get DOWN result
    down_row <- down_df[down_df$Pathway == pathway, ]
    pval_down <- if (nrow(down_row) > 0) down_row$pval[1] else NA
    log10p_down <- if (!is.na(pval_down)) -log10(pval_down) else 0
    
    # Combine (positive for UP, negative for DOWN)
    if (log10p_up > log10p_down) {
      log10p_direction <- log10p_up
      direction <- "UP"
    } else {
      log10p_direction <- -log10p_down
      direction <- "DOWN"
    }
    
    pval_data <- rbind(pval_data, data.frame(
      pathway = pathway,
      contrast = contrast,
      log10p = log10p_direction,
      direction = direction,
      stringsAsFactors = FALSE
    ))
  }
}

# Create plot
p_panel_a <- ggplot(pval_data, aes(x = reorder(pathway, log10p), 
                                    y = log10p, 
                                    fill = direction)) +
  geom_col(alpha = 0.8, position = "identity") +
  facet_wrap(~contrast, nrow = 1) +
  coord_flip() +
  scale_fill_manual(values = c("UP" = "darkred", "DOWN" = "darkblue")) +
  labs(title = "Panel A: Selected Pathways - Enrichment Significance",
       x = "Pathway",
       y = "-log10(P-value)",
       fill = "Direction") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 7),
    strip.text = element_text(size = 8),
    panel.grid.major.x = element_line(color = "gray90")
  )

ggsave(file.path(outdir, "Panel_A_pathway_significance.pdf"),
       p_panel_a, width = 16, height = 8)

cat("✓ Saved: Panel_A_pathway_significance.pdf\n\n")

###############################################################
## 11. FINAL SUMMARY
###############################################################

cat("\n", "="*70, "\n")
cat("✅ ANALYSIS COMPLETE!\n")
cat("="*70, "\n\n")

cat("Output files generated in:", outdir, "\n\n")

cat("📊 DATA FILES:\n")
cat("  1. phosphosite_selected_pathways_all_conditions.csv\n")
cat("  2. pathway_summary_by_condition.csv\n")
cat("  3. top_5_phosphosites_per_pathway.csv\n")
cat("  4. clathrin_phosphosites_detailed.csv\n\n")

cat("📈 FIGURES:\n")
cat("  1. Panel_A_pathway_significance.pdf\n")
cat("  2. heatmaps_all_pathways.pdf (Panels B-E equivalent)\n")
cat("  3. comparison_clathrin_3way.pdf (Panel F equivalent)\n\n")

cat("📌 STATISTICS:\n")
cat("  Total phosphosites analyzed:", nrow(phosphosite_all), "\n")
cat("  Total conditions:", n_distinct(phosphosite_all$condition), "\n")
cat("  Selected pathways:", length(selected_pathways), "\n")
if (length(clathrin_pathways) > 0 && nrow(clathrin_all) > 0) {
  cat("  Clathrin phosphosites:", nrow(clathrin_all), "\n")
  if ("PACS1" %in% clathrin_all$name) {
    cat("  ✅ PACS1 found: YES\n")
  } else {
    cat("  ⚠️ PACS1 found: NO (suggests selective regulation)\n")
  }
}

cat("\n✅ Ready for publication!\n\n")

