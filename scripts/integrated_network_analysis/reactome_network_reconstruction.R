###############################################################
## REACTOME-ONLY NETWORK WITH MULTI-PSITE VISUALIZATION
## Handles: Multiple psites per node, multiple timepoints
###############################################################

library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(patchwork)

## ============================================================
## STEP 1: LOAD YOUR PHOSPHO-DATA
## ============================================================

cat("STEP 1: Loading phosphosite data...\n")

# Your existing phospho_wide from Step 14
# Already has: phosphosite_id, name, PSite, logFC_10/600/1800, PValue_10/600/1800

# Aggregate to protein level (for NODE visualization)
protein_summary <- phospho_wide %>%
  group_by(name) %>%
  summarize(
    n_psites = n(),
    
    # Mean logFC per timepoint (for node color)
    mean_logFC_10 = mean(logFC_10, na.rm = TRUE),
    mean_logFC_600 = mean(logFC_600, na.rm = TRUE),
    mean_logFC_1800 = mean(logFC_1800, na.rm = TRUE),
    
    # Max absolute change (for node size scaling)
    max_abs_logFC = max(abs(c(logFC_10, logFC_600, logFC_1800)), na.rm = TRUE),
    
    # Significance
    min_pvalue = min(c(PValue_10, PValue_600, PValue_1800), na.rm = TRUE),
    n_sig_timepoints = sum(c(PValue_10 < 0.05, PValue_600 < 0.05, PValue_1800 < 0.05), na.rm = TRUE),
    
    # Which timepoint has strongest effect?
    dominant_timepoint = case_when(
      abs(mean_logFC_10) >= abs(mean_logFC_600) & abs(mean_logFC_10) >= abs(mean_logFC_1800) ~ "10s",
      abs(mean_logFC_600) >= abs(mean_logFC_10) & abs(mean_logFC_600) >= abs(mean_logFC_1800) ~ "600s",
      TRUE ~ "1800s"
    ),
    
    # Representative sites (for labels)
    top_3_sites = paste(head(PSite[order(min_pvalue)], 3), collapse = ";"),
    
    .groups = "drop"
  )

# Get pathway genes
pathway_genes <- gene_pathway_map %>%
  filter(pathway_num %in% pathway_summary$pathway_num) %>%
  pull(gene) %>%
  unique()

cat(sprintf("✓ Pathway genes: %d\n", length(pathway_genes)))
cat(sprintf("✓ Proteins with phospho-data: %d\n", nrow(protein_summary)))
cat(sprintf("✓ Total phosphosites: %d\n\n", nrow(phospho_wide)))

## ============================================================
## STEP 2: DOWNLOAD REACTOME INTERACTIONS
## ============================================================
options(timeout = 600)
file.remove("reactome_interactions.txt") 

cat("STEP 2: Getting Reactome interactions...\n")

# Download once and cache
if (!file.exists("reactome_human_interactions.rds")) {
  cat("  Downloading Reactome database...\n")
  
  download.file(
    url = "https://reactome.org/download/current/interactors/reactome.all_species.interactions.tab-delimited.txt",
    destfile = "reactome_interactions.txt"
  )
  
  reactome_all <- read_tsv(
    "reactome_interactions.txt",
    skip = 1,
    col_names = c("uniprot_A", "uniprot_B", "gene_A", "gene_B",
                  "interaction_type", "pubmed_ids", "score",
                  "species_A", "species_B"),
    show_col_types = FALSE
  )
  
  reactome_human <- reactome_all %>%
    filter(species_A == "Homo sapiens" & species_B == "Homo sapiens")
  
  saveRDS(reactome_human, "reactome_human_interactions.rds")
  file.remove("reactome_interactions.txt")
}

reactome <- readRDS("reactome_human_interactions.rds")

# Filter to pathway genes
reactome_pathway <- reactome %>%
  filter(gene_A %in% pathway_genes & gene_B %in% pathway_genes) %>%
  select(
    source = gene_A,
    target = gene_B,
    interaction_type,
    pubmed_ids,
    score
  ) %>%
  distinct(source, target, .keep_all = TRUE)

cat(sprintf("✓ Reactome interactions: %d\n\n", nrow(reactome_pathway)))

## ============================================================
## STEP 3: CREATE NODES (WITH REACTOME + NON-REACTOME)
## ============================================================

cat("STEP 3: Creating node table...\n")

# All nodes from Reactome
reactome_nodes <- unique(c(reactome_pathway$source, reactome_pathway$target))

# Proteins with phospho-data but NOT in Reactome
orphan_proteins <- setdiff(protein_summary$name, reactome_nodes)

cat(sprintf("  Reactome nodes: %d\n", length(reactome_nodes)))
cat(sprintf("  Orphan proteins (not in Reactome): %d\n", length(orphan_proteins)))

# Combine
all_nodes <- data.frame(
  name = unique(c(reactome_nodes, protein_summary$name)),
  stringsAsFactors = FALSE
) %>%
  left_join(protein_summary, by = "name") %>%
  mutate(
    in_reactome = name %in% reactome_nodes,
    has_phospho = !is.na(n_psites),
    n_psites = replace_na(n_psites, 0),
    min_pvalue = replace_na(min_pvalue, 1),
    mean_logFC_10 = replace_na(mean_logFC_10, 0),
    mean_logFC_600 = replace_na(mean_logFC_600, 0),
    mean_logFC_1800 = replace_na(mean_logFC_1800, 0)
  )

cat(sprintf("✓ Total nodes: %d\n", nrow(all_nodes)))
cat(sprintf("  - In Reactome: %d\n", sum(all_nodes$in_reactome)))
cat(sprintf("  - With phospho-data: %d\n", sum(all_nodes$has_phospho)))
cat(sprintf("  - Both: %d\n\n", sum(all_nodes$in_reactome & all_nodes$has_phospho)))

## ============================================================
## STEP 4: BUILD NETWORK WITH TRAFFIC SCORES
## ============================================================

cat("STEP 4: Calculating edge importance...\n")

# Build graph
g <- graph_from_data_frame(reactome_pathway, directed = TRUE, vertices = all_nodes)

# Calculate betweenness (centrality)
betweenness_scores <- betweenness(g, directed = TRUE)

# Add traffic scores to edges
edges_scored <- reactome_pathway %>%
  left_join(protein_summary %>% select(name, pval_source = min_pvalue),
            by = c("source" = "name")) %>%
  left_join(protein_summary %>% select(name, pval_target = min_pvalue),
            by = c("target" = "name")) %>%
  mutate(
    # Both proteins significant?
    both_sig = (pval_source < 0.05 & !is.na(pval_source)) & 
      (pval_target < 0.05 & !is.na(pval_target)),
    
    # Significance score
    sig_score = -log10((pval_source + 0.001) * (pval_target + 0.001)),
    
    # Centrality score
    betweenness_source = betweenness_scores[source],
    betweenness_target = betweenness_scores[target],
    centrality_score = (betweenness_source + betweenness_target) / 2,
    
    # Normalize
    sig_norm = scale(sig_score)[,1],
    cent_norm = scale(centrality_score)[,1],
    
    # Combined traffic
    traffic = 0.7 * sig_norm + 0.3 * cent_norm
  ) %>%
  arrange(desc(traffic))

cat("✓ Top 10 edges:\n")
edges_scored %>%
  select(source, target, both_sig, sig_score, traffic) %>%
  head(10) %>%
  print()

cat("\n")

## ============================================================
## STEP 5: EXTRACT CORE NETWORK
## ============================================================

cat("STEP 5: Extracting core network (top 30%)...\n")

traffic_threshold <- quantile(edges_scored$traffic, probs = 0.70, na.rm = TRUE)

core_edges <- edges_scored %>%
  filter(traffic >= traffic_threshold)

core_node_names <- unique(c(core_edges$source, core_edges$target))
core_nodes <- all_nodes %>%
  filter(name %in% core_node_names)

cat(sprintf("✓ Core edges: %d\n", nrow(core_edges)))
cat(sprintf("✓ Core nodes: %d\n\n", nrow(core_nodes)))

g_core <- graph_from_data_frame(
  core_edges %>% select(from = source, to = target, traffic, interaction_type),
  directed = TRUE,
  vertices = core_nodes
)

## ============================================================
## STEP 6: VISUALIZATION - 3 TEMPORAL SNAPSHOTS
## ============================================================

cat("STEP 6: Creating temporal network plots...\n")

# Plot function
plot_network <- function(graph, nodes_df, logfc_col, timepoint_name) {
  
  tbl <- as_tbl_graph(graph) %>%
    activate(nodes) %>%
    left_join(nodes_df %>% select(name, !!sym(logfc_col), min_pvalue, n_psites, in_reactome),
              by = "name")
  
  ggraph(tbl, layout = "stress") +
    # Edges
    geom_edge_link(
      aes(alpha = traffic),
      arrow = arrow(length = unit(2, 'mm')),
      end_cap = circle(3, 'mm'),
      color = "grey50",
      width = 0.5
    ) +
    # Nodes
    geom_node_point(
      aes(size = n_psites,
          fill = !!sym(logfc_col),
          shape = ifelse(in_reactome, 21, 24),  # circle vs triangle
          color = ifelse(min_pvalue < 0.05, "black", "grey80")),
      stroke = 1.2
    ) +
    # Labels (only for significant multi-site proteins)
    geom_node_text(
      aes(label = ifelse(n_psites >= 4 & min_pvalue < 0.05, name, "")),
      repel = TRUE,
      size = 2.5,
      fontface = "bold",
      max.overlaps = 20
    ) +
    scale_fill_gradient2(
      low = "#0571B0", 
      mid = "white", 
      high = "#CA0020",
      midpoint = 0,
      limits = c(-1.5, 1.5),
      name = "mean logFC"
    ) +
    scale_size_continuous(
      range = c(1, 10),
      name = "# psites"
    ) +
    scale_shape_identity() +
    scale_color_identity() +
    scale_edge_alpha_continuous(
      range = c(0.1, 0.8),
      guide = "none"
    ) +
    theme_graph() +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10)
    ) +
    labs(
      title = sprintf("ACKR3/CXCR7 Network: %s", timepoint_name),
      subtitle = sprintf("%d nodes (%d with phospho) | %d edges | Circle=in Reactome, Triangle=orphan",
                         vcount(graph),
                         sum(V(graph)$has_phospho, na.rm = TRUE),
                         ecount(graph))
    )
}

p1 <- plot_network(g_core, core_nodes, "mean_logFC_10", "10 seconds")
p2 <- plot_network(g_core, core_nodes, "mean_logFC_600", "600 seconds")
p3 <- plot_network(g_core, core_nodes, "mean_logFC_1800", "1800 seconds")

pdf("ACKR3_network_temporal.pdf", width = 22, height = 7)
p1 + p2 + p3
dev.off()

cat("✓ Saved: ACKR3_network_temporal.pdf\n\n")

## ============================================================
## STEP 7: SOLVE MULTI-PSITE VISUALIZATION
## ============================================================

cat("STEP 7: Creating phosphosite-level detail plots...\n")

# Strategy: For each hub protein, show ALL phosphosites as small multiples

# Select top hub proteins (≥4 psites, significant)
hub_proteins <- core_nodes %>%
  filter(n_psites >= 4, min_pvalue < 0.05) %>%
  arrange(desc(n_psites)) %>%
  head(20) %>%  # Top 20 hubs
  pull(name)

# Get all psites for these hubs
hub_psites <- phospho_wide %>%
  filter(name %in% hub_proteins) %>%
  select(name, PSite, phosphosite_id,
         logFC_10, logFC_600, logFC_1800,
         PValue_10, PValue_600, PValue_1800) %>%
  pivot_longer(
    cols = c(logFC_10, logFC_600, logFC_1800),
    names_to = "timepoint",
    values_to = "logFC",
    names_prefix = "logFC_"
  ) %>%
  mutate(
    pvalue = case_when(
      timepoint == "10" ~ PValue_10,
      timepoint == "600" ~ PValue_600,
      timepoint == "1800" ~ PValue_1800
    ),
    significant = pvalue < 0.05,
    timepoint = factor(timepoint, levels = c("10", "600", "1800"),
                       labels = c("10s", "600s", "1800s"))
  )

# Plot: Heatmap grid
p_psites <- ggplot(hub_psites, 
                   aes(x = timepoint, y = phosphosite_id, fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_point(
    data = hub_psites %>% filter(significant),
    aes(size = 3),
    shape = 8,
    color = "black"
  ) +
  facet_wrap(~ name, scales = "free_y", ncol = 4) +
  scale_fill_gradient2(
    low = "#0571B0",
    mid = "white",
    high = "#CA0020",
    midpoint = 0,
    name = "logFC"
  ) +
  scale_size_identity() +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 8, angle = 0),
    strip.text = element_text(face = "bold", size = 9),
    panel.spacing = unit(0.5, "lines")
  ) +
  labs(
    title = "Phosphosite-level dynamics for hub proteins",
    subtitle = "Star (*) = p < 0.05 | Each row = one phosphosite",
    x = "Time",
    y = "Phosphosite"
  )

ggsave("ACKR3_phosphosite_heatmap.pdf", p_psites, 
       width = 16, height = 20, limitsize = FALSE)

cat("✓ Saved: ACKR3_phosphosite_heatmap.pdf\n\n")

## ============================================================
## STEP 8: ALTERNATIVE - SUMMARY HEATMAP (COLLAPSED)
## ============================================================

cat("STEP 8: Creating summary heatmap (protein-level)...\n")

# Matrix: Proteins × Timepoints
protein_matrix <- protein_summary %>%
  filter(name %in% hub_proteins) %>%
  select(name, mean_logFC_10, mean_logFC_600, mean_logFC_1800) %>%
  column_to_rownames("name") %>%
  as.matrix()

colnames(protein_matrix) <- c("10s", "600s", "1800s")

# Significance matrix
protein_sig <- protein_summary %>%
  filter(name %in% hub_proteins) %>%
  left_join(
    phospho_wide %>%
      group_by(name) %>%
      summarize(
        sig_10 = sum(PValue_10 < 0.05, na.rm = TRUE),
        sig_600 = sum(PValue_600 < 0.05, na.rm = TRUE),
        sig_1800 = sum(PValue_1800 < 0.05, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "name"
  ) %>%
  select(name, sig_10, sig_600, sig_1800) %>%
  column_to_rownames("name") %>%
  as.matrix()

# Create annotation
library(pheatmap)

pdf("ACKR3_protein_summary_heatmap.pdf", width = 6, height = 10)
pheatmap(
  protein_matrix,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  display_numbers = protein_sig,
  number_color = "black",
  fontsize_number = 8,
  fontsize_row = 8,
  cellwidth = 40,
  cellheight = 15,
  color = colorRampPalette(c("#0571B0", "#F7F7F7", "#CA0020"))(100),
  main = "Hub proteins: Mean logFC (numbers = # sig psites)",
  border_color = "grey60"
)
dev.off()

cat("✓ Saved: ACKR3_protein_summary_heatmap.pdf\n\n")

## ============================================================
## STEP 9: EXPORT DATA
## ============================================================

cat("STEP 9: Exporting network data...\n")

# Nodes with all psite info
nodes_export <- core_nodes %>%
  left_join(
    phospho_wide %>%
      group_by(name) %>%
      summarize(
        all_psites = paste(PSite, collapse = ";"),
        all_psites_ids = paste(phosphosite_id, collapse = ";"),
        .groups = "drop"
      ),
    by = "name"
  )

write_csv(nodes_export, "ACKR3_network_nodes.csv")

# Edges
write_csv(core_edges, "ACKR3_network_edges.csv")

# Orphan proteins (not in Reactome but significant)
orphan_export <- all_nodes %>%
  filter(!in_reactome & has_phospho & min_pvalue < 0.05) %>%
  arrange(min_pvalue)

write_csv(orphan_export, "ACKR3_orphan_proteins.csv")

cat("✓ Saved: ACKR3_network_nodes.csv\n")
cat("✓ Saved: ACKR3_network_edges.csv\n")
cat("✓ Saved: ACKR3_orphan_proteins.csv\n\n")

## ============================================================
## SUMMARY
## ============================================================

cat(strrep("=", 80), "\n")
cat("REACTOME NETWORK RECONSTRUCTION COMPLETE\n")
cat(strrep("=", 80), "\n\n")
cat(sprintf("Total phosphosites: %d\n", nrow(phospho_wide)))
cat(sprintf("Proteins with phospho-data: %d\n", nrow(protein_summary)))
cat(sprintf("Reactome interactions: %d\n", nrow(reactome_pathway)))
cat(sprintf("Core network nodes: %d\n", nrow(core_nodes)))
cat(sprintf("Core network edges: %d\n", nrow(core_edges)))
cat(sprintf("Orphan proteins (sig but not in Reactome): %d\n", nrow(orphan_export)))
cat("\nOutputs:\n")
cat("  1. ACKR3_network_temporal.pdf - 3 timepoint views\n")
cat("  2. ACKR3_phosphosite_heatmap.pdf - All psites for hub proteins\n")
cat("  3. ACKR3_protein_summary_heatmap.pdf - Collapsed protein-level view\n")
cat("  4. ACKR3_network_nodes.csv\n")
cat("  5. ACKR3_network_edges.csv\n")
cat("  6. ACKR3_orphan_proteins.csv\n\n")

cat("Key features:\n")
cat("  ✓ Reactome-only interactions (no STRING)\n")
cat("  ✓ Multiple psites per protein handled via aggregation\n")
cat("  ✓ Multiple timepoints shown in separate panels\n")
cat("  ✓ Orphan proteins identified (sig but not in Reactome)\n")
cat("  ✓ Hub proteins get detailed psite-level heatmap\n\n")