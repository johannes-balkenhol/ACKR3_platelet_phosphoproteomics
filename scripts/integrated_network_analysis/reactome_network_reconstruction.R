###############################################################
## REACTOME NETWORK RECONSTRUCTION WITH PHOSPHO-DATA OVERLAY
## Strategy: Reactome = backbone, phospho-data = evidence layer
###############################################################

library(tidyverse)
library(igraph)
library(ReactomePA)
library(reactome.db)
library(OmnipathR)
library(ggraph)
library(tidygraph)

## ============================================================
## STEP 1: EXTRACT ALL PHOSPHOSITES FROM PATHWAY HEATMAPS
## ============================================================

cat("STEP 1: Extracting phosphosites from pathway heatmaps...\n")

# Read the pathway summary (from your Step 14)
pathway_summary <- read_csv("analysis/pathway_enrichment/reactome/pathway_heatmaps_FINAL/pathway_summary.csv")

# Load UNCOLLAPSED phosphosite data (from your Step 14 code)
val_datasets <- c("val_10.dmso.vs.cxcr7", 
                  "val_600.dmso.vs.cxcr7", 
                  "val_1800.dmso.vs.cxcr7")

all_phosphosite_data <- lapply(val_datasets, function(time_name) {
  time_pt <- str_extract(time_name, "\\d+")
  dataset <- dfs_new_raw[[time_name]]
  
  dataset %>%
    as.data.frame() %>%
    filter(!is.na(name)) %>%
    mutate(
      phosphosite_id = paste(name, PSite, sep = "_"),
      timepoint = time_pt
    ) %>%
    select(phosphosite_id, name, PSite, uniprot_id, logFC, PValue, timepoint)
}) %>% 
  bind_rows()

# Pivot to wide format
phospho_wide <- all_phosphosite_data %>%
  pivot_wider(
    id_cols = c(phosphosite_id, name, PSite, uniprot_id),
    names_from = timepoint,
    values_from = c(logFC, PValue),
    names_sep = "_"
  ) %>%
  mutate(
    mean_logFC = rowMeans(cbind(logFC_10, logFC_600, logFC_1800), na.rm = TRUE),
    max_abs_logFC = pmax(abs(logFC_10), abs(logFC_600), abs(logFC_1800), na.rm = TRUE),
    min_pvalue = pmin(PValue_10, PValue_600, PValue_1800, na.rm = TRUE),
    n_sig = (PValue_10 < 0.05) + (PValue_600 < 0.05) + (PValue_1800 < 0.05)
  )

# Get genes from pathways
pathway_genes <- gene_pathway_map %>%
  filter(pathway_num %in% pathway_summary$pathway_num) %>%
  pull(gene) %>%
  unique()

cat(sprintf("✓ Total genes in top 30 pathways: %d\n", length(pathway_genes)))
cat(sprintf("✓ Total phosphosites: %d\n\n", nrow(phospho_wide)))

## ============================================================
## STEP 2: GET REACTOME INTERACTIONS (BACKBONE)
## ============================================================

cat("STEP 2: Downloading Reactome interactions...\n")

# Get Reactome pathway IDs for your top 30 pathways
# You need to map pathway names to Reactome IDs
top30_pathway_ids <- c(
  "R-HSA-2029481",  # BRAF/RAF1 fusions
  "R-HSA-1474228",  # Degradation of ECM
  "R-HSA-1474244",  # ECM organization
  "R-HSA-165159",   # mTOR signaling
  "R-HSA-216083",   # Integrin cell surface interactions
  # ... add all 30 pathway IDs
)

# Download Reactome interactions via OmniPath
reactome_interactions <- import_post_translational_interactions(
  organisms = 9606,  # Human
  resources = "Reactome"
) %>%
  filter(
    source_genesymbol %in% pathway_genes | 
      target_genesymbol %in% pathway_genes
  ) %>%
  select(
    source = source_genesymbol,
    target = target_genesymbol,
    interaction_type = type,
    is_directed = is_directed,
    is_stimulation = is_stimulation,
    is_inhibition = is_inhibition,
    sources,
    references
  ) %>%
  distinct()

cat(sprintf("✓ Reactome interactions (initial): %d\n", nrow(reactome_interactions)))

# Filter: keep only if at least ONE protein is in your phospho-data
proteins_with_phospho <- unique(phospho_wide$name)

reactome_filtered <- reactome_interactions %>%
  filter(
    source %in% proteins_with_phospho | 
      target %in% proteins_with_phospho
  )

cat(sprintf("✓ Reactome interactions (with phospho evidence): %d\n\n", 
            nrow(reactome_filtered)))

## ============================================================
## STEP 3: AGGREGATE PHOSPHOSITES PER PROTEIN (NODE LEVEL)
## ============================================================

cat("STEP 3: Aggregating phosphosites to protein level...\n")

# Strategy: Per protein, summarize phosphorylation
protein_summary <- phospho_wide %>%
  group_by(name) %>%
  summarize(
    n_psites = n(),
    mean_logFC_10 = mean(logFC_10, na.rm = TRUE),
    mean_logFC_600 = mean(logFC_600, na.rm = TRUE),
    mean_logFC_1800 = mean(logFC_1800, na.rm = TRUE),
    max_abs_logFC_10 = max(abs(logFC_10), na.rm = TRUE),
    max_abs_logFC_600 = max(abs(logFC_600), na.rm = TRUE),
    max_abs_logFC_1800 = max(abs(logFC_1800), na.rm = TRUE),
    min_pvalue = min(min_pvalue, na.rm = TRUE),
    n_sig_sites = sum(n_sig > 0),
    # Dominant direction (more up or down sites?)
    direction = case_when(
      sum(mean_logFC > 0) > sum(mean_logFC < 0) ~ "up",
      sum(mean_logFC < 0) > sum(mean_logFC > 0) ~ "down",
      TRUE ~ "mixed"
    ),
    # Peak timepoint (when is max effect?)
    peak_timepoint = case_when(
      abs(mean_logFC_10) == pmax(abs(mean_logFC_10), abs(mean_logFC_600), abs(mean_logFC_1800)) ~ "10s",
      abs(mean_logFC_600) == pmax(abs(mean_logFC_10), abs(mean_logFC_600), abs(mean_logFC_1800)) ~ "600s",
      abs(mean_logFC_1800) == pmax(abs(mean_logFC_10), abs(mean_logFC_600), abs(mean_logFC_1800)) ~ "1800s",
      TRUE ~ "unclear"
    ),
    # Representative phosphosites (top 3 by significance)
    top_psites = paste(head(PSite[order(min_pvalue)], 3), collapse = ";"),
    .groups = "drop"
  )

cat(sprintf("✓ Proteins with phospho-data: %d\n\n", nrow(protein_summary)))

## ============================================================
## STEP 4: CREATE NETWORK (NODES + EDGES)
## ============================================================

cat("STEP 4: Building network...\n")

# NODES: All proteins in Reactome interactions
all_nodes <- unique(c(reactome_filtered$source, reactome_filtered$target))

# Create node table
nodes <- data.frame(
  name = all_nodes,
  stringsAsFactors = FALSE
) %>%
  left_join(protein_summary, by = "name") %>%
  mutate(
    # If protein not in phospho-data, mark as "no_data"
    has_phospho = !is.na(n_psites),
    n_psites = replace_na(n_psites, 0),
    mean_logFC_1800 = replace_na(mean_logFC_1800, 0),
    min_pvalue = replace_na(min_pvalue, 1),
    direction = replace_na(direction, "no_data")
  )

# EDGES: Reactome interactions
edges <- reactome_filtered %>%
  select(from = source, to = target, 
         interaction_type, is_stimulation, is_inhibition, 
         references)

cat(sprintf("✓ Nodes: %d (%d with phospho-data)\n", 
            nrow(nodes), sum(nodes$has_phospho)))
cat(sprintf("✓ Edges: %d\n\n", nrow(edges)))

## ============================================================
## STEP 5: CALCULATE "TRAFFIC" (EDGE IMPORTANCE)
## ============================================================

cat("STEP 5: Calculating edge traffic/importance...\n")

# Build initial graph
g <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes)

# TRAFFIC METRIC 1: Both source AND target phosphorylated + significant
edges_with_traffic <- edges %>%
  left_join(protein_summary %>% select(name, min_pvalue_source = min_pvalue), 
            by = c("from" = "name")) %>%
  left_join(protein_summary %>% select(name, min_pvalue_target = min_pvalue), 
            by = c("to" = "name")) %>%
  mutate(
    # Both significant?
    both_significant = (min_pvalue_source < 0.05) & (min_pvalue_target < 0.05),
    # Score based on significance
    traffic_score = -log10(min_pvalue_source * min_pvalue_target + 1e-10),
    # Bonus if activation/inhibition is clear
    traffic_bonus = case_when(
      is_stimulation == 1 ~ 1.5,
      is_inhibition == 1 ~ 1.5,
      TRUE ~ 1.0
    ),
    traffic_total = traffic_score * traffic_bonus
  ) %>%
  arrange(desc(traffic_total))

# TRAFFIC METRIC 2: Betweenness centrality (connectivity)
betweenness <- betweenness(g, directed = TRUE)
edges_with_betweenness <- edges_with_traffic %>%
  mutate(
    betweenness_from = betweenness[from],
    betweenness_to = betweenness[to],
    centrality_score = (betweenness_from + betweenness_to) / 2
  )

# COMBINED TRAFFIC SCORE
edges_final <- edges_with_betweenness %>%
  mutate(
    # Normalize both scores
    traffic_norm = scale(traffic_total)[,1],
    centrality_norm = scale(centrality_score)[,1],
    # Combined score
    combined_traffic = 0.7 * traffic_norm + 0.3 * centrality_norm,
    # Rank
    traffic_rank = rank(-combined_traffic)
  )

cat(sprintf("✓ Top 10 highest-traffic edges:\n"))
edges_final %>%
  select(from, to, interaction_type, traffic_total, centrality_score, combined_traffic) %>%
  head(10) %>%
  print()

cat("\n")

## ============================================================
## STEP 6: IDENTIFY CORE NETWORK (HIGHEST TRAFFIC)
## ============================================================

cat("STEP 6: Extracting core high-traffic network...\n")

# Strategy: Keep top X% of edges by traffic
traffic_percentile <- 0.25  # Top 25% of edges

traffic_threshold <- quantile(edges_final$combined_traffic, 
                              probs = 1 - traffic_percentile, 
                              na.rm = TRUE)

core_edges <- edges_final %>%
  filter(combined_traffic >= traffic_threshold)

cat(sprintf("✓ Core edges (top %.0f%%): %d\n", traffic_percentile * 100, nrow(core_edges)))

# Keep only nodes that are in core edges
core_nodes_names <- unique(c(core_edges$from, core_edges$to))
core_nodes <- nodes %>%
  filter(name %in% core_nodes_names)

cat(sprintf("✓ Core nodes: %d\n\n", nrow(core_nodes)))

# Build core graph
g_core <- graph_from_data_frame(
  core_edges %>% select(from, to, interaction_type, combined_traffic),
  directed = TRUE,
  vertices = core_nodes
)

## ============================================================
## STEP 7: VISUALIZATION - TEMPORAL SNAPSHOTS
## ============================================================

cat("STEP 7: Creating temporal network visualizations...\n")

# Function to create network plot for one timepoint
plot_temporal_network <- function(graph, nodes_df, timepoint_col, title) {
  
  tbl <- as_tbl_graph(graph) %>%
    activate(nodes) %>%
    left_join(nodes_df %>% select(name, !!sym(timepoint_col), min_pvalue, n_psites), 
              by = "name")
  
  ggraph(tbl, layout = "stress") +
    # Edges
    geom_edge_link(
      aes(edge_alpha = combined_traffic),
      arrow = arrow(length = unit(2, 'mm')),
      end_cap = circle(3, 'mm'),
      color = "grey50"
    ) +
    # Nodes
    geom_node_point(
      aes(size = n_psites,
          fill = !!sym(timepoint_col),
          color = ifelse(min_pvalue < 0.05, "black", "grey70")),
      shape = 21,
      stroke = 1.5
    ) +
    # Labels (only for top nodes)
    geom_node_text(
      aes(label = ifelse(n_psites >= 3, name, "")),
      repel = TRUE,
      size = 3,
      fontface = "bold"
    ) +
    # Colors
    scale_fill_gradient2(
      low = "#0571B0", mid = "white", high = "#CA0020",
      midpoint = 0,
      name = "mean logFC",
      limits = c(-1.5, 1.5)
    ) +
    scale_size_continuous(
      range = c(2, 12),
      name = "# phosphosites"
    ) +
    scale_edge_alpha_continuous(
      range = c(0.2, 1),
      name = "traffic"
    ) +
    theme_graph() +
    theme(legend.position = "right") +
    labs(title = title,
         subtitle = sprintf("%d nodes, %d edges", 
                            vcount(graph), ecount(graph)))
}

# Create plots
p1 <- plot_temporal_network(g_core, core_nodes, "mean_logFC_10", 
                            "ACKR3/CXCR7 Network: 10 seconds")
p2 <- plot_temporal_network(g_core, core_nodes, "mean_logFC_600", 
                            "ACKR3/CXCR7 Network: 600 seconds")
p3 <- plot_temporal_network(g_core, core_nodes, "mean_logFC_1800", 
                            "ACKR3/CXCR7 Network: 1800 seconds")

# Save
library(patchwork)
pdf("ACKR3_network_temporal.pdf", width = 18, height = 6)
p1 + p2 + p3
dev.off()

cat("✓ Saved: ACKR3_network_temporal.pdf\n\n")

## ============================================================
## STEP 8: IDENTIFY MISSING LINKS
## ============================================================

cat("STEP 8: Identifying missing links...\n")

# STRATEGY: Find proteins that are:
# 1. Significantly phosphorylated
# 2. NOT connected in Reactome network

# All proteins with significant phosphorylation
sig_proteins <- protein_summary %>%
  filter(min_pvalue < 0.05) %>%
  pull(name)

# Which are NOT in core network?
missing_proteins <- sig_proteins[!sig_proteins %in% core_nodes_names]

cat(sprintf("✓ Proteins significantly phosphorylated: %d\n", length(sig_proteins)))
cat(sprintf("✓ Missing from core network: %d\n", length(missing_proteins)))

if (length(missing_proteins) > 0) {
  cat("\nTop 20 missing proteins:\n")
  protein_summary %>%
    filter(name %in% missing_proteins) %>%
    arrange(min_pvalue) %>%
    select(name, n_psites, min_pvalue, mean_logFC_1800, top_psites) %>%
    head(20) %>%
    print()
}

cat("\n")

## ============================================================
## STEP 9: SOLVE MULTIPLE PHOSPHOSITES PER NODE
## ============================================================

cat("STEP 9: Creating detailed node-level visualization...\n")

# Create a "bubble plot" for key nodes showing all phosphosites

key_nodes <- core_nodes %>%
  filter(n_psites >= 3, min_pvalue < 0.05) %>%
  pull(name)

phosphosite_details <- phospho_wide %>%
  filter(name %in% key_nodes) %>%
  select(name, PSite, phosphosite_id, 
         logFC_10, logFC_600, logFC_1800,
         PValue_10, PValue_600, PValue_1800) %>%
  pivot_longer(
    cols = starts_with("logFC"),
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
    timepoint = factor(timepoint, levels = c("10", "600", "1800"))
  )

# Plot: Heatmap of phosphosites for top nodes
library(ggplot2)

p_detail <- ggplot(phosphosite_details %>% filter(name %in% head(key_nodes, 10)),
                   aes(x = timepoint, y = phosphosite_id, fill = logFC)) +
  geom_tile(color = "white", size = 0.5) +
  geom_point(aes(size = ifelse(significant, 3, NA)), 
             shape = 8, color = "black") +
  facet_wrap(~ name, scales = "free_y", ncol = 2) +
  scale_fill_gradient2(
    low = "#0571B0", mid = "white", high = "#CA0020",
    midpoint = 0,
    name = "logFC"
  ) +
  scale_size_identity() +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 10),
    strip.text = element_text(face = "bold", size = 11)
  ) +
  labs(
    title = "Phosphosite-level detail for key hub proteins",
    subtitle = "Star (*) indicates p < 0.05",
    x = "Time (seconds)",
    y = "Phosphosite"
  )

ggsave("ACKR3_phosphosite_details.pdf", p_detail, width = 12, height = 16)
cat("✓ Saved: ACKR3_phosphosite_details.pdf\n\n")

## ============================================================
## STEP 10: EXPORT NETWORK WITH FULL ANNOTATIONS
## ============================================================

cat("STEP 10: Exporting annotated network...\n")

# Export node table with all phosphosite info
nodes_export <- core_nodes %>%
  left_join(
    phospho_wide %>%
      group_by(name) %>%
      summarize(
        all_psites = paste(PSite, collapse = ";"),
        all_psites_logFC_1800 = paste(round(logFC_1800, 2), collapse = ";"),
        all_psites_pvalue = paste(format(PValue_1800, scientific = TRUE, digits = 2), collapse = ";"),
        .groups = "drop"
      ),
    by = "name"
  )

write_csv(nodes_export, "ACKR3_network_nodes.csv")

# Export edge table with traffic scores
edges_export <- core_edges %>%
  select(from, to, interaction_type, 
         is_stimulation, is_inhibition,
         traffic_score, centrality_score, combined_traffic,
         references)

write_csv(edges_export, "ACKR3_network_edges.csv")

cat("✓ Saved: ACKR3_network_nodes.csv\n")
cat("✓ Saved: ACKR3_network_edges.csv\n\n")

## ============================================================
## STEP 11: CREATE SIMPLIFIED SUMMARY NETWORK
## ============================================================

cat("STEP 11: Creating simplified summary network...\n")

# Super-condensed: Only show proteins with ≥5 phosphosites
super_core_nodes <- core_nodes %>%
  filter(n_psites >= 5)

super_core_edges <- core_edges %>%
  filter(from %in% super_core_nodes$name & to %in% super_core_nodes$name)

g_super <- graph_from_data_frame(
  super_core_edges,
  directed = TRUE,
  vertices = super_core_nodes
)

p_super <- plot_temporal_network(
  g_super, super_core_nodes, "mean_logFC_1800",
  "ACKR3/CXCR7 Core Network (≥5 phosphosites per protein)"
)

ggsave("ACKR3_network_CORE.pdf", p_super, width = 12, height = 10)
cat("✓ Saved: ACKR3_network_CORE.pdf\n\n")

## ============================================================
## SUMMARY
## ============================================================

cat(strrep("=", 80), "\n")
cat("NETWORK RECONSTRUCTION COMPLETE\n")
cat(strrep("=", 80), "\n\n")
cat(sprintf("Total phosphosites: %d\n", nrow(phospho_wide)))
cat(sprintf("Total proteins with phospho-data: %d\n", nrow(protein_summary)))
cat(sprintf("Reactome interactions (full): %d\n", nrow(reactome_filtered)))
cat(sprintf("High-traffic edges: %d\n", nrow(core_edges)))
cat(sprintf("Core network nodes: %d\n", nrow(core_nodes)))
cat(sprintf("Missing proteins (significant but not connected): %d\n", 
            length(missing_proteins)))
cat("\nOutput files:\n")
cat("  - ACKR3_network_temporal.pdf (3 timepoints)\n")
cat("  - ACKR3_network_CORE.pdf (simplified)\n")
cat("  - ACKR3_phosphosite_details.pdf (multi-site nodes)\n")
cat("  - ACKR3_network_nodes.csv (node annotations)\n")
cat("  - ACKR3_network_edges.csv (edge annotations)\n\n")