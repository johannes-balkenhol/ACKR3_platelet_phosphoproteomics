###############################################################
## REACTOME NETWORK RECONSTRUCTION
## Phosphosite-level temporal signaling network
##
## PREREQUISITES: Run pl.phospho.7 (Steps 1-15) first!
## This script uses: phospho_wide_filtered, gene_pathway_map, 
##   top_pathways_master, dfs_new_raw
###############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(patchwork)
  library(pheatmap)
  library(org.Hs.eg.db)
  library(annotate)
})

###############################################################
## STEP 0: BRIDGE - Connect to enrichment pipeline data
###############################################################

cat(strrep("=", 80), "\n")
cat("STEP 0: Connecting to enrichment pipeline data\n")
cat(strrep("=", 80), "\n\n")

# Verify required objects exist from pl.phospho.7
required_objects <- c("phospho_wide_filtered", "gene_pathway_map", 
                      "top_pathways_master", "dfs_new_raw")

for (obj in required_objects) {
  if (!exists(obj)) {
    stop(sprintf("Missing object: %s\nRun pl.phospho.7 (Steps 1-15) first!", obj))
  }
}

# Use uncollapsed phosphosite data from Step 14B
phospho_wide <- phospho_wide_filtered

# Get gene names from top 30 pathways
pathway_genes_vec <- gene_pathway_map %>%
  pull(gene) %>%
  unique()

cat(sprintf("✓ phospho_wide: %d phosphosites (uncollapsed)\n", nrow(phospho_wide)))
cat(sprintf("✓ pathway genes: %d unique genes from top %d pathways\n", 
            length(pathway_genes_vec), n_distinct(gene_pathway_map$pathway_num)))
cat(sprintf("✓ gene_pathway_map: %d gene-pathway links\n\n", nrow(gene_pathway_map)))

###############################################################
## STEP 1: AGGREGATE TO PROTEIN LEVEL (for nodes)
###############################################################

cat(strrep("=", 80), "\n")
cat("STEP 1: Aggregating phosphosite data to protein level\n")
cat(strrep("=", 80), "\n\n")

protein_summary <- phospho_wide %>%
  group_by(name) %>%
  summarize(
    n_psites = n(),
    uniprot_id = first(uniprot_id),
    
    # Mean logFC per timepoint
    mean_logFC_10 = mean(logFC_10, na.rm = TRUE),
    mean_logFC_600 = mean(logFC_600, na.rm = TRUE),
    mean_logFC_1800 = mean(logFC_1800, na.rm = TRUE),
    
    # Max absolute change
    max_abs_logFC = max(abs(c(logFC_10, logFC_600, logFC_1800)), na.rm = TRUE),
    
    # Significance
    min_pvalue = min(c(PValue_10, PValue_600, PValue_1800), na.rm = TRUE),
    n_sig_timepoints = sum(c(
      any(PValue_10 < 0.05, na.rm = TRUE),
      any(PValue_600 < 0.05, na.rm = TRUE),
      any(PValue_1800 < 0.05, na.rm = TRUE)
    )),
    
    # Dominant timepoint
    dominant_timepoint = case_when(
      abs(mean(logFC_10, na.rm = TRUE)) >= abs(mean(logFC_600, na.rm = TRUE)) & 
        abs(mean(logFC_10, na.rm = TRUE)) >= abs(mean(logFC_1800, na.rm = TRUE)) ~ "10s",
      abs(mean(logFC_600, na.rm = TRUE)) >= abs(mean(logFC_10, na.rm = TRUE)) & 
        abs(mean(logFC_600, na.rm = TRUE)) >= abs(mean(logFC_1800, na.rm = TRUE)) ~ "600s",
      TRUE ~ "1800s"
    ),
    
    # Top 3 sites for labels
    top_3_sites = paste(head(PSite, 3), collapse = ";"),
    
    # Pathway memberships
    pathway_nums = ifelse(
      name[1] %in% gene_pathway_map$gene,
      paste(sort(unique(gene_pathway_map$pathway_num[gene_pathway_map$gene == name[1]])), collapse = ","),
      ""
    ),
    
    .groups = "drop"
  )

cat(sprintf("✓ %d proteins with phospho-data\n", nrow(protein_summary)))
cat(sprintf("✓ %d proteins in top 30 pathways\n\n", 
            sum(protein_summary$name %in% pathway_genes_vec)))

###############################################################
## STEP 2: DOWNLOAD REACTOME HUMAN INTERACTIONS
###############################################################

cat(strrep("=", 80), "\n")
cat("STEP 2: Getting Reactome human interactions\n")
cat(strrep("=", 80), "\n\n")

options(timeout = 600)

# Reactome tab-delimited format (human only):
# Col 1: UniProt_A
# Col 2: Ensembl_A (pipe-delimited)
# Col 3: Entrez_A (pipe-delimited)
# Col 4: UniProt_B
# Col 5: Ensembl_B (pipe-delimited)
# Col 6: Entrez_B (pipe-delimited)
# Col 7: interaction_type (direct_complex, indirect_complex, reaction, neighbouring_reaction)
# Col 8: context (Reactome instance ID)
# Col 9: PubMed IDs (comma-delimited)

if (!file.exists("reactome_human_interactions.rds")) {
  
  cat("  Downloading Reactome human interactions...\n")
  
  download.file(
    url = "https://reactome.org/download/current/interactors/reactome.homo_sapiens.interactions.tab-delimited.txt",
    destfile = "reactome_human_raw.txt",
    mode = "wb"
  )
  
  # Read, skipping first line (comment with # count)
  raw_lines <- readLines("reactome_human_raw.txt")
  cat(sprintf("  Raw lines: %d\n", length(raw_lines)))
  cat(sprintf("  Header: %s\n", substr(raw_lines[1], 1, 80)))
  
  # Skip comment line(s) starting with #
  data_start <- min(which(!grepl("^#", raw_lines)))
  cat(sprintf("  Data starts at line: %d\n", data_start))
  
  reactome_human <- read_tsv(
    "reactome_human_raw.txt",
    skip = data_start - 1,
    col_names = c("uniprot_A", "ensembl_A", "entrez_A",
                  "uniprot_B", "ensembl_B", "entrez_B",
                  "interaction_type", "context", "pubmed_ids"),
    show_col_types = FALSE
  )
  
  cat(sprintf("  Parsed rows: %d\n", nrow(reactome_human)))
  cat(sprintf("  Interaction types:\n"))
  print(table(reactome_human$interaction_type))
  
  saveRDS(reactome_human, "reactome_human_interactions.rds")
  file.remove("reactome_human_raw.txt")
  cat("  ✓ Saved: reactome_human_interactions.rds\n\n")
  
} else {
  reactome_human <- readRDS("reactome_human_interactions.rds")
  
  # Check if it's the old empty file
  if (nrow(reactome_human) == 0) {
    cat("  ⚠️ Old empty file detected, re-downloading...\n")
    file.remove("reactome_human_interactions.rds")
    stop("Deleted empty RDS. Please re-run this step.")
  }
}

cat(sprintf("✓ Reactome human interactions loaded: %d rows\n\n", nrow(reactome_human)))

###############################################################
## STEP 3: MAP REACTOME UNIPROT IDs TO GENE SYMBOLS
###############################################################

cat(strrep("=", 80), "\n")
cat("STEP 3: Mapping Reactome UniProt IDs to gene symbols\n")
cat(strrep("=", 80), "\n\n")

# Build UniProt → Gene symbol map from our phospho data
uniprot_to_gene <- protein_summary %>%
  select(uniprot_id, name) %>%
  distinct() %>%
  filter(!is.na(uniprot_id), !is.na(name))

cat(sprintf("  UniProt→Gene map from phospho data: %d entries\n", nrow(uniprot_to_gene)))

# Also use org.Hs.eg.db for broader mapping
# Get all UniProt IDs from Reactome that involve our pathway genes
all_reactome_uniprots <- unique(c(reactome_human$uniprot_A, reactome_human$uniprot_B))
cat(sprintf("  Total unique UniProt IDs in Reactome: %d\n", length(all_reactome_uniprots)))

# Map via org.Hs.eg.db: UniProt → Entrez → Symbol
tryCatch({
  uniprot_entrez <- AnnotationDbi::select(org.Hs.eg.db, 
                                           keys = all_reactome_uniprots,
                                           keytype = "UNIPROT",
                                           columns = c("SYMBOL", "ENTREZID"))
  
  uniprot_to_symbol_db <- uniprot_entrez %>%
    filter(!is.na(SYMBOL)) %>%
    mutate(SYMBOL = toupper(SYMBOL)) %>%
    select(uniprot_id = UNIPROT, name = SYMBOL) %>%
    distinct()
  
  cat(sprintf("  org.Hs.eg.db mapping: %d UniProt → Symbol\n", nrow(uniprot_to_symbol_db)))
  
  # Combine: prefer our phospho data mapping, supplement with org.Hs.eg.db
  uniprot_to_gene <- bind_rows(
    uniprot_to_gene,
    uniprot_to_symbol_db %>% 
      filter(!uniprot_id %in% uniprot_to_gene$uniprot_id)
  ) %>%
    distinct(uniprot_id, .keep_all = TRUE)
  
}, error = function(e) {
  cat("  ⚠️ org.Hs.eg.db mapping failed, using phospho data only\n")
})

cat(sprintf("  ✓ Combined map: %d UniProt → Gene entries\n\n", nrow(uniprot_to_gene)))

# Add gene symbols to Reactome interactions
reactome_mapped <- reactome_human %>%
  left_join(uniprot_to_gene %>% rename(gene_A = name), 
            by = c("uniprot_A" = "uniprot_id")) %>%
  left_join(uniprot_to_gene %>% rename(gene_B = name), 
            by = c("uniprot_B" = "uniprot_id")) %>%
  filter(!is.na(gene_A), !is.na(gene_B)) %>%
  mutate(gene_A = toupper(gene_A), gene_B = toupper(gene_B))

cat(sprintf("✓ Reactome interactions with gene symbols: %d\n", nrow(reactome_mapped)))

###############################################################
## STEP 4: FILTER TO PATHWAY GENES & BUILD EDGE TABLE
###############################################################

cat(strrep("=", 80), "\n")
cat("STEP 4: Filtering to pathway genes\n")
cat(strrep("=", 80), "\n\n")

# Filter: both genes must be in our pathway gene list
reactome_pathway <- reactome_mapped %>%
  filter(gene_A %in% pathway_genes_vec & gene_B %in% pathway_genes_vec) %>%
  filter(gene_A != gene_B) %>%  # Remove self-loops
  select(source = gene_A, target = gene_B, interaction_type, pubmed_ids) %>%
  distinct(source, target, .keep_all = TRUE)

cat(sprintf("✓ Pathway interactions: %d edges\n", nrow(reactome_pathway)))
cat(sprintf("  Interaction types:\n"))
print(table(reactome_pathway$interaction_type))

# Also include ALL proteins with phospho data for broader network
reactome_phospho <- reactome_mapped %>%
  filter(gene_A %in% protein_summary$name & gene_B %in% protein_summary$name) %>%
  filter(gene_A != gene_B) %>%
  select(source = gene_A, target = gene_B, interaction_type, pubmed_ids) %>%
  distinct(source, target, .keep_all = TRUE)

cat(sprintf("\n✓ All phospho-protein interactions: %d edges\n\n", nrow(reactome_phospho)))

###############################################################
## STEP 5: BUILD NETWORK WITH SIGNAL FLOW SCORES
###############################################################

cat(strrep("=", 80), "\n")
cat("STEP 5: Building network with signal flow scores\n")
cat(strrep("=", 80), "\n\n")

# Use pathway-filtered edges for core network
edges_for_network <- reactome_pathway

# Build node table
node_names <- unique(c(edges_for_network$source, edges_for_network$target))

nodes <- data.frame(name = node_names, stringsAsFactors = FALSE) %>%
  left_join(protein_summary, by = "name") %>%
  mutate(
    has_phospho = !is.na(n_psites),
    n_psites = replace_na(n_psites, 0),
    min_pvalue = replace_na(min_pvalue, 1),
    mean_logFC_10 = replace_na(mean_logFC_10, 0),
    mean_logFC_600 = replace_na(mean_logFC_600, 0),
    mean_logFC_1800 = replace_na(mean_logFC_1800, 0),
    max_abs_logFC = replace_na(max_abs_logFC, 0)
  )

cat(sprintf("✓ Network nodes: %d\n", nrow(nodes)))
cat(sprintf("  - With phospho data: %d\n", sum(nodes$has_phospho)))
cat(sprintf("  - Without phospho data: %d\n", sum(!nodes$has_phospho)))

# Build igraph
g <- graph_from_data_frame(
  edges_for_network %>% select(from = source, to = target, interaction_type),
  directed = FALSE,
  vertices = nodes
)

cat(sprintf("✓ Graph: %d nodes, %d edges\n", vcount(g), ecount(g)))

# Calculate centrality metrics
V(g)$degree <- degree(g)
V(g)$betweenness <- betweenness(g, directed = FALSE)

cat(sprintf("  Top 10 by degree:\n"))
deg_order <- order(V(g)$degree, decreasing = TRUE)
for (i in 1:min(10, length(deg_order))) {
  v <- deg_order[i]
  cat(sprintf("    %s: degree=%d, betweenness=%.0f\n", 
              V(g)$name[v], V(g)$degree[v], V(g)$betweenness[v]))
}

###############################################################
## STEP 6: COMPUTE TEMPORAL EDGE WEIGHTS (SIGNAL FLOW)
###############################################################

cat("\n", strrep("=", 80), "\n")
cat("STEP 6: Computing temporal edge weights\n")
cat(strrep("=", 80), "\n\n")

# For each edge, compute a "signal flow" score at each timepoint
# = mean of -log10(pvalue) * |logFC| for source and target proteins

compute_edge_flow <- function(edges_df, protein_data, timepoint_suffix) {
  
  logfc_col <- paste0("mean_logFC_", timepoint_suffix)
  
  edges_df %>%
    left_join(protein_data %>% 
                select(name, min_pvalue, !!sym(logfc_col)) %>%
                rename(source_pval = min_pvalue, source_logfc = !!sym(logfc_col)),
              by = c("source" = "name")) %>%
    left_join(protein_data %>% 
                select(name, min_pvalue, !!sym(logfc_col)) %>%
                rename(target_pval = min_pvalue, target_logfc = !!sym(logfc_col)),
              by = c("target" = "name")) %>%
    mutate(
      source_pval = replace_na(source_pval, 1),
      target_pval = replace_na(target_pval, 1),
      source_logfc = replace_na(source_logfc, 0),
      target_logfc = replace_na(target_logfc, 0),
      
      # Signal intensity: both nodes active
      source_signal = -log10(source_pval + 1e-10) * abs(source_logfc),
      target_signal = -log10(target_pval + 1e-10) * abs(target_logfc),
      
      # Edge flow = geometric mean of endpoint signals
      edge_flow = sqrt(pmax(source_signal * target_signal, 0)),
      
      # Both significant?
      both_sig = (source_pval < 0.05) & (target_pval < 0.05),
      
      # Concordant direction? (both up or both down)
      concordant = sign(source_logfc) == sign(target_logfc)
    ) %>%
    select(source, target, edge_flow, both_sig, concordant, 
           source_logfc, target_logfc, source_pval, target_pval)
}

edge_flow_10 <- compute_edge_flow(edges_for_network, protein_summary, "10")
edge_flow_600 <- compute_edge_flow(edges_for_network, protein_summary, "600")
edge_flow_1800 <- compute_edge_flow(edges_for_network, protein_summary, "1800")

cat("Edge flow statistics:\n")
for (tp in list(list("10s", edge_flow_10), list("600s", edge_flow_600), list("1800s", edge_flow_1800))) {
  ef <- tp[[2]]
  cat(sprintf("  %s: mean=%.2f, max=%.2f, edges with both sig: %d/%d\n",
              tp[[1]], mean(ef$edge_flow, na.rm = TRUE), max(ef$edge_flow, na.rm = TRUE),
              sum(ef$both_sig, na.rm = TRUE), nrow(ef)))
}

###############################################################
## STEP 7: FIND TOP SIGNAL ROUTES (SHORTEST PATHS)
###############################################################

cat("\n", strrep("=", 80), "\n")
cat("STEP 7: Finding major signal flow routes\n")
cat(strrep("=", 80), "\n\n")

find_signal_routes <- function(graph, edge_flow_df, protein_data, 
                               timepoint_name, n_routes = 15) {
  
  cat(sprintf("\n--- %s ---\n", timepoint_name))
  
  # Create weighted graph: invert flow for shortest path (high flow = short distance)
  max_flow <- max(edge_flow_df$edge_flow, na.rm = TRUE) + 1
  
  # Build edge list with inverted weights
  el <- edge_flow_df %>%
    mutate(weight = max_flow - edge_flow) %>%
    filter(!is.na(weight), is.finite(weight))
  
  g_weighted <- graph_from_data_frame(
    el %>% select(from = source, to = target, weight),
    directed = FALSE,
    vertices = data.frame(name = unique(c(el$source, el$target)))
  )
  
  # Find hub proteins (high flow endpoints)
  node_flow <- edge_flow_df %>%
    pivot_longer(cols = c(source, target), values_to = "node") %>%
    group_by(node) %>%
    summarize(total_flow = sum(edge_flow, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total_flow))
  
  # Get top hub pairs
  top_hubs <- head(node_flow$node, 10)
  
  cat(sprintf("  Top 10 hubs by flow: %s\n", paste(top_hubs, collapse = ", ")))
  
  # Find shortest paths between all hub pairs
  routes <- list()
  route_scores <- numeric()
  
  hub_pairs <- combn(top_hubs, 2, simplify = FALSE)
  
  for (pair in hub_pairs) {
    src <- pair[1]
    tgt <- pair[2]
    
    if (src %in% V(g_weighted)$name & tgt %in% V(g_weighted)$name) {
      sp <- tryCatch(
        shortest_paths(g_weighted, from = src, to = tgt, output = "vpath", weights = E(g_weighted)$weight),
        error = function(e) NULL
      )
      
      if (!is.null(sp) && length(sp$vpath[[1]]) >= 3) {
        path_nodes <- V(g_weighted)$name[sp$vpath[[1]]]
        
        # Score = sum of edge flows along path
        path_score <- 0
        for (j in 1:(length(path_nodes) - 1)) {
          ef <- edge_flow_df %>%
            filter((source == path_nodes[j] & target == path_nodes[j + 1]) |
                     (source == path_nodes[j + 1] & target == path_nodes[j]))
          if (nrow(ef) > 0) path_score <- path_score + ef$edge_flow[1]
        }
        
        route_key <- paste(path_nodes, collapse = " → ")
        routes[[route_key]] <- path_nodes
        route_scores[route_key] <- path_score
      }
    }
  }
  
  # Sort by score, take top N
  if (length(routes) > 0) {
    top_idx <- order(route_scores, decreasing = TRUE)[1:min(n_routes, length(routes))]
    
    cat(sprintf("\n  Top %d signal routes:\n", length(top_idx)))
    cat(strrep("-", 80), "\n")
    
    top_routes <- list()
    for (i in seq_along(top_idx)) {
      key <- names(routes)[top_idx[i]]
      score <- route_scores[top_idx[i]]
      path <- routes[[top_idx[i]]]
      
      cat(sprintf("  #%2d (score=%.1f): %s\n", i, score, key))
      top_routes[[i]] <- list(nodes = path, score = score)
    }
    
    return(top_routes)
  } else {
    cat("  ⚠️ No routes found\n")
    return(list())
  }
}

routes_10 <- find_signal_routes(g, edge_flow_10, protein_summary, "10 seconds")
routes_600 <- find_signal_routes(g, edge_flow_600, protein_summary, "600 seconds")
routes_1800 <- find_signal_routes(g, edge_flow_1800, protein_summary, "1800 seconds")

###############################################################
## STEP 8: VISUALIZATION - 3 TEMPORAL NETWORK SNAPSHOTS
###############################################################

cat("\n", strrep("=", 80), "\n")
cat("STEP 8: Creating temporal network plots\n")
cat(strrep("=", 80), "\n\n")

# Extract top route edges for highlighting
get_route_edges <- function(routes, n_top = 5) {
  route_edges <- data.frame(source = character(), target = character())
  for (r in routes[1:min(n_top, length(routes))]) {
    nodes <- r$nodes
    for (j in 1:(length(nodes) - 1)) {
      route_edges <- bind_rows(route_edges, 
                                data.frame(source = nodes[j], target = nodes[j + 1]))
    }
  }
  distinct(route_edges)
}

route_edges_10 <- get_route_edges(routes_10)
route_edges_600 <- get_route_edges(routes_600)
route_edges_1800 <- get_route_edges(routes_1800)

# Compute layout ONCE (shared across timepoints)
set.seed(42)
layout_coords <- create_layout(as_tbl_graph(g), layout = "stress")

plot_temporal_network <- function(graph, nodes_df, edge_flow_df, route_edges,
                                  logfc_col, timepoint_name, layout) {
  
  # Build tbl_graph with temporal data
  tbl <- as_tbl_graph(graph) %>%
    activate(nodes) %>%
    left_join(nodes_df %>% select(name, !!sym(logfc_col), min_pvalue, n_psites, has_phospho),
              by = "name")
  
  # Mark route edges
  edge_df <- tbl %>% activate(edges) %>% as_tibble()
  node_names_map <- tbl %>% activate(nodes) %>% as_tibble() %>% pull(name)
  
  edge_df <- edge_df %>%
    mutate(
      source_name = node_names_map[from],
      target_name = node_names_map[to],
      is_route = paste(source_name, target_name) %in% paste(route_edges$source, route_edges$target) |
        paste(target_name, source_name) %in% paste(route_edges$source, route_edges$target)
    )
  
  # Merge edge flow
  edge_df <- edge_df %>%
    left_join(edge_flow_df %>% select(source, target, edge_flow, both_sig),
              by = c("source_name" = "source", "target_name" = "target"))
  
  tbl <- tbl %>%
    activate(edges) %>%
    mutate(
      is_route = edge_df$is_route,
      edge_flow = replace_na(edge_df$edge_flow, 0),
      both_sig = replace_na(edge_df$both_sig, FALSE)
    )
  
  p <- ggraph(tbl, layout = "stress") +
    # Regular edges (grey, thin)
    geom_edge_link(
      aes(alpha = edge_flow, 
          color = ifelse(is_route, "route", "normal"),
          width = ifelse(is_route, 1.5, 0.3)),
      show.legend = FALSE
    ) +
    scale_edge_color_manual(values = c("normal" = "grey70", "route" = "#FF4500")) +
    scale_edge_width_identity() +
    scale_edge_alpha_continuous(range = c(0.05, 0.6)) +
    
    # Nodes
    geom_node_point(
      aes(size = pmax(n_psites, 1),
          fill = !!sym(logfc_col)),
      shape = 21,
      color = "grey30",
      stroke = 0.5
    ) +
    
    # Labels (significant proteins only)
    geom_node_text(
      aes(label = ifelse(min_pvalue < 0.05 & n_psites >= 2, name, "")),
      repel = TRUE,
      size = 2.2,
      fontface = "bold",
      max.overlaps = 25,
      segment.color = "grey50",
      segment.size = 0.2
    ) +
    
    scale_fill_gradient2(
      low = "#0571B0", mid = "#F7F7F7", high = "#CA0020",
      midpoint = 0, limits = c(-1.5, 1.5), oob = scales::squish,
      name = "mean logFC"
    ) +
    scale_size_continuous(range = c(1.5, 8), name = "# psites") +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40")
    ) +
    labs(
      title = timepoint_name,
      subtitle = sprintf("%d nodes | %d edges | red = top signal routes",
                         vcount(graph), ecount(graph))
    )
  
  return(p)
}

p1 <- plot_temporal_network(g, nodes, edge_flow_10, route_edges_10,
                             "mean_logFC_10", "10 seconds", layout_coords)
p2 <- plot_temporal_network(g, nodes, edge_flow_600, route_edges_600,
                             "mean_logFC_600", "600 seconds", layout_coords)
p3 <- plot_temporal_network(g, nodes, edge_flow_1800, route_edges_1800,
                             "mean_logFC_1800", "1800 seconds", layout_coords)

# Save combined plot
pdf("ACKR3_network_temporal_signalflow.pdf", width = 24, height = 9)
print(p1 + p2 + p3 + plot_layout(guides = "collect") & 
        theme(legend.position = "bottom"))
dev.off()
cat("✓ Saved: ACKR3_network_temporal_signalflow.pdf\n")

# Save individual plots
for (tp in list(list("10s", p1), list("600s", p2), list("1800s", p3))) {
  fname <- sprintf("ACKR3_network_%s.pdf", tp[[1]])
  pdf(fname, width = 10, height = 10)
  print(tp[[2]])
  dev.off()
  cat(sprintf("✓ Saved: %s\n", fname))
}

###############################################################
## STEP 9: PHOSPHOSITE-LEVEL HEATMAP FOR HUB PROTEINS
###############################################################

cat("\n", strrep("=", 80), "\n")
cat("STEP 9: Hub protein phosphosite heatmaps\n")
cat(strrep("=", 80), "\n\n")

# Select top hub proteins (≥3 psites, significant)
hub_proteins <- nodes %>%
  filter(has_phospho, n_psites >= 3, min_pvalue < 0.05) %>%
  arrange(min_pvalue) %>%
  head(25) %>%
  pull(name)

cat(sprintf("Hub proteins for detailed view: %d\n", length(hub_proteins)))

# Get all psites for hubs
hub_psites <- phospho_wide %>%
  filter(name %in% hub_proteins) %>%
  select(name, PSite, phosphosite_id,
         logFC_10, logFC_600, logFC_1800,
         PValue_10, PValue_600, PValue_1800) %>%
  pivot_longer(
    cols = starts_with("logFC_"),
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

p_psites <- ggplot(hub_psites, 
                   aes(x = timepoint, y = phosphosite_id, fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_point(
    data = hub_psites %>% filter(significant),
    shape = 8, color = "black", size = 1
  ) +
  facet_wrap(~ name, scales = "free_y", ncol = 5) +
  scale_fill_gradient2(
    low = "#0571B0", mid = "#F7F7F7", high = "#CA0020",
    midpoint = 0, name = "logFC"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 4),
    axis.text.x = element_text(size = 7, angle = 0),
    strip.text = element_text(face = "bold", size = 8),
    panel.spacing = unit(0.3, "lines")
  ) +
  labs(
    title = "Phosphosite-level dynamics for hub proteins",
    subtitle = "Star = p < 0.05 | Each row = one phosphosite",
    x = "Time", y = NULL
  )

ggsave("ACKR3_phosphosite_heatmap_hubs.pdf", p_psites,
       width = 18, height = 22, limitsize = FALSE)
cat("✓ Saved: ACKR3_phosphosite_heatmap_hubs.pdf\n")

###############################################################
## STEP 10: SIGNAL ROUTE SUMMARY TABLE
###############################################################

cat("\n", strrep("=", 80), "\n")
cat("STEP 10: Exporting signal route summaries\n")
cat(strrep("=", 80), "\n\n")

# Combine routes into table
export_routes <- function(routes, timepoint_name) {
  if (length(routes) == 0) return(data.frame())
  
  bind_rows(lapply(seq_along(routes), function(i) {
    data.frame(
      timepoint = timepoint_name,
      rank = i,
      score = routes[[i]]$score,
      path_length = length(routes[[i]]$nodes),
      route = paste(routes[[i]]$nodes, collapse = " → "),
      stringsAsFactors = FALSE
    )
  }))
}

all_routes <- bind_rows(
  export_routes(routes_10, "10s"),
  export_routes(routes_600, "600s"),
  export_routes(routes_1800, "1800s")
)

write.csv(all_routes, "ACKR3_signal_routes.csv", row.names = FALSE)
cat("✓ Saved: ACKR3_signal_routes.csv\n")

# Export network data
write.csv(nodes, "ACKR3_network_nodes.csv", row.names = FALSE)
write.csv(edges_for_network, "ACKR3_network_edges.csv", row.names = FALSE)
cat("✓ Saved: ACKR3_network_nodes.csv\n")
cat("✓ Saved: ACKR3_network_edges.csv\n")

# Orphan proteins (in pathway, significant, but not in Reactome network)
orphan_proteins <- protein_summary %>%
  filter(name %in% pathway_genes_vec) %>%
  filter(!name %in% node_names) %>%
  filter(min_pvalue < 0.05) %>%
  arrange(min_pvalue)

write.csv(orphan_proteins, "ACKR3_orphan_proteins.csv", row.names = FALSE)
cat(sprintf("✓ Saved: ACKR3_orphan_proteins.csv (%d proteins)\n", nrow(orphan_proteins)))

###############################################################
## SUMMARY
###############################################################

cat("\n", strrep("█", 80), "\n")
cat("NETWORK RECONSTRUCTION COMPLETE\n")
cat(strrep("█", 80), "\n\n")

cat(sprintf("Network: %d nodes, %d edges\n", vcount(g), ecount(g)))
cat(sprintf("  - With phospho data: %d / %d\n", sum(nodes$has_phospho), nrow(nodes)))
cat(sprintf("  - Hub proteins (≥3 psites, sig): %d\n", length(hub_proteins)))
cat(sprintf("\nSignal routes found:\n"))
cat(sprintf("  10s:   %d routes\n", length(routes_10)))
cat(sprintf("  600s:  %d routes\n", length(routes_600)))
cat(sprintf("  1800s: %d routes\n", length(routes_1800)))
cat(sprintf("\nOrphan proteins: %d (sig but not in Reactome)\n", nrow(orphan_proteins)))

cat("\n📁 Output files:\n")
cat("  ACKR3_network_temporal_signalflow.pdf  (3 timepoints side-by-side)\n")
cat("  ACKR3_network_10s/600s/1800s.pdf       (individual timepoints)\n")
cat("  ACKR3_phosphosite_heatmap_hubs.pdf     (phosphosite-level detail)\n")
cat("  ACKR3_signal_routes.csv                (top routes per timepoint)\n")
cat("  ACKR3_network_nodes.csv\n")
cat("  ACKR3_network_edges.csv\n")
cat("  ACKR3_orphan_proteins.csv\n\n")
