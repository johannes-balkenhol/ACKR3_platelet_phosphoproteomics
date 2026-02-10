#############################
## FLEXIBLE QUERY SYSTEM FOR ENRICHMENT RESULTS (FULLY FIXED!)
#############################

cat("\n", strrep("=", 80), "\n")
cat("ENRICHMENT RESULTS QUERY SYSTEM\n")
cat(strrep("=", 80), "\n\n")

#############################
## PARSING FUNCTIONS (ALL VECTORIZED!)
#############################

# Extract batch (val or init) - VECTORIZED
extract_batch <- function(time_str) {
  sapply(time_str, function(x) {
    if (grepl("^val_", x)) return("VAL")
    if (grepl("^init_", x)) return("INIT")
    return(NA)
  }, USE.NAMES = FALSE)
}

# Extract timepoint (10, 30, 60, 300, 600, 900, 1800 seconds) - VECTORIZED
extract_timepoint <- function(time_str) {
  sapply(time_str, function(x) {
    match <- regexpr("[0-9]+(?=\\.|$)", x, perl = TRUE)
    if (match > 0) {
      return(as.numeric(regmatches(x, match)))
    }
    return(NA)
  }, USE.NAMES = FALSE)
}

# Extract comparison type - VECTORIZED
extract_comparison <- function(time_str) {
  sapply(time_str, function(x) {
    if (grepl("dmso\\.vs\\.cxcr7", x, ignore.case = TRUE)) return("CXCR7_vs_DMSO")
    if (grepl("cxcr7\\.vs\\.0s", x, ignore.case = TRUE)) return("CXCR7_vs_0s")
    if (grepl("dmso\\.vs\\.0s", x, ignore.case = TRUE)) return("DMSO_vs_0s")
    return(NA)
  }, USE.NAMES = FALSE)
}

# Add parsed columns to data
annotate_results <- function(df) {
  df %>%
    mutate(
      batch = extract_batch(time),
      timepoint_sec = extract_timepoint(time),
      comparison = extract_comparison(time),
      significant = pvalue < 0.05
    )
}

# Annotate both datasets
cat("Annotating UP-regulated results...\n")
up_anno <- annotate_results(up_summary_long_all)

cat("Annotating DOWN-regulated results...\n")
down_anno <- annotate_results(down_summary_long_all)

cat("✓ Data annotated with batch, timepoint, and comparison type\n\n")

# Quick check
cat("UP annotation sample:\n")
print(head(up_anno %>% select(time, batch, timepoint_sec, comparison)))

cat("\n")

#############################
## SUMMARY FUNCTION
#############################

summarize_enrichment <- function(
    data = NULL,  # up_anno or down_anno
    batch = NULL,  # "VAL", "INIT", or NULL for both
    timepoint = NULL,  # 10, 600, 1800, or NULL for all
    comparison = NULL,  # "CXCR7_vs_DMSO", "CXCR7_vs_0s", "DMSO_vs_0s", or NULL
    min_pvalue = 1,  # Show trends down to this p-value (e.g., 0.1)
    title = NULL
) {
  
  if (is.null(data)) {
    stop("Must provide data (up_anno or down_anno)")
  }
  
  # Filter data
  result <- data
  
  if (!is.null(batch)) {
    result <- result %>% filter(batch %in% batch)
  }
  
  if (!is.null(timepoint)) {
    result <- result %>% filter(timepoint_sec %in% timepoint)
  }
  
  if (!is.null(comparison)) {
    result <- result %>% filter(comparison %in% comparison)
  }
  
  result <- result %>% 
    filter(pvalue <= min_pvalue) %>%
    distinct(pathway, time, .keep_all = TRUE) %>%
    arrange(pvalue)
  
  # Print title
  if (!is.null(title)) {
    cat("\n", strrep("=", 80), "\n")
    cat(title, "\n")
    cat(strrep("=", 80), "\n\n")
  }
  
  # Print summary stats
  sig_count <- sum(result$pvalue < 0.05, na.rm = TRUE)
  trend_count <- sum(result$pvalue >= 0.05 & result$pvalue <= min_pvalue, na.rm = TRUE)
  
  cat(sprintf("Filters: Batch=%s | Timepoint=%s | Comparison=%s | p≤%.3f\n",
              paste(batch %||% "ALL", collapse=","),
              paste(timepoint %||% "ALL", collapse=","),
              paste(comparison %||% "ALL", collapse=","),
              min_pvalue))
  cat(sprintf("Results: %d significant (p<0.05) + %d trends (0.05≤p≤%0.3f) = %d total\n\n",
              sig_count, trend_count, min_pvalue, nrow(result)))
  
  # Print results
  if (nrow(result) > 0) {
    print(result %>% select(pathway, pvalue, ratio, batch, timepoint_sec, comparison))
  } else {
    cat("(No pathways match criteria)\n")
  }
  
  invisible(result)
}

#############################
## QUICK ANALYSES
#############################

cat("\n\n", strrep("=", 80), "\n")
cat("AUTOMATIC ANALYSES (10 key queries)\n")
cat(strrep("=", 80), "\n")

# ============================================================================
# 1. CXCR7 vs DMSO - THE GOLD STANDARD (Pure CXCR7 effect)
# ============================================================================

cat("\n\n1️⃣  CXCR7 vs DMSO (PURE CXCR7 EFFECT - USE THIS!)\n")

result_1 <- summarize_enrichment(
  data = up_anno,
  comparison = "CXCR7_vs_DMSO",
  min_pvalue = 0.05,
  title = "UP-regulated in CXCR7 vs DMSO (Pure CXCR7 activation)"
)

result_1_down <- summarize_enrichment(
  data = down_anno,
  comparison = "CXCR7_vs_DMSO",
  min_pvalue = 0.05,
  title = "DOWN-regulated in CXCR7 vs DMSO (Pure CXCR7 inhibition)"
)

# ============================================================================
# 2. Vehicle control (DMSO vs 0s)
# ============================================================================

cat("\n\n2️⃣  DMSO vs 0s (VEHICLE CONTROL - should be quiet!)\n")

result_2 <- summarize_enrichment(
  data = up_anno,
  comparison = "DMSO_vs_0s",
  min_pvalue = 0.1,
  title = "UP-regulated in DMSO vs 0s (Vehicle effects only)"
)

result_2_down <- summarize_enrichment(
  data = down_anno,
  comparison = "DMSO_vs_0s",
  min_pvalue = 0.1,
  title = "DOWN-regulated in DMSO vs 0s (Vehicle effects only)"
)

if ((nrow(result_2) + nrow(result_2_down)) < 5) {
  cat("\n✅ EXCELLENT: Clean vehicle control!\n")
} else {
  cat("\n⚠️  WARNING: Vehicle causes some pathway changes\n")
}

# ============================================================================
# 3. Temporal dynamics at 600s
# ============================================================================

cat("\n\n3️⃣  PEAK RESPONSE at 600s (CXCR7 vs DMSO)\n")

result_3 <- summarize_enrichment(
  data = up_anno,
  timepoint = 600,
  comparison = "CXCR7_vs_DMSO",
  min_pvalue = 0.05,
  title = "UP-regulated at 600s (CXCR7 vs DMSO)"
)

result_3_down <- summarize_enrichment(
  data = down_anno,
  timepoint = 600,
  comparison = "CXCR7_vs_DMSO",
  min_pvalue = 0.05,
  title = "DOWN-regulated at 600s (CXCR7 vs DMSO)"
)

# ============================================================================
# 4. Batch comparison at 600s
# ============================================================================

cat("\n\n4️⃣  BATCH REPRODUCIBILITY at 600s\n")

cat("\nVAL batch:\n")
result_4a <- summarize_enrichment(
  data = up_anno,
  batch = "VAL",
  timepoint = 600,
  comparison = "CXCR7_vs_DMSO",
  min_pvalue = 0.05,
  title = "VAL: 600s CXCR7 vs DMSO (UP)"
)

cat("\nINIT batch:\n")
result_4b <- summarize_enrichment(
  data = up_anno,
  batch = "INIT",
  timepoint = 600,
  comparison = "CXCR7_vs_DMSO",
  min_pvalue = 0.05,
  title = "INIT: 600s CXCR7 vs DMSO (UP)"
)

# ============================================================================
# 5. Consensus pathways
# ============================================================================

cat("\n\n5️⃣  CONSENSUS UP PATHWAYS (appear in both VAL and INIT)\n")

val_pathways <- up_anno %>%
  filter(batch == "VAL", comparison == "CXCR7_vs_DMSO", pvalue < 0.05) %>%
  distinct(pathway) %>%
  pull(pathway)

init_pathways <- up_anno %>%
  filter(batch == "INIT", comparison == "CXCR7_vs_DMSO", pvalue < 0.05) %>%
  distinct(pathway) %>%
  pull(pathway)

consensus <- intersect(val_pathways, init_pathways)

cat(sprintf("VAL batch: %d pathways\n", length(val_pathways)))
cat(sprintf("INIT batch: %d pathways\n", length(init_pathways)))
cat(sprintf("Consensus (both): %d pathways ✅\n", length(consensus)))

if (length(consensus) > 0) {
  cat("\nConsensus UP pathways:\n")
  consensus_data <- up_anno %>%
    filter(pathway %in% consensus, comparison == "CXCR7_vs_DMSO") %>%
    arrange(pvalue) %>%
    distinct(pathway, .keep_all = TRUE)
  
  print(consensus_data %>% select(pathway, pvalue, ratio))
}

# ============================================================================
# 6. Consensus DOWN pathways
# ============================================================================

cat("\n\n6️⃣  CONSENSUS DOWN PATHWAYS (appear in both VAL and INIT)\n")

val_down_paths <- down_anno %>%
  filter(batch == "VAL", comparison == "CXCR7_vs_DMSO", pvalue < 0.05) %>%
  distinct(pathway) %>%
  pull(pathway)

init_down_paths <- down_anno %>%
  filter(batch == "INIT", comparison == "CXCR7_vs_DMSO", pvalue < 0.05) %>%
  distinct(pathway) %>%
  pull(pathway)

consensus_down <- intersect(val_down_paths, init_down_paths)

cat(sprintf("VAL batch: %d pathways\n", length(val_down_paths)))
cat(sprintf("INIT batch: %d pathways\n", length(init_down_paths)))
cat(sprintf("Consensus (both): %d pathways ✅\n", length(consensus_down)))

if (length(consensus_down) > 0) {
  cat("\nConsensus DOWN pathways:\n")
  consensus_down_data <- down_anno %>%
    filter(pathway %in% consensus_down, comparison == "CXCR7_vs_DMSO") %>%
    arrange(pvalue) %>%
    distinct(pathway, .keep_all = TRUE)
  
  print(consensus_down_data %>% select(pathway, pvalue, ratio))
}

# ============================================================================
# 7. Overall summary table
# ============================================================================

cat("\n\n7️⃣  COMPLETE OVERVIEW (all conditions)\n")

overview <- bind_rows(
  up_anno %>% 
    filter(pvalue < 0.05) %>%
    group_by(batch, timepoint_sec, comparison) %>%
    summarize(
      n_pathways = n_distinct(pathway),
      min_pvalue = min(pvalue),
      direction = "UP",
      .groups = "drop"
    ),
  down_anno %>%
    filter(pvalue < 0.05) %>%
    group_by(batch, timepoint_sec, comparison) %>%
    summarize(
      n_pathways = n_distinct(pathway),
      min_pvalue = min(pvalue),
      direction = "DOWN",
      .groups = "drop"
    )
) %>%
  arrange(batch, comparison, timepoint_sec, desc(direction))

cat("Significant pathways (p<0.05) by all dimensions:\n\n")
print(overview)

# ============================================================================
# 8. Timepoint comparison
# ============================================================================

cat("\n\n8️⃣  TEMPORAL DYNAMICS: Pathways per timepoint (CXCR7 vs DMSO)\n")

temporal_summary <- up_anno %>%
  filter(comparison == "CXCR7_vs_DMSO", pvalue < 0.05) %>%
  group_by(timepoint_sec) %>%
  summarize(
    n_pathways = n_distinct(pathway),
    min_pvalue = min(pvalue),
    .groups = "drop"
  ) %>%
  arrange(timepoint_sec)

cat("UP-regulated pathways by timepoint:\n")
print(temporal_summary)

# ============================================================================
# 9. Key pathways summary
# ============================================================================

cat("\n\n9️⃣  KEY PATHWAYS: MAPK, RHO, PKA, Hemostasis\n")

key_pathways <- c(
  "GRB2:SOS provides linkage to MAPK signaling for Integrins",
  "Signalling to ERKs",
  "Signaling by Rho GTPases",
  "RHOA GTPase cycle",
  "PKA activation",
  "Hemostasis"
)

key_results <- bind_rows(
  up_anno %>% filter(pathway %in% key_pathways),
  down_anno %>% filter(pathway %in% key_pathways)
) %>%
  arrange(pathway, batch, timepoint_sec, comparison) %>%
  filter(pvalue < 0.1)

cat("Key pathway behavior (p<0.1):\n\n")
if (nrow(key_results) > 0) {
  print(key_results %>% 
          select(pathway, pvalue, batch, timepoint_sec, comparison) %>%
          arrange(pathway, pvalue))
} else {
  cat("(No key pathways in results)\n")
}

# ============================================================================
# 10. Final summary
# ============================================================================

cat("\n\n🔟  FINAL SUMMARY\n")
cat(strrep("=", 80), "\n\n")

total_up_cxcr7dmso <- nrow(up_anno %>% 
                             filter(comparison == "CXCR7_vs_DMSO", pvalue < 0.05) %>%
                             distinct(pathway))

total_down_cxcr7dmso <- nrow(down_anno %>% 
                               filter(comparison == "CXCR7_vs_DMSO", pvalue < 0.05) %>%
                               distinct(pathway))

cat(sprintf("CXCR7 vs DMSO (PURE SIGNAL):\n"))
cat(sprintf("  ✅ %d unique UP-regulated pathways\n", total_up_cxcr7dmso))
cat(sprintf("  ✅ %d unique DOWN-regulated pathways\n", total_down_cxcr7dmso))
cat(sprintf("  ✅ %d UP consensus (both VAL & INIT)\n", length(consensus)))
cat(sprintf("  ✅ %d DOWN consensus (both VAL & INIT)\n\n", length(consensus_down)))

cat("Vehicle control quality:\n")
if ((nrow(result_2) + nrow(result_2_down)) < 5) {
  cat("  ✅ EXCELLENT - DMSO has minimal effects\n\n")
} else {
  cat("  ⚠️  Vehicle causes some effects\n\n")
}

cat("✓ QUERY SYSTEM READY!\n")
cat("  Use: summarize_enrichment(up_anno, comparison='CXCR7_vs_DMSO', ...)\n")
cat("  See: QUERY_SYSTEM_GUIDE.md for more examples\n\n")








#############################
## SYSTEMATIC PATHWAY OVERVIEW
## Complete data: all significant pathways p<0.05
## Detailed information per pathway
## VAL and INIT completely separated
#############################

cat("\n", strrep("█", 140), "\n")
cat("COMPREHENSIVE PATHWAY ENRICHMENT ANALYSIS\n")
cat("All significant pathways (p<0.05) with detailed information\n")
cat(strrep("█", 140), "\n")

# Helper function to format pathway output
print_pathway_table <- function(data, direction_name, title) {
  if (nrow(data) == 0) {
    cat(sprintf("%s: (No significant pathways)\n\n", direction_name))
    return()
  }
  
  cat(sprintf("%s (%d pathways):\n", direction_name, nrow(data)))
  cat("┌" %+% strrep("─", 138) %+% "┐\n")
  cat(sprintf("│ %-3s │ %-65s │ p-value  │ ratio  │\n", "Rank", "Pathway Name"))
  cat("├" %+% strrep("─", 138) %+% "┤\n")
  
  for (i in seq_len(nrow(data))) {
    pathway_name <- data$pathway[i]
    if (nchar(pathway_name) > 65) {
      pathway_name <- substr(pathway_name, 1, 62) %+% "..."
    }
    cat(sprintf("│ %3d │ %-65s │ %.5f  │ %.4f  │\n",
                i,
                pathway_name,
                data$pvalue[i],
                data$ratio[i]))
  }
  
  cat("└" %+% strrep("─", 138) %+% "┘\n\n")
}

#############################
## VAL BATCH
#############################

cat("\n")
cat(strrep("▓▓▓", 47), "\n")
cat("BATCH: VAL (Validation cohort)\n")
cat(strrep("▓▓▓", 47), "\n\n")

# Get all comparisons and timepoints for VAL
val_combos <- up_anno %>%
  filter(batch == "VAL") %>%
  distinct(comparison, timepoint_sec) %>%
  arrange(comparison, timepoint_sec)

for (row_idx in seq_len(nrow(val_combos))) {
  comp <- val_combos$comparison[row_idx]
  tp <- val_combos$timepoint_sec[row_idx]
  
  cat("\n")
  cat(strrep("═", 140), "\n")
  cat(sprintf("COMPARISON: %s | TIMEPOINT: %d seconds\n", comp, tp))
  cat(strrep("═", 140), "\n\n")
  
  # UP-regulated
  up_data <- up_anno %>%
    filter(batch == "VAL", comparison == comp, timepoint_sec == tp, pvalue < 0.05) %>%
    distinct(pathway, .keep_all = TRUE) %>%
    arrange(pvalue) %>%
    select(pathway, pvalue, ratio)
  
  print_pathway_table(up_data, "UP-REGULATED", "")
  
  # DOWN-regulated
  down_data <- down_anno %>%
    filter(batch == "VAL", comparison == comp, timepoint_sec == tp, pvalue < 0.05) %>%
    distinct(pathway, .keep_all = TRUE) %>%
    arrange(pvalue) %>%
    select(pathway, pvalue, ratio)
  
  print_pathway_table(down_data, "DOWN-REGULATED", "")
  
  # Summary
  total_sig <- nrow(up_data) + nrow(down_data)
  cat(sprintf("TOTAL: %d pathways | UP: %d | DOWN: %d\n", total_sig, nrow(up_data), nrow(down_data)))
  cat("\n")
}

#############################
## INIT BATCH
#############################

cat("\n\n")
cat(strrep("▓▓▓", 47), "\n")
cat("BATCH: INIT (Initial cohort)\n")
cat(strrep("▓▓▓", 47), "\n\n")

# Get all comparisons and timepoints for INIT
init_combos <- up_anno %>%
  filter(batch == "INIT") %>%
  distinct(comparison, timepoint_sec) %>%
  arrange(comparison, timepoint_sec)

if (nrow(init_combos) == 0) {
  cat("\n(No INIT data available)\n")
} else {
  for (row_idx in seq_len(nrow(init_combos))) {
    comp <- init_combos$comparison[row_idx]
    tp <- init_combos$timepoint_sec[row_idx]
    
    cat("\n")
    cat(strrep("═", 140), "\n")
    cat(sprintf("COMPARISON: %s | TIMEPOINT: %d seconds\n", comp, tp))
    cat(strrep("═", 140), "\n\n")
    
    # UP-regulated
    up_data <- up_anno %>%
      filter(batch == "INIT", comparison == comp, timepoint_sec == tp, pvalue < 0.05) %>%
      distinct(pathway, .keep_all = TRUE) %>%
      arrange(pvalue) %>%
      select(pathway, pvalue, ratio)
    
    print_pathway_table(up_data, "UP-REGULATED", "")
    
    # DOWN-regulated
    down_data <- down_anno %>%
      filter(batch == "INIT", comparison == comp, timepoint_sec == tp, pvalue < 0.05) %>%
      distinct(pathway, .keep_all = TRUE) %>%
      arrange(pvalue) %>%
      select(pathway, pvalue, ratio)
    
    print_pathway_table(down_data, "DOWN-REGULATED", "")
    
    # Summary
    total_sig <- nrow(up_data) + nrow(down_data)
    cat(sprintf("TOTAL: %d pathways | UP: %d | DOWN: %d\n", total_sig, nrow(up_data), nrow(down_data)))
    cat("\n")
  }
}

#############################
## SUMMARY STATISTICS
#############################

cat("\n\n")
cat(strrep("█", 140), "\n")
cat("SUMMARY STATISTICS\n")
cat(strrep("█", 140), "\n\n")

# Count by comparison, batch, timepoint
up_stats <- up_anno %>%
  filter(pvalue < 0.05) %>%
  group_by(comparison, batch, timepoint_sec) %>%
  summarize(
    up_pathways = n_distinct(pathway),
    up_results = n(),
    up_min_p = min(pvalue),
    up_max_p = max(pvalue),
    .groups = "drop"
  )

down_stats <- down_anno %>%
  filter(pvalue < 0.05) %>%
  group_by(comparison, batch, timepoint_sec) %>%
  summarize(
    down_pathways = n_distinct(pathway),
    down_results = n(),
    down_min_p = min(pvalue),
    down_max_p = max(pvalue),
    .groups = "drop"
  )

summary_stats <- full_join(up_stats, down_stats, 
                           by = c("comparison", "batch", "timepoint_sec")) %>%
  replace(is.na(.), 0) %>%
  mutate(
    total_pathways = up_pathways + down_pathways,
    total_results = up_results + down_results
  ) %>%
  arrange(comparison, batch, timepoint_sec)

cat("All Condition Combinations:\n")
cat("┌" %+% strrep("─", 138) %+% "┐\n")
cat(sprintf("│ %-20s │ %-6s │ %-10s │ %-4s ↑ │ %-4s ↓ │ %-4s │ p-val range UP | p-val range DOWN │\n",
            "Comparison", "Batch", "Timepoint", "Up", "Down", "Total"))
cat("├" %+% strrep("─", 138) %+% "┤\n")

for (i in seq_len(nrow(summary_stats))) {
  comp <- summary_stats$comparison[i]
  batch <- summary_stats$batch[i]
  tp <- summary_stats$timepoint_sec[i]
  up_p <- summary_stats$up_pathways[i]
  down_p <- summary_stats$down_pathways[i]
  total_p <- summary_stats$total_pathways[i]
  
  up_range <- if (up_p > 0) {
    sprintf("%.4f-%.4f", summary_stats$up_min_p[i], summary_stats$up_max_p[i])
  } else {
    "─────────────"
  }
  
  down_range <- if (down_p > 0) {
    sprintf("%.4f-%.4f", summary_stats$down_min_p[i], summary_stats$down_max_p[i])
  } else {
    "───────────────"
  }
  
  cat(sprintf("│ %-20s │ %-6s │ %5ds     │ %4d │ %4d │ %4d │ %s │ %s │\n",
              comp, batch, tp, up_p, down_p, total_p, up_range, down_range))
}

cat("└" %+% strrep("─", 138) %+% "┘\n\n")

# Top 5 most significant pathways overall
cat("\nTop 10 Most Significant Pathways (p<0.05) Across All Conditions:\n")
cat("┌" %+% strrep("─", 138) %+% "┐\n")
cat(sprintf("│ %-3s │ Direction │ %-20s │ Batch │ Timepoint │ %-45s │ p-value │\n",
            "Rank", "Comparison", "Pathway"))
cat("├" %+% strrep("─", 138) %+% "┤\n")

top_pathways <- bind_rows(
  up_anno %>%
    filter(pvalue < 0.05) %>%
    select(comparison, batch, timepoint_sec, pathway, pvalue) %>%
    mutate(direction = "UP"),
  down_anno %>%
    filter(pvalue < 0.05) %>%
    select(comparison, batch, timepoint_sec, pathway, pvalue) %>%
    mutate(direction = "DOWN")
) %>%
  arrange(pvalue) %>%
  head(10)

for (i in seq_len(nrow(top_pathways))) {
  comp <- top_pathways$comparison[i]
  batch <- top_pathways$batch[i]
  tp <- top_pathways$timepoint_sec[i]
  pathway <- top_pathways$pathway[i]
  pval <- top_pathways$pvalue[i]
  dir <- top_pathways$direction[i]
  
  if (nchar(pathway) > 45) {
    pathway <- substr(pathway, 1, 42) %+% "..."
  }
  
  cat(sprintf("│ %3d │ %-9s │ %-20s │ %-5s │ %5ds     │ %-45s │ %.5f │\n",
              i, dir, comp, batch, tp, pathway, pval))
}

cat("└" %+% strrep("─", 138) %+% "┘\n\n")

cat(strrep("█", 140), "\n")
cat("✓ COMPREHENSIVE ANALYSIS COMPLETE\n")
cat(strrep("█", 140), "\n\n")











































#############################
## COMPREHENSIVE PATHWAY OVERVIEW
## Including suggestive and borderline pathways
## Multiple p-value thresholds for complete picture
#############################

cat("\n", strrep("█", 150), "\n")
cat("EXTENDED PATHWAY ANALYSIS: ALL SIGNIFICANT & BORDERLINE PATHWAYS\n")
cat("Showing p<0.05 (significant), p<0.10 (suggestive), p<0.20 (borderline)\n")
cat(strrep("█", 150), "\n\n")

#############################
## Helper function for extended table
#############################

print_extended_pathway_table <- function(data, direction_name, p_thresholds = c(0.05, 0.10, 0.20)) {
  if (nrow(data) == 0) {
    cat(sprintf("%s: (No pathways at any threshold)\n\n", direction_name))
    return()
  }
  
  # Arrange by p-value
  data <- data %>% arrange(pvalue)
  
  cat(sprintf("%s:\n", direction_name))
  cat("┌" %+% strrep("─", 148) %+% "┐\n")
  cat(sprintf("│ %-3s │ %-60s │ p-value  │ ratio  │ Sig   │ Threshold       │\n", "Rank", "Pathway Name"))
  cat("├" %+% strrep("─", 148) %+% "┤\n")
  
  for (i in seq_len(nrow(data))) {
    pathway_name <- data$pathway[i]
    if (nchar(pathway_name) > 60) {
      pathway_name <- substr(pathway_name, 1, 57) %+% "..."
    }
    
    p_val <- data$pvalue[i]
    
    # Determine significance level
    if (p_val < 0.05) {
      threshold <- "p<0.05 ⭐⭐⭐"
    } else if (p_val < 0.10) {
      threshold <- "p<0.10 ⭐⭐"
    } else if (p_val < 0.20) {
      threshold <- "p<0.20 ⭐"
    } else {
      threshold <- "p≥0.20"
    }
    
    sig_marker <- if (p_val < 0.05) "YES" else if (p_val < 0.10) "?" else "weak"
    
    cat(sprintf("│ %3d │ %-60s │ %.5f  │ %.4f  │ %-5s │ %-15s │\n",
                i,
                pathway_name,
                p_val,
                data$ratio[i],
                sig_marker,
                threshold))
  }
  
  cat("└" %+% strrep("─", 148) %+% "┘\n\n")
}

#############################
## VAL BATCH: CXCR7_vs_DMSO (PURE EFFECT)
#############################

cat(strrep("▓", 50), "\n")
cat("VAL BATCH: CXCR7_vs_DMSO (PURE CXCR7 EFFECT)\n")
cat("Gold standard - vehicle subtracted\n")
cat(strrep("▓", 50), "\n\n")

# Get all comparisons and timepoints for VAL CXCR7_vs_DMSO
val_cxcr7dmso_combos <- up_anno %>%
  filter(batch == "VAL", comparison == "CXCR7_vs_DMSO") %>%
  distinct(timepoint_sec) %>%
  arrange(timepoint_sec) %>%
  pull(timepoint_sec)

for (tp in val_cxcr7dmso_combos) {
  cat("\n")
  cat(strrep("═", 150), "\n")
  cat(sprintf("TIMEPOINT: %d seconds\n", tp))
  cat(strrep("═", 150), "\n\n")
  
  # UP-regulated - all pathways (p < 0.20)
  up_data_all <- up_anno %>%
    filter(batch == "VAL", comparison == "CXCR7_vs_DMSO", timepoint_sec == tp, pvalue < 0.20) %>%
    distinct(pathway, .keep_all = TRUE) %>%
    arrange(pvalue) %>%
    select(pathway, pvalue, ratio)
  
  print_extended_pathway_table(up_data_all, "UP-REGULATED")
  
  # DOWN-regulated - all pathways (p < 0.20)
  down_data_all <- down_anno %>%
    filter(batch == "VAL", comparison == "CXCR7_vs_DMSO", timepoint_sec == tp, pvalue < 0.20) %>%
    distinct(pathway, .keep_all = TRUE) %>%
    arrange(pvalue) %>%
    select(pathway, pvalue, ratio)
  
  print_extended_pathway_table(down_data_all, "DOWN-REGULATED")
  
  # Summary
  n_up_sig <- nrow(up_anno %>% filter(batch == "VAL", comparison == "CXCR7_vs_DMSO", timepoint_sec == tp, pvalue < 0.05))
  n_up_suggest <- nrow(up_anno %>% filter(batch == "VAL", comparison == "CXCR7_vs_DMSO", timepoint_sec == tp, pvalue < 0.10, pvalue >= 0.05))
  n_up_borderline <- nrow(up_anno %>% filter(batch == "VAL", comparison == "CXCR7_vs_DMSO", timepoint_sec == tp, pvalue < 0.20, pvalue >= 0.10))
  
  n_down_sig <- nrow(down_anno %>% filter(batch == "VAL", comparison == "CXCR7_vs_DMSO", timepoint_sec == tp, pvalue < 0.05))
  n_down_suggest <- nrow(down_anno %>% filter(batch == "VAL", comparison == "CXCR7_vs_DMSO", timepoint_sec == tp, pvalue < 0.10, pvalue >= 0.05))
  n_down_borderline <- nrow(down_anno %>% filter(batch == "VAL", comparison == "CXCR7_vs_DMSO", timepoint_sec == tp, pvalue < 0.20, pvalue >= 0.10))
  
  cat(sprintf("SUMMARY: UP [sig:%d suggest:%d border:%d] | DOWN [sig:%d suggest:%d border:%d]\n\n",
              n_up_sig, n_up_suggest, n_up_borderline,
              n_down_sig, n_down_suggest, n_down_borderline))
}

#############################
## VAL BATCH: CXCR7_vs_0s (TOTAL EFFECT)
#############################

cat("\n\n")
cat(strrep("▓", 50), "\n")
cat("VAL BATCH: CXCR7_vs_0s (TOTAL EFFECT - includes vehicle)\n")
cat("For comparison with INIT batch\n")
cat(strrep("▓", 50), "\n\n")

val_cxcr7_0s_combos <- up_anno %>%
  filter(batch == "VAL", comparison == "CXCR7_vs_0s") %>%
  distinct(timepoint_sec) %>%
  arrange(timepoint_sec) %>%
  pull(timepoint_sec)

for (tp in val_cxcr7_0s_combos) {
  cat("\n")
  cat(strrep("═", 150), "\n")
  cat(sprintf("TIMEPOINT: %d seconds\n", tp))
  cat(strrep("═", 150), "\n\n")
  
  up_data_all <- up_anno %>%
    filter(batch == "VAL", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.20) %>%
    distinct(pathway, .keep_all = TRUE) %>%
    arrange(pvalue) %>%
    select(pathway, pvalue, ratio)
  
  print_extended_pathway_table(up_data_all, "UP-REGULATED")
  
  down_data_all <- down_anno %>%
    filter(batch == "VAL", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.20) %>%
    distinct(pathway, .keep_all = TRUE) %>%
    arrange(pvalue) %>%
    select(pathway, pvalue, ratio)
  
  print_extended_pathway_table(down_data_all, "DOWN-REGULATED")
  
  n_up_sig <- nrow(up_anno %>% filter(batch == "VAL", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.05))
  n_up_suggest <- nrow(up_anno %>% filter(batch == "VAL", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.10, pvalue >= 0.05))
  n_up_borderline <- nrow(up_anno %>% filter(batch == "VAL", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.20, pvalue >= 0.10))
  
  n_down_sig <- nrow(down_anno %>% filter(batch == "VAL", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.05))
  n_down_suggest <- nrow(down_anno %>% filter(batch == "VAL", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.10, pvalue >= 0.05))
  n_down_borderline <- nrow(down_anno %>% filter(batch == "VAL", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.20, pvalue >= 0.10))
  
  cat(sprintf("SUMMARY: UP [sig:%d suggest:%d border:%d] | DOWN [sig:%d suggest:%d border:%d]\n\n",
              n_up_sig, n_up_suggest, n_up_borderline,
              n_down_sig, n_down_suggest, n_down_borderline))
}

#############################
## INIT BATCH: CXCR7_vs_0s (ONLY COMPARISON AVAILABLE)
#############################

cat("\n\n")
cat(strrep("▓", 50), "\n")
cat("INIT BATCH: CXCR7_vs_0s (ONLY AVAILABLE COMPARISON)\n")
cat("⚠️  Missing CXCR7_vs_DMSO comparison\n")
cat(strrep("▓", 50), "\n\n")

init_cxcr7_0s_combos <- up_anno %>%
  filter(batch == "INIT", comparison == "CXCR7_vs_0s") %>%
  distinct(timepoint_sec) %>%
  arrange(timepoint_sec) %>%
  pull(timepoint_sec)

for (tp in init_cxcr7_0s_combos) {
  cat("\n")
  cat(strrep("═", 150), "\n")
  cat(sprintf("TIMEPOINT: %d seconds\n", tp))
  cat(strrep("═", 150), "\n\n")
  
  up_data_all <- up_anno %>%
    filter(batch == "INIT", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.20) %>%
    distinct(pathway, .keep_all = TRUE) %>%
    arrange(pvalue) %>%
    select(pathway, pvalue, ratio)
  
  print_extended_pathway_table(up_data_all, "UP-REGULATED")
  
  down_data_all <- down_anno %>%
    filter(batch == "INIT", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.20) %>%
    distinct(pathway, .keep_all = TRUE) %>%
    arrange(pvalue) %>%
    select(pathway, pvalue, ratio)
  
  print_extended_pathway_table(down_data_all, "DOWN-REGULATED")
  
  n_up_sig <- nrow(up_anno %>% filter(batch == "INIT", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.05))
  n_up_suggest <- nrow(up_anno %>% filter(batch == "INIT", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.10, pvalue >= 0.05))
  n_up_borderline <- nrow(up_anno %>% filter(batch == "INIT", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.20, pvalue >= 0.10))
  
  n_down_sig <- nrow(down_anno %>% filter(batch == "INIT", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.05))
  n_down_suggest <- nrow(down_anno %>% filter(batch == "INIT", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.10, pvalue >= 0.05))
  n_down_borderline <- nrow(down_anno %>% filter(batch == "INIT", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.20, pvalue >= 0.10))
  
  cat(sprintf("SUMMARY: UP [sig:%d suggest:%d border:%d] | DOWN [sig:%d suggest:%d border:%d]\n\n",
              n_up_sig, n_up_suggest, n_up_borderline,
              n_down_sig, n_down_suggest, n_down_borderline))
}

#############################
## COMPARISON: DMSO vs 0s (Vehicle effects)
#############################

cat("\n\n")
cat(strrep("▓", 50), "\n")
cat("VAL BATCH: DMSO_vs_0s (VEHICLE CONTROL EFFECTS)\n")
cat("Understanding baseline platelet effects from DMSO\n")
cat(strrep("▓", 50), "\n\n")

dmso_combos <- up_anno %>%
  filter(batch == "VAL", comparison == "DMSO_vs_0s") %>%
  distinct(timepoint_sec) %>%
  arrange(timepoint_sec) %>%
  pull(timepoint_sec)

for (tp in dmso_combos) {
  cat("\n")
  cat(strrep("═", 150), "\n")
  cat(sprintf("TIMEPOINT: %d seconds\n", tp))
  cat(strrep("═", 150), "\n\n")
  
  up_data_all <- up_anno %>%
    filter(batch == "VAL", comparison == "DMSO_vs_0s", timepoint_sec == tp, pvalue < 0.20) %>%
    distinct(pathway, .keep_all = TRUE) %>%
    arrange(pvalue) %>%
    select(pathway, pvalue, ratio)
  
  print_extended_pathway_table(up_data_all, "UP-REGULATED")
  
  down_data_all <- down_anno %>%
    filter(batch == "VAL", comparison == "DMSO_vs_0s", timepoint_sec == tp, pvalue < 0.20) %>%
    distinct(pathway, .keep_all = TRUE) %>%
    arrange(pvalue) %>%
    select(pathway, pvalue, ratio)
  
  print_extended_pathway_table(down_data_all, "DOWN-REGULATED")
  
  n_up_sig <- nrow(up_anno %>% filter(batch == "VAL", comparison == "DMSO_vs_0s", timepoint_sec == tp, pvalue < 0.05))
  n_up_suggest <- nrow(up_anno %>% filter(batch == "VAL", comparison == "DMSO_vs_0s", timepoint_sec == tp, pvalue < 0.10, pvalue >= 0.05))
  n_up_borderline <- nrow(up_anno %>% filter(batch == "VAL", comparison == "DMSO_vs_0s", timepoint_sec == tp, pvalue < 0.20, pvalue >= 0.10))
  
  n_down_sig <- nrow(down_anno %>% filter(batch == "VAL", comparison == "DMSO_vs_0s", timepoint_sec == tp, pvalue < 0.05))
  n_down_suggest <- nrow(down_anno %>% filter(batch == "VAL", comparison == "DMSO_vs_0s", timepoint_sec == tp, pvalue < 0.10, pvalue >= 0.05))
  n_down_borderline <- nrow(down_anno %>% filter(batch == "VAL", comparison == "DMSO_vs_0s", timepoint_sec == tp, pvalue < 0.20, pvalue >= 0.10))
  
  cat(sprintf("SUMMARY: UP [sig:%d suggest:%d border:%d] | DOWN [sig:%d suggest:%d border:%d]\n\n",
              n_up_sig, n_up_suggest, n_up_borderline,
              n_down_sig, n_down_suggest, n_down_borderline))
}

#############################
## SUMMARY STATISTICS
#############################

cat("\n\n", strrep("█", 150), "\n")
cat("SUMMARY STATISTICS: Pathway counts at different p-value thresholds\n")
cat(strrep("█", 150), "\n\n")

summary_extended <- data.frame(
  Comparison = character(),
  Batch = character(),
  Timepoint = character(),
  Direction = character(),
  p05_count = integer(),
  p10_count = integer(),
  p20_count = integer()
)

# VAL CXCR7_vs_DMSO
for (tp in sort(unique(up_anno$timepoint_sec[up_anno$batch == "VAL" & up_anno$comparison == "CXCR7_vs_DMSO"]))) {
  n_05 <- nrow(up_anno %>% filter(batch == "VAL", comparison == "CXCR7_vs_DMSO", timepoint_sec == tp, pvalue < 0.05))
  n_10 <- nrow(up_anno %>% filter(batch == "VAL", comparison == "CXCR7_vs_DMSO", timepoint_sec == tp, pvalue < 0.10))
  n_20 <- nrow(up_anno %>% filter(batch == "VAL", comparison == "CXCR7_vs_DMSO", timepoint_sec == tp, pvalue < 0.20))
  summary_extended <- rbind(summary_extended, data.frame(
    Comparison = "CXCR7_vs_DMSO", Batch = "VAL", Timepoint = sprintf("%ds", tp), Direction = "UP",
    p05_count = n_05, p10_count = n_10, p20_count = n_20
  ))
  
  n_05 <- nrow(down_anno %>% filter(batch == "VAL", comparison == "CXCR7_vs_DMSO", timepoint_sec == tp, pvalue < 0.05))
  n_10 <- nrow(down_anno %>% filter(batch == "VAL", comparison == "CXCR7_vs_DMSO", timepoint_sec == tp, pvalue < 0.10))
  n_20 <- nrow(down_anno %>% filter(batch == "VAL", comparison == "CXCR7_vs_DMSO", timepoint_sec == tp, pvalue < 0.20))
  summary_extended <- rbind(summary_extended, data.frame(
    Comparison = "CXCR7_vs_DMSO", Batch = "VAL", Timepoint = sprintf("%ds", tp), Direction = "DOWN",
    p05_count = n_05, p10_count = n_10, p20_count = n_20
  ))
}

# VAL CXCR7_vs_0s
for (tp in sort(unique(up_anno$timepoint_sec[up_anno$batch == "VAL" & up_anno$comparison == "CXCR7_vs_0s"]))) {
  n_05 <- nrow(up_anno %>% filter(batch == "VAL", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.05))
  n_10 <- nrow(up_anno %>% filter(batch == "VAL", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.10))
  n_20 <- nrow(up_anno %>% filter(batch == "VAL", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.20))
  summary_extended <- rbind(summary_extended, data.frame(
    Comparison = "CXCR7_vs_0s", Batch = "VAL", Timepoint = sprintf("%ds", tp), Direction = "UP",
    p05_count = n_05, p10_count = n_10, p20_count = n_20
  ))
  
  n_05 <- nrow(down_anno %>% filter(batch == "VAL", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.05))
  n_10 <- nrow(down_anno %>% filter(batch == "VAL", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.10))
  n_20 <- nrow(down_anno %>% filter(batch == "VAL", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.20))
  summary_extended <- rbind(summary_extended, data.frame(
    Comparison = "CXCR7_vs_0s", Batch = "VAL", Timepoint = sprintf("%ds", tp), Direction = "DOWN",
    p05_count = n_05, p10_count = n_10, p20_count = n_20
  ))
}

# INIT CXCR7_vs_0s
for (tp in sort(unique(up_anno$timepoint_sec[up_anno$batch == "INIT" & up_anno$comparison == "CXCR7_vs_0s"]))) {
  n_05 <- nrow(up_anno %>% filter(batch == "INIT", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.05))
  n_10 <- nrow(up_anno %>% filter(batch == "INIT", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.10))
  n_20 <- nrow(up_anno %>% filter(batch == "INIT", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.20))
  summary_extended <- rbind(summary_extended, data.frame(
    Comparison = "CXCR7_vs_0s", Batch = "INIT", Timepoint = sprintf("%ds", tp), Direction = "UP",
    p05_count = n_05, p10_count = n_10, p20_count = n_20
  ))
  
  n_05 <- nrow(down_anno %>% filter(batch == "INIT", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.05))
  n_10 <- nrow(down_anno %>% filter(batch == "INIT", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.10))
  n_20 <- nrow(down_anno %>% filter(batch == "INIT", comparison == "CXCR7_vs_0s", timepoint_sec == tp, pvalue < 0.20))
  summary_extended <- rbind(summary_extended, data.frame(
    Comparison = "CXCR7_vs_0s", Batch = "INIT", Timepoint = sprintf("%ds", tp), Direction = "DOWN",
    p05_count = n_05, p10_count = n_10, p20_count = n_20
  ))
}

cat("Pathway counts at different p-value thresholds:\n\n")
cat("┌──────────────────┬──────────┬──────────┬───────────┬────────┬────────┬────────┐\n")
cat("│ Comparison       │ Batch    │ Time     │ Direction │ p<0.05 │ p<0.10 │ p<0.20 │\n")
cat("├──────────────────┼──────────┼──────────┼───────────┼────────┼────────┼────────┤\n")

for (i in seq_len(nrow(summary_extended))) {
  cat(sprintf("│ %-16s │ %-8s │ %-8s │ %-9s │ %6d │ %6d │ %6d │\n",
              summary_extended$Comparison[i],
              summary_extended$Batch[i],
              summary_extended$Timepoint[i],
              summary_extended$Direction[i],
              summary_extended$p05_count[i],
              summary_extended$p10_count[i],
              summary_extended$p20_count[i]))
}

cat("└──────────────────┴──────────┴──────────┴───────────┴────────┴────────┴────────┘\n\n")

cat("Legend:\n")
cat("  p<0.05 = Statistically significant (gold standard)\n")
cat("  p<0.10 = Suggestive of significance (may be real)\n")
cat("  p<0.20 = Borderline/weak signal (worth noting)\n\n")

cat(strrep("█", 150), "\n")
cat("✓ EXTENDED PATHWAY OVERVIEW COMPLETE\n")
cat(strrep("█", 150), "\n\n")









