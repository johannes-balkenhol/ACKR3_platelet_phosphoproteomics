###############################################################
## COMBINED PIPELINE FOR NEW + OLD PHOSPHOPROTEOMICS DATASET
## Reactome Enrichment • Up/Down Regulation • Combined Plots
###############################################################

#############################
## 0) Refined Phosphoproteomics Enrichment and Visualization Script
#############################
suppressPackageStartupMessages({
  library(remotes)
  library(rlist)
  library(sjmisc)
  library(stringr)
  library(plyr)
  library(dplyr)
  library(tidyr)
  
  library(annotate)
  library(reactome.db)
  library(ReactomePA)
  library(org.Hs.eg.db) 
  
  library(RColorBrewer)
  library(basicPlotteR)
  library(pheatmap)
  library(ggplot2)
  library(calibrate)
  library(enrichplot)
  library(cowplot)
  library(limma)
  
  library(PhosR)
  library(clusterProfiler)
  library(directPA)
})




###############################################################
# 1) HARMONIZATION FUNCTIONS
###############################################################

## ------------------------------------------------------------
## A) Harmonize NEW validation datasets
## ------------------------------------------------------------
harmonize_new <- function(df) {
  df %>%
    dplyr::rename(
      uniprot_id = uniprot,
      name       = symbol,
      Average    = AveExpr,
      PValue     = adj.P.Val,
      PSite      = psite
    ) %>%
    dplyr::select(uniprot_id, name, Average, logFC, PValue, PSite) %>%
    as.data.frame()
}


## ------------------------------------------------------------
## B) Harmonize OLD initial datasets  
##    structure from df$id: "UNIPROT;GENE;PSITE;PEPTIDE;INDEX"
## ------------------------------------------------------------
harmonize_old <- function(df) {
  
  parts <- strsplit(df$id, ";")
  
  df$uniprot_id <- sapply(parts, `[`, 1)
  df$name       <- sapply(parts, `[`, 2)
  df$PSite      <- sapply(parts, `[`, 3)
  df$Peptide    <- sapply(parts, `[`, 4)
  
  df %>%
    dplyr::transmute(
      uniprot_id,
      name,
      Average = AveExpr,
      logFC,
      PValue = P.Value,
      PSite
    ) %>%
    as.data.frame()
}


## ------------------------------------------------------------
## C) Collapse AFTER intersection (not used yet here)
##    - group by uniprot_id
##    - select phosphosite with highest |logFC|
## ------------------------------------------------------------
collapse_uniprot <- function(df) {
  df %>%
    mutate(abs_logFC = abs(logFC)) %>%
    group_by(uniprot_id) %>%
    slice_max(abs_logFC, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(uniprot_id, name, PSite, Average, logFC, PValue)
}



###############################################################
# 2) LOAD & HARMONIZE ALL DATA
###############################################################

## ----------------------------  
## NEW validation datasets
## ----------------------------

dfs_new <- list(
  top.10, top.600, top.1800,
  top.10.dmso.vs.0s, top.600.dmso.vs.0s, top.1800.dmso.vs.0s,
  top.10.cxcr7.vs.0s, top.600.cxcr7.vs.0s, top.1800.cxcr7.vs.0s
)

names(dfs_new) <- c(
  "val_10.dmso.vs.cxcr7", "val_600.dmso.vs.cxcr7", "val_1800.dmso.vs.cxcr7",  #CXCR7 vs DMSO (at each timepoint)
  "val_10.dmso.vs.0s",    "val_600.dmso.vs.0s",    "val_1800.dmso.vs.0s",     #DMSO vs 0s
  "val_10.cxcr7.vs.0s",   "val_600.cxcr7.vs.0s",   "val_1800.cxcr7.vs.0s"     #CXCR7 vs 0s
)

dfs_new_raw <- lapply(dfs_new, harmonize_new)


## ----------------------------  
## OLD initial datasets
## ----------------------------

old_path  <- "SubProjects/CXCR7_initial/phosphoproteomics/data/processed_data"
old_files <- list.files(old_path, pattern = "cxcr7", full.names = TRUE)

old_tables_raw        <- lapply(old_files, read.delim)
old_tables_harmonized <- lapply(old_tables_raw, harmonize_old)

names(old_tables_harmonized) <- c(
  "init_10.cxcr7.vs.0s","init_30.cxcr7.vs.0s","init_60.cxcr7.vs.0s",  #CXCR7 vs 0s
  "init_300.cxcr7.vs.0s","init_600.cxcr7.vs.0s","init_900.cxcr7.vs.0s",  #CXCR7 vs 0s
  "init_1800.cxcr7.vs.0s"  #CXCR7 vs 0s
)



###############################################################
# 3) INTERSECT PHOSPHOSITES ACROSS ALL DATASETS
###############################################################

## ----------------------------  
## A) Add true phosphosite key (uniprot_id + PSite)  
## ----------------------------
add_key <- function(df) {
  df$phospho_id <- paste(df$uniprot_id, df$PSite, sep = "@")
  df
}

new_keyed  <- lapply(dfs_new_raw, add_key)
init_keyed <- lapply(old_tables_harmonized, add_key)


## ----------------------------  
## B) Compute intersection across ALL datasets  
## ----------------------------
common_phosphosites <- Reduce(
  intersect,
  lapply(c(new_keyed, init_keyed), function(df) df$phospho_id)
)

cat("Number of common phosphosites:", length(common_phosphosites), "\n")


## ----------------------------  
## C) Filter each dataset to intersected sites  
## ----------------------------
filter_common <- function(df) {
  df[df$phospho_id %in% common_phosphosites, ]
}

dfs_new_intersect  <- lapply(new_keyed,  filter_common)
init_intersect     <- lapply(init_keyed, filter_common)



###############################################################
# 4) BUILD FINAL all_inputs (INTERSECTED, NOT COLLAPSED)
###############################################################

## Combine validation + initial intersected tables
all_inputs <- c(dfs_new_intersect, init_intersect)

## Sort each dataset by Uniprot + PSite for consistency
all_inputs <- lapply(all_inputs, function(df) {
  df[order(df$uniprot_id, df$PSite), ]
})

## Report
cat("\nFinal all_inputs created:\n")
cat("Validation datasets:", length(dfs_new_intersect), "\n")
cat("Initial datasets:   ", length(init_intersect), "\n")
cat("Total tables:       ", length(all_inputs), "\n\n")

## Optional: show dimensions of each table
for (nm in names(all_inputs)) {
  cat(sprintf("%s → %d rows, %d columns\n",
              nm, nrow(all_inputs[[nm]]), ncol(all_inputs[[nm]])))
}

###############################################################
# 5) COLLAPSE EACH DATASET INDIVIDUALLY BY UNIPROT
#    choose PSite with max |logFC| per protein
###############################################################

collapse_uniprot <- function(df) {
  df %>%
    mutate(abs_logFC = abs(logFC)) %>%
    group_by(uniprot_id) %>%
    slice_max(abs_logFC, with_ties = FALSE) %>%
    ungroup() %>%
    select(uniprot_id, name, PSite, Average, logFC, PValue)
}

# Apply to ALL datasets (validation + initial)
all_inputs_collapsed <- lapply(all_inputs, collapse_uniprot)

# Sort rows
all_inputs_collapsed <- lapply(all_inputs_collapsed, function(df) {
  df[order(df$uniprot_id), ]
})

###############################################################
# Build consistent UniProt set (intersection across all collapsed)
###############################################################

common_uniprot <- Reduce(
  intersect,
  lapply(all_inputs_collapsed, function(df) df$uniprot_id)
)

# Filter each dataset to shared proteins
all_inputs_collapsed <- lapply(all_inputs_collapsed, function(df) {
  df[df$uniprot_id %in% common_uniprot, ]
})

# Final sorting again
all_inputs_collapsed <- lapply(all_inputs_collapsed, function(df) {
  df[order(df$uniprot_id), ]
})

###############################################################
# BUILD LOGFC MATRIX
###############################################################

logfc_matrix <- do.call(cbind, lapply(all_inputs_collapsed, `[[`, "logFC"))
colnames(logfc_matrix) <- names(all_inputs_collapsed)
rownames(logfc_matrix) <- all_inputs_collapsed[[1]]$uniprot_id

dim(logfc_matrix)
head(logfc_matrix)


gene_symbols <- all_inputs_collapsed[[1]]$name


###############################################################
### ----------------------------------------------------------
### Pathway enrichemnt
### ----------------------------------------------------------
###############################################################


## Build a logFC matrix from the collapsed tables.
#Tc.gene <- as.matrix(cbind(input[[1]]$logFC, input[[2]]$logFC, input[[3]]$logFC,
#                           input[[4]]$logFC, input[[5]]$logFC, input[[6]]$logFC,
#                           input[[7]]$logFC, input[[8]]$logFC, input[[9]]$logFC))
#rownames(Tc.gene) <- paste(sapply(strsplit(rownames(input[[1]]), ";"), "[[", 2))
#colnames(Tc.gene) <- names_input



Tc.gene <- logfc_matrix
rownames(Tc.gene) <- gene_symbols
names_input <- colnames(Tc.gene)



#############################
## 1.2 Prepare Reactome pathways (pathway loading and matching to platelet)
#############################

# Step 1: Load Reactome pathways as before
pathways <- as.list(reactomePATHID2EXTID)
path_names <- as.list(reactomePATHID2NAME)
name_id <- match(names(pathways), names(path_names))
names(pathways) <- unlist(path_names)[name_id]

# ✓ KEEP ONLY HOMO SAPIENS
pathways <- pathways[grepl("Homo sapiens", names(pathways), ignore.case = TRUE)]

# Convert to gene symbols
pathways <- lapply(pathways, function(path) {
  gene_name <- unname(getSYMBOL(path, data = "org.Hs.eg"))
  toupper(unique(gene_name))
})

cat("✓ Total Reactome pathways (Homo sapiens):", length(pathways), "\n\n")

# Step 2: Clean pathway names (remove prefix)
pathway_names_clean <- names(pathways)
pathway_names_clean <- gsub("Homo sapiens: ", "", pathway_names_clean, ignore.case = TRUE)
pathway_names_clean <- gsub("Homo sapiens >> ", "", pathway_names_clean, ignore.case = TRUE)
pathway_names_clean <- trimws(pathway_names_clean)
names(pathways) <- pathway_names_clean

cat("Sample cleaned names:\n")
print(head(pathway_names_clean, 10))
cat("\n")

# Step 3: Your curated list
manual_path_refined <- c(
  "Platelet activation",
  "Platelet degranulation",
  "Response to elevated platelet cytosolic Ca2+",
  "Platelet Adhesion to exposed collagen",
  "Hemostasis",
  "Platelet calcium homeostasis",
  "Platelet Aggregation (Plug Formation)",
  "Clathrin-mediated endocytosis",
  "Cargo recognition for clathrin-mediated endocytosis",
  "Integrin cell surface interactions",
  "Integrin signaling",
  "Signaling by Receptor Tyrosine Kinases",
  "Signaling by GPCR",
  "GPCR downstream signalling",
  "G beta:gamma signalling through PI3Kgamma",
  "Activation of G protein gated Potassium channels",
  "G beta:gamma signalling through PLC beta",
  "G beta:gamma signalling through BTK",
  "G beta:gamma signalling through CDC42",
  "G protein gated Potassium channels",
  "GPVI-mediated activation cascade",
  "FCGR activation",
  "Fcgamma receptor (FCGR) dependent phagocytosis",
  "PKA activation",
  "PKA activation in glucagon signalling",
  "Adenylate cyclase activating pathway",
  "Adenylate cyclase inhibitory pathway",
  "CREB1 phosphorylation through the activation of Adenylate Cyclase",
  "PKA-mediated phosphorylation of CREB",
  "PKA-mediated phosphorylation of key metabolic factors",
  "Gastrin-CREB signalling pathway via PKC and MAPK",
  "cGMP effects",
  "Nitric oxide stimulates guanylate cyclase",
  "Energy dependent regulation of mTOR by LKB1-AMPK",
  "Activation of AMPK downstream of NMDARs",
  "AMPK inhibits chREBP transcriptional activation activity",
  "RAF activation",
  "Signalling to ERKs",
  "ERK/MAPK targets",
  "ERKs are inactivated",
  "Signaling by BRAF and RAF1 fusions",
  "Negative regulation of MAPK pathway",
  "GRB2:SOS provides linkage to MAPK signaling for Integrins",
  "Regulation of RAS by GAPs",
  "Signaling by Rho GTPases",
  "RHOA GTPase cycle",
  "RHOB GTPase cycle",
  "RHOC GTPase cycle",
  "RHOV GTPase cycle",
  "RHO GTPase cycle",
  "RHO GTPase Effectors",
  "RHO GTPases Activate WASPs and WAVEs",
  "RHO GTPases Activate NADPH Oxidases",
  "PIP3 activates AKT signaling",
  "AKT phosphorylates targets in the cytosol",
  "Negative regulation of the PI3K/AKT network",
  "PI3K Cascade",
  "Effects of PIP2 hydrolysis",
  "Ca-dependent events",
  "MTOR signalling",
  "mTORC1-mediated signalling",
  "Amino acids regulate mTORC1",
  "Autophagy",
  "Signaling by Non-Receptor Tyrosine Kinases",
  "DAP12 signaling",
  "Extracellular matrix organization",
  "Degradation of the extracellular matrix",
  "Membrane Trafficking",
  "Arachidonic acid metabolism",
  "Synthesis of Prostaglandins (PG) and Thromboxanes (TX)",
  "Thromboxane signalling through TP receptor",
  "Eicosanoids",
  "DAG and IP3 signaling",
  "Inositol phosphate metabolism",
  "Synthesis of IP3 and IP4 in the cytosol",
  "Arachidonate production from DAG",
  "Eicosanoid ligand-binding receptors",
  "Acyl chain remodeling of DAG and TAG",
  "PLC beta mediated events",
  "Role of phospholipids in phagocytosis"
)

# Step 4: BETTER MATCHING FUNCTION
match_pathway <- function(query, all_names) {
  # Try exact match first (case-insensitive)
  idx <- which(tolower(all_names) == tolower(query))
  if (length(idx) > 0) return(all_names[idx[1]])
  
  # If no exact match, try greedy partial match
  # This handles cases like "Platelet activation" → "Platelet activation, signaling and aggregation"
  idx <- grep(tolower(gsub(" ", ".*", query)), tolower(all_names))
  if (length(idx) > 0) return(all_names[idx[1]])
  
  return(NA)
}

# Step 5: Match all 80 pathways
selected_pathways_filtered <- list()
matched_count <- 0
not_found <- character()

cat("Matching your 80 curated pathways...\n\n")

for (path in manual_path_refined) {
  found_name <- match_pathway(path, pathway_names_clean)
  
  if (!is.na(found_name)) {
    selected_pathways_filtered[[path]] <- pathways[[found_name]]
    cat("✓", path, "\n")
    matched_count <- matched_count + 1
  } else {
    cat("✗", path, "\n")
    not_found <- c(not_found, path)
  }
}

# Step 6: SUMMARY
cat("\n", strrep("=", 70), "\n")
cat("FINAL MATCHING SUMMARY\n")
cat(strrep("=", 70), "\n\n")

cat("Total curated pathways:", length(manual_path_refined), "\n")
cat("Successfully matched:", matched_count, "\n")
cat("NOT found:", length(not_found), "\n")
cat("Match rate:", round(100 * matched_count / length(manual_path_refined), 1), "%\n\n")

if (length(not_found) > 0) {
  cat("NOT FOUND:\n")
  for (p in not_found) cat(" -", p, "\n")
}

cat("\n✓ Ready for enrichment with", length(selected_pathways_filtered), "pathways!\n\n")

# Step 7: Save for reference
selected_pathways_df <- data.frame(
  pathway_number = seq_along(selected_pathways_filtered),
  pathway_name = names(selected_pathways_filtered),
  n_genes = sapply(selected_pathways_filtered, length)
)

write.csv(selected_pathways_df, 
          "selected_pathways_matched.csv", 
          row.names = FALSE)

cat("✓ Saved to: selected_pathways_matched.csv\n\n")



#############################
## ENRICHMENT WITH 80 PATHWAYS - FINAL FIX
#############################

cat("\n", strrep("=", 70), "\n")
cat("RUNNING ENRICHMENT WITH 80 PATHWAYS (FIXED!)\n")
cat(strrep("=", 70), "\n\n")

up_results <- list()
down_results <- list()

pathways <- selected_pathways_filtered

# Function to process enrichment results
process_enrichment <- function(res, selected_pathways_filtered) {
  # Convert matrix to data frame
  path3 <- as.data.frame(res)
  path3$pathway <- rownames(path3)
  
  # Create pathway sizes in matching format
  pathways_df <- data.frame(
    pathway = names(selected_pathways_filtered),
    pw.size = sapply(selected_pathways_filtered, length),
    stringsAsFactors = FALSE
  )
  
  # Merge by pathway name
  path4 <- merge(path3, pathways_df, by = "pathway", all = FALSE)
  
  # Clean up columns
  colnames(path4)[2] <- "pvalue"
  colnames(path4)[3] <- "number.substrates"
  colnames(path4)[4] <- "substrates"
  
  # Convert to numeric
  path4$pvalue <- as.numeric(path4$pvalue)
  path4$number.substrates <- as.numeric(path4$number.substrates)
  path4$pw.size <- as.numeric(path4$pw.size)
  path4$ratio <- path4$number.substrates / path4$pw.size
  path4$neg_log10_p <- -log10(path4$pvalue)
  
  # Sort by ratio (descending)
  path4 <- path4[order(path4$ratio, decreasing = TRUE), ]
  
  return(path4)
}

# UP-regulated
for (i in seq_len(ncol(Tc.gene))) {
  res <- pathwayRankBasedEnrichment(Tc.gene[, i], annotation = pathways, alter = "greater")
  path4 <- process_enrichment(res, selected_pathways_filtered)
  
  up_results[[names_input[i]]] <- path4
  
  sig <- sum(path4$pvalue < 0.05)
  cat(sprintf("✓ UP   - %s: %d pathways, %d significant\n", names_input[i], nrow(path4), sig))
}

cat("\n")

# DOWN-regulated
for (i in seq_len(ncol(Tc.gene))) {
  res <- pathwayRankBasedEnrichment(Tc.gene[, i], annotation = pathways, alter = "less")
  path4 <- process_enrichment(res, selected_pathways_filtered)
  
  down_results[[names_input[i]]] <- path4
  
  sig <- sum(path4$pvalue < 0.05)
  cat(sprintf("✓ DOWN - %s: %d pathways, %d significant\n", names_input[i], nrow(path4), sig))
}

#############################
## BUILD SUMMARIES
#############################

cat("\n", strrep("=", 70), "\n")
cat("BUILDING SUMMARIES\n")
cat(strrep("=", 70), "\n\n")

up_summary_long_all <- lapply(names(up_results), function(t) {
  df <- up_results[[t]]
  if (!is.null(df) && nrow(df) > 0) {
    df %>% dplyr::select(pathway, pvalue, ratio, substrates, neg_log10_p) %>%
      mutate(time = t)
  }
}) %>% bind_rows()

down_summary_long_all <- lapply(names(down_results), function(t) {
  df <- down_results[[t]]
  if (!is.null(df) && nrow(df) > 0) {
    df %>% dplyr::select(pathway, pvalue, ratio, substrates, neg_log10_p) %>%
      mutate(time = t)
  }
}) %>% bind_rows()

cat("UP-regulated:  ", nrow(up_summary_long_all), "results |", 
    n_distinct(up_summary_long_all$pathway), "unique pathways |",
    sum(up_summary_long_all$pvalue < 0.05), "significant (p<0.05)\n\n")

cat("DOWN-regulated:", nrow(down_summary_long_all), "results |", 
    n_distinct(down_summary_long_all$pathway), "unique pathways |",
    sum(down_summary_long_all$pvalue < 0.05), "significant (p<0.05)\n\n")


#############################
## DISPLAY RESULTS
#############################

cat(strrep("=", 70), "\n")
cat("TOP 20 UP-REGULATED PATHWAYS (by p-value)\n")
cat(strrep("=", 70), "\n\n")

top_up <- up_summary_long_all %>%
  arrange(pvalue) %>%
  head(20) %>%
  select(pathway, pvalue, ratio, time)

print(top_up)

cat("\n\n", strrep("=", 70), "\n")
cat("TOP 20 DOWN-REGULATED PATHWAYS (by p-value)\n")
cat(strrep("=", 70), "\n\n")

top_down <- down_summary_long_all %>%
  arrange(pvalue) %>%
  head(20) %>%
  select(pathway, pvalue, ratio, time)

print(top_down)

#############################
## SHOW PATHWAYS BY SIGNIFICANCE
#############################

cat("\n\n", strrep("=", 70), "\n")
cat("SIGNIFICANTLY ENRICHED PATHWAYS (p < 0.05)\n")
cat(strrep("=", 70), "\n\n")

sig_up <- up_summary_long_all %>% 
  filter(pvalue < 0.05) %>%
  arrange(pvalue)

sig_down <- down_summary_long_all %>%
  filter(pvalue < 0.05) %>%
  arrange(pvalue)

cat("UP-regulated (", nrow(sig_up), " total):\n\n")
if (nrow(sig_up) > 0) {
  print(sig_up %>% select(pathway, pvalue, ratio, time))
} else {
  cat("  (None)\n")
}

cat("\n\nDOWN-regulated (", nrow(sig_down), " total):\n\n")
if (nrow(sig_down) > 0) {
  print(sig_down %>% select(pathway, pvalue, ratio, time))
} else {
  cat("  (None)\n")
}

cat("\n✓ ENRICHMENT ANALYSIS COMPLETE!\n\n")
cat("Objects ready:\n")
cat("  - up_summary_long_all\n")
cat("  - down_summary_long_all\n")
cat("  - up_results (by timepoint)\n")
cat("  - down_results (by timepoint)\n\n")





#############################
## 3. Save Enrichment Results
#############################

# Create output directory
output_dir <- "/mnt/user-data/outputs/enrichment_results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Save UP-regulated enrichment tables
for (time in names(up_results)) {
  write.table(up_results[[time]],
              file = file.path(output_dir, paste0("greater_", time, "_pathways.txt")),
              sep = "\t", row.names = FALSE)
}

# Save DOWN-regulated enrichment tables
for (time in names(down_results)) {
  write.table(down_results[[time]],
              file = file.path(output_dir, paste0("less_", time, "_pathways.txt")),
              sep = "\t", row.names = FALSE)
}




#############################
## 4. (optional) or FILTER TO SELECTED PATHWAYS
#############################

manual_filter <- c(
  # === CORE PLATELET ACTIVATION & ADHESION (7) ===
  "Platelet activation",
  "Platelet degranulation",
  "Response to elevated platelet cytosolic Ca2+",
  "Platelet Adhesion to exposed collagen",
  "Hemostasis",
  "Platelet calcium homeostasis",
  "Platelet Aggregation (Plug Formation)",
  
  # === ENDOCYTOSIS - KEY FOR ACKR3! (2) ===
  "Clathrin-mediated endocytosis",
  "Cargo recognition for clathrin-mediated endocytosis",
  
  # === INTEGRIN & CELL SIGNALING (3) ===
  "Integrin cell surface interactions",
  "Integrin signaling",
  "Signaling by Receptor Tyrosine Kinases",
  
  # === GPCR & G-PROTEIN SIGNALING (8) ⭐ CRITICAL ===
  "Signaling by GPCR",
  "GPCR downstream signalling",
  "G beta:gamma signalling through PI3Kgamma",
  "Activation of G protein gated Potassium channels",
  "G beta:gamma signalling through PLC beta",
  "G beta:gamma signalling through BTK",
  "G beta:gamma signalling through CDC42",
  "G protein gated Potassium channels",
  
  # === COLLAGEN/GPVI - SYK-MEDIATED ACTIVATION (1) ===
  "GPVI-mediated activation cascade",
  
  # === FC GAMMA RECEPTOR - IMMUNE (2) ===
  "FCGR activation",
  "Fcgamma receptor (FCGR) dependent phagocytosis",
  
  # === PKA/cAMP SIGNALING (7) ⭐ NEW ===
  "PKA activation",
  "PKA activation in glucagon signalling",
  "Adenylate cyclase activating pathway",
  "Adenylate cyclase inhibitory pathway",
  "CREB1 phosphorylation through the activation of Adenylate Cyclase",
  "PKA-mediated phosphorylation of CREB",
  "PKA-mediated phosphorylation of key metabolic factors",
  
  # === PKC SIGNALING (1) ⭐ NEW ===
  "Gastrin-CREB signalling pathway via PKC and MAPK",
  
  # === cGMP/PKG SIGNALING (2) ⭐ ENHANCED ===
  "cGMP effects",
  "Nitric oxide stimulates guanylate cyclase",
  
  # === AMPK SIGNALING (3) ⭐ NEW ===
  "Energy dependent regulation of mTOR by LKB1-AMPK",
  "Activation of AMPK downstream of NMDARs",
  "AMPK inhibits chREBP transcriptional activation activity",
  
  # === RAF/MAPK CASCADE (6) ===
  "RAF activation",
  "Signalling to ERKs",
  "ERK/MAPK targets",
  "ERKs are inactivated",
  "Signaling by BRAF and RAF1 fusions",
  "Negative regulation of MAPK pathway",
  
  # === RAS/GRB2 (2) ===
  "GRB2:SOS provides linkage to MAPK signaling for Integrins",
  "Regulation of RAS by GAPs",
  
  # === RHO GTPASE SIGNALING (8) ===
  "Signaling by Rho GTPases",
  "RHOA GTPase cycle",
  "RHOB GTPase cycle",
  "RHOC GTPase cycle",
  "RHOV GTPase cycle",
  "RHO GTPase cycle",
  "RHO GTPase Effectors",
  "RHO GTPases Activate WASPs and WAVEs",
  
  # === ROS/NOX PATHWAY (1) ===
  "RHO GTPases Activate NADPH Oxidases",
  
  # === PI3K/AKT SIGNALING (4) ===
  "PIP3 activates AKT signaling",
  "AKT phosphorylates targets in the cytosol",
  "Negative regulation of the PI3K/AKT network",
  "PI3K Cascade",
  
  # === PHOSPHOLIPID SIGNALING (1) ===
  "Effects of PIP2 hydrolysis",
  
  # === CALCIUM SIGNALING (1) ===
  "Ca-dependent events",
  
  # === METABOLIC SIGNALING (3) ===
  "MTOR signalling",
  "mTORC1-mediated signalling",
  "Amino acids regulate mTORC1",
  
  # === AUTOPHAGY (1) ===
  "Autophagy",
  
  # === NON-RECEPTOR TYROSINE KINASES (1) ===
  "Signaling by Non-Receptor Tyrosine Kinases",
  
  # === PLATELET SPREADING/CYTOSKELETON (1) ===
  "DAP12 signaling",
  
  # === EXTRACELLULAR MATRIX (2) ===
  "Extracellular matrix organization",
  "Degradation of the extracellular matrix",
  
  # === MEMBRANE TRAFFICKING (1) ===
  "Membrane Trafficking",
  
  # === CRITICAL THROMBOXANE PATHWAY (4) ⭐⭐⭐ ===
  "Arachidonic acid metabolism",
  "Synthesis of Prostaglandins (PG) and Thromboxanes (TX)",
  "Thromboxane signalling through TP receptor",
  "Eicosanoids",
  
  # === DAG/IP3 SIGNALING & METABOLISM (4) ⭐⭐⭐ ===
  "DAG and IP3 signaling",
  "Inositol phosphate metabolism",
  "Synthesis of IP3 and IP4 in the cytosol",
  "Arachidonate production from DAG",
  
  # === SPHINGOLIPID/S1P PATHWAY (3) ⭐⭐⭐ CRITICAL! ===
  #"Lysosphingolipid and LPA receptors",
  #"Ceramide signalling",
  #"Sphingolipid metabolism",
  
  # === ADDITIONAL LIPID PATHWAYS (4) ⭐ ===
  "Eicosanoid ligand-binding receptors",
  "Acyl chain remodeling of DAG and TAG",
  "PLC beta mediated events",
  "Role of phospholipids in phagocytosis"
)


manual_filter <- c(
  "Platelet activation",
  "Platelet degranulation",
  "Response to elevated platelet cytosolic Ca2+",
  "Platelet Adhesion to exposed collagen",
  "Hemostasis",
  "Platelet calcium homeostasis",
  "Platelet Aggregation (Plug Formation)",
  "Clathrin-mediated endocytosis",
  "Cargo recognition for clathrin-mediated endocytosis",
  "Integrin cell surface interactions",
  "Integrin signaling",
  "Signaling by Receptor Tyrosine Kinases",
  "Signaling by GPCR",
  "GPCR downstream signalling",
  "G beta:gamma signalling through PI3Kgamma",
  "Activation of G protein gated Potassium channels",
  "G beta:gamma signalling through PLC beta",
  "G beta:gamma signalling through BTK",
  "G beta:gamma signalling through CDC42",
  "G protein gated Potassium channels",
  "GPVI-mediated activation cascade",
  "FCGR activation",
  "Fcgamma receptor (FCGR) dependent phagocytosis",
  "PKA activation",
  "PKA activation in glucagon signalling",
  "Adenylate cyclase activating pathway",
  "Adenylate cyclase inhibitory pathway",
  "CREB1 phosphorylation through the activation of Adenylate Cyclase",
  "PKA-mediated phosphorylation of CREB",
  "PKA-mediated phosphorylation of key metabolic factors",
  "Gastrin-CREB signalling pathway via PKC and MAPK",
  "cGMP effects",
  "Nitric oxide stimulates guanylate cyclase",
  "Energy dependent regulation of mTOR by LKB1-AMPK",
  "Activation of AMPK downstream of NMDARs",
  "AMPK inhibits chREBP transcriptional activation activity",
  "RAF activation",
  "Signalling to ERKs",
  "ERK/MAPK targets",
  "ERKs are inactivated",
  "Signaling by BRAF and RAF1 fusions",
  "Negative regulation of MAPK pathway",
  "GRB2:SOS provides linkage to MAPK signaling for Integrins",
  "Regulation of RAS by GAPs",
  "Signaling by Rho GTPases",
  "RHOA GTPase cycle",
  "RHOB GTPase cycle",
  "RHOC GTPase cycle",
  "RHOV GTPase cycle",
  "RHO GTPase cycle",
  "RHO GTPase Effectors",
  "RHO GTPases Activate WASPs and WAVEs",
  "RHO GTPases Activate NADPH Oxidases",
  "PIP3 activates AKT signaling",
  "AKT phosphorylates targets in the cytosol",
  "Negative regulation of the PI3K/AKT network",
  "PI3K Cascade",
  "Effects of PIP2 hydrolysis",
  "Ca-dependent events",
  "MTOR signalling",
  "mTORC1-mediated signalling",
  "Amino acids regulate mTORC1",
  "Autophagy",
  "Signaling by Non-Receptor Tyrosine Kinases",
  "DAP12 signaling",
  "Extracellular matrix organization",
  "Degradation of the extracellular matrix",
  "Membrane Trafficking",
  "Arachidonic acid metabolism",
  "Synthesis of Prostaglandins (PG) and Thromboxanes (TX)",
  "Thromboxane signalling through TP receptor",
  "Eicosanoids",
  "DAG and IP3 signaling",
  "Inositol phosphate metabolism",
  "Synthesis of IP3 and IP4 in the cytosol",
  "Arachidonate production from DAG",
  "Eicosanoid ligand-binding receptors",
  "Acyl chain remodeling of DAG and TAG",
  "PLC beta mediated events",
  "Role of phospholipids in phagocytosis"
)


#############################
## PARSE COMPARISONS CORRECTLY
#############################

cat("\n", strrep("=", 70), "\n")
cat("PARSING PATHWAY COMPARISONS\n")
cat(strrep("=", 70), "\n\n")

# Use your curated pathway list
selected_pathways <- unique(manual_filter)

# Filter UP-regulated pathways
up_summary_long_filtered <- up_summary_long_all %>% 
  filter(pathway %in% selected_pathways) %>% 
  mutate(pvalue = as.numeric(pvalue),
         ratio = as.numeric(ratio),
         neg_log10_p = -log10(pvalue)) %>%
  # Simple parsing: extract batch (first 3-4 chars), timepoint (first digits after batch), comparison (rest)
  mutate(
    batch = substr(time, 1, regexpr("_", time) - 1),  # before first underscore
    rest = substr(time, regexpr("_", time) + 1, nchar(time)),  # after first underscore
    timepoint_sec = as.numeric(gsub("[^0-9]", "", substr(rest, 1, 10))),  # first digits
    comparison = sub("^[0-9]+\\.", "", rest)  # remove leading digits and dot
  ) %>%
  select(-rest) %>%
  arrange(pvalue)

# Filter DOWN-regulated pathways
down_summary_long_filtered <- down_summary_long_all %>% 
  filter(pathway %in% selected_pathways) %>% 
  mutate(pvalue = as.numeric(pvalue),
         ratio = as.numeric(ratio),
         neg_log10_p = -log10(pvalue)) %>%
  mutate(
    batch = substr(time, 1, regexpr("_", time) - 1),
    rest = substr(time, regexpr("_", time) + 1, nchar(time)),
    timepoint_sec = as.numeric(gsub("[^0-9]", "", substr(rest, 1, 10))),
    comparison = sub("^[0-9]+\\.", "", rest)
  ) %>%
  select(-rest) %>%
  arrange(pvalue)

cat("UP-regulated combinations:", nrow(up_summary_long_filtered), "\n")
cat("DOWN-regulated combinations:", nrow(down_summary_long_filtered), "\n\n")

cat("Unique batches:\n")
print(unique(up_summary_long_filtered$batch))
cat("\nUnique timepoints:\n")
print(unique(up_summary_long_filtered$timepoint_sec))
cat("\nUnique comparisons:\n")
print(unique(up_summary_long_filtered$comparison))

cat("\n", strrep("=", 70), "\n")
cat("Sample UP-regulated:\n")
print(head(up_summary_long_filtered %>% select(pathway, pvalue, batch, timepoint_sec, comparison)))
cat("\nSample DOWN-regulated:\n")
print(head(down_summary_long_filtered %>% select(pathway, pvalue, batch, timepoint_sec, comparison)))



#################################
## 5. ggplot mirrored barplot (SPECIFIC BATCH & COMPARISON)
#################################


# Choose your batch and comparison:
select_batch <- "val"          
select_comparison <- "dmso.vs.cxcr7"

# Get top 20 most significant pathways
top_20_paths <- bind_rows(
  up_summary_long_filtered %>% 
    filter(batch == select_batch, comparison == select_comparison) %>%
    select(pathway, pvalue),
  down_summary_long_filtered %>%
    filter(batch == select_batch, comparison == select_comparison) %>%
    select(pathway, pvalue)
) %>%
  group_by(pathway) %>%
  summarize(mean_pvalue = mean(pvalue), .groups = 'drop') %>%
  arrange(mean_pvalue) %>%
  head(30) %>%
  pull(pathway)

# Filter combined data for top 20 + specific batch/comparison + ALL TIMEPOINTS
combined_top20 <- bind_rows(
  up_summary_long_filtered %>% 
    filter(batch == select_batch, comparison == select_comparison, pathway %in% top_20_paths) %>%
    mutate(value = neg_log10_p, regulation = "up"),
  down_summary_long_filtered %>%
    filter(batch == select_batch, comparison == select_comparison, pathway %in% top_20_paths) %>%
    mutate(value = -neg_log10_p, regulation = "down")
) %>%
  mutate(pathway = factor(pathway, levels = top_20_paths))

# Create timepoint groups (separate bars, not clustered)
combined_top20 <- combined_top20 %>%
  mutate(
    timepoint_group = paste0(timepoint_sec, " sec"),
    timepoint_group = factor(timepoint_group, levels = c("10 sec", "600 sec", "1800 sec"))
  )

# Color palette: light → medium → dark purple (matching image 2)
colors_timepoint <- c(
  "10 sec" = "#E8D5F0",      # Light purple
  "600 sec" = "#B366D9",     # Medium purple
  "1800 sec" = "#4B0082"     # Dark purple
)

# Create plot - bars side by side, NOT dodged
top20_plot <- ggplot(combined_top20, aes(x = pathway, y = value, fill = timepoint_group)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  scale_fill_manual(
    name = "Timepoint",
    values = colors_timepoint,
    guide = guide_legend(reverse = FALSE)
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "-log10(p-value)",
    title = sprintf("Top 20 Pathways: %s | %s", toupper(select_batch), select_comparison)
  ) +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
  theme_cowplot(font_size = 11) +
  theme(
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.background = element_rect(fill = "white", color = NA)
  )

print(top20_plot)

# Save
ggsave("top20_dmso_vs_cxcr7.pdf", plot = top20_plot, width = 10, height = 8, dpi = 300)
ggsave("top20_dmso_vs_cxcr7.png", plot = top20_plot, width = 10, height = 8, dpi = 300)

cat(sprintf("✓ Top 20 pathways: %s | %s\n", select_batch, select_comparison))
cat(sprintf("✓ Timepoints: %s\n", paste(sort(unique(combined_top20$timepoint_sec)), collapse = ", ")))
cat("✓ Saved: top20_dmso_vs_cxcr7.pdf/png\n")






#############################################################
## 6.1 Filter for the 6 desired contrasts (val + init)
#############################################################

target_conditions <- c(
  "val_10.cxcr7.vs.0s",
  "val_600.cxcr7.vs.0s",
  "val_1800.cxcr7.vs.0s",
  "init_10.cxcr7.vs.0s",
  "init_600.cxcr7.vs.0s",
  "init_1800.cxcr7.vs.0s"
)



target_conditions <- c(
  "val_10.dmso.vs.0s",
  "val_600.dmso.vs.0s",
  "val_1800.dmso.vs.0s"
)



# Filter up_long and down_long
up_long_sub <- up_long %>% filter(time %in% target_conditions)
down_long_sub <- down_long %>% filter(time %in% target_conditions)

#############################################################
## 6.2 Build mirrored combined frame
#############################################################

up_long2 <- up_long_sub %>% 
  mutate(value = neg_log10_p,
         regulation = "up")

down_long2 <- down_long_sub %>% 
  mutate(value = -neg_log10_p,
         regulation = "down")

combined <- bind_rows(up_long2, down_long2)

# Pathway order
combined <- combined %>%
  mutate(pathway = factor(pathway, levels = selected_pathways))

# Time order: val first, then init
combined$time <- factor(
  combined$time,
  levels = target_conditions
)

#############################################################
## 6.3 Plot
#############################################################

colors6 <- c(
  "#a6cee3", "#1f78b4", "#08306b",   # val 10,600,1800
  "#fb9a99", "#e31a1c", "#99000d"    # init 10,600,1800
)

mirrored_plot <- ggplot(combined, aes(x = pathway, y = value, fill = time)) +
  geom_col(position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = colors6) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 5)) +
  labs(
    x = "Pathway",
    y = "-log10(p-value)",
    title = "Mirrored Up/Down Pathway Enrichment (CXCR7 vs 0s)"
  ) +
  theme_cowplot() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.title  = element_text(size = 10),
    plot.title  = element_text(size = 12, face = "bold")
  )

mirrored_plot




# Publication-ready mirrored pathway enrichment plot
library(ggplot2)
library(cowplot)
library(stringr)

# Color palette (blues for validation, reds for initial)
colors6 <- c(
  "#a6cee3", "#1f78b4", "#08306b",   # val 10,600,1800
  "#fb9a99", "#e31a1c", "#99000d"    # init 10,600,1800
)

# Create publication-ready plot with larger text
mirrored_plot <- ggplot(combined, aes(x = pathway, y = value, fill = time)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(
    values = colors6,
    name = "Time Point"
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) +
  labs(
    x = "Pathway",
    y = "-log10(p-value)",
    title = "Mirrored Up/Down Pathway Enrichment (CXCR7 vs 0s)"
  ) +
  theme_cowplot(font_size = 14) +
  theme(
    # Axis text
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    
    # Axis titles
    axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
    
    # Plot title
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    
    # Legend
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    
    # Panel and grid
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    
    # Margins
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

# Display the plot
print(mirrored_plot)

# Save as high-resolution TIFF (publication standard)
ggsave(
  filename = "D:/Research/CXCR7_platelet_analysis/figures/Supplementary_Figures/mirrored_pathway_enrichment.tiff",
  plot = mirrored_plot,
  width = 10,
  height = 7,
  dpi = 300,
  compression = "lzw"
)

# Also save as PDF (vector format, scalable)
ggsave(
  filename = "D:/Research/CXCR7_platelet_analysis/figures/Supplementary_Figures/mirrored_pathway_enrichment.pdf",
  plot = mirrored_plot,
  width = 10,
  height = 7,
  device = "pdf"
)

# Optional: save as PNG for presentations
ggsave(
  filename = "D:/Research/CXCR7_platelet_analysis/figures/Supplementary_Figures/mirrored_pathway_enrichment.png",
  plot = mirrored_plot,
  width = 10,
  height = 7,
  dpi = 300
)

cat("Plots saved successfully to Supplementary_Figures folder!\n")
cat("Formats: TIFF (publication), PDF (vector), PNG (presentation)\n")













###############################################################
# 7. PHOSPHOSITE-LEVEL HEATMAP ANALYSIS - Extract relevant phosphoproteins and psites
###############################################################

## Create pathway numbering lookup table
pathway_numbers <- data.frame(
  pathway = selected_pathways,
  pathway_num = 1:length(selected_pathways)
)


cat("Starting phosphosite-level clustering analysis...\n")

## ------------------------------------------------------------
## A) Extract genes from selected pathways with numbering
## ------------------------------------------------------------

# Get substrate lists from enrichment results
pathway_genes <- lapply(names(up_results), function(time_pt) {
  up_results[[time_pt]] %>%
    filter(pathway %in% selected_pathways) %>%
    select(pathway, substrates)
}) %>% 
  bind_rows() %>%
  distinct(pathway, substrates)

# Expand substrates (semicolon-separated) to individual genes
gene_pathway_map <- pathway_genes %>%
  separate_rows(substrates, sep = ";") %>%
  mutate(gene = trimws(substrates)) %>%
  left_join(pathway_numbers, by = "pathway") %>%
  select(gene, pathway_num) %>%
  distinct()

# Create pathway membership as numbers (e.g., "1, 5, 12")
gene_pathway_summary <- gene_pathway_map %>%
  group_by(gene) %>%
  summarise(pathway_nums = paste(sort(unique(pathway_num)), collapse = ", "), 
            .groups = "drop")

cat("  Extracted", nrow(gene_pathway_summary), "genes from selected pathways\n")

## ------------------------------------------------------------
## B) Build phosphosite-level matrix from validation data
## ------------------------------------------------------------

# Use validation CXCR7 vs 0s data (indices 7-9 in dfs_new_harm)
phosphosite_data <- lapply(names(all_datasets_filt)[4:6], function(time_name) {
  time_pt <- str_extract(time_name, "\\d+")
  
  all_datasets_filt[[time_name]] %>%
    filter(toupper(name) %in% toupper(gene_pathway_summary$gene)) %>%
    mutate(
      phosphosite_id = paste(name, PSite, sep = "_"),
      timepoint = time_pt
    ) %>%
    select(phosphosite_id, name, PSite, logFC, PValue, timepoint)
}) %>% bind_rows()

# Convert to wide format
phospho_wide <- phosphosite_data %>%
  pivot_wider(
    id_cols = c(phosphosite_id, name),
    names_from = timepoint,
    values_from = c(logFC, PValue),
    names_sep = "."
  )

# Merge with pathway numbering
phospho_annotated <- phospho_wide %>%
  left_join(gene_pathway_summary, by = c("name" = "gene"))

cat("  Created matrix with", nrow(phospho_annotated), "phosphosites\n")

## ------------------------------------------------------------
## C) Create master heatmap matrix
## ------------------------------------------------------------

# Extract logFC columns and create matrix
logfc_cols <- grep("^logFC\\.", colnames(phospho_annotated), value = TRUE)
logFC_mat <- as.matrix(phospho_annotated[, logfc_cols])
rownames(logFC_mat) <- phospho_annotated$phosphosite_id
colnames(logFC_mat) <- c("10s", "600s", "1800s")

# Remove rows with any NAs
logFC_mat <- logFC_mat[complete.cases(logFC_mat), ]

# Z-score normalization (row-wise)
logFC_z <- t(scale(t(logFC_mat)))

cat("  Z-scored matrix:", nrow(logFC_z), "sites ×", ncol(logFC_z), "timepoints\n")

## ------------------------------------------------------------
## D) Create significance star matrix
## ------------------------------------------------------------

create_sig_matrix <- function(phospho_data, logFC_matrix) {
  sig_mat <- matrix("", 
                    nrow = nrow(logFC_matrix), 
                    ncol = ncol(logFC_matrix),
                    dimnames = dimnames(logFC_matrix))
  
  for (i in 1:nrow(sig_mat)) {
    site_id <- rownames(logFC_matrix)[i]
    site_data <- phospho_data[phospho_data$phosphosite_id == site_id, ]
    
    if (nrow(site_data) > 0) {
      if (!is.na(site_data$PValue.10) && site_data$PValue.10 < 0.05) sig_mat[i, "10s"] <- "*"
      if (!is.na(site_data$PValue.600) && site_data$PValue.600 < 0.05) sig_mat[i, "600s"] <- "*"
      if (!is.na(site_data$PValue.1800) && site_data$PValue.1800 < 0.05) sig_mat[i, "1800s"] <- "*"
    }
  }
  sig_mat
}

sig_stars <- create_sig_matrix(phospho_annotated, logFC_z)

## ------------------------------------------------------------
## E) Master heatmap - all phosphosites
## ------------------------------------------------------------

cat("\nCreating master heatmap...\n")

pheatmap(
  logFC_z,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  display_numbers = sig_stars,
  number_color = "black",
  cellwidth = 20,   # Make cells more square
  cellheight = 3,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Master Heatmap: All Phosphosites (Z-scored logFC)",
  fontsize = 12,
  fontsize_number = 10
)

## ------------------------------------------------------------
## F) Hierarchical clustering into subclusters
## ------------------------------------------------------------

cat("\nPerforming hierarchical clustering...\n")

dist_rows <- dist(logFC_z)
hc_rows <- hclust(dist_rows, method = "complete")
subclusters <- cutree(hc_rows, k = 4)

cluster_counts <- table(subclusters)
cat("  Cluster sizes:", paste(names(cluster_counts), "=", cluster_counts, collapse = ", "), "\n")

###############################################################
# 8. SUBCLUSTER VISUALIZATION WITH PATHWAY NUMBERS ----
###############################################################

## Helper function to plot subcluster with compact labeling
plot_subcluster <- function(cluster_num, logFC_matrix, sig_matrix, 
                            phospho_data, subclusters, 
                            top_n = NULL, sort_by = "none",
                            save_pdf = FALSE, out_dir = NULL) {
  
  # Get rows in this cluster
  cluster_rows <- names(subclusters[subclusters == cluster_num])
  sub_mat <- logFC_matrix[cluster_rows, , drop = FALSE]
  sub_data <- phospho_data[phospho_data$phosphosite_id %in% cluster_rows, ]
  
  # Sorting options
  if (!is.null(top_n) && sort_by != "none") {
    if (sort_by == "pvalue") {
      sub_data$min_pval <- pmin(sub_data$PValue.10, sub_data$PValue.600, 
                                sub_data$PValue.1800, na.rm = TRUE)
      sub_data <- sub_data %>% arrange(min_pval) %>% head(top_n)
    } else if (sort_by == "logfc") {
      sub_data$max_absFC <- pmax(abs(sub_data$logFC.10), abs(sub_data$logFC.600),
                                 abs(sub_data$logFC.1800), na.rm = TRUE)
      sub_data <- sub_data %>% arrange(desc(max_absFC)) %>% head(top_n)
    }
    
    sub_mat <- sub_mat[match(sub_data$phosphosite_id, rownames(sub_mat)), , drop = FALSE]
  }
  
  # Create compact row labels: "GENE_SITE | pathway_nums"
  row_labels <- paste0(sub_data$phosphosite_id, " | ", sub_data$pathway_nums)
  
  # Create significance stars for subset
  sub_sig <- sig_matrix[rownames(sub_mat), , drop = FALSE]
  
  # Title
  title_text <- paste0("Cluster ", cluster_num, " (n=", nrow(sub_mat), ")")
  if (!is.null(top_n)) {
    title_text <- paste0(title_text, " - Top ", top_n, " by ", 
                         ifelse(sort_by == "pvalue", "P-value", "|logFC|"))
  }
  
  # Calculate dimensions for square-ish cells
  cell_width <- 20
  cell_height <- 12
  plot_width <- (ncol(sub_mat) * cell_width + 150) / 25.4  # Convert mm to inches
  plot_height <- (nrow(sub_mat) * cell_height + 100) / 25.4
  
  # Constrain dimensions
  plot_width <- max(6, min(plot_width, 12))
  plot_height <- max(4, min(plot_height, 20))
  
  # Plot
  if (save_pdf && !is.null(out_dir)) {
    suffix <- ifelse(!is.null(top_n), 
                     paste0("_top", top_n, "_by_", sort_by), 
                     "_full")
    pdf(file.path(out_dir, paste0("Subcluster_", cluster_num, suffix, ".pdf")), 
        width = plot_width, height = plot_height)
  }
  
  pheatmap(
    sub_mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    display_numbers = sub_sig,
    number_color = "black",
    labels_row = row_labels,
    fontsize_row = 7,
    fontsize_col = 10,
    fontsize_number = 8,
    cellwidth = cell_width,
    cellheight = cell_height,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    main = title_text
  )
  
  if (save_pdf && !is.null(out_dir)) {
    dev.off()
    cat("  Saved:", basename(file.path(out_dir, paste0("Subcluster_", cluster_num, suffix, ".pdf"))), "\n")
  }
}

## ------------------------------------------------------------
## Generate all subcluster plots
## ------------------------------------------------------------

cat("\nGenerating subcluster heatmaps...\n")

# Interactive display: top 20 by p-value for each cluster
for (i in 1:4) {
  plot_subcluster(i, logFC_z, sig_stars, phospho_annotated, subclusters,
                  top_n = 20, sort_by = "pvalue")
}

## Save PDFs
cat("\nSaving subcluster PDFs...\n")
out_dir <- file.path(fig_dir, "Subclusters")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

for (i in 1:4) {
  # Full cluster
  plot_subcluster(i, logFC_z, sig_stars, phospho_annotated, subclusters,
                  save_pdf = TRUE, out_dir = out_dir)
  
  # Top 20 by p-value
  plot_subcluster(i, logFC_z, sig_stars, phospho_annotated, subclusters,
                  top_n = 20, sort_by = "pvalue", save_pdf = TRUE, out_dir = out_dir)
  
  # Top 20 by logFC
  plot_subcluster(i, logFC_z, sig_stars, phospho_annotated, subclusters,
                  top_n = 20, sort_by = "logfc", save_pdf = TRUE, out_dir = out_dir)
}

cat("\n✓ Phosphosite clustering analysis complete!\n")
cat("  PDFs saved to:", out_dir, "\n")
cat("  - Full subclusters (4 PDFs)\n")
cat("  - Top 20 by p-value (4 PDFs)\n")
cat("  - Top 20 by |logFC| (4 PDFs)\n\n")

cat(rep("=", 60), "\n", sep = "")
cat("PIPELINE COMPLETE\n")
cat(rep("=", 60), "\n", sep = "")











































###############################################################
# 10. PHOSPHOSITE-LEVEL HEATMAP ANALYSIS ----
###############################################################

fig_dir <- "D:/Research/CXCR7_platelet_analysis/figures/Supplementary_Figures"

cat("Starting phosphosite-level clustering analysis...\n")

## ------------------------------------------------------------
## A) Extract genes from selected pathways
## ------------------------------------------------------------

# Get substrate lists from enrichment results (using first available result)
pathway_genes <- lapply(names(up_results), function(time_pt) {
  up_results[[time_pt]] %>%
    filter(pathway %in% selected_pathways) %>%
    select(pathway, substrates)
}) %>% 
  bind_rows() %>%
  distinct(pathway, substrates)

# Expand substrates (semicolon-separated) to individual genes
gene_pathway_map <- pathway_genes %>%
  separate_rows(substrates, sep = ";") %>%
  mutate(
    gene = trimws(substrates),
    pathway_clean = gsub("\\.", " ", pathway)
  ) %>%
  select(gene, pathway_clean) %>%
  distinct()

# Create pathway membership string for each gene
gene_pathway_summary <- gene_pathway_map %>%
  group_by(gene) %>%
  summarise(pathway_membership = paste(unique(pathway_clean), collapse = "; "), 
            .groups = "drop")

cat("  Extracted", nrow(gene_pathway_summary), "genes from selected pathways\n")

## ------------------------------------------------------------
## B) Build phosphosite-level matrix from NON-collapsed data
## ------------------------------------------------------------

# Use the original harmonized data (before collapse) for phosphosite-level info
# Filter to genes in selected pathways
phosphosite_data <- lapply(names(dfs_new_intersect)[4:6], function(time_name) {
  # Extract timepoint number from name (e.g., "val_10.cxcr7.vs.0s" -> "10")
  time_pt <- str_extract(time_name, "\\d+")
  
  dfs_new_intersect[[time_name]] %>%
    filter(toupper(name) %in% toupper(gene_pathway_summary$gene)) %>%
    mutate(
      phosphosite_id = paste(name, PSite, sep = "_"),
      timepoint = time_pt
    ) %>%
    select(phosphosite_id, name, PSite, logFC, PValue, timepoint)
}) %>% bind_rows()

# Convert to wide format
phospho_wide <- phosphosite_data %>%
  pivot_wider(
    id_cols = c(phosphosite_id, name),
    names_from = timepoint,
    values_from = c(logFC, PValue),
    names_sep = "."
  )

# Merge with pathway membership
phospho_annotated <- phospho_wide %>%
  left_join(gene_pathway_summary, by = c("name" = "gene"))

cat("  Created matrix with", nrow(phospho_annotated), "phosphosites\n")

## ------------------------------------------------------------
## C) Create master heatmap matrix
## ------------------------------------------------------------

# Extract logFC columns and create matrix
logfc_cols <- grep("^logFC\\.", colnames(phospho_annotated), value = TRUE)
logFC_mat <- as.matrix(phospho_annotated[, logfc_cols])
rownames(logFC_mat) <- phospho_annotated$phosphosite_id
colnames(logFC_mat) <- paste0(gsub("logFC\\.", "", colnames(logFC_mat)), "s")

# Remove rows with all NAs
logFC_mat <- logFC_mat[complete.cases(logFC_mat), ]

# Z-score normalization (row-wise)
logFC_z <- t(scale(t(logFC_mat)))

cat("  Z-scored matrix:", nrow(logFC_z), "sites ×", ncol(logFC_z), "timepoints\n")

## ------------------------------------------------------------
## D) Create significance star matrix
## ------------------------------------------------------------

create_sig_matrix <- function(phospho_data, logFC_matrix) {
  sig_mat <- matrix("", 
                    nrow = nrow(logFC_matrix), 
                    ncol = ncol(logFC_matrix),
                    dimnames = dimnames(logFC_matrix))
  
  for (i in 1:nrow(sig_mat)) {
    site_id <- rownames(logFC_matrix)[i]
    site_data <- phospho_data[phospho_data$phosphosite_id == site_id, ]
    
    if (nrow(site_data) > 0) {
      if (!is.na(site_data$PValue.10) && site_data$PValue.10 < 0.05) sig_mat[i, "10s"] <- "*"
      if (!is.na(site_data$PValue.600) && site_data$PValue.600 < 0.05) sig_mat[i, "600s"] <- "*"
      if (!is.na(site_data$PValue.1800) && site_data$PValue.1800 < 0.05) sig_mat[i, "1800s"] <- "*"
    }
  }
  sig_mat
}

sig_stars <- create_sig_matrix(phospho_annotated, logFC_z)

## ------------------------------------------------------------
## E) Master heatmap - all phosphosites
## ------------------------------------------------------------

cat("\nCreating master heatmap...\n")

pheatmap(
  logFC_z,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  display_numbers = sig_stars,
  number_color = "black",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Master Heatmap: All Phosphosites (Z-scored logFC)",
  fontsize = 12
)

## ------------------------------------------------------------
## F) Hierarchical clustering into subclusters
## ------------------------------------------------------------

cat("\nPerforming hierarchical clustering...\n")

# Compute distance and cluster
dist_rows <- dist(logFC_z)
hc_rows <- hclust(dist_rows, method = "complete")
subclusters <- cutree(hc_rows, k = 4)

# Count sites per cluster
cluster_counts <- table(subclusters)
cat("  Cluster sizes:", paste(names(cluster_counts), "=", cluster_counts, collapse = ", "), "\n")

###############################################################
# 11. SUBCLUSTER VISUALIZATION ----
###############################################################

## ------------------------------------------------------------
## Helper function to plot subcluster
## ------------------------------------------------------------

plot_subcluster <- function(cluster_num, logFC_matrix, sig_matrix, 
                            phospho_data, subclusters, 
                            top_n = NULL, sort_by = "none",
                            save_pdf = FALSE, out_dir = NULL) {
  
  # Get rows in this cluster
  cluster_rows <- names(subclusters[subclusters == cluster_num])
  sub_mat <- logFC_matrix[cluster_rows, , drop = FALSE]
  sub_data <- phospho_data[phospho_data$phosphosite_id %in% cluster_rows, ]
  
  # Sorting options
  if (!is.null(top_n) && sort_by != "none") {
    if (sort_by == "pvalue") {
      sub_data$min_pval <- pmin(sub_data$PValue.10, sub_data$PValue.600, 
                                sub_data$PValue.1800, na.rm = TRUE)
      sub_data <- sub_data %>% arrange(min_pval) %>% head(top_n)
    } else if (sort_by == "logfc") {
      sub_data$max_absFC <- pmax(abs(sub_data$logFC.10), abs(sub_data$logFC.600),
                                 abs(sub_data$logFC.1800), na.rm = TRUE)
      sub_data <- sub_data %>% arrange(desc(max_absFC)) %>% head(top_n)
    }
    
    # Reorder matrix
    sub_mat <- sub_mat[match(sub_data$phosphosite_id, rownames(sub_mat)), , drop = FALSE]
  }
  
  # Create row labels with pathway info
  row_labels <- paste(sub_data$phosphosite_id, 
                      substr(sub_data$pathway_membership, 1, 40), 
                      sep = " | ")
  
  # Create significance stars for subset
  sub_sig <- sig_matrix[rownames(sub_mat), , drop = FALSE]
  
  # Title
  title_text <- paste0("Cluster ", cluster_num, " (n=", nrow(sub_mat), ")")
  if (!is.null(top_n)) {
    title_text <- paste0(title_text, " - Top ", top_n, " by ", 
                         ifelse(sort_by == "pvalue", "P-value", "logFC"))
  }
  
  # Plot
  if (save_pdf && !is.null(out_dir)) {
    pdf(file.path(out_dir, paste0("Subcluster_", cluster_num, 
                                  ifelse(!is.null(top_n), paste0("_top", top_n), ""),
                                  ".pdf")), 
        width = 8, height = max(6, nrow(sub_mat) * 0.15))
  }
  
  pheatmap(
    sub_mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    display_numbers = sub_sig,
    number_color = "black",
    labels_row = row_labels,
    fontsize_row = 8,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    main = title_text
  )
  
  if (save_pdf && !is.null(out_dir)) {
    dev.off()
    cat("  Saved:", file.path(out_dir, paste0("Subcluster_", cluster_num, ".pdf")), "\n")
  }
}

## ------------------------------------------------------------
## A) Plot all subclusters (full)
## ------------------------------------------------------------

cat("\nGenerating subcluster heatmaps (all sites)...\n")
for (i in 1:4) {
  plot_subcluster(i, logFC_z, sig_stars, phospho_annotated, subclusters)
}

## ------------------------------------------------------------
## B) Plot top 20 by p-value
## ------------------------------------------------------------

cat("\nGenerating top 20 by p-value for each cluster...\n")
for (i in 1:4) {
  plot_subcluster(i, logFC_z, sig_stars, phospho_annotated, subclusters,
                  top_n = 20, sort_by = "pvalue")
}

## ------------------------------------------------------------
## C) Plot top 20 by logFC
## ------------------------------------------------------------

cat("\nGenerating top 20 by logFC for each cluster...\n")
for (i in 1:4) {
  plot_subcluster(i, logFC_z, sig_stars, phospho_annotated, subclusters,
                  top_n = 20, sort_by = "logfc")
}

## ------------------------------------------------------------
## D) Save all subclusters to PDF
## ------------------------------------------------------------

cat("\nSaving subcluster PDFs...\n")
out_dir <- file.path(fig_dir, "Subclusters")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

for (i in 1:4) {
  # Full cluster
  plot_subcluster(i, logFC_z, sig_stars, phospho_annotated, subclusters,
                  save_pdf = TRUE, out_dir = out_dir)
  
  # Top 20 by p-value
  plot_subcluster(i, logFC_z, sig_stars, phospho_annotated, subclusters,
                  top_n = 20, sort_by = "pvalue", save_pdf = TRUE, out_dir = out_dir)
  
  # Top 20 by logFC
  plot_subcluster(i, logFC_z, sig_stars, phospho_annotated, subclusters,
                  top_n = 20, sort_by = "logfc", save_pdf = TRUE, out_dir = out_dir)
}

cat("\n✓ Phosphosite clustering analysis complete!\n")
cat("  PDFs saved to:", out_dir, "\n")
cat("  - Full subclusters (4 PDFs)\n")
cat("  - Top 20 by p-value (4 PDFs)\n")
cat("  - Top 20 by logFC (4 PDFs)\n\n")

cat("="*60, "\n")
cat("PIPELINE COMPLETE\n")
cat("="*60, "\n")











