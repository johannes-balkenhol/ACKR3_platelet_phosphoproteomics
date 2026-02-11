###############################################################
## PHOSPHOPROTEOMICS PATHWAY ENRICHMENT PIPELINE
## NEW Validation Dataset (Default) + Optional OLD Initial Dataset
###############################################################

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
## 1) HARMONIZATION FUNCTIONS
###############################################################

cat("STEP 1: Setting up harmonization functions\n")
cat(strrep("─", 80), "\n\n")

## Harmonize NEW validation datasets
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

## Harmonize OLD initial datasets  
## structure from df$id: "UNIPROT;GENE;PSITE;PEPTIDE;INDEX"
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

cat("✓ harmonize_new() - for NEW validation datasets\n")
cat("✓ harmonize_old() - for OLD initial datasets\n\n")


###############################################################
## 2) LOAD & HARMONIZE NEW VALIDATION DATASETS
###############################################################

cat("\nSTEP 2: Loading NEW validation datasets\n")
cat(strrep("─", 80), "\n\n")


## select all comparisons
dfs_new <- list(
  top.10, top.600, top.1800,
  top.10.dmso.vs.0s, top.600.dmso.vs.0s, top.1800.dmso.vs.0s,
  top.10.cxcr7.vs.0s, top.600.cxcr7.vs.0s, top.1800.cxcr7.vs.0s
)

names(dfs_new) <- c(
  "val_10.dmso.vs.cxcr7", "val_600.dmso.vs.cxcr7", "val_1800.dmso.vs.cxcr7",
  "val_10.dmso.vs.0s",    "val_600.dmso.vs.0s",    "val_1800.dmso.vs.0s",
  "val_10.cxcr7.vs.0s",   "val_600.cxcr7.vs.0s",   "val_1800.cxcr7.vs.0s"
)


## select only the important comparison
dfs_new <- list(
  top.10, top.600, top.1800
)

names(dfs_new) <- c(
  "val_10.dmso.vs.cxcr7", "val_600.dmso.vs.cxcr7", "val_1800.dmso.vs.cxcr7"
)




dfs_new_raw <- lapply(dfs_new, harmonize_new)

cat(sprintf("✓ Loaded %d NEW validation datasets\n\n", length(dfs_new_raw)))
for (nm in names(dfs_new_raw)) {
  cat(sprintf("  • %s (%d phosphosites)\n", nm, nrow(dfs_new_raw[[nm]])))
}

###############################################################
## 3) OPTIONAL: LOAD & HARMONIZE OLD INITIAL DATASETS
###############################################################

cat("\n")
cat(strrep("▓", 80), "\n")
cat("OPTIONAL: Load OLD initial datasets?\n")
cat(strrep("▓", 80), "\n\n")

cat("Set: use_old_data <- TRUE  (to include OLD datasets)\n")
cat("     use_old_data <- FALSE (NEW datasets only) [DEFAULT]\n\n")

use_old_data <- FALSE  # ← CHANGE THIS TO TRUE IF YOU WANT OLD DATASETS

if (use_old_data) {
  
  cat("Loading OLD initial datasets...\n\n")
  
  old_path  <- "SubProjects/CXCR7_initial/phosphoproteomics/data/processed_data"
  old_files <- list.files(old_path, pattern = "cxcr7", full.names = TRUE)
  
  old_tables_raw        <- lapply(old_files, read.delim)
  old_tables_harmonized <- lapply(old_tables_raw, harmonize_old)
  
  names(old_tables_harmonized) <- c(
    "init_10.cxcr7.vs.0s","init_30.cxcr7.vs.0s","init_60.cxcr7.vs.0s",
    "init_300.cxcr7.vs.0s","init_600.cxcr7.vs.0s","init_900.cxcr7.vs.0s",
    "init_1800.cxcr7.vs.0s"
  )
  
  cat(sprintf("✓ Loaded %d OLD initial datasets\n\n", length(old_tables_harmonized)))
  for (nm in names(old_tables_harmonized)) {
    cat(sprintf("  • %s (%d phosphosites)\n", nm, nrow(old_tables_harmonized[[nm]])))
  }
  
  # Combine NEW + OLD
  all_inputs_raw <- c(dfs_new_raw, old_tables_harmonized)
  
} else {
  
  cat("✓ Using NEW validation datasets ONLY (recommended)\n\n")
  all_inputs_raw <- dfs_new_raw
  
}


###############################################################
## 4) COLLAPSE BY UNIPROT (CHOOSE TOP PHOSPHOSITE PER PROTEIN)
###############################################################


cat("\nSTEP 3: Collapsing phosphosites by protein\n")
cat(strrep("─", 80), "\n\n")

# ← CHOOSE COLLAPSE METHOD
collapse_by <- "individual"  # Options: "mean_logfc", "max_logfc", "pvalue", "individual"

## Combine dmso.vs.cxcr7 timepoints to find best phosphosite per protein
combined_timepoints <- bind_rows(
  dfs_new_raw$`val_10.dmso.vs.cxcr7` %>% mutate(timepoint = "10"),
  dfs_new_raw$`val_600.dmso.vs.cxcr7` %>% mutate(timepoint = "600"),
  dfs_new_raw$`val_1800.dmso.vs.cxcr7` %>% mutate(timepoint = "1800")
)

## ------------------------------------------------------------
## Conditional: Collapse globally or individually per timepoint
## ------------------------------------------------------------

if (collapse_by == "individual") {
  
  cat("✓ Mode: INDIVIDUAL COLLAPSE per timepoint\n")
  cat("  (Each timepoint picks its own best phosphosite)\n\n")
  
  ## Collapse function for individual datasets
  collapse_individual <- function(df, method = "pvalue") {
    if (method == "pvalue") {
      df %>%
        group_by(uniprot_id) %>%
        slice_min(PValue, n = 1, with_ties = FALSE) %>%
        ungroup()
    } else if (method == "max_logfc") {
      df %>%
        mutate(abs_logFC = abs(logFC)) %>%
        group_by(uniprot_id) %>%
        slice_max(abs_logFC, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        select(-abs_logFC)
    } else {  # mean_logfc - not applicable for single timepoint
      df %>%
        mutate(abs_logFC = abs(logFC)) %>%
        group_by(uniprot_id) %>%
        slice_max(abs_logFC, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        select(-abs_logFC)
    }
  }
  
  # Collapse each dataset individually
  all_inputs_collapsed <- lapply(dfs_new_raw, function(df) {
    collapse_individual(df, method = "pvalue") %>%
      arrange(uniprot_id) %>%
      select(uniprot_id, name, PSite, Average, logFC, PValue)
  })
  
  # Track which phosphosites were selected at each timepoint
  selection_tracking <- bind_rows(
    all_inputs_collapsed$`val_10.dmso.vs.cxcr7` %>% 
      select(uniprot_id, name, PSite) %>% 
      mutate(timepoint = "10s"),
    all_inputs_collapsed$`val_600.dmso.vs.cxcr7` %>% 
      select(uniprot_id, name, PSite) %>% 
      mutate(timepoint = "600s"),
    all_inputs_collapsed$`val_1800.dmso.vs.cxcr7` %>% 
      select(uniprot_id, name, PSite) %>% 
      mutate(timepoint = "1800s")
  )
  
  # Reshape to wide format to see differences
  selection_wide <- selection_tracking %>%
    pivot_wider(
      id_cols = c(uniprot_id, name),
      names_from = timepoint,
      values_from = PSite
    )
  
  # Save tracking
  write.csv(selection_wide, 
            "selected_phosphosites_individual_per_timepoint.csv", 
            row.names = FALSE)
  
  cat("✓ Method: Individual collapse (best p-value per timepoint)\n")
  cat("✓ Saved: selected_phosphosites_individual_per_timepoint.csv\n\n")
  
  for (nm in names(all_inputs_collapsed)) {
    cat(sprintf("  %s: %d proteins\n", nm, nrow(all_inputs_collapsed[[nm]])))
  }
  
  # Show key kinases - which sites selected at each timepoint?
  cat("\nKey kinases - phosphosite selection per timepoint:\n")
  cat(strrep("-", 80), "\n")
  
  for (kinase in c("SHC1", "SRC", "AKT1", "PRKCA", "PRKCD")) {
    kinase_selection <- selection_wide %>% filter(name == kinase)
    
    if (nrow(kinase_selection) > 0) {
      cat(sprintf("  %s:\n", kinase))
      cat(sprintf("    10s:   %s\n", 
                  ifelse(!is.na(kinase_selection$`10s`), kinase_selection$`10s`, "N/A")))
      cat(sprintf("    600s:  %s\n", 
                  ifelse(!is.na(kinase_selection$`600s`), kinase_selection$`600s`, "N/A")))
      cat(sprintf("    1800s: %s\n", 
                  ifelse(!is.na(kinase_selection$`1800s`), kinase_selection$`1800s`, "N/A")))
      
      # Check if different
      sites <- unique(c(kinase_selection$`10s`, 
                        kinase_selection$`600s`, 
                        kinase_selection$`1800s`))
      sites <- sites[!is.na(sites)]
      
      if (length(sites) > 1) {
        cat(sprintf("    ⚠️  DIFFERENT sites selected! (%s)\n", 
                    paste(sites, collapse = ", ")))
      } else {
        cat(sprintf("    ✓  Same site across all timepoints\n"))
      }
    }
  }
  
  # Count how many proteins have different sites
  different_sites <- selection_wide %>%
    rowwise() %>%
    mutate(n_unique = n_distinct(c(`10s`, `600s`, `1800s`), na.rm = TRUE)) %>%
    ungroup()
  
  cat("\n", strrep("-", 80), "\n")
  cat(sprintf("Proteins with SAME site across timepoints: %d\n", 
              sum(different_sites$n_unique == 1, na.rm = TRUE)))
  cat(sprintf("Proteins with DIFFERENT sites: %d\n", 
              sum(different_sites$n_unique > 1, na.rm = TRUE)))
  
} else {
  
  ## Global collapse (same phosphosite across all timepoints)
  best_psites <- combined_timepoints %>%
    group_by(uniprot_id, name, PSite) %>%
    summarise(
      mean_abs_logFC = mean(abs(logFC), na.rm = TRUE),
      max_abs_logFC = max(abs(logFC), na.rm = TRUE),
      min_pvalue = min(PValue, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(uniprot_id) %>%
    {
      if (collapse_by == "mean_logfc") {
        slice_max(., mean_abs_logFC, n = 1, with_ties = FALSE)
      } else if (collapse_by == "max_logfc") {
        slice_max(., max_abs_logFC, n = 1, with_ties = FALSE)
      } else {  # pvalue
        slice_min(., min_pvalue, n = 1, with_ties = FALSE)
      }
    } %>%
    ungroup() %>%
    select(uniprot_id, PSite)
  
  cat(sprintf("✓ Selected %d phosphosites (one per protein)\n", nrow(best_psites)))
  
  all_inputs_collapsed <- lapply(dfs_new_raw, function(df) {
    df %>%
      inner_join(best_psites, by = c("uniprot_id", "PSite")) %>%
      arrange(uniprot_id) %>%
      select(uniprot_id, name, PSite, Average, logFC, PValue)
  })
  
  method_name <- case_when(
    collapse_by == "mean_logfc" ~ "Mean |logFC| (consistent effects)",
    collapse_by == "max_logfc" ~ "Max |logFC| (peak effects)",
    collapse_by == "pvalue" ~ "Min p-value (significance)",
    TRUE ~ "Unknown"
  )
  
  cat(sprintf("✓ Method: %s\n", method_name))
  
  for (nm in names(all_inputs_collapsed)) {
    cat(sprintf("  %s: %d proteins\n", nm, nrow(all_inputs_collapsed[[nm]])))
  }
  
  best_psites_annotated <- best_psites %>%
    left_join(combined_timepoints %>% distinct(uniprot_id, name), by = "uniprot_id")
  
  write.csv(best_psites_annotated, 
            "selected_phosphosites_for_enrichment.csv", 
            row.names = FALSE)
  
  cat("\n✓ Saved: selected_phosphosites_for_enrichment.csv\n")
  cat("\n✓ Same phosphosite used across all timepoints!\n")
}

cat("\n")

###############################################################
## 4) BUILD LOG2FC MATRIX
###############################################################


cat("\nSTEP 4: Building logFC matrix\n")
cat(strrep("─", 80), "\n\n")

logfc_matrix <- do.call(cbind, lapply(all_inputs_collapsed, `[[`, "logFC"))
colnames(logfc_matrix) <- names(all_inputs_collapsed)
rownames(logfc_matrix) <- all_inputs_collapsed[[1]]$uniprot_id

gene_symbols <- all_inputs_collapsed[[1]]$name

cat(sprintf("✓ Matrix: %d proteins × %d comparisons\n", 
            nrow(logfc_matrix), ncol(logfc_matrix)))
cat("✓ Rownames: UniProt IDs\n")
cat("✓ Gene symbols: stored separately\n\n")

cat("Matrix preview:\n")
print(head(logfc_matrix[, 1:3]))

cat("\n")


###############################################################
## 5) LOAD REACTOME PATHWAYS & MATCH TO CURATED LIST
###############################################################

cat("\n")
cat(strrep("▓", 80), "\n")
cat("STEP 5: Loading Reactome pathways & curated selection\n")
cat(strrep("▓", 80), "\n\n")

# Load Reactome pathways (ONCE)
pathways_all <- as.list(reactomePATHID2EXTID)
path_names <- as.list(reactomePATHID2NAME)
name_id <- match(names(pathways_all), names(path_names))
names(pathways_all) <- unlist(path_names)[name_id]

# Keep only HOMO SAPIENS
pathways_all <- pathways_all[grepl("Homo sapiens", names(pathways_all), ignore.case = TRUE)]

# Convert Entrez IDs → Gene symbols
pathways_all <- lapply(pathways_all, function(path) {
  gene_name <- unname(getSYMBOL(path, data = "org.Hs.eg"))
  toupper(unique(gene_name))
})

# Clean pathway names
pathway_names_clean <- names(pathways_all)
pathway_names_clean <- gsub("Homo sapiens: ", "", pathway_names_clean, ignore.case = TRUE)
pathway_names_clean <- gsub("Homo sapiens >> ", "", pathway_names_clean, ignore.case = TRUE)
pathway_names_clean <- trimws(pathway_names_clean)
names(pathways_all) <- pathway_names_clean

cat(sprintf("✓ Total Reactome pathways (Homo sapiens): %d\n\n", length(pathways_all)))

###############################################################
## 6) MATCH CURATED PATHWAYS TO REACTOME DATABASE
###############################################################

cat("Matching curated pathways to Reactome database...\n\n")

manual_path_refined <- c(
  # Platelet pathways
  "Platelet activation",
  "Platelet degranulation",
  "Response to elevated platelet cytosolic Ca2+",
  "Platelet Adhesion to exposed collagen",
  "Hemostasis",
  "Platelet calcium homeostasis",
  "Platelet Aggregation (Plug Formation)",
  
  # Clathrin & Endocytosis
  "Clathrin-mediated endocytosis",
  "Cargo recognition for clathrin-mediated endocytosis",
  
  # Integrin signaling
  "Integrin cell surface interactions",
  "Integrin signaling",
  
  # RTK signaling
  "Signaling by Receptor Tyrosine Kinases",
  #"SHC1 events in EGFR signaling",                    # ADDED - from original figure
  #"SHC1 events in ERBB2 signaling",                   # ADDED - from original figure
  #"SHC1 events in ERBB4 signaling",                   # ADDED - from original figure
  
  # GPCR pathways
  "Signaling by GPCR",
  "GPCR downstream signalling",
  "G beta:gamma signalling through PI3Kgamma",
  "Activation of G protein gated Potassium channels",
  "G beta:gamma signalling through PLC beta",
  "G beta:gamma signalling through BTK",
  "G beta:gamma signalling through CDC42",
  "G protein gated Potassium channels",
  "GPVI-mediated activation cascade",
  
  # FCGR
  "FCGR activation",
  "Fcgamma receptor (FCGR) dependent phagocytosis",
  
  # PKA pathways
  "PKA activation",
  "PKA activation in glucagon signalling",
  "Adenylate cyclase activating pathway",
  "Adenylate cyclase inhibitory pathway",
  "CREB1 phosphorylation through the activation of Adenylate Cyclase",
  "PKA-mediated phosphorylation of CREB",
  "PKA-mediated phosphorylation of key metabolic factors",
  "Gastrin-CREB signalling pathway via PKC and MAPK",
  
  # cGMP/PKG
  "cGMP effects",
  "Nitric oxide stimulates guanylate cyclase",
  
  # AMPK
  "Energy dependent regulation of mTOR by LKB1-AMPK",
  "Activation of AMPK downstream of NMDARs",
  "AMPK inhibits chREBP transcriptional activation activity",
  
  # MAPK/ERK pathways
  "RAF activation",
  "Signalling to ERKs",                               # Already included
  "ERK/MAPK targets",
  "ERKs are inactivated",
  "Signaling by BRAF and RAF1 fusions",
  "Negative regulation of MAPK pathway",
  "GRB2:SOS provides linkage to MAPK signaling for Integrins",
  "Signal attenuation",                               # ADDED - from original figure
  
  # RAS pathways
  "Regulation of RAS by GAPs",
  
  # RHO GTPases
  "Signaling by Rho GTPases",
  "RHOA GTPase cycle",
  #"RHOB GTPase cycle",
  #"RHOC GTPase cycle",
  #"RHOV GTPase cycle",
  "RHO GTPase cycle",
  "RHO GTPase Effectors",
  "RHO GTPases Activate WASPs and WAVEs",
  "RHO GTPases Activate NADPH Oxidases",
  
  # PI3K/AKT/mTOR
  "PIP3 activates AKT signaling",
  "AKT phosphorylates targets in the cytosol",
  "Negative regulation of the PI3K/AKT network",      # Already included
  "PI3K Cascade",
  "Effects of PIP2 hydrolysis",
  "Ca-dependent events",
  "MTOR signalling",
  "mTORC1-mediated signalling",
  "Amino acids regulate mTORC1",
  
  # Autophagy
  "Autophagy",
  
  # Non-receptor TKs
  "Signaling by Non-Receptor Tyrosine Kinases",
  "DAP12 signaling",
  
  # ECM
  "Extracellular matrix organization",
  "Degradation of the extracellular matrix",
  
  # Trafficking
  "Membrane Trafficking",
  
  # Lipid signaling
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
  "Role of phospholipids in phagocytosis",
  
  # Core platelet function (MUST ADD)
  "Platelet activation, signaling and aggregation",
  "Platelet homeostasis",
  "P2Y receptors",
  "ADP signalling through P2Y purinoceptor 12",
  "ADP signalling through P2Y purinoceptor 1",
  "Thrombin signalling through proteinase activated receptors (PARs)",
  
  # Platelet adhesion to collagen (IMPORTANT)
  "Defects of platelet adhesion to exposed collagen",
  
  # Collagen pathways (relevant for platelet adhesion)
  "Collagen degradation",
  "Collagen formation",
  "Collagen biosynthesis and modifying enzymes",
  
  # PKC/Phospholipase C (MISSING from your list!)
  "PLCG1 events in ERBB2 signaling",
  "PLC-gamma1 signalling",
  "EGFR interacts with phospholipase C-gamma",
  "Phospholipase C-mediated cascade: FGFR1",
  
  # Additional platelet-relevant
  "Platelet sensitization by LDL",
  #"RUNX1 regulates genes involved in megakaryocyte differentiation and platelet function",
  #"Factors involved in megakaryocyte development and platelet production",
  
  # Prostanoid receptors (IMPORTANT - includes TXA2 receptor)
  "Prostanoid ligand receptors",
  
  # Serotonin receptors (dense granule release, platelet aggregation)
  "Serotonin receptors"
)



# Use your curated pathway list
#selected_pathways <- unique(manual_filter)
selected_pathways <- unique(manual_path_refined) 


# Matching function
match_pathway <- function(query, all_names) {
  # Try exact match first
  idx <- which(tolower(all_names) == tolower(query))
  if (length(idx) > 0) return(all_names[idx[1]])
  
  # Try partial match if no exact match
  idx <- grep(tolower(gsub(" ", ".*", query)), tolower(all_names))
  if (length(idx) > 0) return(all_names[idx[1]])
  
  return(NA)
}

# Match all 80 curated pathways
pathways <- list()
matched_count <- 0
not_found <- character()

for (path in manual_path_refined) {
  found_name <- match_pathway(path, pathway_names_clean)
  
  if (!is.na(found_name)) {
    pathways[[path]] <- pathways_all[[found_name]]
    cat("✓", path, "\n")
    matched_count <- matched_count + 1
  } else {
    cat("✗", path, "\n")
    not_found <- c(not_found, path)
  }
}

cat("\n", strrep("=", 80), "\n")
cat("PATHWAY MATCHING SUMMARY\n")
cat(strrep("=", 80), "\n\n")

cat(sprintf("Total curated pathways: %d\n", length(manual_path_refined)))
cat(sprintf("Successfully matched:  %d\n", matched_count))
cat(sprintf("NOT found:             %d\n", length(not_found)))
cat(sprintf("Match rate:            %.1f%%\n\n", 100 * matched_count / length(manual_path_refined)))

if (length(not_found) > 0) {
  cat("NOT FOUND:\n")
  for (p in not_found) cat(" -", p, "\n")
}

cat(sprintf("\n✓ Ready for enrichment with %d pathways!\n\n", length(pathways)))

# Save reference
selected_pathways_df <- data.frame(
  pathway_number = seq_along(pathways),
  pathway_name = names(pathways),
  n_genes = sapply(pathways, length)
)

write.csv(selected_pathways_df, 
          "selected_pathways_matched.csv", 
          row.names = FALSE)

cat("✓ Saved reference: selected_pathways_matched.csv\n\n")

Tc.gene <- logfc_matrix
rownames(Tc.gene) <- toupper(gene_symbols)  # ← FIX: Use gene symbols, uppercase
names_input <- colnames(Tc.gene)

cat(strrep("█", 150), "\n")
cat("✓ READY FOR PATHWAY ENRICHMENT\n")
cat(strrep("█", 150), "\n\n")




#############################
## 9) ENRICHMENT WITH all PATHWAYS
#############################

cat("\n", strrep("=", 70), "\n")
cat("RUNNING PATHWAY ENRICHMENT\n")
cat(strrep("=", 70), "\n\n")

up_results <- list()
down_results <- list()

# Function to process enrichment results
process_enrichment <- function(res, selected_pathways) {
  # Convert to data frame
  path3 <- as.data.frame(res)
  path3$pathway <- rownames(path3)
  
  # Create pathway sizes
  pathways_df <- data.frame(
    pathway = names(selected_pathways),
    pw.size = sapply(selected_pathways, length),
    stringsAsFactors = FALSE
  )
  
  # Merge
  path4 <- merge(path3, pathways_df, by = "pathway", all = FALSE)
  
  # Clean columns
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

# ============================================================
# UP-REGULATED PATHWAYS
# ============================================================
cat("UP-REGULATED (alter='greater'):\n")
cat(strrep("─", 70), "\n\n")

for (i in seq_len(ncol(Tc.gene))) {
  res <- pathwayRankBasedEnrichment(Tc.gene[, i], annotation = pathways, alter = "greater")
  path4 <- process_enrichment(res, pathways)
  
  up_results[[names_input[i]]] <- path4
  
  sig <- sum(path4$pvalue < 0.05)
  cat(sprintf("✓ %s: %d pathways, %d significant (p<0.05)\n", 
              names_input[i], nrow(path4), sig))
}

# ============================================================
# DOWN-REGULATED PATHWAYS
# ============================================================
cat("\n")
cat("DOWN-REGULATED (alter='less'):\n")
cat(strrep("─", 70), "\n\n")

for (i in seq_len(ncol(Tc.gene))) {
  res <- pathwayRankBasedEnrichment(Tc.gene[, i], annotation = pathways, alter = "less")
  path4 <- process_enrichment(res, pathways)
  
  down_results[[names_input[i]]] <- path4
  
  sig <- sum(path4$pvalue < 0.05)
  cat(sprintf("✓ %s: %d pathways, %d significant (p<0.05)\n", 
              names_input[i], nrow(path4), sig))
}



#############################
## 10) BUILD SUMMARIES
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
## 11. Save Enrichment Results
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
## PARSE COMPARISONS (time, batch, comparison)
#############################

cat("\n", strrep("=", 70), "\n")
cat("PARSING PATHWAY COMPARISONS\n")
cat(strrep("=", 70), "\n\n")


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
    timepoint_sec = as.numeric(sub("\\..*", "", rest)),  # first digits
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
    timepoint_sec = as.numeric(sub("\\..*", "", rest)),
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
print(head(up_summary_long_filtered %>% select(pathway, pvalue, batch, timepoint_sec, comparison), n = 40))

cat("\n", strrep("=", 70), "\n")
cat("Sample DOWN-regulated:\n")
print(head(down_summary_long_filtered %>% select(pathway, pvalue, batch, timepoint_sec, comparison), n = 40))
cat("\n")





###############################################################
## STEP 13.3: Create enrichment plot (UPDATED FUNCTION)
###############################################################
cat("\n", strrep("=", 80), "\n")
cat("CLEAN RESTART: STEPS 13 & 14\n")
cat(strrep("=", 80), "\n\n")

# Remove old pathway_heatmaps directory to avoid confusion
if (dir.exists("pathway_heatmaps")) {
  cat("Removing old pathway_heatmaps directory...\n")
  unlink("pathway_heatmaps", recursive = TRUE)
  cat("✓ Old files removed\n\n")
}


plot_pathway_enrichment <- function(
    batch = "val",
    comparison = "dmso.vs.cxcr7",
    n_pathways = 10,
    sort_by = "pvalue",
    save_plot = TRUE,
    output_format = c("pdf", "png")
) {
  
  cat("\n", strrep("▓", 80), "\n")
  cat(sprintf("STEP 13: Plotting enrichment - %s | %s | Top %d pathways\n", 
              toupper(batch), comparison, n_pathways))
  cat(strrep("▓", 80), "\n\n")
  
  # Validate inputs
  if (!batch %in% c("val", "init")) {
    stop("batch must be 'val' or 'init'")
  }
  if (!comparison %in% c("cxcr7.vs.0s", "dmso.vs.cxcr7", "dmso.vs.0s")) {
    stop("comparison must be one of: cxcr7.vs.0s, dmso.vs.cxcr7, dmso.vs.0s")
  }
  
  # Combine UP and DOWN data (fix scoping bug with !!)
  all_data_filtered <- bind_rows(
    up_summary_long_filtered %>% 
      filter(batch == !!batch, comparison == !!comparison) %>%
      mutate(regulation = "up") %>%
      group_by(pathway, timepoint_sec, regulation) %>%
      slice_min(pvalue, n = 1, with_ties = FALSE) %>%
      ungroup(),
    down_summary_long_filtered %>%
      filter(batch == !!batch, comparison == !!comparison) %>%
      mutate(regulation = "down") %>%
      group_by(pathway, timepoint_sec, regulation) %>%
      slice_min(pvalue, n = 1, with_ties = FALSE) %>%
      ungroup()
  )
  
  # Rank pathways
  pathway_ranking <- all_data_filtered %>%
    filter(!is.na(pvalue), !is.infinite(pvalue)) %>%
    group_by(pathway) %>%
    summarize(
      min_pvalue = min(pvalue, na.rm = TRUE),
      max_neg_log10_p = max(abs(neg_log10_p), na.rm = TRUE),
      best_timepoint = timepoint_sec[which.min(pvalue)],
      best_regulation = regulation[which.min(pvalue)],
      .groups = 'drop'
    ) %>%
    arrange(min_pvalue) %>%
    mutate(pathway_num = row_number())
  
  # Select top N
  top_pathways_master <- pathway_ranking %>%
    head(n_pathways)
  
  top_paths <- top_pathways_master$pathway
  pathway_order <- top_paths
  
  cat("\n📊 TOP PATHWAYS BY MINIMUM P-VALUE:\n")
  cat(strrep("-", 80), "\n")
  ranking_display <- top_pathways_master %>% 
    select(pathway_num, pathway, min_pvalue, best_timepoint, best_regulation)
  print(ranking_display, n = Inf)
  cat("\n")
  
  # Filter for plotting
  combined_data <- all_data_filtered %>%
    filter(pathway %in% top_paths) %>%
    mutate(
      value = ifelse(regulation == "up", neg_log10_p, -neg_log10_p),
      pathway = factor(pathway, levels = pathway_order),
      dodge_group = factor(
        paste0(timepoint_sec),
        levels = c("10", "600", "1800")
      )
    ) %>%
    filter(!is.na(value), !is.infinite(value))
  
  if (nrow(combined_data) == 0) {
    cat("⚠️  No data found!\n")
    return(NULL)
  }
  
  # Y-axis limits
  max_abs_value <- max(abs(combined_data$value), na.rm = TRUE)
  y_limit <- ceiling(max_abs_value * 1.1)
  y_limit <- max(y_limit, 3)
  
  # Timepoint groups
  combined_data <- combined_data %>%
    mutate(
      timepoint_group = factor(
        paste0(timepoint_sec, " sec"),
        levels = c("10 sec", "600 sec", "1800 sec")
      )
    )
  
  # Colors
  colors_timepoint <- c(
    "10 sec" = "#E8D5F0",
    "600 sec" = "#B366D9",
    "1800 sec" = "#4B0082"
  )
  
  # Create plot
  plot <- ggplot(combined_data, aes(x = pathway, y = value, fill = timepoint_group, group = dodge_group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    scale_fill_manual(
      name = "Timepoint",
      values = colors_timepoint,
      guide = guide_legend(reverse = FALSE)
    ) +
    labs(
      x = NULL,
      y = "-log10(p-value) [UP | DOWN]",
      title = sprintf("Top %d Pathways: %s | %s\n(sorted by minimum p-value)", 
                      n_pathways, toupper(batch), comparison)
    ) +
    scale_y_continuous(
      limits = c(-y_limit, y_limit),
      breaks = seq(-y_limit, y_limit, by = 1)
    ) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
    theme_cowplot(font_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.y = element_text(size = 11, face = "bold"),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  print(plot)
  
  # Save
  if (save_plot) {
    filename_base <- sprintf("pathway_enrichment_%s_%s_top%d_FINAL", batch, comparison, n_pathways)
    
    plot_width <- max(12, n_pathways * 0.5)
    plot_height <- 8
    
    if ("pdf" %in% output_format) {
      pdf_file <- paste0(filename_base, ".pdf")
      ggsave(pdf_file, plot = plot, width = plot_width, height = plot_height, dpi = 300)
      cat(sprintf("✓ Saved: %s\n", pdf_file))
    }
    
    if ("png" %in% output_format) {
      png_file <- paste0(filename_base, ".png")
      ggsave(png_file, plot = plot, width = plot_width, height = plot_height, dpi = 300)
      cat(sprintf("✓ Saved: %s\n", png_file))
    }
  }
  
  cat(sprintf("\n✓ Pathways plotted: %d\n", n_pathways))
  cat(sprintf("✓ Total bars: %d\n", nrow(combined_data)))
  
  # RETURN pathway ranking
  return(list(
    plot = plot,
    pathway_ranking = top_pathways_master,
    batch = batch,
    comparison = comparison
  ))
}

## RUN STEP 13
result <- plot_pathway_enrichment(
  batch = "val",
  comparison = "dmso.vs.cxcr7",
  n_pathways = 30,
  save_plot = TRUE,
  output_format = c("pdf", "png")
)

# Save pathway ranking
top_pathways_master <- result$pathway_ranking
analysis_batch <- result$batch
analysis_comparison <- result$comparison

write.csv(top_pathways_master, 
          "pathway_ranking_master_FINAL.csv", 
          row.names = FALSE)

cat("\n", strrep("=", 80), "\n")
cat("✓ STEP 13 COMPLETE\n")
cat(sprintf("✓ Pathway list saved: top_pathways_master (%d pathways)\n", nrow(top_pathways_master)))
cat("✓ CSV saved: pathway_ranking_master_FINAL.csv\n")
cat(strrep("=", 80), "\n\n")

cat("⏭️  NOW RUN STEP 14 to create heatmaps using this EXACT pathway list\n\n")





###############################################################
## STEP 14: CREATE INDIVIDUAL PATHWAY HEATMAPS
## Using UNCOLLAPSED data, ranked by 1800s p-value
###############################################################

cat("\n", strrep("=", 80), "\n")
cat("STEP 14: CREATING INDIVIDUAL PATHWAY HEATMAPS\n")
cat(strrep("=", 80), "\n\n")

cat("✓ Using pathway ranking from STEP 13\n")
cat(sprintf("  Pathways: %d\n", nrow(top_pathways_master)))
cat(sprintf("  Batch: %s\n", analysis_batch))
cat(sprintf("  Comparison: %s\n\n", analysis_comparison))

# Display the pathway list we're using
cat("Pathway list from Step 13:\n")
cat(strrep("-", 80), "\n")
print(top_pathways_master %>% select(pathway_num, pathway), n = 30)
cat("\n")

## ------------------------------------------------------------
## A) Map genes to master pathways
## ------------------------------------------------------------

cat("Step A: Mapping genes to pathways...\n")

pathway_genes <- bind_rows(
  up_summary_long_filtered %>%
    filter(batch == analysis_batch, comparison == analysis_comparison) %>%
    select(pathway, substrates),
  down_summary_long_filtered %>%
    filter(batch == analysis_batch, comparison == analysis_comparison) %>%
    select(pathway, substrates)
) %>%
  distinct(pathway, substrates)

gene_pathway_map <- pathway_genes %>%
  separate_rows(substrates, sep = ";") %>%
  mutate(gene = trimws(substrates)) %>%
  inner_join(top_pathways_master %>% select(pathway, pathway_num), 
             by = "pathway") %>%
  select(pathway_num, pathway, gene) %>%
  distinct()

cat(sprintf("  ✓ Mapped %d genes to %d pathways\n\n", 
            n_distinct(gene_pathway_map$gene), 
            n_distinct(gene_pathway_map$pathway)))

## ------------------------------------------------------------
## B) Load UNCOLLAPSED phosphosite data
## ------------------------------------------------------------

cat("Step B: Loading UNCOLLAPSED phosphosite data...\n")

val_datasets <- c("val_10.dmso.vs.cxcr7", 
                  "val_600.dmso.vs.cxcr7", 
                  "val_1800.dmso.vs.cxcr7")

# Use RAW data (dfs_new_raw) instead of collapsed data
all_phosphosite_data <- lapply(val_datasets, function(time_name) {
  time_pt <- str_extract(time_name, "\\d+")
  dataset <- dfs_new_raw[[time_name]]  # ← UNCOLLAPSED!
  
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

phospho_wide_all <- all_phosphosite_data %>%
  pivot_wider(
    id_cols = c(phosphosite_id, name, PSite, uniprot_id),
    names_from = timepoint,
    values_from = c(logFC, PValue),
    names_sep = "_"
  ) %>%
  mutate(
    mean_abs_logFC = rowMeans(abs(cbind(logFC_10, logFC_600, logFC_1800)), na.rm = TRUE),
    min_pvalue = pmin(PValue_10, PValue_600, PValue_1800, na.rm = TRUE),
    pvalue_1800 = PValue_1800  # ← For ranking!
  )

cat(sprintf("  ✓ Loaded %d total phosphosites (UNCOLLAPSED)\n", nrow(phospho_wide_all)))

# Filter: At least 2/3 timepoints with data
phospho_wide_filtered <- phospho_wide_all %>%
  mutate(
    n_complete = rowSums(!is.na(cbind(logFC_10, logFC_600, logFC_1800)))
  ) %>%
  filter(n_complete >= 2)

cat(sprintf("  ✓ After filtering (≥2 timepoints): %d phosphosites\n\n", 
            nrow(phospho_wide_filtered)))

## ------------------------------------------------------------
## C) Function to create pathway heatmap (ranked by 1800s)
## ------------------------------------------------------------

create_pathway_heatmap_RANKED <- function(pathway_num, 
                                          pathway_name,
                                          gene_map,
                                          phospho_data,
                                          batch,
                                          comparison,
                                          max_sites = 20,
                                          output_dir = "pathway_heatmaps_FINAL") {
  
  cat(sprintf("  #%02d: %-65s | ", pathway_num, 
              substr(pathway_name, 1, 65)))
  
  pathway_genes <- gene_map %>%
    filter(pathway_num == !!pathway_num) %>%
    pull(gene) %>%
    unique()
  
  cat(sprintf("G:%2d | ", length(pathway_genes)))
  
  # Get ALL phosphosites for pathway genes
  pathway_psites <- phospho_data %>%
    filter(name %in% pathway_genes) %>%
    arrange(pvalue_1800)  # RANK by 1800s p-value
  
  cat(sprintf("P:%3d | ", nrow(pathway_psites)))
  
  if (nrow(pathway_psites) < 3) {
    cat("⚠️ SKIP (too few)\n")
    return(NULL)
  }
  
  # Take top N sites by 1800s p-value
  if (nrow(pathway_psites) > max_sites) {
    pathway_psites <- pathway_psites %>% head(max_sites)
    cat(sprintf("→%2d | ", max_sites))
  }
  
  # Create logFC matrix
  logfc_mat <- as.matrix(pathway_psites[, c("logFC_10", "logFC_600", "logFC_1800")])
  rownames(logfc_mat) <- pathway_psites$phosphosite_id
  colnames(logfc_mat) <- c("10 sec", "600 sec", "1800 sec")
  
  # Impute missing values with row mean
  logfc_mat_imputed <- t(apply(logfc_mat, 1, function(row) {
    if (any(is.na(row))) {
      row[is.na(row)] <- mean(row, na.rm = TRUE)
    }
    return(row)
  }))
  
  cat(sprintf("S:%2d | ", nrow(logfc_mat_imputed)))
  
  # Z-score normalization
  logfc_z <- t(scale(t(logfc_mat_imputed)))
  
  # Significance stars
  sig_mat <- matrix("", 
                    nrow = nrow(logfc_z), 
                    ncol = ncol(logfc_z),
                    dimnames = dimnames(logfc_z))
  
  for (i in 1:nrow(sig_mat)) {
    site_id <- rownames(logfc_z)[i]
    site_data <- pathway_psites[pathway_psites$phosphosite_id == site_id, ]
    
    if (!is.na(site_data$PValue_10) && site_data$PValue_10 < 0.05) {
      sig_mat[i, "10 sec"] <- "*"
    }
    if (!is.na(site_data$PValue_600) && site_data$PValue_600 < 0.05) {
      sig_mat[i, "600 sec"] <- "*"
    }
    if (!is.na(site_data$PValue_1800) && site_data$PValue_1800 < 0.05) {
      sig_mat[i, "1800 sec"] <- "*"
    }
  }
  
  n_sig <- sum(sig_mat == "*")
  cat(sprintf("★:%2d | ", n_sig))
  
  # Row labels with gene_psite
  row_labels <- paste0(pathway_psites$name, "_", pathway_psites$PSite)
  
  # Clean pathway name for filename
  clean_name <- gsub("[^A-Za-z0-9_]", "_", pathway_name)
  clean_name <- substr(clean_name, 1, 60)
  
  # Plot dimensions
  plot_height <- max(5, nrow(logfc_z) * 0.18 + 2)
  plot_width <- 9
  
  # Title
  title_text <- sprintf("Pathway #%d: %s\n%d sites | %d proteins | %s | %s",
                        pathway_num,
                        pathway_name,
                        nrow(logfc_z),
                        n_distinct(pathway_psites$name),
                        batch,
                        comparison)
  
  # Create PDF
  filename <- sprintf("%s/pathway_%02d_%s.pdf", output_dir, pathway_num, clean_name)
  pdf(filename, width = plot_width, height = plot_height)
  
  library(pheatmap)
  
  pheatmap(
    logfc_z,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    labels_row = row_labels,
    display_numbers = sig_mat,
    number_color = "black",
    # annotation_row = protein_annot,  # ← REMOVED!
    show_colnames = TRUE,
    fontsize_row = 7,
    fontsize_col = 11,
    fontsize_number = 9,
    cellwidth = 35,
    cellheight = 11,
    color = colorRampPalette(c("#0571B0", "#92C5DE", "#F7F7F7", "#F4A582", "#CA0020"))(100),
    main = title_text,
    border_color = "grey70"
  )
  
  dev.off()
  
  cat("✓\n")
  
  return(data.frame(
    pathway_num = pathway_num,
    pathway = pathway_name,
    n_sites = nrow(logfc_z),
    n_proteins = n_distinct(pathway_psites$name),
    n_significant = n_sig
  ))
}

## ------------------------------------------------------------
## D) Generate all heatmaps
## ------------------------------------------------------------

cat("Step C: Generating heatmaps (max 20 sites, ranked by 1800s)...\n")
cat(strrep("-", 80), "\n")

dir.create("pathway_heatmaps_FINAL", showWarnings = FALSE)

pathway_summaries <- list()

for (i in 1:nrow(top_pathways_master)) {
  pathway_info <- top_pathways_master[i, ]
  
  summary <- create_pathway_heatmap_RANKED(
    pathway_num = pathway_info$pathway_num,
    pathway_name = pathway_info$pathway,
    gene_map = gene_pathway_map,
    phospho_data = phospho_wide_filtered,
    batch = analysis_batch,
    comparison = analysis_comparison,
    max_sites = 20,  # ← LIMIT TO 20 SITES
    output_dir = "pathway_heatmaps_FINAL"
  )
  
  if (!is.null(summary)) {
    pathway_summaries[[i]] <- summary
  }
}

## ------------------------------------------------------------
## E) Summary report
## ------------------------------------------------------------

cat("\n", strrep("=", 80), "\n")
cat("SUMMARY REPORT\n")
cat(strrep("=", 80), "\n\n")

pathway_summaries_clean <- pathway_summaries[!sapply(pathway_summaries, is.null)]
summary_df <- bind_rows(pathway_summaries_clean)

cat("Statistics:\n")
cat(strrep("-", 80), "\n")
print(summary_df, row.names = FALSE)

write.csv(summary_df, "pathway_heatmaps_FINAL/pathway_summary.csv", row.names = FALSE)
cat("\n✓ Saved: pathway_heatmaps_FINAL/pathway_summary.csv\n")

cat("\n", strrep("-", 80), "\n")
cat("OVERALL STATISTICS\n")
cat(strrep("-", 80), "\n")
cat(sprintf("Step 13 pathways: %d\n", nrow(top_pathways_master)))
cat(sprintf("Step 14 heatmaps: %d\n", nrow(summary_df)))
cat(sprintf("Skipped: %d\n", nrow(top_pathways_master) - nrow(summary_df)))
cat(sprintf("Total sites: %d\n", sum(summary_df$n_sites)))
cat(sprintf("Total proteins: %d\n", sum(summary_df$n_proteins)))
cat(sprintf("Total significant: %d\n", sum(summary_df$n_significant)))

cat("\n", strrep("=", 80), "\n")
cat("✅ STEP 14 COMPLETE - UNCOLLAPSED DATA, RANKED BY 1800s!\n")
cat(strrep("=", 80), "\n\n")
cat("📁 Output: pathway_heatmaps_FINAL/\n")
cat(sprintf("   %d PDF files\n", nrow(summary_df)))
cat("   pathway_summary.csv\n\n")








###############################################################
## STEP 15: COMPREHENSIVE CANDIDATE CLUSTERING ANALYSIS
## Top 40 per timepoint (10s, 600s, 1800s) for temporal resolution
###############################################################

cat("\n", strrep("=", 80), "\n")
cat("STEP 16: CANDIDATE PHOSPHOSITE CLUSTERING\n")
cat(strrep("=", 80), "\n\n")

## ------------------------------------------------------------
## A) Prepare UNCOLLAPSED phosphosite matrix from top 30 pathways
## ------------------------------------------------------------

cat("Preparing UNCOLLAPSED phosphosite matrix...\n")

# Get all genes from the top 30 pathways
genes_from_top_pathways <- gene_pathway_map %>%
  pull(gene) %>%
  unique()

cat(sprintf("  Genes from top 30 pathways: %d\n", length(genes_from_top_pathways)))

# Extract UNCOLLAPSED phosphosite data
phospho_for_clustering <- phospho_wide_filtered %>%  # ← Already filtered (≥2 timepoints)
  filter(name %in% genes_from_top_pathways) %>%
  left_join(
    gene_pathway_map %>%
      group_by(gene) %>%
      summarise(
        pathway_nums = paste(sort(unique(pathway_num)), collapse = ","),
        n_pathways = n_distinct(pathway_num),
        .groups = "drop"
      ),
    by = c("name" = "gene")
  )

cat(sprintf("  Total phosphosites (uncollapsed): %d\n\n", nrow(phospho_for_clustering)))

## ------------------------------------------------------------
## B) Select top 40 per timepoint
## ------------------------------------------------------------

cat("Selecting top 40 phosphosites per timepoint...\n")

# Top 40 by 10s p-value
top40_10s <- phospho_for_clustering %>%
  arrange(PValue_10) %>%
  head(60) %>%
  mutate(selected_by = "10s")

cat(sprintf("  Top 40 at 10s:   p-value range %.2e - %.2e\n", 
            min(top40_10s$PValue_10, na.rm = TRUE),
            max(top40_10s$PValue_10, na.rm = TRUE)))

# Top 40 by 600s p-value
top40_600s <- phospho_for_clustering %>%
  arrange(PValue_600) %>%
  head(40) %>%
  mutate(selected_by = "600s")

cat(sprintf("  Top 40 at 600s:  p-value range %.2e - %.2e\n", 
            min(top40_600s$PValue_600, na.rm = TRUE),
            max(top40_600s$PValue_600, na.rm = TRUE)))

# Top 40 by 1800s p-value
top40_1800s <- phospho_for_clustering %>%
  arrange(PValue_1800) %>%
  head(40) %>%
  mutate(selected_by = "1800s")

cat(sprintf("  Top 40 at 1800s: p-value range %.2e - %.2e\n\n", 
            min(top40_1800s$PValue_1800, na.rm = TRUE),
            max(top40_1800s$PValue_1800, na.rm = TRUE)))

# Combine and deduplicate
# Combine and deduplicate
phospho_combined <- bind_rows(top40_10s, top40_600s, top40_1800s) %>%
  group_by(phosphosite_id) %>%
  slice(1) %>%  # Take first occurrence
  mutate(selected_by = paste(sort(unique(c(
    if_else(phosphosite_id %in% top40_10s$phosphosite_id, "10s", NA_character_),
    if_else(phosphosite_id %in% top40_600s$phosphosite_id, "600s", NA_character_),
    if_else(phosphosite_id %in% top40_1800s$phosphosite_id, "1800s", NA_character_)
  ))), collapse = ",")) %>%
  ungroup()

cat(sprintf("  Combined unique sites: %d\n", nrow(phospho_combined)))
cat("\n  Selection breakdown:\n")
selection_counts <- table(phospho_combined$selected_by)
for (sel in names(selection_counts)) {
  cat(sprintf("    %s: %d sites\n", sel, selection_counts[sel]))
}
cat("\n")

## ------------------------------------------------------------
## C) Create logFC matrix
## ------------------------------------------------------------

cat("Creating logFC matrix...\n")

# Create logFC matrix
logFC_mat_master <- as.matrix(phospho_combined[, c("logFC_10", "logFC_600", "logFC_1800")])
rownames(logFC_mat_master) <- paste0(phospho_combined$name, "_", phospho_combined$PSite)
colnames(logFC_mat_master) <- c("10s", "600s", "1800s")

# Impute missing values (for sites with 2/3 timepoints)
logFC_mat_master <- t(apply(logFC_mat_master, 1, function(row) {
  if (any(is.na(row))) {
    row[is.na(row)] <- mean(row, na.rm = TRUE)
  }
  return(row)
}))

# Z-score normalization (row-wise)
logFC_z_master <- t(scale(t(logFC_mat_master)))

cat(sprintf("  ✓ Z-scored matrix: %d sites × %d timepoints\n\n", 
            nrow(logFC_z_master), ncol(logFC_z_master)))

## ------------------------------------------------------------
## D) Create significance star matrix
## ------------------------------------------------------------

cat("Creating significance annotations...\n")

sig_stars_master <- matrix("", 
                           nrow = nrow(logFC_z_master), 
                           ncol = ncol(logFC_z_master),
                           dimnames = dimnames(logFC_z_master))

for (i in 1:nrow(sig_stars_master)) {
  site_data <- phospho_combined[i, ]
  
  if (!is.na(site_data$PValue_10) && site_data$PValue_10 < 0.05) {
    sig_stars_master[i, "10s"] <- "*"
  }
  if (!is.na(site_data$PValue_600) && site_data$PValue_600 < 0.05) {
    sig_stars_master[i, "600s"] <- "*"
  }
  if (!is.na(site_data$PValue_1800) && site_data$PValue_1800 < 0.05) {
    sig_stars_master[i, "1800s"] <- "*"
  }
}

n_sig <- sum(sig_stars_master == "*")
cat(sprintf("  ✓ Significant measurements: %d\n\n", n_sig))

## ------------------------------------------------------------
## E) Create timepoint annotation (FIXED - DYNAMIC)
## ------------------------------------------------------------

# Annotation showing which timepoint(s) selected each site
timepoint_annot <- data.frame(
  Selected = phospho_combined$selected_by,
  row.names = rownames(logFC_z_master)
)

# Get all unique selection combinations
unique_selections <- unique(phospho_combined$selected_by)

cat("Selection combinations found:\n")
print(table(phospho_combined$selected_by))
cat("\n")

# Define ALL possible colors
all_selection_colors <- c(
  "10s" = "#E41A1C",
  "600s" = "#377EB8", 
  "1800s" = "#4DAF4A",
  "10s,600s" = "#984EA3",
  "10s,1800s" = "#FF7F00",
  "600s,1800s" = "#A65628",
  "10s,600s,1800s" = "#F781BF"
)

# Only keep colors that actually exist in the data
annot_colors <- list(
  Selected = all_selection_colors[names(all_selection_colors) %in% unique_selections]
)

cat("Annotation colors defined for:\n")
print(names(annot_colors$Selected))
cat("\n")

## ------------------------------------------------------------
## F) Master heatmap - WITHOUT ANNOTATION (cleaner)
## ------------------------------------------------------------

cat("Creating master heatmap...\n")

library(pheatmap)

pdf("candidates_MASTER_heatmap_top40_per_timepoint.pdf", 
    width = 11, 
    height = max(14, nrow(logFC_z_master) * 0.12))

pheatmap(
  logFC_z_master,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  labels_row = rownames(logFC_z_master),
  display_numbers = sig_stars_master,
  number_color = "black",
  # annotation_row = timepoint_annot,  # ← REMOVED - causes errors
  # annotation_colors = annot_colors,   # ← REMOVED
  cellwidth = 35,
  cellheight = 8,
  fontsize_row = 6,
  fontsize_col = 11,
  fontsize_number = 7,
  color = colorRampPalette(c("#2166AC", "#67A9CF", "#D1E5F0", "#FFFFFF", 
                             "#FDDBC7", "#EF8A62", "#B2182B"))(100),
  main = sprintf("Top 40 Candidates per Timepoint (n=%d unique)\nZ-scored logFC | * = p<0.05",
                 nrow(logFC_z_master)),
  border_color = "grey70"
)

dev.off()
cat("  ✓ Saved: candidates_MASTER_heatmap_top40_per_timepoint.pdf\n\n")


## ------------------------------------------------------------
## G) Hierarchical clustering into 4 subclusters
## ------------------------------------------------------------

cat("Performing hierarchical clustering...\n")

dist_rows <- dist(logFC_z_master, method = "euclidean")
hc_rows <- hclust(dist_rows, method = "ward.D2")
subclusters <- cutree(hc_rows, k = 4)

cluster_counts <- table(subclusters)
cat("  Cluster sizes:\n")
for (i in 1:4) {
  cat(sprintf("    Cluster %d: %d sites\n", i, cluster_counts[i]))
}
cat("\n")

# Add cluster assignment to phospho_combined
phospho_combined$cluster <- subclusters

cat("✓ Cluster assignments added to data\n\n")


## ------------------------------------------------------------
## H) UPDATED Function to plot subclusters (NO ANNOTATION)
## ------------------------------------------------------------

plot_subcluster <- function(cluster_num, 
                            logFC_matrix, 
                            sig_matrix, 
                            phospho_data, 
                            subclusters,
                            top_n = NULL, 
                            sort_by = "pvalue",
                            save_pdf = FALSE, 
                            filename = NULL) {  # ← REMOVED annot parameters
  
  # Get rows in this cluster
  cluster_idx <- which(subclusters == cluster_num)
  sub_mat <- logFC_matrix[cluster_idx, , drop = FALSE]
  sub_data <- phospho_data[cluster_idx, ]
  
  # Sort and select top N
  if (!is.null(top_n) && nrow(sub_data) > top_n) {
    if (sort_by == "pvalue") {
      order_idx <- order(sub_data$min_pvalue)
      sub_data <- sub_data[order_idx[1:top_n], ]
      sub_mat <- sub_mat[order_idx[1:top_n], , drop = FALSE]
    } else if (sort_by == "logfc") {
      order_idx <- order(-sub_data$mean_abs_logFC)
      sub_data <- sub_data[order_idx[1:top_n], ]
      sub_mat <- sub_mat[order_idx[1:top_n], , drop = FALSE]
    }
  }
  
  # Create row labels with pathway numbers and selection info
  row_labels <- paste0(rownames(sub_mat), " [", sub_data$pathway_nums, "] (", sub_data$selected_by, ")")
  
  # Subset significance matrix
  sub_sig <- sig_matrix[rownames(sub_mat), , drop = FALSE]
  
  # Title
  title_text <- sprintf("Cluster %d (n=%d)", cluster_num, nrow(sub_mat))
  if (!is.null(top_n)) {
    title_text <- sprintf("%s - Top %d by %s", 
                          title_text, 
                          min(top_n, nrow(sub_mat)), 
                          ifelse(sort_by == "pvalue", "min p-value", "|logFC|"))
  }
  
  # Dimensions
  plot_height <- max(6, nrow(sub_mat) * 0.18 + 2)
  plot_width <- 12  # ← Wider for longer labels
  
  # Create plot
  if (save_pdf && !is.null(filename)) {
    pdf(filename, width = plot_width, height = plot_height)
  }
  
  pheatmap(
    sub_mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    labels_row = row_labels,
    display_numbers = sub_sig,
    number_color = "black",
    # annotation_row = sub_annot,  # ← REMOVED
    # annotation_colors = annot_cols,  # ← REMOVED
    fontsize_row = 6,  # ← Smaller font for longer labels
    fontsize_col = 11,
    fontsize_number = 8,
    cellwidth = 40,
    cellheight = 12,
    color = colorRampPalette(c("#2166AC", "#67A9CF", "#D1E5F0", "#FFFFFF", 
                               "#FDDBC7", "#EF8A62", "#B2182B"))(100),
    main = title_text,
    border_color = "grey70"
  )
  
  if (save_pdf && !is.null(filename)) {
    dev.off()
    cat(sprintf("  ✓ Saved: %s\n", basename(filename)))
  }
  
  invisible(sub_data)
}

## ------------------------------------------------------------
## I) Generate subcluster plots (UPDATED CALLS)
## ------------------------------------------------------------

cat("Generating subcluster heatmaps...\n")

# Create output directory
dir.create("candidate_subclusters", showWarnings = FALSE)
# For each cluster: full + top 20
for (i in 1:4) {
  cat(sprintf("\nCluster %d:\n", i))
  
  # Full cluster
  plot_subcluster(
    i, logFC_z_master, sig_stars_master, phospho_combined, subclusters,
    save_pdf = TRUE,
    filename = sprintf("candidate_subclusters/cluster_%d_full.pdf", i)
  )
  
  # Top 20 by min p-value
  top_sites <- plot_subcluster(
    i, logFC_z_master, sig_stars_master, phospho_combined, subclusters,
    top_n = 20, sort_by = "pvalue",
    save_pdf = TRUE,
    filename = sprintf("candidate_subclusters/cluster_%d_top20.pdf", i)
  )
  
  # Save top candidates to CSV (add cluster number manually)
  write.csv(top_sites %>% 
              mutate(cluster = i) %>%  # ← ADD cluster number
              select(name, PSite, phosphosite_id, 
                     logFC_10, logFC_600, logFC_1800,
                     PValue_10, PValue_600, PValue_1800,
                     min_pvalue, selected_by, pathway_nums, n_pathways, cluster), 
            sprintf("candidate_subclusters/cluster_%d_top20_candidates.csv", i),
            row.names = FALSE)

}

## ------------------------------------------------------------
## J) Characterize temporal patterns per cluster
## ------------------------------------------------------------

cat("\n\nCharacterizing temporal patterns...\n")

cluster_patterns <- phospho_combined %>%
  group_by(cluster) %>%
  summarise(
    n_sites = n(),
    n_early = sum(grepl("10s", selected_by)),
    n_mid = sum(grepl("600s", selected_by)),
    n_late = sum(grepl("1800s", selected_by)),
    mean_logFC_10 = mean(logFC_10, na.rm = TRUE),
    mean_logFC_600 = mean(logFC_600, na.rm = TRUE),
    mean_logFC_1800 = mean(logFC_1800, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    dominant_time = case_when(
      n_early > n_mid & n_early > n_late ~ "Early (10s)",
      n_late > n_mid & n_late > n_early ~ "Late (1800s)",
      n_mid > n_early & n_mid > n_late ~ "Mid (600s)",
      TRUE ~ "Mixed"
    )
  )

cat("\nCluster characterization:\n")
cat(strrep("-", 80), "\n")
print(cluster_patterns, row.names = FALSE)

write.csv(cluster_patterns, 
          "candidate_subclusters/cluster_characterization.csv",
          row.names = FALSE)

## ------------------------------------------------------------
## K) Final summary
## ------------------------------------------------------------

cat("\n", strrep("=", 80), "\n")
cat("CLUSTERING SUMMARY\n")
cat(strrep("=", 80), "\n\n")

cat(sprintf("Total unique phosphosites: %d\n", nrow(phospho_combined)))
cat(sprintf("  - Selected only at 10s: %d\n", sum(phospho_combined$selected_by == "10s")))
cat(sprintf("  - Selected only at 600s: %d\n", sum(phospho_combined$selected_by == "600s")))
cat(sprintf("  - Selected only at 1800s: %d\n", sum(phospho_combined$selected_by == "1800s")))
cat(sprintf("  - Selected at multiple times: %d\n", 
            sum(grepl(",", phospho_combined$selected_by))))
cat(sprintf("Number of clusters: 4\n"))

cat("\n", strrep("=", 80), "\n")
cat("✓ STEP 16 COMPLETE!\n")
cat(strrep("=", 80), "\n\n")

cat("📁 Output files:\n")
cat("   candidates_MASTER_heatmap_top40_per_timepoint.pdf\n")
cat("   candidate_subclusters/\n")
cat("     cluster_1-4_full.pdf\n")
cat("     cluster_1-4_top20.pdf\n")
cat("     cluster_1-4_top20_candidates.csv\n")
cat("     cluster_characterization.csv\n\n")

































###############################################################
## 15. COMPREHENSIVE CANDIDATE CLUSTERING ANALYSIS
## Master heatmap + subclusters with top candidates
###############################################################

cat("\n", strrep("=", 80), "\n")
cat("STEP 16: CANDIDATE PHOSPHOSITE CLUSTERING\n")
cat(strrep("=", 80), "\n\n")

## ------------------------------------------------------------
## A) Prepare phosphosite matrix from top 30 pathways
## ------------------------------------------------------------

cat("Preparing phosphosite matrix...\n")

# Get all genes from the top 30 pathways
genes_from_top_pathways <- gene_pathway_map %>%
  pull(gene) %>%
  unique()

cat(sprintf("  Genes from top 30 pathways: %d\n", length(genes_from_top_pathways)))

# Extract phosphosite data with pathway annotations
phospho_for_clustering <- phospho_wide_all %>%
  filter(name %in% genes_from_top_pathways) %>%
  left_join(
    gene_pathway_map %>%
      group_by(gene) %>%
      summarise(
        pathway_nums = paste(sort(unique(pathway_num)), collapse = ", "),
        n_pathways = n_distinct(pathway_num),
        .groups = "drop"
      ),
    by = c("name" = "gene")
  )

cat(sprintf("  Total phosphosites: %d\n", nrow(phospho_for_clustering)))

# Create logFC matrix (only complete cases)
logFC_mat_master <- as.matrix(phospho_for_clustering[, c("logFC_10", "logFC_600", "logFC_1800")])
rownames(logFC_mat_master) <- paste0(phospho_for_clustering$name, "_", phospho_for_clustering$PSite)
colnames(logFC_mat_master) <- c("10s", "600s", "1800s")

# Keep only complete cases
complete_idx <- complete.cases(logFC_mat_master)
logFC_mat_master <- logFC_mat_master[complete_idx, ]
phospho_complete <- phospho_for_clustering[complete_idx, ]

cat(sprintf("  Complete cases: %d phosphosites\n", nrow(logFC_mat_master)))

# Z-score normalization (row-wise)
logFC_z_master <- t(scale(t(logFC_mat_master)))

cat(sprintf("  ✓ Z-scored matrix: %d sites × %d timepoints\n\n", 
            nrow(logFC_z_master), ncol(logFC_z_master)))

## ------------------------------------------------------------
## B) Create significance star matrix
## ------------------------------------------------------------

cat("Creating significance annotations...\n")

sig_stars_master <- matrix("", 
                           nrow = nrow(logFC_z_master), 
                           ncol = ncol(logFC_z_master),
                           dimnames = dimnames(logFC_z_master))

for (i in 1:nrow(sig_stars_master)) {
  site_id <- rownames(logFC_z_master)[i]
  site_data <- phospho_complete[i, ]
  
  if (!is.na(site_data$PValue_10) && site_data$PValue_10 < 0.05) {
    sig_stars_master[i, "10s"] <- "*"
  }
  if (!is.na(site_data$PValue_600) && site_data$PValue_600 < 0.05) {
    sig_stars_master[i, "600s"] <- "*"
  }
  if (!is.na(site_data$PValue_1800) && site_data$PValue_1800 < 0.05) {
    sig_stars_master[i, "1800s"] <- "*"
  }
}

n_sig <- sum(sig_stars_master == "*")
cat(sprintf("  ✓ Significant measurements: %d\n\n", n_sig))

## ------------------------------------------------------------
## C) Master heatmap - all phosphosites
## ------------------------------------------------------------

cat("Creating master heatmap...\n")

library(pheatmap)

pdf("candidates_MASTER_heatmap_all_sites.pdf", 
    width = 10, 
    height = max(8, nrow(logFC_z_master) * 0.05))

pheatmap(
  logFC_z_master,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  display_numbers = sig_stars_master,
  number_color = "black",
  cellwidth = 25,
  cellheight = 3,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = sprintf("Master Heatmap: All Candidate Phosphosites (Z-scored logFC)\nn=%d sites | * = p<0.05",
                 nrow(logFC_z_master)),
  fontsize = 10,
  fontsize_number = 6,
  border_color = NA
)

dev.off()
cat("  ✓ Saved: candidates_MASTER_heatmap_all_sites.pdf\n\n")

## ------------------------------------------------------------
## D) Hierarchical clustering into 4 subclusters
## ------------------------------------------------------------

cat("Performing hierarchical clustering...\n")

dist_rows <- dist(logFC_z_master, method = "euclidean")
hc_rows <- hclust(dist_rows, method = "ward.D2")
subclusters <- cutree(hc_rows, k = 4)

cluster_counts <- table(subclusters)
cat("  Cluster sizes:\n")
for (i in 1:4) {
  cat(sprintf("    Cluster %d: %d sites\n", i, cluster_counts[i]))
}
cat("\n")

# Add cluster assignment to data
phospho_complete$cluster <- subclusters

## ------------------------------------------------------------
## E) Function to plot subclusters
## ------------------------------------------------------------

plot_subcluster <- function(cluster_num, 
                            logFC_matrix, 
                            sig_matrix, 
                            phospho_data, 
                            subclusters,
                            top_n = NULL, 
                            sort_by = "pvalue",
                            save_pdf = FALSE, 
                            filename = NULL) {
  
  # Get rows in this cluster
  cluster_rows <- names(subclusters[subclusters == cluster_num])
  sub_mat <- logFC_matrix[cluster_rows, , drop = FALSE]
  sub_data <- phospho_data[phospho_data$phosphosite_id %in% 
                             paste0(phospho_data$name, "_", phospho_data$PSite) &
                             subclusters == cluster_num, ]
  
  # Sort and select top N
  if (!is.null(top_n)) {
    if (sort_by == "pvalue") {
      sub_data <- sub_data %>%
        arrange(min_pvalue) %>%
        head(top_n)
    } else if (sort_by == "logfc") {
      sub_data <- sub_data %>%
        arrange(desc(mean_abs_logFC)) %>%
        head(top_n)
    }
    
    # Match rows
    row_ids <- paste0(sub_data$name, "_", sub_data$PSite)
    sub_mat <- sub_mat[row_ids, , drop = FALSE]
  }
  
  # Create row labels with pathway numbers
  row_labels <- paste0(rownames(sub_mat), " [", sub_data$pathway_nums, "]")
  
  # Subset significance matrix
  sub_sig <- sig_matrix[rownames(sub_mat), , drop = FALSE]
  
  # Title
  title_text <- sprintf("Cluster %d (n=%d)", cluster_num, nrow(sub_mat))
  if (!is.null(top_n)) {
    title_text <- sprintf("%s - Top %d by %s", 
                          title_text, top_n, 
                          ifelse(sort_by == "pvalue", "p-value", "|logFC|"))
  }
  
  # Dimensions
  plot_height <- max(6, nrow(sub_mat) * 0.18 + 2)
  plot_width <- 10
  
  # Create plot
  if (save_pdf && !is.null(filename)) {
    pdf(filename, width = plot_width, height = plot_height)
  }
  
  pheatmap(
    sub_mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    labels_row = row_labels,
    display_numbers = sub_sig,
    number_color = "black",
    fontsize_row = 7,
    fontsize_col = 11,
    fontsize_number = 8,
    cellwidth = 30,
    cellheight = 12,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    main = title_text,
    border_color = "grey70"
  )
  
  if (save_pdf && !is.null(filename)) {
    dev.off()
    cat(sprintf("  ✓ Saved: %s\n", basename(filename)))
  }
  
  invisible(sub_data)
}

## ------------------------------------------------------------
## F) Generate subcluster plots
## ------------------------------------------------------------

cat("Generating subcluster heatmaps...\n")

# Create output directory
dir.create("candidate_subclusters", showWarnings = FALSE)

# For each cluster: full + top 20 by p-value + top 20 by logFC
for (i in 1:4) {
  cat(sprintf("\nCluster %d:\n", i))
  
  # Full cluster
  plot_subcluster(
    i, logFC_z_master, sig_stars_master, phospho_complete, subclusters,
    save_pdf = TRUE,
    filename = sprintf("candidate_subclusters/cluster_%d_full.pdf", i)
  )
  
  # Top 20 by p-value
  top_pval <- plot_subcluster(
    i, logFC_z_master, sig_stars_master, phospho_complete, subclusters,
    top_n = 20, sort_by = "pvalue",
    save_pdf = TRUE,
    filename = sprintf("candidate_subclusters/cluster_%d_top20_pvalue.pdf", i)
  )
  
  # Top 20 by logFC
  top_logfc <- plot_subcluster(
    i, logFC_z_master, sig_stars_master, phospho_complete, subclusters,
    top_n = 20, sort_by = "logfc",
    save_pdf = TRUE,
    filename = sprintf("candidate_subclusters/cluster_%d_top20_logFC.pdf", i)
  )
  
  # Save top candidates to CSV
  write.csv(top_pval, 
            sprintf("candidate_subclusters/cluster_%d_top20_pvalue.csv", i),
            row.names = FALSE)
}

## ------------------------------------------------------------
## G) Characterize temporal patterns per cluster
## ------------------------------------------------------------

cat("\n\nCharacterizing temporal patterns...\n")

cluster_patterns <- phospho_complete %>%
  group_by(cluster) %>%
  summarise(
    n_sites = n(),
    n_significant = sum(min_pvalue < 0.05),
    mean_logFC_10 = mean(logFC_10, na.rm = TRUE),
    mean_logFC_600 = mean(logFC_600, na.rm = TRUE),
    mean_logFC_1800 = mean(logFC_1800, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pattern = case_when(
      abs(mean_logFC_1800) > abs(mean_logFC_10) ~ "Late response",
      abs(mean_logFC_10) > abs(mean_logFC_1800) ~ "Early response",
      TRUE ~ "Sustained"
    )
  )

cat("\nCluster characterization:\n")
cat(strrep("-", 80), "\n")
print(cluster_patterns, row.names = FALSE)

write.csv(cluster_patterns, 
          "candidate_subclusters/cluster_characterization.csv",
          row.names = FALSE)

## ------------------------------------------------------------
## H) Final summary
## ------------------------------------------------------------

cat("\n", strrep("=", 80), "\n")
cat("CLUSTERING SUMMARY\n")
cat(strrep("=", 80), "\n\n")

cat(sprintf("Total phosphosites analyzed: %d\n", nrow(logFC_z_master)))
cat(sprintf("Number of clusters: 4\n"))
cat(sprintf("Significant sites: %d (%.1f%%)\n", 
            sum(phospho_complete$min_pvalue < 0.05),
            100 * sum(phospho_complete$min_pvalue < 0.05) / nrow(phospho_complete)))

cat("\n", strrep("=", 80), "\n")
cat("✓ STEP 16 COMPLETE!\n")
cat(strrep("=", 80), "\n\n")

cat("📁 Output files:\n")
cat("   candidates_MASTER_heatmap_all_sites.pdf\n")
cat("   candidate_subclusters/\n")
cat("     cluster_1-4_full.pdf\n")
cat("     cluster_1-4_top20_pvalue.pdf\n")
cat("     cluster_1-4_top20_logFC.pdf\n")
cat("     cluster_1-4_top20_pvalue.csv\n")
cat("     cluster_characterization.csv\n\n")


















###############################################################
## 15. PHOSPHOSITE-LEVEL HEATMAP ANALYSIS - TOP PATHWAYS ONLY
## Extract phosphosites from TOP VISUALIZED pathways
###############################################################

cat("\n", strrep("=", 80), "\n")
cat("PHOSPHOSITE ANALYSIS: TOP PATHWAYS ONLY\n")
cat(strrep("=", 80), "\n\n")

## ------------------------------------------------------------
## A) Get TOP PATHWAYS from enrichment plot (e.g., top 30)
## ------------------------------------------------------------

cat("Step 1: Identifying top pathways from enrichment analysis...\n")

# Define how many top pathways to use (match your enrichment plot)
n_top_pathways <- 30

# Get pathway ranking (same as enrichment plot)
all_data_filtered <- bind_rows(
  up_summary_long_filtered %>% 
    filter(batch == "val", comparison == "dmso.vs.cxcr7") %>%
    mutate(regulation = "up") %>%
    group_by(pathway, timepoint_sec, regulation) %>%
    slice_min(pvalue, n = 1, with_ties = FALSE) %>%
    ungroup(),
  down_summary_long_filtered %>%
    filter(batch == "val", comparison == "dmso.vs.cxcr7") %>%
    mutate(regulation = "down") %>%
    group_by(pathway, timepoint_sec, regulation) %>%
    slice_min(pvalue, n = 1, with_ties = FALSE) %>%
    ungroup()
)

pathway_ranking <- all_data_filtered %>%
  filter(!is.na(pvalue), !is.infinite(pvalue)) %>%
  group_by(pathway) %>%
  summarize(
    min_pvalue = min(pvalue, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(min_pvalue) %>%
  mutate(pathway_num = row_number())

# Select TOP N pathways
top_pathways_list <- pathway_ranking %>%
  head(n_top_pathways) %>%
  pull(pathway)

# Create numbered lookup
pathway_numbers <- pathway_ranking %>%
  head(n_top_pathways) %>%
  select(pathway, pathway_num)

cat(sprintf("  ✓ Selected top %d pathways\n", n_top_pathways))
cat("\n  Top 10 pathways:\n")
print(pathway_ranking %>% head(10) %>% select(pathway_num, pathway, min_pvalue))

## ------------------------------------------------------------
## B) Extract genes ONLY from TOP pathways
## ------------------------------------------------------------

cat("\nStep 2: Extracting genes from top pathways only...\n")

# Get substrate lists ONLY from top pathways
pathway_genes_top <- bind_rows(
  up_summary_long_filtered %>%
    filter(batch == "val", comparison == "dmso.vs.cxcr7", 
           pathway %in% top_pathways_list) %>%
    select(pathway, substrates),
  down_summary_long_filtered %>%
    filter(batch == "val", comparison == "dmso.vs.cxcr7",
           pathway %in% top_pathways_list) %>%
    select(pathway, substrates)
) %>%
  distinct(pathway, substrates)

# Expand substrates to individual genes
gene_pathway_map_top <- pathway_genes_top %>%
  separate_rows(substrates, sep = ";") %>%
  mutate(gene = trimws(substrates)) %>%
  left_join(pathway_numbers, by = "pathway") %>%
  filter(!is.na(pathway_num)) %>%
  select(gene, pathway_num, pathway) %>%
  distinct()

# Create pathway membership with numbers
gene_pathway_summary_top <- gene_pathway_map_top %>%
  group_by(gene) %>%
  summarise(
    pathway_nums = paste(sort(unique(pathway_num)), collapse = ","),
    pathway_names = paste(unique(pathway), collapse = " | "),
    n_pathways = n_distinct(pathway_num),
    .groups = "drop"
  )

cat(sprintf("  ✓ Extracted %d unique genes from top %d pathways\n", 
            nrow(gene_pathway_summary_top), n_top_pathways))
cat("  Genes per pathway distribution:\n")
print(table(gene_pathway_summary_top$n_pathways))

## ------------------------------------------------------------
## C) Extract phosphosite data for these genes
## ------------------------------------------------------------

cat("\nStep 3: Extracting phosphosite data...\n")

val_datasets <- c("val_10.dmso.vs.cxcr7", 
                  "val_600.dmso.vs.cxcr7", 
                  "val_1800.dmso.vs.cxcr7")

phosphosite_data_top <- lapply(val_datasets, function(time_name) {
  
  time_pt <- str_extract(time_name, "(?<=val_)\\d+")
  dataset <- all_inputs_collapsed[[time_name]]
  
  cat(sprintf("  Processing %s...\n", time_name))
  
  result <- dataset %>%
    as.data.frame() %>%
    filter(!is.na(name), toupper(name) %in% toupper(gene_pathway_summary_top$gene)) %>%
    mutate(
      phosphosite_id = paste(name, PSite, sep = "_"),
      timepoint = time_pt
    ) %>%
    select(phosphosite_id, name, PSite, uniprot_id, logFC, PValue, timepoint)
  
  cat(sprintf("    Found %d phosphosites\n", nrow(result)))
  return(result)
}) %>% 
  bind_rows()

cat(sprintf("\n  ✓ Total: %d phosphosite measurements\n", nrow(phosphosite_data_top)))
cat(sprintf("    From %d proteins\n", n_distinct(phosphosite_data_top$name)))

# Convert to wide format
phospho_wide_top <- phosphosite_data_top %>%
  pivot_wider(
    id_cols = c(phosphosite_id, name, PSite, uniprot_id),
    names_from = timepoint,
    values_from = c(logFC, PValue),
    names_sep = "_"
  )

# Merge with pathway info
phospho_annotated_top <- phospho_wide_top %>%
  left_join(gene_pathway_summary_top, by = c("name" = "gene")) %>%
  mutate(
    mean_abs_logFC = rowMeans(abs(cbind(logFC_10, logFC_600, logFC_1800)), na.rm = TRUE),
    min_pvalue = pmin(PValue_10, PValue_600, PValue_1800, na.rm = TRUE)
  )

## ------------------------------------------------------------
## D) Select best phosphosites per protein
## ------------------------------------------------------------

cat("\nStep 4: Selecting top phosphosites per protein...\n")

# Select top 2 phosphosites per protein (fewer for cleaner visualization)
top_psites_top_pathways <- phospho_annotated_top %>%
  filter(!is.na(pathway_nums)) %>%
  group_by(name) %>%
  arrange(min_pvalue) %>%
  slice_head(n = 2) %>%
  ungroup()

cat(sprintf("  ✓ Selected %d phosphosites (top 2 per protein)\n", 
            nrow(top_psites_top_pathways)))
cat(sprintf("    From %d proteins\n", n_distinct(top_psites_top_pathways$name)))

## ------------------------------------------------------------
## E) Create matrix and filter complete cases
## ------------------------------------------------------------

cat("\nStep 5: Creating heatmap matrix...\n")

logFC_mat_top <- as.matrix(top_psites_top_pathways[, c("logFC_10", "logFC_600", "logFC_1800")])
rownames(logFC_mat_top) <- top_psites_top_pathways$phosphosite_id
colnames(logFC_mat_top) <- c("10 sec", "600 sec", "1800 sec")

# Keep only complete cases
complete_rows <- complete.cases(logFC_mat_top)
logFC_mat_top <- logFC_mat_top[complete_rows, ]
top_psites_complete <- top_psites_top_pathways[complete_rows, ]

cat(sprintf("  ✓ Matrix: %d sites × %d timepoints\n", 
            nrow(logFC_mat_top), ncol(logFC_mat_top)))

# Z-score normalization
logFC_z_top <- t(scale(t(logFC_mat_top)))

## ------------------------------------------------------------
## F) Create significance markers
## ------------------------------------------------------------

cat("\nStep 6: Creating significance markers...\n")

sig_mat_top <- matrix("", 
                      nrow = nrow(logFC_z_top), 
                      ncol = ncol(logFC_z_top),
                      dimnames = dimnames(logFC_z_top))

for (i in 1:nrow(sig_mat_top)) {
  site_id <- rownames(logFC_z_top)[i]
  site_data <- top_psites_complete[top_psites_complete$phosphosite_id == site_id, ]
  
  if (nrow(site_data) > 0) {
    if (!is.na(site_data$PValue_10) && site_data$PValue_10 < 0.05) {
      sig_mat_top[i, "10 sec"] <- "*"
    }
    if (!is.na(site_data$PValue_600) && site_data$PValue_600 < 0.05) {
      sig_mat_top[i, "600 sec"] <- "*"
    }
    if (!is.na(site_data$PValue_1800) && site_data$PValue_1800 < 0.05) {
      sig_mat_top[i, "1800 sec"] <- "*"
    }
  }
}

n_sig <- sum(sig_mat_top == "*")
cat(sprintf("  ✓ Found %d significant measurements (p < 0.05)\n", n_sig))

## ------------------------------------------------------------
## G) MASTER HEATMAP - All phosphosites
## ------------------------------------------------------------

cat("\nStep 7: Creating master heatmap...\n")

row_labels_master <- paste0(
  top_psites_complete$name, "_", top_psites_complete$PSite,
  " [", top_psites_complete$pathway_nums, "]"
)

library(pheatmap)

pdf("phosphosite_heatmap_top_pathways_MASTER.pdf", 
    width = 10, 
    height = max(10, nrow(logFC_z_top) * 0.15))

pheatmap(
  logFC_z_top,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  labels_row = row_labels_master,
  display_numbers = sig_mat_top,
  number_color = "black",
  fontsize_row = 6,
  fontsize_col = 12,
  fontsize_number = 8,
  cellwidth = 35,
  cellheight = 10,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = sprintf("Phosphosites from Top %d Pathways (Z-scored logFC)\n* = p<0.05 | [Pathway Numbers]", 
                 n_top_pathways),
  border_color = "grey80"
)

dev.off()
cat("  ✓ Saved: phosphosite_heatmap_top_pathways_MASTER.pdf\n")

## ------------------------------------------------------------
## H) HIERARCHICAL CLUSTERING into 4 subclusters
## ------------------------------------------------------------

cat("\nStep 8: Hierarchical clustering...\n")

dist_rows <- dist(logFC_z_top, method = "euclidean")
hc_rows <- hclust(dist_rows, method = "ward.D2")
subclusters <- cutree(hc_rows, k = 4)

cluster_counts <- table(subclusters)
cat("  Cluster sizes:\n")
for (i in 1:4) {
  cat(sprintf("    Cluster %d: %d sites\n", i, cluster_counts[i]))
}

## Continue in next message with subcluster plotting...




































###############################################################
## 8. PHOSPHOSITE-LEVEL HEATMAP ANALYSIS
## Extract top phosphosites and create clustered heatmaps
###############################################################

cat("\n", strrep("=", 80), "\n")
cat("PHOSPHOSITE-LEVEL HEATMAP ANALYSIS\n")
cat(strrep("=", 80), "\n\n")

## ------------------------------------------------------------
## A) Create pathway numbering from ranked pathways
## ------------------------------------------------------------

cat("Step 1: Creating pathway numbering based on significance ranking...\n")

# Get pathway ranking (reuse from enrichment plot)
all_data_filtered <- bind_rows(
  up_summary_long_filtered %>% 
    filter(batch == "val", comparison == "dmso.vs.cxcr7") %>%
    mutate(regulation = "up") %>%
    group_by(pathway, timepoint_sec, regulation) %>%
    slice_min(pvalue, n = 1, with_ties = FALSE) %>%
    ungroup(),
  down_summary_long_filtered %>%
    filter(batch == "val", comparison == "dmso.vs.cxcr7") %>%
    mutate(regulation = "down") %>%
    group_by(pathway, timepoint_sec, regulation) %>%
    slice_min(pvalue, n = 1, with_ties = FALSE) %>%
    ungroup()
)

pathway_ranking <- all_data_filtered %>%
  filter(!is.na(pvalue), !is.infinite(pvalue)) %>%
  group_by(pathway) %>%
  summarize(
    min_pvalue = min(pvalue, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(min_pvalue) %>%
  mutate(pathway_num = row_number())

# Create lookup table
pathway_numbers <- pathway_ranking %>%
  select(pathway, pathway_num)

cat("  ✓ Created numbering for", nrow(pathway_numbers), "pathways\n")
cat("  Top 5 pathways:\n")
print(pathway_ranking %>% head(5) %>% select(pathway_num, pathway, min_pvalue))

## ------------------------------------------------------------
## B) Extract genes from pathways
## ------------------------------------------------------------

cat("\nStep 2: Extracting genes from pathways...\n")

# Get substrate lists from enrichment results
pathway_genes <- bind_rows(
  up_summary_long_filtered %>%
    filter(batch == "val", comparison == "dmso.vs.cxcr7") %>%
    select(pathway, substrates),
  down_summary_long_filtered %>%
    filter(batch == "val", comparison == "dmso.vs.cxcr7") %>%
    select(pathway, substrates)
) %>%
  distinct(pathway, substrates)

# Expand substrates (semicolon-separated) to individual genes
gene_pathway_map <- pathway_genes %>%
  separate_rows(substrates, sep = ";") %>%
  mutate(gene = trimws(substrates)) %>%
  left_join(pathway_numbers, by = "pathway") %>%
  filter(!is.na(pathway_num)) %>%
  select(gene, pathway_num, pathway) %>%
  distinct()

# Create pathway membership as numbers (e.g., "1,5,12")
gene_pathway_summary <- gene_pathway_map %>%
  group_by(gene) %>%
  summarise(
    pathway_nums = paste(sort(unique(pathway_num)), collapse = ","),
    n_pathways = n_distinct(pathway_num),
    .groups = "drop"
  )

cat("  ✓ Extracted", nrow(gene_pathway_summary), "unique genes\n")
cat("  Genes per pathway count:\n")
print(table(gene_pathway_summary$n_pathways))



## ------------------------------------------------------------
## C) Build phosphosite-level matrix from validation data
## ------------------------------------------------------------

cat("\nStep 3: Building phosphosite matrix from validation data...\n")

# Check the actual column structure
cat("\n  Checking column structure...\n")
test_dataset <- all_inputs_collapsed[["val_10.dmso.vs.cxcr7"]]

cat("  Columns:", paste(colnames(test_dataset), collapse = ", "), "\n")
cat("  First 3 rows:\n")
print(head(test_dataset, 3))

# Extract val dmso.vs.cxcr7 data for all timepoints
val_datasets <- c("val_10.dmso.vs.cxcr7", 
                  "val_600.dmso.vs.cxcr7", 
                  "val_1800.dmso.vs.cxcr7")

phosphosite_data <- lapply(val_datasets, function(time_name) {
  
  time_pt <- str_extract(time_name, "(?<=val_)\\d+")
  
  # Get the dataset
  dataset <- all_inputs_collapsed[[time_name]]
  
  if (is.null(dataset)) {
    cat("  ⚠️  Warning: Dataset", time_name, "not found\n")
    return(NULL)
  }
  
  cat(sprintf("  Processing %s (n=%d rows)...\n", time_name, nrow(dataset)))
  
  # Filter for genes in pathways and create phosphosite_id
  result <- dataset %>%
    as.data.frame() %>%
    filter(!is.na(name), toupper(name) %in% toupper(gene_pathway_summary$gene)) %>%
    mutate(
      phosphosite_id = paste(name, PSite, sep = "_"),
      timepoint = time_pt
    ) %>%
    select(phosphosite_id, name, PSite, uniprot_id, logFC, PValue, timepoint)
  
  cat(sprintf("    Found %d matching phosphosites\n", nrow(result)))
  
  return(result)
}) %>% 
  bind_rows()

cat("\n  ✓ Found", nrow(phosphosite_data), "phosphosite measurements\n")
cat("    From", n_distinct(phosphosite_data$name), "proteins\n")

# Show sample
cat("\n  Sample data:\n")
print(head(phosphosite_data, 5))

# Check if we got data
if (nrow(phosphosite_data) == 0) {
  cat("\n⚠️  ERROR: No phosphosite data found!\n")
  stop("No phosphosite data found")
}

# Convert to wide format
phospho_wide <- phosphosite_data %>%
  pivot_wider(
    id_cols = c(phosphosite_id, name, PSite, uniprot_id),
    names_from = timepoint,
    values_from = c(logFC, PValue),
    names_sep = "_"
  )

# Merge with pathway numbering
phospho_annotated <- phospho_wide %>%
  left_join(gene_pathway_summary, by = c("name" = "gene"))

cat("  ✓ Created matrix with", nrow(phospho_annotated), "phosphosites\n")
cat("  ✓ With pathway annotations:", sum(!is.na(phospho_annotated$pathway_nums)), "\n")

## ------------------------------------------------------------
## D) Select top phosphosites per protein
## ------------------------------------------------------------

cat("\nStep 4: Selecting top phosphosites per protein...\n")

# Calculate average absolute logFC and minimum p-value across timepoints
phospho_annotated <- phospho_annotated %>%
  mutate(
    mean_abs_logFC = rowMeans(abs(cbind(logFC_10, logFC_600, logFC_1800)), na.rm = TRUE),
    min_pvalue = pmin(PValue_10, PValue_600, PValue_1800, na.rm = TRUE)
  )

# Select top 3 phosphosites per protein by significance
top_psites_per_protein <- phospho_annotated %>%
  filter(!is.na(pathway_nums)) %>%  # Only keep those with pathway annotation
  group_by(name) %>%
  arrange(min_pvalue) %>%
  slice_head(n = 3) %>%  # Top 3 per protein
  ungroup()

cat("  ✓ Selected", nrow(top_psites_per_protein), "top phosphosites\n")
cat("    From", n_distinct(top_psites_per_protein$name), "proteins\n")

## ------------------------------------------------------------
## E) Create logFC matrix for heatmap
## ------------------------------------------------------------

cat("\nStep 5: Creating heatmap matrix...\n")

# Extract logFC columns
logFC_mat <- as.matrix(top_psites_per_protein[, c("logFC_10", "logFC_600", "logFC_1800")])
rownames(logFC_mat) <- top_psites_per_protein$phosphosite_id
colnames(logFC_mat) <- c("10 sec", "600 sec", "1800 sec")

# Remove rows with any NAs
complete_rows <- complete.cases(logFC_mat)
logFC_mat <- logFC_mat[complete_rows, ]
top_psites_complete <- top_psites_per_protein[complete_rows, ]

cat("  ✓ Matrix dimensions:", nrow(logFC_mat), "sites ×", ncol(logFC_mat), "timepoints\n")

if (nrow(logFC_mat) == 0) {
  cat("  ⚠️  ERROR: No complete cases found!\n")
  stop("No phosphosites with complete data across all timepoints")
}

# Z-score normalization (row-wise)
logFC_z <- t(scale(t(logFC_mat)))

cat("  ✓ Applied Z-score normalization\n")

## ------------------------------------------------------------
## F) Create significance annotation matrix
## ------------------------------------------------------------

cat("\nStep 6: Creating significance annotations...\n")

sig_mat <- matrix("", 
                  nrow = nrow(logFC_z), 
                  ncol = ncol(logFC_z),
                  dimnames = dimnames(logFC_z))

for (i in 1:nrow(sig_mat)) {
  site_id <- rownames(logFC_z)[i]
  site_data <- top_psites_complete[top_psites_complete$phosphosite_id == site_id, ]
  
  if (nrow(site_data) > 0) {
    if (!is.na(site_data$PValue_10) && site_data$PValue_10 < 0.05) {
      sig_mat[i, "10 sec"] <- "*"
    }
    if (!is.na(site_data$PValue_600) && site_data$PValue_600 < 0.05) {
      sig_mat[i, "600 sec"] <- "*"
    }
    if (!is.na(site_data$PValue_1800) && site_data$PValue_1800 < 0.05) {
      sig_mat[i, "1800 sec"] <- "*"
    }
  }
}

n_sig <- sum(sig_mat == "*")
cat("  ✓ Created significance markers:", n_sig, "significant measurements\n")

## ------------------------------------------------------------
## G) MASTER OVERVIEW HEATMAP
## ------------------------------------------------------------

cat("\nStep 7: Creating master overview heatmap...\n")

# Create compact labels with pathway numbers
row_labels_overview <- paste0(
  top_psites_complete$name, "_", top_psites_complete$PSite,
  " [", top_psites_complete$pathway_nums, "]"
)

library(pheatmap)

pdf("phosphosite_heatmap_overview.pdf", width = 10, height = max(8, nrow(logFC_z) * 0.12))

pheatmap(
  logFC_z,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  labels_row = row_labels_overview,
  display_numbers = sig_mat,
  number_color = "black",
  fontsize_row = 6,
  fontsize_col = 12,
  fontsize_number = 8,
  cellwidth = 30,
  cellheight = 8,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Overview: Top Phosphosites by Protein (Z-scored logFC)\nNumbers = Pathway Ranks",
  border_color = "grey80"
)

dev.off()
cat("  ✓ Saved: phosphosite_heatmap_overview.pdf\n")

cat("\n", strrep("=", 80), "\n")
cat("✓ PHOSPHOSITE HEATMAP CREATED!\n")
cat(strrep("=", 80), "\n\n")





















































###############################################################
# 6. PHOSPHOSITE-LEVEL HEATMAP ANALYSIS - Extract relevant phosphoproteins and psites
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











