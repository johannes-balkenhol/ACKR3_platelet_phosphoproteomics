###############################################################
## PHOSPHOPROTEOMICS PATHWAY ENRICHMENT PIPELINE
## NEW Validation Dataset (Default) + Optional OLD Initial Dataset
###############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(remotes)
  library(rlist)
  library(sjmisc)
  library(stringr)
  
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
cat("\nSTEP 4: Collapsing phosphosites by protein\n")
cat(strrep("─", 80), "\n\n")

# ← CHOOSE COLLAPSE METHOD
collapse_by <- "individual"  # Options: "mean_logfc", "max_logfc", "pvalue", "individual"

## Combine dmso.vs.cxcr7 timepoints to find best phosphosite per protein
combined_timepoints <- dplyr::bind_rows(
  dfs_new_raw$`val_10.dmso.vs.cxcr7` %>% dplyr::mutate(timepoint = "10"),
  dfs_new_raw$`val_600.dmso.vs.cxcr7` %>% dplyr::mutate(timepoint = "600"),
  dfs_new_raw$`val_1800.dmso.vs.cxcr7` %>% dplyr::mutate(timepoint = "1800")
)

## ────────────────────────────────────────────────────────────
## Conditional: Collapse globally or individually per timepoint
## ────────────────────────────────────────────────────────────

if (collapse_by == "individual") {
  
  cat("✓ Mode: INDIVIDUAL COLLAPSE per timepoint\n")
  cat("  (Each timepoint picks its own best phosphosite)\n\n")
  
  ## Collapse function for individual datasets
  collapse_individual <- function(df, method = "pvalue") {
    if (method == "pvalue") {
      df %>%
        dplyr::group_by(uniprot_id) %>%
        dplyr::slice_min(PValue, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup()
    } else if (method == "max_logfc") {
      df %>%
        dplyr::mutate(abs_logFC = abs(logFC)) %>%
        dplyr::group_by(uniprot_id) %>%
        dplyr::slice_max(abs_logFC, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::select(-abs_logFC)
    } else {  # mean_logfc - not applicable for single timepoint
      df %>%
        dplyr::mutate(abs_logFC = abs(logFC)) %>%
        dplyr::group_by(uniprot_id) %>%
        dplyr::slice_max(abs_logFC, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::select(-abs_logFC)
    }
  }
  
  # Collapse each dataset individually
  all_inputs_collapsed <- lapply(dfs_new_raw, function(df) {
    collapse_individual(df, method = "pvalue") %>%
      dplyr::arrange(uniprot_id) %>%
      dplyr::select(uniprot_id, name, PSite, Average, logFC, PValue)
  })
  
  # Track which phosphosites were selected at each timepoint
  selection_tracking <- dplyr::bind_rows(
    all_inputs_collapsed$`val_10.dmso.vs.cxcr7` %>% 
      dplyr::select(uniprot_id, name, PSite) %>% 
      dplyr::mutate(timepoint = "10s"),
    all_inputs_collapsed$`val_600.dmso.vs.cxcr7` %>% 
      dplyr::select(uniprot_id, name, PSite) %>% 
      dplyr::mutate(timepoint = "600s"),
    all_inputs_collapsed$`val_1800.dmso.vs.cxcr7` %>% 
      dplyr::select(uniprot_id, name, PSite) %>% 
      dplyr::mutate(timepoint = "1800s")
  )
  
  # Reshape to wide format to see differences
  selection_wide <- selection_tracking %>%
    tidyr::pivot_wider(
      id_cols = c(uniprot_id, name),
      names_from = timepoint,
      values_from = PSite
    )
  
  # Save tracking
  utils::write.csv(selection_wide, 
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
    kinase_selection <- selection_wide %>% dplyr::filter(name == kinase)
    
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
    dplyr::rowwise() %>%
    dplyr::mutate(n_unique = dplyr::n_distinct(c(`10s`, `600s`, `1800s`), na.rm = TRUE)) %>%
    dplyr::ungroup()
  
  cat("\n", strrep("-", 80), "\n")
  cat(sprintf("Proteins with SAME site across timepoints: %d\n", 
              sum(different_sites$n_unique == 1, na.rm = TRUE)))
  cat(sprintf("Proteins with DIFFERENT sites: %d\n", 
              sum(different_sites$n_unique > 1, na.rm = TRUE)))
  
} else {
  
  ## Global collapse (same phosphosite across all timepoints)
  best_psites <- combined_timepoints %>%
    dplyr::group_by(uniprot_id, name, PSite) %>%
    dplyr::summarise(
      mean_abs_logFC = mean(abs(logFC), na.rm = TRUE),
      max_abs_logFC = max(abs(logFC), na.rm = TRUE),
      min_pvalue = min(PValue, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(uniprot_id) %>%
    {
      if (collapse_by == "mean_logfc") {
        dplyr::slice_max(., mean_abs_logFC, n = 1, with_ties = FALSE)
      } else if (collapse_by == "max_logfc") {
        dplyr::slice_max(., max_abs_logFC, n = 1, with_ties = FALSE)
      } else {  # pvalue
        dplyr::slice_min(., min_pvalue, n = 1, with_ties = FALSE)
      }
    } %>%
    dplyr::ungroup() %>%
    dplyr::select(uniprot_id, PSite)
  
  cat(sprintf("✓ Selected %d phosphosites (one per protein)\n", nrow(best_psites)))
  
  all_inputs_collapsed <- lapply(dfs_new_raw, function(df) {
    df %>%
      dplyr::inner_join(best_psites, by = c("uniprot_id", "PSite")) %>%
      dplyr::arrange(uniprot_id) %>%
      dplyr::select(uniprot_id, name, PSite, Average, logFC, PValue)
  })
  
  method_name <- dplyr::case_when(
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
    dplyr::left_join(combined_timepoints %>% dplyr::distinct(uniprot_id, name), by = "uniprot_id")
  
  utils::write.csv(best_psites_annotated, 
                   "selected_phosphosites_for_enrichment.csv", 
                   row.names = FALSE)
  
  cat("\n✓ Saved: selected_phosphosites_for_enrichment.csv\n")
  cat("\n✓ Same phosphosite used across all timepoints!\n")
}

cat("\n")

###############################################################
## 5) BUILD LOG2FC MATRIX
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
## 6) LOAD REACTOME PATHWAYS & MATCH TO CURATED LIST
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
## 7) MATCH CURATED PATHWAYS TO REACTOME DATABASE
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
                                          max_sites = 999,
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
    max_sites = 999,  # ← LIMIT TO 20 SITES
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






## ------------------------------------------------------------
## F) Master heatmap - WITHOUT ANNOTATION (cleaner)
## Save as PNG - NO BORDERS
## ------------------------------------------------------------

cat("Creating master heatmap...\n")
library(pheatmap)

# Save as PNG
png("candidates_MASTER_heatmap_top40_per_timepoint.png", 
    width = 11 * 200, 
    height = max(14, nrow(logFC_z_master) * 0.12) * 1000,
    res = 1000)

pheatmap(
  logFC_z_master,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  display_numbers = sig_stars_master,
  number_color = "black",
  cellwidth = 10,
  cellheight = 7,
  fontsize_row = 6,
  fontsize_col = 10,
  fontsize_number = 7,
  color = colorRampPalette(c("#0571B0", "#92C5DE", "#F7F7F7", "#F4A582", "#CA0020"))(100),
  main = sprintf("Top 40 Candidates per Timepoint (n=%d unique)\nZ-scored logFC | * = p<0.05",
                 nrow(logFC_z_master)),
  border_color = NA
)

dev.off()

cat("  ✓ Saved: candidates_MASTER_heatmap_top40_per_timepoint.png\n\n")


## ------------------------------------------------------------
## G) Hierarchical clustering into 3 subclusters
## ------------------------------------------------------------

cat("Performing hierarchical clustering...\n")

dist_rows <- dist(logFC_z_master, method = "euclidean")
hc_rows <- hclust(dist_rows, method = "ward.D2")
subclusters <- cutree(hc_rows, k = 3)  # ← Already correct: 3 clusters

cluster_counts <- table(subclusters)
cat("  Cluster sizes:\n")
for (i in 1:3) {  # ← FIXED: Changed from 4 to 3
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
                            filename = NULL) {
  
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
  row_labels <- paste0(rownames(sub_mat), " | ", sub_data$pathway_nums)
  
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
  plot_width <- 12
  
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
    fontsize_row = 6,
    fontsize_col = 11,
    fontsize_number = 8,
    cellwidth = 10,
    cellheight = 10,
    color = colorRampPalette(c("#0571B0", "#92C5DE", "#F7F7F7", "#F4A582", "#CA0020"))(100),
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

# For each cluster: full + top 30
for (i in 1:3) {  # ← FIXED: Changed from 4 to 3
  cat(sprintf("\nCluster %d:\n", i))
  
  # Full cluster
  plot_subcluster(
    i, logFC_z_master, sig_stars_master, phospho_combined, subclusters,
    save_pdf = TRUE,
    filename = sprintf("candidate_subclusters/cluster_%d_full.pdf", i)
  )
  
  # Top 30 by min p-value
  top_sites <- plot_subcluster(
    i, logFC_z_master, sig_stars_master, phospho_combined, subclusters,
    top_n = 32, sort_by = "pvalue",
    save_pdf = TRUE,
    filename = sprintf("candidate_subclusters/cluster_%d_top30.pdf", i)  # ← Renamed to top30
  )
  
  # Save top candidates to CSV
  write.csv(top_sites %>% 
              mutate(cluster = i) %>%
              select(name, PSite, phosphosite_id, 
                     logFC_10, logFC_600, logFC_1800,
                     PValue_10, PValue_600, PValue_1800,
                     min_pvalue, selected_by, pathway_nums, n_pathways, cluster), 
            sprintf("candidate_subclusters/cluster_%d_top30_candidates.csv", i),  # ← Renamed to top30
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
cat(sprintf("Number of clusters: 3\n"))  # ← FIXED: Changed from 4 to 3

cat("\n", strrep("=", 80), "\n")
cat("✓ STEP 16 COMPLETE!\n")
cat(strrep("=", 80), "\n\n")

cat("📁 Output files:\n")
cat("   candidates_MASTER_heatmap_top40_per_timepoint.pdf\n")
cat("   candidate_subclusters/\n")
cat("     cluster_1-3_full.pdf\n")  # ← FIXED: Changed to 1-3
cat("     cluster_1-3_top30.pdf\n")  # ← FIXED: Changed to 1-3
cat("     cluster_1-3_top30_candidates.csv\n")  # ← FIXED: Changed to 1-3
cat("     cluster_characterization.csv\n\n")



