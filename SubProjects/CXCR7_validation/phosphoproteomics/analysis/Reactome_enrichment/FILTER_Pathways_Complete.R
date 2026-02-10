###############################################################
## FILTER REACTOME PATHWAYS TO SELECTED LIST
## Before running enrichment
## Uses your 59 curated pathways only
###############################################################

cat("\n", "="*70, "\n")
cat("FILTERING REACTOME PATHWAYS TO CURATED LIST\n")
cat("="*70, "\n\n")

# Your curated pathway list (59 pathways)
manual_path_refined <- c(
  
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
  
  # === GPCR & G-PROTEIN SIGNALING (8) ===
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
  
  # === PKA/cAMP SIGNALING (7) ===
  "PKA activation",
  "PKA activation in glucagon signalling",
  "Adenylate cyclase activating pathway",
  "Adenylate cyclase inhibitory pathway",
  "CREB1 phosphorylation through the activation of Adenylate Cyclase",
  "PKA-mediated phosphorylation of CREB",
  "PKA-mediated phosphorylation of key metabolic factors",
  
  # === PKC SIGNALING (1) ===
  "Gastrin-CREB signalling pathway via PKC and MAPK",
  
  # === cGMP/PKG SIGNALING (2) ===
  "cGMP effects",
  "Nitric oxide stimulates guanylate cyclase",
  
  # === AMPK SIGNALING (3) ===
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
  "Membrane Trafficking"
)

cat("Selected pathways:", length(manual_path_refined), "\n\n")

###############################################################
## STEP 1: Load FULL Reactome pathways (all 2,200+)
###############################################################

cat("Loading Reactome pathways...\n")

pathways <- as.list(reactomePATHID2EXTID)
path_names <- as.list(reactomePATHID2NAME)
name_id <- match(names(pathways), names(path_names))
names(pathways) <- unlist(path_names)[name_id]

# Filter to Homo sapiens only
pathways <- pathways[grepl("Homo sapiens", names(pathways), ignore.case = TRUE)]

# Convert gene IDs to symbols
pathways <- lapply(pathways, function(path) {
  gene_name <- unname(getSYMBOL(path, data = "org.Hs.eg"))
  toupper(unique(gene_name[!is.na(gene_name)]))
})

cat("✓ Loaded", length(pathways), "Homo sapiens pathways\n\n")

###############################################################
## STEP 2: CLEAN PATHWAY NAMES FOR MATCHING
###############################################################

cat("Cleaning pathway names for matching...\n\n")

# Function to clean pathway name
clean_pathway_name <- function(name) {
  name <- gsub("Homo sapiens: ", "", name, ignore.case = TRUE)
  name <- gsub("Homo sapiens >> ", "", name, ignore.case = TRUE)
  name <- trimws(name)
  return(name)
}

# Clean all pathway names
pathway_names_clean <- sapply(names(pathways), clean_pathway_name)
names(pathways) <- pathway_names_clean

###############################################################
## STEP 3: MATCH YOUR CURATED PATHWAYS
###############################################################

cat("Matching your curated pathways...\n\n")

# Function to find matching pathway
find_pathway <- function(query_name, all_names) {
  # Exact match
  exact <- grep(paste0("^", query_name, "$"), all_names, ignore.case = TRUE)
  if (length(exact) > 0) return(all_names[exact[1]])
  
  # Partial match
  partial <- grep(query_name, all_names, ignore.case = TRUE)
  if (length(partial) > 0) return(all_names[partial[1]])
  
  return(NA)
}

# Try to find each curated pathway
matched_pathways <- list()
not_found <- character()

for (path in manual_path_refined) {
  found_name <- find_pathway(path, pathway_names_clean)
  
  if (!is.na(found_name)) {
    matched_pathways[[path]] <- pathways[[found_name]]
    cat(sprintf("✓ %s\n", path))
  } else {
    not_found <- c(not_found, path)
    cat(sprintf("✗ %s\n", path))
  }
}

cat("\n", "="*70, "\n")
cat("MATCHING SUMMARY\n")
cat("="*70, "\n\n")

cat("✓ Successfully matched:", length(matched_pathways), "pathways\n")
cat("✗ Not found:", length(not_found), "pathways\n\n")

if (length(not_found) > 0) {
  cat("Not found (will be excluded):\n")
  for (p in not_found) {
    cat(sprintf("  - %s\n", p))
  }
}

###############################################################
## STEP 4: FINAL FILTERED PATHWAYS
###############################################################

cat("\n", "="*70, "\n")
cat("FINAL FILTERED PATHWAY SET\n")
cat("="*70, "\n\n")

selected_pathways_filtered <- matched_pathways

cat("Final filtered pathways:", length(selected_pathways_filtered), "\n")
cat("Selected from your list:", length(manual_path_refined), "\n")
cat("Match rate:", 
    round(100 * length(selected_pathways_filtered) / length(manual_path_refined), 1), 
    "%\n\n")

# Save info
selected_pathways_df <- data.frame(
  pathway_number = seq_along(selected_pathways_filtered),
  pathway_name = names(selected_pathways_filtered),
  n_genes = sapply(selected_pathways_filtered, length)
)

write.csv(selected_pathways_df, 
          "selected_pathways_for_enrichment.csv", 
          row.names = FALSE)

cat("✓ Saved selected pathways to: selected_pathways_for_enrichment.csv\n\n")

cat("="*70, "\n")
cat("✓ READY FOR ENRICHMENT!\n")
cat("Use: selected_pathways_filtered\n")
cat("="*70, "\n\n")

