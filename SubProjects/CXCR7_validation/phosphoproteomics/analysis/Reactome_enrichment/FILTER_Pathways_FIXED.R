###############################################################
## FIXED FILTER PATHWAYS SCRIPT
## Handles case sensitivity and special characters properly
###############################################################

cat("\n", strrep("=", 70), "\n")
cat("FILTERING REACTOME PATHWAYS TO CURATED LIST (FIXED)\n")
cat(strrep("=", 70), "\n\n")

# Your 74 curated pathways
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
  "Synthesis of PIPs at the plasma membrane",
  "Phospholipid metabolism",
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
  "Role of phospholipids in phagocytosis",
  "Lysosphingolipid and LPA receptors",
  "Ceramide signalling",
  "Sphingolipid metabolism"
)

cat("Selected pathways:", length(manual_path_refined), "\n\n")

###############################################################
## IMPROVED MATCHING FUNCTION
## Handles case, special characters, and whitespace
###############################################################

# Improved function for better matching
find_pathway_improved <- function(query_name, all_names) {
  # Normalize both query and all names for comparison
  query_normalized <- tolower(trimws(query_name))
  
  # Try exact match (case-insensitive)
  for (i in seq_along(all_names)) {
    if (tolower(all_names[i]) == query_normalized) {
      return(all_names[i])
    }
  }
  
  # Try partial match with good scoring
  matches <- list()
  for (i in seq_along(all_names)) {
    name_normalized <- tolower(trimws(all_names[i]))
    
    # Calculate similarity (how many words match)
    query_words <- strsplit(query_normalized, " ")[[1]]
    name_words <- strsplit(name_normalized, " ")[[1]]
    
    common_words <- sum(query_words %in% name_words)
    similarity <- common_words / max(length(query_words), length(name_words))
    
    if (similarity > 0.6) {  # At least 60% word match
      matches[[length(matches) + 1]] <- list(
        name = all_names[i],
        similarity = similarity
      )
    }
  }
  
  # Return best match if found
  if (length(matches) > 0) {
    best_match <- matches[[which.max(sapply(matches, function(x) x$similarity))]]
    return(best_match$name)
  }
  
  return(NA)
}

###############################################################
## LOAD AND CLEAN PATHWAYS (same as before)
###############################################################

cat("Loading Reactome pathways...\n")

pathways <- as.list(reactomePATHID2EXTID)
path_names <- as.list(reactomePATHID2NAME)
name_id <- match(names(pathways), names(path_names))
names(pathways) <- unlist(path_names)[name_id]

pathways <- pathways[grepl("Homo sapiens", names(pathways), ignore.case = TRUE)]

pathways <- lapply(pathways, function(path) {
  gene_name <- unname(getSYMBOL(path, data = "org.Hs.eg"))
  toupper(unique(gene_name[!is.na(gene_name)]))
})

cat("✓ Loaded", length(pathways), "Homo sapiens pathways\n\n")

cat("Cleaning pathway names...\n\n")

# Clean pathway names
clean_pathway_name <- function(name) {
  name <- gsub("Homo sapiens: ", "", name, ignore.case = TRUE)
  name <- gsub("Homo sapiens >> ", "", name, ignore.case = TRUE)
  name <- trimws(name)
  return(name)
}

pathway_names_clean <- sapply(names(pathways), clean_pathway_name)
names(pathways) <- pathway_names_clean

###############################################################
## MATCH CURATED PATHWAYS WITH IMPROVED FUNCTION
###############################################################

cat("Matching your curated pathways with improved matching...\n\n")

matched_pathways <- list()
not_found <- character()

for (path in manual_path_refined) {
  found_name <- find_pathway_improved(path, pathway_names_clean)
  
  if (!is.na(found_name)) {
    matched_pathways[[path]] <- pathways[[found_name]]
    cat(sprintf("✓ %s\n", path))
  } else {
    not_found <- c(not_found, path)
    cat(sprintf("✗ %s\n", path))
  }
}

###############################################################
## SUMMARY
###############################################################

cat("\n", strrep("=", 70), "\n")
cat("MATCHING SUMMARY\n")
cat(strrep("=", 70), "\n\n")

cat("Successfully matched:", length(matched_pathways), "pathways\n")
cat("NOT found:", length(not_found), "pathways\n\n")

if (length(not_found) > 0) {
  cat("NOT FOUND (will be excluded):\n")
  for (p in not_found) {
    cat(sprintf("  - %s\n", p))
  }
  cat("\n")
}

###############################################################
## FINAL FILTERED PATHWAYS
###############################################################

selected_pathways_filtered <- matched_pathways

cat(strrep("=", 70), "\n")
cat("FILTERED PATHWAY SET READY FOR ENRICHMENT\n")
cat(strrep("=", 70), "\n\n")

cat("Final count:", length(selected_pathways_filtered), "pathways\n")
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
          "selected_pathways_for_enrichment_FIXED.csv", 
          row.names = FALSE)

cat("✓ Saved to: selected_pathways_for_enrichment_FIXED.csv\n\n")

cat(strrep("=", 70), "\n")
cat("✓ READY FOR ENRICHMENT!\n")
cat("Use: selected_pathways_filtered\n")
cat(strrep("=", 70), "\n\n")

