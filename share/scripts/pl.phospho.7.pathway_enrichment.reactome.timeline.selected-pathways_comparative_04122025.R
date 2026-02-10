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
  "val_10.dmso.vs.cxcr7", "val_600.dmso.vs.cxcr7", "val_1800.dmso.vs.cxcr7",
  "val_10.dmso.vs.0s",    "val_600.dmso.vs.0s",    "val_1800.dmso.vs.0s",
  "val_10.cxcr7.vs.0s",   "val_600.cxcr7.vs.0s",   "val_1800.cxcr7.vs.0s"
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
  "init_10.cxcr7.vs.0s","init_30.cxcr7.vs.0s","init_60.cxcr7.vs.0s",
  "init_300.cxcr7.vs.0s","init_600.cxcr7.vs.0s","init_900.cxcr7.vs.0s",
  "init_1800.cxcr7.vs.0s"
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
## 1.2 Prepare Reactome pathways (only Homo sapiens)
#############################


pathways <- as.list(reactomePATHID2EXTID)
path_names <- as.list(reactomePATHID2NAME)
name_id <- match(names(pathways), names(path_names))
names(pathways) <- unlist(path_names)[name_id]
pathways <- pathways[grepl("Homo sapiens", names(pathways), ignore.case = TRUE)]
pathways <- lapply(pathways, function(path) {
  gene_name <- unname(getSYMBOL(path, data = "org.Hs.eg"))
  toupper(unique(gene_name))
})

#############################
## 2. Rank-Based Enrichment per Timepoint: Up vs. Down
#############################

# We'll process only the three main timepoints (first three elements)
up_results <- list()
down_results <- list()


# Loop over timepoints (i=1:3) for UP-regulated (alter = "greater")
#for (i in 1:3) {
for (i in seq_len(ncol(Tc.gene))) {
  res <- pathwayRankBasedEnrichment(Tc.gene[, i], annotation = pathways, alter = "greater")#[TcP.gene[,i] < 0.05]
  
  # Merge enrichment result with pathway sizes
  pathways2 <- as.data.frame(t(data.frame(lapply(pathways, length))))
  path3 <- as.data.frame(res)
  path4 <- merge(path3, pathways2, by = "row.names", all = FALSE)
  colnames(path4) <- c("pathway", "pvalue", "number.substrates", "substrates", "pw.size")
  
  # Clean and convert numeric columns
  path4 <- na.omit(path4)
  path4$number.substrates <- as.numeric(path4$number.substrates)
  path4$pw.size <- as.numeric(path4$pw.size)
  path4$pvalue <- as.numeric(path4$pvalue)
  path4$ratio <- path4$number.substrates / path4$pw.size
  
  # Clean pathway names
  path4$pathway <- gsub("_", " ", gsub("REACTOME_|Homo.sapiens..", "", path4$pathway))
  
  # Filter significant pathways (pvalue < 0.05) and sort by ratio (descending)
  #path4 <- path4[path4$pvalue < 0.05, ]
  path4 <- path4[order(path4$ratio, decreasing = TRUE), ]
  
  up_results[[ names_input[i] ]] <- path4
}

# Loop over timepoints (i=1:3) for DOWN-regulated (alter = "less")
#for (i in 1:3) {
for (i in seq_len(ncol(Tc.gene))) {
  res <- pathwayRankBasedEnrichment(Tc.gene[, i], annotation = pathways, alter = "less")
  
  pathways2 <- as.data.frame(t(data.frame(lapply(pathways, length))))
  path3 <- as.data.frame(res)
  path4 <- merge(path3, pathways2, by = "row.names", all = FALSE)
  colnames(path4) <- c("pathway", "pvalue", "number.substrates", "substrates", "pw.size")
  
  path4 <- na.omit(path4)
  path4$number.substrates <- as.numeric(path4$number.substrates)
  path4$pw.size <- as.numeric(path4$pw.size)
  path4$pvalue <- as.numeric(path4$pvalue)
  path4$ratio <- path4$number.substrates / path4$pw.size
  
  path4$pathway <- gsub("_", " ", gsub("REACTOME_|Homo.sapiens..", "", path4$pathway))
  
  #path4 <- path4[path4$pvalue < 0.05, ]
  path4 <- path4[order(path4$ratio, decreasing = TRUE), ]
  
  down_results[[ names_input[i] ]] <- path4
}



#############################
## 3. Save Enrichment Results
#############################

# Save UP-regulated enrichment tables
for (time in names(up_results)) {
  write.table(up_results[[time]],
              file = paste0("D:/Research/CXCR7_platelet_analysis/SubProjects/CXCR7_validation/phosphoproteomics/analysis/Reactome_enrichment/greater_", time, "_pathways.txt"),
              sep = "\t", row.names = FALSE)
}

# Save DOWN-regulated enrichment tables
for (time in names(down_results)) {
  write.table(down_results[[time]],
              file = paste0("D:/Research/CXCR7_platelet_analysis/SubProjects/CXCR7_validation/phosphoproteomics/analysis/Reactome_enrichment/less_", time, "_pathways.txt"),
              sep = "\t", row.names = FALSE)
}


#############################
## 4. Extract Top 5 Pathways (by ratio) per timepoint for each condition and Visualize
#############################
# For up-regulated:
top_up <- lapply(names(up_results)[1:16], function(time) {
  df <- up_results[[time]]
  if(nrow(df) > 0) {
    df_top <- df %>% arrange((pvalue)) %>% head(5)
    df_top$time <- time
    df_top$regulation <- "up"
    return(df_top)
  } else {
    return(NULL)
  }
})
top_up <- bind_rows(top_up)


# For down-regulated:
top_down <- lapply(names(down_results)[1:16], function(time) {
  df <- down_results[[time]]
  if(nrow(df) > 0) {
    df_top <- df %>% arrange((pvalue)) %>% head(5)
    df_top$time <- time
    df_top$regulation <- "down"
    # For plotting, convert ratio to negative values:
    df_top$ratio <- -df_top$ratio
    return(df_top)
  } else {
    return(NULL)
  }
})
top_down <- bind_rows(top_down)


# --- Combine top pathway data for each condition ---
# For plotting, we want a common grouping by pathway. If the same pathway is found at multiple timepoints, we'll show them together.


## First, extract the top pathways (from your previous step):
up_summary <- top_up %>% 
  dplyr::select(pathway, ratio, pvalue, time, regulation)
down_summary <- top_down %>% 
  dplyr::select(pathway, ratio, pvalue, time, regulation)



## Build a long-format summary for all up/down-regulated results from timepoints 10, 600, 1800.
up_summary_long_all <- lapply(names(up_results), function(t) {
  df <- up_results[[t]]
  if (!is.null(df) && nrow(df) > 0) {
    df %>% 
      dplyr::select(pathway, pvalue, ratio, substrates) %>%
      mutate(time = t)
  } else {
    NULL
  }
}) %>% bind_rows()


down_summary_long_all <- lapply(names(down_results), function(t) {
  df <- down_results[[t]]
  if (!is.null(df) && nrow(df) > 0) {
    df %>% 
      dplyr::select(pathway, pvalue, ratio, substrates) %>%
      mutate(time = t)
  } else {
    NULL
  }
}) %>% bind_rows()




## Filter long-format up/down data to include only slected pathways.
### Define your manual list of pathway keywords
library(dplyr)
library(stringr)

# Define your manual list of pathway keywords in natural language
manual_filter <- c(
  "Grb2 SOS provides linkage to MAPK signaling for Integrins",
  "Negative regulators of MAPK signaling",
  "Negative regulators of MAPK pathway",
  "Negative regulation of the PI3K AKT network",
  "Signalling by ERKs",
  "Integrin cell surface interactions",
  "RAF activation",
  "Regulation of RAS by GAPs",
  "Interleukin 37 signaling",
  "Effects of PIP2 hydrolysis",
  "Signal attenuation",
  "Sphingolipid de novo synthesis",
  "Sphingolipid metabolism",
  "Plasma Lipoprotein assembly remodeling and clearance",
  "AKT phosphorylates targets in cytosol",
  "DAP12 signaling",
  "cGMP signaling",
  "Activated Adenylyl cyclase synthesizes cyclic AMP",
  "cAMP degradation by Phosphodiesterases",
  "Adenylate cyclase produces cAMP",
  "RhoA GTPase cycle",
  "Mtor signalling",
  "Ca dependent events",
  "Clathrin mediated endocytosis",
  "ABC family proteins mediated transport",
  "Signaling by Hedgehog",
  "Autophagy",
  "Signaling by non receptor tyrosine kinases",
  "Extracellular Matrix Organization",
  "Degradation of extracellular matrix"
)

# Convert the manual filter to the dot-separated format used in your pathway names.
# Each space is replaced by a literal dot.
manual_filter_dots <- gsub(" ", "\\.", manual_filter)

# Also include pathways that mention "cAMP" or "cyclic AMP"
cAMP_filter <- "cAMP|cyclic AMP"

selected_pathways_cAMP <- unique(c(
  names(pathways)[str_detect(up_summary_long_all$pathway, regex("cGMP|cAMP|cyclic AMP|cyclic", ignore_case = TRUE)) ]
))

# Create the union of all pathway names in up_summary_long_all that match any of these keywords.
selected_pathways <- unique(up_summary_long_all$pathway[
  sapply(up_summary_long_all$pathway, function(x) {
    any(sapply(manual_filter_dots, function(keyword) {
      str_detect(x, regex(keyword, ignore_case = TRUE))
    })) || str_detect(x, regex(cAMP_filter, ignore_case = TRUE))
  })
])

# Now, filter the long-format data using the selected pathways:
up_summary_long <- up_summary_long_all %>% 
  filter(pathway %in% selected_pathways) %>% 
  mutate(pvalue = as.numeric(pvalue),
         ratio = as.numeric(ratio),
         neg_log10_p = -log10(pvalue))

down_summary_long <- down_summary_long_all %>% 
  filter(pathway %in% selected_pathways) %>% 
  mutate(pvalue = as.numeric(pvalue),
         ratio = as.numeric(ratio),
         neg_log10_p = -log10(pvalue))


## after fuzzy match of pathways, I have a selected set of pathways for further analysis 
manual_path_refined <- c(
  "RAF.activation",
  "Signal.attenuation",
  "DAP12.signaling",
  "Regulation.of.RAS.by.GAPs",
  "GRB2.SOS.provides.linkage.to.MAPK.signaling.for.Integrins.",
  "Integrin.cell.surface.interactions",
  "Clathrin.mediated.endocytosis",
  "Signalling.to.ERKs",
  "AKT.phosphorylates.targets.in.the.cytosol",
  "Negative.regulation.of.the.PI3K.AKT.network",
  "Negative.regulation.of.MAPK.pathway",
  "Signaling.by.Non.Receptor.Tyrosine.Kinases",
  "Interleukin.37.signaling",
  "Sphingolipid.metabolism",
  "Sphingolipid.de.novo.biosynthesis",
  "Plasma.lipoprotein.assembly..remodeling..and.clearance",
  "Effects.of.PIP2.hydrolysis",
  "cGMP.effects",
  "Ca.dependent.events",
  "MTOR.signalling",
  "RHOA.GTPase.cycle",
  "Autophagy",
  "ABC.family.proteins.mediated.transport",
  "Extracellular.matrix.organization",
  "Degradation.of.the.extracellular.matrix"
)

manual_path_refined <- c(
  "SHC1.events.in.ERBB2.signaling",
  "Signaling.by.ERBB2",
  "Constitutive.Signaling.by.Overexpressed.ERBB2",
  "Signaling.by.PDGFR.in.disease",
  "Signaling.by.FGFR1",
  "Signaling.by.MET",
  "Signaling.by.BRAF.and.RAF1.fusions",
  "Paradoxical.activation.of.RAF.signaling.by.kinase.inactive.BRAF",
  "Signaling.by.moderate.kinase.activity.BRAF.mutants",
  "GP1b.IX.V.activation.signalling",
  "RHOB.GTPase.cycle",
  "RHOC.GTPase.cycle",
  "RHOV.GTPase.cycle",
  "Cargo.recognition.for.clathrin.mediated.endocytosis",
  "Signaling.by.Hedgehog"
)




manual_path_refined <- c(
  "RAF.activation",
  "Signal.attenuation",
  "Negative.regulation.of.MAPK.pathway",
  "Signaling.by.BRAF.and.RAF1.fusions",
  "Negative.regulation.of.the.PI3K.AKT.network",
  "RHOA.GTPase.cycle",
  "DAP12.signaling",
  "Signalling.to.ERKs",
  "Degradation.of.the.extracellular.matrix",
  "Extracellular.matrix.organization",
  "MTOR.signalling",
  "Ca.dependent.events",
  "Clathrin.mediated.endocytosis",
  "Autophagy",
  "GRB2.SOS.provides.linkage.to.MAPK.signaling.for.Integrins.",
  "Regulation.of.RAS.by.GAPs",
  "Integrin.cell.surface.interactions",
  "Signaling.by.Non.Receptor.Tyrosine.Kinases",
  "AKT.phosphorylates.targets.in.the.cytosol",
  "Effects.of.PIP2.hydrolysis",
  "cGMP.effects",
  "Sphingolipid.metabolism",
  "SHC1.events.in.ERBB2.signaling",
  "Interleukin.37.signaling",
  "Sphingolipid.de.novo.biosynthesis"
)




# Remove duplicates:
selected_pathways <- unique(manual_path_refined )
print(selected_pathways)


up_summary_long <- up_summary_long_all %>% 
  filter(pathway %in% selected_pathways) %>% 
  mutate(pvalue = as.numeric(pvalue),
         ratio = as.numeric(ratio),
         neg_log10_p = -log10(pvalue))


down_summary_long <- down_summary_long_all %>% 
  filter(pathway %in% selected_pathways) %>% 
  mutate(pvalue = as.numeric(pvalue),
         ratio = as.numeric(ratio),
         neg_log10_p = -log10(pvalue))


# Transform to wide format:
up_summary_long <- up_summary_long %>% 
  mutate(ratio = as.numeric(ratio),
         neg_log10_p = as.numeric(neg_log10_p))

down_summary_long <- down_summary_long %>% 
  mutate(ratio = as.numeric(ratio),
         neg_log10_p = as.numeric(neg_log10_p))


up_wide <- dplyr::select(up_summary_long, time, pathway, neg_log10_p) %>% 
  pivot_wider(names_from = pathway, 
              values_from = neg_log10_p,
              values_fill = list(neg_log10_p = 0))

down_wide <- dplyr::select(down_summary_long, time, pathway, neg_log10_p) %>% 
  pivot_wider(names_from = pathway, 
              values_from = neg_log10_p,
              values_fill = list(neg_log10_p = 0))


common_paths <- selected_pathways[selected_pathways %in% colnames(up_wide)[-1]]
up_wide <- dplyr::select(up_wide, time, all_of(common_paths))


common_paths <- selected_pathways[selected_pathways %in% colnames(down_wide)[-1]]
down_wide <- dplyr::select(down_wide, time, all_of(common_paths))

# For the up-regulated data:
# Convert the wide data frame to long format.
common_paths <- selected_pathways[selected_pathways %in% colnames(up_wide)[-1]]
up_wide <- dplyr::select(up_wide, time, all_of(common_paths))

up_long <- up_wide %>% 
  pivot_longer(cols = -time, names_to = "pathway", values_to = "neg_log10_p")

# For the down-regulated data:
# Convert the wide data frame to long format.
common_paths <- selected_pathways[selected_pathways %in% colnames(down_wide)[-1]]
down_wide <- dplyr::select(down_wide, time, all_of(common_paths))

down_long <- down_wide %>% 
  pivot_longer(cols = -time, names_to = "pathway", values_to = "neg_log10_p")





#################################
## ggplot mirrored barplot (ALL CONDITIONS)
#################################

# UP and DOWN long tables already exist: up_long, down_long

# Mirror down-regulated values
up_long2 <- up_long %>% 
  mutate(value = neg_log10_p,
         regulation = "up")

down_long2 <- down_long %>% 
  mutate(value = -neg_log10_p,   # mirrored
         regulation = "down")

# Combine datasets
combined <- bind_rows(up_long2, down_long2)

# Order pathways according to your manual list
combined <- combined %>% 
  mutate(pathway = factor(pathway, levels = selected_pathways))

# Order contrasts using the exact order in your dataset
combined$time <- factor(combined$time, 
                        levels = names(up_results))   # <— IMPORTANT

# Colors (same as before)
violet_colors <- c("#E6E6FA", "#BA55D3", "#4B0082",
                   "#9370DB", "#8A2BE2", "#9400D3", 
                   "#9932CC", "#DA70D6", "#DDA0DD",
                   rep("grey70", 7))  # fallback if >3 contrasts

# Create plot
mirrored_plot <- ggplot(combined, aes(x = pathway, y = value, fill = time)) +
  geom_col(position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = violet_colors) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) +
  labs(x = "Pathway", y = "-log10(p-value)",
       title = "Mirrored Up- and Down-regulated Pathway Enrichment Across Conditions") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.title   = element_text(size = 10),
        plot.title   = element_text(size = 12, face = "bold"))

# Show
dev.new()
print(mirrored_plot)






#############################################################
## Filter for the 6 desired contrasts (val + init)
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
## Build mirrored combined frame
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
## Plot
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
# 10. PHOSPHOSITE-LEVEL HEATMAP ANALYSIS ----
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
# 11. SUBCLUSTER VISUALIZATION WITH PATHWAY NUMBERS ----
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











