#############################
## Refined Phosphoproteomics Enrichment and Visualization Script
#############################

## --- INSTALL & LOAD PACKAGES ---
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("annotate", "org.Hs.eg.db", "reactome.db", "clusterProfiler", "ReactomePA", "limma"))

#Generat Reactome DB
#https://gitlab.com/42analytics1/public/reactome.db-r-package/-/blob/main/generate_reactome_package.R?ref_type=heads

install.packages(c("ggplot2", "cowplot", "pheatmap", "dplyr", "tidyr", "RColorBrewer", "sjmisc", "rlist", "directPA", "RMySQL"))

library(devtools)
devtools::install_github("PYangLab/PhosR")
devtools::install_github("JosephCrispell/basicPlotteR")


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


#############################
## 1. Prepare Collapsed Input & Pathway Annotation
#############################

## Define input2 (collapsed tables)
input <- list(top.collapse.10, top.collapse.600, top.collapse.1800, 
              top.collapse.10.dmso.vs.0s, top.collapse.600.dmso.vs.0s, top.collapse.1800.dmso.vs.0s,
              top.collapse.10.cxcr7.vs.0s, top.collapse.600.cxcr7.vs.0s, top.collapse.1800.cxcr7.vs.0s)

names_input <- c("10", "600", "1800", 
                 "10.dmso.vs.0s", "600.dmso.vs.0s", "1800.dmso.vs.0s",
                 "10.cxcr7.vs.0s", "600.cxcr7.vs.0s", "1800.cxcr7.vs.0s")

for (i in seq_along(input)) {
  input[[i]] <- input[[i]][order(rownames(input[[i]])), ]
}
names(input) <- names_input

## Build a logFC matrix from the collapsed tables.
Tc.gene <- as.matrix(cbind(input[[1]]$logFC, input[[2]]$logFC, input[[3]]$logFC,
                           input[[4]]$logFC, input[[5]]$logFC, input[[6]]$logFC,
                           input[[7]]$logFC, input[[8]]$logFC, input[[9]]$logFC))
rownames(Tc.gene) <- paste(sapply(strsplit(rownames(input[[1]]), ";"), "[[", 2))
colnames(Tc.gene) <- names_input

## Prepare Reactome pathways (only Homo sapiens)
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
for (i in 1:3) {
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
for (i in 1:3) {
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
              file = paste0("D:/Research/CXCR7_platelet_analysis/CXCR7_validation/phosphoproteomics/analysis/Reactome_enrichment/greater_", time, "_pathways.txt"),
              sep = "\t", row.names = FALSE)
}

# Save DOWN-regulated enrichment tables
for (time in names(down_results)) {
  write.table(down_results[[time]],
              file = paste0("D:/Research/CXCR7_platelet_analysis/CXCR7_validation/phosphoproteomics/analysis/Reactome_enrichment/less_", time, "_pathways.txt"),
              sep = "\t", row.names = FALSE)
}




#############################
## 4. Extract Top 5 Pathways (by ratio) per timepoint for each condition and Visualize
#############################

# For up-regulated:
top_up <- lapply(names(up_results)[1:3], function(time) {
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
top_down <- lapply(names(down_results)[1:3], function(time) {
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
up_summary_long_all <- lapply(c("10", "600", "1800"), function(t) {
  if (!is.null(up_results[[t]])) {
    df <- up_results[[t]] %>% 
      dplyr::select(pathway, pvalue, ratio, substrates) %>% 
      mutate(time = t)
    return(df)
  } else {
    return(NULL)
  }
}) %>% bind_rows()

down_summary_long_all <- lapply(c("10", "600", "1800"), function(t) {
  if (!is.null(down_results[[t]])) {
    df <- down_results[[t]] %>% 
      dplyr::select(pathway, pvalue, ratio, substrates) %>% 
      mutate(time = t)
    return(df)
  } else {
    return(NULL)
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
## ggplot barplot mirror
#################################

# Suppose up_long and down_long are already created from up_wide and down_wide respectively.
# up_long: contains columns "time", "pathway", and "neg_log10_p" for up-regulated pathways.
# down_long: contains columns "time", "pathway", and "neg_log10_p" for down-regulated pathways.
# Here we mirror the down-regulated values by multiplying by -1.

up_long2 <- up_long %>% 
  mutate(value = neg_log10_p, 
         regulation = "up")
down_long2 <- down_long %>% 
  mutate(value = -neg_log10_p,  # mirror the down values
         regulation = "down")

# Combine the two datasets.
combined <- bind_rows(up_long2, down_long2)

# (Optional) If you want the x-axis (pathways) to appear in a specific order, 
# you can convert pathway to a factor with your desired ordering.
# For example, if you have a selected_pathways vector (in the desired order):
combined <- combined %>% 
  mutate(pathway = factor(pathway, levels = selected_pathways))

# 1. Convert 'time' into a factor with the desired level order
combined$time <- factor(combined$time, levels = c("10", "600", "1800"))

# 2. Plot again
violet_colors <- c("#E6E6FA", "#BA55D3", "#4B0082")

mirrored_plot <- ggplot(combined, aes(x = pathway, y = value, fill = time)) +
  geom_col(position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = violet_colors) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) +
  labs(x = "Pathway", y = "-log10(pvalue)",
       title = "Mirrored Up- and Down-regulated Pathway Enrichment") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.title   = element_text(size = 10),
        plot.title   = element_text(size = 12, face = "bold"))

dev.new()
print(mirrored_plot)








################
### prepare heatmaps
################

### get the genes of the seleced pathways
filtered_df <- as.data.frame(path4) %>% 
  dplyr::select(pathway, substrates)

filtered_df2 <- filtered_df %>% 
  filter(pathway %in% manual_path_refined)

filtered_df2$pathway_number <- seq_len(nrow(filtered_df2))
View(filtered_df2)


### deaggregate the aggregated genes
gene_df <- filtered_df2 %>% 
  separate_rows(substrates, sep = ";")

gene_df <- gene_df %>% 
  mutate(substrates = trimws(substrates))

View(gene_df)


### group by genes and concatenate the pathways/pathways number
gene_df <- filtered_df2 %>%
  separate_rows(substrates, sep = ";") %>%
  mutate(
    substrates = trimws(substrates),
    substrates = as.character(substrates),
    pathway = as.character(pathway)
  )

# Now aggregate pathway_number
gene_summary <- aggregate(
  pathway_number ~ substrates,
  data = gene_df,
  FUN = function(x) paste(unique(x), collapse = "; ")
)

View(gene_summary)
print(dim(gene_summary))


gene_summary <- gene_summary %>%
  rename(pathway = pathway_number)

gene_summary <- gene_summary %>%
  rename(symbol = substrates)

### get top logFC values of each phosphosite
library(tibble)  # for column_to_rownames()

make_merged_df <- function(top_collapse_df, top_df) {
  step1 <- top_collapse_df %>%
    inner_join(
      top_df,
      by = c("logFC" = "logFC", "name" = "symbol")
      # (We’re not matching on PValue=adj.P.Val here; remove comment if needed)
    )
  
  step2 <- step1 %>%
    mutate(name_psite = paste(name, psite, sep = "_"))
  
  step3 <- step2 %>%
    select(
      Average,
      logFC,
      PValue       = adj.P.Val,   # rename adj.P.Val -> PValue
      uniprot_id,
      name         = name_psite,  # rename name_psite -> name
      symbol       = name         # rename 'name' (from top.collapse) -> symbol
    )
  
  step4 <- step3 %>%
    mutate(row_id = paste(uniprot_id, symbol, sep = ";"))
  
  final_df <- step4 %>%
    column_to_rownames("row_id")
  
  final_df
}

# --- Now call the helper function for each pair --- #

top.collapse.10.2 <- make_merged_df(top.collapse.10, top.10)
top.collapse.600.2 <- make_merged_df(top.collapse.600, top.600)
top.collapse.1800.2 <- make_merged_df(top.collapse.1800, top.1800)



### merge the new annotated top tables collpased by the top site
merged_10_600 <- merge(top.collapse.10.2, top.collapse.600.2, by = 0, all = TRUE)
rownames(merged_10_600) <- merged_10_600$Row.names
merged_10_600$Row.names <- NULL

merged_all <- merge(merged_10_600, top.collapse.1800.2, by = 0, all = TRUE)
rownames(merged_all) <- merged_all$Row.names
merged_all$Row.names <- NULL

colnames(merged_all)[colnames(merged_all) == "logFC"]   <- "logFC.1800"
colnames(merged_all)[colnames(merged_all) == "PValue"] <- "PValue.1800"
colnames(merged_all)[colnames(merged_all) == "logFC.x"]   <- "logFC.10"
colnames(merged_all)[colnames(merged_all) == "PValue.x"] <- "PValue.10"
colnames(merged_all)[colnames(merged_all) == "logFC.y"]   <- "logFC.600"
colnames(merged_all)[colnames(merged_all) == "PValue.y"] <- "PValue.600"

filtered_cols <- merged_all[, c("uniprot_id",   
                                "name",
                                "symbol",    
                                "logFC.10", 
                                "PValue.10",
                                "logFC.600", 
                                "PValue.600",
                                "logFC.1800",
                                "PValue.1800")]



merged_final <- merge(filtered_cols, gene_summary, by = "symbol")




##### make a heatmap of genes in the pathway
# If needed, install pheatmap:
# install.packages("pheatmap")

library(pheatmap)

# -- 1) Subset your merged_final data to the three logFC columns --
# We assume merged_final has columns: "logFC.10", "logFC.600", "logFC.1800"
# Also "PValue.10", "PValue.600", "PValue.1800" for significance
# and "name" for row identification

logFC_mat <- as.matrix(merged_final[, c("logFC.10", "logFC.600", "logFC.1800")])

# Set row names so each row (gene/protein) is identifiable internally
# For example, use the 'name' column
rownames(logFC_mat) <- merged_final$name

# -- 2) (Optional) Z-score each row --
# Scales each row to mean=0, sd=1
logFC_z <- t(scale(t(logFC_mat)))

# -- 3) Build a matrix of significance symbols (stars) --
# Initialize a matrix with the same shape, filled with ""
sig_stars <- matrix("",
                    nrow = nrow(logFC_z), 
                    ncol = ncol(logFC_z),
                    dimnames = list(rownames(logFC_z), colnames(logFC_z)))

# If p-value < 0.05, put a star in the corresponding cell
for (i in seq_len(nrow(sig_stars))) {
  if (merged_final$PValue.10[i] < 0.05) {
    sig_stars[i, "logFC.10"] <- "*"
  }
  if (merged_final$PValue.600[i] < 0.05) {
    sig_stars[i, "logFC.600"] <- "*"
  }
  if (merged_final$PValue.1800[i] < 0.05) {
    sig_stars[i, "logFC.1800"] <- "*"
  }
}

dev.new()
#logFC_z_t       <- t(logFC_z)
#sig_stars_t     <- t(sig_stars)


# -- 4) Create the heatmap with pheatmap --
pheatmap(
  logFC_z,                # The row-scaled logFC matrix
  cluster_rows  = TRUE,   # Hierarchical clustering on rows
  cluster_cols  = FALSE,   # Hierarchical clustering on columns
  show_rownames = FALSE,  # Don't show row labels (we have 230 rows)
  display_numbers = sig_stars,   # Show the star if p < 0.05
  number_color   = "black",      # Color of the text overlay
  main           = "Heatmap of logFC (Z-scored) with Significance Stars"
)




##### make heatmap of subclusters
# Suppose you have already run something like:

library(pheatmap)

logFC_mat <- as.matrix(merged_final[, c("logFC.10", "logFC.600", "logFC.1800")])
rownames(logFC_mat) <- merged_final$name

# Z-score by row
logFC_z <- t(scale(t(logFC_mat)))


# 1. Compute distance among rows
dist_rows <- dist(logFC_z)

# 2. Hierarchical clustering
hc_rows <- hclust(dist_rows, method = "complete")  # or "ward.D2", "average", etc.

# 3. Cut tree into 4 clusters
subclusters <- cutree(hc_rows, k = 4)

# subclusters is an integer vector (1..4) that tells you which cluster each row belongs to.
# Make sure it has names identical to rownames(logFC_z).




## take all from a cluster
for (i in 1:4) {
  # Identify rows in subcluster i
  cluster_i_rows <- names(subclusters[subclusters == i])
  sub_mat <- logFC_z[cluster_i_rows, , drop = FALSE]
  
  sub_df <- merged_final[match(cluster_i_rows, merged_final$name), ]
  
  # Create row labels
  new_labels <- paste(sub_df$name, sub_df$pathway, sep = " | ")
  
  # Build significance star matrix for just these rows
  sig_stars <- matrix("",
                      nrow = nrow(sub_mat), 
                      ncol = ncol(sub_mat),
                      dimnames = list(rownames(sub_mat), colnames(sub_mat)))
  
  # For each row in the subcluster, if PValue.10 < 0.05, place a star in [row, "logFC.10"], etc.
  for (r in rownames(sub_mat)) {
    idx <- which(merged_final$name == r)
    if (merged_final$PValue.10[idx] < 0.05) {
      sig_stars[r, "logFC.10"] <- "*"
    }
    if (merged_final$PValue.600[idx] < 0.05) {
      sig_stars[r, "logFC.600"] <- "*"
    }
    if (merged_final$PValue.1800[idx] < 0.05) {
      sig_stars[r, "logFC.1800"] <- "*"
    }
  }
  
  # Open a new graphics window
  dev.new()
  
  pheatmap(
    sub_mat,
    cluster_rows    = TRUE,
    cluster_cols    = FALSE,
    display_numbers = sig_stars,
    number_color    = "black",
    labels_row      = new_labels,
    main            = paste("Subcluster", i, "Heatmap (with stars)"),
    scale           = "none"
  )
}




## subset pvalue
library(pheatmap)

for (i in 4:4) {
  # Identify rows in subcluster i
  cluster_i_rows <- names(subclusters[subclusters == i])
  
  # Subset the Z-scored matrix to only these rows
  sub_mat <- logFC_z[cluster_i_rows, , drop = FALSE]
  
  # Subset your main data frame to these rows
  sub_df <- merged_final[match(cluster_i_rows, merged_final$name), ]
  
  # Create row labels
  new_labels <- paste(sub_df$name, sub_df$pathway, sep = " | ")
  
  # Calculate the minimum p-value across the three time points
  # (or pick whichever method you prefer)
  sub_df$minPValue <- pmin(sub_df$PValue.10, sub_df$PValue.600, sub_df$PValue.1800,
                           na.rm = TRUE)
  
  # Sort by minPValue ascending (smallest = most significant)
  sub_df <- sub_df[order(sub_df$minPValue), ]
  
  # Keep only the top 20 rows
  top_n <- min(35, nrow(sub_df))
  sub_df_top20 <- sub_df[seq_len(top_n), ]
  
  # Now figure out which rows those correspond to
  top20_rows <- sub_df_top20$name
  
  # Subset the expression matrix accordingly, in the same order
  # Use match(...) to preserve the order in sub_df_top20
  sub_mat_top20 <- sub_mat[match(top20_rows, rownames(sub_mat)), , drop = FALSE]
  
  # Build significance star matrix for these top-20 rows
  sig_stars <- matrix("",
                      nrow = nrow(sub_mat_top20), 
                      ncol = ncol(sub_mat_top20),
                      dimnames = list(rownames(sub_mat_top20), colnames(sub_mat_top20)))
  
  for (r in rownames(sub_mat_top20)) {
    idx <- which(merged_final$name == r)
    if (merged_final$PValue.10[idx] < 0.05) {
      sig_stars[r, "logFC.10"] <- "*"
    }
    if (merged_final$PValue.600[idx] < 0.05) {
      sig_stars[r, "logFC.600"] <- "*"
    }
    if (merged_final$PValue.1800[idx] < 0.05) {
      sig_stars[r, "logFC.1800"] <- "*"
    }
  }
  
  # Subset the labels to the top-20
  # 'new_labels' must match the same row order
  new_labels_top20 <- paste(
    sub_df_top20$name,
    sub_df_top20$pathway,
    sep = " | "
  )
  
  # Optional: open a new dev window (on some systems dev.new() may not work)
  # dev.new()
  
  pheatmap(
    sub_mat_top20,
    cluster_rows    = TRUE,
    cluster_cols    = FALSE,
    display_numbers = sig_stars,
    number_color    = "black",
    labels_row      = new_labels_top20,
    main            = paste("Subcluster", i, "Heatmap (Top 20 by min P-value)"),
    scale           = "none"
  )
}




## subset logFC
library(pheatmap)

for (i in 1:1) {
  # Identify rows in subcluster i
  cluster_i_rows <- names(subclusters[subclusters == i])
  
  # Subset the Z-scored matrix for these rows
  sub_mat <- logFC_z[cluster_i_rows, , drop = FALSE]
  
  # Subset your main data frame
  sub_df <- merged_final[match(cluster_i_rows, merged_final$name), ]
  
  # Create row labels
  new_labels <- paste(sub_df$name, sub_df$pathway, sep = " | ")
  
  # 1) Compute maximum absolute logFC across the three timepoints
  sub_df$maxAbsFC <- pmax(
    abs(sub_df$logFC.10),
    abs(sub_df$logFC.600),
    abs(sub_df$logFC.1800),
    na.rm = TRUE
  )
  
  # 2) Sort descending by 'maxAbsFC'
  sub_df <- sub_df[order(sub_df$maxAbsFC, decreasing = TRUE), ]
  
  # 3) Keep the top 20 rows (or fewer if cluster has <20)
  top_n <- min(20, nrow(sub_df))
  sub_df_top20 <- sub_df[seq_len(top_n), ]
  
  # Figure out which row names those correspond to:
  top20_rows <- sub_df_top20$name
  
  # Subset the expression matrix accordingly, in the new order
  sub_mat_top20 <- sub_mat[match(top20_rows, rownames(sub_mat)), , drop = FALSE]
  
  # Build significance star matrix for these top-20 rows
  sig_stars <- matrix("",
    nrow = nrow(sub_mat_top20), 
    ncol = ncol(sub_mat_top20),
    dimnames = list(rownames(sub_mat_top20), colnames(sub_mat_top20))
  )
  
  for (r in rownames(sub_mat_top20)) {
    idx <- which(merged_final$name == r)
    if (merged_final$PValue.10[idx] < 0.05) {
      sig_stars[r, "logFC.10"] <- "*"
    }
    if (merged_final$PValue.600[idx] < 0.05) {
      sig_stars[r, "logFC.600"] <- "*"
    }
    if (merged_final$PValue.1800[idx] < 0.05) {
      sig_stars[r, "logFC.1800"] <- "*"
    }
  }
  
  # Subset the labels in the same order
  new_labels_top20 <- paste(
    sub_df_top20$name,
    sub_df_top20$pathway,
    sep = " | "
  )
  
  # Plot the heatmap with only the top-20 rows by maxAbsFC
  pheatmap(
    sub_mat_top20,
    cluster_rows    = TRUE,
    cluster_cols    = FALSE,
    display_numbers = sig_stars,
    number_color    = "black",
    labels_row      = new_labels_top20,
    main            = paste("Subcluster", i, "Heatmap (Top 20 by maxAbsFC)"),
    scale           = "none"
  )
}






### save subcluster to file


head(subclusters)

outdir <- "D:/Research/SubProjects/CXCR7_platelet_analysis/CXCR7_validation/phosphoproteomics/analysis/Reactome_enrichment"

for (i in 1:4) {
  # Identify rows in subcluster i
  cluster_i_rows <- names(subclusters[subclusters == i])
  sub_mat <- logFC_z[cluster_i_rows, , drop = FALSE]
  
  sub_df <- merged_final[match(cluster_i_rows, merged_final$name), ]
  
  # Create row labels
  new_labels <- paste(sub_df$name, sub_df$pathway, sep = " | ")
  
  # Build significance star matrix for just these rows
  sig_stars <- matrix("",
                      nrow = nrow(sub_mat),
                      ncol = ncol(sub_mat),
                      dimnames = list(rownames(sub_mat), colnames(sub_mat)))
  
  for (r in rownames(sub_mat)) {
    idx <- which(merged_final$name == r)
    if (merged_final$PValue.10[idx] < 0.05) {
      sig_stars[r, "logFC.10"] <- "*"
    }
    if (merged_final$PValue.600[idx] < 0.05) {
      sig_stars[r, "logFC.600"] <- "*"
    }
    if (merged_final$PValue.1800[idx] < 0.05) {
      sig_stars[r, "logFC.1800"] <- "*"
    }
  }

  # Create a unique file name for each subcluster
  outfile <- file.path(outdir, paste0("Subcluster_", i, "_heatmap.pdf"))
  
  pdf(outfile, width = 6, height = 8)  # Adjust width/height as needed
  pheatmap(
    sub_mat,
    cluster_rows    = TRUE,
    cluster_cols    = FALSE,
    display_numbers = sig_stars,
    number_color    = "black",
    labels_row      = new_labels,
    main            = paste("Subcluster", i, "Heatmap (with stars)"),
    scale           = "none"
  )
  dev.off()
  
  message("Saved: ", outfile)
}






###################
## ggplot barplot
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(stringr)

# For the up-regulated data:
# Convert the wide data frame to long format.
common_paths <- selected_pathways[selected_pathways %in% colnames(up_wide)[-1]]
up_wide <- dplyr::select(up_wide, time, all_of(common_paths))

up_long <- up_wide %>% 
  pivot_longer(cols = -time, names_to = "pathway", values_to = "neg_log10_p")

# Create the upregulated grouped bar plot.
up_gg <- ggplot(up_long, aes(x = pathway, y = neg_log10_p, fill = time)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("#E6E6FA", "#BA55D3", "#4B0082")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) + 
  labs(x = "Pathway", y = "-log10(pvalue)",
       title = "Up-regulated Pathway Enrichment (-log10(pvalue))") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold"))
  
# Open a new graphics window and display the plot.
dev.new()
print(up_gg)


# For the down-regulated data:
# Convert the wide data frame to long format.
common_paths <- selected_pathways[selected_pathways %in% colnames(down_wide)[-1]]
down_wide <- dplyr::select(down_wide, time, all_of(common_paths))

down_long <- down_wide %>% 
  pivot_longer(cols = -time, names_to = "pathway", values_to = "neg_log10_p")

# Create the down-regulated grouped bar plot.
down_gg <- ggplot(down_long, aes(x = pathway, y = neg_log10_p, fill = time)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("#E6E6FA", "#BA55D3", "#4B0082")) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 5)) + 
  labs(x = "Pathway", y = "-log10(pvalue)",
       title = "Down-regulated Pathway Enrichment (-log10(pvalue))") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold"))

      
# Open another new graphics window and display the down-regulated plot.
dev.new()
print(down_gg)




###################
## old plot
# Convert to matrix: Remove the "time" column and then assign rownames from it
up_mat <- as.matrix(up_wide[,-1])
rownames(up_mat) <- up_wide$time

# Define a color vector for each timepoint (adjust colors as needed)
time_colors <- c("#E6E6FA", "#BA55D3", "#4B0082")

dev.new()
# Create grouped barplot
barplot(up_mat, 
         beside = TRUE, 
         col = time_colors, 
         border = "white", 
         legend.text = rownames(up_mat), 
         xlab = "Pathway",
         main = "Up-regulated Pathway Enrichment (-log10(pvalue))",
         font.axis = 6,
         font.lab = 6,
         las = 2,         # makes axis labels vertical
         cex.names = 0.8) # reduces the font size of axis names

# Convert to matrix: Remove the "time" column and then assign rownames from it
down_mat <- as.matrix(down_wide[,-1])
rownames(down_mat) <- down_wide$time

# Define a color vector for each timepoint (adjust colors as needed)
time_colors <- c("#E6E6FA", "#BA55D3", "#4B0082")

dev.new()
# Create grouped barplot
# Create grouped barplot
barplot(down_mat, 
         beside = TRUE, 
         col = time_colors, 
         border = "white", 
         legend.text = rownames(up_mat), 
         xlab = "Pathway",
         main = "Down-regulated Pathway Enrichment (-log10(pvalue))",
         font.axis = 6,
         font.lab = 6,
         las = 2,         # makes axis labels vertical
         cex.names = 0.8) # reduces the font size of axis names