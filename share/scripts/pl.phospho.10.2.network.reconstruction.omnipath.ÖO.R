# Load packages ----
suppressPackageStartupMessages({ 
  library(clusterProfiler)
  library(dplyr)
  library(OmnipathR)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggVennDiagram)
  library(org.Hs.eg.db)
})

# Download protein interactions ----
#interactions <- import_omnipath_interactions(c("SignaLink3", "SIGNOR", "PhosphoSite"))

#post_translational i.e. physical interactions of proteins, protein-protein interactions (or PPIs)
interactions_PPI <- import_post_translational_interactions( 
  organism = 9606
)


# Perform quality control ----
qc_interactions <- interactions_PPI %>%
  filter(curation_effort > 2) %>%
  mutate(
    type = case_when(
      is_stimulation == 1 ~ "activation",
      is_inhibition == 1 ~ "inhibition",
      consensus_stimulation == 1 ~ "activation",
      consensus_inhibition == 1 ~ "inhibition",
      TRUE ~ NA_character_
    )
  ) %>%
  drop_na(type)


EGFR_interactions <- interactions_PPI %>%
  filter(source_genesymbol == "EGFR") %>%
  mutate(
    type = case_when(
      is_stimulation == 1 ~ "activation",
      is_inhibition == 1 ~ "inhibition",
      consensus_stimulation == 1 ~ "activation",
      consensus_inhibition == 1 ~ "inhibition",
      TRUE ~ NA_character_
    )
  ) %>%
  drop_na(type)




## load proteomics data ----
proteomicsData1 <- read.table("/Users/ozgeosmanoglu/Nextcloud/CXCR7_platelet_analysis/CXCR7_initial/proteomics/data/processed_data/top.all.txt", header = TRUE, sep = "\t")
proteomicsData2 <- read.table("/Users/ozgeosmanoglu/Nextcloud/CXCR7_platelet_analysis/CXCR7_validation/proteomics/data/processed_data/top.all.txt", header = TRUE, sep = "\t")


proteins1 <- rownames(proteomicsData1) %>%
  mapIds(org.Hs.eg.db,
         keys = .,
         keytype="UNIPROT",
         column="SYMBOL", multiVals = "first") %>%
  .[!is.na(.)] %>%
  as_tibble(rownames = "UNIPROT_ID")  %>%
  pull(value)

proteins2 <- rownames(proteomicsData2) %>%
  mapIds(org.Hs.eg.db,
         keys = .,
         keytype="UNIPROT",
         column="SYMBOL", multiVals = "first") %>%
  .[!is.na(.)] %>%
  as_tibble(rownames = "UNIPROT_ID")  %>%
  pull(value)

proteins <- unique(c(proteins1, proteins2))

#phosphodata

phosphoData1 <- read.table("/Users/ozgeosmanoglu/Nextcloud/CXCR7_platelet_analysis/CXCR7_initial/phosphoproteomics/data/processed_data/top.all.txt", header = TRUE, sep = "\t")
phosphoData2 <- read.table("/Users/ozgeosmanoglu/Nextcloud/CXCR7_platelet_analysis/CXCR7_validation/phosphoproteomics/data/processed_data/top.all.txt", header = TRUE, sep = "\t")

phosphoproteins1 <- unique(phosphoData1$symbol)
phosphoproteins2 <- unique(phosphoData2$symbol)

phosphoproteins <- unique(c(phosphoproteins1, phosphoproteins2))

ggVennDiagram(list(proteins1, proteins2, phosphoproteins1, phosphoproteins2))


#all together
plProteins <- unique(c(phosphoproteins, proteins)) 

## get platelet network ----
PlNetwork <- qc_interactions %>%
  filter(source_genesymbol %in% plProteins & 
           target_genesymbol %in% plProteins)

PlNetwork <- rbind(PlNetwork, EGFR_interactions)


# add edge weights

phosphoDataVal <- read.table("/Users/ozgeosmanoglu/Nextcloud/CXCR7_platelet_analysis/CXCR7_validation/phosphoproteomics/data/processed_data/top.all.collapse.txt", header = TRUE, sep = "\t",
                             dec = ",")

min(abs(as.numeric(phosphoDataVal$X10)))
min(abs(as.numeric(phosphoDataVal$X600)))
min(abs(as.numeric(phosphoDataVal$X1800)))



PlNetwork <- PlNetwork %>%
  mutate(
    source_logFC = as.numeric(phosphoDataVal$X10[match(source, phosphoDataVal$ID)]),
    target_logFC = as.numeric(phosphoDataVal$X10[match(target, phosphoDataVal$ID)])) %>%
  mutate(
    source_logFC = ifelse(is.na(source_logFC), 10^-6, source_logFC),
    target_logFC = ifelse(is.na(target_logFC), 10^-6, target_logFC),
    weight_10 = 1/abs((source_logFC * target_logFC)) #edge weights
  ) %>%
  select(-source_logFC, -target_logFC) %>%
  mutate(
    source_logFC = as.numeric(phosphoDataVal$X600[match(source, phosphoDataVal$ID)]),
    target_logFC = as.numeric(phosphoDataVal$X600[match(target, phosphoDataVal$ID)])) %>%
  mutate(
    source_logFC = ifelse(is.na(source_logFC), 10^-6, source_logFC),
    target_logFC = ifelse(is.na(target_logFC), 10^-6, target_logFC),
    weight_600 = 1/abs((source_logFC * target_logFC)) #edge weights
  ) %>%
  select(-source_logFC, -target_logFC) %>%
  mutate(
    source_logFC = as.numeric(phosphoDataVal$X1800[match(source, phosphoDataVal$ID)]),
    target_logFC = as.numeric(phosphoDataVal$X1800[match(target, phosphoDataVal$ID)])) %>%
  mutate(
    source_logFC = ifelse(is.na(source_logFC), 10^-6, source_logFC),
    target_logFC = ifelse(is.na(target_logFC), 10^-6, target_logFC),
    weight_1800 = 1/abs((source_logFC * target_logFC)) #edge weights
  ) %>%
  select(-source_logFC, -target_logFC)


options(scipen = 999)

write.table(PlNetwork, "/Users/ozgeosmanoglu/Nextcloud/CXCR7_platelet_analysis/CXCR7_validation/phosphoproteomics/analysis/intact/PLNetworkÖO.txt", 
            row.names = F, sep = "\t")


options(scipen = 0)
