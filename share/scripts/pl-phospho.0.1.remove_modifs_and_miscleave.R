###removing miscleaved peptides and other modifications####
### a script by Özge Osmanoglu###


library(tidyr)
library(stringr)
library(dplyr)




peptide_info_split <- peptide_info %>%
  separate_rows(Modifications, sep = "];") %>%
  mutate(Modifications = str_trim(Modifications)) %>%
  separate(Modifications, into = c("Count", "Modification", "Location"), sep = "(?<=[0-9])x| \\[")

counts <- data.frame(table(peptide_info_split$peptide_id))

only_phospho_index <- counts[counts$Freq == 1, "Var1" ]

# first remove miscleaved ones 8291 --> 5697
peptide_info_filtered <- peptide_info_split[peptide_info_split$missed_cleavages == 0,]

#then remove anything that does not only have phosphorylation 5697 --Y 4824
peptide_info_filtered <- peptide_info_filtered[peptide_info_filtered$peptide_id %in% only_phospho_index,]


raw_abundance2$pepid <- sapply(strsplit(rownames(raw_abundance2), ";"), "[[", 5)

raw_abundance2_filtered <- raw_abundance2[raw_abundance2$pepid %in% peptide_info_filtered$peptide_id,]
raw_abundance2_filtered <-  raw_abundance2_filtered %>% select(-pepid)
raw_abundance2 <-  raw_abundance2 %>% select(-pepid)
