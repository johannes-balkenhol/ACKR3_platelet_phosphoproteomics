library(dplyr)

### for 

## prepare Proteomics file
selected_cols <- top.10 %>%
  select(symbol, psite, logFC) %>%
  mutate(effect = "", ID = rownames(.))

rearranged_cols <- selected_cols %>%
  select(ID, symbol, psite, effect, logFC)

colnames(rearranged_cols) <- c("ID", "Symbols", "Sites", "Effect", "Value")


# Transpose the data frame
#transposed_df <- t(rearranged_cols)

#transposed_df2 <- cbind(rownames(transposed_df), transposed_df)
#colnames(transposed_df2)[1] <- "column_header"

rearranged_cols2 <- rbind(colnames(rearranged_cols), rearranged_cols)

# Export as tab-delimited file
write.table(rearranged_cols2, file = "../analysis/Causalpath/phospho_data.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

## prepare Parameter file

parameter-name = parameter-value

proteomics-values-file

proteomics-repeat-values-file

proteomics-platform-file

id-column

symbols-column

sites-column

feature-column

effect-column

value-transformation

value-column




## prepare Proteomics file
#norm_intensity
#grps


mean_by_sample <- t(aggregate(t(norm_intensity), by = list(groups = grps), FUN = mean))
colnames(mean_by_sample) <- mean_by_sample[1,]
mean_by_sample <- as.data.frame(mean_by_sample[-1,])
#mean_by_sample$ID <- rownames(mean_by_sample)

effect <- vector("numeric", nrow(mean_by_sample))
mean_by_sample <- cbind(Effect = effect, mean_by_sample)
mean_by_sample <- cbind(Sites = sapply(strsplit(rownames(mean_by_sample), ";"), "[[", 3), mean_by_sample)
mean_by_sample <- cbind(Symbols = sapply(strsplit(rownames(mean_by_sample), ";"), "[[", 2), mean_by_sample)
mean_by_sample <- cbind(Uniprot_id = sapply(strsplit(rownames(mean_by_sample), ";"), "[[", 1), mean_by_sample)
mean_by_sample <- cbind(ID = rownames(mean_by_sample), mean_by_sample)
mean_by_sample$ID <- gsub("\\|", "-", mean_by_sample$ID)
mean_by_sample2 <- rbind(colnames = colnames(mean_by_sample), mean_by_sample)
mean_by_sample2 <- mean_by_sample2[, !(colnames(mean_by_sample2) == "x0sek_Ctrl")]

mean_by_sample2$ID <- sapply(strsplit(mean_by_sample2$ID, split=";"), function(x) paste(x[1:3], collapse="_"))



## mapping ID (uniprot to HGNC)

# https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
head(datasets)


ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
attributes = listAttributes(ensembl)
attributes[1:5,]

uniprot_id <- sapply(strsplit(rownames(mean_by_sample), ";"), "[[", 1)

uniprot2hgnc <- getBM(attributes = c('uniprotswissprot', 'hgnc_symbol'),
      filters = 'uniprotswissprot',
      values = uniprot_id, 
      mart = ensembl)

uniprot2hgnc <- getBM(attributes = c('uniprot_gn_id', 'hgnc_symbol'),
      filters = 'uniprot_gn_id',
      values = uniprot_id, 
      mart = ensembl)

colnames(uniprot2hgnc) = c("Uniprot_id", "hgnc_symbol")


# Merge the data frames by uniprot_id and update symbols with hgnc_symbol
df_merged <- merge(mean_by_sample2, uniprot2hgnc, by = "Uniprot_id", all.x = TRUE)
df_merged$symbols <- ""
df_merged$symbols <- df_merged$hgnc_symbol
df_merged <- subset(df_merged, select = -Uniprot_id)
df_merged$Effect <- ""
df_merged <- df_merged[!is.na(df_merged$symbols) & df_merged$symbols != "", ]


write.table(df_merged, file = "../analysis/Causalpath/phospho_data_v2.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)




## search data
pattern <- "JNK"
matched_rows <- df_merged[grepl(pattern, df_merged$Symbols), ]
View(matched_rows)





# BioPax
# https://www.bioconductor.org/packages/devel/bioc/vignettes/paxtoolsr/inst/doc/using_paxtoolsr.html







