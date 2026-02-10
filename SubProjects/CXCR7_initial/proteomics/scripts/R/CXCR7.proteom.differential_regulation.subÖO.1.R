## script by Johannes Balkenhol and Özge Osmanoglu
## differential analysis

# Load packages ----
suppressPackageStartupMessages({
  library(data.table)
  library(calibrate)
  library(annotate)
  library(clusterProfiler)
  library(sqldf)
  library(limma)
  library(gt)
  library(dplyr)
  library(tibble)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# Choose dataset ----
dataset <-  scaled_data
# Decide on the design matrix ----
design <- model.matrix(~0 + grps) #(~0 + grps + grps2) for including batch
# Fit model ----
fit <- lmFit(dataset, design)

# Choose contrasts and get statistics ----

## tt10 vs. tt00 ----
contrast.matrix <- makeContrasts((grps0010-grps0000), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.10 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

## tt30 vs. tt00 ----
contrast.matrix <- makeContrasts((grps0030-grps0000), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.30 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

## tt60 vs. tt00 ----
contrast.matrix <- makeContrasts((grps0060-grps0000), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.60 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

## tt300 vs. tt00 ----
contrast.matrix <- makeContrasts((grps0300-grps0000), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.300 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

##  tt600 vs. tt00 ----
contrast.matrix <- makeContrasts((grps0600-grps0000), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.600 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

## tt900 vs. tt00 ----
contrast.matrix <- makeContrasts((grps0900-grps0000), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.900 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

## tt1800 vs. tt00 ----
contrast.matrix <- makeContrasts(grps1800-grps0000, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.1800 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

# significant proteins ----
top.10.sign <- top.10[top.10[, "adj.P.Val"] <0.05,]
top.30.sign <- top.30[top.30[, "adj.P.Val"] <0.05,]
top.60.sign <- top.60[top.60[, "adj.P.Val"] <0.05,]
top.300.sign <- top.300[top.300[, "adj.P.Val"] <0.05,]
top.600.sign <-  top.600[top.600[, "adj.P.Val"] <0.05,]
top.900.sign <- top.900[top.900[, "adj.P.Val"] <0.05,]
top.1800.sign <-  top.1800[top.1800[, "adj.P.Val"] <0.05,]

length(unique(c(rownames(top.10.sign), rownames(top.30.sign), rownames(top.60.sign),
    rownames(top.300.sign), rownames(top.600.sign), rownames(top.900.sign),
    rownames(top.1800.sign))))

# Differential and significant proteins ----

# Function to filter significant differentially expressed genes (DEGs)
filter_significant <- function(df, adj_p_val_col = "adj.P.Val", 
                               logfc_col = "logFC", p_val_thresh = 0.05, logfc_thresh = 0.5) {
  df %>%
    filter(!!sym(adj_p_val_col) < p_val_thresh & abs(!!sym(logfc_col)) > logfc_thresh)
}

# Apply the filter to each dataset
datasets <- list(
  t10 = top.10, t30 = top.30, t60 = top.60, t300= top.300,
  t600 = top.600, t900 = top.900, t1800 = top.1800)

# Get the significant and differential proteins
significantD_datasets <- lapply(datasets, function(df) 
  filter_significant(df,adj_p_val_col = "P.Value", logfc_thresh = 0.5))
significant_datasets <- lapply(datasets, function(df) 
  filter_significant(df, adj_p_val_col = "P.Value", logfc_thresh = 0))


# Extract row names for each set of significant DEGs
extract_rownames <- function(dfs) unique(unlist(lapply(dfs, rownames)))

signDnames <- extract_rownames(significantD_datasets)

# Prepare publication-ready tables ----

SignDiff_tb <- top.all[signDnames, c(1, 4, 7, 10, 13, 16, 19)] %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  arrange(desc(rowMeans(dplyr::across(2:8, abs)))) %>%
  mutate(symbol = mapIds(org.Hs.eg.db, keys = gene, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first"),
         description = mapIds(org.Hs.eg.db, keys = gene, column = "GENENAME", keytype = "UNIPROT", multiVals = "first")
  ) %>%
  dplyr::select(gene, symbol, description, 2:8)

# Manually updating specific rows
custom_updates <- tibble(
  gene = c("A6NIZ1", "A0A0C4DH25", "A0A0B4J1V0", "A0A0B4J1U7"),
  symbol = c("RAP1BL", "IGKV3D-20", "IGHV3-15", "IGHV6-1"),
  description = c("Ras-related protein Rap-1b-like protein", "Immunoglobulin kappa variable 3D-20", "Immunoglobulin heavy variable 3-15", "Immunoglobulin heavy variable 6-1")
)

SignDiff_tb <- SignDiff_tb %>%
  dplyr::left_join(custom_updates, by = "gene") %>%
  mutate(symbol = dplyr::coalesce(symbol.y, symbol.x),
         description = dplyr::coalesce(description.y, description.x)) %>%
  dplyr::select(-symbol.x, -symbol.y, -description.x, -description.y) %>%
  dplyr::select(gene, symbol, description, 2:8)


write.table(SignDiff_tb, "../data/processed_data/allSignDEGs.txt", row.names = F, sep = "\t")

number = 10

gt(SignDiff_tb[1:number,2:10]) %>% 
  tab_options(table.font.size = 14) %>%
  fmt_number(columns=c(3:9), decimals = 1) %>%
  cols_width(
    symbol ~ px(200),
    description ~ px(800)) %>%
  cols_label(#gene = md("**Uniprot**"),
             symbol = md("**Symbol**"),
             description = md("**Name**"),
             logFC.10 = md("**t10**"),
             logFC.30 = md("**t30**"),
             logFC.60 = md("**t60**"),
             logFC.300 = md("**t300**"),
             logFC.600 = md("**t600**"),
             logFC.900 = md("**t900**"),
             logFC.1800 = md("**t1800**")) %>%
  cols_align("center",
             columns = 3:9) %>%
  tab_header(title = md(paste0("Top ", number,  " differentially regulated proteins")),
             #subtitle = md("**to choose fit type***")
  ) %>%
  gtsave(paste0("../data/processed_data/", number, "signDEGs.pdf"))




# Statistics ----
### number of significant dergulated targets
DE2.RUV2 <- c(length(rownames(top.10.sign)),length(rownames(top.30.sign)),length(rownames(top.60.sign)),
              length(rownames(top.300.sign)),length(rownames(top.600.sign)),length(rownames(top.900.sign)),
              length(rownames(top.1800.sign)))
DE2.RUV2 #[1] 0 0 0 0 0 0 0


# Prepare top tables for export ----

top.10$id <- rownames(top.10)
top.30$id <- rownames(top.30)
top.60$id <- rownames(top.60)
top.300$id <- rownames(top.300)
top.600$id <- rownames(top.600)
top.900$id <- rownames(top.900)
top.1800$id <- rownames(top.1800)


top.10$symbol <- mapIds(org.Hs.eg.db, keys=rownames(top.10), column="SYMBOL", keytype="UNIPROT", multiVals="first")
top.30$symbol <- mapIds(org.Hs.eg.db, keys=rownames(top.30), column="SYMBOL", keytype="UNIPROT", multiVals="first")
top.60$symbol <- mapIds(org.Hs.eg.db, keys=rownames(top.60), column="SYMBOL", keytype="UNIPROT", multiVals="first")
top.300$symbol <- mapIds(org.Hs.eg.db, keys=rownames(top.300), column="SYMBOL", keytype="UNIPROT", multiVals="first")
top.600$symbol <- mapIds(org.Hs.eg.db, keys=rownames(top.600), column="SYMBOL", keytype="UNIPROT", multiVals="first")
top.900$symbol <- mapIds(org.Hs.eg.db, keys=rownames(top.900), column="SYMBOL", keytype="UNIPROT", multiVals="first")
top.1800$symbol <- mapIds(org.Hs.eg.db, keys=rownames(top.1800), column="SYMBOL", keytype="UNIPROT", multiVals="first")


new_10 <- top.10[c("logFC","P.Value","adj.P.Val","id")][ order(row.names(top.10[c("logFC","P.Value","adj.P.Val")])), ]
new_30 <- top.30[c("logFC","P.Value","adj.P.Val","id")][ order(row.names(top.30[c("logFC","P.Value","adj.P.Val")])), ]
new_60 <- top.60[c("logFC","P.Value","adj.P.Val","id")][ order(row.names(top.60[c("logFC","P.Value","adj.P.Val")])), ]
new_300 <- top.300[c("logFC","P.Value","adj.P.Val","id")][ order(row.names(top.300[c("logFC","P.Value","adj.P.Val")])), ]
new_600 <- top.600[c("logFC","P.Value","adj.P.Val","id")][ order(row.names(top.600[c("logFC","P.Value","adj.P.Val")])), ]
new_900 <- top.900[c("logFC","P.Value","adj.P.Val","id")][ order(row.names(top.900[c("logFC","P.Value","adj.P.Val")])), ]
new_1800 <- top.1800[c("logFC","P.Value","adj.P.Val","id")][ order(row.names(top.1800[c("logFC","P.Value","adj.P.Val")])), ]


top.all <- as.data.frame(merge(new_900[,1:3], new_1800[,1:3], by="row.names"))
rownames(top.all) <- top.all$Row.names
top.all <- top.all[,-1]


top.all <- merge(new_600[,1:3], top.all, by="row.names", all = TRUE)
colnames(top.all) <- c("Row.names", "logFC.600","P.Value.600","adj.P.Val.600","logFC.900","P.Value.900","adj.P.Val.900", "logFC.1800","P.Value.1800","adj.P.Val.1800")
rownames(top.all) <- top.all$Row.names
top.all <- top.all[,-1]


top.all <- merge(new_300[,1:3], top.all, by="row.names", all = TRUE)
colnames(top.all) <- c("Row.names", "logFC.300","P.Value.300","adj.P.Val.300", "logFC.600","P.Value.600","adj.P.Val.600","logFC.900","P.Value.900","adj.P.Val.900","logFC.1800","P.Value.1800","adj.P.Val.1800")
rownames(top.all) <- top.all$Row.names
top.all <- top.all[,-1]


top.all <- merge(new_60[,1:3], top.all, by="row.names", all = TRUE)
colnames(top.all) <- c("Row.names", "logFC.60","P.Value.60","adj.P.Val.60", "logFC.300","P.Value.300","adj.P.Val.300", "logFC.600","P.Value.600","adj.P.Val.600","logFC.900","P.Value.900","adj.P.Val.900", "logFC.1800","P.Value.1800","adj.P.Val.1800")
rownames(top.all) <- top.all$Row.names
top.all <- top.all[,-1]


top.all <- merge(new_30[,1:3], top.all, by="row.names", all = TRUE)
colnames(top.all) <- c("Row.names", "logFC.30","P.Value.30","adj.P.Val.30", "logFC.60","P.Value.60","adj.P.Val.60", "logFC.300","P.Value.300","adj.P.Val.300", "logFC.600","P.Value.600","adj.P.Val.600","logFC.900","P.Value.900","adj.P.Val.900", "logFC.1800","P.Value.1800","adj.P.Val.1800")
rownames(top.all) <- top.all$Row.names
top.all <- top.all[,-1]


top.all <- merge(new_10[,c(1:3)], top.all, by="row.names", all = TRUE)
rownames(top.all) <- top.all$Row.names
top.all <- top.all[,-1]
colnames(top.all) <- c("logFC.10","P.Value.10","adj.P.Val.10", "logFC.30","P.Value.30","adj.P.Val.30", "logFC.60","P.Value.60","adj.P.Val.60", "logFC.300","P.Value.300","adj.P.Val.300", "logFC.600","P.Value.600","adj.P.Val.600","logFC.900","P.Value.900","adj.P.Val.900", "logFC.1800","P.Value.1800","adj.P.Val.1800")

## Save top.all ----
write.table(top.all, "../data/processed_data/top.all.txt", sep="\t", , row.names=TRUE)

## Save individial toptables ----
write.table(top.10, "../data/processed_data/top.10.txt", sep="\t", , row.names=FALSE)
write.table(top.30, "../data/processed_data/top.30.txt", sep="\t", , row.names=FALSE)
write.table(top.60, "../data/processed_data/top.60.txt", sep="\t", , row.names=FALSE)
write.table(top.300, "../data/processed_data/top.300.txt", sep="\t", , row.names=FALSE)
write.table(top.600, "../data/processed_data/top.600.txt", sep="\t", , row.names=FALSE)
write.table(top.900, "../data/processed_data/top.900.txt", sep="\t", , row.names=FALSE)
write.table(top.1800, "../data/processed_data/top.1800.txt", sep="\t", , row.names=FALSE)

# Make log2FC histogram ----

log2data <- top.all[c(1,4,7,10,13,16,19)]

tiff(filename = "../analysis/Volcano_plots/log2FCHistograms.tiff",
     width = 12* 300, 
     height = 6 * 300,
     res = 300,
     compression = "lzw")
par(mfrow=c(2,4))
for(i in 1:ncol(log2data)){
  hist(log2data[,i], na.rm=T, main=names(log2data)[i],
       breaks=seq(-10,10,0.1), xlim= c(-3,3), ylim = c(0,1000))
}
dev.off()
  

