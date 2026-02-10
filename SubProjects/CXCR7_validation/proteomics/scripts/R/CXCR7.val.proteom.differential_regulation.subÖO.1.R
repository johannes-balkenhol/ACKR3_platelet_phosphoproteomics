####################################################################################
##### differential analysis
### define the groups and replicates
### literature
### https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2021/RNAseq/Markdowns/09_Linear_Models.html
### http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
### http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments


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
design <- model.matrix(~0 + grps) #(~0 + grps + grps4) for including batch
# Fit model ----
fit <- lmFit(dataset, design)

# Choose contrasts and get statistics ----

## CXCR7 vs. DMSO ----
### tt10 ----
contrast.matrix <- makeContrasts((grps0010_CXCR7 - grps0010_DMSO), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.10 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

### tt600 ----
contrast.matrix <- makeContrasts((grps0600_CXCR7 - grps0600_DMSO), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.600 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

### tt1800 ----
contrast.matrix <- makeContrasts(grps1800_CXCR7 - grps1800_DMSO, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.1800 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

## DMSO vs. Ctrl (0s) ----
### tt10 ----
contrast.matrix <- makeContrasts((grps0010_DMSO - grps0000_Ctrl), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.10.dmso.vs.0s <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

###  tt600 ----
contrast.matrix <- makeContrasts((grps0600_DMSO - grps0000_Ctrl), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.600.dmso.vs.0s <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

###  tt1800 ----
contrast.matrix <- makeContrasts(grps1800_DMSO - grps0000_Ctrl, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.1800.dmso.vs.0s <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

## CXCR7 vs. Ctrl (0s) ----
### tt10 ----
contrast.matrix <- makeContrasts((grps0010_CXCR7 - grps0000_Ctrl), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.10.cxcr7.vs.0s <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

###  tt600 ----
contrast.matrix <- makeContrasts((grps0600_CXCR7 - grps0000_Ctrl), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.600.cxcr7.vs.0s <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

###  tt1800 ----
contrast.matrix <- makeContrasts(grps1800_CXCR7 - grps0000_Ctrl, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.1800.cxcr7.vs.0s <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

##### Add gene names to top tables ----
top.10$id <- rownames(top.10)
top.600$id <- rownames(top.600)
top.1800$id <- rownames(top.1800)
top.10.dmso.vs.0s$id <- rownames(top.10.dmso.vs.0s)
top.600.dmso.vs.0s$id <- rownames(top.600.dmso.vs.0s)
top.1800.dmso.vs.0s$id <- rownames(top.1800.dmso.vs.0s)
top.10.cxcr7.vs.0s$id <- rownames(top.10.cxcr7.vs.0s)
top.600.cxcr7.vs.0s$id <- rownames(top.600.cxcr7.vs.0s)
top.1800.cxcr7.vs.0s$id <- rownames(top.1800.cxcr7.vs.0s)

top.10$symbol <- mapIds(org.Hs.eg.db, keys=rownames(top.10), column="SYMBOL", keytype="UNIPROT", multiVals="first")
top.600$symbol <- mapIds(org.Hs.eg.db, keys=rownames(top.600), column="SYMBOL", keytype="UNIPROT", multiVals="first")
top.1800$symbol <- mapIds(org.Hs.eg.db, keys=rownames(top.1800), column="SYMBOL", keytype="UNIPROT", multiVals="first")
top.10.cxcr7.vs.0s$symbol <- mapIds(org.Hs.eg.db, keys=rownames(top.10.cxcr7.vs.0s), column="SYMBOL", keytype="UNIPROT", multiVals="first")
top.600.cxcr7.vs.0s$symbol <- mapIds(org.Hs.eg.db, keys=rownames(top.600.cxcr7.vs.0s), column="SYMBOL", keytype="UNIPROT", multiVals="first")
top.1800.cxcr7.vs.0s$symbol <- mapIds(org.Hs.eg.db, keys=rownames(top.1800.cxcr7.vs.0s), column="SYMBOL", keytype="UNIPROT", multiVals="first")
top.10.dmso.vs.0s$symbol <- mapIds(org.Hs.eg.db, keys=rownames(top.10.dmso.vs.0s), column="SYMBOL", keytype="UNIPROT", multiVals="first")
top.600.dmso.vs.0s$symbol <- mapIds(org.Hs.eg.db, keys=rownames(top.600.dmso.vs.0s), column="SYMBOL", keytype="UNIPROT", multiVals="first")
top.1800.dmso.vs.0s$symbol <- mapIds(org.Hs.eg.db, keys=rownames(top.1800.dmso.vs.0s), column="SYMBOL", keytype="UNIPROT", multiVals="first")

# Differential and significant tables ----

# Function to filter significant differentially expressed genes (DEGs)
filter_significant <- function(df, adj_p_val_col = "adj.P.Val", logfc_col = "logFC", p_val_thresh = 0.05, logfc_thresh = 0.5) {
  df %>%
    filter(!!sym(adj_p_val_col) < p_val_thresh & abs(!!sym(logfc_col)) > logfc_thresh)
}

# Apply the filter to each dataset
datasets <- list(
  t10 = top.10, t600 = top.600, t1800 = top.1800,
  t10.cxcr7.vs.0s = top.10.cxcr7.vs.0s, t600.cxcr7.vs.0s = top.600.cxcr7.vs.0s, t1800.cxcr7.vs.0s = top.1800.cxcr7.vs.0s,
  t10.dmso.vs.0s = top.10.dmso.vs.0s, t600.dmso.vs.0s = top.600.dmso.vs.0s, t1800.dmso.vs.0s = top.1800.dmso.vs.0s
)

# Get the significant and differential proteins
significantD_datasets <- lapply(datasets, function(df) 
  filter_significant(df,adj_p_val_col = "P.Value", logfc_thresh = 0.5))
significant_datasets <- lapply(datasets, function(df) 
  filter_significant(df, adj_p_val_col = "P.Value", logfc_thresh = 0))


# Extract row names for each set of significant DEGs
extract_rownames <- function(dfs) unique(unlist(lapply(dfs, rownames)))

signDnames <- extract_rownames(significantD_datasets[1:3])
signDnames.cxcr7.vs.0s <- extract_rownames(significantD_datasets[4:6])
signDnames.dmso.vs.0s <- extract_rownames(significantD_datasets[7:9])

signnames <- extract_rownames(significant_datasets[1:3])
signnames.cxcr7.vs.0s <- extract_rownames(significant_datasets[4:6])
signnames.dmso.vs.0s <- extract_rownames(significant_datasets[7:9])

# Create input list for each comparison
inputList <- list(
  CXCR7.vs.DMSO = list(signDnames, top.all),
  CXCR7.vs.0 = list(signDnames.cxcr7.vs.0s, top.all.cxcr7.vs.0s),
  DMSO.vs.0 = list(signDnames.dmso.vs.0s, top.all.dmso.vs.0s)
)
# Function to process and save each comparison
process_comparison <- function(signDnames, top_all, comparison_name, number) {
  SignDiff_tb <- top_all[signDnames[1:number], c(1, 4, 7)] %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    arrange(desc(rowMeans(dplyr::across(2:4, abs)))) %>%
    mutate(
      symbol = mapIds(org.Hs.eg.db, keys = gene, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first"),
      description = mapIds(org.Hs.eg.db, keys = gene, column = "GENENAME", keytype = "UNIPROT", multiVals = "first")
    ) %>%
    dplyr::select(gene, symbol, description, 2:4)
  
  # Write to file
  write.table(SignDiff_tb, paste0("../data/processed_data/", number, "signDEGs_", comparison_name, ".txt"), row.names = FALSE, sep = "\t")
  
  # Create and save the GT table
  number <- nrow(SignDiff_tb)
  gt(SignDiff_tb[1:number, 2:6]) %>%
    tab_options(table.font.size = 14) %>%
    fmt_number(columns = 3:5, decimals = 1) %>%
    cols_label(
      symbol = md("**Symbol**"),
      description = md("**Name**"),
      logFC.10 = md("**t10**"),
      logFC.600 = md("**t600**"),
      logFC.1800 = md("**t1800**")
    ) %>%
    cols_align("center", columns = 3:5) %>%
    tab_header(
      title = md(paste0("Top ", number, " proteins (", gsub("\\.", " ", comparison_name), ")"))
    ) %>%
    gtsave(paste0("../data/processed_data/", number, "signDEGs_", comparison_name, ".pdf"))
}

# Process each comparison
for (i in seq_along(inputList)) {
  comparison_name <- names(inputList)[i]
  signDnames <- inputList[[i]][[1]]
  top_all <- inputList[[i]][[2]]
  process_comparison(signDnames, top_all, comparison_name, 10)
  process_comparison(signDnames, top_all, comparison_name, length(signDnames))
}


# statistics ----
### number of significant dergulated targets
input <- significant_datasets
DE2.RUV <- c(nrow(input[[1]]),nrow(input[[2]]),nrow(input[[3]]),
             nrow(input[[4]]),nrow(input[[5]]),nrow(input[[6]]),
             nrow(input[[7]]),nrow(input[[8]]),nrow(input[[9]]))

input <- significantD_datasets
DE2.RUV2 <- c(nrow(input[[1]]),nrow(input[[2]]),nrow(input[[3]]),
             nrow(input[[4]]),nrow(input[[5]]),nrow(input[[6]]),
             nrow(input[[7]]),nrow(input[[8]]),nrow(input[[9]]))

DE2.RUV    # [1] 36 25 29 28 27 33 25 14 17
DE2.RUV2   # [1] 19 15 21 16 18 18 17  6 12



# Prepare top tables for export ----

new_10 <- top.10[c("logFC","P.Value","adj.P.Val","id")][ order(row.names(top.10[c("logFC","P.Value","adj.P.Val")])), ]
new_600 <- top.600[c("logFC","P.Value","adj.P.Val","id")][ order(row.names(top.600[c("logFC","P.Value","adj.P.Val")])), ]
new_1800 <- top.1800[c("logFC","P.Value","adj.P.Val","id")][ order(row.names(top.1800[c("logFC","P.Value","adj.P.Val")])), ]
new_10.dmso.vs.0s <- top.10.dmso.vs.0s[c("logFC","P.Value","adj.P.Val","id")][ order(row.names(top.10.dmso.vs.0s[c("logFC","P.Value","adj.P.Val")])), ]
new_600.dmso.vs.0s <- top.600.dmso.vs.0s[c("logFC","P.Value","adj.P.Val","id")][ order(row.names(top.600.dmso.vs.0s[c("logFC","P.Value","adj.P.Val")])), ]
new_1800.dmso.vs.0s <- top.1800.dmso.vs.0s[c("logFC","P.Value","adj.P.Val","id")][ order(row.names(top.1800.dmso.vs.0s[c("logFC","P.Value","adj.P.Val")])), ]
new_10.cxcr7.vs.0s <- top.10.cxcr7.vs.0s[c("logFC","P.Value","adj.P.Val","id")][ order(row.names(top.10.cxcr7.vs.0s[c("logFC","P.Value","adj.P.Val")])), ]
new_600.cxcr7.vs.0s <- top.600.cxcr7.vs.0s[c("logFC","P.Value","adj.P.Val","id")][ order(row.names(top.600.cxcr7.vs.0s[c("logFC","P.Value","adj.P.Val")])), ]
new_1800.cxcr7.vs.0s <- top.1800.cxcr7.vs.0s[c("logFC","P.Value","adj.P.Val","id")][ order(row.names(top.1800.cxcr7.vs.0s[c("logFC","P.Value","adj.P.Val")])), ]

## CXCR7 vs. DMSO (0s) ----
top.all <- as.data.frame(merge(new_600[,1:3], new_1800[,1:3], by="row.names"))
rownames(top.all) <- top.all$Row.names
top.all <- top.all[,-1]
top.all <- merge(new_10[,1:3], top.all, by="row.names", all = TRUE)
rownames(top.all) <- top.all$Row.names
top.all <- top.all[,-1]
colnames(top.all) <- c("logFC.10","P.Value.10","adj.P.Val.10","logFC.600","P.Value.600","adj.P.Val.600","logFC.1800","P.Value.1800","adj.P.Val.1800")

write.table(top.all, "../data/processed_data/top.all.txt", sep="\t", , row.names=TRUE)

## DMSO vs. Ctrl (0s) ----
top.all.dmso.vs.0s <- as.data.frame(merge(new_600.dmso.vs.0s[,1:3], new_1800.dmso.vs.0s[,1:3], by="row.names"))
rownames(top.all.dmso.vs.0s) <- top.all.dmso.vs.0s$Row.names
top.all.dmso.vs.0s <- top.all.dmso.vs.0s[,-1]
top.all.dmso.vs.0s <- merge(new_10.dmso.vs.0s[,1:3], top.all.dmso.vs.0s, by="row.names", all = TRUE)
rownames(top.all.dmso.vs.0s) <- top.all.dmso.vs.0s$Row.names
top.all.dmso.vs.0s <- top.all.dmso.vs.0s[,-1]
colnames(top.all.dmso.vs.0s) <- c("logFC.10","P.Value.10","adj.P.Val.10","logFC.600","P.Value.600","adj.P.Val.600","logFC.1800","P.Value.1800","adj.P.Val.1800")

write.table(top.all.dmso.vs.0s, "../data/processed_data/top.all.dmso.vs.0s.txt", sep="\t", , row.names=TRUE)

## CXCR7 vs. Ctrl (0s) ----
top.all.cxcr7.vs.0s <- as.data.frame(merge(new_600.cxcr7.vs.0s[,1:3], new_1800.cxcr7.vs.0s[,1:3], by="row.names"))
rownames(top.all.cxcr7.vs.0s) <- top.all.cxcr7.vs.0s$Row.names
top.all.cxcr7.vs.0s <- top.all.cxcr7.vs.0s[,-1]
top.all.cxcr7.vs.0s <- merge(new_10.cxcr7.vs.0s[,1:3], top.all.cxcr7.vs.0s, by="row.names", all = TRUE)
rownames(top.all.cxcr7.vs.0s) <- top.all.cxcr7.vs.0s$Row.names
top.all.cxcr7.vs.0s <- top.all.cxcr7.vs.0s[,-1]
colnames(top.all.cxcr7.vs.0s) <- c("logFC.10","P.Value.10","adj.P.Val.10","logFC.600","P.Value.600","adj.P.Val.600","logFC.1800","P.Value.1800","adj.P.Val.1800")

write.table(top.all.cxcr7.vs.0s, "../data/processed_data/top.all.cxcr7.vs.0s.txt", sep="\t", , row.names=TRUE)


# Save individial toptables ----

## CXCR7 vs. DMSO (0s) ----
write.table(top.10, "../data/processed_data/top.10.txt", sep="\t", , row.names=TRUE)
write.table(top.600, "../data/processed_data/top.600.txt", sep="\t", , row.names=TRUE)
write.table(top.1800, "../data/processed_data/top.1800.txt", sep="\t", , row.names=TRUE)

## DMSO vs. Ctrl (0s) ----
write.table(top.10.dmso.vs.0s, "../data/processed_data/top.10.dmso.vs.0s.txt", sep="\t", , row.names=TRUE)
write.table(top.600.dmso.vs.0s, "../data/processed_data/top.600.dmso.vs.0s.txt", sep="\t", , row.names=TRUE)
write.table(top.1800.dmso.vs.0s, "../data/processed_data/top.1800.dmso.vs.0s.txt", sep="\t", , row.names=TRUE)

## CXCR7 vs. Ctrl (0s) ----
write.table(top.10.cxcr7.vs.0s, "../data/processed_data/top.10.cxcr7.vs.0s.txt", sep="\t", , row.names=TRUE)
write.table(top.600.cxcr7.vs.0s, "../data/processed_data/top.600.cxcr7.vs.0s.txt", sep="\t", , row.names=TRUE)
write.table(top.1800.cxcr7.vs.0s, "../data/processed_data/top.1800.cxcr7.vs.0s.txt", sep="\t", , row.names=TRUE)

# Make log2FC histogram ----

log2data <- top.all[c(1,4,7)]

tiff(filename = "../analysis/Volcano_plots/log2FCHist_CXCR7vsDMSO.tiff",
     width = 9* 300, 
     height = 3 * 300,
     res = 300,
     compression = "lzw")
par(mfrow=c(1,3))
for(i in 1:ncol(log2data)){
  hist(log2data[,i], na.rm=T, main=names(log2data)[i],
       breaks=seq(-10,10,0.1), xlim= c(-3,3), ylim = c(0, 1000))
}
dev.off()
  
