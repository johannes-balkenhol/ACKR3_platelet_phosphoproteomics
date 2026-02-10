####################################################################################
##### differential analysis
#### fitting
#dataset <- SummarizedExperiment::assay(ppe_imputed_tmp, "Quantification")
dataset <- SummarizedExperiment::assay(ppe, "normalised")
design <- model.matrix(~ grps - 1)
#design_interaction <- model.matrix(~0+time+treatment)
#f <- gsub("_\\d", "", colnames(ppe))
#X <- model.matrix(~ f - 1)
fit <- lmFit(dataset, design)


#treatment <- sapply(strsplit(colnames(ppe), "_"), "[[",2)
time <- sapply(strsplit(colnames(ppe), "_"), "[[",1)


#(treat10sec - ctrl0sec) - (ctrl10sec - ctrl0sec)

## tt10
#contrast.matrix <- makeContrasts(((grpsx10sek_CXCR7-grpsx0sek_Ctrl)-(grpsx10sek_DMSO-grpsx0sek_Ctrl)), levels=design)
contrast.matrix <- makeContrasts((grpsx10sek_CXCR7-grpsx10sek_DMSO), levels=design)
#contrast.matrix <- makeContrasts((grpsx10sek_DMSO-grpsx0sek_Ctrl), levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.10 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

##  tt600
contrast.matrix <- makeContrasts((grpsx600sek_CXCR7-grpsx600sek_DMSO), levels=design)
#contrast.matrix <- makeContrasts((grpsx600sek_DMSO-grpsx0sek_Ctrl), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.600 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

## tt1800
contrast.matrix <- makeContrasts(grpsx1800sek_CXCR7-grpsx1800sek_DMSO, levels=design)
#contrast.matrix <- makeContrasts(grpsx1800sek_DMSO-grpsx0sek_Ctrl, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.1800 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)



### annotation
uniprot <- sapply(strsplit(rownames(top.10), ";"), "[[", 1)
top.10 <- cbind(top.10, uniprot)

uniprot <- sapply(strsplit(rownames(top.600), ";"), "[[", 1)
top.600 <- cbind(top.600, uniprot)

uniprot <- sapply(strsplit(rownames(top.1800), ";"), "[[", 1)
top.1800 <- cbind(top.1800, uniprot)



#### significant table of different norm methods
top.10.sign <- top.10[top.10[, "adj.P.Val"] <0.05,]
top.600.sign <-  top.600[top.600[, "adj.P.Val"] <0.05,]
top.1800.sign <-  top.1800[top.1800[, "adj.P.Val"] <0.05,]

# top.10.sign.RUVIII <- top.10[top.10[, "adj.P.Val"] <0.05,]
# top.600.sign.RUVIII <-  top.600[top.600[, "adj.P.Val"] <0.05,]
# top.1800.sign.RUVIII <-  top.1800[top.1800[, "adj.P.Val"] <0.05,]
# 
# top.10.sign.RUVg <- top.10[top.10[, "adj.P.Val"] <0.05,]
# top.600.sign.RUVg <-  top.600[top.600[, "adj.P.Val"] <0.05,]
# top.1800.sign.RUVg <-  top.1800[top.1800[, "adj.P.Val"] <0.05,]
# 
# top.10.sign.com <- top.10[top.10[, "adj.P.Val"] <0.05,]
# top.600.sign.com <-  top.600[top.600[, "adj.P.Val"] <0.05,]
# top.1800.sign.com <-  top.1800[top.1800[, "adj.P.Val"] <0.05,]


#### statistics
### number of significant dergulated targets
DE2.RUV <- c(length(rownames(top.10.sign)),length(rownames(top.600.sign)),length(rownames(top.1800.sign)))
DE2.RUV


# DE2.RUV <- c(length(rownames(top.10.sign.RUVIII)),length(rownames(top.600.sign.RUVIII)),length(rownames(top.1800.sign.RUVIII)))
# DE2.RUV
# 
# 
# DE2.RUV <- c(length(rownames(top.10.sign.com )),length(rownames(top.600.sign.com )),length(rownames(top.1800.sign.com )))
# DE2.RUV



#### investigate specific proteins

top.10[grep(pattern = "*PLNM*", x = rownames(top.10)),]
top.10[grep(pattern = "VASP", x = rownames(top.10)),]
top.600["P50552;VASP;S239;KVSKQEEASGGPTAPK;9984",]
top.600.vasp <- top.600[top.900$uniprot=="P50552",]
top.10.vasp <- top.10[top.900$uniprot=="P50552",]
top.1800.vasp <- top.1800[top.1800$uniprot=="P50552",]


#### overlap of different normalization techniques
# common_differences <- Reduce(intersect, list(rownames(top.10.sign), rownames(top.600.sign), rownames(top.1800.sign)))
# common_differences <- Reduce(intersect, list(rownames(top.600.sign), rownames(top.1800.sign)))
# common_differences <- Reduce(intersect, list(rownames(top.10.sign.RUVIII), rownames(top.600.sign.RUVIII), rownames(top.1800.sign.RUVIII)))
# common_differences <- Reduce(intersect, list(rownames(top.600.sign.RUVIII), rownames(top.1800.sign.RUVIII)))
# common_differences <- Reduce(intersect, list(rownames(top.10.sign), rownames(top.10.sign.RUVIII)))
# common_differences <- Reduce(intersect, list(rownames(top.600.sign), rownames(top.600.sign.RUVIII)))
# common_differences <- Reduce(intersect, list(rownames(top.1800.sign), rownames(top.1800.sign.RUVIII)))
# common_differences <- Reduce(intersect, list(rownames(top.10.sign), rownames(top.10.sign.com)))
# common_differences <- Reduce(intersect, list(rownames(top.600.sign), rownames(top.600.sign.com)))
# common_differences <- Reduce(intersect, list(rownames(top.1800.sign), rownames(top.1800.sign.com)))
# common_differences <- Reduce(intersect, list(rownames(top.10.sign.com), rownames(top.10.sign.RUVIII)))


##### prepare export for excel and mysql

top.10$id <- rownames(top.10)
top.600$id <- rownames(top.600)
top.1800$id <- rownames(top.1800)


top.10$symbol <- sapply(strsplit(rownames(top.10), ";"), "[[", 2)
top.600$symbol <- sapply(strsplit(rownames(top.600), ";"), "[[", 2)
top.1800$symbol <- sapply(strsplit(rownames(top.1800), ";"), "[[", 2)


new_10 <- top.10[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.10[c("uniprot","logFC","P.Value","adj.P.Val")])), ]
new_600 <- top.600[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.600[c("uniprot","logFC","P.Value","adj.P.Val")])), ]
new_1800 <- top.1800[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.1800[c("uniprot","logFC","P.Value","adj.P.Val")])), ]


top.all <- as.data.frame(merge(new_600[,2:4], new_1800[,2:4], by="row.names"))
rownames(top.all) <- top.all$Row.names
top.all <- top.all[,-1]
top.all <- merge(new_10[,c(1,5,2,3,4)], top.all, by="row.names", all = TRUE)
rownames(top.all) <- top.all$Row.names
top.all <- top.all[,-1]
colnames(top.all) <- c("uniprot","symbol","logFC.10","P.Value.10","adj.P.Val.10","logFC.600","P.Value.600","adj.P.Val.600","logFC.1800","P.Value.1800","adj.P.Val.1800")
## read tables locally or
## RS: connect to mysql server and wrtie the table from the mysql database
#write.table(top.all, "../top.all.txt", sep="\t", , row.names=TRUE)
write.table(top.all, "../top.all_DMSOtoctrl0_ÖO_151122.txt", sep="\t", , row.names=TRUE)


#### transfer mysql database
#### export main tables 

## read tables locally or
## RS: connect to mysql server and wrtie the table from the mysql database
write.table(top.10, "../data/top.10.txt", sep="\t", , row.names=FALSE)
write.table(top.600, "../data/top.600.txt", sep="\t", , row.names=FALSE)
write.table(top.1800, "../data/top.1800.txt", sep="\t", , row.names=FALSE)