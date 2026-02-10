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
#time <- sapply(strsplit(colnames(ppe), "_"), "[[",1)


#(treat10sec - ctrl0sec) - (ctrl10sec - ctrl0sec)

## tt10 vs. tt00
#contrast.matrix <- makeContrasts(((grpsx10sek_CXCR7-grpsx0sek_Ctrl)-(grpsx10sek_DMSO-grpsx0sek_Ctrl)), levels=design)
contrast.matrix <- makeContrasts((grpsX0010-grpsX0000), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.10 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)


## tt30 vs. tt00
#contrast.matrix <- makeContrasts(((grpsx10sek_CXCR7-grpsx0sek_Ctrl)-(grpsx10sek_DMSO-grpsx0sek_Ctrl)), levels=design)
contrast.matrix <- makeContrasts((grpsX0030-grpsX0000), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.30 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

## tt60 vs. tt00
#contrast.matrix <- makeContrasts(((grpsx10sek_CXCR7-grpsx0sek_Ctrl)-(grpsx10sek_DMSO-grpsx0sek_Ctrl)), levels=design)
contrast.matrix <- makeContrasts((grpsX0060-grpsX0000), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.60 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

## tt300 vs. tt00
#contrast.matrix <- makeContrasts(((grpsx10sek_CXCR7-grpsx0sek_Ctrl)-(grpsx10sek_DMSO-grpsx0sek_Ctrl)), levels=design)
contrast.matrix <- makeContrasts((grpsX0300-grpsX0000), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.300 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

##  tt600 vs. tt00
contrast.matrix <- makeContrasts((grpsX0600-grpsX0000), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.600 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

## tt900 vs. tt00
#contrast.matrix <- makeContrasts(((grpsx10sek_CXCR7-grpsx0sek_Ctrl)-(grpsx10sek_DMSO-grpsx0sek_Ctrl)), levels=design)
contrast.matrix <- makeContrasts((grpsX0900-grpsX0000), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.900 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

## tt1800 vs. tt00
contrast.matrix <- makeContrasts(grpsX1800-grpsX0000, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.1800 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)



### annotation
uniprot <- sapply(strsplit(rownames(top.10), ";"), "[[", 1)
top.10 <- cbind(top.10, uniprot)

uniprot <- sapply(strsplit(rownames(top.30), ";"), "[[", 1)
top.30 <- cbind(top.30, uniprot)

uniprot <- sapply(strsplit(rownames(top.60), ";"), "[[", 1)
top.60 <- cbind(top.60, uniprot)

uniprot <- sapply(strsplit(rownames(top.300), ";"), "[[", 1)
top.300 <- cbind(top.300, uniprot)

uniprot <- sapply(strsplit(rownames(top.600), ";"), "[[", 1)
top.600 <- cbind(top.600, uniprot)

uniprot <- sapply(strsplit(rownames(top.900), ";"), "[[", 1)
top.900 <- cbind(top.900, uniprot)

uniprot <- sapply(strsplit(rownames(top.1800), ";"), "[[", 1)
top.1800 <- cbind(top.1800, uniprot)



#### significant table of different norm methods
top.10.sign <- top.10[top.10[, "adj.P.Val"] <0.05,]
top.30.sign <- top.30[top.30[, "adj.P.Val"] <0.05,]
top.60.sign <- top.60[top.60[, "adj.P.Val"] <0.05,]
top.300.sign <- top.300[top.300[, "adj.P.Val"] <0.05,]
top.600.sign <-  top.600[top.600[, "adj.P.Val"] <0.05,]
top.900.sign <- top.900[top.900[, "adj.P.Val"] <0.05,]
top.1800.sign <-  top.1800[top.1800[, "adj.P.Val"] <0.05,]


#### statistics
### number of significant dergulated targets
DE2.RUV <- c(length(rownames(top.10.sign)),length(rownames(top.30.sign)),length(rownames(top.60.sign)),
             length(rownames(top.300.sign)),length(rownames(top.600.sign)),length(rownames(top.900.sign)),
             length(rownames(top.1800.sign)))
DE2.RUV




#### investigate specific proteins
# 
# top.10[grep(pattern = "*PLNM*", x = rownames(top.10)),]
# top.10[grep(pattern = "VASP", x = rownames(top.10)),]
# top.600["P50552;VASP;S239;KVSKQEEASGGPTAPK;9984",]
# top.600.vasp <- top.600[top.900$uniprot=="P50552",]
# top.10.vasp <- top.10[top.900$uniprot=="P50552",]
# top.1800.vasp <- top.1800[top.1800$uniprot=="P50552",]



##### prepare export for excel and mysql

top.10$id <- rownames(top.10)
top.30$id <- rownames(top.30)
top.60$id <- rownames(top.60)
top.300$id <- rownames(top.300)
top.600$id <- rownames(top.600)
top.900$id <- rownames(top.900)
top.1800$id <- rownames(top.1800)


top.10$symbol <- sapply(strsplit(rownames(top.10), ";"), "[[", 2)
top.30$symbol <- sapply(strsplit(rownames(top.30), ";"), "[[", 2)
top.60$symbol <- sapply(strsplit(rownames(top.60), ";"), "[[", 2)
top.300$symbol <- sapply(strsplit(rownames(top.300), ";"), "[[", 2)
top.600$symbol <- sapply(strsplit(rownames(top.600), ";"), "[[", 2)
top.900$symbol <- sapply(strsplit(rownames(top.900), ";"), "[[", 2)
top.1800$symbol <- sapply(strsplit(rownames(top.1800), ";"), "[[", 2)


new_10 <- top.10[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.10[c("uniprot","logFC","P.Value","adj.P.Val")])), ]
new_30 <- top.30[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.30[c("uniprot","logFC","P.Value","adj.P.Val")])), ]
new_60 <- top.60[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.60[c("uniprot","logFC","P.Value","adj.P.Val")])), ]
new_300 <- top.300[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.300[c("uniprot","logFC","P.Value","adj.P.Val")])), ]
new_600 <- top.600[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.600[c("uniprot","logFC","P.Value","adj.P.Val")])), ]
new_900 <- top.900[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.900[c("uniprot","logFC","P.Value","adj.P.Val")])), ]
new_1800 <- top.1800[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.1800[c("uniprot","logFC","P.Value","adj.P.Val")])), ]


top.all <- as.data.frame(merge(new_900[,2:4], new_1800[,2:4], by="row.names"))
rownames(top.all) <- top.all$Row.names
top.all <- top.all[,-1]


top.all <- merge(new_600[,2:4], top.all, by="row.names", all = TRUE)
colnames(top.all) <- c("Row.names", "logFC.600","P.Value.600","adj.P.Val.600","logFC.900","P.Value.900","adj.P.Val.900", "logFC.1800","P.Value.1800","adj.P.Val.1800")
rownames(top.all) <- top.all$Row.names
top.all <- top.all[,-1]


top.all <- merge(new_300[,2:4], top.all, by="row.names", all = TRUE)
colnames(top.all) <- c("Row.names", "logFC.300","P.Value.300","adj.P.Val.300", "logFC.600","P.Value.600","adj.P.Val.600","logFC.900","P.Value.900","adj.P.Val.900","logFC.1800","P.Value.1800","adj.P.Val.1800")
rownames(top.all) <- top.all$Row.names
top.all <- top.all[,-1]


top.all <- merge(new_60[,2:4], top.all, by="row.names", all = TRUE)
colnames(top.all) <- c("Row.names", "logFC.60","P.Value.60","adj.P.Val.60", "logFC.300","P.Value.300","adj.P.Val.300", "logFC.600","P.Value.600","adj.P.Val.600","logFC.900","P.Value.900","adj.P.Val.900", "logFC.1800","P.Value.1800","adj.P.Val.1800")
rownames(top.all) <- top.all$Row.names
top.all <- top.all[,-1]


top.all <- merge(new_30[,2:4], top.all, by="row.names", all = TRUE)
colnames(top.all) <- c("Row.names", "logFC.30","P.Value.30","adj.P.Val.30", "logFC.60","P.Value.60","adj.P.Val.60", "logFC.300","P.Value.300","adj.P.Val.300", "logFC.600","P.Value.600","adj.P.Val.600","logFC.900","P.Value.900","adj.P.Val.900", "logFC.1800","P.Value.1800","adj.P.Val.1800")
rownames(top.all) <- top.all$Row.names
top.all <- top.all[,-1]


top.all <- merge(new_10[,c(1,5,2,3,4)], top.all, by="row.names", all = TRUE)
rownames(top.all) <- top.all$Row.names
top.all <- top.all[,-1]
colnames(top.all) <- c("uniprot","symbol", "logFC.10","P.Value.10","adj.P.Val.10", "logFC.30","P.Value.30","adj.P.Val.30", "logFC.60","P.Value.60","adj.P.Val.60", "logFC.300","P.Value.300","adj.P.Val.300", "logFC.600","P.Value.600","adj.P.Val.600","logFC.900","P.Value.900","adj.P.Val.900", "logFC.1800","P.Value.1800","adj.P.Val.1800")



## read tables locally or
## RS: connect to mysql server and wrtie the table from the mysql database
#write.table(top.all, "../top.all.txt", sep="\t", , row.names=TRUE)
write.table(top.all, "data/processed_data/top.all.txt", sep="\t", , row.names=TRUE)




#### transfer mysql database
#### export main tables 

## read tables locally or
## RS: connect to mysql server and wrtie the table from the mysql database
write.table(top.10, "data/processed_data/top.10.txt", sep="\t", , row.names=FALSE)
write.table(top.30, "data/processed_data/top.30.txt", sep="\t", , row.names=FALSE)
write.table(top.60, "data/processed_data/top.60.txt", sep="\t", , row.names=FALSE)
write.table(top.300, "data/processed_data/top.300.txt", sep="\t", , row.names=FALSE)
write.table(top.600, "data/processed_data/top.600.txt", sep="\t", , row.names=FALSE)
write.table(top.900, "data/processed_data/top.900.txt", sep="\t", , row.names=FALSE)
write.table(top.1800, "data/processed_data/top.1800.txt", sep="\t", , row.names=FALSE)

