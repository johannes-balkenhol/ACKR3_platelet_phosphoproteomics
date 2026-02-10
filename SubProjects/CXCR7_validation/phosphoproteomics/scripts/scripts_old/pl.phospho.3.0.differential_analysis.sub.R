####################################################################################
##### differential analysis
### define the groups and replicates
grps = gsub("_[0-9][0-9]", "", colnames(ppe0))
grps


## define header
donor_nr <- gsub("Donor", "",sapply(strsplit(colnames(ppe0), "_"), "[[", 3))
time_point <- sapply(strsplit(colnames(ppe0), "_"), "[[", 1)
#time_point <- paste("x", time_point, sep = "")
timepoint_fac <- gsub("x", "",time_point)
timepoint_fac <- as.numeric(gsub("sek", "",timepoint_fac))
condition <- sapply(strsplit(colnames(ppe0), "_"), "[[", 2)

#### fitting
#dataset <- SummarizedExperiment::assay(ppe_imputed_tmp, "Quantification")
dataset <- SummarizedExperiment::assay(ppe, "normalised")
# donor_nr, condition, time_point, timpoint_fac
#design <- model.matrix(~0 + grps + timepoint_fac:condition)
#design <- model.matrix(~0 + grps + condition)
design <- model.matrix(~ grps - 1)
colnames(design) <- make.names(colnames(design))
#design <- subset(design, select = -(timepoint_fac.conditionCtrl))

#design_interaction <- model.matrix(~0+time+treatment)
#f <- gsub("_\\d", "", colnames(ppe))
#X <- model.matrix(~ f - 1)
fit <- lmFit(dataset, design)


## timecourse comparison example
#times <- rep(c(0, 1, 3, 5, 7), 2)
#geno <- rep(c("WT", "KO"), each=5)
#design2 <- model.matrix(~0 + geno + times:geno)
#colnames(design) <- make.names(colnames(design)) # make names syntactically vald
#con <- makeContrasts(genoKO.times - genoWT.times, levels=design)


#treatment <- sapply(strsplit(colnames(ppe), "_"), "[[",2)
#time <- sapply(strsplit(colnames(ppe), "_"), "[[",1)


#(treat10sec - ctrl0sec) - (ctrl10sec - ctrl0sec)
#######
## tt10
contrast.matrix <- makeContrasts((grpsx10sek_CXCR7 - grpsx10sek_DMSO), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.10 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

##  tt600
contrast.matrix <- makeContrasts((grpsx600sek_CXCR7 - grpsx600sek_DMSO), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.600 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

## tt1800
contrast.matrix <- makeContrasts(grpsx1800sek_CXCR7 - grpsx1800sek_DMSO, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.1800 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

#######
## tt10.dmso vs tt0 ctrl
contrast.matrix <- makeContrasts((grpsx10sek_DMSO - grpsx0sek_Ctrl), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.10.dmso.vs.0s <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

##  tt600.dmso vs tt0 ctrl
contrast.matrix <- makeContrasts((grpsx600sek_DMSO - grpsx0sek_Ctrl), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.600.dmso.vs.0s <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

## tt1800.dmso vs tt0 ctrl
contrast.matrix <- makeContrasts(grpsx1800sek_DMSO - grpsx0sek_Ctrl, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.1800.dmso.vs.0s <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

#######
## tt10.cxcr7 vs tt0 ctrl
contrast.matrix <- makeContrasts((grpsx10sek_CXCR7 - grpsx0sek_Ctrl), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.10.cxcr7.vs.0s <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

##  tt600 vs tt0 ctrl
contrast.matrix <- makeContrasts((grpsx600sek_CXCR7 - grpsx0sek_Ctrl), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.600.cxcr7.vs.0s <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)


## tt1800 vs tt0 ctrl
contrast.matrix <- makeContrasts(grpsx1800sek_CXCR7 - grpsx0sek_Ctrl, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.1800.cxcr7.vs.0s <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

## timecourse comparison
#contrast.matrix <- makeContrasts(timepoint_fac.conditionCXCR7 - timepoint_fac.conditionDMSO, levels=design)
#fit2 <- contrasts.fit(fit, contrast.matrix)
#fit2 <- eBayes(fit2, trend=TRUE)
#top.timecourse <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)


#contrast_matrix <- makeContrasts(
#  contrast = (grpsx10sek_CXCR7 - grpsx10sek_DMSO) + (grpsx600sek_CXCR7 - grpsx600sek_DMSO) +
#    (grpsx1800sek_CXCR7 - grpsx1800sek_DMSO),
#  levels = colnames(design)
#)


### annotation
top.10$uniprot <- sapply(strsplit(rownames(top.10), ";"), "[[", 1)
top.600$uniprot <- sapply(strsplit(rownames(top.600), ";"), "[[", 1)
top.1800$uniprot <- sapply(strsplit(rownames(top.1800), ";"), "[[", 1)
top.10.dmso.vs.0s$uniprot <- sapply(strsplit(rownames(top.10.dmso.vs.0s), ";"), "[[", 1)
top.600.dmso.vs.0s$uniprot <- sapply(strsplit(rownames(top.600.dmso.vs.0s), ";"), "[[", 1)
top.1800.dmso.vs.0s$uniprot <- sapply(strsplit(rownames(top.1800.dmso.vs.0s), ";"), "[[", 1)
top.10.cxcr7.vs.0s$uniprot <- sapply(strsplit(rownames(top.10.cxcr7.vs.0s), ";"), "[[", 1)
top.600.cxcr7.vs.0s$uniprot <- sapply(strsplit(rownames(top.600.cxcr7.vs.0s), ";"), "[[", 1)
top.1800.cxcr7.vs.0s$uniprot <- sapply(strsplit(rownames(top.1800.cxcr7.vs.0s), ";"), "[[", 1)

top.10$symbol <- sapply(strsplit(rownames(top.10), ";"), "[[", 2)
top.600$symbol <- sapply(strsplit(rownames(top.600), ";"), "[[", 2)
top.1800$symbol <- sapply(strsplit(rownames(top.1800), ";"), "[[", 2)
top.10.dmso.vs.0s$symbol <- sapply(strsplit(rownames(top.10.dmso.vs.0s), ";"), "[[", 2)
top.600.dmso.vs.0s$symbol <- sapply(strsplit(rownames(top.600.dmso.vs.0s), ";"), "[[", 2)
top.1800.dmso.vs.0s$symbol <- sapply(strsplit(rownames(top.1800.dmso.vs.0s), ";"), "[[", 2)
top.10.cxcr7.vs.0s$symbol <- sapply(strsplit(rownames(top.10.cxcr7.vs.0s), ";"), "[[", 2)
top.600.cxcr7.vs.0s$symbol <- sapply(strsplit(rownames(top.600.cxcr7.vs.0s), ";"), "[[", 2)
top.1800.cxcr7.vs.0s$symbol <- sapply(strsplit(rownames(top.1800.cxcr7.vs.0s), ";"), "[[", 2)

top.10$psite <- sapply(strsplit(rownames(top.10), ";"), "[[", 3)
top.600$psite <- sapply(strsplit(rownames(top.600), ";"), "[[", 3)
top.1800$psite <- sapply(strsplit(rownames(top.1800), ";"), "[[", 3)
top.10.dmso.vs.0s$psite <- sapply(strsplit(rownames(top.10.dmso.vs.0s), ";"), "[[", 3)
top.600.dmso.vs.0s$psite <- sapply(strsplit(rownames(top.600.dmso.vs.0s), ";"), "[[", 3)
top.1800.dmso.vs.0s$psite <- sapply(strsplit(rownames(top.1800.dmso.vs.0s), ";"), "[[", 3)
top.10.cxcr7.vs.0s$psite <- sapply(strsplit(rownames(top.10.cxcr7.vs.0s), ";"), "[[", 3)
top.600.cxcr7.vs.0s$psite <- sapply(strsplit(rownames(top.600.cxcr7.vs.0s), ";"), "[[", 3)
top.1800.cxcr7.vs.0s$psite <- sapply(strsplit(rownames(top.1800.cxcr7.vs.0s), ";"), "[[", 3)


top.10$id <- rownames(top.10)
top.600$id <- rownames(top.600)
top.1800$id <- rownames(top.1800)
top.10.dmso.vs.0s$id <- rownames(top.10.dmso.vs.0s)
top.600.dmso.vs.0s$id <- rownames(top.600.dmso.vs.0s)
top.1800.dmso.vs.0s$id <- rownames(top.1800.dmso.vs.0s)
top.10.cxcr7.vs.0s$id <- rownames(top.10.cxcr7.vs.0s)
top.600.cxcr7.vs.0s$id <- rownames(top.600.cxcr7.vs.0s)
top.1800.cxcr7.vs.0s$id <- rownames(top.1800.cxcr7.vs.0s)



#### significant table
top.10.sign <- top.10[top.10[, "adj.P.Val"] <0.05,]
top.600.sign <-  top.600[top.600[, "adj.P.Val"] <0.05,]
top.1800.sign <-  top.1800[top.1800[, "adj.P.Val"] <0.05,]
top.10.dmso.vs.0s.sign <- top.10.dmso.vs.0s[top.10.dmso.vs.0s[, "adj.P.Val"] <0.05,]
top.600.dmso.vs.0s.sign <-  top.600.dmso.vs.0s[top.600.dmso.vs.0s[, "adj.P.Val"] <0.05,]
top.1800.dmso.vs.0s.sign <-  top.1800.dmso.vs.0s[top.1800.dmso.vs.0s[, "adj.P.Val"] <0.05,]
top.10.cxcr7.vs.0s.sign <- top.10.cxcr7.vs.0s[top.10.cxcr7.vs.0s[, "adj.P.Val"] <0.05,]
top.600.cxcr7.vs.0s.sign <-  top.600.cxcr7.vs.0s[top.600.cxcr7.vs.0s[, "adj.P.Val"] <0.05,]
top.1800.cxcr7.vs.0s.sign <-  top.1800.cxcr7.vs.0s[top.1800.cxcr7.vs.0s[, "adj.P.Val"] <0.05,]


#### statistics
### number of significant dergulated targets
DE2.RUV <- c(length(rownames(top.10.sign)),length(rownames(top.600.sign)),length(rownames(top.1800.sign)),
length(rownames(top.10.dmso.vs.0s.sign)),length(rownames(top.600.dmso.vs.0s.sign)),length(rownames(top.1800.dmso.vs.0s.sign)),
length(rownames(top.10.cxcr7.vs.0s.sign)),length(rownames(top.600.cxcr7.vs.0s.sign)),length(rownames(top.1800.cxcr7.vs.0s.sign))
)
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



##### prepare export for excel and mysql

new_10 <- top.10[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.10[c("uniprot","logFC","P.Value","adj.P.Val")])), ]
new_600 <- top.600[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.600[c("uniprot","logFC","P.Value","adj.P.Val")])), ]
new_1800 <- top.1800[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.1800[c("uniprot","logFC","P.Value","adj.P.Val")])), ]
new_10.dmso.vs.0s <- top.10.dmso.vs.0s[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.10.dmso.vs.0s[c("uniprot","logFC","P.Value","adj.P.Val")])), ]
new_600.dmso.vs.0s <- top.600.dmso.vs.0s[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.600.dmso.vs.0s[c("uniprot","logFC","P.Value","adj.P.Val")])), ]
new_1800.dmso.vs.0s <- top.1800.dmso.vs.0s[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.1800.dmso.vs.0s[c("uniprot","logFC","P.Value","adj.P.Val")])), ]
new_10.cxcr7.vs.0s <- top.10.cxcr7.vs.0s[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.10.cxcr7.vs.0s[c("uniprot","logFC","P.Value","adj.P.Val")])), ]
new_600.cxcr7.vs.0s <- top.600.cxcr7.vs.0s[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.600.cxcr7.vs.0s[c("uniprot","logFC","P.Value","adj.P.Val")])), ]
new_1800.cxcr7.vs.0s <- top.1800.cxcr7.vs.0s[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.1800.cxcr7.vs.0s[c("uniprot","logFC","P.Value","adj.P.Val")])), ]




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
write.table(top.all, "data/processed_data/top.all.txt", sep="\t", , row.names=TRUE)

top.all.dmso.vs.0s <- as.data.frame(merge(new_600.dmso.vs.0s[,2:4], new_1800.dmso.vs.0s[,2:4], by="row.names"))
rownames(top.all.dmso.vs.0s) <- top.all.dmso.vs.0s$Row.names
top.all.dmso.vs.0s <- top.all.dmso.vs.0s[,-1]
top.all.dmso.vs.0s <- merge(new_10.dmso.vs.0s[,c(1,5,2,3,4)], top.all.dmso.vs.0s, by="row.names", all = TRUE)
rownames(top.all.dmso.vs.0s) <- top.all.dmso.vs.0s$Row.names
top.all.dmso.vs.0s <- top.all.dmso.vs.0s[,-1]
colnames(top.all.dmso.vs.0s) <- c("uniprot","symbol","logFC.10","P.Value.10","adj.P.Val.10","logFC.600","P.Value.600","adj.P.Val.600","logFC.1800","P.Value.1800","adj.P.Val.1800")
## read tables locally or
## RS: connect to mysql server and wrtie the table from the mysql database
#write.table(top.all, "../top.all.txt", sep="\t", , row.names=TRUE)
write.table(top.all.dmso.vs.0s, "data/processed_data/top.all.dmso.vs.0s.txt", sep="\t", , row.names=TRUE)


top.all.cxcr7.vs.0s <- as.data.frame(merge(new_600.cxcr7.vs.0s[,2:4], new_1800.cxcr7.vs.0s[,2:4], by="row.names"))
rownames(top.all.cxcr7.vs.0s) <- top.all.cxcr7.vs.0s$Row.names
top.all.cxcr7.vs.0s <- top.all.cxcr7.vs.0s[,-1]
top.all.cxcr7.vs.0s <- merge(new_10.cxcr7.vs.0s[,c(1,5,2,3,4)], top.all.cxcr7.vs.0s, by="row.names", all = TRUE)
rownames(top.all.cxcr7.vs.0s) <- top.all.cxcr7.vs.0s$Row.names
top.all.cxcr7.vs.0s <- top.all.cxcr7.vs.0s[,-1]
colnames(top.all.cxcr7.vs.0s) <- c("uniprot","symbol","logFC.10","P.Value.10","adj.P.Val.10","logFC.600","P.Value.600","adj.P.Val.600","logFC.1800","P.Value.1800","adj.P.Val.1800")
## read tables locally or
## RS: connect to mysql server and wrtie the table from the mysql database
#write.table(top.all, "../top.all.txt", sep="\t", , row.names=TRUE)
write.table(top.all.cxcr7.vs.0s, "data/processed_data/top.all.cxcr7.vs.0s.txt", sep="\t", , row.names=TRUE)


#### transfer mysql database
#### export main tables 

## read tables locally or
## RS: connect to mysql server and wrtie the table from the mysql database
write.table(top.10, "data/top.10.txt", sep="\t", , row.names=FALSE)
write.table(top.600, "data/top.600.txt", sep="\t", , row.names=FALSE)
write.table(top.1800, "data/top.1800.txt", sep="\t", , row.names=FALSE)

write.table(top.10.dmso.vs.0s, "data/top.10.dmso.vs.0s.txt", sep="\t", , row.names=FALSE)
write.table(top.600.dmso.vs.0s, "data/top.600.dmso.vs.0s.txt", sep="\t", , row.names=FALSE)
write.table(top.1800.dmso.vs.0s, "data/top.1800.dmso.vs.0s.txt", sep="\t", , row.names=FALSE)

write.table(top.10.cxcr7.vs.0s, "data/top.10.cxcr7.vs.0s.txt", sep="\t", , row.names=FALSE)
write.table(top.600.cxcr7.vs.0s, "data/top.600.cxcr7.vs.0s.txt", sep="\t", , row.names=FALSE)
write.table(top.1800.cxcr7.vs.0s, "data/top.1800.cxcr7.vs.0s.txt", sep="\t", , row.names=FALSE)


##### check some candidiates
top.1800[grep("P25098", rownames(top.1800)),]

which(rownames(top.600) 

grep("prkd2", rownames(top.600) )

top.600[grep("PRK", rownames(top.600)),]