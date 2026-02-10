####################################################################################
##### Diagnosing and Correcting for batch effect
## by SPSs from phosR study
## by Genes with low variance
## by Genes with low Rank in differential regulation in this study  --> this i do
## by generating SPSs empirically
## by a combination of all

## empirical stable phosphosites
## perform deferential regulation analysis and 
## determine the peptides that are least deregulated
## define new expression set and look at quality 

setwd("F:/Masterarbeit_BioWis/Phosphoproteomics/Phosphoproteom validation/scripts/")

library("PhosR")
library(matrixStats)


##### define input and visualize

x_tt <- as.factor(grps)
x_tt2 <- as.factor(grps2)
x_tt3 <- as.factor(grps3)
x_tt4 <- as.factor(grps4)


dataset <- 2^SummarizedExperiment::assay(ppe_imputed_scaled,"scaled")
sum(is.na(dataset))
#dataset <- na.omit(dataset)

test4 <- SummarizedExperiment::assay(ppe_imputed_scaled,"scaled")
dev.new()
plotQC(test4, grps=grps4, 
       labels = sapply(strsplit(colnames(ppe0), "_"), "[[",1), panel="pca")


#### Prenormalize make a betweenLane Normalization
dataset_norm <- betweenLaneNormalization(dataset, which="upper")

colors <- brewer.pal(8, "Paired")

plotRLE(dataset_norm, outline=FALSE, ylim=c(-1, 1), col=colors[x_tt])
plotPCA(dataset_norm, col=colors[x_tt], cex=1, ylim=c(-0.3, 0.3), xlim=c(-0.3, 0.3))

plotQC(log2(dataset_norm), grps=grps4, 
       labels = colnames(dataset_norm), panel = "pca") +
  ggplot2::ggtitle("back transformed")

plotQC(log2(dataset_norm), grps=grps, 
       labels = colnames(dataset_norm), panel = "pca") +
  ggplot2::ggtitle("back transformed")




make_all_contrasts <- function (group, delim="_vs_", design_matrix){
  
  suppressMessages(require(limma))
  
  #/ ensure that group levels are unique
  group <- sort(unique(as.character(group)))
  
  #/ make all combinations
  cb   <- combn(group, 2, FUN = function(x){paste0(x[1], "-", x[2])})
  
  
  #/ make contrasts
  contrasts<- limma::makeContrasts(contrasts=cb, levels=design_matrix)
  colnames(contrasts) <- gsub("-", delim, colnames(contrasts))
  
  return(contrasts)
}



design_new <- model.matrix(~0+x_tt)
colnames(design_new) <- gsub("x_tt", "", colnames(design_new))
v_new <- voom(dataset_norm, design_new)
#v_new <- log2(dataset_norm)
fit_new <- lmFit(v_new, design_new)

contrast.matrix_new <- make_all_contrasts(x_tt, delim= "_vs_", design_new)

fit2 <- contrasts.fit(fit_new, contrast.matrix_new)
fit2 <- eBayes(fit2, trend=TRUE)


top1 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=1, sort.by = "none")
top2 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=2, sort.by = "none")
top3 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=3, sort.by = "none")
top4 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=4, sort.by = "none")
top5 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=5, sort.by = "none")
top6 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=6, sort.by = "none")
top7 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=7, sort.by = "none")
top8 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=8, sort.by = "none")
top9 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=9, sort.by = "none")
top10 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=10, sort.by = "none")
top11 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=11, sort.by = "none")
top12 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=12, sort.by = "none")
top13 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=13, sort.by = "none")
top14 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=14, sort.by = "none")
top15 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=15, sort.by = "none")
top16 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=16, sort.by = "none")
top17 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=17, sort.by = "none")
top18 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=18, sort.by = "none")
top19 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=19, sort.by = "none")
top20 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=20, sort.by = "none")
top21 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=21, sort.by = "none")


p.value <- cbind(top1[,c(5)],top2[,c(5)],top3[,c(5)],top4[,c(5)],top5[,c(5)],top6[,c(5)],top7[,c(5)],top8[,c(5)],top9[,c(5)],top10[,c(5)]
,top11[,c(5)],top12[,c(5)],top13[,c(5)],top14[,c(5)],top15[,c(5)],top16[,c(5)],top17[,c(5)],top18[,c(5)],top19[,c(5)],top20[,c(5)],top21[,c(5)])
rownames(p.value) <- rownames(fit2)

logfc <- cbind(top1[,c(1)],top2[,c(1)],top3[,c(1)],top4[,c(1)],top5[,c(1)],top6[,c(1)],top7[,c(1)],top8[,c(1)],top9[,c(1)],top10[,c(1)]
,top11[,c(1)],top12[,c(1)],top13[,c(1)],top14[,c(1)],top15[,c(1)],top16[,c(1)],top17[,c(1)],top18[,c(1)],top19[,c(1)],top20[,c(1)],top21[,c(1)])
rownames(logfc) <- rownames(fit2)

basemean <- cbind(top1[,c(2)],top2[,c(2)],top3[,c(2)],top4[,c(2)],top5[,c(2)],top6[,c(2)],top7[,c(2)],top8[,c(2)],top9[,c(2)],top10[,c(2)]
,top11[,c(2)],top12[,c(2)],top13[,c(2)],top14[,c(2)],top15[,c(2)],top16[,c(2)],top17[,c(2)],top18[,c(2)],top19[,c(2)],top20[,c(2)],top21[,c(2)])
rownames(basemean) <- rownames(fit2)


#p.value <- cbind(top7[,c(5)],top16[,c(5)],top21[,c(5)])
#rownames(p.value) <- rownames(fit2)

#logfc <- cbind(top7[,c(1)],top16[,c(1)],top21[,c(1)])
#rownames(logfc) <- rownames(fit2)

#basemean <- cbind(top7[,c(2)],top16[,c(2)],top21[,c(2)])
#rownames(basemean) <- rownames(fit2)



p_all <- as.data.frame(rowMedians(p.value))
colnames(p_all) <- c("P")
rownames(p_all) <- rownames(p.value)

logfc_all <- as.data.frame(rowMedians(logfc))
colnames(logfc_all) <- c("logfc")
rownames(logfc_all) <- rownames(logfc)

basemean_all <- as.data.frame(rowMedians(basemean))
colnames(basemean_all) <- c("basemean")
rownames(basemean_all) <- rownames(basemean)

p.value.count <- as.data.frame(sapply(1:nrow(p.value), function(i) sum(p.value[i,] < 0.05)))
colnames(p.value.count) <- c("P")
rownames(p.value.count) <- rownames(p.value)



empirical_p <- rownames(subset(p_all, P >= 0.9))
length(empirical_p)

empirical_logfc <- rownames(subset(logfc_all, abs(logfc) <= 0.1))
length(empirical_logfc)

#quantile(basemean)
empirical_basemean <- rownames(subset(basemean_all, basemean >= as.numeric(quantile(basemean, probs = seq(0, 1, 0.05))['25%'])))
length(empirical_basemean)

empirical_p_count <- rownames(subset(p.value.count, P < 1))
length(empirical_p_count)

p_rank <- as.data.frame(colMedians(colRanks(p.value)))
rownames(p_rank) <- rownames(p.value)
colnames(p_rank) <- "V1"
empirical_p_rank <- rownames(subset(p_rank, V1 >3000))
length(empirical_p_rank)


#Reduce(intersect, list(a,b,c))
empirical_topall <- Reduce(intersect, list(empirical_p,empirical_p_count,empirical_p_rank,empirical_logfc))
#empirical_topall <- empirical_p
length(empirical_topall)


## write tables locally or
## RS: connect to mysql server and wrtie the table from the mysql database
write.table(empirical_topall, "../data/processed_data/empirical_topall.txt", sep = "\t")


#use empirical_control_proteins3.R

empirical_top <- Reduce(intersect, list(empirical_p_rank,empirical_basemean))
#empirical_topall <- empirical_p
length(empirical_top)

## write tables locally or
## RS: connect to mysql server and wrtie the table from the mysql database
write.table(empirical_top, "../data/processed_data/empirical_top.txt", sep = "\t")



##### RUVphospho (RUVIII) normalization
design = model.matrix(~ grps - 1)
ctl2 = which(rownames(ppe_imputed_scaled) %in% empirical_topall)
ctl2 = which(rownames(ppe_imputed_scaled) %in% new_data_empirical)

#### RUVphospho normalization
## RUVphospho shows a good tendency for statistics, but the intensity value are curiously shifted 
## but still 
ppe_norm = RUVphospho(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), M = design, k = 16, ctl = ctl2)

ppe <- PhosphoExperiment(assays = list(normalised = as.matrix(ppe_norm)), 
                          UniprotID = sapply(strsplit(rownames(ppe_imputed_scaled), ";"), "[[", 1),
                          Site = as.numeric(gsub("[A-Z]","", sapply(strsplit(rownames(ppe_imputed_scaled), ";"), "[[", 3))),
                          GeneSymbol = sapply(strsplit(rownames(ppe_imputed_scaled), ";"), "[[", 2),
                          Residue = gsub("[0-9]","", sapply(strsplit(rownames(ppe_imputed_scaled), ";"), "[[", 3)),
                          Sequence = sapply(strsplit(rownames(ppe_imputed_scaled), ";"), "[[", 4))




## plot clusteirng PCA
test2 <- SummarizedExperiment::assay(ppe_imputed_scaled,"scaled")
sum(is.na(test2))
test4 <- SummarizedExperiment::assay(ppe,"normalised")
sum(is.na(test4))
dev.new()
plotQC(test4, grps=grps, 
       labels = sapply(strsplit(colnames(ppe), "_"), "[[",1), panel="pca")

###DIFFERENT NORMALIZATION METHODS: we dont do this part anymore#####
#### RUVg normalization 
abundance <- 2^SummarizedExperiment::assay(ppe_imputed_scaled,"scaled")
abundance2 = apply(abundance, 1:2, function(x) round(x))
set <- newSeqExpressionSet(as.matrix(abundance2),phenoData = data.frame(x_tt, row.names=colnames(abundance2)))

set2 <- RUVg(set, empirical_topall, k=14)
SummarizedExperiment::assay(ppe,"normalised") <- log2(set2@assayData[["normalizedCounts"]])



## plot clusteirng PCA
test4 <- log2(set2@assayData[["normalizedCounts"]])
dev.new()
plotQC(test4, grps=grps, 
       labels = sapply(strsplit(colnames(ppe0), "_"), "[[",1), panel="pca")


## write tables locally or
## RS: connect to mysql server and wrtie the table from the mysql database
write.table(SummarizedExperiment::assay(ppe,"normalised"), "../data/normalized_intensities.txt", sep = "\t")



#### RUVg visualization of the influence of the value k on normalization 
plot_list = list()
for (j in c(10, 600, 1800)) {
  for (i in c(1, seq(5, 65, 2), 67)) {
    set2 <- RUVg(set, empirical_topall, k=i) #how many dimensions (replicates) to delete to get better clustering?
    SummarizedExperiment::assay(ppe,"normalised") <- log2(set2@assayData[["normalizedCounts"]])
    test4 <-SummarizedExperiment::assay(ppe,"normalised")
    p1 <- plotQC(test4, grps=grps, 
                 labels = colnames(ppe), panel="pca")
    p2<-plotQC(SummarizedExperiment::assay(ppe, "normalised"), grps=grps, 
               labels = colnames(ppe), panel="dendrogram")+
      ggtitle("batch corrected with empirical top genes")
    plot_list[[i]] <- ggpubr::ggarrange(p1, p2, nrow = 2)
    ggsave(plot_list[[i]], file=paste0("../analysis/ruv_and_k/", "plot_k_", i,".tiff"), width = 44.45, height = 27.78, units = "cm", dpi=300)
    #plot_list<-c(plot_list, figure)
    #assign(paste('figure_k', i, sep="_"), figure)
    dataset <- SummarizedExperiment::assay(ppe, "normalised")
    design <- model.matrix(~0 + grps)
    fit0 <- lmFit(dataset, design)
    
    #print(j)
    contra<- paste("grpsx", j, "sek_CXCR7-grpsx", j, "sek_DMSO", sep = "")
	print(contra)
    #10sek
    contrast.matrix <- makeContrasts(contra, levels=design)
    fit1 <- contrasts.fit(fit0, contrast.matrix)
    fit1 <- eBayes(fit1, trend=TRUE)
    top <- topTable(fit1, number=nrow(fit1), adjust.method="fdr", sort.by="p", p.value=1)
    assign(paste('top',j, '_k', i, sep="_"), top)
    a <- paste(j, "sec(k=", i, "):",
               dim(top[top$adj.P.Val < 0.05, ])[1], sep = "")
    print(a)
  }
}


#### combat
library(bladderbatch)
data(bladderdata)
dat <- bladderEset[1:50,]

pheno = pData(dat)
edata = exprs(dat)
batch = pheno$batch
mod = model.matrix(~as.factor(cancer), data=pheno)

# parametric adjustment
combat_edata1 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

# non-parametric adjustment, mean-only version
combat_edata2 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=FALSE, mean.only=TRUE)

# reference-batch version, with covariates
combat_edata3 = ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE, ref.batch=3)


batch2  <- as.numeric(sapply(strsplit(colnames(abundance2), "_"), "[[", 3))

pheno2 <- data.frame(sample = 1:length(colnames(abundance2)),
						time = sapply(strsplit(colnames(abundance2), "_"), "[[", 1),
						treatment = sapply(strsplit(colnames(abundance2), "_"), "[[", 2),
						batch = sapply(strsplit(colnames(abundance2), "_"), "[[", 3)
						)
rownames(pheno2) <- colnames(abundance2)
#colnames(pheno) <- c("sample", "time", "treatmnt", "batch")


abndance3 <- log2(abundance)
combat_edata3 = ComBat(dat=as.matrix(log2(abundance)), batch=batch2, mod=NULL, par.prior=TRUE, prior.plots=FALSE)



####################################################################################
#### statistics of normalization and quality control
### get different datasets to compare to
### make firgure for publication

test1 <- SummarizedExperiment::assay(ppe0,"Quantification")
test2 <- SummarizedExperiment::assay(ppe_imputed_omit,"imputed")
test3 <- SummarizedExperiment::assay(ppe_imputed_scaled,"scaled")
test4 <- SummarizedExperiment::assay(ppe,"normalised")
test5 <- combat_edata3

### define the groups and replicates
grps = gsub("_[0-9][0-9]", "", colnames(ppe0))
grps

grps2 <- sapply(strsplit(colnames(ppe0), "_"), "[[",2)
grps2

grps3 <- sapply(strsplit(colnames(ppe0), "_"), "[[", 1)
grps3

grps4 <- sapply(strsplit(colnames(ppe0), "_"), "[[", 3)
grps4


x_tt <- as.factor(grps)
x_tt2 <- as.factor(grps2)
x_tt3 <- as.factor(grps3)
x_tt4 <- as.factor(grps4)

##### plot distributions
dev.new()
dd <- melt(test3, variable.name = "sample")
ggplot(dd, aes(value, colour = Var2)) + geom_density(bw = "SJ") 
#+ xlim(0,1e+7)

dev.new()
dd <- melt(test4, variable.name = "sample")
ggplot(dd, aes(value, colour = Var2)) + geom_density(bw = "SJ") 
#+ xlim(0,1e+7)


###### plot RLE and PCA
colors <- brewer.pal(8, "Paired")
dev.new()
plotRLE(test4, outline=FALSE, col=colors[x_tt], cex.axis=0.8, las = 2)
dev.new()
plotPCA(test4, col=colors[x_tt], cex=1, ylim=c(-0.3, 0.3), xlim=c(-0.3, 0.3))



####### plot clusteirng PCA, dendogram
dev.new()
plotQC(test3, grps=grps4, 
       labels = sapply(strsplit(colnames(ppe0), "_"), "[[",1), panel="pca")
dev.new()
plotQC(test4, grps=grps, 
       labels = sapply(strsplit(colnames(ppe0), "_"), "[[",1), panel="pca")




dev.new()
p1 = plotQC(SummarizedExperiment::assay(ppe0, "Quantification"), grps=grps, 
            labels = colnames(ppe0), panel = "dendrogram" )+
  ggtitle("log2 transformed")
p2 = plotQC(SummarizedExperiment::assay(ppe, "scaled"), grps=grps, 
            labels = colnames(ppe), panel="dendrogram")+
  ggtitle("median scaled")
p3 = plotQC(SummarizedExperiment::assay(ppe, "normalised"), grps=grps, 
            labels = colnames(ppe), panel="dendrogram")+
  ggtitle("batch corrected with empirical top genes")
ggpubr::ggarrange(p1, p2, p3, nrow = 3)


dev.new()
# plot after batch correction #NA ERROR
p1 = plotQC(SummarizedExperiment::assay(ppe, "scaled"), panel = "pca", 
            grps=grps3, labels = colnames(ppe)) +
  ggplot2::ggtitle("Before Batch correction")
p2 = plotQC(SummarizedExperiment::assay(ppe, "normalised"), panel="pca", 
            grps=grps, labels = colnames(ppe)) +
  ggplot2::ggtitle("After Batch correction")
ggpubr::ggarrange(p1, p2, nrow = 2)

