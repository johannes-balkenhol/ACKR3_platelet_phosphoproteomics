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
## build a vector accoridng to the sample
#vect1 <- as.data.frame(table(grps))$Freq

setwd("F:/Masterarbeit_BioWis/Phosphoproteomics/Phosphoproteom validation/scripts/")



library("PhosR")


## hierachical cluster
plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), panel = "dendrogram", grps=grps, labels = colnames(ppe_imputed_scaled)) +
  ggplot2::ggtitle("Before batch correction")

## PCA: grps redefined for visualization
plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), grps=grps, 
       labels = colnames(ppe_imputed_scaled), panel = "pca") +
  ggplot2::ggtitle("Before batch correction")


x_tt <- as.factor(grps)
x_tt2 <- as.factor(grps2)
x_tt3 <- as.factor(grps3)


dataset <- 2^SummarizedExperiment::assay(ppe_imputed_scaled,"scaled")
#sum(is.na(dataset))
#dataset <- na.omit(dataset)


colors <- brewer.pal(8, "Paired")
plotRLE(dataset, outline=FALSE, ylim=c(-1, 1), col=colors[x_tt], cex.axis=0.8, las = 2)
plotPCA(dataset, col=colors[x_tt], cex=1, ylim=c(-0.3, 0.3), xlim=c(-0.3, 0.3))

plotQC(dataset, grps=grps3, 
       labels = colnames(dataset), panel = "pca") +
  ggplot2::ggtitle("back transformed")

plotQC(dataset, grps=grps2, 
       labels = colnames(dataset), panel = "pca") +
  ggplot2::ggtitle("back transformed")

#### Prenormalize make a betweenLane Normalization
dataset_norm <- betweenLaneNormalization(dataset, which="upper")


plotRLE(dataset_norm, outline=FALSE, ylim=c(-1, 1), col=colors[x_tt])
plotPCA(dataset_norm, col=colors[x_tt], cex=1, ylim=c(-0.3, 0.3), xlim=c(-0.3, 0.3))

plotQC(dataset_norm, grps=grps, 
       labels = colnames(dataset_norm), panel = "pca") +
  ggplot2::ggtitle("back transformed")

plotQC(dataset_norm, grps=grps2, 
       labels = colnames(dataset_norm), panel = "pca") +
  ggplot2::ggtitle("back transformed")

#### use Differential Expression Analysis for determining empirical control genes
### differential analysis all vs. all  to find empirical control genes
design_new <- model.matrix(~0+x_tt)
colnames(design_new) <- gsub("x_tt", "", colnames(design_new))
#v_new <- voom(dataset_norm, design_new)
v_new <- log2(dataset_norm)
fit_new <- lmFit(v_new, design_new)



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

contrast.matrix_new <- make_all_contrasts(x_tt, delim= "_vs_", design_new)

fit2_new <- contrasts.fit(fit_new, contrast.matrix_new)
fit2_new <- eBayes(fit2_new, trend=TRUE)


top_all_newmethod <- topTable(fit2_new, number=nrow(fit2_new), adjust.method="BH", sort.by="F", p.value=1)


empirical_topall <- rownames(dataset_norm)[which((!rownames(dataset_norm) %in% rownames(top_all_newmethod)[1:4085]))]
length(empirical_topall)

## write tables locally or
## RS: connect to mysql server and wrtie the table from the mysql database
write.table(empirical_topall, "../data/empirical_topall.txt", sep = "\t")



##### RUVphospho (RUVIII) normalization
design = model.matrix(~ grps - 1)
ctl2 = which(rownames(ppe_imputed_scaled) %in% empirical_topall)


#### RUVphospho normalization
## RUVphospho shows a good tendency for statistics, but the intensity value are curiously shifted 
## but still 
ppe = RUVphospho(ppe_imputed_scaled, M = design, k = 14, ctl = ctl2)



#### RUVg normalization
abundance <- 2^SummarizedExperiment::assay(ppe_imputed_scaled,"scaled")
abundance2 = apply(abundance, 1:2, function(x) round(x))
set <- newSeqExpressionSet(as.matrix(abundance2),phenoData = data.frame(x_tt, row.names=colnames(abundance2)))

set2 <- RUVg(set, empirical_topall, k=19)
SummarizedExperiment::assay(ppe,"normalised") <- log2(set2@assayData[["normalizedCounts"]])


## write tables locally or
## RS: connect to mysql server and wrtie the table from the mysql database
write.table(SummarizedExperiment::assay(ppe_imputed_scaled,"normalised"), "../data/normalized_intensities.txt", sep = "\t")



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