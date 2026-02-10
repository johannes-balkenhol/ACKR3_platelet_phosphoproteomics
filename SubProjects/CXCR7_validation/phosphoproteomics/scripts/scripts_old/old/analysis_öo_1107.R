########################################################
########################################################
########################################################
####### phosphoproteom analysis
# load data
# filter data
# determine empirical control genes
# reduction of unwanted variance
# determine log2foldchanges of the samples
# downstream analysis
#

###install packages####


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("RUVSeq")
BiocManager::install("RUV") #NOT AVAILABLE, DO we need this?
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("org.Rn.eg.db")
BiocManager::install("EDASeq")
BiocManager::install("PhosR")
install.packages("remotes")

remotes::install_github("biobenkj/ATACseeker") #DO WE NEED THIS?
BiocManager::install("directPA")
BiocManager::install("reactome.db")
BiocManager::install("MSnbase")

###load packages####



suppressPackageStartupMessages({
  library("BiocManager")
  library("RColorBrewer")
  library("DESeq2")
  library("ruv")
  library("limma")
  library("dplyr")
  library(PhosR)
  library("edgeR")
  library(ggplot2)
  library("reshape2")
  library(calibrate)
  library(limma)
  library(directPA)
  library(annotate)
  library(PhosR)
  library(stringr)
  library(ggplot2)
  library("RUVSeq")
  library("remotes")
  library(org.Rn.eg.db)
  library("MSnbase")
  library("EDASeq")
  library(psych) #for describe function
})



#library("ATACseeker")


### load data #### 
#We dont use peptide info for anything?

########## before the data are delivered by ISAS (TRR240_A08_01_Phosphoproteome.exe)
########## in excel (A08_pp_info.xlxs) relevant columns are selected and the file is seprated in
########## A08_phosR.txt and A08_pp_info.txt
########## A08_phosR.txt is prprocessed in perl (parserPhosR.sql) to cope with duplicated Uniport_IDs
########## the resulting table A08_phosR_v2.txt is loaded to mysql (phosR_gene_symbol_mapping) to get the corresponding protein symbol 
########## and the correspondign protein name and description
########## an idea is build by uniport_id;smbol;phosphosite;sequence;peptide_id
########## and that the table A08_phosR_v3.txt is generated for R analysis

#F:\Masterarbeit_BioWis\Proteomics\Irina Pleines\data
setwd("C:/Users/Ozge/Nextcloud/Phosphoproteom/phosphoproteom validation") #folder structure should fit.

raw_abundance <- read.table("data/A08_val_phosphoR_v3.txt", sep="\t", header=TRUE, dec=",")
#peptide_info <- read.table("A08_pp_info.txt", sep="\t", header=TRUE, dec=",")
#pka_target <- read.table("pka_target.csv", sep="\t", header=TRUE, dec=",")


###check & filter data########################################

raw_abundance2 <- raw_abundance[,c(-1)]
rownames(raw_abundance2) <- raw_abundance[,1]

raw_abundance2 <- as.data.frame(raw_abundance2)


donor_nr <- gsub("Donor", "",sapply(strsplit(colnames(raw_abundance2), "_"), "[[", 2))
time_point <- sapply(strsplit(colnames(raw_abundance2), "_"), "[[", 3)
time_point <- paste("x", time_point, sep = "")
condition <- sapply(strsplit(colnames(raw_abundance2), "_"), "[[", 4)

colnames(raw_abundance2) <- paste(time_point, condition, donor_nr, sep = "_")


raw_abundance3 <- raw_abundance2 %>% dplyr::select(sort(names(raw_abundance2))) #sorts the columns
raw_abundance3_log <- log2(raw_abundance3) 





test2 <- as.matrix(raw_abundance3)
test3 <- as.matrix(raw_abundance3_log)

#check distribution ####

# ## plot distributions
dd <- melt(test2, variable.name = "sample")
dd_log <- melt(test3, variable.name = "sample")



ggplot(dd, aes(value, colour = Var2)) + geom_density(bw = "SJ") +
  theme(legend.position="none")
#+ xlim(0,1e+7)

ggplot(dd_log, aes(value, colour = Var2)) + geom_density(bw = "SJ") +
  theme(legend.position="none")
#+ xlim(0,1e+7)


desc_data<- describe(raw_abundance3)
c(min(desc_data$skew), max(desc_data$skew))  #skewed to right, skewness all positive [15.74092 31.60951]
c(min(desc_data$kurtosis), max(desc_data$kurtosis)) #heavy tails [292.8265 1293.2373]

desc_logdata <- describe(raw_abundance3_log)
c(min(desc_logdata$skew), max(desc_logdata$skew))  #almost no skewness [0.5383390 0.7913213]
c(min(desc_logdata$kurtosis), max(desc_logdata$kurtosis)) #almost no kurtosis [0.5881861 0.8762273]


#####imputation#



########## imputation as there are still NA and RUV does not work with NA
########## use PhoR 
########## prepare table in mysql


### Setting up the PhosphoExperiment object ####


ppe0 <- PhosphoExperiment(assays = list(Quantification = as.matrix(raw_abundance3_log)))

ppe0 <- PhosphoExperiment(assays = list(Quantification = as.matrix(raw_abundance3_log)), 
                          UniprotID = sapply(strsplit(rownames(ppe0), ";"), "[[", 1),
                          Site = as.numeric(gsub("[A-Z]","", sapply(strsplit(rownames(ppe0), ";"), "[[", 3))),
                          GeneSymbol = sapply(strsplit(rownames(ppe0), ";"), "[[", 2),
                          Residue = gsub("[0-9]","", sapply(strsplit(rownames(ppe0), ";"), "[[", 3)),
                          Sequence = sapply(strsplit(rownames(ppe0), ";"), "[[", 4))


grps = gsub("_[0-9][0-9]", "", colnames(ppe0))
grps



grps2 <- sapply(strsplit(colnames(ppe0), "_"), "[[",3)
grps2

grps3 <- sapply(strsplit(colnames(ppe0), "_"), "[[", 1)
grps3


###Part A. Preprocessing####
####A.1. Filtering of phosphosites ####


#filtering of phosphosites so that only phosphosites with quantification for
#at least 50% of the replicates in at least eight of the conditions are retained.

#how do we decide n here? we have one condition but 8 timepoints  --> here filtering already
#removes half, and the imputations do not actually impute anything?

ppe_filtered <- selectGrps(ppe0, grps, 0.3, n=7) #same method or 0.3 

#1266 phosphosites have no NAs.

dim(ppe_filtered)
dim(ppe0)


####A.2. Imputation of phosphosites ####

#####A.2.1. Site- and condition-specific imputation

set.seed(123)
ppe_imputed_tmp <- scImpute(ppe_filtered, 0.5, grps)[,colnames(ppe_filtered)]



test0 <- SummarizedExperiment::assay(ppe0,"Quantification")
test <- SummarizedExperiment::assay(ppe_filtered,"Quantification")
test1 <- SummarizedExperiment::assay(ppe_imputed_tmp,"Quantification")
test2 <- SummarizedExperiment::assay(ppe_imputed_tmp,"imputed")

sum(is.na(test0))
sum(is.na(test))
sum(is.na(test1))
sum(is.na(test2))



#####A.2.2. Paired tail-based imputation (second imputing is optional) #DONT####
 
 #vect1 <- c(7, 6, 7, 7, 7, 6, 7, 6)

set.seed(123)
ppe_imputed <- ppe_imputed_tmp

SummarizedExperiment::assay(ppe_imputed[,seq(1,10)],"imputed") <- PhosR::ptImpute(SummarizedExperiment::assay(ppe_imputed[,seq(11,70)],"imputed"),
                                                                                   SummarizedExperiment::assay(ppe_imputed[,seq(1,10)],"imputed"),
                                                                                   percent1 = 1, percent2 = 0.3, paired = FALSE)

SummarizedExperiment::assay(ppe_imputed[,seq(11,20)],"imputed") <- PhosR::ptImpute(SummarizedExperiment::assay(ppe_imputed[,seq(21,30)],"imputed"),
                                                                                 SummarizedExperiment::assay(ppe_imputed[,seq(11,20)],"imputed"),
                                                                                 percent1 = 1, percent2 = 0.3, paired = FALSE)


SummarizedExperiment::assay(ppe_imputed[,seq(31,40)],"imputed") <- PhosR::ptImpute(SummarizedExperiment::assay(ppe_imputed[,seq(41,50)],"imputed"),
                                                                                  SummarizedExperiment::assay(ppe_imputed[,seq(31,40)],"imputed"),
                                                                                  percent1 = 1, percent2 = 0.3, paired = FALSE)

SummarizedExperiment::assay(ppe_imputed[,seq(51,60)],"imputed") <- PhosR::ptImpute(SummarizedExperiment::assay(ppe_imputed[,seq(61,70)],"imputed"),
                                                                                  SummarizedExperiment::assay(ppe_imputed[,seq(51,60)],"imputed"),
                                                                                  percent1 = 1, percent2 = 0.3, paired = FALSE)



####what is this for?
#ppe_imputed <- ppe_imputed_tmp


######### remove NA
abundance_no.na <- na.omit(SummarizedExperiment::assay(ppe_imputed, "imputed"))

ppe_imputed_omit <- PhosphoExperiment(assays = list(imputed = as.matrix(abundance_no.na)))

UniprotID(ppe_imputed_omit)  <- sapply(strsplit(rownames(ppe_imputed_omit), ";"), "[[", 1)
GeneSymbol(ppe_imputed_omit) <- sapply(strsplit(rownames(ppe_imputed_omit), ";"), "[[", 2)
Residue(ppe_imputed_omit) <- gsub("[0-9]","", sapply(strsplit(rownames(ppe_imputed_omit), ";"), "[[", 3))
Site(ppe_imputed_omit) <- gsub("[A-Z]","", sapply(strsplit(rownames(ppe_imputed_omit), ";"), "[[", 3))
Sequence(ppe_imputed_omit) <- sapply(strsplit(rownames(ppe_imputed_omit), ";"), "[[", 4)

ppe_imputed <- ppe_imputed_omit


#######################scaling imputation#######################################
######### scaling imputation
#ppe_imputed_scaled <- medianScaling(ppe_imputed_omit, scale = FALSE, assay = "imputed")
ppe_imputed_scaled <- medianScaling(ppe_imputed, scale = FALSE, assay = "imputed")

test3 <- SummarizedExperiment::assay(ppe_imputed,"imputed")
test3 <- ppe_imputed_scaled@assays@data@listData[["imputed"]]
test4 <- SummarizedExperiment::assay(ppe_imputed_scaled,"scaled")


sum(is.na(test))
sum(is.na(test2))
sum(is.na(test3))
sum(is.na(test4))
dim(test4)





#PCA # grps redefined for visualization
plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), grps=grps, 
       labels = colnames(ppe_imputed_scaled), panel = "pca") +
  ggplot2::ggtitle("Before batch correction")



########################quantification before and after impuation######################################
######### quantification before and after impuation
p0 = plotQC(SummarizedExperiment::assay(ppe0,"Quantification"), 
            labels=colnames(ppe0), 
            panel = "quantify", grps = grps)
p1 = plotQC(SummarizedExperiment::assay(ppe_filtered,"Quantification"), 
            labels=colnames(ppe_filtered), 
            panel = "quantify", grps = grps)
p2 = plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), 
            labels=colnames(ppe_imputed_scaled), panel = "quantify", grps = grps)
ggpubr::ggarrange(p0, p1, p2, nrow = 3)

#hierachical clustering before and after impuation

p0 = plotQC(SummarizedExperiment::assay(ppe0,"Quantification"), 
            labels=colnames(ppe_filtered), panel = "dendrogram", 
            grps = grps2)

p1 = plotQC(SummarizedExperiment::assay(ppe_filtered,"Quantification"), 
            labels=colnames(ppe_filtered), panel = "dendrogram", 
            grps = grps2)
p2 = plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), 
            labels=colnames(ppe_imputed_scaled), 
            panel = "dendrogram", grps = grps2)

ggpubr::ggarrange(p0, p1, p2, nrow = 3)

p1 = plotQC(SummarizedExperiment::assay(ppe_filtered,"Quantification"), 
            labels=colnames(ppe_filtered), panel = "dendrogram", 
            grps = grps2)
p2 = plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), 
            labels=colnames(ppe_imputed_scaled), 
            panel = "dendrogram", grps = grps2)
ggpubr::ggarrange(p1, p2, nrow = 1)


######### diagnosis for batch correction####


######## Diagnosing batch effect
#hierachical cluster
plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), panel = "dendrogram", grps=grps, labels = colnames(ppe_imputed_scaled)) +
  ggplot2::ggtitle("Before batch correction")

#PCA # grps redefined for visualization
plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), grps=grps, 
       labels = colnames(ppe_imputed_scaled), panel = "pca") +
  ggplot2::ggtitle("Before batch correction")




######## Correcting batch effect####
### by SPSs from phosR study
### by Genes with low variance
### by Genes with low Rank in differential regulation in this study  --> this i do
### by generating SPSs empirically
### by a combination of all


##############################################################
########## empirical stable phosphosites
########## perform deferential regulation analysis and 
########## determine the peptides that are least deregulated
### define new expression set and look at quality 
##### build a vector accoridng to the sample
#vect1 <- as.data.frame(table(grps))$Freq
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



########## make a betweenLaneNorm 
dataset_norm <- betweenLaneNormalization(dataset, which="upper")


plotRLE(dataset_norm, outline=FALSE, ylim=c(-1, 1), col=colors[x_tt])
plotPCA(dataset_norm, col=colors[x_tt], cex=1, ylim=c(-0.3, 0.3), xlim=c(-0.3, 0.3))

plotQC(dataset_norm, grps=grps, 
       labels = colnames(dataset_norm), panel = "pca") +
  ggplot2::ggtitle("back transformed")

plotQC(dataset_norm, grps=grps2, 
       labels = colnames(dataset_norm), panel = "pca") +
  ggplot2::ggtitle("back transformed")

########## use Differential Expression Analysis for determining empirical control genes ####
#### fitting

design_new <- model.matrix(~0+x_tt)
colnames(design_new) <- gsub("x_tt", "", colnames(design_new))
v_new <- voom(dataset_norm, design_new)
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

#test_fit <- topTable(fit2_new, number=nrow(fit2_new), coef = 1, adjust.method="BH", sort.by = "M", p.value=1)

#write.table(dataset_norm, "dataset_norm.txt", sep = "\t")

top_all_newmethod <- topTable(fit2_new, number=nrow(fit2_new), adjust.method="BH", sort.by="F", p.value=1)

#write.table(top_all_newmethod, "top_all_new.txt", sep = "\t")
empirical_topall <- rownames(dataset_norm)[which((!rownames(dataset_norm) %in% rownames(top_all_newmethod)[1:3100]))]
#ncol(contrast.matrix_new)


#a <- dataset_norm[!rownames(dataset_norm) %in% rownames(top_all_newmethod)[1:3100],]


#pairwise empirical#####

bottom_list = c()

x=matrix(data=NA, nrow=nrow(fit2_new), ncol=ncol(contrast.matrix_new))


for (i in 1:ncol(contrast.matrix_new)) {
    df<- topTable(fit2_new, coef = i, number=nrow(fit2_new),  
                  adjust.method="BH", sort.by="M", p.value=1)
    df2 <- df[df$adj.P.Val <= 0.05,]
    new_df <- df[ order(row.names(df)), ]
    x[,i]= new_df$adj.P.Val
    bottom <- rownames(dataset_norm)[which(!(rownames(dataset_norm) %in% rownames(df2)[1:3100]))]
    assign(paste('top',colnames(contrast.matrix_new)[i],sep='_'),df)
    assign(paste('newdf', i, sep="_"), new_df)
    #assign(paste('pvals', i, sep="_"), pvals)
    assign(paste('unsig',colnames(contrast.matrix_new)[i],sep='_'),df2)
    assign(paste('unsig',i,sep='_'),bottom)
    bottom_list <- c(bottom_list, paste('unsig', i,sep='_' ))
    
}

x <- rowMedians(x)
x<- as.data.frame(x)
rownames(x) <- rownames(top_x0sek_Ctrl_vs_x10sek_CXCR7[order(row.names(top_x0sek_Ctrl_vs_x10sek_CXCR7)), ])
#rownames(x) <- rownames(top_01_vs_02[order(row.names(top_01_vs_02)), ])

colnames(x) <- "medianp"
x$name <- rownames(x)


x_filt<-x[x$medianp > 0.1, ]
empirical_x_filt <- rownames(x_filt)

emp_intersect_filt <- intersect(empirical_topall, empirical_x_filt)


#idea: if we're seeing donors cluster together, can't we find the genes that do not change between the donors
#and use those to normalize?

empirical_pairwise_new <- Reduce(intersect, list(unsig_1, unsig_2, unsig_3, unsig_4, unsig_5, unsig_6, unsig_7, unsig_8, unsig_9, 
                       unsig_10, unsig_11, unsig_12, unsig_13, unsig_14, unsig_15, unsig_16, unsig_17, unsig_18,
                       unsig_19, unsig_20, unsig_21, unsig_22, unsig_23, unsig_24, unsig_25, unsig_26, unsig_27, unsig_28, 
                       unsig_29, unsig_30, unsig_31, unsig_32, unsig_33, unsig_34, unsig_35, unsig_36, 
                       unsig_37, unsig_38, unsig_39, unsig_40, unsig_41, unsig_42, unsig_43, unsig_44, unsig_45))


empirical_pair_all <- intersect(empirical_topall, empirical_pairwise_new)
empirical_x_pair <- intersect(empirical_x_filt, empirical_pairwise_new)
emp_x_pair_all <- Reduce(intersect, list(empirical_topall, empirical_pairwise_new, empirical_x_filt))
# old_diff <- setdiff(empirical, empirical_new)
# new_diff <- setdiff(empirical_new, empirical)
# common <- intersect(empirical, empirical_new)
# old_diff_norm <- dataset_norm[rownames(dataset_norm) %in% old_diff, ]
# new_diff_norm <- dataset_norm[rownames(dataset_norm) %in% new_diff, ]
# common_norm <- dataset_norm[rownames(dataset_norm) %in% common, ]
# 
# 
# heatmap(old_diff_norm)
# heatmap(new_diff_norm)
# heatmap(common_norm)

#empirical_new_uniprot<- sapply(strsplit(empirical_new, ";"), "[[", 2)

#length(empirical_new)
#134
#length(unique(empirical_new_uniprot))
#113

#write.table(empirical_new, "empirical_new.txt", sep="\t", , row.names=FALSE)
#write.table(empirical_new_uniprot, "empirical_uniprot_new.txt", sep="\t", , row.names=FALSE)





######### prepare RUVIII norm####
 design = model.matrix(~ grps - 1) #why do we define model matrix like this? it's the same (removes intercept)
 #design

 
# #ctl0 = which(sites2 %in% SPSs)
# ctl0 = which(sites %in% SPSs)
# #ctl0 = which(sites %in% inhouse_SPSs)
# specified_control <- rownames(ppe_imputed_scaled[ctl0,])

 
empiricalcontrol <-  which(rownames(ppe_imputed_scaled) %in% dual)
#length(empiricalcontrol)


ctlvar = which(rownames(ppe_imputed_scaled) %in% ControlGenes3)
# #used 27.09.2021
ctl2 = which(rownames(ppe_imputed_scaled) %in% empirical_topall)
# ctl3 = which(rownames(ppe_imputed_scaled) %in% empiricalcontrol)
# ctl4 <- c(ctl,ctl2)
# empirical2 <- c(empirical,specified_control)


### RUVphospho shows a good tendency for statistics, but the intensity value are curiously shifted 
#and with increasing k the number of diff. regulated sites increase continuously
# thats why we here create norm data but dont use those, but replace those by RUVg
ppe = RUVphospho(ppe_imputed_scaled, M = design, k = 7, ctl = ctl2)

### therefore we need to translate to RUVg usage
abundance <- 2^SummarizedExperiment::assay(ppe_imputed_scaled,"scaled")
abundance2 = apply(abundance, 1:2, function(x) round(x))
set <- newSeqExpressionSet(as.matrix(abundance2),phenoData = data.frame(x_tt, row.names=colnames(abundance2)))

#set2_mk <- RUVg(set0_mk, Control_mk, k=7)
set2 <- RUVg(set, empirical_topall, k=15) #how many dimensions (replicates) to delete to get better clustering?
SummarizedExperiment::assay(ppe,"normalised") <- log2(set2@assayData[["normalizedCounts"]])
###
##deciding on k####

plot_list = list()
for (j in c(10, 600, 1800)) {
  for (i in c(1, seq(5, 65, 5), 67)) {
    set2 <- RUVg(set, empirical_topall, k=i) #how many dimensions (replicates) to delete to get better clustering?
    SummarizedExperiment::assay(ppe,"normalised") <- log2(set2@assayData[["normalizedCounts"]])
    test4 <-SummarizedExperiment::assay(ppe,"normalised")
    p1 <- plotQC(test4, grps=grps, 
                 labels = colnames(ppe), panel="pca")
    p2<-plotQC(SummarizedExperiment::assay(ppe, "normalised"), grps=grps, 
               labels = colnames(ppe), panel="dendrogram")+
      ggtitle("batch corrected with empirical top genes")
    plot_list[[i]] <- ggpubr::ggarrange(p1, p2, nrow = 2)
    ggsave(plot_list[[i]], file=paste0("plot_k_", i,".tiff"), width = 44.45, height = 27.78, units = "cm", dpi=300)
    #plot_list<-c(plot_list, figure)
    #assign(paste('figure_k', i, sep="_"), figure)
    dataset <- SummarizedExperiment::assay(ppe, "normalised")
    design <- model.matrix(~0 + grps)
    fit0 <- lmFit(dataset, design)
    
    #print(j)
    contra<- paste("grpsx", j, "sek_CXCR7-grpsx", j, "sek_DMSO", sep = "")
    #10sek
    contrast.matrix <- makeContrasts(contra, levels=design)
    fit1 <- contrasts.fit(fit0, contrast.matrix)
    fit1 <- eBayes(fit1, trend=TRUE)
    top <- topTable(fit1, number=nrow(fit1), adjust.method="fdr", sort.by="p", p.value=1)
    assign(paste('top',j, '_k', i, sep="_"), top)
    a <- paste(j, "sec(k=", i, "):",
               dim(top[top$adj.P.Val < 0.05, ])[1], sep = "")
    print(a)
  }}

#####

#
test2 <- SummarizedExperiment::assay(ppe0,"Quantification")
test3 <-SummarizedExperiment::assay(ppe,"scaled")
test4 <-SummarizedExperiment::assay(ppe,"normalised")


# #iterative second emp####
# 
# dataset <- 2^SummarizedExperiment::assay(ppe,"normalised")
# dataset_norm <- betweenLaneNormalization(dataset, which="upper")
# 
# design_new <- model.matrix(~0+x_tt)
# colnames(design_new) <- gsub("x_tt", "", colnames(design_new))
# v_new <- voom(dataset_norm, design_new)
# fit_new <- lmFit(v_new, design_new)
# 
# contrast.matrix_new <- make_all_contrasts(x_tt, delim= "_vs_", design_new)
# 
# fit2_new <- contrasts.fit(fit_new, contrast.matrix_new)
# fit2_new <- eBayes(fit2_new, trend=TRUE)
# 
# 
# 
# top_all_newmethod <- topTable(fit2_new, number=nrow(fit2_new), adjust.method="BH", sort.by="F", p.value=1)
# 
# empirical_topall4 <- rownames(dataset_norm)[which((!rownames(dataset_norm) %in% rownames(top_all_newmethod)[1:3100]))]
# 
# empirical_iter4 <- intersect(empirical_iter3, empirical_topall4)
# 
# 
# ctl2 = which(rownames(ppe_imputed_scaled) %in% empirical_topall3)
# 
# ppe = RUVphospho(ppe_imputed_scaled, M = design, k = 7, ctl = ctl2)
# 
# 
# abundance <- 2^SummarizedExperiment::assay(ppe_imputed_scaled,"scaled")
# abundance2 = apply(abundance, 1:2, function(x) round(x))
# set <- newSeqExpressionSet(as.matrix(abundance2),phenoData = data.frame(x_tt, row.names=colnames(abundance2)))
# 
# set2 <- RUVg(set, empirical_topall3, k=7) #how many dimensions (replicates) to delete to get better clustering?
# SummarizedExperiment::assay(ppe,"normalised") <- log2(set2@assayData[["normalizedCounts"]])
# ###
# ##
# #
# test2 <- SummarizedExperiment::assay(ppe,"Quantification")
# test3 <-SummarizedExperiment::assay(ppe,"scaled")
# test4 <-SummarizedExperiment::assay(ppe,"normalised")
# 
# 
######

## plot distributions
dev.new()
dd <- melt(test2, variable.name = "sample")
ggplot(dd, aes(value, colour = Var2)) + geom_density(bw = "SJ") 
#+ xlim(0,1e+7)

dev.new()
dd <- melt(test4, variable.name = "sample")
ggplot(dd, aes(value, colour = Var2)) + geom_density(bw = "SJ") 
#+ xlim(0,1e+7)


dev.new()
plotQC(test4, grps=grps, 
       labels = colnames(ppe), panel="pca")

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






###############################################################
##############################################################
##############################################################
###############################################################
######### differential regulationanalysis v2 (source Balkenhol)

#### fitting
#dataset <- SummarizedExperiment::assay(ppe_imputed_tmp, "Quantification")
dataset <- SummarizedExperiment::assay(ppe, "normalised")
design <- model.matrix(~0 + grps)
#f <- gsub("_\\d", "", colnames(ppe))
#X <- model.matrix(~ f - 1)
fit <- lmFit(dataset, design)


treatment <- sapply(strsplit(colnames(ppe0), "_"), "[[",2)
time <- sapply(strsplit(colnames(ppe0), "_"), "[[",1)

datesetgrps3 <- 
design_interaction <- model.matrix(~0+time+treatment)

#

(treat10sec - ctrl0sec) - (ctrl10sec - ctrl0sec)

##### tt00 vs. tt10
contrast.matrix <- makeContrasts(((grpsx10sek_CXCR7-grpsx0sek_Ctrl)-(grpsx10sek_DMSO-grpsx0sek_Ctrl)), levels=design)
contrast.matrix <- makeContrasts((grpsx10sek_CXCR7-grpsx10sek_DMSO), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.10_2 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

##### tt00 vs. tt1800
contrast.matrix <- makeContrasts(grpsx1800sek_CXCR7-grpsx1800sek_DMSO, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.1800 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

##### tt00 vs. tt600
contrast.matrix <- makeContrasts((grpsx600sek_CXCR7-grpsx600sek_DMSO), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.600 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

##### tt00 vs. tt300
contrast.matrix <- makeContrasts(grpsx600sek_CXCR7-grpsx0sek_Ctrl, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)

top.600toctrl <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)


##### tt00 vs. tt600
contrast.matrix <- makeContrasts(x_tttt600-x_tttt00, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.600 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)


##### tt00 vs. tt900
contrast.matrix <- makeContrasts(x_tttt900-x_tttt00, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.900 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)


##### tt00 vs. tt1800
contrast.matrix <- makeContrasts(x_tttt1800-x_tttt00, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.1800 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)



##############################################################
##############################################################
######### map PKA targets
######### and filter top deregulated #WE DONT NEED

UniprotID<- sapply(strsplit(rownames(ppe), ";"), "[[", 1)
GeneSymbol <- sapply(strsplit(rownames(ppe), ";"), "[[", 2)

pka_target <- read.table("pka_target.csv", sep="\t", header=TRUE, dec=",")
pka <- as.list(pka_target)

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


top.10.pka <- top.10[top.10$uniprot %in% pka_target$uniprot,] 
top.30.pka <- top.30[top.30$uniprot %in% pka_target$uniprot,]
top.60.pka <- top.60[top.60$uniprot %in% pka_target$uniprot,]
top.300.pka <- top.300[top.300$uniprot %in% pka_target$uniprot,]
top.600.pka <- top.600[top.600$uniprot %in% pka_target$uniprot,]
top.900.pka <- top.900[top.900$uniprot %in% pka_target$uniprot,]
top.1800.pka <- top.1800[top.1800$uniprot %in% pka_target$uniprot,]


top.10.pka.sign <- top.10.pka[top.10.pka[, "adj.P.Val"] <0.05,]
top.30.pka.sign <- top.30.pka[top.30.pka[, "adj.P.Val"] <0.05,]
top.60.pka.sign <- top.60.pka[top.60.pka[, "adj.P.Val"] <0.05,]
top.300.pka.sign <- top.300.pka[top.300.pka[, "adj.P.Val"] <0.05,]
top.600.pka.sign <- top.600.pka[top.600.pka[, "adj.P.Val"] <0.05,]
top.900.pka.sign <- top.900.pka[top.900.pka[, "adj.P.Val"] <0.05,]
top.1800.pka.sign <- top.1800.pka[top.1800.pka[, "adj.P.Val"] <0.05,]


top.10.sign <- top.10[top.10[, "adj.P.Val"] <0.05,]
top.30.sign <- top.30[top.30[, "adj.P.Val"] <0.05,]
top.60.sign <- top.60[top.60[, "adj.P.Val"] <0.05,]
top.300.sign <- top.300[top.300[, "adj.P.Val"] <0.05,]
top.600.sign <-  top.600[top.600[, "adj.P.Val"] <0.05,]
top.900.sign <-  top.900[top.900[, "adj.P.Val"] <0.05,]
top.1800.sign <-  top.1800[top.1800[, "adj.P.Val"] <0.05,]


##############################################################
##############################################################
######### statistics

######### number of significant dergulated targets
DE2.RUV <- c(length(rownames(top.10.pka.sign)),length(rownames(top.30.pka.sign)),length(rownames(top.60.pka.sign)),length(rownames(top.300.pka.sign)),length(rownames(top.600.pka.sign)),length(rownames(top.900.pka.sign)),length(rownames(top.1800.pka.sign)))
DE2.RUV

DE2.RUV <- c(length(rownames(top.10.sign)),length(rownames(top.30.sign)),length(rownames(top.60.sign)),length(rownames(top.300.sign)),length(rownames(top.600.sign)),length(rownames(top.900.sign)),length(rownames(top.1800.sign)))
DE2.RUV

top.900["P50552;VASP;S239;KVSKQEEASGGPTAPK;9984",]
top.300.vasp <- top.900[top.900$uniprot=="P50552",]
top.900.vasp <- top.900[top.900$uniprot=="P50552",]
top.1800.vasp <- top.1800[top.1800$uniprot=="P50552",]

######### PKA Target Analysis
UniprotID<- sapply(strsplit(rownames(ppe), ";"), "[[", 1)
GeneSymbol <- sapply(strsplit(rownames(ppe), ";"), "[[", 2)

pka <- as.list(pka_target)

pka_id <- which(UniprotID(ppe) %in% pka[["uniprot_id"]])
length(which(UniprotID(ppe) %in% pka[["uniprot_id"]]))


######### analyze PKA target phosphosites
######### filtered by Masterprotein

# PLot for Discussion (Quantification vs, Normalised) PCA
p1 = plotQC(SummarizedExperiment::assay(ppe[pka_id,], "Quantification"), panel = "pca", 
            grps=grps, labels = colnames(ppe[pka_id,])) +
  ggplot2::ggtitle("Before Batch correction")
p2 = plotQC(SummarizedExperiment::assay(ppe[pka_id,], "normalised"), grps=grps, 
            labels = colnames(ppe[pka_id,]), panel="pca") +
  ggplot2::ggtitle("After Batch correction")
ggpubr::ggarrange(p1, p2, nrow = 2)

# PLot for Discussion (Quantification vs, Normalised) Dendrogram
p1 = plotQC(SummarizedExperiment::assay(ppe[pka_id,], "Quantification"), grps=grps, 
            labels = colnames(ppe[pka_id,]), panel = "dendrogram" )
p2 = plotQC(SummarizedExperiment::assay(ppe[pka_id,], "normalised"), grps=grps, 
            labels = colnames(ppe[pka_id,]), panel="dendrogram")
ggpubr::ggarrange(p1, p2, nrow = 1,  widths = 1, heights = 1)



################################################################
################################################################
############ prepare tables for
############ export for excel and mysql


################################################################
############ export PKA targets merged
top.10.pka$id <- rownames(top.10.pka)
top.30.pka$id <- rownames(top.30.pka)
top.60.pka$id <- rownames(top.60.pka)
top.300.pka$id <- rownames(top.300.pka)
top.600.pka$id <- rownames(top.600.pka)
top.900.pka$id <- rownames(top.900.pka)
top.1800.pka$id <- rownames(top.1800.pka)


top.10.pka$symbol <- sapply(strsplit(rownames(top.10.pka), ";"), "[[", 2)
top.30.pka$symbol <- sapply(strsplit(rownames(top.30.pka), ";"), "[[", 2)
top.60.pka$symbol <- sapply(strsplit(rownames(top.60.pka), ";"), "[[", 2)
top.300.pka$symbol <- sapply(strsplit(rownames(top.300.pka), ";"), "[[", 2)
top.600.pka$symbol <- sapply(strsplit(rownames(top.600.pka), ";"), "[[", 2)
top.900.pka$symbol <- sapply(strsplit(rownames(top.900.pka), ";"), "[[", 2)
top.1800.pka$symbol <- sapply(strsplit(rownames(top.1800.pka), ";"), "[[", 2)


new_10 <- top.10.pka[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.10.pka[c("uniprot","logFC","P.Value","adj.P.Val")])), ]
new_30 <- top.30.pka[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.30.pka[c("uniprot","logFC","P.Value","adj.P.Val")])), ]
new_60 <- top.60.pka[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.60.pka[c("uniprot","logFC","P.Value","adj.P.Val")])), ]
new_300 <- top.300.pka[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.300.pka[c("uniprot","logFC","P.Value","adj.P.Val")])), ]
new_600 <- top.600.pka[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.600.pka[c("uniprot","logFC","P.Value","adj.P.Val")])), ]
new_900 <- top.900.pka[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.900.pka[c("uniprot","logFC","P.Value","adj.P.Val")])), ]
new_1800 <- top.1800.pka[c("uniprot","logFC","P.Value","adj.P.Val","symbol")][ order(row.names(top.1800.pka[c("uniprot","logFC","P.Value","adj.P.Val")])), ]

#
top.pka <- cbind(rownames(new_10),new_10, new_30, new_60, new_300, new_600, new_900, new_1800)

write.table(top.pka, "top.pka.txt", sep="\t", , row.names=FALSE)


################################################################
############ export all targets

#### export PKA targets merged
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

#
top <- cbind(rownames(new_10),new_10, new_30, new_60, new_300, new_600, new_900, new_1800)

write.table(top, "top.txt", sep="\t", , row.names=FALSE)



################################################################
#### export main tables 
write.table(top.10.pka, "top.10.pka.txt", sep="\t", , row.names=FALSE)
write.table(top.30.pka, "top.30.pka.txt", sep="\t", , row.names=FALSE)
write.table(top.60.pka, "top.60.pka.txt", sep="\t", , row.names=FALSE)
write.table(top.300.pka, "top.300.pka.txt", sep="\t", , row.names=FALSE)
write.table(top.600.pka, "top.600.pka.txt", sep="\t", , row.names=FALSE)
write.table(top.900.pka, "top.900.pka.txt", sep="\t", , row.names=FALSE)
write.table(top.1800.pka, "top.1800.pka.txt", sep="\t", , row.names=FALSE)


write.table(top.10, "top.10.txt", sep="\t", , row.names=FALSE)
write.table(top.30, "top.30.txt", sep="\t", , row.names=FALSE)
write.table(top.60, "top.60.txt", sep="\t", , row.names=FALSE)
write.table(top.300, "top.300.txt", sep="\t", , row.names=FALSE)
write.table(top.600, "top.600.txt", sep="\t", , row.names=FALSE)
write.table(top.900, "top.900.txt", sep="\t", , row.names=FALSE)
write.table(top.1800, "top.1800.txt", sep="\t", , row.names=FALSE)


