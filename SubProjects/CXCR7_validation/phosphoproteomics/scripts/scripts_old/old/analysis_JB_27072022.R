############################################
##### script for phosphoproteom analysis
## R pipeline for analyzing pproteom data
## by ÖO and JB



##### nomenclature
## comments
##JB comments by JB
##ÖO comments by ÖO
# outcommented code
##### header 1
#### header 2
### header 3


##### structure of the script
## load data
## filter data
## determine empirical control genes
## reduction of unwanted variance
## determine log2foldchanges of the samples
## downstream analysis

##### to DO/implement 
## mysql intersection



##### install packages####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("RUVSeq")
BiocManager::install("rub") ##ÖO DO we need this? ##JB optional normailization
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("org.Rn.eg.db")
BiocManager::install("EDASeq")
BiocManager::install("PhosR")
install.packages("remotes")
remotes::install_github("biobenkj/ATACseeker") ##ÖO DO WE NEED THIS? ##JB I think not
BiocManager::install("directPA")
BiocManager::install("reactome.db") ##only needed in downstream analysis
BiocManager::install("MSnbase")
BiocManager::install("psych")
BiocManager::install("MSPrep")
install.packages("harmony")


##### load packages
suppressPackageStartupMessages({
  library("BiocManager")
  library("RColorBrewer")
  library("DESeq2")
  library("ruv")
  library("limma")
  library("dplyr")
  library("PhosR")
  library("edgeR")
  library("ggplot2")
  library("reshape2")
  library("calibrate")
  library("limma")
  library("directPA")
  library("annotate")
  library("stringr")
  library("ggplot2")
  library("RUVSeq")
  library("remotes")
  library("org.Rn.eg.db")
  library("MSnbase")
  library("EDASeq")
  library("psych") #for describe function
  library("MSPrep")
  library("harmony")
})



###### load data and preprocess
##OÖ We dont use peptide info for anything? ##JB its important for quality control, backtracking, and for creating final tables 

## before the data are delivered by ISAS (TRR240_A08_01_Phosphoproteome.exe)
## in excel (A08_pp_info.xlxs) relevant columns are selected and the file is seprated in
## A08_phosR.txt and A08_pp_info.txt
## A08_phosR.txt is prprocessed in perl (parserPhosR.sql) to cope with duplicated Uniport_IDs
## the resulting table A08_phosR_v2.txt is loaded to mysql (phosR_gene_symbol_mapping) to get the corresponding protein symbol 
## and the correspondign protein name and description
## an idea is build by uniport_id;smbol;phosphosite;sequence;peptide_id
## and that the table A08_phosR_v3.txt is generated for R analysis


### go to script folder
setwd("C:/Users/Ozge/Nextcloud/Phosphoproteom/phosphoproteom validation/scripts") ##ÖO folder structure should fit.
setwd("D:/Eigene Datein/Phosphoproteom/phosphoproteom validation/scripts/")
setwd("F:/Masterarbeit_BioWis/Phosphoproteomics/Phosphoproteom validation/scripts/")

raw_abundance <- read.table("../data/A08_val_phosphoR_v3.txt", sep="\t", header=TRUE, dec=",")
peptide_info <- read.table("../data/A8_val_pp_info.txt", sep="\t", header=TRUE, dec=",")



### check data
raw_abundance2 <- raw_abundance[,c(-1)]
rownames(raw_abundance2) <- raw_abundance[,1]

raw_abundance2 <- as.data.frame(raw_abundance2)


donor_nr <- gsub("Donor", "",sapply(strsplit(colnames(raw_abundance2), "_"), "[[", 2))
time_point <- sapply(strsplit(colnames(raw_abundance2), "_"), "[[", 3)
time_point <- paste("x", time_point, sep = "")
condition <- sapply(strsplit(colnames(raw_abundance2), "_"), "[[", 4)

colnames(raw_abundance2) <- paste(time_point, condition, donor_nr, sep = "_")

raw_abundance3 <- raw_abundance2 %>% dplyr::select(sort(names(raw_abundance2))) ## sorts the columns
raw_abundance3_log <- log2(raw_abundance3) 


test2 <- as.matrix(raw_abundance3)
test3 <- as.matrix(raw_abundance3_log)

### check distribution and plot distributions
dd <- melt(test2, variable.name = "sample")
dd_log <- melt(test3, variable.name = "sample")

ggplot(dd, aes(value, colour = Var2)) + geom_density(bw = "SJ") +
  theme(legend.position="none") + xlim(0,1e+7)

ggplot(dd_log, aes(value, colour = Var2)) + geom_density(bw = "SJ") +
  theme(legend.position="none") 
#+ xlim(0,1e+2)


desc_data<- describe(raw_abundance3)
c(min(desc_data$skew), max(desc_data$skew))  #skewed to right, skewness all positive [15.74092 31.60951]
c(min(desc_data$kurtosis), max(desc_data$kurtosis)) #heavy tails [292.8265 1293.2373]

desc_logdata <- describe(raw_abundance3_log)
c(min(desc_logdata$skew), max(desc_logdata$skew))  #almost no skewness [0.5383390 0.7913213]
c(min(desc_logdata$kurtosis), max(desc_logdata$kurtosis)) #almost no kurtosis [0.5881861 0.8762273]



##### filter data & imputation of NA values
##### Part A. Preprocessing
##### A.1. Filtering of phosphosites
## imputation as there are still NA and RUV does not work with NA
## use PhoR 
## prepare table in mysql

## Setting up the PhosphoExperiment object 

ppe0 <- PhosphoExperiment(assays = list(Quantification = as.matrix(raw_abundance3_log)))

ppe0 <- PhosphoExperiment(assays = list(Quantification = as.matrix(raw_abundance3_log)), 
                          UniprotID = sapply(strsplit(rownames(ppe0), ";"), "[[", 1),
                          Site = as.numeric(gsub("[A-Z]","", sapply(strsplit(rownames(ppe0), ";"), "[[", 3))),
                          GeneSymbol = sapply(strsplit(rownames(ppe0), ";"), "[[", 2),
                          Residue = gsub("[0-9]","", sapply(strsplit(rownames(ppe0), ";"), "[[", 3)),
                          Sequence = sapply(strsplit(rownames(ppe0), ";"), "[[", 4))

### define the groups and replicates
grps = gsub("_[0-9][0-9]", "", colnames(ppe0))
grps

grps2 <- sapply(strsplit(colnames(ppe0), "_"), "[[",3)
grps2

grps3 <- sapply(strsplit(colnames(ppe0), "_"), "[[", 1)
grps3


##ÖO filtering of phosphosites so that only phosphosites with quantification for
##ÖO at least 50% of the replicates in at least eight of the conditions are retained.

##ÖO how do we decide n here? we have one condition but 8 timepoints  --> here filtering already
##ÖO removes half, and the imputations do not actually impute anything?

ppe_filtered <- selectGrps(ppe0, grps, 0.3, n=7) #same method or 0.3 
## 1266 phosphosites have no NAs.
dim(ppe_filtered)
dim(ppe0)


##### A.2. Imputation of phosphosites ####
#### A.2.1. Site- and condition-specific imputation

set.seed(123)
ppe_imputed_tmp <- scImpute(ppe_filtered, 0.3, grps)[,colnames(ppe_filtered)]


## stats
test0 <- SummarizedExperiment::assay(ppe0,"Quantification")
test <- SummarizedExperiment::assay(ppe_filtered,"Quantification")
test1 <- SummarizedExperiment::assay(ppe_imputed_tmp,"Quantification")
test2 <- SummarizedExperiment::assay(ppe_imputed_tmp,"imputed")

sum(is.na(test0))
sum(is.na(test))
sum(is.na(test1))
sum(is.na(test2))


dim(ppe0)
dim(ppe_filtered)
dim(ppe_imputed_tmp)


#### A.2.2. Paired tail-based imputation (second imputing is optional) ##JB isnt used here so jump over it
 
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



##ÖO what is this for? ##JB just need if we jump over pair tailed imputation
ppe_imputed <- ppe_imputed_tmp


### remove NA if any left
#SummarizedExperiment::assay(ppe_imputed, "imputed") <- na.omit(SummarizedExperiment::assay(ppe_imputed, "imputed"))
abundance_no.na <- na.omit(SummarizedExperiment::assay(ppe_imputed, "imputed"))


### generate phosphoproteom object again? ##JB why here again?
ppe_imputed_omit <- PhosphoExperiment(assays = list(imputed = as.matrix(abundance_no.na)))

UniprotID(ppe_imputed_omit)  <- sapply(strsplit(rownames(ppe_imputed_omit), ";"), "[[", 1)
GeneSymbol(ppe_imputed_omit) <- sapply(strsplit(rownames(ppe_imputed_omit), ";"), "[[", 2)
Residue(ppe_imputed_omit) <- gsub("[0-9]","", sapply(strsplit(rownames(ppe_imputed_omit), ";"), "[[", 3))
Site(ppe_imputed_omit) <- gsub("[A-Z]","", sapply(strsplit(rownames(ppe_imputed_omit), ";"), "[[", 3))
Sequence(ppe_imputed_omit) <- sapply(strsplit(rownames(ppe_imputed_omit), ";"), "[[", 4)

ppe_imputed <- ppe_imputed_omit
#SummarizedExperiment::assay(ppe_imputed, "Quantification") <- SummarizedExperiment::assay(ppe_filtered, "Quantification")


### scaling imputation
#ppe_imputed_scaled <- medianScaling(ppe_imputed_omit, scale = FALSE, assay = "imputed")
ppe_imputed_scaled <- medianScaling(ppe_imputed_omit, scale = FALSE, assay = "imputed")

test3 <- SummarizedExperiment::assay(ppe_imputed,"imputed")
test3 <- ppe_imputed_scaled@assays@data@listData[["imputed"]]
test4 <- SummarizedExperiment::assay(ppe_imputed_scaled,"scaled")


sum(is.na(test))
sum(is.na(test2))
sum(is.na(test3))
sum(is.na(test4))
dim(test4)





## PCA: grps redefined for visualization
plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), grps=grps, 
       labels = colnames(ppe_imputed_scaled), panel = "pca") +
  ggplot2::ggtitle("Before batch correction")

plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), grps=grps2, 
       labels = colnames(ppe_imputed_scaled), panel = "pca") +
  ggplot2::ggtitle("Before batch correction")


## quantification before and after impuation

p0 = plotQC(SummarizedExperiment::assay(ppe0,"Quantification"), 
            labels=colnames(ppe0), 
            panel = "quantify", grps = grps)
p1 = plotQC(SummarizedExperiment::assay(ppe_filtered,"Quantification"), 
            labels=colnames(ppe_filtered), 
            panel = "quantify", grps = grps)
p2 = plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), 
            labels=colnames(ppe_imputed_scaled), panel = "quantify", grps = grps)
ggpubr::ggarrange(p0, p1, p2, nrow = 3)

## hierachical clustering before and after impuation
p0 = plotQC(SummarizedExperiment::assay(ppe0,"Quantification"), 
            labels=colnames(ppe_filtered), panel = "dendrogram", 
            grps = grps)

p1 = plotQC(SummarizedExperiment::assay(ppe_filtered,"Quantification"), 
            labels=colnames(ppe_filtered), panel = "dendrogram", 
            grps = grps)
p2 = plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), 
            labels=colnames(ppe_imputed_scaled), 
            panel = "dendrogram", grps = grps)

ggpubr::ggarrange(p0, p1, p2, nrow = 3)

p1 = plotQC(SummarizedExperiment::assay(ppe_filtered,"Quantification"), 
            labels=colnames(ppe_filtered), panel = "dendrogram", 
            grps = grps2)
p2 = plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), 
            labels=colnames(ppe_imputed_scaled), 
            panel = "dendrogram", grps = grps2)
ggpubr::ggarrange(p1, p2, nrow = 1)



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

#### make a betweenLane Normalization
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
empirical_topall <- rownames(dataset_norm)[which((!rownames(dataset_norm) %in% rownames(top_all_newmethod)[1:4085]))]
#ncol(contrast.matrix_new)


a <- rownames(dataset_norm)[which(!(rownames(dataset_norm) %in% rownames(top_all_newmethod)[1:3100]))]
#rownames(dataset_norm)[which(!(rownames(dataset_norm) %in% rownames(top_0_1800)[1:2000]))]

### differential analysis pairwise to find empirical control genes

bottom_list = c()

x=matrix(data=NA, nrow=nrow(fit2_new), ncol=ncol(contrast.matrix_new))

for (i in 1:ncol(contrast.matrix_new)) {
    df<- topTable(fit2_new, coef = i, number=nrow(fit2_new),  
                  adjust.method="BH", sort.by="M", p.value=1)
    #df2 <- df[df$adj.P.Val <= 0.05,]
    new_df <- df[ order(row.names(df)), ]
    x[,i]= new_df$adj.P.Val
    bottom <- rownames(dataset_norm)[which(!(rownames(dataset_norm) %in% rownames(df)[1:2500]))]
    assign(paste('top',colnames(contrast.matrix_new)[i],sep='_'),df)
    assign(paste('newdf', i, sep="_"), new_df)
    #assign(paste('pvals', i, sep="_"), pvals)
    assign(paste('unsig',colnames(contrast.matrix_new)[i],sep='_'),bottom)
    assign(paste('unsig',i,sep='_'),bottom)
    bottom_list <- c(bottom_list, paste('unsig', i,sep='_' ))
}

x <- rowMedians(x)
x<- as.data.frame(x)
rownames(x) <- rownames(top_x0sek_Ctrl_vs_x10sek_CXCR7[order(row.names(top_x0sek_Ctrl_vs_x10sek_CXCR7)), ])
#rownames(x) <- rownames(top_01_vs_02[order(row.names(top_01_vs_02)), ])

colnames(x) <- "medianp"
x$name <- rownames(x)


x_filt<-x[x$medianp > 0.995, ]
empirical_x_filt <- rownames(x_filt)
length(empirical_x_filt)

emp_intersect_filt <- intersect(empirical_topall, empirical_x_filt)
length(emp_intersect_filt)

#idea: if we're seeing donors cluster together, can't we find the genes that do not change between the donors
#and use those to normalize?
##JB why here 45 comparisons cause the contrast matrix has 21 comparisons?
##ÖO answer: here i tried to group them by donor and make comparisons between donors, so 10 donors, 45 pairwise combinations, but it didn't work better

#empirical_pairwise_new <- Reduce(intersect, list(unsig_1, unsig_2, unsig_3, unsig_4, unsig_5, unsig_6, unsig_7, unsig_8, unsig_9, 
#                       unsig_10, unsig_11, unsig_12, unsig_13, unsig_14, unsig_15, unsig_16, unsig_17, unsig_18,
#                       unsig_19, unsig_20, unsig_21, unsig_22, unsig_23, unsig_24, unsig_25, unsig_26, unsig_27, unsig_28, 
#                       unsig_29, unsig_30, unsig_31, unsig_32, unsig_33, unsig_34, unsig_35, unsig_36, 
#                       unsig_37, unsig_38, unsig_39, unsig_40, unsig_41, unsig_42, unsig_43, unsig_44, unsig_45))
					
empirical_pairwise_new <- Reduce(intersect, list(unsig_1, unsig_2, unsig_3, unsig_4, unsig_5, unsig_6, unsig_7, unsig_8, unsig_9, 
                       unsig_10, unsig_11, unsig_12, unsig_13, unsig_14, unsig_15, unsig_16, unsig_17, unsig_18,
                       unsig_19, unsig_20, unsig_21))
length(empirical_pairwise_new)


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

#Spss

data("SPSs")

SPSs2 <- gsub(".{1}$", "", SPSs)


##### RUV normalization
design = model.matrix(~ grps - 1) ##ÖO why do we define model matrix like this? it's the same (removes intercept) ##JB ja thats double. sometime model matrix differs either using limma or edgeR. we can remove that.
#design

#ctl0 = which(sites2 %in% SPSs)
#ctl0 = which(sites %in% SPSs)
#ctl0 = which(sites %in% inhouse_SPSs)
#specified_control <- rownames(ppe_imputed_scaled[ctl0,])

 
#empiricalcontrol <-  which(rownames(ppe_imputed_scaled) %in% dual)
#length(empiricalcontrol)


#ctlvar = which(rownames(ppe_imputed_scaled) %in% ControlGenes3)

ctl2 = which(rownames(ppe_imputed_scaled) %in% empirical_topall)
ctl3 <- which(paste(sapply(strsplit(rownames(ppe_imputed_scaled), ";"), "[[", 2), 
                    sapply(strsplit(rownames(ppe_imputed_scaled), ";"),
                           "[[", 3), sep = ";") %in%SPSs2) #ä#ÖO only 9 spss in our data, just wanted to check what it does but doesnt help.
ctl <- c(ctl2, ctl3)  #ÖO doesnt help neither.
#ctl3 = which(rownames(ppe_imputed_scaled) %in% empiricalcontrol)
#ctl4 <- c(ctl,ctl2)
#empirical2 <- c(empirical,specified_control)


#### RUVphospho normalization
## RUVphospho shows a good tendency for statistics, but the intensity value are curiously shifted 
## but still 
ppe = RUVphospho(ppe_imputed_scaled, M = design, k = 23, ctl = ctl2)


#### RUVg normalization
abundance <- 2^SummarizedExperiment::assay(ppe_imputed_scaled,"scaled")
abundance2 = apply(abundance, 1:2, function(x) round(x))
set <- newSeqExpressionSet(as.matrix(abundance2),phenoData = data.frame(x_tt, row.names=colnames(abundance2)))

#set2_mk <- RUVg(set0_mk, Control_mk, k=7)
set2 <- RUVg(set, empirical_topall, k=23) ##ÖO how many dimensions (replicates) to delete to get better clustering?
SummarizedExperiment::assay(ppe_imputed_scaled,"normalised") <- log2(set2@assayData[["normalizedCounts"]])

ppe <- ppe_imputed_scaled


#### visualization of the influence of the value k on normalization
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



#### statistics of normalization and quality control
### get different datasets to compare to
test2 <- SummarizedExperiment::assay(ppe0,"Quantification")
test3 <-SummarizedExperiment::assay(ppe,"scaled")
test4 <-SummarizedExperiment::assay(ppe,"normalised")

plotQC(test4, grps=grps, 
       labels = colnames(ppe), panel="pca")
plotQC(test4, grps=grps, 
       labels = colnames(ppe), panel="dendrogram")+
  ggtitle("ruvg normalized")



##### msNormalize normalization


#medianNormalizedDF <- msNormalize(hmImputedDF, normalizeMethod = "median",
#                                  compVars = c("mz", "rt"),
#                                  sampleVars = c("spike", "batch", 
#                                  "subject_id"),
#                                  separator = "_")

test5 <- as.data.frame(test4)
test5$id2 <- rownames(test4)


medianNormalizedDF <- msNormalize(test5, normalizeMethod = "median",
                                  compVars = "id2",
								  #nComp = 2,
                                  sampleVars = c("subject_id"),
                                  separator = c("time","batch","subject_id"))


SummarizedExperiment::assay(ppe_imputed_scaled,"normalised") <- medianNormalizedDF[,-1]

test_msnorm <-SummarizedExperiment::assay(ppe_imputed_scaled,"normalised")

plotQC(test_msnorm, grps=grps, 
             labels = colnames(ppe), panel="pca")
msg <- plotQC(test_msnorm, grps=grps, 
           labels = colnames(ppe), panel="dendrogram")+
  ggtitle("ms normalized")

ggpubr::ggarrange(t4g, msg, nrow = 2)


ppe <- ppe_imputed_scaled

##### VSN normalization

test5 <- normalizeVSN(test3)
SummarizedExperiment::assay(ppe_imputed_scaled,"normalised") <- test5
ppe <- ppe_imputed_scaled


####
## iterative second emp
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
####

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



##### differential analysis
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


#(treat10sec - ctrl0sec) - (ctrl10sec - ctrl0sec)

## tt10
#contrast.matrix <- makeContrasts(((grpsx10sek_CXCR7-grpsx0sek_Ctrl)-(grpsx10sek_DMSO-grpsx0sek_Ctrl)), levels=design)
contrast.matrix <- makeContrasts((grpsx10sek_CXCR7-grpsx10sek_DMSO), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.10 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)

##  tt600
contrast.matrix <- makeContrasts((grpsx600sek_CXCR7-grpsx600sek_DMSO), levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend=TRUE)
top.600 <- topTable(fit2, number=nrow(fit2), adjust.method="fdr", sort.by="p", p.value=1)


## tt1800
contrast.matrix <- makeContrasts(grpsx1800sek_CXCR7-grpsx1800sek_DMSO, levels=design)
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


top.10.sign <- top.10[top.10[, "adj.P.Val"] <0.05,]
top.600.sign <-  top.600[top.600[, "adj.P.Val"] <0.05,]
top.1800.sign <-  top.1800[top.1800[, "adj.P.Val"] <0.05,]


#### statistics
### number of significant dergulated targets
DE2.RUV <- c(length(rownames(top.10.sign)),length(rownames(top.600.sign)),length(rownames(top.1800.sign)))
DE2.RUV


top.600["P50552;VASP;S239;KVSKQEEASGGPTAPK;9984",]
top.600.vasp <- top.600[top.900$uniprot=="P50552",]
top.10.vasp <- top.10[top.900$uniprot=="P50552",]
top.1800.vasp <- top.1800[top.1800$uniprot=="P50552",]


common_differences <- Reduce(intersect, list(rownames(top.10.sign), rownames(top.600.sign), rownames(top.1800.sign)))


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

#
top <- cbind(rownames(new_10),new_10, new_30, new_60, new_300, new_600, new_900, new_1800)

write.table(top, "top.txt", sep="\t", , row.names=FALSE)


#### export main tables 
write.table(top.10, "top.10.txt", sep="\t", , row.names=FALSE)
write.table(top.30, "top.30.txt", sep="\t", , row.names=FALSE)
write.table(top.60, "top.60.txt", sep="\t", , row.names=FALSE)
write.table(top.300, "top.300.txt", sep="\t", , row.names=FALSE)
write.table(top.600, "top.600.txt", sep="\t", , row.names=FALSE)
write.table(top.900, "top.900.txt", sep="\t", , row.names=FALSE)
write.table(top.1800, "top.1800.txt", sep="\t", , row.names=FALSE)


