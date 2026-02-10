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
## psiquic network 
## GO enrichment
## kinase targets


####################################################################################
##### install packages
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
BiocManager::install("PSICQUIC")


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
  library("PSICQUIC")
  library("sva")  #package SVA
})


####################################################################################
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
#setwd("C:/Users/Ozge/Nextcloud/Phosphoproteom/phosphoproteom validation/scripts") ##ÖO folder structure should fit.
#setwd("D:/Eigene Datein/Phosphoproteom/phosphoproteom validation/scripts/")
setwd("F:/Masterarbeit_BioWis/Phosphoproteomics/Phosphoproteom validation/scripts/")


## read tables locally or
## RS: connect to mysql server and load the table from the mysql database
raw_abundance <- read.table("../data/raw_data/A08_val_phosphoR_v3.txt", sep="\t", header=TRUE, dec=",")
peptide_info <- read.table("../data/raw_data/A8_val_pp_info.txt", sep="\t", header=TRUE, dec=",")

## format
raw_abundance2 <- raw_abundance[,c(-1)]
rownames(raw_abundance2) <- raw_abundance[,1]
raw_abundance2 <- as.data.frame(raw_abundance2)

## define header
donor_nr <- gsub("Donor", "",sapply(strsplit(colnames(raw_abundance2), "_"), "[[", 2))
time_point <- sapply(strsplit(colnames(raw_abundance2), "_"), "[[", 3)
time_point <- paste("x", time_point, sep = "")
condition <- sapply(strsplit(colnames(raw_abundance2), "_"), "[[", 4)
colnames(raw_abundance2) <- paste(time_point, condition, donor_nr, sep = "_")

## sort
raw_abundance3 <- raw_abundance2 %>% dplyr::select(sort(names(raw_abundance2))) ## sorts the columns
raw_abundance3_log <- log2(raw_abundance3)

write.table(raw_abundance3, "../data/processed_data/raw_intensities.txt", sep = "\t")

## check distribution and plot distributions
test2 <- as.matrix(raw_abundance3)
test3 <- as.matrix(raw_abundance3_log)

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



## Setting up the PhosphoExperiment object 

ppe0 <- PhosphoExperiment(assays = list(Quantification = as.matrix(raw_abundance3_log)))

ppe0 <- PhosphoExperiment(assays = list(Quantification = as.matrix(raw_abundance3_log)), 
                          UniprotID = sapply(strsplit(rownames(ppe0), ";"), "[[", 1),
                          Site = as.numeric(gsub("[A-Z]","", sapply(strsplit(rownames(ppe0), ";"), "[[", 3))),
                          GeneSymbol = sapply(strsplit(rownames(ppe0), ";"), "[[", 2),
                          Residue = gsub("[0-9]","", sapply(strsplit(rownames(ppe0), ";"), "[[", 3)),
                          Sequence = sapply(strsplit(rownames(ppe0), ";"), "[[", 4))

UniprotID(ppe_imputed_omit)  <- sapply(strsplit(rownames(ppe_imputed_omit), ";"), "[[", 1)
GeneSymbol(ppe_imputed_omit) <- sapply(strsplit(rownames(ppe_imputed_omit), ";"), "[[", 2)
Residue(ppe_imputed_omit) <- gsub("[0-9]","", sapply(strsplit(rownames(ppe_imputed_omit), ";"), "[[", 3))
Site(ppe_imputed_omit) <- gsub("[A-Z]","", sapply(strsplit(rownames(ppe_imputed_omit), ";"), "[[", 3))
Sequence(ppe_imputed_omit) <- sapply(strsplit(rownames(ppe_imputed_omit), ";"), "[[", 4)


### define the groups and replicates
grps = gsub("_[0-9][0-9]", "", colnames(ppe0))
grps

grps2 <- sapply(strsplit(colnames(ppe0), "_"), "[[",2)
grps2

grps3 <- sapply(strsplit(colnames(ppe0), "_"), "[[", 1)
grps3

grps4 <- sapply(strsplit(colnames(ppe0), "_"), "[[", 3)
grps4


## filter data & imputation of NA values
ppe_filtered <- selectGrps(ppe0, grps, 0.3, n=7) #same method or 0.3 
# 1266 phosphosites have no NAs.
dim(ppe_filtered)

set.seed(123)
ppe_imputed <- scImpute(ppe_filtered, 0.3, grps)
dim(ppe_imputed)

ppe_imputed_tmp <- PhosphoExperiment(assays = list(imputed = as.matrix(SummarizedExperiment::assay(ppe_imputed,"imputed")))

ppe_imputed_scaled <- medianScaling(ppe_imputed_tmp, scale = FALSE, assay = "imputed")
dim(ppe_imputed_scaled)



write.table(2^ppe_imputed_scaled@assays@data@listData[["scaled"]], "../data/processed_data/intensities_imputed_scaled.txt", sep = "\t")


## stats
test0 <- SummarizedExperiment::assay(ppe0,"Quantification")
test <- SummarizedExperiment::assay(ppe_filtered,"Quantification")
test1 <- SummarizedExperiment::assay(ppe_imputed,"imputed")
test2 <- SummarizedExperiment::assay(ppe_imputed_scaled,"scaled")

sum(is.na(test0))
sum(is.na(test))
sum(is.na(test1))
sum(is.na(test2))


#### Visualization

## hierachical cluster
plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), panel = "dendrogram", grps=grps, labels = colnames(ppe_imputed_scaled)) +
  ggplot2::ggtitle("Before batch correction")


## PCA: grps redefined for visualization
plotQC(SummarizedExperiment::assay(ppe_imputed,"imputed"), grps=grps4,
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
