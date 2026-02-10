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

# install packages ----
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
BiocManager::install("RUVSeq")
BiocManager::install("clusterProfiler")


# load packages ----
suppressPackageStartupMessages({
  library("BiocManager")
  library("PhosR")
  library("RColorBrewer")
  library("DESeq2")
  library("ruv")
  library("RUVSeq")
  library("limma")
  library("dplyr")
  library("ggplot2")
  library("reshape2")
  library("limma")
  library("directPA")
  library("annotate")
  library("stringr")
  library("ggplot2")
  library("remotes")
  library("psych") #for describe function
  library("plyr")
  library(clusterProfiler)
  library(ggvenn)
  library("RMySQL")
  library(dplyr)
  library(cowplot)
})


# load data and preprocess ----
##OÖ We dont use peptide info for anything? ##JB its important for quality control, backtracking, and for creating final tables 

## before the data are delivered by ISAS (TRR240_A08_01_Phosphoproteome.exe)
## in excel (A08_pp_info.xlxs) relevant columns are selected and the file is seprated in
## A08_phosR.txt and A08_pp_info.txt
## A08_phosR.txt is prprocessed in perl (parserPhosR.sql) to cope with duplicated Uniport_IDs
## the resulting table A08_phosR_v2.txt is loaded to mysql (phosR_gene_symbol_mapping) to get the corresponding protein symbol 
## and the correspondign protein name and description
## an idea is build by uniport_id;smbol;phosphosite;sequence;peptide_id
## and that the table A08_phosR_v3.txt is generated for R analysis

## map the gene name
## load organism database
organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)


## go to script folder
setwd("/Users/ozgeosmanoglu/Nextcloud/Phosphoproteom/phosphoproteom analysis 2/scripts") ##ÖO folder structure
setwd("F:/Masterarbeit_BioWis/Phosphoproteomics/phosphoproteom analysis 2") ## JB Laptop privat
setwd("C:/Users/job37yv/Research/Phospho/phosphoproteom analysis 2")  ## JB Laptop RVZ
setwd("D:/Eigene Datein/Phospho/Phospho.analysis.A08.v2") ## JB Computer RVZ

## read tables locally or
## RS: connect to mysql server and load the table from the mysql database
raw_data <- read.table("../data/raw_data/A08_phosR.txt", sep="\t", header=TRUE, dec=",")
peptide_info <- read.table("../data/raw_data/A08_pp_info.txt", sep="\t", header=TRUE, dec=",")

raw_abundance <- read.table("../data/raw_data/A08_phosR_v2.txt", sep="\t", header=TRUE, dec=",", fileEncoding = "UCS-2LE")

## filter psites that are not annotated
colnames(raw_abundance)[1:4] <- c("peptide_id", "uniprot_id", "sequence", "psite")
raw_abundance = raw_abundance[-which(raw_abundance$psite == ""),]


uniprot_id <- raw_abundance$uniprot
psite <- raw_abundance$psite
sequence <- raw_abundance$sequence
peptide_id <- raw_abundance$peptide_id
#peptide_id <- sapply(strsplit(rownames(raw_abundance2), ";"), "[[", 4)

keytypes(org.Hs.eg.db)
uniprot2symbol<-bitr(uniprot_id, fromType = "UNIPROT", toType = "SYMBOL", OrgDb=organism, )

uniprot_id2 <- as.data.frame(uniprot_id)
colnames(uniprot_id2) <- "UNIPROT"
uniprot_id2$SYMBOL  <- "NA"



### if update funtion is like this in R maybe we dont need mysql
### drop the mysql tables any time

con = dbConnect(RMySQL::MySQL(),
                dbname='phospho_A08_v2',
                host='localhost',
                port=3306,
                user='lumpi',
                password='ipmul')

dbListTables(con)

dbRemoveTable(con,"uniprot_id2")
dbRemoveTable(con,"uniprot2symbol")
dbWriteTable(con, value = uniprot_id2, name = "uniprot_id2", append = TRUE)
dbWriteTable(con, value = uniprot2symbol, name = "uniprot2symbol", append = TRUE)

result = dbSendQuery(con, "select * from uniprot_id2")
data.frame = fetch(result, n=5)
print(data.frame)
dbClearResult(result)

result2 = dbSendQuery(con, "UPDATE uniprot_id2
       inner join uniprot2symbol
       ON uniprot_id2.UNIPROT = uniprot2symbol.UNIPROT
       set uniprot_id2.symbol = uniprot2symbol.symbol")
data.frame = fetch(result2, n=5)
print(data.frame)
dbClearResult(result2)

result = dbSendQuery(con, "select * from uniprot_id2")
uniprot2symbol2 = fetch(result, n = -1 )
print(uniprot2symbol2)
dbClearResult(result)

### or similar with likewise dplyr but different row number

### or similar with likewise dplyr but different row number

#updated_df <- uniprot_id2 %>%
#  inner_join(uniprot2symbol, by = c("UNIPROT")) %>%
#  mutate(SYMBOL.x = SYMBOL.y) %>%
#  select(-SYMBOL.y) %>%
#  rename(SYMBOL = SYMBOL.x)



## raw_abundance2 <- raw_abundance
raw_abundance2 <- raw_abundance[,-c(1,2,3,4)]
rownames(raw_abundance2) <- paste(uniprot2symbol2$UNIPROT, uniprot2symbol2$SYMBOL, psite, sequence, peptide_id, sep = ";")

## format
# raw_abundance2 <- raw_abundance[,c(-1)]
# rownames(raw_abundance2) <- raw_abundance[,1]
raw_abundance2 <- as.data.frame(raw_abundance2)

## further filter (sites with several uniport IDs)
not.uniq = read.delim(file = "../data/processed_data/not_unique_ids.txt", sep = "\t") #ÖO did it manually

raw_abundance2 <- raw_abundance2[!(row.names(raw_abundance2) %in% not.uniq[,1]), ]

write.table(raw_abundance2, "../data/processed_data/raw_abundance2.txt", sep = "\t", row.names = TRUE)

## define header
donor_nr <- gsub("Donor", "",sapply(strsplit(colnames(raw_abundance2), "_"), "[[", 2))
time_point <- sapply(strsplit(colnames(raw_abundance2), "_"), "[[", 1)
#time_point <- paste("x", time_point, sep = "")
timepoint_fac <- as.factor(gsub("X", "",time_point))
#timepoint_fac <- gsub("sek", "",timepoint_fac)
condition <- gsub("0000", "ctrl",timepoint_fac)
condition <- ifelse(grepl("ctrl", condition), "ctrl", "CXCR7")
#colnames(raw_abundance2) <- paste(time_point, condition, donor_nr, sep = "_")


sample_info <- data.frame(donor_nr, timepoint_fac, condition, colnames(raw_abundance2))
colnames(sample_info) <- c("donor", "time", "condition", "sample")
rownames(sample_info) <- colnames(raw_abundance2)

## sort
input <- raw_abundance2

raw_abundance3 <- input %>% dplyr::select(sort(names(input))) ## sorts the columns
raw_abundance3_log <- log2(raw_abundance3)

############ filter psites that are also in the old dataset

#top.all_old <- read.table("../phosphoproteom analysis 2/data/processed_data/top.all.txt", 
#                        sep="\t", header=TRUE, dec=",")

#top.all_old_names <- paste(sapply(strsplit(rownames(top.all_old), ";"), "[[", 2),
#                    sapply(strsplit(rownames(top.all_old), ";"), "[[", 3), sep = ";")

#raw_abundance3_log_names <- paste(sapply(strsplit(rownames(raw_abundance3_log), ";"), "[[", 2),
#                    sapply(strsplit(rownames(raw_abundance3_log), ";"), "[[", 3), sep = ";")

#which(raw_abundance3_log_names %in% top.all_old_names)
#length(which(raw_abundance3_log_names %in% top.all_old_names))
#length(raw_abundance3_log_names)
#length(top.all_old_names)

#sel <- which(raw_abundance3_log_names %in% top.all_old_names)

#raw_abundance3_log <- raw_abundance3_log[sel,]


## write
write.table(raw_abundance3, "../data/processed_data/raw_abundance3.txt", sep = "\t")
write.table(raw_abundance3_log, "../data/processed_data/raw_intensities.txt", sep = "\t")

# check correlations ----
library(ComplexHeatmap)
library(corrplot)

colorBatch = function(batch){
  batchCol =  unlist(ifelse(length(unique(batch)) == 1, rainbow(1, start = 0.5), list(rainbow(length(unique(batch))))))
  names(batchCol) = sort(unique(batch))
  return(batchCol)
}
colorGroup = function(groups){
  groupCol = viridis::viridis(length(unique(groups)))
  names(groupCol) = sort(unique(groups))
  return(groupCol)
}
heatmapCorr = function(data, groups, batch, sampleLabels){
  cor_mat <- cor(data, use = "pairwise.complete.obs", method = "spearman")
  
  ColAnn <- HeatmapAnnotation(Sample = groups, 
                              col = list(Sample = colorGroup(groups)), 
                              annotation_legend_param = list(
                                Sample = list(
                                  title = "Sample",
                                  at = sort(unique(groups)),
                                  labels = paste("Group", sort(unique(groups)))
                                )
                              ))
  RowAnn <- rowAnnotation(Batch = batch, 
                          col = list(Batch = colorBatch(batch)), 
                          annotation_legend_param = list(
                            Batch = list(
                              title = "Donor",
                              at = sort(unique(batch)),
                              labels = paste("Batch", sort(unique(batch)))
                            )
                          ))
  hm_corr = Heatmap(cor_mat, 
                    col = circlize::colorRamp2(seq(min(cor_mat), 1, ((1 - min(cor_mat))/7)),RColorBrewer::brewer.pal(8, "Reds")), 
                    heatmap_legend_param = list(color_bar = "continuous", 
                                                legend_direction = "horizontal", 
                                                legend_width = unit(5, "cm"), 
                                                title_position = "topcenter"), 
                    name = "Spearman correlation", 
                    column_names_gp = gpar(fontsize = 12), 
                    row_names_gp = gpar(fontsize = 12), 
                    top_annotation = ColAnn, 
                    left_annotation = RowAnn,
                    column_labels = sampleLabels, row_labels = sampleLabels,
                    cluster_columns = F,
                    cluster_rows = F
  )
  
  draw(hm_corr, heatmap_legend_side = "top")
}

##prep data
data <- as.data.frame(raw_abundance3)
##or
data <- as.data.frame(norm_intensity)

#continue
donor_nr <- sapply(strsplit(colnames(data), "_"), "[[", 2)
time_point <- sapply(strsplit(colnames(data), "_"), "[[", 1)
timepoint_fac <- as.factor(time_point)
colnames(data) <- paste(time_point, donor_nr, sep = "_")

#prep metadata
metadata <-  data.frame(donor_nr, timepoint_fac, colnames(data))
colnames(metadata) <- c("donor", "time", "sample")
metadata$group <- gsub("X", "", metadata$time)
rownames(metadata) <- colnames(data)

#prep heatmap inputs
groups <- metadata$group
batch <- metadata$donor
sampleLabels <- metadata$sample

tiff("../analysis/PCA/correlation_norm_style2.tiff", 
     width = 20, height = 12, units = 'in', res = 600, compression = "lzw")
heatmapCorr(data, groups, batch, sampleLabels)
dev.off()

# Stefans script for corr. ----
pearsonInput <-  data

#Calculate Pearson CorrelationMatrix on log10 intensities
pearson <- matrix(ncol=53,nrow=53)
for(i in 1:ncol(pearson))	{
  for(j in 1:nrow(pearson))	{
    pearson[i,j] <- cor(pearsonInput[,i],pearsonInput[,j],use = "pairwise.complete.obs",method="pearson")
    print(paste("i=",i,"j"=j))
    flush.console()
  }
}

cols1 <-  colorRampPalette(c("yellow","firebrick1"))(200)

tiff("../analysis/PCA/correlation_norm_stefanStyle.tiff", 
     width = 12, height = 10, units = 'in', res = 600, compression = "lzw")
corrplot(pearson,is.corr=F,col=cols1,tl.col="black",diag=T,addCoef.col="black",method="color",number.cex=0.5)
dev.off()


# check distribution and plot distributions ----
test2 <- as.matrix(raw_abundance3)
test3 <- as.matrix(raw_abundance3_log)

dd <- reshape2::melt(test2, variable.name = "sample")
dd_log <- reshape2::melt(test3, variable.name = "sample")

ggplot(dd, aes(value, colour = Var2)) + geom_density(bw = "sj") +
  theme_cowplot() +
  theme(legend.position="none")
#+ xlim(0,1e+7)

ggplot(dd_log, aes(value, colour = Var2)) + geom_density(bw = "sj") +
  theme_cowplot() +
  theme(legend.position="none")
#+ xlim(0,1e+7)

ggsave("../analysis/PCA/DensityPlotLog.tiff",
      device = "tiff",
      bg = "white",
      dpi = 300,
      width = 6,
      height = 5)

desc_data<- describe(raw_abundance3)
c(min(desc_data$skew), max(desc_data$skew))  #skewed to right, skewness all positive [17.07702 28.59698]
c(min(desc_data$kurtosis), max(desc_data$kurtosis)) #heavy tails [363.3107 1079.3505]

desc_logdata <- describe(raw_abundance3_log)
c(min(desc_logdata$skew), max(desc_logdata$skew))  #almost no skewness [0.5449089 0.7728712]
c(min(desc_logdata$kurtosis), max(desc_logdata$kurtosis)) #almost no kurtosis [0.6261776 1.1108691]


## Setting up the PhosphoExperiment object 
ppe0 <- PhosphoExperiment(assays = list(Quantification = as.matrix(raw_abundance3_log)))

ppe0 <- PhosphoExperiment(assays = list(Quantification = as.matrix(raw_abundance3_log)), 
                          UniprotID = sapply(strsplit(rownames(ppe0), ";"), "[[", 1),
                          Site = as.numeric(gsub("[A-Z]","", sapply(strsplit(rownames(ppe0), ";"), "[[", 3))),
                          GeneSymbol = sapply(strsplit(rownames(ppe0), ";"), "[[", 2),
                          Residue = gsub("[0-9]","", sapply(strsplit(rownames(ppe0), ";"), "[[", 3)),
                          Sequence = sapply(strsplit(rownames(ppe0), ";"), "[[", 4))


### define the groups and replicates
grps = gsub("_[0-9]{1}", "", colnames(ppe0))
grps


grps2 <- sapply(strsplit(colnames(ppe0), "_"), "[[",2)
grps2

#grps3 <- sapply(strsplit(colnames(ppe0), "_"), "[[", 1)
#grps3

#grps4 <- sapply(strsplit(colnames(ppe0), "_"), "[[", 3)
#grps4


##ÖO filtering of phosphosites so that only phosphosites with quantification for
##ÖO at least 50% of the replicates in at least eight of the conditions are retained.

##ÖO how do we decide n here? we have one condition but 8 timepoints  --> here filtering already
##ÖO removes half, and the imputations do not actually impute anything?

ppe_filtered <- selectGrps(ppe0, grps, 0.35, n=8) #similar methods, 

dim(ppe_filtered)
dim(ppe0)


##### A.2. Imputation of phosphosites ####
#### A.2.1. Site- and condition-specific imputation

set.seed(123)
ppe_imputed <- scImpute(ppe_filtered, 0.35, grps)[,colnames(ppe_filtered)]
dim(ppe_imputed)


#### A.2.2. Paired tail-based imputation (second imputing is optional) ##JB isnt used here so jump over it

#vect1 <- c(7, 6, 7, 7, 7, 6, 7, 6)

# set.seed(123)
# ppe_imputed <- ppe_imputed_tmp
# 
# SummarizedExperiment::assay(ppe_imputed[,seq(1,10)],"imputed") <- PhosR::ptImpute(SummarizedExperiment::assay(ppe_imputed[,seq(11,70)],"imputed"),
#                                                                                   SummarizedExperiment::assay(ppe_imputed[,seq(1,10)],"imputed"),
#                                                                                   percent1 = 1, percent2 = 0.3, paired = FALSE)
# 
# SummarizedExperiment::assay(ppe_imputed[,seq(11,20)],"imputed") <- PhosR::ptImpute(SummarizedExperiment::assay(ppe_imputed[,seq(21,30)],"imputed"),
#                                                                                    SummarizedExperiment::assay(ppe_imputed[,seq(11,20)],"imputed"),
#                                                                                    percent1 = 1, percent2 = 0.3, paired = FALSE)
# 
# 
# SummarizedExperiment::assay(ppe_imputed[,seq(31,40)],"imputed") <- PhosR::ptImpute(SummarizedExperiment::assay(ppe_imputed[,seq(41,50)],"imputed"),
#                                                                                    SummarizedExperiment::assay(ppe_imputed[,seq(31,40)],"imputed"),
#                                                                                    percent1 = 1, percent2 = 0.3, paired = FALSE)
# 
# SummarizedExperiment::assay(ppe_imputed[,seq(51,60)],"imputed") <- PhosR::ptImpute(SummarizedExperiment::assay(ppe_imputed[,seq(61,70)],"imputed"),
#                                                                                    SummarizedExperiment::assay(ppe_imputed[,seq(51,60)],"imputed"),
#                                                                                    percent1 = 1, percent2 = 0.3, paired = FALSE)
# 

ppe_imputed_tmp <- PhosphoExperiment(assays = list(imputed = as.matrix(SummarizedExperiment::assay(ppe_imputed,"imputed"))))


## remove NA if any left
#abundance_no.na <- na.omit(SummarizedExperiment::assay(ppe_imputed_tmp, "imputed"))
#ppe_imputed_omit <- PhosphoExperiment(assays = list(imputed = as.matrix(abundance_no.na)))
#ppe_imputed_tmp <- ppe_imputed_omit



## scaling imputation
ppe_imputed_scaled <- medianScaling(ppe_imputed_tmp, scale = FALSE, assay = "imputed")
scaled_data <- SummarizedExperiment::assay(ppe_imputed_scaled,"scaled")

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

dim(ppe0)
dim(ppe_filtered)
dim(ppe_imputed_tmp)



#### Visualization

## hierachical cluster
plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), panel = "dendrogram", grps=grps, labels = colnames(ppe_imputed_scaled)) +
  ggplot2::ggtitle("Before batch correction")


## PCA: grps redefined for visualization
p1 <-plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), grps=grps,
            labels = colnames(ppe_imputed_scaled), panel = "pca") +
  ggplot2::ggtitle("Before batch correction")

p2 <-plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), grps=grps2, 
            labels = colnames(ppe_imputed_scaled), panel = "pca") +
  ggplot2::ggtitle("Before batch correction")


tiff(filename = paste0("../analysis/PCA/pca_donor_effect.tiff"),
     width = 14 * 300, 
     height = 10 * 300,
     res = 300,
     compression = "lzw")
ggpubr::ggarrange(p1,p2, nrow = 2, labels= c("grps","donors"), 
                  font.label = list(size = 12), label.x = c(0.3,0.3))
dev.off()



## quantification before and after impuation

p0 = plotQC(SummarizedExperiment::assay(ppe0,"Quantification"), 
            labels=colnames(ppe0), 
            panel = "quantify", grps = grps) + 
  theme_cowplot() +
  scale_fill_brewer(palette = "RdPu", direction = 1)+
  theme(legend.position = "none",title = element_blank(),
        axis.text.x = element_text(angle = 90, size =10))

p1 = plotQC(SummarizedExperiment::assay(ppe_filtered,"Quantification"), 
            labels=colnames(ppe_filtered), 
            panel = "quantify", grps = grps) + theme_cowplot() +
  scale_fill_brewer(palette = "RdPu", direction = 1)+
  theme(legend.position = "none",title = element_blank(),
        axis.text.x = element_text(angle = 90, size = 10))

p2 = plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), 
            labels=colnames(ppe_imputed_scaled), 
            panel = "quantify", grps = grps) + theme_cowplot() +
  scale_fill_brewer(palette = "RdPu", direction = 1)+
  theme(legend.position = "none", title = element_blank(),
        axis.text.x = element_text(angle = 90, size = 10))


tiff(filename = paste0("../analysis/PCA/quantification2.tiff"),
     width = 10 * 300, 
     height = 10 * 300,
     res = 300,
     compression = "lzw")

ggpubr::ggarrange(p0, p1, p2, nrow = 3, labels= c("raw(log-transformed)","filtered","imputed and scaled"), 
                  font.label = list(size = 12), label.x = c(0.5,0.5,0.5),
                  hjust = 0.5, vjust = 1)

dev.off()

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

ggpubr::ggarrange(p0, p1, p2, nrow = 3, labels= c("raw(log-transformed)","filtered","imputed and scaled"), 
                  font.label = list(size = 12), label.x = c(0.3,0.35,0.3))

tiff(filename = paste0("../analysis/PCA/dendrogram.tiff"),
     width = 14 * 300, 
     height = 10 * 300,
     res = 300,
     compression = "lzw")
ggpubr::ggarrange(p0, p1, p2, nrow = 3, labels= c("raw(log-transformed)","filtered","imputed and scaled"), 
                  font.label = list(size = 12), label.x = c(0.3,0.35,0.3))
dev.off()

p1 = plotQC(SummarizedExperiment::assay(ppe_filtered,"Quantification"), 
            labels=colnames(ppe_filtered), panel = "dendrogram", 
            grps = grps)
p2 = plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), 
            labels=colnames(ppe_imputed_scaled), 
            panel = "dendrogram", grps = grps)
ggpubr::ggarrange(p1, p2, nrow = 1)


