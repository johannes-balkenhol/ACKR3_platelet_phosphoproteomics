## script by Johannes Balkenhol & Özge Osmanoglu
## load data
## explore data
## filter data
## impute data
## median center data

# Install packages ----

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("RUVSeq")
BiocManager::install("RUV")
BiocManager::install("limma")
install.packages("remotes")

# Load packages ----
suppressPackageStartupMessages({
  library(readxl)
  library(matrixStats)
  library(data.table)
  library(dplyr)
  library("PhosR")
  library("psych")
  library("ggplot2")
  library(limma)
  library(plyr)
  library(ggplot2)
  library(corrplot) 
  library(MatrixGenerics)
})



# DATA EXPLORATION ----

## Load and prep data ----
setwd("/Users/ozgeosmanoglu/Nextcloud/CXCR7_platelet_analysis/global_proteomics_CXCR7/A08_1/scripts")

raw_data <- read.csv("../data/raw_data/TR240_A08_01_Global_1stAnalysis_Repl_1_2_3_4_5_6_7_8_forJB.csv",
                    sep = ";", dec = ".")

rownames(raw_data) <- raw_data[,2]
colnames(raw_data) <-  gsub("Abundance_","", colnames(raw_data))
peptide_info <- raw_data[c(1:16, 81:82)]
raw_abundance <- raw_data[17:80]

## Define sample names ----
donor_nr <- gsub("Donor", "",sapply(strsplit(colnames(raw_abundance), "_"), "[[", 2))
time_point <- sapply(strsplit(colnames(raw_abundance), "_"), "[[", 1)
timepoint_fac <- as.factor(time_point)
condition <- gsub("0000", "ctrl",timepoint_fac)
condition <- ifelse(grepl("ctrl", condition), "ctrl", "CXCR7")
colnames(raw_abundance) <- paste(time_point, condition, donor_nr, sep = "_")

## Check sample sizes ----
num_proteins <- colSums(!is.na(raw_abundance))
summary(num_proteins)

# Make boxplot
tiff(filename = paste0("../analysis/PCA/protein_numbers_raw.tiff"),
     width = 4* 300, 
     height = 8 * 300,
     res = 300,
     compression = "lzw")
par(mfrow = c(3,1))
boxplot(num_proteins, horizontal = TRUE,
        xlab = "Number of Proteins Detected")
plot(num_proteins, pch = 19, xlab = "Sample Index", 
     ylab = "Number of Proteins Detected")
plot(density(num_proteins, na.rm = TRUE), 
     xlim = range(num_proteins, na.rm = TRUE),
     main = "Density Plot", xlab = "Number of Proteins Detected")
dev.off()

tiff(filename = paste0("../analysis/PCA/Boxplot_Raw.tiff"),
     width = 20* 300, 
     height = 5 * 300,
     res = 300,
     compression = "lzw")
boxplot(raw_abundance, main = "Protein Abundances in Raw Data", col = timepoint_fac)
dev.off()


## Check and remove contamination ----

### Blood contamination ----
blood_prots <- c("P02088","P01942","Q8VCM7","E9PV24","Q8K0E8","P07724","Q62261",
                 "P16546","P15508","P08032") #hemoglobin, albumin, fibrinogen, spectrin


idx <- raw_data$Description %like% "Hemoglobin" # indexing matches
hemoglobin <- raw_data[idx, ]
idx <- raw_data$Description %like% "Fibrinogen"  # indexing matches
fibrinogen <- raw_data[idx, ]
idx <- raw_data$Description %like% "albumin"  # indexing matches
albumin <- raw_data[idx, ]
idx <- raw_data$Description %like% "Spectrin"  # indexing matches
spectrin <- raw_data[idx, ]
blood_prots <- rbind(hemoglobin,fibrinogen,albumin,spectrin)

# Plot blood contamination 
tiff(filename = paste0("../analysis/PCA/raw_abundance_blood_contamination.tiff"),
     width = 10* 300, 
     height = 5 * 300,
     res = 300,
     compression = "lzw")
boxplot(blood_prots[17:80], main = "Blood contamination", col = timepoint_fac)
dev.off()

### Keratin contamination ----
idx <- raw_data$Description %like% "Keratin" # indexing matches
keratins <- raw_data[idx, ]

# Plot keratin contamination
tiff(filename = paste0("../analysis/PCA/raw_abundance_keratin_contamination.tiff"),
     width = 10* 300, 
     height = 5 * 300,
     res = 300,
     compression = "lzw")
boxplot(raw_data[17:80][idx, ], main = "Keratin contamination",col = timepoint_fac)
dev.off()

keratin <- colMeans(raw_data[17:80][idx, ], na.rm = TRUE) # sample means
plot(density(keratin, na.rm = TRUE), 
     xlim = range(keratin, na.rm = TRUE),
     main = "Keratin", xlab = "Average Protein Abundance")

### Immunoglobulin contamination ----
idx <- raw_data$Description %like% "Immunoglobulin" # indexing matches
immunoglobulins1 <- raw_data[idx, ]
idx <- raw_data$Description %like% "Ig"  # indexing matches
immunoglobulins2 <- raw_data[idx, ]
immunoglobulins1 <-  immunoglobulins1[!rownames(immunoglobulins1) %in% rownames(immunoglobulins2),]
immunoglobulins <- rbind(immunoglobulins1,immunoglobulins2)

# Plot Ig contamination
tiff(filename = paste0("../analysis/PCA/raw_abundance_Ig_contamination.tiff"),
     width = 10* 300, 
     height = 7 * 300,
     res = 300,
     compression = "lzw")
boxplot(immunoglobulins[17:80], main = "Immunoglobulin contamination", col = timepoint_fac)
dev.off()

### Remove contamination (OPTIONAL) ----
#HMM DO WE REMOVE?
contamination <- rbind(blood_prots, keratins, immunoglobulins)
#raw_abundance <- raw_abundance[!rownames(raw_abundance) %in% rownames(contamination),]
# tiff(filename = paste0("../analysis/PCA/Boxplot_rawContamiRemoved.tiff"),
#      width = 20* 300, 
#      height = 5 * 300,
#      res = 300,
#      compression = "lzw")
# boxplot(raw_abundance, main = "Contamination removed", col = timepoint_fac)
# dev.off()



## Check for all NA samples ----

all_miss <- apply(raw_abundance, 2, function(x) all(is.na(x)))

# display columns with all missing values and remove if any columns is all NA
na_columns <- names(all_miss[all_miss>0])   
na_col_nums <- c(which(colnames(raw_abundance) == na_columns[1]), which(colnames(raw_abundance) == na_columns[2]))
raw_abundance2 <- raw_abundance[-na_col_nums]

# remove samples with num_proteins a lot of NAs : 1800_CXCR7_1
raw_abundance2 <- raw_abundance2[-which(colnames(raw_abundance2) == "1800_CXCR7_1")]

# check if any duplications #should be true if no duplications
length(unique(rownames(raw_abundance2))) == nrow(raw_abundance2)

## Define new sample names ----
donor_nr <- sapply(strsplit(colnames(raw_abundance2), "_"), "[[", 3)
time_point <- sapply(strsplit(colnames(raw_abundance2), "_"), "[[", 1)
timepoint_fac <- as.factor(time_point)
condition <- sapply(strsplit(colnames(raw_abundance2), "_"), "[[", 2)

## Prep metadata ----
sample_info <- data.frame(donor_nr, timepoint_fac, condition, colnames(raw_abundance2))
colnames(sample_info) <- c("donor", "time", "condition", "sample")
rownames(sample_info) <- colnames(raw_abundance2)

# DATA FILTERING ----

## filter NAs ----
### Check how many NAs ----
data_na_content <- apply(X = is.na(raw_abundance2), MARGIN = 1, FUN = sum)
names(data_na_content) <- rownames(raw_abundance2)
data_na_content <- data.frame(data_na_content)
data_na_content$name <- rownames(data_na_content)
 
### NA plot ----
tiff("../analysis/PCA/hist_of_NAs_inraw.tiff", width = 8, height = 6, units = 'in', res = 600)
a <- hist(data_na_content$data_na_content, 
          main = "Raw Data Initial Check",
          xlab = "number of NAs",
          col = "#EA8331", 
          breaks = 20) 
text(a$mids,a$counts,labels=a$counts, adj=c(0.5, -0.2))
dev.off()

### Remove lines with all NAs ----
raw_abundance2_filt<-raw_abundance2[!rownames(raw_abundance2) %in% data_na_content[data_na_content$data_na_content == ncol(raw_abundance2),"name"],]

### Sort columns ----
raw_abundance3 <- raw_abundance2_filt %>% dplyr::select(sort(names(raw_abundance2_filt))) ## sorts the columns

### Log2-transformation ----
raw_abundance3_log <- log2(raw_abundance3)

# write
write.table(raw_abundance3_log, "../data/processed_data/rawAbundance3Log.txt", sep = "\t")

## Plot distribution ----
test2 <- as.matrix(raw_abundance3)
test3 <- as.matrix(raw_abundance3_log)

dd <- reshape2::melt(test2, variable.name = "sample")
dd_log <- reshape2::melt(test3, variable.name = "sample")

ggplot(dd, aes(value, colour = Var2)) + geom_density(bw = "sj") +
  theme(legend.position="none")
#+ xlim(0,1e+7)

### Density plot ----

tiff("../analysis/PCA/densityPlotRA3Log2.tiff", 
     width = 6, height = 4, units = 'in', 
     res = 600, compression = "lzw")
ggplot(dd_log, aes(value, colour = Var2)) + geom_density(bw = "sj") +
  theme(legend.position="none") 
#+ xlim(0,1e+7)
dev.off()

### Violin Plot ----
#data <- as.data.frame(norm_intensity)
#data$ProteinID <- rownames(data)

data <- raw_abundance3_log
data$ProteinID <- rownames(raw_abundance3_log)

# Melting the data
# Adding protein IDs as a column if not already present
# Melting the data with ProteinID as the id variable
long_data <- reshape2::melt(data, id.vars = "ProteinID", variable.name = "Group", value.name = "Intensity")
#long_data <- melt(data, variable.name = "Group", value.name = "Intensity")
long_data$timepoint <- sapply(strsplit(as.character(long_data$Group), "_"), "[[", 1)

tiff("../analysis/PCA/violin_plot_raw.tiff", 
     width = 10, height = 3, units = 'in', 
     res = 600, compression = "lzw")
ggplot(long_data, aes(x = Group, y = Intensity, fill = timepoint)) +
  geom_violin(trim = FALSE) +
  #geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none")+
  labs(title = "Raw data", x = "Group", y = "Intensity")
dev.off()

## Get summary statistics  ----
desc_data<- describe(raw_abundance3)
c(min(desc_data$skew), max(desc_data$skew))  #skewed to right, skewness all positive [21.75514 32.45829]
c(min(desc_data$kurtosis), max(desc_data$kurtosis)) #heavy tails [577.0445 1279.6645]

desc_logdata <- describe(raw_abundance3_log)
c(min(desc_logdata$skew), max(desc_logdata$skew))  #almost no skewness [0.3406330 0.4650117]
c(min(desc_logdata$kurtosis), max(desc_logdata$kurtosis)) #almost no kurtosis [-0.09287762  0.09680030]

## Define the groups ----
grps = timepoint_fac
grps

grps2 <- as.factor(donor_nr)
grps2

## Select by treatment groups (replicate block) ----
data <- raw_abundance3_log
ppe <- PhosphoExperiment(assays = list(Quantification = as.matrix(data)), 
                         UniprotID = rownames(data))
ppe_filtered <- selectGrps(ppe, grps, 0.35, n=8) #similar methods, 

dim(ppe_filtered)
dim(ppe)

## Impute missing values ----
set.seed(123)
ppe_imputed <- scImpute(ppe_filtered, 0.35, grps)

imputed_filtered_intensity <- as.data.frame(ppe_imputed@assays@data@listData[["imputed"]])

## remove remaining NAs ----
abundance_no.na <- na.omit(imputed_filtered_intensity)
ppe_imputed_omit <- PhosphoExperiment(assays = list(imputed = as.matrix(abundance_no.na)))
ppe_imputed_tmp <- ppe_imputed_omit

imputed_data <- as.data.frame(ppe_imputed_tmp@assays@data@listData[["imputed"]])
colnames(imputed_data) <- colnames(imputed_data)


# DATA NORMALIZATION (median centering) ----
ppe_imputed_scaled <- medianScaling(ppe_imputed_tmp, scale = FALSE, assay = "imputed")

filtered2 <- 2^SummarizedExperiment::assay(ppe_imputed_scaled, "scaled")
scaled_data <-as.data.frame(SummarizedExperiment::assay(ppe_imputed_scaled, "scaled"))

write.table(scaled_data, "../data/processed_data/scaled_prot_data.txt", row.names = T, sep = "\t")

## Get stats and summary statistics ----
test0 <- SummarizedExperiment::assay(ppe,"Quantification")
test <- SummarizedExperiment::assay(ppe_filtered,"Quantification")
test1 <- SummarizedExperiment::assay(ppe_imputed,"imputed")
test2 <- SummarizedExperiment::assay(ppe_imputed_scaled,"scaled")

sum(is.na(test0))
sum(is.na(test))
sum(is.na(test1))
sum(is.na(test2))

dim(ppe)
dim(ppe_filtered)
dim(ppe_imputed_tmp)
dim(ppe_imputed_scaled)

desc_scaleddata <- describe(scaled_data, na.rm = TRUE)
c(min(desc_scaleddata$skew), max(desc_scaleddata$skew))  #skewed to left[0.3689168 0.5171845]
c(min(desc_scaleddata$kurtosis), max(desc_scaleddata$kurtosis)) #almost no kurtosis [-0.05348012  0.11985387]

## Correlation in scaled data ----
CorrFunction <- function(pearsonInput, digitN) {
  #Calculate Pearson CorrelationMatrix on log2 intensities
  nSamples <- ncol(pearsonInput)
  pearson <- matrix(ncol=nSamples,nrow=nSamples)
  for(i in 1:ncol(pearson))	{
    for(j in 1:nrow(pearson))	{
      pearson[i,j] <- cor(pearsonInput[,i],pearsonInput[,j],use = "pairwise.complete.obs",method="pearson")
      print(paste("i=",i,"j=", j))
      flush.console()
    }
  }
  
  cols1 <-  colorRampPalette(c("yellow","firebrick1"))(200)
  cp <- corrplot(pearson,is.corr=F,col=cols1,tl.col="black",
                 number.digits = digitN,
                 diag=T,addCoef.col="black",method="color",number.cex=0.5)
  return(cp)
}

DataList <- list(Raw = raw_abundance3_log, 
                 Scaled = scaled_data)
### All Samples ----
for (i in 1:length(DataList)) {
  tiff(paste0("../analysis/PCA/Corr", names(DataList)[i], ".tiff"), 
       width = 12, height = 10, units = 'in', res = 600, compression = "lzw")
  CorrFunction(DataList[[i]], digitN = 2)
  dev.off()
}
### Medians ----
for (i in 1:length(DataList)) {
  Input <- as.matrix(DataList[[i]])
  
  t00 <- rowMedians(Input[,1:8])
  t10 <- rowMedians(Input[,9:15])
  t30 <- rowMedians(Input[,16:23])
  t60 <- rowMedians(Input[,24:31])
  t300 <- rowMedians(Input[,32:39])
  t600 <- rowMedians(Input[,40:47])
  t900 <- rowMedians(Input[,48:55])
  t1800 <- rowMedians(Input[,56:61])
  
  Inputmedian <- cbind(t00,t10,t30, t60, t300, t600,t900, t1800)
  rownames(Inputmedian) <- rownames(input)
  colnames(Inputmedian) <- c("X0000","X0010","X0030", "X0060","X0300",
                             "X0600","X0900", "X1800")
  #plot
  tiff(paste0("../analysis/PCA/Corr", names(DataList)[i], "Median.tiff"), 
       width = 4, height = 3, units = 'in', res = 600, compression = "lzw")
  CorrFunction(Inputmedian, digitN = 3)
  dev.off()
}


## Visualization of scaled data ----

#### determine groups
vect1 <- c(8,7,8,8,8,8,8,6)
x <- grps

### Check no of proteins ----
num_proteins <- colSums(!is.na(scaled_data))
summary(num_proteins)
boxplot(num_proteins, horizontal = TRUE,
        xlab = "Number of Proteins Detected")

plot(num_proteins, pch = 19, xlab = "Sample Index", 
     ylab = "Number of Proteins Detected")

plot(density(num_proteins, na.rm = TRUE), 
     xlim = range(num_proteins, na.rm = TRUE),
     main = "Density Plot", xlab = "Number of Proteins Detected")

tiff(filename = paste0("../analysis/PCA/scaled_data_boxplot.tiff"),
     width = 20* 300, 
     height = 5 * 300,
     res = 300,
     compression = "lzw")
boxplot(norm_intensity, main = "After Scaling ", col = timepoint_fac)
dev.off()

### Violin Plot ----
#data <- as.data.frame(norm_intensity)
#data$ProteinID <- rownames(data)

data <- scaled_data
data$ProteinID <- rownames(scaled_data)

# Melting the data
# Adding protein IDs as a column if not already present
# Melting the data with ProteinID as the id variable
long_data <- reshape2::melt(data, id.vars = "ProteinID", variable.name = "Group", value.name = "Intensity")
#long_data <- melt(data, variable.name = "Group", value.name = "Intensity")
long_data$timepoint <- sapply(strsplit(as.character(long_data$Group), "_"), "[[", 1)

tiff("../analysis/PCA/violin_plot_scaled.tiff", 
     width = 10, height = 3, units = 'in', 
     res = 600, compression = "lzw")
ggplot(long_data, aes(x = Group, y = Intensity, fill = timepoint)) +
  geom_violin(trim = FALSE) +
  #geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none")+
  labs(title = "Scaled data", x = "Group", y = "Intensity")
dev.off()

# VISUALIZATION ----
inputList <- list(raw =  as.matrix(raw_abundance2),     #NAs filtered out
                  logT = as.matrix(raw_abundance3_log), #log_transformed
                  Imp = as.matrix(imputed_data),        #imputed&scaled
                  ImpSca = as.matrix(scaled_data),      #imputed
                  RUV = norm_intensity)                 #ruv
                  

## Density Plots ----
##change input!
i = 4
dd2 <- reshape2::melt(inputList[[i]], variable.name = "sample")

vect2 <- c( nrow(dd2)/length(x)*vect1[1],  nrow(dd2)/length(x)*vect1[2],  
            nrow(dd2)/length(x)*vect1[3],  nrow(dd2)/length(x)*vect1[4],
            nrow(dd2)/length(x)*vect1[5],  nrow(dd2)/length(x)*vect1[6],
            nrow(dd2)/length(x)*vect1[7],  nrow(dd2)/length(x)*vect1[8]) 
dd2$Var3 <- as.factor(c(rep(c("0000", "0010", "0030", "0060", "0300", "0600", "0900", "1800"), vect2)))
colnames(dd2) <- c("ID", "samples", "value", "condition")

#another way to separate also the legends
dd2_ctrl <- dd2[dd2$samples %in% sample_info[sample_info$time == "0000",]$sample,]
dd2_10 <- dd2[dd2$samples %in% sample_info[sample_info$time == "0010",]$sample,]
dd2_30<- dd2[dd2$samples %in% sample_info[sample_info$time == "0030",]$sample,]
dd2_60 <- dd2[dd2$samples %in% sample_info[sample_info$time == "0060",]$sample,]
dd2_300 <- dd2[dd2$samples %in% sample_info[sample_info$time == "0300",]$sample,]
dd2_600<- dd2[dd2$samples %in% sample_info[sample_info$time == "0600",]$sample,]
dd2_900 <- dd2[dd2$samples %in% sample_info[sample_info$time == "0900",]$sample,]
dd2_1800 <- dd2[dd2$samples %in% sample_info[sample_info$time == "1800",]$sample,]


p0 <-  ggplot(dd2, aes(value, colour = samples)) + geom_density(bw = "sj") +
  theme(legend.position="right", legend.title = element_text(color = "white"))
p1 <- ggplot(dd2_ctrl, aes(value, colour = samples)) + geom_density(bw = "sj") +
  theme(legend.position="right", legend.title = element_text(color = "white"))
  #scale_color_manual(values = c("#F8766D", "#EA8331", "#D89000", "#C09B00", "#A3A500"))
p2 <-  ggplot(dd2_10, aes(value, colour = samples)) + geom_density(bw = "sj") +
  theme(legend.position="right", legend.title = element_text(color = "white"))
  #scale_color_manual(values = c("#7CAE00", "#39B600", "#00BB4E", "#00BF7D", "#00C1A3"))
p3 <-  ggplot(dd2_30, aes(value, colour = samples)) + geom_density(bw = "sj") +
  theme(legend.position="right", legend.title = element_text(color = "white"))
#scale_color_manual(values = c("#7CAE00", "#39B600", "#00BB4E", "#00BF7D", "#00C1A3"))
p4 <-  ggplot(dd2_60, aes(value, colour = samples)) + geom_density(bw = "sj") +
  theme(legend.position="right", legend.title = element_text(color = "white"))
#scale_color_manual(values = c("#7CAE00", "#39B600", "#00BB4E", "#00BF7D", "#00C1A3"))
p5 <-  ggplot(dd2_300, aes(value, colour = samples)) + geom_density(bw = "sj") +
  theme(legend.position="right", legend.title = element_text(color = "white"))
#scale_color_manual(values = c("#7CAE00", "#39B600", "#00BB4E", "#00BF7D", "#00C1A3"))
p6 <-  ggplot(dd2_600, aes(value, colour = samples)) + geom_density(bw = "sj") +
  theme(legend.position="right", legend.title = element_text(color = "white"))
#scale_color_manual(values = c("#7CAE00", "#39B600", "#00BB4E", "#00BF7D", "#00C1A3"))
p7 <-  ggplot(dd2_900, aes(value, colour = samples)) + geom_density(bw = "sj") +
  theme(legend.position="right", legend.title = element_text(color = "white"))
#scale_color_manual(values = c("#7CAE00", "#39B600", "#00BB4E", "#00BF7D", "#00C1A3"))
p8 <-  ggplot(dd2_1800, aes(value, colour = samples)) + geom_density(bw = "sj") +
  theme(legend.position="right", legend.title = element_text(color = "white"))
#scale_color_manual(values = c("#7CAE00", "#39B600", "#00BB4E", "#00BF7D", "#00C1A3"))


sep<-ggpubr::ggarrange(p1, p2, p3, p4,p5,p6,p7,p8, nrow = 2, ncol=4)

tiff(filename = paste0("../analysis/PCA/density_plots_", names(inputList)[i], ".tiff"),
     width = 36* 300, 
     height = 10 * 300,
     res = 300,
     compression = "lzw")
ggpubr::ggarrange(p0, sep, nrow = 1, ncol=2)
dev.off()



## Quantification ----

ppe0 <- PhosphoExperiment(assays = list(Quantification = as.matrix(raw_abundance)), 
                         UniprotID = rownames(raw_abundance))
ppe_filtered <- PhosphoExperiment(assays = list(Quantification = as.matrix(raw_abundance2_filt)), 
                                  UniprotID = rownames(raw_abundance2_filt))
vect <- c(8,8,8,8,8,8,8,8)
x_init <- as.factor(c(rep(c("0000", "0010", "0030", "0060", "0300", "0600", "0900", "1800"), vect)))


p0 = plotQC(SummarizedExperiment::assay(ppe0,"Quantification"), 
            labels=colnames(ppe0), 
            panel = "quantify", grps = x_init) + labs(tag = "A")
p1 = plotQC(SummarizedExperiment::assay(ppe_filtered,"Quantification"), 
            labels=colnames(ppe_filtered), 
            panel = "quantify", grps = x)+ labs(tag = "B")
p2 = plotQC(SummarizedExperiment::assay(ppe_imputed_omit,"imputed"), 
            labels=colnames(ppe_imputed_omit), panel = "quantify", grps = x)+ labs(tag = "C")
p3 = plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), 
            labels=colnames(ppe_imputed_scaled), panel = "quantify", grps = x)+ labs(tag = "D")

tiff(filename = paste0("../analysis/PCA/quantification.tiff"),
     width = 6 * 300, 
     height = 12 * 300,
     res = 300,
     compression = "lzw")

ggpubr::ggarrange(p0, p1, p2,p3, nrow = 4, labels= c("raw","filtered","imputed","imputed and scaled"), 
                  font.label = list(size = 12), label.x = c(0.3,0.25,0.25, 0.15))

dev.off()

## Hierachical clustering ----

grps_temp <- grps

p0 = plotQC(as.matrix(raw_abundance3_log), 
            labels=colnames(raw_abundance3_log), panel = "dendrogram", 
            grps = grps_temp)+ labs(tag = "A", title = "") 
p1 = plotQC(SummarizedExperiment::assay(ppe_filtered,"Quantification"), 
            labels=colnames(ppe_filtered), panel = "dendrogram", 
            grps = grps_temp)+ labs(tag = "B", title = "")
p2 = plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"imputed"), 
            labels=colnames(ppe_imputed_scaled),
            panel = "dendrogram", grps = grps_temp)+ labs(tag = "C",title = "")
p3 = plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), 
            labels=colnames(ppe_imputed_scaled), 
            panel = "dendrogram", grps = grps_temp)+ labs(tag = "C", title = "")
p4 = plotQC(norm_intensity, 
            labels = colnames(norm_intensity), 
            panel = "dendrogram", grps = grps_temp)+ labs(tag = "D", title = "")


tiff(filename = paste0("../analysis/PCA/dendrograms.tiff"),
     width = 12* 300, 
     height = 10 * 300,
     res = 300,
     compression = "lzw")
ggpubr::ggarrange(p0, p2,p3,p4, nrow = 4, labels= c("raw(log-transformed)","imputed", "imputed and scaled", "ruvphospho"), 
                  font.label = list(size = 12), label.x = c(0.3,0.35,0.3))
dev.off()


## PCA ----

### Alltogether ----
grps_temp <- grps

p0 <- plotQC(as.matrix(raw_abundance3_log),grps=grps_temp, labels = "", panel="pca") +
  #scale_color_manual(values=c( "#A6CEE3","#1F78B4","#085304", "#33A02C")) + 
  #geom_text(aes(label =colnames(as.matrix(raw_abundance3_log))), size = 3, color = "black")+
  geom_point(size=5,alpha=0.5) + theme(legend.position="right")
p1 <-  plotQC(SummarizedExperiment::assay(ppe_filtered,"Quantification"),grps=grps_temp, labels = "", panel="pca") +
  #scale_color_manual(values=c( "#A6CEE3","#1F78B4","#085304", "#33A02C")) + 
  #geom_text(aes(label =colnames(ppe_filtered)), size = 3, color = "black")+
  geom_point(size=5,alpha=0.5) + theme(legend.position="right")
p2 <- plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"imputed"),grps=grps_temp, labels = "", panel="pca") +
  #scale_color_manual(values=c( "#A6CEE3","#1F78B4","#085304", "#33A02C")) + 
  #geom_text(aes(label =colnames(ppe_imputed_scaled)), size = 3, color = "black")+
  geom_point(size=5,alpha=0.5) + theme(legend.position="right")
p3 <- plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"),grps=grps_temp, labels = "", panel="pca") +
  #scale_color_manual(values=c( "#A6CEE3","#1F78B4","#085304", "#33A02C")) + 
  #geom_text(aes(label =colnames(ppe_imputed_scaled)), size = 3, color = "black")+
  geom_point(size=5,alpha=0.5) + theme(legend.position="right")
p4 <- plotQC(norm_intensity,grps=grps_temp, labels = "", panel="pca") +
  #scale_color_manual(values=c( "#A6CEE3","#1F78B4","#085304", "#33A02C")) + 
  #geom_text(aes(label =colnames(norm_intensity)), size = 3, color = "black")+
  geom_point(size=5,alpha=0.5) + theme(legend.position="right")

tiff(filename = paste0("../analysis/PCA/pcas.tiff"),
     width = 15* 300, 
     height = 10 * 300,
     res = 300,
     compression = "lzw")
ggpubr::ggarrange(p0, p2,p3,p4, nrow = 2, ncol = 2, labels= c("raw(log-transformed)","imputed", "imputed and scaled", "ruvphospho"), 
                  font.label = list(size = 12), label.x = c(0.3,0.35,0.3))
dev.off()


### Individual ----
i = 4

tiff(filename = paste0("../analysis/PCA/PCA", names(inputList)[i], "tiff"),
     width = 7 * 300, 
     height = 4  * 300,
     res = 300,
     compression = "lzw")
plotQC(inputList[[i]],grps=grps, labels = "", panel="pca") +
  #scale_color_manual(values=c( "#A6CEE3","#1F78B4","#085304", "#33A02C")) + 
  geom_text(aes(label =colnames(inputList[[i]])), size = 3, color = "black")+
  geom_point(size=5,alpha=0.5) + theme(legend.position="right")+ 
  labs(title = names(inputList)[i])
dev.off()


