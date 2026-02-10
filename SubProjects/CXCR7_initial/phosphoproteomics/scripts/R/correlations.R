library(ggvenn)
library(RColorBrewer)
library(corrplot)
library(matrixStats)
library(pheatmap)


# COMPARE EMPIRICALS ----
empirical_topall_old <- read.table("../../../CXCR7_initial/phosphoproteomics/data/processed_data/empirical_topall.txt", 
                                   sep="\t", header=TRUE, dec=",")


empirical.old <- paste(sapply(strsplit(empirical_topall_old$x, ";"), "[[", 2),
                       sapply(strsplit(empirical_topall_old$x, ";"), "[[", 3), sep = ";")

#empirical.new <- paste(sapply(strsplit(empirical_topall, ";"), "[[", 2),
#                       sapply(strsplit(empirical_topall, ";"), "[[", 3), sep = ";")

sel <- which(empirical.new %in% empirical.old)
empirical_match <- empirical_topall[sel]

sel <- which(empirical.old %in% empirical.new)
empirical_topall_old <- as.data.frame(empirical_topall_old)
empirical_match <- empirical_topall_old[sel,]

write.table(empirical_match, "data/processed_data/empirical_match.txt", sep="\t", , row.names=FALSE)


# COMPARE PHOSPHO LOGFCs ----

## load the data and change rownames to prot;site ----
top.all_val<- read.table("../../../CXCR7_validation/phosphoproteomics/data/processed_data/top.all.cxcr7.vs.0s.txt", 
                         sep="\t", header=TRUE, dec=".")
top.all_old <- top.all

rownames(top.all_val) <- paste(sapply(strsplit(rownames(top.all_val), ";"),"[[", 1),
                               sapply(strsplit(rownames(top.all_val), ";"),"[[", 3), sep = ";")

rownames(top.all_old) <- paste(sapply(strsplit(rownames(top.all_old), ";"),"[[", 1),
                               sapply(strsplit(rownames(top.all_old), ";"),"[[", 3), sep = ";")

#3323 dim(top.all_old)
#3150 dim(top.all_val)
#1889 overlap

## check overlaps ----

tiff("../analysis/PCA/venn_oldnew.tiff",
     width = 5 * 300, 
     height = 5 * 300,
     res = 300,
     compression = "lzw")
ggvenn(list("validation data" = rownames(top.all_val),
            "initial data" = rownames(top.all_old)), 
       fill_color = c("white", "white"), set_name_size = 4)
dev.off()


mergedtops <- merge(top.all_old, top.all_val, by = "row.names", all = T)

#if you only want to use common phosphopeptides, but it does not really change the correlations.
mergedtops <- merge(top.all_old, top.all_val, by = "row.names", all = F)

colnames(mergedtops) = gsub(pattern = "\\.x", ".old", colnames(mergedtops))
colnames(mergedtops) = gsub(pattern = "\\.y", ".val", colnames(mergedtops))

#plot(mergedtops$logFC.10.old, mergedtops$logFC.10.val)      


#input <- mergedtops[mergedtops$Row.names %in% intersect(rownames(top.all_val), rownames(top.all_old)),]
input <- mergedtops[c(4, 16, 22, 27, 30,33)]

## make corr plot ----
correlation_matrix <- cor(input, use = "pairwise.complete.obs") #for logfcs
#correlation_matrix <- cor(mergedtops[c(5, 17, 23, 28, 31,34)], use = "pairwise.complete.obs") #for pvals
#correlation_matrix <- cor(mergedtops[c(6, 18, 24, 29, 32,35)], use = "pairwise.complete.obs") #for pvals

corrplot(correlation_matrix, method = "circle", 
         order = "hclust", col = COL2('BrBG'),
         addrect = 2, addCoef.col = 'black')


tiff("../analysis/PCA/correlation_oldnew_commons.tiff",
     width = 6 * 300, 
     height = 6 * 300,
     res = 300,
     compression = "lzw")
corrplot(correlation_matrix, method = "circle", 
         order = "hclust", col = COL2('BrBG'),number.digits = 3, 
         addrect = 2, addCoef.col = 'black')
dev.off()



# COMPARE PHOSPHO INTENSITIES ----

# Load data efficiently into a list to avoid redundant code
files <- list(
  norm_intensity_val = "../../../CXCR7_validation/phosphoproteomics/data/processed_data/norm_intensity.txt",
  norm_intensity_old = "../../../CXCR7_initial/phosphoproteomics/data/processed_data/norm_intensity.txt",
  raw_intensity_val  = "../../../CXCR7_validation/phosphoproteomics/data/processed_data/raw_intensities.txt",
  raw_intensity_old  = "../../../CXCR7_initial/phosphoproteomics/data/processed_data/raw_intensities.txt"
)

# Read data into matrices
intensity_data <- lapply(files, function(file) as.matrix(read.table(file, sep="\t", header=TRUE, dec=".")))

# Assign individual variables from list
norm_intensity_val <- intensity_data[["norm_intensity_val"]]
norm_intensity_old <- intensity_data[["norm_intensity_old"]]
raw_intensity_val  <- intensity_data[["raw_intensity_val"]]
raw_intensity_old  <- intensity_data[["raw_intensity_old"]]

# Process row names for both datasets using a function to avoid duplication
process_row_names <- function(data) {
  rownames(data) <- paste(sapply(strsplit(rownames(data), ";"), "[[", 1),
                          sapply(strsplit(rownames(data), ";"), "[[", 3), sep = ";")
  data
}

inputOldNorm <- process_row_names(norm_intensity_old)
inputValNorm <- process_row_names(norm_intensity_val)
inputOldRaw <- process_row_names(raw_intensity_old)
inputValRaw <- process_row_names(raw_intensity_val)

## get medians ----
# Define a function to compute medians for specified column ranges
compute_medians <- function(data, ranges) {
  sapply(ranges, function(cols) rowMedians(data[, cols]))
}

# Column ranges and names for old data (both raw and normalized)
old_ranges <- list(1:7, 8:13, 14:20, 21:27, 28:34, 35:40, 41:47, 48:53)  # Adjust as needed for desired timepoints
old_colnames <- c("X0000", "X0010", "X0030", "X0060", 
                  "X0300", "X0600", "X0900",  "X1800")

phospho_norm_old_median <- compute_medians(inputOldNorm, old_ranges)
rownames(phospho_norm_old_median) <- rownames(inputOldNorm)
colnames(phospho_norm_old_median) <- old_colnames

phospho_raw_old_median <- compute_medians(inputOldRaw, old_ranges)
rownames(phospho_raw_old_median) <- rownames(inputOldRaw)
colnames(phospho_raw_old_median) <- old_colnames

# Column ranges and names for validation data (both raw and normalized)
val_ranges <- list(1:10, 11:20, 31:40, 51:60, 21:30, 41:50, 61:70)
val_colnames <- c("X0000v", "X0010v", "X0600v", "X1800v", "X0010vDMSO", "X0600vDMSO", "X1800vDMSO")

phospho_norm_val_median <- compute_medians(inputValNorm, val_ranges)
rownames(phospho_norm_val_median) <- rownames(inputValNorm)
colnames(phospho_norm_val_median) <- val_colnames

phospho_raw_val_median <- compute_medians(inputValRaw, val_ranges)
rownames(phospho_raw_val_median) <- rownames(inputValRaw)
colnames(phospho_raw_val_median) <- val_colnames

## get medians for combined figure ----

# decide on input
inputOld <- raw_intensity_old # change input accordingly
inputVal <-raw_intensity_val


rownames(inputVal) <- paste(sapply(strsplit(rownames(inputVal), ";"),"[[", 1),
                            sapply(strsplit(rownames(inputVal), ";"),"[[", 3), sep = ";")

rownames(inputOld) <- paste(sapply(strsplit(rownames(inputOld), ";"),"[[", 1),
                            sapply(strsplit(rownames(inputOld), ";"),"[[", 3), sep = ";")

### get medians for old data ----

t00 <- rowMedians(inputOld[,1:7])
t10 <- rowMedians(inputOld[,8:13])
t30 <- rowMedians(inputOld[,14:20])
t60 <- rowMedians(inputOld[,21:27])
t300 <- rowMedians(inputOld[,28:34])
t600 <- rowMedians(inputOld[,35:40])
t900 <- rowMedians(inputOld[,41:47])
t1800 <- rowMedians(inputOld[,48:53])

inputOld_median <- cbind(t00,t10,t600,t1800)
rownames(inputOld_median) <- rownames(inputOld)
colnames(inputOld_median) <- c("X0000","X0010","X0600","X1800")

### get medians for val data ----

t00 <- rowMedians(inputVal[,1:10])
t10 <- rowMedians(inputVal[,11:20])
t10wt <- rowMedians(inputVal[,21:30])
t600 <- rowMedians(inputVal[,31:40])
t600wt <- rowMedians(inputVal[,41:50])
t1800 <- rowMedians(inputVal[,51:60])
t1800wt <- rowMedians(inputVal[,61:70])

inputVal_median <- cbind(t00,t10,t600,t1800,t10wt,t600wt,t1800wt)
rownames(inputVal_median) <- rownames(inputVal)
colnames(inputVal_median) <- c("X0000v",
                               "X0010v", "X0600v","X1800v",
                               "X0010vDMSO","X0600vDMSO","X1800vDMSO")

## check overlaps ----

ggvenn(list("validation data" = rownames(inputVal_median),
            "initial data" = rownames(inputOld_median)), 
       fill_color = c("white", "white"), set_name_size = 4)


mergedtopsNormInt <- merge(inputOld_median, inputVal_median, by = "row.names", all = F)

input <- mergedtopsNormInt[2:12]
rownames(input) <- mergedtopsNormInt$Row.names
colnames(input) <- c("0s", "10s", "600s", "1800s","0s.val", 
                     "10s.val","600s.val", "1800s.val",
                     "10s.DMSO.val","600s.DMSO.val","1800s.DMSO.val")

## make corr plot ----
correlation_matrix <- cor(input, use = "pairwise.complete.obs", method = "pearson") #for logfcs

tiff("../../../analysis/correlations/correlation_normMed_phospho.tiff", 
     width = 8, height = 6, units = 'in', res = 600, compression = "lzw")
corrplot(correlation_matrix, method = "number", 
         order = "original", col = COL2('RdBu'),number.digits = 2,
         addCoef.col = 'black', type = "lower")
dev.off()



## Stefans script for corr. ----

### all samples ----
mergedtopsNormInt <- merge(inputOld, inputVal, by = "row.names", all = F)
pearsonInput <-  mergedtopsNormInt[2:ncol(mergedtopsNormInt)]
rownames(pearsonInput) <- mergedtopsNormInt$Row.names

#Calculate Pearson Correlation Matrix on intensities
nsamples <- ncol(pearsonInput)
pearson <- matrix(ncol=nsamples,nrow=nsamples)
for(i in 1:ncol(pearson))	{
  for(j in 1:nrow(pearson))	{
    pearson[i,j] <- cor(pearsonInput[,i],pearsonInput[,j],use = "pairwise.complete.obs",method="pearson")
    print(paste("i=",i,"j"=j))
    flush.console()
  }
}

cols1 <-  colorRampPalette(c("yellow","firebrick1"))(200)
tiff("../../../analysis/correlations/correlation_rawAll_phospho_stefanStyle.tiff", 
     width = 25, height = 20, units = 'in', res = 600, compression = "lzw")
corrplot(pearson,is.corr=F,col=cols1,
         number.digits = 2, tl.col="black",diag=T,addCoef.col="black",method="color",number.cex=0.5)
dev.off()

### median samples ----

mergedtopsNormInt <- merge(inputOld_median, inputVal_median, by = "row.names", all = F)
pearsonInput <-  mergedtopsNormInt[2:ncol(mergedtopsNormInt)]
rownames(pearsonInput) <- mergedtopsNormInt$Row.names


#Calculate Pearson Correlation Matrix on intensities
nsamples <- ncol(pearsonInput)
pearson <- matrix(ncol=nsamples,nrow=nsamples)
for(i in 1:ncol(pearson))	{
  for(j in 1:nrow(pearson))	{
    pearson[i,j] <- cor(pearsonInput[,i],pearsonInput[,j],use = "pairwise.complete.obs",method="pearson")
    print(paste("i=",i,"j"=j))
    flush.console()
  }
}
colnames(pearson) <- c("0s", "10s", "600s", "1800s","0s.val", 
                       "10s.val","600s.val", "1800s.val",
                       "10s.DMSO.val","600s.DMSO.val","1800s.DMSO.val")
rownames(pearson) <- c("0s", "10s", "600s", "1800s","0s.val", 
                       "10s.val","600s.val", "1800s.val",
                       "10s.DMSO.val","600s.DMSO.val","1800s.DMSO.val")


tiff("../../../analysis/correlations/correlation_rawMed_phospho_stefanStyle.tiff", 
     width = 8, height = 6, units = 'in', res = 600, compression = "lzw")
corrplot(pearson,is.corr=F,col=cols1,tl.col="black",diag=T,
         number.digits = 3, addCoef.col="black",method="color",number.cex=0.8)
#mtext("Correlation of Proteome Datasets (median raw intensities)", at=7, line=3, cex=1)
dev.off()


# COMPARE PROT INTENSITIES ----
## norm intensities ----
prot_norm_old <- as.matrix(read.table("../../../CXCR7_initial/proteomics/data/processed_data/norm_intensity.txt", 
                                      sep="\t", header=TRUE, dec="."))
prot_norm_val <- as.matrix(read.table("../../../CXCR7_validation/proteomics/data/processed_data/norm_intensity.txt", 
                                      sep="\t", header=TRUE, dec="."))
## raw intensities ----
prot_raw_old <- as.matrix(read.table("../../../CXCR7_initial/proteomics/data/processed_data/rawAbundance3Log.txt",
                                     sep="\t", header=TRUE, dec="."))
prot_raw_val <-  as.matrix(read.table("../../../CXCR7_validation/proteomics/data/processed_data/rawAbundance3Log.txt",
                                      sep="\t", header=TRUE, dec="."))

# Load data efficiently into a list to avoid redundant code
files <- list(
  pnorm_intensity_val = "../../../CXCR7_validation/proteomics/data/processed_data/norm_intensity.txt",
  pnorm_intensity_old = "../../../CXCR7_initial/proteomics/data/processed_data/norm_intensity.txt",
  praw_intensity_val  = "../../../CXCR7_validation/proteomics/data/processed_data/rawAbundance3Log.txt",
  praw_intensity_old  = "../../../CXCR7_initial/proteomics/data/processed_data/rawAbundance3Log.txt"
)

# Read data into matrices
intensity_data <- lapply(files, function(file) as.matrix(read.table(file, sep="\t", header=TRUE, dec=".")))

# Assign individual variables from list
inputValNorm <- intensity_data[["pnorm_intensity_val"]]
inputOldNorm <- intensity_data[["pnorm_intensity_old"]]
inputValRaw  <- intensity_data[["praw_intensity_val"]]
inputOldRaw  <- intensity_data[["praw_intensity_old"]]

## get medians ----
# Define a function to compute medians for specified column ranges
compute_medians <- function(data, ranges) {
  sapply(ranges, function(cols) rowMedians(data[, cols]))
}

# Column ranges and names for old data (both raw and normalized)
old_ranges <- list(1:8, 9:15, 16:23, 24:31, 32:39, 40:47, 48:55, 56:61)  # Adjust as needed for desired timepoints
old_colnames <- c("X0000", "X0010", "X0030", "X0060", 
                  "X0300", "X0600", "X0900",  "X1800")


prot_norm_old_median <- compute_medians(inputOldNorm, old_ranges)
rownames(prot_norm_old_median) <- rownames(inputOldNorm)
colnames(prot_norm_old_median) <- old_colnames

prot_raw_old_median <- compute_medians(inputOldRaw, old_ranges)
rownames(prot_raw_old_median) <- rownames(inputOldRaw)
colnames(prot_raw_old_median) <- old_colnames

# Column ranges and names for validation data (both raw and normalized)
val_ranges <- list(1:10, 11:19, 20:29,  30:39, 40:49, 50:59, 60:69)
val_colnames <- c("X0000v", "X0010v", "X0600v", "X1800v", "X0010vDMSO", "X0600vDMSO", "X1800vDMSO")

prot_norm_val_median <- compute_medians(inputValNorm, val_ranges)
rownames(prot_norm_val_median) <- rownames(inputValNorm)
colnames(prot_norm_val_median) <- val_colnames

prot_raw_val_median <- compute_medians(inputValRaw, val_ranges)
rownames(prot_raw_val_median) <- rownames(inputValRaw)
colnames(prot_raw_val_median) <- val_colnames



## get medians for combined figure ----
### decide on input ----
inputOld <- prot_norm_old # change input accordingly
inputVal <- prot_norm_val

### get medians for old data ----

t00 <- rowMedians(inputOld[,1:8])
t10 <- rowMedians(inputOld[,9:15])
t30 <- rowMedians(inputOld[,16:23])
t60 <- rowMedians(inputOld[,24:31])
t300 <- rowMedians(inputOld[,32:39])
t600 <- rowMedians(inputOld[,40:47])
t900 <- rowMedians(inputOld[,48:55])
t1800 <- rowMedians(inputOld[,56:61])

inputOld_median <- cbind(t00,t10,t600,t1800)
rownames(inputOld_median) <- rownames(inputOld)
colnames(inputOld_median) <- c("X0000","X0010","X0600","X1800")

### get medians for val data ----

t00 <- rowMedians(inputVal[,1:10])
t10 <- rowMedians(inputVal[,11:19])
t10wt <- rowMedians(inputVal[,20:29])
t600 <- rowMedians(inputVal[,30:39])
t600wt <- rowMedians(inputVal[,40:49])
t1800 <- rowMedians(inputVal[,50:59])
t1800wt <- rowMedians(inputVal[,60:69])

inputVal_median <- cbind(t00,t10,t600,t1800, t10wt, t600wt, t1800wt)
rownames(inputVal_median) <- rownames(inputVal)
colnames(inputVal_median) <- c("X0000v",
                               "X0010v", "X0600v","X1800v",
                               "X0010vDMSO","X0600vDMSO","X1800vDMSO")


## check overlaps ----

ggvenn(list("validation data" = rownames(inputVal_median),
            "initial data" = rownames(inputOld_median)), 
       fill_color = c("white", "white"), set_name_size = 4)


mergedtopsNormInt <- merge(inputOld_median, inputVal_median, by = "row.names", all = F)

#plot(mergedtopsNormInt$X0000, mergedtopsNormInt$X0000v)   

input <- mergedtopsNormInt[2:12]
rownames(input) <- mergedtopsNormInt$Row.names
colnames(input) <- c("0s", "10s", "600s", "1800s","0s.val", 
                     "10s.val","600s.val", "1800s.val",
                     "10s.DMSO.val","600s.DMSO.val","1800s.DMSO.val")


## make corr plot ----
correlation_matrix <- cor(input, use = "pairwise.complete.obs", method = "pearson") #for logfcs

tiff("../../../analysis/correlations/correlation_normMed_prot.tiff", 
     width = 8, height = 6, units = 'in', res = 600, compression = "lzw")
corrplot(correlation_matrix, method = "number", 
         order = "original", col = COL2('RdBu'),
         addCoef.col = 'black', type = "lower")
dev.off()

## Stefans script for corr. ----

### all samples ----
mergedtopsNormInt <- merge(inputOld, inputVal, by = "row.names", all = F)
pearsonInput <-  mergedtopsNormInt[2:ncol(mergedtopsNormInt)]
rownames(pearsonInput) <- mergedtopsNormInt$Row.names

#Calculate Pearson Correlation Matrix on intensities
nsamples <- ncol(pearsonInput)
pearson <- matrix(ncol=nsamples,nrow=nsamples)
for(i in 1:ncol(pearson))	{
  for(j in 1:nrow(pearson))	{
    pearson[i,j] <- cor(pearsonInput[,i],pearsonInput[,j],use = "pairwise.complete.obs",method="pearson")
    print(paste("i=",i,"j"=j))
    flush.console()
  }
}

cols1 <-  colorRampPalette(c("yellow","firebrick1"))(200)
tiff("../../../analysis/correlations/correlation_rawAll_prot_stefanStyle.tiff", 
     width = 25, height = 20, units = 'in', res = 600, compression = "lzw")
corrplot(pearson,is.corr=F,col=cols1,tl.col="black",diag=T,addCoef.col="black",method="color",number.cex=0.5)
dev.off()

### median samples ----

mergedtopsNormInt <- merge(inputOld_median, inputVal_median, by = "row.names", all = F)
pearsonInput <-  mergedtopsNormInt[2:ncol(mergedtopsNormInt)]
rownames(pearsonInput) <- mergedtopsNormInt$Row.names


#Calculate Pearson Correlation Matrix on intensities
nsamples <- ncol(pearsonInput)
pearson <- matrix(ncol=nsamples,nrow=nsamples)
for(i in 1:ncol(pearson))	{
  for(j in 1:nrow(pearson))	{
    pearson[i,j] <- cor(pearsonInput[,i],pearsonInput[,j],use = "pairwise.complete.obs",method="pearson")
    print(paste("i=",i,"j"=j))
    flush.console()
  }
}
colnames(pearson) <- c("0s", "10s", "600s", "1800s","0s.val", 
                       "10s.val","600s.val", "1800s.val",
                       "10s.DMSO.val","600s.DMSO.val","1800s.DMSO.val")
rownames(pearson) <- c("0s", "10s", "600s", "1800s","0s.val", 
                       "10s.val","600s.val", "1800s.val",
                       "10s.DMSO.val","600s.DMSO.val","1800s.DMSO.val")


tiff("../../../analysis/correlations/correlation_rawMed_prot_stefanStyle.tiff", 
     width = 8, height = 6, units = 'in', res = 600, compression = "lzw")
corrplot(pearson,is.corr=F,col=cols1,tl.col="black",diag=T,
         number.digits = 3, addCoef.col="black",method="color",number.cex=0.8)
#mtext("Correlation of Proteome Datasets (median raw intensities)", at=7, line=3, cex=1)
dev.off()


# INDIVIDUALS CORRPLOTS FOR DATASETS

# ---- Correlation Plots ----
## correlation plots

OldDataList <- list(ProtRawOld = prot_raw_old, 
                    ProtRawOldMedian = prot_raw_old_median,
                    PhosphoRawOld = raw_intensity_old,
                    PhosphoRawOldMedian = phospho_raw_old_median,
                    ProtNormOld = prot_norm_old, 
                    ProtNormOldMedian = prot_norm_old_median,
                    PhosphoNormOld = norm_intensity_old,
                    PhosphoNormOldMedian = phospho_norm_old_median
)

ValDataList <- list(ProtRawVal = prot_raw_val, 
                    ProtRawValMedian = prot_raw_val_median,
                    PhosphoRawVal = raw_intensity_val,
                    PhosphoRawValMedian = phospho_raw_val_median,
                    ProtNormVal = prot_norm_val,
                    ProtNormValMedian = prot_norm_val_median,
                    PhosphoNormVal = norm_intensity_val, 
                    PhosphoNormValMedian = phospho_norm_val_median)
                        

# Stefans script for corr. ----

CorrFunction <- function(pearsonInput, description, path) {
  #Calculate Pearson CorrelationMatrix on log2 intensities
  nSamples <- ncol(pearsonInput)
  pearson <- matrix(ncol=nSamples,nrow=nSamples)
  for(i in 1:ncol(pearson))	{
    for(j in 1:nrow(pearson))	{
      pearson[i,j] <- cor(pearsonInput[,i],pearsonInput[,j],
                          use = "pairwise.complete.obs",method="pearson")
      print(paste("i=",i,"j=", j))
      flush.console()
    }
    }

  cols1 <-  colorRampPalette(c("yellow","firebrick1"))(200)
  size_factor <- ifelse(grepl("Median$", description), 0.3, 1)
  digit_factor <- ifelse(grepl("Median$", description), 1.5, 1)
  nm <-   ifelse(grepl("Median$", description), 1.67, 1)

  if (path == "initial") {
    tiff(paste0("../../../CXCR7_initial/Correlations/Corr", description, ".tiff"), 
         width = 12 * size_factor, height = 10 * size_factor, units = 'in', res = 600, compression = "lzw")
    corrplot(pearson, is.corr=FALSE, col=cols1, tl.col="black", diag=TRUE, addCoef.col="black", method="color", number.cex=0.3 * nm, number.digits = 2 * digit_factor)
    dev.off()
  } else if (path == "validation") {
    tiff(paste0("../../../CXCR7_validation/Correlations/Corr", description, ".tiff"), 
         width = 12 * size_factor, height = 10 * size_factor, units = 'in', res = 600, compression = "lzw")
    corrplot(pearson, is.corr=FALSE, col=cols1, tl.col="black", diag=TRUE, addCoef.col="black", method="color", number.cex=0.5, number.digits = 2 * digit_factor)
    dev.off()
  } else {
    print("Wrong Path Input! Please give either -initial- or -validation-")
  }
}

ListInput <- OldDataList
path <- "initial"
#i = 3
for (i in 1:length(ListInput)) {
  CorrFunction(ListInput[[i]], names(ListInput)[i], path)
}
 