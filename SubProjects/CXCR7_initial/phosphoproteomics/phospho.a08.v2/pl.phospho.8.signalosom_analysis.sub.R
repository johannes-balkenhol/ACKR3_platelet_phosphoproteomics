#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#########Signalome Analysis####
### use mysql database
### get norm abundance
### use parserNorm.pl
### get norm abundance with single phosphosite
suppressPackageStartupMessages({
  library(PhosR)
  library(dplyr)
  library(ggplot2)
  library(GGally)
  library(ggpubr)
  library(calibrate)
  library(network)
  library(tidyr)
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
  library(reshape2)
  library(matrixStats)
  library("stringr")      
})


setwd("/Users/ozgeosmanoglu/Nextcloud/Phosphoproteom/phosphoproteom analysis 2/")


### load data
data("KinaseMotifs")
data("KinaseFamily")
data("phospho_L6_ratio_pe")
data("SPSs")
data("PhosphoSitePlus")


### table of all kinase subunits
pka_targets_psplus <- PhosphoSite.human[["PRKACA"]]
pkg_targets_psplus <- PhosphoSite.human[["PRKG1"]] #prkcg is pkc subunit, this should be changed

write.table(pka_targets_psplus, "pka_targets_psiteplus.txt", sep="\t", , row.names=TRUE)
write.table(pkg_targets_psplus, "pkg_targets_psiteplus.txt", sep="\t", , row.names=TRUE)


### load all differential sites

input = list(top.filter.10, top.filter.30, top.filter.60, top.filter.300, top.filter.600,
             top.filter.900, top.filter.1800)



names_input = c("10", "30", "60", "300", "600", "900", "1800")

input2 <- input[c(1:7)]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########analysis for all time points together###########
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (i in 1:length(input2)) {
  #top.names <- append(top.names, rownames(input2[[i]][input2[[i]]$PValue<0.05,]))
  #assign(paste0(".", names_input[[i]]), )
  if (i == 1 & i <= 2) { 
    top.names <- rownames(input2[[i]][input2[[i]]$PValue<0.05 & abs(input2[[i]]$logFC)>0.5,])
  } else {
    print(as.list(test[[i]]))
    top.names <- append(top.names, rownames(input2[[i]][input2[[i]]$PValue<0.05 & abs(input2[[i]]$logFC)>0.5,]))
  }
}

top.names <- unique(top.names) #duplicates should be removed!
norm_intensity2 = norm_intensity[!(row.names(norm_intensity) %in% not.uniq[,1]), ]
norm_intensity2 <- as.data.frame(norm_intensity2)

#norm_intensity2 <- meanAbundance(norm_intensity2, grps)
top.norm_intensity <- norm_intensity2[top.names,]  ##hmmm here it creates duplicates, why?
top.norm_intensity <- na.omit(top.norm_intensity)


### determine kinase substrate interaction

mat.reg <- as.matrix(top.norm_intensity)

### standardize matrix
mat.std <- PhosR::standardise(mat.reg)

rownames(mat.std) <- sapply(strsplit(rownames(mat.std), ";"), function(x) { 
  gsub(" ", "", paste(toupper(x[2]), x[3], "", sep=";"))
})

### Run PhosR kinase-substrate prediction with the default parameters
data("KinaseMotifs")

seqs <- sapply(strsplit(rownames(mat.reg), ";"), function(x) { gsub(" ", "", x[4]) })

kssMat <- kinaseSubstrateScore(substrate.list = PhosphoSite.human, 
                               mat = mat.std, seqs = seqs,
                               numMotif = 5, numSub = 1, verbose = FALSE)

write.table(kssMat$combinedScoreMatrix, "analysis/signalome/kssMat_global.txt", sep = "\t", row.names = TRUE)
write.table(kssMat$ksActivityMatrix, "analysis/signalome/kinaseactivities_global.txt", sep = "\t", row.names = TRUE)



####### print the heatmap of the top targets
tiff(filename = paste0("analysis/signalome/KSS", "all_top3.tiff"),
     width = 8 * 300, 
     height = 8* 300,
     res = 300,
     compression = "lzw")

print(kinaseSubstrateHeatmap(kssMat,top = 3))

dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####function for separate analyses for each time point####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

separate_kss <- function(i) {
  #top.names <- append(top.names, rownames(input2[[i]][input2[[i]]$PValue<0.05,]))
  #assign(paste0(".", names_input[[i]]), )
  top.names <- rownames(input2[[i]][input2[[i]]$PValue<0.05 & abs(input2[[i]]$logFC)>0.5,])
  
  #if we don't want a loop but use all timepoints together, use the for loop above with if&else
  #(so close this loop before) and continue from here downwards.
  
  
  norm_intensity2 = norm_intensity[!(row.names(norm_intensity) %in% not.uniq[,1]), ]
  norm_intensity2 <- as.data.frame(norm_intensity2)
  
  #norm_intensity2 <- meanAbundance(norm_intensity2, grps)
  top.norm_intensity <- norm_intensity2[top.names,]
  top.norm_intensity <- na.omit(top.norm_intensity)
  
  ### determine kinase substrate interaction
  
  mat.reg <- as.matrix(top.norm_intensity)
  
  ### standardize matrix
  mat.std <- PhosR::standardise(mat.reg)
  
  rownames(mat.std) <- sapply(strsplit(rownames(mat.std), ";"), function(x) { 
    gsub(" ", "", paste(toupper(x[2]), x[3], "", sep=";"))
  })
  
  ### Run PhosR kinase-substrate prediction with the default parameters
  data("KinaseMotifs")
  #seqs2 <- sapply(strsplit(rownames(mat.reg), ";"), "[[", 4)
  seqs <- sapply(strsplit(rownames(mat.reg), ";"), function(x) { gsub(" ", "", x[4]) })
  
  kssMat <- kinaseSubstrateScore(substrate.list = PhosphoSite.human, 
                                 mat = mat.std, seqs = seqs,
                                 numMotif = 5, numSub = 1, verbose = FALSE)
  
  assign(paste0("kssMat", names_input[[i]]), kssMat, globalenv())
  
  
  
  
  ####### print the heatmap of the top targets
  tiff(filename = paste0("analysis/signalome/KSS",names_input[i], ".tiff"),
       width = 8 * 300, 
       height = 8* 300,
       res = 300,
       compression = "lzw")
  
  print(kinaseSubstrateHeatmap(kssMat,top = 3))
  
  dev.off()
  
}

for (j in 1:length(input2)) {
  separate_kss(j)
}

#manual saving to tiff to change the dimensions
tiff(filename = paste0("analysis/signalome/KSS","60.tiff"),
     width = 6 * 300, 
     height = 4* 300,
     res = 300,
     compression = "lzw")

print(kinaseSubstrateHeatmap(kssMat60,top = 3))

dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####extract predicted targets of specific kinases and make heatmap of significants####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data(PhosphoSitePlus)

dataset <- SummarizedExperiment::assay(ppe, "normalised")
dataset_df <- as.data.frame(dataset)
dataset_df <- dataset_df[!(row.names(dataset_df) %in% not.uniq[,1]), ]

#get the list of significant ids in all
union_sig <- Reduce(union, list(rownames(top.10.sign), rownames(top.30.sign), 
                                rownames(top.60.sign), row.names(top.300.sign),
                                rownames(top.600.sign), rownames(top.900.sign), 
                                rownames(top.1800.sign)))
dataset_df <- dataset_df[rownames(dataset_df) %in% union_sig, ]
name <- sapply(strsplit(rownames(dataset_df), ";"), "[[", 2)
site <- sapply(strsplit(rownames(dataset_df), ";"), "[[", 3)
dataset_df$namesite <- paste0(name,";",site,";")

##specific kinases####
kinase_list <- list("PRKACA", "CDK2", "MAPK3", "PRKG1", "PRKG2", "MTOR", "SRC", "CDK1", "CDK5")
kinase <- kinase_list[[9]]
kinase_scores <- data.frame(kssMat$combinedScoreMatrix[,kinase])
kinase_scores$name <- rownames(kssMat$combinedScoreMatrix)
colnames(kinase_scores) <- c("score", "id")
kinase_scores <- kinase_scores[ kinase_scores$score > 0.75,] 
predicted_kinase_targets <- unique(kinase_scores$id)
kinasename  <- paste0(kinase, "_predicted")
kinase_targets <- predicted_kinase_targets
targets_norm_abundance <- dataset_df[dataset_df$namesite %in% kinase_targets,1:53]
###




targets_norm_abundance <-targets_norm_abundance[order(row.names(targets_norm_abundance)),]
targets_norm_abundance <- as.matrix(targets_norm_abundance)

targets_pvals <- top.all[rownames(top.all) %in% rownames(targets_norm_abundance), c(5,8,11,14,17,20,23)]
targets_pvals<-targets_pvals[order(row.names(targets_pvals)), ]

targets_pvals <- as.matrix(targets_pvals)
rownames(targets_pvals) <- sapply(strsplit(rownames(targets_pvals), ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})
targets_pvals <- as.data.frame(targets_pvals)
targets_pvals[targets_pvals <= 0.05] <- "*"
targets_pvals[targets_pvals> 0.05] <- ""
targets_pvals$adj.P.Val.00 <- ""
targets_pvals <- targets_pvals[c(8,1:7)]



targets_logfcs <- top.all[rownames(top.all) %in% rownames(targets_norm_abundance), c(3,6,9,12,15,18,21)]
targets_logfcs<-targets_logfcs[order(row.names(targets_logfcs)), ]
targets_logfcs <- as.matrix(targets_logfcs)
rownames(targets_logfcs) <- sapply(strsplit(rownames(targets_logfcs), ";"),  
                                   function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})
targets_logfcs <- as.data.frame(targets_logfcs)

targets_logfcs<-round(targets_logfcs, digits = 3)

for (j in 1:(length(colnames(targets_pvals))-1)) {
  print(j)
  for (i in 1:length(rownames(targets_logfcs))) {
    if (targets_pvals[i, j+1] == "") {
      targets_logfcs[i, j] <- ""}}}

targets_logfcs$logFC.00<- ""
targets_logfcs <- targets_logfcs[c(8,1:7)]
targets_logfcs<-targets_logfcs[order(row.names(targets_logfcs)), ]



rownames(targets_norm_abundance) <- sapply(strsplit(rownames(targets_norm_abundance), ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})
x_tt <- as.factor(grps)


t00 <- rowMedians(targets_norm_abundance[,1:7])
t10 <- rowMedians(targets_norm_abundance[,8:13])
t30 <- rowMedians(targets_norm_abundance[,14:20])
t60 <- rowMedians(targets_norm_abundance[,21:27])
t300 <- rowMedians(targets_norm_abundance[,28:34])
t600 <- rowMedians(targets_norm_abundance[,35:40])
t900 <- rowMedians(targets_norm_abundance[,41:47])
t1800 <- rowMedians(targets_norm_abundance[,48:53])


targets_norm_abundance.collapse <- cbind(t00,t10,t30,t60,t300,t600,t900,t1800)
rownames(targets_norm_abundance.collapse) <- rownames(targets_norm_abundance)
colnames(targets_norm_abundance.collapse) <- c("t00","t10","t30","t60","t300","t600","t900", "t1800")
targets_norm_abundance.collapse <- as.data.frame(targets_norm_abundance.collapse)

targets_norm_abundance2 <- targets_norm_abundance.collapse


#dd <- melt(as.matrix(targets_norm_abundance2), variable.name = "normalised")
#ggplot(dd, aes(value, colour = Var2)) + geom_density(bw = "SJ")



## use median, as more robust than average, for the Z values in heatmap
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

########### coloring according to sample names
col_groups <- substr(colnames(targets_norm_abundance2), 1, 5)
table(col_groups)


mat_col <- data.frame(time = col_groups)
rownames(mat_col) <- colnames(targets_norm_abundance2)

mat_colors <- list(time = brewer.pal(8, "BuPu"))
#mat_colors <- list(group = brewer.pal(length(table(col_groups)), "Set3"))

names(mat_colors$time) <- unique(col_groups)



############ create a simple heatmap
targets_data_subset_norm <- t(apply(targets_norm_abundance2, 1, cal_z_score))  ##applies to rows
#rownames(data_subset_norm) <- c(sapply(strsplit(rownames(data_subset_norm), ";"), "[[", c(1,2,3)), ";", sapply(strsplit(rownames(data_subset_norm), ";"), "[[", 2), ";" , sapply(strsplit(rownames(data_subset_norm), ";"), "[[", 3))
#rownames(data_subset_norm) <- sapply(strsplit(rownames(data_subset_norm), ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})


############ include breaks in the heatmap 
############ for better visualization in tailed data
############ as we use a color code cutoff regarding quantiles

mat_breaks <- seq(min(targets_norm_abundance2, na.rm=TRUE), max(targets_norm_abundance2, na.rm=TRUE), length.out = 20)


# define function
quantile_breaks <- function(xs, n) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm=TRUE)
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(targets_data_subset_norm, n = 40)


test_labels <- as.matrix(targets_logfcs) 



tiff(filename = paste0("analysis/signalome/Heatmap_", kinasename, "_sign.tiff", sep =""),
     width = 10 * 300, 
     height = 13 * 300,
     res = 300,
     compression = "lzw")

pheatmap(mat = targets_data_subset_norm, 
         color = colorRampPalette(c("navy", "white", "red"))(length(mat_breaks) - 1),
         gaps_row=c(5,10,15,20,25,30,35,40,45,50, 55),
         scale="row",
         na_col = "grey",
         breaks = mat_breaks,
         border_color = "white", 
         show_colnames = TRUE, 
         show_rownames = TRUE, 
         annotation_col = mat_col, 
         annotation_colors = mat_colors, 
         drop_levels = TRUE, 
         fontsize = 10, 
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         cex=1,
         clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean",
         clustering_method="complete",
         main = paste(kinase, "predicted targets"),
         display_numbers = test_labels,
         number_color = "green", 
         fontsize_number = 8,
         cellheight=15, cellwidth = 30
)

dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############ create signalome ########################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### PhosR uses the ‘kinaseSubstratePred’ function to synthesise the scores generated from ‘kinaseSubstrateScore’ 
#to predict the kinase-substrate relationships using an adaptive sampling-based positive-unlabeled learning method (Yang et al., 2018)

set.seed(1)
predMat <- kinaseSubstratePred(kssMat, top = 30, cs = 0.8, inclusion = 20, iter = 5, verbose = TRUE)
#predMat <- kinaseSubstratePred(kssMat, top=30, verbose = FALSE) 
#assign(paste0("predMat.", names_input[[i]]), predMat)
pheatmap(mat = predMat)

# tiff(filename = paste0("analysis/signalome/heatmap_predictions", names_input[[7]] ,".tiff"),
#       width = 12 * 300, 
#       height = 12 * 300,
#       res = 300,
#       compression = "lzw")
#  pheatmap(mat = predMat)
#  dev.off()


### Constructing Signalling Networks (Signalomes)
kinaseOI = c("CDK2","MAPK3")
signalomesRes <- Signalomes(KSR = kssMat, 
                            predMatrix = predMat, 
                            exprsMat = mat.std, 
                            module_res = 10,
                            KOI = kinaseOI)

### generate palette
my_color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))
kinase_all_color <- my_color_palette(ncol(kssMat$combinedScoreMatrix))
names(kinase_all_color) <- colnames(kssMat$combinedScoreMatrix)
kinase_signalome_color <- kinase_all_color[colnames(predMat)]

plotSignalomeMap(signalomes = signalomesRes,
                 color = kinase_signalome_color)


### plot the signalome network that illustrates the connectivity between kinase signalome networks.


plotKinaseNetwork(KSR = kssMat, 
                  predMatrix = predMat, 
                  threshold = 0.9, 
                  color = kinase_all_color)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### plot the abundance of interesting sites####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## Use grep to extract rows that match the patterns
patterns <- c("CAVIN2;S192", "CLASP1;S559","BIN;280|S285","GAS2L1;S297|S306")
matched_rows <- Tc[grep(paste(patterns, collapse = "|"), rownames(Tc)), ]
matched_rows <- top.all[grep(paste(patterns, collapse = "|"), rownames(top.all)), ]


kss <- kssMat[["combinedScoreMatrix"]]


matched_rows <- Tc[rownames(kss1),]
rownames(matched_rows) <- rownames(kss)


df <- matched_rows[,1:3]
df$rowname <- rownames(df)

# Get the number of columns in the dataframe
num_rows <- nrow(df)

# Set color palette for lines based on number of columns
colors <- viridis::viridis(num_rows)

# Convert dataframe to long format
df_long <- tidyr::pivot_longer(df, cols = -rowname, names_to = "variable", values_to = "value")

# Filter to keep only the last column
df_long_last_col <- df_long[df_long$variable == colnames(df)[ncol(df)-1], ]


# Convert variable to factor with custom levels to specify x-axis order
df_long$variable <- factor(df_long$variable, levels = colnames(df))

# Create the plot
plot <- ggplot(df_long, aes(x = variable, y = value, color = rowname, group = rowname)) +
  geom_line(size = 1.5) +
  geom_label(data = df_long_last_col, aes(label = rowname), hjust = -0.2, vjust = 0.5, size = 3, show.legend = FALSE) +
  scale_color_manual(values = colors) +
  labs(title = "Line Plot of DataFrame Rows",
       x = "X Axis Label",
       y = "Y Axis Label")

# Print the plot
print(plot)

