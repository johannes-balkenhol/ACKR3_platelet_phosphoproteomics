##############################################################
#########Signalome Analysis
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
})


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

###################################
### load all differential sites

input = list(top.filter.10, top.filter.30, top.filter.60, top.filter.300, top.filter.600,
             top.filter.900, top.filter.1800)



names_input = c("10", "30", "60", "300", "600", "900", "1800")


#input_tables_names <- list("top.10", "top.600", "top.1800")

#top.rownames <- c(rownames(top.filter.10),rownames(top.600[1:20,]),rownames(top.1800[1:20,]))
#top.norm_intensity <- norm_intensity[top.rownames,]
#rownames(top.norm_intensity) <- sapply(strsplit(rownames(top.norm_intensity), ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})

###################################
### select top differential sites
#i <- 0

# Loop with increment of 3
#for (iter in 1:10) {  # Replace 10 with the desired number of iterations
#  i <- i + 3  # Increment i by 3
#  print(i)  # Print the updated value of i
#}




cr <- 1
## for loop okay, but to many errors still
#for (cr in 1:length(input)) {



#input2 <- input[c(1,5,7)]
input2 <- input[c(1:7)]

####analysis for all time points together####

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
tiff(filename = paste0("analysis/signalome/KSS", "all_top3.tiff"),
     width = 8 * 300, 
     height = 8* 300,
     res = 300,
     compression = "lzw")

print(kinaseSubstrateHeatmap(kssMat,top = 3))

dev.off()
 


####loop for separate analyses for each time point####


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

  ##############################################################
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
tiff(filename = paste0("analysis/signalome/KSS","300.tiff"),
     width = 6 * 300, 
     height = 6* 300,
     res = 300,
     compression = "lzw")

print(kinaseSubstrateHeatmap(kssMat300,top = 3))

dev.off()

#kssMat$combinedScoreMatrix[, which(colnames(kssMat$combinedScoreMatrix)== "PRKACA")] #gives scores for PKA targets

  ##############################################################
  ### export tables 

  ### table of the scoring matrix
  kss <- kssMat[["combinedScoreMatrix"]]
  kks_act <- kssMat[["ksActivityMatrix"]]
  # kss1[kks1>0.8,]
  kss1 <- kss
  rownames(kss1) <- rownames(mat.reg)
  assign(paste0("kss1.", names_input[[i]]), kss1)
  #write.table(kss1, paste0("analysis/signalome/kss_all_combinedscore",".txt"), sep="\t", , row.names=TRUE)
  write.table(kss1, paste0("analysis/signalome/kss", names_input[[i]] ,".txt"), sep="\t", , row.names=TRUE)

  #filter <- apply(kss1, 1, function(x) length(x[x>0.7])>=1)
  #kss1_filt <- kss1[filter,]
  #kss1_filt2 <- as.data.frame(kss1_filt[kss1_filt[,"ROCK1"]>0.7,"ROCK1"])
  #colnames(kss1_filt2) <- "ROCK1"

  ## get the target of the kinases
  #kss_logfc_time <- logfc3[rownames(logfc3) %in% rownames(kss1_filt2),]

  #write.table(kss_logfc_time, "ROCK1_target_intermediate_logfc.txt", sep="\t", , row.names=TRUE)

  ## get the kinases
  #logfc_kss <- logfc3 %>% slice(grep("PRK", row.names(.)))

#}

  ##############################################################
  ### create signalome

  ### PhosR uses the ‘kinaseSubstratePred’ function to synthesise the scores generated from ‘kinaseSubstrateScore’ 
  #to predict the kinase-substrate relationships using an adaptive sampling-based positive-unlabeled learning method (Yang et al., 2018)
  set.seed(1)
  predMat <- kinaseSubstratePred(kssMat, top = 30, cs = 0.8, inclusion = 20, iter = 5, verbose = TRUE)
  #predMat <- kinaseSubstratePred(kssMat, top=30, verbose = FALSE) 
  assign(paste0("predMat.", names_input[[i]]), predMat)
  pheatmap(mat = predMat)

  # tiff(filename = paste0("analysis/signalome/heatmap_predictions", names_input[[7]] ,".tiff"),
  #      width = 12 * 300, 
  #      height = 12 * 300,
  #      res = 300,
  #      compression = "lzw")
  # pheatmap(mat = predMat)
  # dev.off()
}


  ##############################################################
  ### plot the abundance of interesting sites


## Use grep to extract rows that match the patterns
patterns <- c("CAVIN2;S192", "CLASP1;S559","BIN;280|S285","GAS2L1;S297|S306")
matched_rows <- Tc[grep(paste(patterns, collapse = "|"), rownames(Tc)), ]
matched_rows <- top.all[grep(paste(patterns, collapse = "|"), rownames(top.all)), ]

matched_rows <- Tc[rownames(kss1),]
rownames(matched_rows) <. rownames(kss)


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



### Constructing Signalling Networks (Signalomes)
kinaseOI = c("PRKACA")
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
                  threshold = 0.95, 
                  color = kinase_all_color)
