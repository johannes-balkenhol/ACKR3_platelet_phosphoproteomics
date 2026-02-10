
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(reshape2)
library(ggplot2)
library(PhosR)
library(dplyr)
library(matrixStats)
library(stringi)
library(stringr)

data("KinaseMotifs")
data("KinaseFamily")
data("phospho_L6_ratio_pe")
data("SPSs")
data("PhosphoSitePlus")

setwd("D:/Research/CXCR7_platelet_analysis/SubProjects/CXCR7_validation/phosphoproteomics/scripts/R")

# write.table(PhosphoSite.human[["SRC"]], "../../data/processed_data/SRC_targets_psplus.txt", row.names = F, col.names = F)
# write.table(PhosphoSite.human[["CSNK2A1"]], "../../data/processed_data/CSNK2A1_targets_psplus.txt", row.names = F, col.names = F)
# write.table(PhosphoSite.human[["LYN"]], "../../data/processed_data/LYN_targets_psplus.txt", row.names = F, col.names = F)
# write.table(PhosphoSite.human[["CDK1"]], "../../data/processed_data/CDK1_targets_psplus.txt", row.names = F, col.names = F)
# write.table(PhosphoSite.human[["CDK2"]], "../../data/processed_data/CDK2_targets_psplus.txt", row.names = F, col.names = F)
# write.table(PhosphoSite.human[["MAPK13"]], "../../data/processed_data/MAPK13_targets_psplus.txt", row.names = F, col.names = F)
# write.table(PhosphoSite.human[["PRKACA"]], "../../data/processed_data/PRKACA_targets_psplus.txt", row.names = F, col.names = F)
# write.table(PhosphoSite.human[["PRKG1"]], "../../data/processed_data/PRKG1_targets_psplus.txt", row.names = F, col.names = F)
# write.table(PhosphoSite.human[["LCK"]], "../../data/processed_data/kinaseTargets/LCK_targets_psplus.txt", row.names = F, col.names = F)


#pkg_targets_psplus <- PhosphoSite.human[["PRKG1"]]
#pka_targets_psplus <- PhosphoSite.human[["PRKACA"]]

dataset <- norm_intensity_filter
dataset_df <- as.data.frame(dataset)

# ### define input 
# input = list(top.filter.10, top.filter.600, top.filter.1800, 
#              top.filter.10.dmso.vs.0s, top.filter.600.dmso.vs.0s, top.filter.1800.dmso.vs.0s,
#              top.filter.10.cxcr7.vs.0s, top.filter.600.cxcr7.vs.0s, top.filter.1800.cxcr7.vs.0s)
# 
# names_input = c("10", "600", "1800", 
#                 "10.dmso.vs.0s", "600.dmso.vs.0s", "1800.dmso.vs.0s",
#                 "10.cxcr7.vs.0s", "600.cxcr7.vs.0s", "1800.cxcr7.vs.0s")
# 
# input <-  input[c(1:3)]
# for (i in 1:length(input)) {
#   #df_pl = top.10.cxcr7.vs.0s
#   df_pl2 = input[[i]]
#   df_pl2 <- df_pl2[df_pl2[, "PValue"] <0.05,]
#   assign(paste0("top.filter.", names_input[[i]], ".sign"), df_pl2)
# }
# union_sig <- Reduce(union, list(rownames(top.filter.10.sign),rownames(top.filter.600.sign),rownames(top.filter.1800.sign)))
# signames <- unique(sapply(strsplit(union_sig, ";"), "[[", 1))


#all #use norm_intensity_filter as dataset_df directly.

#continue
name <- sapply(strsplit(rownames(dataset_df), ";"), "[[", 2)
site <- sapply(strsplit(rownames(dataset_df), ";"), "[[", 3)

dataset_df$namesite <- paste0(name,";",site,";")

##specific kinases
kinases = c("PRKACA", "PRKG1", "PRKAA1", "CDK1", "MAPKAPK2", "AKT1", "MTOR", "SRC",  "CDK2", "CSNK2A1", "LCK", "GSK3B", "PRKD1", "MAPK3", "CDK5", "PRKCA",  "MAPK14")

# PLATELET CC ----

CC_proteins <- c("RGS18","P2RY1","PIK3R1","F2R","ADCY5","ADCY3","ADCY1","HTR2A","ADCY8","AKT2",
                 "ARHGEF1","COL1A2","COL3A1","GP1BA","TBXA2R","ITPR1","GP1BB","P2RX1","ORAI1","ITGA2B",
                 "PIK3R5","PIK3CB","PIK3R6","PLA2G4A","RAC1","RAP1A","RAP1B","RHOA","AKT1","STIM1",
                 "F2","TRPM7","RASGRP2","TRPC6","GP5","ADCY7","COL1A1","F2RL3","MAPK1","ROCK1",
                 "PIK3R2","PIK3CD","PLCB4","ADCY6","PLCB3","PIK3CA","P2RY12","PRKG1","PIK3CG","GUCY1A1",
                 "ADCY9","PLCG2","PIK3R3","PRKACA","PLCB1","ADCY4","SRC","PLCB2","GUCY1A2","PTGIR",
                 "GP9","PRKCA","ITGB3","GP6","ADCY2","GUCY1B1","VWF","TLN1","VASP")

CC_proteins_sub <- c("GP1BB", "ITGB3", "ADCY3", "GP5", "TLN1", "STIM1", "VWF", "TBXA2R", "ITGA2B",
                     "PIK3R3", "PRKCA", "RAC1")

CC_proteins_sub2 <- c("F2R", "RGS18", "P2RY12", "RAP1A", "ROCK1", "P2RY1", "HTR2A", "PIK3R1", "ITPR1", "PIK3CG")

proteinlist <- paste(CC_proteins_sub2, collapse = ";|")


#continue ----
#kinases = c(paste0("GRK", "1":"7"))
top.all.input.list <- list("DMSOvs0" = top.all.dmso.vs.0s, 
                           "ACKR3vs0" = top.all.cxcr7.vs.0s, 
                           "ACKR3vsDMSO" = top.all)


x = "SRC"
ExpCondition <- "ACKR3vsDMSO"


kinase_function <- function(x, ExpCondition) {
  top.all.input <- top.all.input.list[[ExpCondition]]
  kinase <- PhosphoSite.human[[x]]
  #targets_norm_abundance <- dataset_df[rownames(dataset_df) %in% proteins, 1:70]
  #targets_norm_abundance <- dataset_df[str_detect(row.names(dataset_df), proteinlist), 1:70]
  
  targets_norm_abundance <- dataset_df[dataset_df$namesite %in% kinase,1:70]
  targets_norm_abundance <-targets_norm_abundance[order(row.names(targets_norm_abundance)), ]
  targets_norm_abundance <- as.matrix(targets_norm_abundance)
  targets_pvals <- top.all.input[rownames(top.all.input) %in% rownames(targets_norm_abundance), c(5,8,11)]
  targets_pvals<-targets_pvals[order(row.names(targets_pvals)), ]
  targets_pvals <- as.matrix(targets_pvals)
  rownames(targets_pvals) <- sapply(strsplit(rownames(targets_pvals), ";"),  function(x){paste(x[[2]], " (",  x[[3]], ")", sep = "")})
  targets_pvals <- as.data.frame(targets_pvals)
  targets_pvals[targets_pvals <= 0.05] <- "*"
  targets_pvals[targets_pvals> 0.05] <- ""
  targets_pvals$adj.P.Val.00 <- ""
  targets_pvals <- targets_pvals[c(4,1:3)]
  
  
  targets_logfcs <- top.all.input[rownames(top.all.input) %in% rownames(targets_norm_abundance), c(3,6,9)]
  targets_logfcs<-targets_logfcs[order(row.names(targets_logfcs)), ]
  targets_logfcs <- as.matrix(targets_logfcs)
  rownames(targets_logfcs) <- sapply(strsplit(rownames(targets_logfcs), ";"),  
                                     function(x){paste(x[[2]], " (",  x[[3]], ")", sep = "")})
  targets_logfcs <- as.data.frame(targets_logfcs)
  
  targets_logfcs<-round(targets_logfcs, digits = 2)
  
  
  #to plot log2FCs
  logfc_input <- targets_logfcs
  logfc_input$logFC.00<- 0
  logfc_input <- logfc_input[c(4,1,2,3)]
  colnames(logfc_input) <- c("t00","t10_CvsD","t600_CvsD","t1800_CvsD")
  
  
  for (j in 1:(length(colnames(targets_pvals))-1)) {
    print(j)
    for (i in 1:length(rownames(targets_logfcs))) {
      if (targets_pvals[i, j+1] == "") {
        targets_logfcs[i, j] <- ""}}}
  
  #targets_pvals$adj.P.Val.10.DMSO <- ""
  #targets_pvals$adj.P.Val.600.DMSO <- ""
  #targets_pvals$adj.P.Val.1800.DMSO <- ""
  #targets_pvals <-  targets_pvals[c(1,5,2,6,3,7,4)]
  #targets_pvals <-  targets_pvals[c(1,2,5,3,6,4,7)]#for dmso
  
  
  targets_logfcs$logFC.00<- ""
  targets_logfcs <- targets_logfcs[c(4,1,2,3)]
  #targets_logfcs$logFC.10.DMSO<- ""
  #targets_logfcs$logFC.600.DMSO<- ""
  #targets_logfcs$logFC.1800.DMSO<- ""
  #targets_logfcs <- targets_logfcs[c(4,5,1,6,2,7,3)]
  #targets_logfcs <- targets_logfcs[c(4,1,5,2,6,3,7)] #for dmso
  
  targets_logfcs<-targets_logfcs[order(row.names(targets_logfcs)), ]
  
  rownames(targets_norm_abundance) <- sapply(strsplit(rownames(targets_norm_abundance), ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})
  
  x_tt <- as.factor(grps)
  t00 <- rowMedians(targets_norm_abundance[,1:10])
  t10_C <- rowMedians(targets_norm_abundance[,11:20])
  t10_D <- rowMedians(targets_norm_abundance[,21:30])
  t600_C <- rowMedians(targets_norm_abundance[,51:60])
  t600_D <- rowMedians(targets_norm_abundance[,61:70])
  t1800_C <- rowMedians(targets_norm_abundance[,31:40])
  t1800_D <- rowMedians(targets_norm_abundance[,41:50])
  
  
  targets_norm_abundance.collapse <- cbind(t00,t10_D,t10_C,t600_D,t600_C,t1800_D,t1800_C)
  rownames(targets_norm_abundance.collapse) <- rownames(targets_norm_abundance)
  colnames(targets_norm_abundance.collapse) <- c("t00","t10_DMSO","t10_ACKR3","t600_DMSO","t600_ACKR3","t1800_DMSO","t1800_ACKR3")
  targets_norm_abundance.collapse <- as.data.frame(targets_norm_abundance.collapse)
  
  if (ExpCondition =="ACKR3vsDMSO") {
    targets_norm_abundance2 = logfc_input
  } else {
    targets_norm_abundance2 <- targets_norm_abundance.collapse[ , 
                                                                grepl(gsub("vs0", "", ExpCondition), names(targets_norm_abundance.collapse))]
    
    targets_norm_abundance2 <- cbind(t00, targets_norm_abundance2)
  }
  
  colnames(targets_norm_abundance2) <- sapply(strsplit(colnames(targets_norm_abundance2), "_"), "[[", 1)
  
  
  
  #dd <- melt(as.matrix(targets_norm_abundance2), variable.name = "normalised")
  #ggplot(dd, aes(value, colour = Var2)) + geom_density(bw = "SJ")
  ###############################################################################
  ## use median, as more robust than average, for the Z values in heatmap
  cal_z_score <- function(x){
    (x - mean(x)) / sd(x)
  }
  
  ########### coloring according to sample names
  col_groups <- substr(colnames(targets_norm_abundance2), 1, 11)

  #table(col_groups)
  
  
  mat_col <- data.frame(time = col_groups)
  rownames(mat_col) <- colnames(targets_norm_abundance2)
  
  mat_colors <- list(time = brewer.pal(7, "BuPu"))
  mat_colors$time <- mat_colors$time[c(1,2,4,6)] #for DMSO
  #mat_colors <- list(group = brewer.pal(length(table(col_groups)), "Set3"))
  
  names(mat_colors$time) <- unique(col_groups)
  
  
  
  ############ create a simple heatmap
  #targets_data_subset_norm <- t(apply(targets_norm_abundance2, 1, cal_z_score))
  targets_data_subset_norm <- targets_norm_abundance2
  #rownames(data_subset_norm) <- c(sapply(strsplit(rownames(data_subset_norm), ";"), "[[", c(1,2,3)), ";", sapply(strsplit(rownames(data_subset_norm), ";"), "[[", 2), ";" , sapply(strsplit(rownames(data_subset_norm), ";"), "[[", 3))
  #rownames(data_subset_norm) <- sapply(strsplit(rownames(data_subset_norm), ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})
  
  
  ############ include breaks in the heatmap 
  ############ for better visualization in tailed data
  ############ as we use a color code cutoff regarding quantiles
  
  #mat_breaks <- seq(min(targets_norm_abundance2, na.rm=TRUE), max(targets_norm_abundance2, na.rm=TRUE), length.out = 20)
  mat_breaks <- seq(-2, 2, length.out = 40)
  # define function
  
  quantile_breaks <- function(xs, n) {
    breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm=TRUE)
    breaks[!duplicated(breaks)]
  }
  
  #mat_breaks <- quantile_breaks(targets_data_subset_norm, n = 40)
  
  # Replace non-empty cells with "*"
  
  targets_logfcs[targets_logfcs != ""] <- stri_unescape_unicode("\u2217")
  
  # assign to the labels
  test_labels <- as.matrix(targets_logfcs) 
  
  
  # Define the threshold for dynamic text colors
  threshold <- mat_breaks[10]
  # Generate dynamic text colors based on cell values
  number_colors <- ifelse(targets_data_subset_norm < threshold, "white", "black")
  
  ph <- pheatmap(
    mat = targets_data_subset_norm,
    color = colorRampPalette(c("navy", "white", "red"))(length(mat_breaks) - 1),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    breaks = mat_breaks,
    silent = TRUE # Suppress immediate plotting to adjust data
  )
  
  # Reorder `number_colors` to match clustering
  reordered_number_colors <- number_colors[ph$tree_row$order, ]
  
  tiff(filename = paste0("../analysis/Heatmap/", ExpCondition, "/Shortened/Heatmap_", x, "_targets.tiff", sep =""), #for kinases
        #tiff(filename = paste0("analysis/Heatmap/Heatmap_", protein, "_sign.tiff", sep =""), #for proteins
        width = 10 * 300, 
        height = 10 * 300,
        res = 300,
        compression = "lzw")
 
 #for the asterisk to be visible, we need to save it with cairo_pdf.
  cairo_pdf(filename = paste0("../analysis/Heatmap/", ExpCondition, "/Shortened/nonstandardized-log2FCs/Heatmap_", x, "_targetsRaw2.pdf", sep =""), #for kinases
 #  #        #tiff(filename = paste0("analysis/Heatmap/Heatmap_", protein, "_sign.tiff", sep =""), #for proteins)
  )
 #  
#   cairo_pdf(filename = paste0("../analysis/Clusters/gene_centric/Diff_all/Standardized/", 
 #                              gsub(" |/", "_", description), "_Heatmapsub.pdf", sep =""),
  #)
  pheatmap(mat = targets_data_subset_norm, 
           color = colorRampPalette(c("navy", "white", "red"))(length(mat_breaks) - 1),
           #gaps_row=c(5,10,15,20,25,30,35,40,45,50, 55),
           scale="none",
           na_col = "grey",
           breaks = mat_breaks,
           border_color = "gray30", 
           show_colnames = TRUE, 
           show_rownames = TRUE, 
           #annotation_col = mat_col, 
           #annotation_colors = mat_colors, 
           drop_levels = TRUE, 
           fontsize = 10, 
           cluster_cols = FALSE,
           cluster_rows = T,
           cex=1,
           clustering_distance_rows="euclidean",
           clustering_distance_cols="euclidean",
           clustering_method="complete",
           main = paste(x, "targets"),
           display_numbers = test_labels[rownames(targets_data_subset_norm),],
           number_color = reordered_number_colors,
           fontsize_number = 20,
           cellheight=10, cellwidth = 25
  )
  
  dev.off()
  
 }


###
i=1
j = 1
for (j in 1:length(names(top.all.input.list))) {
  condition = names(top.all.input.list)[j]
  for (i in 1:length(kinases)) {
    kinasename <-  kinases[[i]]
    kinase_function(kinasename, condition)
  }
}


names(top.all.input.list)

##specific proteins####

#PDE3A;Q14432
#VASP; P50552
#ENSA; O43768
#AMPK subunits/isoforms
##AMPK; Q9Y478
##AMPK1; Q13131
##PRKAB2: O43741
##AMPKG3; Q9UGI9
##PRKAG1; P54619
##PRKAG2; Q9UGJ0
##AMPK2; P54646
#ACACA; Q13085
#ACACB; O00763

protein_list <- c("Q14432", "P50552", "O43768", "Q9Y478", "Q13131", "O43741", 
                  "Q9UGI9", "P54619", "Q9UGJ0", "P54646", "Q13085", "O00763")
protein <- paste(protein_list, collapse = ";|")

protein <- gsub(";", ";|", CalciumPath)

#protein <-  protein_list[[1]]
#protein <- "GAPDH"
#dataset_df <- top.all[,c(3,6,9)]
#name <- sapply(strsplit(rownames(dataset_df), ";"), "[[", 2)
#site <- sapply(strsplit(rownames(dataset_df), ";"), "[[", 3)
#top.all.input <- top.all.dmso.vs.0s
top.all.input <- top.all
#dataset_df$namesite <- paste0(name,";",site,";")

protein_function <- function(description, proteins) {
  #uncomment for specific prot names:
  #targets_norm_abundance <- dataset_df[str_detect(row.names(dataset_df), protein), 1:70]
  #for cluster phosphosites use following directly:
  # Get normalized intensities for proteins of interest
  targets_norm_abundance <- dataset_df[rownames(dataset_df) %in% proteins, 1:70]
  targets_norm_abundance <-targets_norm_abundance[order(row.names(targets_norm_abundance)), ]
  targets_norm_abundance <- as.matrix(targets_norm_abundance)
  
  # Get p-values for their expression at 3 main timepoints
  targets_pvals <- top.all.input[rownames(top.all.input) %in% rownames(targets_norm_abundance), c(5,8,11)]
  targets_pvals<-targets_pvals[order(row.names(targets_pvals)), ]
  targets_pvals <- as.matrix(targets_pvals)
  rownames(targets_pvals) <- sapply(strsplit(rownames(targets_pvals), ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})
  targets_pvals <- as.data.frame(targets_pvals)
  targets_pvals[targets_pvals <= 0.05] <- "*"
  targets_pvals[targets_pvals> 0.05] <- ""
  targets_pvals$adj.P.Val.00 <- ""
  targets_pvals <- targets_pvals[c(4,1:3)]
  
  # Get logFCs for their expression at 3 main timepoints
  targets_logfcs <- top.all.input[rownames(top.all.input) %in% rownames(targets_norm_abundance), c(3,6,9)]
  targets_logfcs<-targets_logfcs[order(row.names(targets_logfcs)), ]
  targets_logfcs <- as.matrix(targets_logfcs)
  rownames(targets_logfcs) <- sapply(strsplit(rownames(targets_logfcs), ";"),  
                                     function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})
  targets_logfcs <- as.data.frame(targets_logfcs)
  targets_logfcs<-round(targets_logfcs, digits = 2)
  
  # Filter log fold changes by p-values
  for (j in 1:(length(colnames(targets_pvals))-1)) {
    #print(j)
    for (i in 1:length(rownames(targets_logfcs))) {
      if (targets_pvals[i, j+1] == "") {
        targets_logfcs[i, j] <- ""}}}
  
  targets_pvals$adj.P.Val.10.DMSO <- ""
  targets_pvals$adj.P.Val.600.DMSO <- ""
  targets_pvals$adj.P.Val.1800.DMSO <- ""
  targets_pvals <-  targets_pvals[c(1,5,2,6,3,7,4)]
  #targets_pvals <-  targets_pvals[c(1,2,5,3,6,4,7)] #for dmso pvals
  #targets_pvals <-  targets_pvals[c(3,5,7)]
  
  targets_logfcs$logFC.00<- ""
  targets_logfcs$logFC.10.DMSO<- ""
  targets_logfcs$logFC.600.DMSO<- ""
  targets_logfcs$logFC.1800.DMSO<- ""
  targets_logfcs <- targets_logfcs[c(4,5,1,6,2,7,3)]
  #targets_logfcs <- targets_logfcs[c(4,1,5,2,6,3,7)] #for dmso logfcs
  #targets_logfcs <- targets_logfcs[c(3,5,7)]
  targets_logfcs<-targets_logfcs[order(row.names(targets_logfcs)), ]
  
  # prepare rownames for the heatmap
  rownames(targets_norm_abundance) <- sapply(strsplit(rownames(targets_norm_abundance), ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})
  x_tt <- as.factor(grps)
  
  # Calculate median values
  t00 <- rowMedians(targets_norm_abundance[,1:10])
  t10 <- rowMedians(targets_norm_abundance[,11:20])
  t10wt <- rowMedians(targets_norm_abundance[,21:30])
  t600 <- rowMedians(targets_norm_abundance[,51:60])
  t600wt <- rowMedians(targets_norm_abundance[,61:70])
  t1800 <- rowMedians(targets_norm_abundance[,31:40])
  t1800wt <- rowMedians(targets_norm_abundance[,41:50])
  
  targets_norm_abundance.collapse <- cbind(t00,t10wt,t10,t600wt,t600,t1800wt,t1800)
  rownames(targets_norm_abundance.collapse) <- rownames(targets_norm_abundance)
  colnames(targets_norm_abundance.collapse) <- c("t00","t10wt","t10","t600wt","t600","t1800wt","t1800")
  targets_norm_abundance.collapse <- as.data.frame(targets_norm_abundance.collapse)
  
  targets_norm_abundance2 <- targets_norm_abundance.collapse
  #targets_norm_abundance2 <-  targets_norm_abundance
  
  # Prepare annotation for heatmap
  col_groups <- substr(colnames(targets_norm_abundance2), 1, 11)
  mat_col <- data.frame(time = col_groups)
  rownames(mat_col) <- colnames(targets_norm_abundance2)
  mat_colors <- list(time = brewer.pal(7, "BuPu"))
  names(mat_colors$time) <- unique(col_groups)
  test_labels <- as.matrix(targets_logfcs) 
  
  # Calculate z-scores
  cal_z_score <- function(x){
    (x - mean(x)) / sd(x)
  }
  targets_data_subset_norm <- t(apply(targets_norm_abundance2, 1, cal_z_score))
  
  # Create quantile breaks for heatmap for better visualization in tailed data
  mat_breaks <- seq(min(targets_norm_abundance2, na.rm=TRUE), max(targets_norm_abundance2, na.rm=TRUE), length.out = 20)
  quantile_breaks <- function(xs, n) {
    breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm=TRUE)
    breaks[!duplicated(breaks)]
  }
  mat_breaks <- quantile_breaks(targets_data_subset_norm, n = 40)
  
  # Plot heatmap
  tiff(filename = paste0("../analysis/Heatmap/Heatmap_", description, ".tiff", sep =""), #for proteins
       width = 10 * 300, 
       height = 15 * 300,
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
           annotation_names_row = F,
           annotation_names_col = F,
           drop_levels = TRUE, 
           fontsize = 10, 
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           cex=1,
           clustering_distance_rows="euclidean",
           clustering_distance_cols="euclidean",
           clustering_method="complete",
           main = gsub("\\.", " ", description),
           display_numbers = test_labels,
           number_color = "green", 
           fontsize_number = 8,
           cellheight=12, cellwidth = 30)
  
  dev.off()
  
}

protein_function(description, proteins)

