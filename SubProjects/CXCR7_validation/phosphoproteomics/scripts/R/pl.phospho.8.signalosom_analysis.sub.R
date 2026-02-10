#########Signalome Analysis
### use mysql database
### get norm abundance
### use parserNorm.pl
### get norm abundance with single phosphosite
suppressPackageStartupMessages({
  library(calibrate)
  library(dplyr)
  library(GGally)
  library(ggplot2)
  library(ggpubr)
  library(grDevices)
  library(matrixStats)
  library(network)
  library(pheatmap)
  library(PhosR)
  library(RColorBrewer)
  library(reshape2)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(viridis)
})


# 1. Load data ----
data("KinaseMotifs")
data("KinaseFamily")
data("phospho_L6_ratio_pe")
data("SPSs")
data("PhosphoSitePlus")

## get known targets of kinases
kinase_list = list("PRKACA", "PRKG1", "PRKAA1", "CDK1", "MAPKAPK2", "AKT1", "MTOR", "SRC", "CDK2", "CSNK2A1", "CDK5", "LCK")
knownTargetsList <- list()
for (i in 1:length(kinase_list)) {
  print(i)
  knownTargetsList[[i]] <- PhosphoSite.human[[kinase_list[[i]]]]
  names(knownTargetsList)[i] <- kinase_list[[i]]
}

# 2. Load all sites ----

input = list(top.filter.10, top.filter.600, top.filter.1800, 
top.filter.10.dmso.vs.0s, top.filter.600.dmso.vs.0s, top.filter.1800.dmso.vs.0s,
top.filter.10.cxcr7.vs.0s, top.filter.600.cxcr7.vs.0s, top.filter.1800.cxcr7.vs.0s)

names_input = c("10", "600", "1800", 
"10.dmso.vs.0s", "600.dmso.vs.0s", "1800.dmso.vs.0s",
"10.cxcr7.vs.0s", "600.cxcr7.vs.0s", "1800.cxcr7.vs.0s")

input2 <- input[c(1:3)]


# 3. Analysis for all time points differential sites together ----
for (i in 1:length(input2)) {
  if (i == 1 & i <= 2) { 
    top.names <- rownames(input2[[i]][input2[[i]]$PValue<0.05 & 
                                        abs(input2[[i]]$logFC)>0,]) #change for stricter filtering
  } else {
    print(as.list(test[[i]]))
    top.names <- append(top.names, rownames(input2[[i]][input2[[i]]$PValue<0.05 & 
                                                          abs(input2[[i]]$logFC)>0,]))
  }
}

top.names <- unique(top.names) #duplicates should be removed!
norm_intensity2 <- as.data.frame(norm_intensity_filter)
top.norm_intensity <- norm_intensity2[top.names,] 
top.norm_intensity <- na.omit(top.norm_intensity) #there shouldn't be any NAs, if there are, there is an incompatibility between your top names and the norm intensity, check if youre using both filtered or not.

# if you want to use only significant & differential targets:
mat.reg <- as.matrix(top.norm_intensity)

# if you want to use all the targets instead:
#mat.reg <- norm_intensity_filter

# standardize matrix
mat.std <- PhosR::standardise(mat.reg)
rownames(mat.std) <- sapply(strsplit(rownames(mat.std), ";"), function(x) { 
  gsub(" ", "", paste(toupper(x[2]), x[3], "", sep=";"))
})

## 3.1. Run PhosR kinase-substrate scoring with the default parameters ----
seqs <- sapply(strsplit(rownames(mat.reg), ";"), function(x) { gsub(" ", "", x[4]) })
kssMat <- kinaseSubstrateScore(substrate.list = PhosphoSite.human, 
                               mat = mat.std, seqs = seqs,
                               species = "human",
                               numMotif = 5, numSub = 1, verbose = FALSE)
# Write results in tables
write.table(kssMat$combinedScoreMatrix, "../analysis/signalome/kssMat_global.txt", sep = "\t", row.names = TRUE)
write.table(kssMat$ksActivityMatrix, "../analysis/signalome/kinaseactivities_global.txt", sep = "\t", row.names = TRUE)

## 3.2. Print the heatmap of the top targets ----
#choose number
number = 5
tiff(filename = paste0("../analysis/signalome/KSS", "all_top", number,".tiff"),
     width = 8 * 300, 
     height = 7* 300,
     res = 300,
     compression = "lzw")
print(kinaseSubstrateHeatmap(kssMat,top = number))
dev.off()

# 4. Separate analyses for each time point----
separate_kss <- function(i) {
  top.names <- rownames(input2[[i]][input2[[i]]$PValue<0.05 & abs(input2[[i]]$logFC)>0.5,])
  norm_intensity2 <- as.data.frame(norm_intensity_filter)
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
  seqs <- sapply(strsplit(rownames(mat.reg), ";"), function(x) { gsub(" ", "", x[4]) })
  
  kssMat <- kinaseSubstrateScore(substrate.list = PhosphoSite.human, 
                                 mat = mat.std, seqs = seqs,
                                 numMotif = 5, numSub = 1, verbose = FALSE)
  
  assign(paste0("kssMat", names_input[[i]]), kssMat, globalenv())
  
  # print the heatmap of the top targets
  tiff(filename = paste0("../analysis/signalome/KSS",names_input[i], ".tiff"),
       width = 8 * 300, 
       height = 8* 300,
       res = 300,
       compression = "lzw")
  
  print(kinaseSubstrateHeatmap(kssMat,top = 3))
  
  dev.off()
  
}
#run the function on timepoints
for (j in 1:length(input2)) {
  separate_kss(j)
}

## 4.1. Manual saving to tiff to change the dimensions ----
tiff(filename = paste0("../analysis/signalome/KSS","60.tiff"),
     width = 6 * 300, 
     height = 4* 300,
     res = 300,
     compression = "lzw")
print(kinaseSubstrateHeatmap(kssMat60,top = 3))
dev.off()


# 5. Create kinase substrate prediction ----
# PhosR uses the ‘kinaseSubstratePred’ function to synthesise the scores ‘kinaseSubstrateScore’ 
#to predict the kinase-substrate relationships using an adaptive sampling-based positive-unlabeled 
#learning method (Yang et al., 2018)

set.seed(123)
predMat <- kinaseSubstratePred(kssMat, ensembleSize = 10, top = 50, cs = 0.6, 
                               inclusion = 8, iter = 20, verbose = TRUE)

#it first sorts the kssMat$combinedScore matrix decreasing, then it checks how many
# of the top 50 (n) phosphosites under that column(kinase) have a score higher than 0.6 (cs).
# If this number is bigger than 10 (inclusion), than it makes a substrate list with 
# these substrate for each kinase. Of course we are always limited by the input substrates
# in our case either the significant ones, or the significantly differential ones. 

#defaults: 
#ensembleSize = 10 An ensemble size
#top = 50 #a number to select top kinase substrates.
#cs = 0.8 #Score threshold - for the kss matrix - to get the substrate list.
#inclusion = 20 #A minimal number of substrates required for a kinase to be selected 
#iter = 5  #A number of iterations for adaSampling

predMatSub <- as.data.frame(predMat) %>% filter_all(any_vars(. > 0.8))

tiff(filename = "../analysis/signalome/KSPredSub.tiff",
     width = 7 * 300, 
     height = 8* 300,
     res = 300,
     compression = "lzw")
pheatmap(mat = predMatSub, main = "Predicted phosphosite(s) for each kinase",
         fontsize = 7, cellwidth = 10, cellheight = 6)
dev.off()


# 6. Extract predicted targets of specific kinases ----
## 6.1. Get all significant phosphoproteins ----
for (i in 1:length(input2)) {
  df_pl2 = input[[i]]
  df_pl2 <- df_pl2[df_pl2[, "PValue"] <0.05,]
  assign(paste0("top.filter.", names_input[[i]], ".sign"), df_pl2)
}
#get the list of significant ids in all
union_sig <- Reduce(union, list(rownames(top.filter.10.sign), 
                                rownames(top.filter.600.sign),
                                rownames(top.filter.1800.sign)))
#prepare the normalized intensity data for significant phosphoproteins
dataset_df <- as.data.frame(norm_intensity_filter)
#if we want all targets and not only the significant ones, skip this line.
dataset_df <- dataset_df[rownames(dataset_df) %in% union_sig, ]
#prepare the feature to match with PhosphositePlus
name <- sapply(strsplit(rownames(dataset_df), ";"), "[[", 2)
site <- sapply(strsplit(rownames(dataset_df), ";"), "[[", 3)
dataset_df$namesite <- paste0(name,";",site,";")

## 6.2. Make the heatmaps ----
##give names of kinases

process_kinase <- function(kinase, score.t, inputMat, class){
  # Check if kinase is in the data
  if (!kinase %in% colnames(inputMat)) {
    message(paste("Kinase", kinase, "not found in the data. Skipping..."))
    return(NULL)
  }
  
  # Filter and prepare kinase scores
  kinase_scores <- data.frame(inputMat[, kinase])
  kinase_scores$name <- rownames(inputMat)
  colnames(kinase_scores) <- c("score", "id")
  kinase_scores <- kinase_scores[ kinase_scores$score > score.t,]
  kinase_scores$colors <- "black"
  
  #add known targets to kinase scores
  knownKs <- dataset_df[dataset_df$namesite %in% knownTargetsList[[kinase]],"namesite"]
  knownKdf <-  data.frame(matrix(NA, nrow =length(knownKs) , ncol = 3))
  rownames(knownKdf) <- knownKs
  colnames(knownKdf) <- c("score", "id", "colors")
  knownKdf$score <- 1
  knownKdf$colors <- "#D73027"
  knownKdf$id <- rownames(knownKdf)
   
  kinase_scores <- rbind(kinase_scores, knownKdf)
  
  # Get predicted kinase targets
  predicted_kinase_targets <- unique(kinase_scores$id)
  kinasename  <- paste0(kinase, "_", class)
  kinase_targets <- predicted_kinase_targets
  
  # Prepare target normalized abundance
  targets_norm_abundance <- dataset_df[dataset_df$namesite %in% kinase_targets,1:70]
  targets_norm_abundance <-targets_norm_abundance[order(row.names(targets_norm_abundance)),]
  targets_norm_abundance <- as.matrix(targets_norm_abundance)
  
  # Prepare targets p-values
  targets_pvals <- top.all[rownames(top.all) %in% rownames(targets_norm_abundance), c(5,8,11)]
  targets_pvals<-targets_pvals[order(row.names(targets_pvals)), ]
  targets_pvals <- as.matrix(targets_pvals)
  rownames(targets_pvals) <- sapply(strsplit(rownames(targets_pvals), ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})
  targets_pvals <- as.data.frame(targets_pvals)
  targets_pvals[targets_pvals <= 0.05] <- "*"
  targets_pvals[targets_pvals> 0.05] <- ""
  targets_pvals$adj.P.Val.00 <- ""
  targets_pvals <- targets_pvals[c(4,1:3)]
  
  # Prepare targets log fold changes
  targets_logfcs <- top.all[rownames(top.all) %in% rownames(targets_norm_abundance), c(3,6,9)]
  targets_logfcs<-targets_logfcs[order(row.names(targets_logfcs)), ]
  targets_logfcs <- as.matrix(targets_logfcs)
  rownames(targets_logfcs) <- sapply(strsplit(rownames(targets_logfcs), ";"),  
                                     function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})
  targets_logfcs <- as.data.frame(targets_logfcs)
  targets_logfcs<-round(targets_logfcs, digits = 3)
  
  # Filter log fold changes by p-values
  for (j in 1:(length(colnames(targets_pvals))-1)) {
    print(j)
    for (i in 1:length(rownames(targets_logfcs))) {
      if (targets_pvals[i, j+1] == "") {
        targets_logfcs[i, j] <- ""}}}
  
  targets_pvals$adj.P.Val.10.DMSO <- ""
  targets_pvals$adj.P.Val.600.DMSO <- ""
  targets_pvals$adj.P.Val.1800.DMSO <- ""
  targets_pvals <-  targets_pvals[c(1,5,2,6,3,7,4)]
  
  targets_logfcs$logFC.00<- ""
  targets_logfcs$logFC.10.DMSO<- ""
  targets_logfcs$logFC.600.DMSO<- ""
  targets_logfcs$logFC.1800.DMSO<- ""
  targets_logfcs <- targets_logfcs[c(4,5,1,6,2,7,3)]
  targets_logfcs<-targets_logfcs[order(row.names(targets_logfcs)), ]
  
  # prepare rownames for the heatmap
  rownames(targets_norm_abundance) <- sapply(strsplit(rownames(targets_norm_abundance), ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})
  x_tt <- as.factor(grps)
  
  # Check if more than two targets found for the kinase
  if (nrow(targets_norm_abundance) < 2) {
    message(paste("Not enough targets for kinase", kinase, "Skipping..."))
    return(NULL)
  }
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
  
  # Prepare column annotation for heatmap
  col_groups <- substr(colnames(targets_norm_abundance2), 1, 10)
  table(col_groups)
  mat_col <- data.frame(time = col_groups)
  rownames(mat_col) <- colnames(targets_norm_abundance2)
  mat_colors <- list(time = brewer.pal(7, "BuPu"))
  names(mat_colors$time) <- unique(col_groups)
  
  
  # Prepare row annotation for heatmap
  mat_rows <- kinase_scores %>%
    filter(rownames(kinase_scores) %in% dataset_df$namesite) %>%
    mutate(namesite = id) %>%
    dplyr::select(namesite, score) %>%
    inner_join(dataset_df %>% mutate(id = rownames(.)) %>% dplyr::select(id, namesite), by = "namesite") %>%
    mutate(new_id = sapply(strsplit(id, ";"), function(x) { paste(x[[1]], x[[2]], x[[3]], sep = "_") })) %>%
    dplyr::select(new_id, score) %>%
    tibble::column_to_rownames("new_id")
  #row_colors <- list(score = c("#FFFAE5","#FFF6CC","#FFE566",
   #                            "#FFBF00","#FF8000","#A45511"))
  row_colors <- list(score =colorRampPalette(rev(brewer.pal(n = 7, 
                                                        name = "RdYlBu")))(50)[c(25, 30, 35,40, 50)])
  # Prepare row text annotation
  mat_text_colors <- kinase_scores %>%
    filter(rownames(kinase_scores) %in% dataset_df$namesite) %>%
    mutate(namesite = id) %>%
    dplyr::select(namesite, colors) %>%
    inner_join(dataset_df %>% mutate(id = rownames(.)) %>% dplyr::select(id, namesite), by = "namesite") %>%
    mutate(new_id = sapply(strsplit(id, ";"), function(x) { paste(x[[1]], x[[2]], x[[3]], sep = "_") })) %>%
    dplyr::select(new_id, colors) %>%
    tibble::column_to_rownames("new_id")

  #get the labels
  test_labels <- as.matrix(targets_logfcs) 
  
  # Calculate z-scores
  cal_z_score <- function(x){
    (x - mean(x)) / sd(x)
  }
  targets_data_subset_norm <- t(apply(targets_norm_abundance2, 1, cal_z_score))  ##applies to rows

  # Create quantile breaks for heatmap for better visualization in tailed data
  mat_breaks <- seq(min(targets_norm_abundance2, na.rm=TRUE), max(targets_norm_abundance2, na.rm=TRUE), length.out = 20)
  quantile_breaks <- function(xs, n) {
    breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm=TRUE)
    breaks[!duplicated(breaks)]
  }
  mat_breaks <- quantile_breaks(targets_data_subset_norm, n = 40)
  
  e <- pheatmap(mat = targets_data_subset_norm, 
                color = colorRampPalette(c("navy", "white", "red"))(length(mat_breaks) - 1),
                gaps_row=c(5,10,15,20,25,30,35,40,45,50, 55),
                scale="row",
                na_col = "grey",
                breaks = mat_breaks,
                border_color = "white", 
                show_colnames = TRUE, 
                show_rownames = TRUE, 
                annotation_row = mat_rows,
                annotation_col = mat_col, 
                annotation_colors = c(mat_colors,row_colors),
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
                main = paste0(kinase, " - ", class, " targets (st:", score_threshold, ")"),
                display_numbers = test_labels,
                number_color = "green", 
                fontsize_number = 8,
                cellheight=12, cellwidth = 30
  )
  
  cols=mat_text_colors[order(match(rownames(mat_text_colors), e$gtable$grobs[[5]]$label)), ]
  
  e$gtable$grobs[[5]]$gp=gpar(col=cols)
  
  # Plot heatmap
  tiff(filename = paste0("../analysis/signalome/Heatmap_", kinasename, "_sign.tiff", sep =""),
       width = 10 * 300, 
       height = 20 * 300,
       res = 300,
       compression = "lzw")
  
  print(e)
  dev.off()

}

# Process each kinase in the list --> directly creates the heatmaps
kinase_list = list("PRKACA", "PRKG1", "PRKAA1", "CDK1", "MAPKAPK2", "AKT1", "MTOR", "SRC", "CDK2", "CSNK2A1", "CDK5")
inputMats <- list(kssMat$combinedScoreMatrix, predMat)  # KSScomb or predMat
classes <-  c("KSS", "PRED") #KSS or PRED
score_threshold <- 0.8

lapply(kinase_list, process_kinase, score_threshold, inputMats[[2]], classes[2])


# 7. Constructing Signalling Networks (Signalomes) ----
kinaseOI = c("PRKACA", "PRKG1")
signalomesRes <- Signalomes(KSR = kssMat, 
                            predMatrix = predMat, 
                            exprsMat = mat.std, 
                            module_res = 10,
                            KOI = kinaseOI)

# generate palette
my_color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))
kinase_all_color <- my_color_palette(ncol(inputMat))
names(kinase_all_color) <- colnames(inputMat)
kinase_signalome_color <- kinase_all_color[colnames(predMat)]

plotSignalomeMap(signalomes = signalomesRes,
                 color = kinase_signalome_color)


# plot the signalome network that illustrates the connectivity between kinase signalome networks.
plotKinaseNetwork(KSR = kssMat, 
                  predMatrix = predMat, 
                  threshold = 0.9, 
                  color = kinase_all_color)




