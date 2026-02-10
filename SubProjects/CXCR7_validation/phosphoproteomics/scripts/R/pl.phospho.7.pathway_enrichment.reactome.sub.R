############################################
#### script for annotation of phosphoproteom data 
## GO and KEGG enrichtment

BiocManager::install("org.Hs.eg.db")
BiocManager::install("reactome.db")
install.packages("sjmisc")
remotes::install_github("JosephCrispell/basicPlotteR")

# Load packges ----
suppressPackageStartupMessages({
  library(annotate)
  library("basicPlotteR")
  library(calibrate)
  library(clusterProfiler)
  library(cowplot)
  library(directPA)
  require(dplyr)
  library(enrichplot)
  library(ggplot2)
  library(limma)
  library(org.Hs.eg.db)
  library(pheatmap)
  library(PhosR)
  library(plyr)
  library(RColorBrewer)
  library(reactome.db)
  library(ReactomePA)
  library(remotes)
  library(rlist)
  library('sjmisc')
  library(stringr)
  library(tidyr)
})



##################################################################
# Rectome enrichment ----
## accroding to https://pyanglab.github.io/PhosR/articles/PhosR.html
#input = list(top.collapse.10, top.collapse.30, top.collapse.60, top.collapse.300, top.collpase.600, 
# top.collapse.900, top.collapse.1800)

## define input  ----
input = list(top.filter.10, top.filter.600, top.filter.1800, 
             top.filter.10.dmso.vs.0s, top.filter.600.dmso.vs.0s, top.filter.1800.dmso.vs.0s,
             top.filter.10.cxcr7.vs.0s, top.filter.600.cxcr7.vs.0s, top.filter.1800.cxcr7.vs.0s)

names_input = c("10", "600", "1800", 
                "10.dmso.vs.0s", "600.dmso.vs.0s", "1800.dmso.vs.0s",
                "10.cxcr7.vs.0s", "600.cxcr7.vs.0s", "1800.cxcr7.vs.0s")

for (i in 1:length(input)) {
  input[[i]] <- input[[i]][order(row.names(input[[i]])),]
}

Tc <- as.matrix(cbind(input[[1]]$logFC, input[[2]]$logFC, input[[3]]$logFC,
                           input[[4]]$logFC, input[[5]]$logFC, input[[6]]$logFC,
                           input[[7]]$logFC, input[[8]]$logFC, input[[9]]$logFC))

rownames(Tc) <- rownames(input[[1]])
colnames(Tc) <- names_input


## define input2 ----
input = list(top.collapse.10, top.collapse.600, top.collapse.1800, 
             top.collapse.10.dmso.vs.0s, top.collapse.600.dmso.vs.0s, top.collapse.1800.dmso.vs.0s,
             top.collapse.10.cxcr7.vs.0s, top.collapse.600.cxcr7.vs.0s, top.collapse.1800.cxcr7.vs.0s)

names_input = c("10", "600", "1800", 
                "10.dmso.vs.0s", "600.dmso.vs.0s", "1800.dmso.vs.0s",
                "10.cxcr7.vs.0s", "600.cxcr7.vs.0s", "1800.cxcr7.vs.0s")


for (i in 1:length(input)) {
  input[[i]] <- input[[i]][order(row.names(input[[i]])),]
}

names(input) <- names_input


# Tc.gene <- as.matrix(cbind(input[[1]]$logFC, input[[1]]$PValue, 
#                            input[[2]]$logFC,input[[2]]$PValue,
#                            input[[4]]$logFC, input[[3]]$PValue))
# colnames(Tc.gene) <- c("logFC.10","adjp.10","logFC.600","adjp.600","logFC.1800","adjp.1800" )
# rownames(Tc.gene) <- paste(sapply(strsplit(rownames(input[[1]]), ";"), "[[", 1))
# write.table(Tc.gene, file = "../data/processed_data/top.all.collapse.txt", 
#             sep = '\t', row.names = TRUE)


Tc.gene <- as.matrix(cbind(input[[1]]$logFC, input[[2]]$logFC, input[[3]]$logFC,
                          input[[4]]$logFC, input[[5]]$logFC, input[[6]]$logFC,
                          input[[7]]$logFC, input[[8]]$logFC, input[[9]]$logFC))

rownames(Tc.gene) <- paste(sapply(strsplit(rownames(input[[1]]), ";"), "[[", 2))
colnames(Tc.gene) <- names_input

#rownames(Tc.gene) <- paste(sapply(strsplit(rownames(input[[1]]), ";"), "[[", 1))
#write.table(Tc.gene, file = "../data/processed_data/top.all.collapse.txt", 
#            sep = '\t', row.names = TRUE)

## select which dataset to choose ----
pathways = as.list(reactomePATHID2EXTID)
path_names = as.list(reactomePATHID2NAME)
name_id = match(names(pathways), names(path_names))
names(pathways) = unlist(path_names)[name_id]

pathways = pathways[which(grepl("Homo sapiens", names(pathways), ignore.case = TRUE))]

pathways = lapply(pathways, function(path) {
  gene_name = unname(getSYMBOL(path, data = "org.Hs.eg"))
  toupper(unique(gene_name))
})

## select up or downregulation ----
class = "UP/more"
order = "greater"
colors <- c("#FDAE61", "red")

class = "DOWN/less"
order = "less"
colors <- c("lightblue2", "blue")

i = 2 # 7,8,9
##  ENRICHED REACTOME PATHWAYS; MAKE FIGURES; SAVE TO TIFF ----
for (i in 1:length(input)) {
  # geneSet <- names(sort(Tc.gene[,i],  #decreasing False for enrichment of downregulated
  #                       decreasing = T))[seq(round(nrow(Tc.gene) * 0.05))]
  #geneSet <- names(sort(TcP.gene[,i], 
  #                       decreasing = FALSE))[TcP.gene[,i]<0.05]
  #geneSet <- input[[1]]$name[input[[4]]$PValue<0.05]
  
  # head(geneSet)
  
  
  # path1 <- pathwayOverrepresent(geneSet, annotation=pathways, 
  #                               universe = rownames(Tc.gene), alter = order)
  path2 <- pathwayRankBasedEnrichment(Tc.gene[,i], 
                                      annotation=pathways, 
                                      alter = order)
  
  # lp1 <- -log10(as.numeric(path1[names(pathways),1]))
  lp2 <- -log10(as.numeric(path2[names(pathways),1]))
  
  
  ######### Visualization1
  #dev.new()
  # tiff(filename = paste0("../analysis/Reactome_enrichment/", class,  names_input[[i]] , "_plot.tiff"),
  #      width = 5 * 300,
  #      height = 5 * 300,
  #      res = 300,
  #      compression = "lzw")
  # #png(file=paste0("../analysis/reactome_enrichment/", "reactome_top.10",".png"), width = 5, height = 5, units = 'in', res = 1200)
  # #par(mfrow=c(1,1))
  # plot(lp1, lp2, xlab="Overrepresentation (-log10 pvalue)", ylab="Rank-based enrichment (-log10 pvalue)", main="", xlim = c(0, 5), 
  #      ylim = c(0, 2.6), col = "grey", bg = "grey", )
  # 
  # # select highly enriched pathways
  # sel <- which(lp1 > 1.3 & lp2 > 1.3)
  # #textxy(lp1[sel], lp2[sel], gsub("_", " ", gsub("REACTOME_|Mus musculus: ", "", names(pathways)))[sel], cex= 0.8,srt=0)
  # label <- gsub("_", " ", gsub("REACTOME_|Homo sapiens: ", "", names(pathways)))[sel]
  # label <- word_wrap(label, 45, linesep = NULL)
  # addTextLabels(lp1[sel], lp2[sel], label, cex.label=0.9, col.label="black", lty=2, col.line=rgb(0,0,0, 0.5),keepLabelsInside = TRUE)
  # dev.off()
  # 
  
  ######### Visualization2
  pathways2 <- as.data.frame(t(data.frame(lapply(pathways, function(x) length(x)))))
  path3 <- as.data.frame(path2)
  path4 <- merge(path3, pathways2, by = "row.names", all = FALSE)
  colnames(path4) <- c("pathway", "pvalue", "number.substrates", "substrates", "pw.size")
  path4 <- na.omit(path4)
  path4$number.substrates <- as.numeric(path4$number.substrates)
  path4$pw.size <- as.numeric(path4$pw.size)
  path4$pvalue <- as.numeric(path4$pvalue)
  path4$ratio <- path4$number.substrates/path4$pw.size
  path4$pathway <- gsub("_", " ", gsub("REACTOME_|Homo.sapiens..", "", path4$pathway))
  #path4 <- path4[order(path4$pvalue),]
  path4 <- path4[path4$pvalue<0.05,]
  path4 <- path4[order(path4$ratio, decreasing = T),]
  write.table(path4, file = paste0("../analysis/Reactome_enrichment/", class, names_input[[i]] ,"_pathways.txt"), 
              sep = '\t', row.names = FALSE)
  assign(paste0("path4.", order, names_input[[i]]), path4)
  #sel = c("2082", "798","758", "1254","910", "1209", "2262","2462", "2084", "1094")
  #path4 <- path4[sel,]
  #path4 <- path4[order(path4$ratio, decreasing = T),]
  #path4 <- path4[path4$number.substrates > 9,]
  path4 <- path4[c(2:4, 6:8,14,20,22,35),] # DOWN10: c(1,3,6,8,12,16,19,20,28,44), UP600:c(7,9,11,14,15,18,21,27,29,34), DOWN600:c(2:5,8,9,12,16,17,20), UP1800:c(1,2,6,15,21,25,30,33,34,37)
  #path4 <- path4[1:10,]
  path4$pvalue <- round(path4$pvalue, digits=4)
  path4 <- path4[order(path4$pvalue),]
  path4$pathway <- gsub("\\.", " ", path4$pathway)
  
  #dev.new()
  
  ggp <- ggplot(path4[10:1,], aes(x=1:10,y=ratio,fill=pvalue)) + 
    coord_flip() +
    geom_bar(stat = "identity") + 
    scale_fill_gradient(low = colors[1], high = colors[2] ,
      limits = c(min(path4[1:10,]$pvalue), max(path4[1:10,]$pvalue)), 
      name = "pvalue", 
      guide = guide_colorbar(barwidth = , 
                             barheight = 10, 
                             title.position = "top", 
                             title.hjust = 0.5,
                             reverse = TRUE))  +
    geom_text(aes(label = pathway, y = 0.005),
              vjust = 0, colour = "black",
              position = position_dodge(0.1), size =8, hjust = 'left') +
    labs(title = paste0(names_input[i]," sec"), size = 25, x = "pathways", y = "gene ratio") +
    theme_cowplot() +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          axis.text.x=element_text(size=20),
          axis.title.x = element_text(size =25),
          axis.title.y = element_text(size = 25)
    )
  tiff(filename = paste0("../analysis/Reactome_enrichment/", class, names_input[[i]] ,"_barplot.tiff"),
       width = 12 * 300,
       height = 8 * 300,
       res = 300,
       compression = "lzw")

  print(ggp)
  dev.off()
}


# Heatmaps of chosen pathways
protein_function <- function(description, proteins, timepoint) {
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
  tiff(filename = paste0("../analysis/Reactome_enrichment/UP/Heatmap.", timepoint, ".", description, ".tiff", sep =""),
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

pathway.list <- list("2262","2138", "1976", "2082", "426", "1209", "899", "1190", "1177")

timepointOI <- "10"
dfI <- path4.greater10
pathway.list <- as.list(rownames(dfI))

for (i in 1:length(pathway.list)) {
  pathOI <- dfI[pathway.list[[i]],]
  ids <- pathOI$substrates %>% strsplit(.,";")
  protein <- paste(ids[[1]], collapse = ";|")

  pm <- median(Tc.gene[ids[[1]],timepointOI])
  nm <- median(Tc.gene[!rownames(Tc.gene) %in% ids[[1]], timepointOI])
  
  if (pm > nm) {
    print(paste0("Median logFC in ",gsub("\\.", " ", pathOI$pathway), " is ", round(pm, digits = 2), 
                 " and higher than the median of the rest (", round(nm,digits = 2), ")"))
  } else {
    print("Lower pathway median, double check!")
  }

  # get ids
  Tc1 <- data.frame(Tc[str_detect(row.names(Tc), protein),])
  Tc2 <- data.frame(Tc.gene)
  col = paste0("X", timepointOI)
  proteinList <- rownames(Tc1[Tc1[col] %in% Tc2[col]])

  # make heatmap
  protein_function(pathOI$pathway, proteinList, timepointOI)
}


# easier and manual heatmaps of chosen pathways ----
library(pheatmap)
names(pathways) <- gsub("Homo sapiens: ", "", names(pathways))

HeatmapOfReactome <- function(description, proteins,
                              cellW = 24, cellH = 12, 
                              showRowNames = TRUE, showlogFCs = test_labels,
                              border_color = "white")  {
  # Get normalized intensities for proteins of interest
  targets_norm_abundance <- dataset_df[rownames(dataset_df) %in% proteins, 1:70]
  targets_norm_abundance <-targets_norm_abundance[order(row.names(targets_norm_abundance)), ]
  targets_norm_abundance <- as.matrix(targets_norm_abundance)
  
  # Get p-values for their expression at 3 main timepoints
  targets_pvals <- top.all.input[rownames(top.all.input) %in% rownames(targets_norm_abundance), c(5,8,11)] #5,8,11 for adj.
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
  colnames(targets_norm_abundance.collapse) <- c("DMSO00","DMSO10","ACKR310","DMSO600","ACKR3600","DMSO1800","ACKR31800")
  targets_norm_abundance.collapse <- as.data.frame(targets_norm_abundance.collapse)
  
  targets_norm_abundance2 <- targets_norm_abundance.collapse
  # Prepare annotation for heatmap
  col_groups <- substr(colnames(targets_norm_abundance2), 1, 11)
  mat_col <- data.frame(time = col_groups)
  rownames(mat_col) <- colnames(targets_norm_abundance2)
  mat_colors <- list(time = brewer.pal(7, "BuPu"))
  names(mat_colors$time) <- unique(col_groups)
  test_labels <- as.matrix(targets_logfcs) 
  
  # Calculate z-scores
  # cal_z_score <- function(x){
  #   (x - mean(x)) / sd(x)
  # }
  #targets_data_subset_norm <- t(apply(targets_norm_abundance2, 1, cal_z_score))
  targets_data_subset_norm <- targets_norm_abundance2 #pheatmap does already z-score calculation
  
  # Create quantile breaks for heatmap for better visualization in tailed data #creates too much color difference.
  # mat_breaks <- seq(min(targets_norm_abundance2, na.rm=TRUE), max(targets_norm_abundance2, na.rm=TRUE), length.out = 20)
  # quantile_breaks <- function(xs, n) {
  #   breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm=TRUE)
  #   breaks[!duplicated(breaks)]
  # }
  # mat_breaks <- quantile_breaks(targets_data_subset_norm, n = 40)
  
  heatmap <- pheatmap(mat = targets_data_subset_norm, 
                      color = colorRampPalette(c("blue", "white", "red"))(100),
                      #gaps_row=c(5,10,15,20,25,30,35,40,45,50, 55),
                      scale="row",
                      na_col = "grey",
                      #breaks = mat_breaks,
                      border_color = border_color, 
                      show_colnames = TRUE, 
                      show_rownames = showRowNames, 
                      annotation_col = mat_col,
                      annotation_colors = mat_colors, 
                      annotation_names_row = F,
                      annotation_names_col = F,
                      drop_levels = TRUE, 
                      fontsize = 10, 
                      cluster_cols = FALSE,
                      cluster_rows = TRUE,
                      cex=1,
                      clustering_distance_rows="correlation",
                      clustering_distance_cols="euclidean",
                      clustering_method="complete",
                      main = gsub("\\.", " ", description),
                      display_numbers = showlogFCs,
                      number_color = "green", 
                      fontsize_number = 8,
                      cellheight=cellH, cellwidth = cellW)
  return(heatmap)
  
  
  
}

# or reactome terms
UPtermsOI <-  list("Signaling by ERBB2","cGMP effects", "Interleukin-15 signaling", 
                 "DAP12 signaling", "Intrinsic Pathway for Apoptosis",
                 "Sphingolipid metabolism", "GP1b-IX-V activation signalling",
                 "RHOA GTPase cycle", "RHOV GTPase cycle",
                 "Signalling to ERKs", "Negative regulation of MAPK pathway",
                 "Cell death signalling via NRAGE, NRIF and NADE",
                 "Smooth Muscle Contraction", "COPII-mediated vesicle transport",
                 "Negative regulation of the PI3K/AKT network",
                 "Amyloid fiber formation",
                 "Signaling by Nuclear Receptors", 
                 "FLT3 signaling in disease",
                 "Sensory processing of sound",
                 "Signaling by NTRK3 (TRKC)",
                 "Mitophagy",
                 "Factors involved in megakaryocyte development and platelet production",
                 "Peptide hormone metabolism",
                 "RAF activation",
                 "Constitutive Signaling by Overexpressed ERBB2",
                 "Signaling by MET",
                 "Signaling by Hippo"
                 )

DOWNtermsOI <- list("AKT phosphorylates targets in the cytosol", "VEGFA-VEGFR2 Pathway",
                    "Generation of second messenger molecules", "Effects of PIP2 hydrolysis",
                    "L1CAM interactions", "Signaling by PTK6", 
                    "Plasma lipoprotein assembly, remodeling, and clearance",
                    "GRB2:SOS provides linkage to MAPK signaling for Integrins",
                    "Signaling by MET", "RND2 GTPase cycle", "Regulation of RAS by GAPs",
                    "Signaling by Hedgehog", "Degradation of GLI1 by the proteasome",
                    "Regulation of Apoptosis", "Integrin cell surface interactions",
                    "ABC-family proteins mediated transport", 
                    "GSK3B and BTRC:CUL1-mediated-degradation of NFE2L2",
                    "Signaling by PDGFR in disease", "MTOR signalling",
                    "Rab regulation of trafficking", "Ca-dependent events",
                    "Energy dependent regulation of mTOR by LKB1-AMPK",
                    "Clathrin-mediated endocytosis", "Signaling by FGFR1 in disease",
                    "Autophagy", "TAK1-dependent IKK and NF-kappa-B activation",
                    "G alpha (q) signalling events")

#names(pathways)[grep("alpha", names(pathways))]

termsOI <- DOWNtermsOI
i = 1
for (i in seq_along(termsOI)) {
  symbols <- pathways[[which(names(pathways) == termsOI[i])]]
  protein <- paste(symbols, collapse = ";|;")
  description <- termsOI[[i]]
  
  # Get all phosphopeptides 
  proteins <- rownames(dataset_df[str_detect(row.names(dataset_df), protein),])
  
  tiff(filename = paste0("../analysis/Reactome_enrichment/Heatmaps/all/", gsub("[/ -:]", "_",description), "_all.tiff"),
       width = 10 * 300, 
       height = 14 * 300,
       res = 300,
       compression = "lzw")
  HeatmapOfReactome(description, proteins, cellH = 12)
  dev.off()
  
  tiff(filename = paste0("../analysis/Reactome_enrichment/Heatmaps/allCondensed/", gsub("[/ -:]", "_",description), "_allC.tiff"),
       width = 10 * 300, 
       height = 14 * 300,
       res = 300,
       compression = "lzw")
  HeatmapOfReactome(description, proteins, cellW = 12, cellH = 4, showRowNames = F, showlogFCs = F, 
                    border_color = NA)
  dev.off()
  
  # Or get significant proteins
  proteins <- union_sig_df[union_sig_df$protein %in% symbols, "ids"] # merged2 for site-centric
  
  # Check if 'proteins' is empty, if so, skip the rest of the iteration
  if (length(proteins) < 2) {
    next  # Skip to the next iteration if 'proteins' is empty
  }
  
  # Create tiff file only if 'proteins' is not empty
  tiff(filename = paste0("../analysis/Reactome_enrichment/Heatmaps/significant/", gsub("[/ -:]", "_",description), ".tiff"),
       width = 10 * 300, 
       height = 12 * 300,
       res = 300,
       compression = "lzw")
  HeatmapOfReactome(description, proteins)
  dev.off()
}

# get enriched reactome for up or down at timepoints (ORA) ----
library(ReactomePA)
library(dplyr)
library(tidyr)


description <-  "1800_DOWN"
proteins <- rownames(top.1800.sign[top.1800.sign$logFC > 1, ])
HeatmapOfReactome(description, proteins, cellH = 6, showlogFCs = F,
                  border_color = NA)


proteins_entrez <- bitr(unique(sapply(strsplit(proteins, ";"), "[[", 1)), 
                        fromType = "UNIPROT",
                        toType = "ENTREZID", 
                        OrgDb =  "org.Hs.eg.db")
proteins_entrez <- proteins_entrez[!duplicated(proteins_entrez$UNIPROT), ]
proteins_entrez <- unique(proteins_entrez$ENTREZID)

# get universe
universe_entrez <- bitr(unique(top.collapse$uniprot_id), 
                        fromType = "UNIPROT",
                        toType = "ENTREZID", 
                        OrgDb =  "org.Hs.eg.db")
universe_entrez <- universe_entrez[!duplicated(universe_entrez$UNIPROT), ]
universe_entrez <- unique(universe_entrez$ENTREZID)

# reactome ORA
Reactome <- enrichPathway(proteins_entrez, 
                          universe=universe_entrez,
                          organism = "human",
                          pvalueCutoff = 0.05, 
                          pAdjustMethod = "BH",
                          minGSSize = 10,
                          maxGSSize = 500,
                          readable=TRUE)

rea.tab = Reactome@result

rea.tab <- rea.tab[rea.tab$Count < 11,] # to get too general pathways out
rea.tab = rea.tab[rea.tab$pvalue < 0.05,]
# Step 1: Split the geneID column into individual genes
expanded_rea <- rea.tab %>%
  separate_rows(geneID, sep = "/")

# Step 2: For each gene, select the pathway with the lowest p-value
Final_rea_df <- expanded_rea %>%
  group_by(geneID) %>%
  slice_min(pvalue, with_ties = FALSE) %>%
  ungroup() %>%
  select(geneID, Description, pvalue) # Adjust selection as per requirement


# MAKE HEATMAP----

# Get normalized intensities for proteins of interest and convert to matrix
targets_norm_abundance <- dataset_df[rownames(dataset_df) %in% proteins, 1:70] %>%
  .[order(rownames(.)), ] %>%
  as.matrix()

# Get p-values for their expression at 3 main timepoints
targets_pvals <- top.all.input[rownames(top.all.input) %in% rownames(targets_norm_abundance), c(5,8,11)] #5,8,11 for adj.
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


targets_logfcs$logFC.00<- ""
targets_logfcs$logFC.10.DMSO<- ""
targets_logfcs$logFC.600.DMSO<- ""
targets_logfcs$logFC.1800.DMSO<- ""
targets_logfcs <- targets_logfcs[c(4,5,1,6,2,7,3)]

#order
targets_logfcs<-targets_logfcs[order(row.names(targets_logfcs)), ]
test_labels <- as.matrix(targets_logfcs) 
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
colnames(targets_norm_abundance.collapse) <- c("DMSO00","DMSO10","ACKR310","DMSO600","ACKR3600","DMSO1800","ACKR31800")
targets_norm_abundance.collapse <- as.data.frame(targets_norm_abundance.collapse)

targets_norm_abundance2 <- targets_norm_abundance.collapse

# Prepare col annotation for heatmap
col_groups <- substr(colnames(targets_norm_abundance2), 1, 11)
mat_col <- data.frame(time = col_groups)
rownames(mat_col) <- colnames(targets_norm_abundance2)
mat_colors <- list(time = brewer.pal(7, "BuPu"))
names(mat_colors$time) <- unique(col_groups)


# Prepare row annotation for heatmap
id_df <- data.frame(id = rownames(targets_norm_abundance2), 
                    name = sapply(strsplit(rownames(targets_norm_abundance2), "_"),"[[", 2))

pathwaydf <- Final_rea_df[Final_rea_df$geneID %in% id_df$name,1:2]
colnames(pathwaydf) <- c("name", "pathway")

id_df <- left_join(id_df, pathwaydf)
rownames(id_df) <- id_df$id
id_df <- id_df[,3,drop=FALSE]
id_df <- id_df[order(id_df$pathway),1,drop = F]


targets_data_subset_norm <- targets_norm_abundance2 #pheatmap does already z-score calculation

heatmap <- pheatmap(mat = targets_data_subset_norm[rownames(id_df),], 
                    color = colorRampPalette(c("darkslategray", "white", "darkorange"))(100),
                    #gaps_row=c(5,10,15,20,25,30,35,40,45,50, 55),
                    scale="row",
                    na_col = "grey",
                    #breaks = mat_breaks,
                    border_color = NA, 
                    show_colnames = TRUE, 
                    show_rownames = T, 
                    annotation_col = mat_col,
                    annotation_row = id_df,
                    annotation_colors = mat_colors, 
                    annotation_names_row = F,
                    annotation_names_col = F,
                    drop_levels = TRUE, 
                    fontsize = 10, 
                    cluster_cols = FALSE,
                    cluster_rows = FALSE,
                    cex=1,
                    #clustering_distance_rows="correlation",
                    #clustering_distance_cols="euclidean",
                    #clustering_method="complete",
                    main = gsub("\\.", " ", description),
                    display_numbers = F,
                    number_color = "green", 
                    fontsize_number = 8,
                    cellheight=10, cellwidth = 20)


# save to image ----
tiff(filename = paste0("../analysis/Reactome_enrichment/ORA_", description, ".tiff"),
     width = 10 * 300, 
     height = 14 * 300,
     res = 300,
     compression = "lzw")
print(heatmap)
dev.off()









#SNIPPETS####
#### extraxt certain pathways
rea1 <- pathways["Homo sapiens: Signaling by ROBO receptors"]
path2["Homo sapiens: Signaling by ROBO receptors"]
rea1 <- rea1[rea1 %in% geneSet]
write.table(rea1, "../../analysis/Reactome_enrichment/rea1.txt", sep="\t",  row.names=FALSE, quote=FALSE)

#### extract certain proteins
top.all[grep("ARHGA", rownames(top.all)),]



barplot(path3$number ~ path3$pvalue, horiz = TRUE)


pl_bar <- barplot(go_enrich_pl, 
                  drop = TRUE, 
                  x= "GeneRatio",
                  showCategory = 10, 
                  title = paste(desc, "GO Overrepresentation", sep=" "),
                  font.size = 18,
)



###### 2D direction site-centric kinase activity analyses
data("PhosphoSitePlus")
data("phospho_L6_ratio_pe")
data("SPSs")

Tc1 <- Tc
rownames(Tc1) <- paste(sapply(strsplit(rownames(Tc), ";"), "[[", 2), sapply(strsplit(rownames(Tc), ";"), "[[", 3), "", sep = ";")
rownames(Tc1) <- lapply(rownames(Tc1), function(x){gsub(";[STY]", ";", x)})

tiff(filename = paste0("../analysis/DPA/DPA_cxcr7vsdmso.tiff"),
     width = 15 * 300, 
     height = 5 * 300,
     res = 300,
     compression = "lzw")

par(mfrow=c(1,3))
dpa1 <- directPA(Tc1[,c(1,2)], direction=0, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dpa2 <- directPA(Tc1[,c(1,3)], direction=0, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dpa3 <- directPA(Tc1[,c(2,3)], direction=0, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dev.off()

tiff(filename = paste0("../analysis/DPA/DPA_dmsovs0.tiff"),
     width = 15 * 300, 
     height =5 * 300,
     res = 300,
     compression = "lzw")

par(mfrow=c(1,3))
dpa4 <- directPA(Tc1[,c(4,5)], direction=0, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dpa5 <- directPA(Tc1[,c(4,6)], direction=0, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dpa6 <- directPA(Tc1[,c(5,6)], direction=0, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dev.off()


tiff(filename = paste0("../analysis/DPA/DPA_cxcr7vs0.tiff"),
     width = 15 * 300, 
     height = 5 * 300,
     res = 300,
     compression = "lzw")

par(mfrow=c(1,3))
dpa7 <- directPA(Tc1[,c(7,8)], direction=0, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dpa8 <- directPA(Tc1[,c(7,9)], direction=0, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dpa9 <- directPA(Tc1[,c(8,9)], direction=0, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dev.off()




##### match psites in data and in phophositeplus
library(rlist)
test <- lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)})
for (i in 1:length(test)) {
  if (i == 1 & i <= 2) { 
    test2 = as.list(test[[i]])
  } else {
    print(as.list(test[[i]]))
    test2 <- append(test2, as.list(test[[i]]))
  }
}

which(rownames(Tc1) %in% test2)

sites.psite <- as.data.frame(t(as.data.frame(test2)))
sites.data <- as.data.frame(rownames(Tc1))
colnames(sites.psite) <- "psite"
colnames(sites.data) <- "psite"

which(sites.psite$psite %in% sites.data$psite)

sites.data$psite <- str_trim(sites.data$psite, side = c("both"))
sites.psite$psite <- str_trim(sites.psite$psite, side = c("both"))

########### top activated kinases
dpa1$pathways[1:5,]

dpa2$pathways[1:5,]

dpa3$pathways[1:5,]




tiff(filename = paste0("../analysis/DPA/Kinase_pert_cxcr7vsdmso.tiff"),
     width = 15 * 300, 
     height = 5 * 300,
     res = 300,
     compression = "lzw")

par(mfrow=c(1,3))
z1 <- perturbPlot2d(Tc1[,c(1,2)], cex=1, xlim=c(-4, 4), ylim=c(-4, 4),  
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                    main="Kinase perturbation analysis")
z2 <- perturbPlot2d(Tc1[,c(1,3)], cex=1, xlim=c(-4, 4), ylim=c(-4, 4),  
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                    main="Kinase perturbation analysis")
z3 <- perturbPlot2d(Tc1[,c(2,3)], cex=1, xlim=c(-4, 4), ylim=c(-4, 4),  
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                    main="Kinase perturbation analysis")
dev.off()



tiff(filename = paste0("../analysis/DPA/Kinase_pert_dmsovs0.tiff"),
     width = 15 * 300, 
     height = 5 * 300,
     res = 300,
     compression = "lzw")

par(mfrow=c(1,3))
z4 <- perturbPlot2d(Tc1[,c(4,5)], cex=1, xlim=c(-4, 4), ylim=c(-4, 4),  
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                    main="Kinase perturbation analysis")
z5 <- perturbPlot2d(Tc1[,c(4,6)], cex=1, xlim=c(-4, 4), ylim=c(-4, 4),  
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                    main="Kinase perturbation analysis")
z6 <- perturbPlot2d(Tc1[,c(5,6)], cex=1, xlim=c(-4, 4), ylim=c(-4, 4),  
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                    main="Kinase perturbation analysis")
dev.off()

tiff(filename = paste0("../analysis/DPA/Kinase_pert_cxcr7vs0.tiff"),
     width = 15 * 300, 
     height = 5 * 300,
     res = 300,
     compression = "lzw")

par(mfrow=c(1,3))
z7 <- perturbPlot2d(Tc1[,c(7,8)], cex=1, xlim=c(-4, 4), ylim=c(-4, 4),  
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                    main="Kinase perturbation analysis")
z8 <- perturbPlot2d(Tc1[,c(7,9)], cex=1, xlim=c(-4, 4), ylim=c(-4, 4),  
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                    main="Kinase perturbation analysis")
z9 <- perturbPlot2d(Tc1[,c(8,9)], cex=1, xlim=c(-4, 4), ylim=c(-4, 4),  
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                    main="Kinase perturbation analysis")
dev.off()

tiff(filename = paste0("../analysis/DPA/Kinase_pert_dmsovscxcr7fcs.tiff"),
     width = 15 * 300, 
     height = 5 * 300,
     res = 300,
     compression = "lzw")

par(mfrow=c(1,3))
z7 <- perturbPlot2d(Tc1[,c(4,7)], cex=1, xlim=c(-4, 4), ylim=c(-4, 4),  
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                    main="Kinase perturbation analysis")
z8 <- perturbPlot2d(Tc1[,c(5,8)], cex=1, xlim=c(-4, 4), ylim=c(-4, 4),  
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                    main="Kinase perturbation analysis")
z9 <- perturbPlot2d(Tc1[,c(6,9)], cex=1, xlim=c(-4, 4), ylim=c(-4, 4),  
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                    main="Kinase perturbation analysis")
dev.off()

###### 3D direction site-centric kinase activity analyses
z10 <- perturbPlot3d(Tc1[,c(1,2,3)], annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                     main="Kinase perturbation analysis", cex=1, xlim=c(-2, 4), ylim=c(-2, 4), zlim=c(-2, 4))


###### 2D direction gene-centric pathway analyses
data("PhosphoSitePlus")
data("phospho_L6_ratio_pe")
data("SPSs")

# load the proteomics dataset
data(PM)

# load pathway annotations
data(Pathways)


# direction pathway analysis in 3-dimensional space. Implemnted as rotating by contrast
# (1) test combined effect of all 3 treatments (stimulation and inhibitions) vs control (basal) 
# on the original direction.
dPA <- directPA(Tc=Tc.gene, direction=c(1,1,1), annotation=Pathways.reactome)

dPA$gst[order(unlist(dPA$gst[,1])),][1:20,]
# rank substrates on the direciton of interest
sort(dPA$gene.pvalues)[1:20]

# (2) test combined effect of all 3 treatments vs controls on direction c(1,-1, 0)
# this rotates Ins by 0 degree, Wmn by 90 degree, and MK by 45 degree.
dPA <- directPA(Tc=Tc.gene, direction=c(1,-1,0), annotation=Pathways.reactome)
dPA$gst[order(unlist(dPA$gst[,1])),][1:20,]
# rank substrates on the direciton of interest
sort(dPA$gene.pvalues)[1:20]

# (3) test combined effect of all 3 perturbations vs controls on direction c(1,-1, 1)
# this rotates Ins by 0 degree, Wmn by 90 degree, and MK by 0 degree.
dPA <- directPA(Tc=Tc.gene, direction=c(1,-1,1), annotation=Pathways.reactome)
dPA$gst[order(unlist(dPA$gst[,1])),][1:20,]
# rank substrates on the direciton of interest
sort(dPA$gene.pvalues)[1:20]