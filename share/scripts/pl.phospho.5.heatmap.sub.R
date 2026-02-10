
################################################################
################################################################
### show heatmap


## install packages
install.packages("tidyverse")
install.packages("pheatmap")
install.packages("viridis")
install.packages("reshape2")
install.packages("matrixStats")


## load packages
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(matrixStats)
library(org.Hs.eg.db)
library(pheatmap)
library(PhosR)
library(RColorBrewer)
library(reshape2)
library(viridis)






###############################################################################
############# select subset to present
############# plot distributions for all proteins
#setwd("F:/Masterarbeit_BioWis/Phosphoproteomics/Phosphoproteom validation/scripts/")
input = list(top.filter.10, top.filter.600, top.filter.1800, 
             top.filter.10.dmso.vs.0s, top.filter.600.dmso.vs.0s, top.filter.1800.dmso.vs.0s,
             top.filter.10.cxcr7.vs.0s, top.filter.600.cxcr7.vs.0s, top.filter.1800.cxcr7.vs.0s)

names_input = c("10", "600", "1800", 
                "10.dmso.vs.0s", "600.dmso.vs.0s", "1800.dmso.vs.0s",
                "10.cxcr7.vs.0s", "600.cxcr7.vs.0s", "1800.cxcr7.vs.0s")

for (i in 1:length(input)) {
  #df_pl = top.10.cxcr7.vs.0s
  df_pl2 = input[[i]]
  df_pl2 <- df_pl2[order(df_pl2$PValue),]
  assign(paste0("top.filter.ordered.", names_input[[i]]), df_pl2)
}


top.rownames <- c(rownames(top.filter.ordered.10[1:20,]),rownames(top.filter.ordered.600[1:20,]),
                  rownames(top.filter.ordered.1800[1:20,]))

#top.norm_intensity <- norm_intensity_filter[top.rownames,]
top.norm_intensity <- norm_intensity_filter #for heatmapall

rownames(top.norm_intensity) <- sapply(strsplit(rownames(top.norm_intensity), ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})


########### aggregate by mean per group 
x_tt <- as.factor(grps)

t00 <- rowMedians(top.norm_intensity[,1:10])
t10 <- rowMedians(top.norm_intensity[,11:20])
t10wt <- rowMedians(top.norm_intensity[,21:30])
t600 <- rowMedians(top.norm_intensity[,31:40])
t600wt <- rowMedians(top.norm_intensity[,41:50])
t1800 <- rowMedians(top.norm_intensity[,51:60])
t1800wt <- rowMedians(top.norm_intensity[,61:70])


top.norm_intensity.collapse <- cbind(t00,t10,t10wt,t600,t600wt,t1800,t1800wt)
rownames(top.norm_intensity.collapse) <- rownames(top.norm_intensity)
colnames(top.norm_intensity.collapse) <- c("t00_ctrl","t10_cxcr7","t10_dmso","t600_cxcr7","t600_dmso","t1800_cxcr7","t1800_dmso")
top.norm_intensity.collapse <- as.data.frame(top.norm_intensity.collapse)

norm_abundance2 <- top.norm_intensity.collapse
##or
#norm_abundance2 <- top.norm_intensity
### plot the distribution of the raw/filtered/imputed/sclaed/normalised intensities




dd <- melt(as.matrix(norm_abundance2), variable.name = "normalised")
ggplot(dd, aes(value, colour = Var2)) + geom_density(bw = "sj")
# + xlim(0,1e+10)



###############################################################################
## use median, as more robust than average, for the Z valus in heatmap
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

########### coloring accoridng to sample names
col_groups <- substr(colnames(norm_abundance2), 1, 10)
table(col_groups)
col_groups <- colnames(norm_abundance2)


mat_col <- data.frame(time = col_groups)
rownames(mat_col) <- colnames(norm_abundance2)

mat_colors <- list(time = brewer.pal(7, "BuPu"))
#mat_colors <- list(group = brewer.pal(length(table(col_groups)), "Set3"))

names(mat_colors$time) <- unique(col_groups)



############ create a simple heatmap
data_subset_norm <- t(apply(norm_abundance2, 1, cal_z_score))
#rownames(data_subset_norm) <- c(sapply(strsplit(rownames(data_subset_norm), ";"), "[[", c(1,2,3)), ";", sapply(strsplit(rownames(data_subset_norm), ";"), "[[", 2), ";" , sapply(strsplit(rownames(data_subset_norm), ";"), "[[", 3))
#rownames(data_subset_norm) <- sapply(strsplit(rownames(data_subset_norm), ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})


############ include breaks in the heatmap 
############ for better visualization in tailed data
############ as we use a color code cutoff regarding quantiles

mat_breaks <- seq(min(norm_abundance2, na.rm=TRUE), max(norm_intensity, na.rm=TRUE), length.out = 20)

# define function
quantile_breaks <- function(xs, n) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm=TRUE)
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(data_subset_norm, n = 40)

tiff(filename = "../analysis/Heatmap/Heatmap_top20all.tiff",
     width = 8 * 300, 
     height = 12 * 300,
     res = 300,
     compression = "lzw")
pheatmap(mat = data_subset_norm, 
         color = colorRampPalette(c("navy", "white", "red"))(length(mat_breaks) - 1),
         gaps_row=c(5,10,15,20,25,30,35,40,45,50,55),
         scale="row",
         na_col = "grey",
         breaks = mat_breaks,
         border_color = "white", 
         show_colnames = TRUE, 
         show_rownames = TRUE, 
         annotation_col = mat_col, 
         annotation_colors = mat_colors, 
         drop_levels = TRUE, 
         fontsize = 8, 
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         cex=1,
         clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean",
         clustering_method="complete",
         main = ""
)

dev.off()


# Hierarchical Clustering ----
set.seed(6)
clustered_data <- pheatmap(mat = data_subset_norm, 
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
                           fontsize = 12, 
                           kmeans_k = 5,
                           cluster_cols =FALSE,
                           cluster_rows = TRUE,
                           cex=1,
                           clustering_distance_rows="correlation",
                           #clustering_distance_cols="euclidean",
                           clustering_method="complete",
                           legend = T,
                           annotation_legend = F,
                           main = ""
)


data_subset_norm.clust <- cbind(data_subset_norm, 
                                cluster = clustered_data$kmeans$cluster)
data_subset_norm.clust <- as.data.frame(data_subset_norm.clust)
table(data_subset_norm.clust$cluster)


data_subset_norm.clust.ordered <- data_subset_norm.clust[order(data_subset_norm.clust$cluster, decreasing = F),]
row_annots <- as.data.frame(paste0("cluster", data_subset_norm.clust.ordered$cluster))
rownames(row_annots) <- rownames(data_subset_norm.clust.ordered)
colnames(row_annots) = "cluster"
row_annots$cluster <- factor(row_annots$cluster, levels = c("cluster1","cluster2","cluster3","cluster4","cluster5",
                                "cluster6","cluster7"))
                                #,"cluster8","cluster9", "cluster10"))

annot_colors <- list(cluster = brewer.pal(7, "PRGn"))
names(annot_colors$cluster) <- unique(row_annots$cluster)

pheatmap(mat = data_subset_norm.clust.ordered[1:7], 
         color = colorRampPalette(c("navy", "white", "red"))(length(mat_breaks) - 1),
         gaps_row=cumsum(as.numeric(table(row_annots$cluster))),
         scale="row",
         #na_col = "grey",
         breaks = mat_breaks,
         #border_color = "white", 
         show_colnames = TRUE, 
         show_rownames = F, 
         annotation_col = mat_col, 
         annotation_colors = c(mat_colors,annot_colors),
         annotation_row = row_annots,
         #drop_levels = TRUE, 
         fontsize = 12, 
         cluster_cols = F,
         cluster_rows = F,
         cex=1,
         #clustering_distance_rows="euclidean",
         #clustering_distance_cols="euclidean",
         #clustering_method="complete",
         main = "")


##get clusters with significant changes ----

kclusters <- cbind(rownames =rownames(top.all),
                   data_subset_norm.clust[8])

for (i in 1:length(unique(data_subset_norm.clust$cluster))) {
  sign <- intersect(kclusters[kclusters$cluster == i, "rownames"], union_sig)
  assign(paste0("kcluster", i,"_sign"), sign)
  signIDs <- sapply(strsplit(sign, ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})
  assign(paste0("kcluster", i, "_signIDs"), signIDs)
  diffSign <- intersect(sign, top.names)
  assign(paste0("kcluster", i, "_diffsign"), diffSign)
  genes <- sapply(strsplit(sign, ";"), "[[", 1)
  assign(paste0("kcluster", i, "_genes"), genes)
  
}


##GO enrichement for clusters ---

library(clusterProfiler)
library(org.Hs.eg.db)


i = 1

go_enrich_pl <- enrichGO(gene = get(paste0("kcluster", i, "_genes")),
                         universe = top.collapse$uniprot_id,
                         OrgDb = org.Hs.eg.db, 
                         keyType = 'UNIPROT',
                         readable = T,
                         ont = "BP",
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 1,
                         pAdjustMethod = "none")

#assign(paste0("go_enrich_pl", names_input[[i]]), go_enrich_pl)
#go overrep plots: save result to file####

pl.tab = go_enrich_pl@result


write.table(pl.tab, file = paste0("../analysis/Heatmap/Clusters/cluster", i, "_GOBP_CUSTOM_all.txt"), sep = "\t", quote = F, 
            row.names = F, col.names = T)



#go overrep plots: custom background


tiff(filename = paste0("../analysis/Heatmap/Clusters/cluster", i, "dotplot.tiff"),
     width = 8 * 300, 
     height = 6 * 300,
     res = 300,
     compression = "lzw")

dotplot(go_enrich_pl, title = paste("Cluster", i, sep=" "), font.size=14,
        showCategory = 10)

dev.off()


tiff(filename = paste0("../analysis/Heatmap/Clusters/cluster", i, "barplot.tiff"),
     width = 8 * 300, 
     height = 6 * 300,
     res = 300,
     compression = "lzw")

barplot(go_enrich_pl, 
        drop = TRUE, 
        x= "GeneRatio",
        showCategory = 10, 
        title = paste("Cluster", i, sep=" "),
        font.size = 14,
)

dev.off()


inputforhm <- kcluster1_diffsign[str_detect(kcluster1_diffsign,
                                            gsub("/", "|", pl.tab[pl.tab$Description == 
                                                                    "second-messenger-mediated signaling","geneID"]))]


## Make heatmaps ----
#for specific go terms after running enrichment
targets_norm_abundance <- dataset_df[inputforhm,1:53] #for direct subsets

#for directly the clusters
#targets_norm_abundance <- dataset_df[kcluster1_sign,1:53] #for direct subsets
#continue from here
targets_norm_abundance <-targets_norm_abundance[order(row.names(targets_norm_abundance)), ]
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

targets_logfcs<-round(targets_logfcs, digits = 2)

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



###############################################################################
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
targets_data_subset_norm <- t(apply(targets_norm_abundance2, 1, cal_z_score))
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

mainOI <- "cluster1 - second-messenger-mediated signaling (GO:0019932)"

tiff(filename = "../analysis/Heatmap/Clusters/Heatmap_kcluster1SC_MSN.tiff",
     width = 10 * 300, 
     height = 12 * 300,
     res = 300,
     compression = "lzw")

pheatmap(mat = targets_data_subset_norm, 
         color = colorRampPalette(c("navy", "white", "red"))(length(mat_breaks) - 1),
         scale="row",
         na_col = "grey",
         breaks = mat_breaks,
         border_color = "white", 
         show_colnames = TRUE, 
         show_rownames = T, 
         annotation_col = mat_col, 
         annotation_colors = mat_colors,
         #annotation_row = row_annots,
         drop_levels = TRUE, 
         fontsize = 8, 
         cluster_cols = F,
         cluster_rows = F,
         cex=1,
         display_numbers = test_labels,
         number_color = "green", 
         fontsize_number = 8,
         cellheight=10, cellwidth = 30,
         main = mainOI)
dev.off()
