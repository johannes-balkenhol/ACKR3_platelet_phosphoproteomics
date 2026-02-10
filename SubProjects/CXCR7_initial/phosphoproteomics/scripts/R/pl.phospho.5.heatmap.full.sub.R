
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
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(reshape2)
library(ggplot2)
library(PhosR)
library(dplyr)
library(matrixStats)


###############################################################################
############# select subset to present
############# plot distributions for all proteins
#setwd("C:/Users/Ozge/Nextcloud/Phosphoproteom/phosphoproteom analysis 2/scripts")

input = list(top.filter.10, top.filter.30, top.filter.60, top.filter.300, top.filter.600, top.filter.900, top.filter.1800)
#top.10.dmso.vs.0s, top.600.dmso.vs.0s, top.1800.dmso.vs.0s,
#top.10.cxcr7.vs.0s, top.600.cxcr7.vs.0s, top.1800.cxcr7.vs.0s)

names_input = c("10", "30", "60", "300", "600", "900", "1800")
#"10.dmso.vs.0s", "600.dmso.vs.0s", "1800.dmso.vs.0s",
#"10.cxcr7.vs.0s", "600.cxcr7.vs.0s", "1800.cxcr7.vs.0s")



for (i in 1:length(input)) {
  #df_pl = top.10.cxcr7.vs.0s
  df_pl2 = input[[i]]
  df_pl2 <- df_pl2[order(df_pl2$PValue),]
  assign(paste0("top.filter.ordered.", names_input[[i]]), df_pl2)
}


top.rownames <- c(rownames(top.filter.ordered.10[1:20,]),rownames(top.filter.ordered.30[1:20,]),
                  rownames(top.filter.ordered.60[1:20,]),rownames(top.filter.ordered.300[1:20,]),
                  rownames(top.filter.ordered.600[1:20,]),rownames(top.filter.ordered.900[1:20,]),
                  rownames(top.filter.ordered.1800[1:20,]))

top.norm_intensity <- norm_intensity_filter[top.rownames,]
#top.norm_intensity <- norm_intensity_filter #for heatmap all.

top.norm_intensity<- top.norm_intensity[!duplicated(top.norm_intensity), ]

rownames(top.norm_intensity) <- sapply(strsplit(rownames(top.norm_intensity), ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})


########### aggregate by mean per group 
x_tt <- as.factor(grps)

t00 <- rowMedians(top.norm_intensity[,1:7])
t10 <- rowMedians(top.norm_intensity[,8:13])
t30 <- rowMedians(top.norm_intensity[,14:20])
t60 <- rowMedians(top.norm_intensity[,21:27])
t300 <- rowMedians(top.norm_intensity[,28:34])
t600 <- rowMedians(top.norm_intensity[,35:40])
t900 <- rowMedians(top.norm_intensity[,41:47])
t1800 <- rowMedians(top.norm_intensity[,48:53])


top.norm_intensity.collapse <- cbind(t00,t10,t30,t60,t300,t600,t900,t1800)
rownames(top.norm_intensity.collapse) <- rownames(top.norm_intensity)
colnames(top.norm_intensity.collapse) <- c("t00","t10","t30","t60","t300","t600","t900", "t1800")
top.norm_intensity.collapse <- as.data.frame(top.norm_intensity.collapse)

norm_abundance2 <- top.norm_intensity.collapse
##or
#norm_abundance2 <- top.norm_intensity
### plot the distribution of the raw/filtered/imputed/sclaed/normalised intensities



#dev.new()
dd <- melt(as.matrix(norm_abundance2), variable.name = "normalised")
ggplot(dd, aes(value, colour = Var2)) + geom_density(bw = "sj")
# + xlim(0,1e+10)



###############################################################################
## use median, as more robust than average, for the Z values in heatmap
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

########### coloring according to sample names
col_groups <- substr(colnames(norm_abundance2), 1, 5)
table(col_groups)


mat_col <- data.frame(time = col_groups)
rownames(mat_col) <- colnames(norm_abundance2)

mat_colors <- list(time = brewer.pal(8, "BuPu"))
#mat_colors <- list(group = brewer.pal(length(table(col_groups)), "Set3"))

names(mat_colors$time) <- unique(col_groups)



############ create a simple heatmap
data_subset_norm <- t(apply(norm_abundance2, 1, cal_z_score))
#rownames(data_subset_norm) <- c(sapply(strsplit(rownames(data_subset_norm), ";"), "[[", c(1,2,3)), ";", sapply(strsplit(rownames(data_subset_norm), ";"), "[[", 2), ";" , sapply(strsplit(rownames(data_subset_norm), ";"), "[[", 3))
#rownames(data_subset_norm) <- sapply(strsplit(rownames(data_subset_norm), ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})


############ include breaks in the heatmap 
############ for better visualization in tailed data
############ as we use a color code cutoff regarding quantiles

mat_breaks <- seq(min(norm_abundance2, na.rm=TRUE), max(norm_abundance2, na.rm=TRUE), length.out = 20)


# define function
quantile_breaks <- function(xs, n) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm=TRUE)
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(data_subset_norm, n = 40)

tiff(filename = "../analysis/Heatmap/Heatmap_top20all.tiff",
     width = 10 * 300, 
     height = 20 * 300,
     res = 300,
     compression = "lzw")

pheatmap(mat = data_subset_norm, 
         color = colorRampPalette(c("navy", "white", "red"))(length(mat_breaks) - 1),
         #gaps_row=c(5,10,15,20,25,30,35,40,45,50, 55),
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
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         cex=1,
         clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean",
         clustering_method="complete",
         main = ""
)

dev.off()





tiff(filename = "../analysis/Heatmap/Heatmap_all_clusters.tiff",
     width = 6 * 300, 
     height = 3 * 300,
     res = 300,
     compression = "lzw")

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
                           kmeans_k = 8,
                           cluster_cols =FALSE,
                           cluster_rows = TRUE,
                           cex=1,
                           clustering_distance_rows="euclidean",
                           clustering_distance_cols="euclidean",
                           clustering_method="complete",
                           legend = T,
                           annotation_legend = F,
                           main = ""
)

dev.off()




data_subset_norm.clust <- cbind(data_subset_norm, 
                                cluster = clustered_data$kmeans$cluster)
data_subset_norm.clust <- as.data.frame(data_subset_norm.clust)
table(data_subset_norm.clust$cluster)
#  1   2   3   4   5   6   7   8 
#328 490 471 375 422 502 330 405




data_subset_norm.clust.ordered <- data_subset_norm.clust[order(data_subset_norm.clust$cluster, decreasing = F),]
row_annots <- as.data.frame(paste0("cluster", data_subset_norm.clust.ordered$cluster))
rownames(row_annots) <- rownames(data_subset_norm.clust.ordered)
colnames(row_annots) = "cluster"

annot_colors <- list(cluster = brewer.pal(8, "PRGn"))
names(annot_colors$cluster) <- unique(row_annots$cluster)


##use clusters from RUNCLUE 
input <- data_subset_norm.clust
input$ID <- rownames(norm_intensity_filter)

clustobjswithIDs <- data.frame()
for (i in 1:length(unique(clustobjs$cluster))) {
  prot_spec_anova_cluster <- as.data.frame(clustobjs[clustobjs$cluster == i,]) #CHANGE CLUSTER NAME
  cluster_prots <- unique(unlist(strsplit(prot_spec_anova_cluster$protein, split = "\\|")))
  proteins <- merged[merged$names %in% cluster_prots, "id"] #merged2 for site-centric
  clustobjsadd <- data.frame(ID  = proteins, cluster = i)
  clustobjswithIDs <- rbind(clustobjswithIDs, clustobjsadd)
}


data_subset_norm.clust2 <- merge(input,clustobjswithIDs, by = "ID",all = TRUE)
rownames(data_subset_norm.clust2) <-rownames(input)
data_subset_norm.clust2 <- data_subset_norm.clust2 %>% select(-ID)
data_subset_norm.clust2[is.na(data_subset_norm.clust2$cluster.y),"cluster.y"] <- 9
data_subset_norm.clust2.ordered <-  data_subset_norm.clust2[order(data_subset_norm.clust2$cluster.y, decreasing = F),]

data_subset_norm.clust2.ordered_noNAS <-  data_subset_norm.clust2.ordered[data_subset_norm.clust2.ordered$cluster.y != "9",]

row_annots2 <- as.data.frame(paste0("cluster", data_subset_norm.clust2.ordered$cluster.y))
rownames(row_annots2) <- rownames(data_subset_norm.clust2.ordered)
colnames(row_annots2) = "cluster"
annot_colors2 <- list(cluster = brewer.pal(9, "PRGn"))
names(annot_colors2$cluster) <- unique(row_annots2$cluster)


tiff(filename = "../analysis/Heatmap/Heatmap_ClueR_all_clusters.tiff",
     width = 10 * 300, 
     height = 10 * 300,
     res = 300,
     compression = "lzw")

pheatmap(mat = data_subset_norm.clust.ordered[1:8], 
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

dev.off()



##get clusters with sign. 

kclusters <- cbind(rownames =rownames(top.all),
                   data_subset_norm.clust[9])

for (i in 1:length(unique(data_subset_norm.clust$cluster))) {
  sign <- intersect(kclusters[kclusters$cluster == i, "rownames"], union_sig)
  assign(paste0("kcluster", i,"_sign"), sign)
  signIDs <- sapply(strsplit(sign, ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})
  assign(paste0("kcluster", i, "_signIDs"), signIDs)
  diffSign <- intersect(sign, top.names)
  assign(paste0("kcluster", i, "_diffsign"), diffSign)
  genes <- sapply(strsplit(diffSign, ";"), "[[", 1)
  assign(paste0("kcluster", i, "_genes"), genes)
  
}



##prep input for heatmap
#for specific go terms after running enrichment
#targets_norm_abundance <- dataset_df[inputforhm,1:53] #for direct subsets

#for directly the clusters
targets_norm_abundance <- dataset_df[kcluster1_sign,1:53] #for direct subsets
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



###GO enrichement for clusters

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



#go overrep plots: custom background####


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


