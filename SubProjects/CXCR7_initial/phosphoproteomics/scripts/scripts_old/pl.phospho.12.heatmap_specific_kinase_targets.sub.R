
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(reshape2)
library(ggplot2)
library(PhosR)
library(dplyr)
library(matrixStats)



pka_targets_psplus <- PhosphoSite.human[["PRKACA"]]
pkg_targets_psplus <- PhosphoSite.human[["PRKG1"]]





dataset <- SummarizedExperiment::assay(ppe, "normalised")
dataset_df <- as.data.frame(dataset)
dataset_df = dataset_df[!(row.names(dataset_df) %in% not.uniq[,1]), ]


name <- sapply(strsplit(rownames(dataset_df), ";"), "[[", 2)
site <- sapply(strsplit(rownames(dataset_df), ";"), "[[", 3)


dataset_df$namesite <- paste0(name,";",site,";")

kinase <- PhosphoSite.human[["ROCK1"]]

targets_norm_abundance <- dataset_df[dataset_df$namesite %in% kinase,1:53]
targets_norm_abundance <- as.matrix(targets_norm_abundance)


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


dd <- melt(as.matrix(targets_norm_abundance2), variable.name = "normalised")
ggplot(dd, aes(value, colour = Var2)) + geom_density(bw = "SJ")



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

tiff(filename = "analysis/Heatmap/Heatmap_rock1.tiff",
     width = 6 * 300, 
     height = 4 * 300,
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



#######

