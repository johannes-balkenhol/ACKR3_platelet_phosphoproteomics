
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

#names_input <- list("10", "600", "1800") 
                    
names_input <-  list("10.cxcr7.vs.0s", "600.cxcr7.vs.0s", "1800.cxcr7.vs.0s",
                    "10.dmso.vs.0s", "600.dmso.vs.0s", "1800.dmso.vs.0s")

#input <- list(top.10, top.600, top.1800)

input <-  list(top.10.cxcr7.vs.0s, top.600.cxcr7.vs.0s, top.1800.cxcr7.vs.0s,
                  top.10.dmso.vs.0s, top.600.dmso.vs.0s, top.1800.dmso.vs.0s)


for (i in 1:length(input)) {
  #df_pl = top.10.cxcr7.vs.0s
  df_pl2 = input[[i]]
  df_pl2 <- df_pl2[order(df_pl2$adj.P.Val),]
  rownames(df_pl2) <- paste0(rownames(df_pl2), ";", df_pl2$symbol)
  assign(paste0("top.ordered.", names_input[[i]]), df_pl2)
}


top.rownames <- c(rownames(top.ordered.10[1:20,]),
                  rownames(top.ordered.600[1:20,]),
                  rownames(top.ordered.1800[1:20,]))

# top.rownames <- c(rownames(top.ordered.10.cxcr7.vs.0s[1:20,]),
#                   rownames(top.ordered.600.cxcr7.vs.0s[1:20,]),
#                   rownames(top.ordered.1800.cxcr7.vs.0s[1:20,]))
# 
# top.rownames <- c(rownames(top.ordered.10.dmso.vs.0s[1:20,]),
#                   rownames(top.ordered.600.dmso.vs.0s[1:20,]),
#                   rownames(top.ordered.1800.dmso.vs.0s[1:20,]))

df_norm_intensity <-as.data.frame(norm_intensity)
df_norm_intensity$symbol <- mapIds(org.Hs.eg.db, keys=rownames(df_norm_intensity), 
                                   column="SYMBOL", keytype="UNIPROT", multiVals="first")
rownames(df_norm_intensity) <- paste0(rownames(df_norm_intensity), ";", df_norm_intensity$symbol)
df_norm_intensity <- df_norm_intensity %>% dplyr::select(-symbol)

top.norm_intensity <- df_norm_intensity[top.rownames,]
top.norm_intensity<- top.norm_intensity[!duplicated(top.norm_intensity), ]
rownames(top.norm_intensity) <- sapply(strsplit(rownames(top.norm_intensity), ";"),  function(x){paste(x[[1]], x[[2]], sep="_")})


########### aggregate by mean per group 
x_tt <- as.factor(grps)
top.norm_intensity <-  as.matrix(top.norm_intensity)

t00 <- rowMedians(top.norm_intensity[,1:10])
t10 <- rowMedians(top.norm_intensity[,11:19])
t10wt <- rowMedians(top.norm_intensity[,20:29])
t600 <- rowMedians(top.norm_intensity[,30:39])
t600wt <- rowMedians(top.norm_intensity[,40:49])
t1800 <- rowMedians(top.norm_intensity[,50:59])
t1800wt <- rowMedians(top.norm_intensity[,60:69])

top.norm_intensity.collapse <- cbind(t00,t10,t10wt,t600,t600wt,t1800,t1800wt)
rownames(top.norm_intensity.collapse) <- rownames(top.norm_intensity)
colnames(top.norm_intensity.collapse) <- c("t00_ctrl","t10_cxcr7","t10_dmso","t600_cxcr7","t600_dmso","t1800_cxcr7","t1800_dmso")
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

mat_breaks <- seq(min(norm_abundance2, na.rm=TRUE), max(norm_abundance2, na.rm=TRUE), length.out = 20)


# define function
quantile_breaks <- function(xs, n) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm=TRUE)
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(data_subset_norm, n = 40)

tiff(filename = "../analysis/Heatmap/Heatmap_top20_CXCR7vsDMSO.tiff",
     width = 10 * 300, 
     height = 14 * 300,
     res = 300,
     compression = "lzw")

pheatmap(mat = data_subset_norm, 
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
  