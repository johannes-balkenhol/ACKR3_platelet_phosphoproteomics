####################################################################################
##### Diagnosing and Correcting for batch effect
## by SPSs from phosR study
## by Genes with low variance
## by Genes with low Rank in differential regulation in this study  --> this i do
## by generating SPSs empirically
## by a combination of all

## empirical stable phosphosites
## perform deferential regulation analysis and 
## determine the peptides that are least deregulated
## define new expression set and look at quality 
## build a vector accoridng to the sample
#vect1 <- as.data.frame(table(grps))$Freq

setwd("/Users/ozgeosmanoglu/Nextcloud/Phosphoproteom/phosphoproteom analysis 2/scripts")



library("PhosR")
library(matrixStats)
library("EDASeq")
library("tidyr")
library("limma")
library("ggplot2")
library("dplyr")
library(RColorBrewer)
library(VennDiagram)
library(plotly)
library(cowplot)




##### define input and visualize

vect1 <- c(7, 6, 7, 7, 7, 6, 7, 6)
x_tt <- as.factor(
  c(rep(c("tt00", "tt10", "tt30", "tt60", "tt300","tt600", "tt900", "tt1800"), vect1)))




dataset <- 2^SummarizedExperiment::assay(ppe_imputed_scaled,"scaled")
#sum(is.na(dataset))
#dataset <- na.omit(dataset)


#### Prenormalize make a betweenLane Normalization
dataset_norm <- betweenLaneNormalization(dataset, which="upper")

colors <- brewer.pal(8, "Paired")

tiff(filename = paste0("../analysis/PCA/upperquantnorm.tiff"),
     width = 14 * 300, 
     height = 10 * 300,
     res = 300,
     compression = "lzw")
par(mfrow = c(2,1))
plotRLE(dataset, outline=FALSE, ylim=c(-1, 1), col=colors[x_tt])
plotRLE(dataset_norm, outline=FALSE, ylim=c(-1, 1), col=colors[x_tt]) #upper quantile normalization doesn't change much? 
dev.off()






n1<-plotQC(log2(dataset_norm), grps=grps, 
           labels = colnames(dataset_norm), panel = "pca") +
  ggplot2::ggtitle("normalized, back transformed")

n2<-plotQC(log2(dataset_norm), grps=grps2, 
           labels = colnames(dataset_norm), panel = "pca") +
  ggplot2::ggtitle("normalized, back transformed")

ggpubr::ggarrange(n1, n2, nrow = 2)


make_all_contrasts <- function (group, delim="_vs_", design_matrix){
  
  suppressMessages(require(limma))
  
  #/ ensure that group levels are unique
  group <- sort(unique(as.character(group)))
  
  #/ make all combinations
  cb   <- combn(group, 2, FUN = function(x){paste0(x[1], "-", x[2])})
  
  
  #/ make contrasts
  contrasts<- limma::makeContrasts(contrasts=cb, levels=design_matrix)
  colnames(contrasts) <- gsub("-", delim, colnames(contrasts))
  
  return(contrasts)
}



design_new <- model.matrix(~0+x_tt)
colnames(design_new) <- gsub("x_tt", "", colnames(design_new))
v_new <- voom(dataset_norm, design_new)
#v_new <- log2(dataset_norm)
fit_new <- lmFit(v_new, design_new)

contrast.matrix_new <- make_all_contrasts(x_tt, delim= "_vs_", design_new)

fit2 <- contrasts.fit(fit_new, contrast.matrix_new)
fit2 <- eBayes(fit2, trend=TRUE)

top.general <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, sort.by = "none")
length(rownames(top.general[top.general$adj.P.Val > 0.9,]))

top1 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=1, sort.by = "none")
top2 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=2, sort.by = "none")
top3 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=3, sort.by = "none")
top4 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=4, sort.by = "none")
top5 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=5, sort.by = "none")
top6 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=6, sort.by = "none")
top7 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=7, sort.by = "none")
top8 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=8, sort.by = "none")
top9 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=9, sort.by = "none")
top10 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=10, sort.by = "none")
top11 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=11, sort.by = "none")
top12 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=12, sort.by = "none")
top13 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=13, sort.by = "none")
top14 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=14, sort.by = "none")
top15 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=15, sort.by = "none")
top16 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=16, sort.by = "none")
top17 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=17, sort.by = "none")
top18 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=18, sort.by = "none")
top19 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=19, sort.by = "none")
top20 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=20, sort.by = "none")
top21 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=21, sort.by = "none")
top22 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=22, sort.by = "none")
top23 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=23, sort.by = "none")
top24 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=24, sort.by = "none")
top25 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=25, sort.by = "none")
top26 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=26, sort.by = "none")
top27 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=27, sort.by = "none")
top28 <- topTable(fit2, number=nrow(fit2), adjust.method="BH", p.value=1, coef=28, sort.by = "none")


p.value <- cbind(top1[,c(5)],top2[,c(5)],top3[,c(5)],top4[,c(5)],top5[,c(5)],top6[,c(5)],top7[,c(5)],top8[,c(5)],top9[,c(5)],top10[,c(5)]
                 ,top11[,c(5)],top12[,c(5)],top13[,c(5)],top14[,c(5)],top15[,c(5)],top16[,c(5)],top17[,c(5)],top18[,c(5)],top19[,c(5)],top20[,c(5)],top21[,c(5)],
                 top22[,c(5)],top23[,c(5)],top24[,c(5)],top25[,c(5)],top26[,c(5)],top27[,c(5)],top28[,c(5)])
rownames(p.value) <- rownames(fit2)

logfc <- cbind(top1[,c(1)],top2[,c(1)],top3[,c(1)],top4[,c(1)],top5[,c(1)],top6[,c(1)],top7[,c(1)],top8[,c(1)],top9[,c(1)],top10[,c(1)]
               ,top11[,c(1)],top12[,c(1)],top13[,c(1)],top14[,c(1)],top15[,c(1)],top16[,c(1)],top17[,c(1)],top18[,c(1)],top19[,c(1)],top20[,c(1)],top21[,c(1)],
               top22[,c(1)],top23[,c(1)],top24[,c(1)],top25[,c(1)],top26[,c(1)],top27[,c(1)],top28[,c(1)])
rownames(logfc) <- rownames(fit2)

basemean <- cbind(top1[,c(2)],top2[,c(2)],top3[,c(2)],top4[,c(2)],top5[,c(2)],top6[,c(2)],top7[,c(2)],top8[,c(2)],top9[,c(2)],top10[,c(2)]
                  ,top11[,c(2)],top12[,c(2)],top13[,c(2)],top14[,c(2)],top15[,c(2)],top16[,c(2)],top17[,c(2)],top18[,c(2)],top19[,c(2)],top20[,c(2)],top21[,c(2)],
                  top22[,c(2)],top23[,c(2)],top24[,c(2)],top25[,c(2)],top26[,c(2)],top27[,c(2)],top28[,c(2)])
rownames(basemean) <- rownames(fit2)


#p.value <- cbind(top7[,c(5)],top16[,c(5)],top21[,c(5)])
#rownames(p.value) <- rownames(fit2)

#logfc <- cbind(top7[,c(1)],top16[,c(1)],top21[,c(1)])
#rownames(logfc) <- rownames(fit2)

#basemean <- cbind(top7[,c(2)],top16[,c(2)],top21[,c(2)])
#rownames(basemean) <- rownames(fit2)



p_all <- as.data.frame(rowMedians(p.value))
colnames(p_all) <- c("P")
rownames(p_all) <- rownames(p.value)

logfc_all <- as.data.frame(rowMedians(logfc))
colnames(logfc_all) <- c("logfc")
rownames(logfc_all) <- rownames(logfc)

basemean_all <- as.data.frame(rowMedians(basemean))
colnames(basemean_all) <- c("basemean")
rownames(basemean_all) <- rownames(basemean)

p.value.count <- as.data.frame(sapply(1:nrow(p.value), function(i) sum(p.value[i,] < 0.05)))
colnames(p.value.count) <- c("P")
rownames(p.value.count) <- rownames(p.value)


#quantile(basemean)
empirical_basemean <- rownames(subset(basemean_all, basemean >= as.numeric(quantile(basemean, probs = seq(0, 1, 0.05))['25%'])))
length(empirical_basemean)




#final
empirical_p <- rownames(subset(p_all, P >= 0.95))
length(empirical_p)

empirical_p_count <- rownames(subset(p.value.count, P < 1))
length(empirical_p_count)

p_rank <- as.data.frame(colMedians(colRanks(p.value)))
rownames(p_rank) <- rownames(p.value)
colnames(p_rank) <- "V1"
empirical_p_rank <- rownames(subset(p_rank, V1 > 3000)) 
length(empirical_p_rank)

empirical_logfc <- rownames(subset(logfc_all, abs(logfc) <= 0.05))
length(empirical_logfc)



#Reduce(intersect, list(a,b,c))
empirical_topall<- Reduce(intersect, list(empirical_p,empirical_p_count,empirical_p_rank,empirical_logfc))

#empirical_topall <- empirical_p
length(empirical_topall)

write.table(empirical_topall, "data/processed_data/empirical_topall.txt", sep = "\t")



#INTERSECT WITH NEW DATA's empirical topall

## use script pl.phospho.3.2.compare.old.new
empirical_match<- read.table("../data/processed_data/empirical_match.txt", sep="\t", header=TRUE, dec=",")
#length(Reduce(intersect, list(empirical_topall, empirical_match)))


#b<- read.delim("reanalysis_141122/intersect_empirical_topall_oldreanalyzed_andnew.txt")
#final_empirical_intersect<- b[,1]

## write tables locally or
## RS: connect to mysql server and wrtie the table from the mysql database
#write.table(empirical_topall, "../data/processed_data/empirical_topall.txt", sep = "\t")


#use empirical_control_proteins3.R

empirical_top <- Reduce(intersect, list(empirical_p_rank,empirical_basemean))
#empirical_topall <- empirical_p
length(empirical_top)

## write tables locally or
## RS: connect to mysql server and wrtie the table from the mysql database
#write.table(empirical_top, "../data/processed_data/empirical_top.txt", sep = "\t")



##### RUVphospho (RUVIII) normalization
design = model.matrix(~ grps - 1)
ctl2 = which(rownames(ppe_imputed_scaled) %in% empirical_topall)
#ctl2 = which(rownames(ppe_imputed_scaled) %in% empirical_match$x)



#### RUVphospho normalization
## RUVphospho shows a good tendency for statistics, but the intensity value are curiously shifted 
## but still 
#ppe = RUVphospho(ppe_imputed_scaled, M = design, k = 16, ctl = ctl2)

ppe_norm = RUVphospho(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), M = design, k = 16, ctl = ctl2)

ppe <- PhosphoExperiment(assays = list(normalised = as.matrix(ppe_norm)), 
                         UniprotID = sapply(strsplit(rownames(ppe_imputed_scaled), ";"), "[[", 1),
                         Site = as.numeric(gsub("[A-Z]","", sapply(strsplit(rownames(ppe_imputed_scaled), ";"), "[[", 3))),
                         GeneSymbol = sapply(strsplit(rownames(ppe_imputed_scaled), ";"), "[[", 2),
                         Residue = gsub("[0-9]","", sapply(strsplit(rownames(ppe_imputed_scaled), ";"), "[[", 3)),
                         Sequence = sapply(strsplit(rownames(ppe_imputed_scaled), ";"), "[[", 4))

norm_intensity <- SummarizedExperiment::assay(ppe, "normalised")



write.table(norm_intensity, "../data/processed_data/norm_intensity.txt", sep = "\t")


############################
### Visualization  for the PCA graphic of in the Power Point

### PCA by group of selected psites
## plot clusteirng PCA
test2 <- SummarizedExperiment::assay(ppe_imputed_scaled,"scaled")
sum(is.na(test2))
test4 <- SummarizedExperiment::assay(ppe,"normalised")
test5 <- norm_intensity_filter




tiff(filename = paste0("../analysis/PCA/PCA_after_normalization.tiff"),
     width = 5 * 300, 
     height = 4 * 300,
     res = 300,
     compression = "lzw")

plotQC(test4, grps=grps, 
       labels = "", panel="pca") + theme_cowplot()

dev.off()



tiff(filename = paste0("../analysis/PCA/dendrogram_after_normalization.tiff"),
     width = 8 * 300, 
     height = 3 * 300,
     res = 300,
     compression = "lzw")

plotQC(test4, grps=grps, 
       labels = sapply(strsplit(colnames(ppe), "_"), "[[",1), panel="dendrogram")


  
dev.off()


tiff(filename = paste0("../analysis/PCA/dendrogram_imputedscaled.tiff"),
     width = 8 * 300, 
     height = 3 * 300,
     res = 300,
     compression = "lzw")
plotQC(SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"), 
       labels=colnames(ppe_imputed_scaled), 
       panel = "dendrogram", grps = grps)
dev.off()




## input
#input <- top.filter.600[top.filter.600[, "adj.P.Val"] <=1,]   ## generated in 3.0
#input <- cbind(top.600.sign)
input <- test2
filter <- test2

## PCA slected proteins
#dev.new()

plotQC(input[row.names(filter),], grps = grps, 
       panel="pca", labels = sapply(strsplit(colnames(filter), "_"), "[[",1))



### reverse PCA by psites 
grps_psite = row.names(input)
up <- sapply(strsplit(grps_psite, ";"), "[[",1)
name <- sapply(strsplit(grps_psite, ";"), "[[",2)
sitex <- sapply(strsplit(grps_psite, ";"), "[[",3)
grps_psite <- paste(up, name, sitex, sep=";")

test5 <- t(input[row.names(input),])

plotQC(test5, grps = "", panel="pca", labels = name)


### PCA tidyverse (figure for publication)

# Create a matrix from our table of counts
#test7 <- test4[row.names(input),]

pca_matrix <- input %>% 
  # coerce to a matrix
  as.matrix() %>% 
  # transpose the matrix so that rows = samples and columns = variables
  t()

# Perform the PCA
sample_pca <- prcomp(pca_matrix)

pca_matrix[1:10, 1:5]
as_tibble(pca_matrix)
as_tibble(pca_matrix, rownames = "sample")


pc_eigenvalues <- sample_pca$sdev^2

# and a variable with the variances
pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues

pc_eigenvalues %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")


pc_scores <- sample_pca$x


pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

pc_scores %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

variance_percent <- sample_pca$sdev^2 / sum(sample_pca$sdev^2) * 100

shape_mapping <- c("0000" = 18 , "0010" = 20, "0030" = 20, "0060" = 20,
                   "0300" = 20, "0600" = 20, "0900" = 20, "1800" = 20)

color_palette_own <- c("0000" = "#FF0000", "0010" = "#FFBF00", "0030"="#80FF00", 
                       "0060"= "#00FF40", "0300"="#00FFFF","0600"=  "#0040FF",
                       "0900"= "#8000FF", "1800"="#FF00BF")

color_palette_own_donor <- c("1" = "#1B9E77", "2" =  "#D95F02", "3"= "#7570B3", 
                             "4"= "#E7298A", "5"= "#66A61E" ,"6"=  "#E6AB02" ,
                             "7"=  "#A6761D", "8"="#666666")


sample_info_IN <- sample_info 
sample_info_IN$time <- as.factor(sample_info_IN$time)
sample_info_IN$donor <- as.factor(sample_info_IN$donor)
sample_info_IN[sample_info_IN$condition == "ctrl","condition"] <- "DMSO"

tiff(filename = paste0("../analysis/PCA/pca_jbstylenew_scaled3.tiff"),
     width = 6 * 300, 
     height = 4* 300,
     res = 300,
     compression = "lzw")
pca_plot <- sample_pca$x %>% # extract the loadings from prcomp
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample") %>% 
  # join with "sample_info" table 
  full_join(sample_info_IN, by = "sample") %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = time, shape=condition)) +
  geom_point(size=3) +
  theme_cowplot() +
  scale_shape_manual(values = c(17,19)) +
  #scale_fill_manual(values = shape_mapping) +
  scale_color_manual(values = color_palette_own)+
  labs(x = paste("PC1 (", round(variance_percent[1], 2), "%)"),
       y = paste("PC2 (", round(variance_percent[2], 2), "%)")) +
  coord_fixed()

pca_plot
dev.off()


##3d pcas

library(plotly)

pcs <- sample_pca$x %>% # extract the loadings from prcomp
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample") %>% 
  # join with "sample_info" table 
  full_join(sample_info_IN, by = "sample") 
  # create the plot

  
plot_ly(x=pcs$PC1, y=pcs$PC2, z=pcs$PC3, type="scatter3d", mode="markers", color=substr(pcs$sample, 2, 5))

  
  

 

##--------------------------------------

   
tiff(filename = paste0("../analysis/PCA/pca_jbstylenew_scaled_donors3.tiff"),
     width = 6 * 300,
     height = 4* 300,
     res = 300,
     compression = "lzw")
pca_plot <- sample_pca$x %>% # extract the loadings from prcomp
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample") %>% 
  # join with "sample_info" table 
  full_join(sample_info_IN, by = "sample") %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, color = donor, shape=condition)) +
  geom_point(size=3) +
  theme_cowplot() +
  scale_shape_manual(values = c(17,19)) +
  #scale_fill_manual(values = shape_mapping) +
  scale_color_manual(values = color_palette_own_donor)+
  labs(x = paste("PC1 (", round(variance_percent[1], 2), "%)"),
       y = paste("PC2 (", round(variance_percent[2], 2), "%)"))+
  coord_fixed()

pca_plot
dev.off()



pc_loadings <- sample_pca$rotation

pc_loadings <- pc_loadings %>% 
  as_tibble(rownames = "gene")

# print the result
pc_loadings


top_genes <- pc_loadings %>% 
  # select only the PCs we are interested in
  select(gene, PC1, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(variable, desc(abs(value))) %>% 
  # take the 10 top rows
  slice(1:10) %>% 
  # pull the gene column as a vector
  pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

topn_loadings <- pc_loadings %>%
  # Select only the PCs we are interested in
  select(gene, PC1, PC2, PC3, PC4) %>%
  # Convert to long format
  pivot_longer(cols = starts_with("PC"), names_to = "variable", values_to = "value") %>%
  # For each PC, sort by the absolute value of the loadings in descending order
  group_by(variable) %>%
  arrange(variable, desc(abs(value))) %>%
  # Get the top 10 entries for each PC
  slice_head(n = 5) %>%
  ungroup()

head(topn_loadings)



topn_loadings$gene <- paste0(sapply(strsplit(topn_loadings$gene, ";"), "[[",2), ";",
                             sapply(strsplit(topn_loadings$gene, ";"), "[[",3))
ggplot(topn_loadings,
      aes(x=gene,y=value)) +  
  geom_bar(stat='identity') + 
  theme_cowplot() +
  theme(axis.text.x = element_text(angle=90, hjust = 1)) +
  facet_grid(~variable, scales = "free_x")


pc_colors <- c("PC1" = "#804862","PC2" = "#B996A5", "PC3" = "#D5BFC8", "PC4" = "#F1E9EC")

ggplot(topn_loadings, 
       aes(x = reorder(gene, abs(value)), y = value)) +  
  geom_bar(stat = 'identity', color = "black", fill ="#D5BFC8") + 
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position = "top") + 
  facet_wrap(~variable, scales = "free_x") +
  #scale_fill_manual(values = pc_colors) +  # Customize colors for 4 PCs
  labs(x = "", y = "value", fill = "Principal Component")

ggsave("../analysis/PCA/top5_Loadings.tiff",
       device = "tiff",
       dpi = 300,
       width = 6,
       height = 10)

top_loadings <- pc_loadings %>% 
  filter(gene %in% top_genes)

loadings_plot <- ggplot(data = top_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
loadings_plot

tiff(filename = paste0("../analysis/PCA/loadings_jbstyle.tiff"),
     width = 14 * 300, 
     height = 4 * 300,
     res = 300,
     compression = "lzw")
ggplot(data = top_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))
dev.off()
