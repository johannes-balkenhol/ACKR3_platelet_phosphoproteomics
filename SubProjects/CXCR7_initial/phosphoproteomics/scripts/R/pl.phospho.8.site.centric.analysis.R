
suppressPackageStartupMessages({
  library(parallel)
  library(ggplot2)
  library(ClueR)
  library(reactome.db)
  library(org.Hs.eg.db)
  library(annotate)
  library(PhosR)
  library(dplyr)
  library(tidyr)
  library(plyr)
  library(cowplot)
  library("ggsci")
})

data(PhosphoSitePlus)

#setwd("/Users/ozgeosmanoglu/Nextcloud/Phosphoproteom/phosphoproteom analysis 2/scripts")
#setwd("/Users/ozgeosmanoglu/Nextcloud/CXCR7_platelet_analysis/phosphoproteom analysis 2/scripts")

background_list <-unique(sapply(strsplit(rownames(norm_intensity_collapse), ";"), "[[", 2))
background_list_sites <-unique(paste0(sapply(strsplit(rownames(norm_intensity_filter), ";"), "[[", 2), 
                                      ";", sapply(strsplit(rownames(norm_intensity_filter), ";"), "[[", 3), ";"))

#OUR WAY - DEGS####

#~~~~~~~~~~~~~~~~~~~~
## Gene-Centric#######
#~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### prepare pathways ----------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pathways = as.list(reactomePATHID2EXTID)
path_names = as.list(reactomePATHID2NAME)
name_id = match(names(pathways), names(path_names))
names(pathways) = unlist(path_names)[name_id]

pathways = pathways[which(grepl("Homo sapiens", names(pathways), ignore.case = TRUE))]

pathways = lapply(pathways, function(path) {
  gene_name = unname(getSYMBOL(path, data = "org.Hs.eg"))
  toupper(unique(gene_name))
})# perform ClueR to identify optimal number of clusters
names(pathways) <- gsub("_", " ", gsub("REACTOME_|Homo sapiens: ", "", names(pathways)))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### our method all differential sites ---------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
input = list(top.filter.10, top.filter.30, top.filter.60, top.filter.300, top.filter.600,
             top.filter.900, top.filter.1800)
names_input = c("10", "30", "60", "300", "600", "900", "1800")

### select top differential sites
input2 <- input[c(1:7)]

for (i in 1:length(input2)) {
  #top.names <- append(top.names, rownames(input2[[i]][input2[[i]]$PValue<0.05,]))
  #assign(paste0(".", names_input[[i]]), )
  if (i == 1 & i <= 2) { 
    top.names <- rownames(input2[[i]][input2[[i]]$PValue<0.05 & abs(input2[[i]]$logFC)>0,])
  } else {
    #print(as.list(test[[i]]))
    top.names <- append(top.names, rownames(input2[[i]][input2[[i]]$PValue<0.05 & abs(input2[[i]]$logFC)>0,]))
  }
}

top.names <- top.names[!duplicated(top.names)] 


norm_intensity2 <- as.data.frame(norm_intensity_filter)
top.norm_intensity <- norm_intensity2[top.names,]
top.norm_intensity$ID <-  rownames(top.norm_intensity)
top.norm_intensity<- top.norm_intensity[!duplicated(top.norm_intensity), ]
top.norm_intensity <- as.matrix(top.norm_intensity[1:53])

########### aggregate by mean per group 

t00 <- rowMedians(top.norm_intensity[,1:7])
t10 <- rowMedians(top.norm_intensity[,8:13])
t30 <- rowMedians(top.norm_intensity[,14:20])
t60 <- rowMedians(top.norm_intensity[,21:27])
t300 <- rowMedians(top.norm_intensity[,28:34])
t600 <- rowMedians(top.norm_intensity[,35:40])
t900 <- rowMedians(top.norm_intensity[,41:47])
t1800 <- rowMedians(top.norm_intensity[,48:53])

top.norm_intensity.median <- cbind(t00,t10,t30,t60,t300,t600,t900,t1800)
rownames(top.norm_intensity.median) <- rownames(top.norm_intensity)
colnames(top.norm_intensity.median) <- c("X0000","X0010","X0030","X0060","X0300","X0600","X0900","X1800")

ppm_gene_degs <- phosCollapse(top.norm_intensity.median, 
                              sapply(strsplit(rownames(top.norm_intensity.median), ";"), "[[", 2), 
                              stat = apply(abs(top.norm_intensity.median), 1, max), by = "max")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### perform ClueR to identify optimal number of clusters ------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ppm_gene <- ppm_gene_degs

#RNGkind("L'Ecuyer-CMRG")

set.seed(123) 
c1 <- runClue(ppm_gene, annotation=pathways, 
              kRange = seq(2,10), rep = 5, effectiveSize = c(5, 100), 
              pvalueCutoff = 0.05, alpha = 0.5,
              universe = background_list)
c1[["Tc"]] <- as.matrix(c1[["Tc"]]) #without this, sometimes you get empty plots in the best graph

# Visualise the evaluation results
data <- data.frame(Success=as.numeric(c1$evlMat), Freq=rep(seq(2,10), each=5))
myplot <- ggplot(data, aes(x=Freq, y=Success)) + 
  geom_boxplot(aes(x = factor(Freq), fill="gray")) +
  stat_smooth(method="loess", colour="red", linewidth=3, span = 0.5) +
  xlab("# of cluster") + 
  ylab("Enrichment score") + 
  #scale_x_discrete(labels= c(10:20))+
  theme_classic()
myplot

supp = "deg"
tiff(filename = paste0("../analysis/Clusters/gene_centric/DEG/gene_centric_boxplot_",supp, ".tiff"),
     width = 5* 300, 
     height = 3* 300,
     res = 300,
     compression = "lzw")
myplot
dev.off()

tiff(filename = paste0("../analysis/Clusters/gene_centric/DEG/gene_centric_best_",supp, ".tiff"),
     width = 15 * 300, 
     height = 6* 300,
     res = 300,
     compression = "lzw")
#dev.new()
set.seed(123)
best <- clustOptimal(c1, rep=5, mfrow=c(2, 5), visualize = T,
                     universe = background_list)
dev.off()

#scaled x-axis
tiff(filename = paste0("../analysis/Clusters/gene_centric/DEG/gene_centric_best_",supp, "_scaled.tiff"),
     width = 15 * 300, 
     height = 6* 300,
     res = 300,
     compression = "lzw")
set.seed(123)
best <- clustOptimalÖO(c1, rep=5, mfrow=c(2, 5), visualize = "scaled", 
                       #user.maxK = 10, #if you want to decide number of clusters yourself.
                       universe = background_list)
dev.off()

clusters_df <-  data.frame()
for (i in 1:length(best$enrichList)) {
  clusters_df <- rbind(clusters_df, 
                       data.frame(best$enrichList[[i]], cluster = names(best[["enrichList"]])[i]))
}

clusters_df$size <-  as.numeric(clusters_df$size)
clusters_df$pvalue <-  as.numeric(clusters_df$pvalue)

theme_set(theme_cowplot())
colors <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF",
            "#8491B4FF", "#91D1C2FF","#DC0000FF", "#7E6148FF", "#B09C85FF")

for (i in 1:length(best$enrichList)) {
  ggdata <-  clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],]
  ggdata <- ggdata[order(ggdata$pvalue,decreasing = TRUE),]
  ordered_kinases <- factor(ggdata$kinase,levels = ggdata$kinase)
  
  g <- ggplot(ggdata, 
              aes(x = ordered_kinases, y = size))+
    geom_col(aes(fill = pvalue),width = 0.75)+
    geom_text(aes(label = kinase, y = 0.2), hjust = 0) +
    scale_fill_gradient(low = "#91D1C2FF",
                        high = "#3C5488FF")+
    labs(title = names(best[["enrichList"]])[i])+
    ylab("pathway")+
    theme(axis.text.y = element_blank())+
    coord_flip()
  tiff(filename = paste0("../analysis/Clusters/gene_centric/DEG/gene_centric_best_",supp, "_", names(best[["enrichList"]])[i], ".tiff"),
       width = 6 * 300, 
       height = 0.75 * nrow(clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],]) * 300,
       res = 300,
       compression = "lzw")
  print(g)
  dev.off()
  
  assign(paste0("g", i), g)
  
}



clustobjs <-  as.data.frame(best[["clustObj"]][["cluster"]])
colnames(clustobjs) <-  "cluster"
clustobjs$protein <- rownames(clustobjs)
rownames(clustobjs) <-  NULL
clustobjs <-  clustobjs[c(2,1)]

write.table(clusters_df, paste0("../analysis/Clusters/gene_centric/DEG/gene_centric_enrichlist_", supp, ".txt"), sep = "\t", row.names = FALSE)
write.table(clustobjs, paste0("../analysis/Clusters/gene_centric/DEG/gene_centric_clustobjs_", supp, ".txt"), sep = "\t", row.names = FALSE)

###individual visualizations----
i=8
ggdata <-  clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],]
ggdata <- ggdata[order(ggdata$pvalue,decreasing = TRUE),]
ordered_kinases <- factor(ggdata$kinase,levels = ggdata$kinase)

g <- ggplot(ggdata, 
            aes(x = ordered_kinases, y = size))+
  geom_col(aes(fill = pvalue),width = 0.75)+
  geom_text(aes(label = kinase, y = 0.2), hjust = 0) +
  scale_fill_gradient(low = "#91D1C2FF",
                      high = "#3C5488FF")+
  labs(title = names(best[["enrichList"]])[i])+
  ylab("pathway")+
  theme(axis.text.y = element_blank())+
  coord_flip()
tiff(filename = paste0("../analysis/Clusters/gene_centric/DEG/gene_centric_best_",supp, "_", names(best[["enrichList"]])[i], ".tiff"),
     width = 6 * 300, 
     height = 1.5 * nrow(clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],]) * 300,
     res = 300,
     compression = "lzw")
print(g)
dev.off()

### all in one figure ----
clusters_df$cluster <- clusters_df$cluster %>%
  gsub("cluster ", "", .) %>%
  as.numeric(.)
clusters_df <- clusters_df[order(clusters_df$cluster,decreasing = FALSE),]
clusters_df$cluster <- as.factor(clusters_df$cluster)

tiff(filename = paste0("../analysis/Clusters/gene_centric/DEG/gene_centric_best_",supp, "_ALL",".tiff"),
     width = 14 * 300, 
     height = 15 * 300,
     res = 300,
     compression = "lzw")
ggplot(clusters_df, aes(x = cluster, y = kinase)) +
  geom_point(aes(size = size,
                 fill = cluster, alpha = -log10(pvalue)),
             shape = 21) +
  # scale_fill_manual(values =c("#c95640",
  #                             "#4baf90",
  #                             "#d54795",
  #                             "#71b249",
  #                             "#9a64ca",
  #                             "#ce9944",
  #                             "#6788cc",
  #                             "#7b7f39",
  #                             "#9e4769"), guide = "none")+
  scale_alpha(range = c(0.3, 0.9)) +
  scale_fill_discrete(guide="none")+
  scale_size(range = c(3,8))+
  theme_bw()+
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 0),
        axis.title.x = element_blank(), axis.title.y = element_blank())
dev.off()


#~~~~~~~~~~~~~~~~~~~~
# Site-CENTRIC#######
#~~~~~~~~~~~~~~~~~~~~

#RNGkind("L'Ecuyer-CMRG")

# PhosphoSite.human2 = mapply(function(kinase) {
#   gsub("(.*)(;[A-Z])([0-9]+;)", "\\1;\\3", kinase)
# }, PhosphoSite.human)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## our method all differential sites -----------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df_pl <-as.data.frame(top.norm_intensity.median)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## perform ClueR to identify optimal number of clusters ------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_pl$id <-  rownames(df_pl)

df_pl <- df_pl %>%
  separate(id, c("Protein", "Gene", "Residue.Both", "Peptide", "delete"), ";")

df_pl = df_pl[, c(9,10,1:8,11,13)]
colnames(df_pl)<- c("uniprot_id", "name", "X0000","X0010","X0030","X0060","X0300",
                    "X0600","X0900","X1800", "PSite","ID")
df_pl$ID <- rownames(df_pl)

## collpase same site by mean (sequence sometime slightly different)
df_pl2 <- plyr::ddply(df_pl, .(uniprot_id, PSite), dplyr::summarise,
                      name = max(name),
                      X0000 = mean(X0000),
                      X0010 = mean(X0010),
                      X0030 = mean(X0030),
                      X0060 = mean(X0060),
                      X0300 = mean(X0300),
                      X0600 = mean(X0600),
                      X0900 = mean(X0900),
                      X1800 = mean(X1800),
                      ID = max(ID)
)

rownames(df_pl2) <- paste(df_pl2$name, df_pl2$PSite, "", sep = ";")
df_pl2 <- df_pl2[c(4:11)]
df_pl2 <-  as.matrix(df_pl2)

# perform ClueR to identify optimal number of clusters
set.seed(321)
c3 <- runClue(df_pl2, annotation=PhosphoSite.human, kRange = 2:10, 
              rep = 5, effectiveSize = c(5, 100), pvalueCutoff = 0.05, alpha = 0.5,
              universe = background_list_sites)
c3[["Tc"]] <- as.matrix(c3[["Tc"]]) #without this, sometimes you get empty plots in the best graph

# Visualise the evaluation results
data <- data.frame(Success=as.numeric(c3$evlMat), Freq=rep(2:10, each=5))

myplot <- ggplot(data, aes(x=Freq, y=Success)) + 
  geom_boxplot(aes(x = factor(Freq), fill="gray"))+
  stat_smooth(method="loess", colour="red", linewidth=3, span = 0.5) + 
  xlab("# of cluster")+ ylab("Enrichment score")+theme_classic()
#scale_x_discrete(labels= c(2:20))


myplot

tiff(filename = paste0("../analysis/Clusters/site_centric/DEG/site_centric_boxplot_",supp, ".tiff"),    
     width = 5 * 300, 
     height = 3* 300,
     res = 300,
     compression = "lzw")
myplot
dev.off()

tiff(filename = paste0("../analysis/Clusters/site_centric/DEG/site_centric_best_",supp, ".tiff"), 
     width = 9 * 300, 
     height = 3* 300,
     res = 300,
     compression = "lzw")
set.seed(321)
best <- clustOptimal(c3, rep=5, mfrow=c(1,3), visualize = T,
                     universe = background_list_sites)
dev.off()

#scaled x-axis
tiff(filename = paste0("../analysis/Clusters/site_centric/DEG/site_centric_best_",supp, "_scaled.tiff"), 
     width = 9 * 300, 
     height = 3* 300,
     res = 300,
     compression = "lzw")
set.seed(321)
best <- clustOptimalÖO(c3, rep=5, mfrow=c(1, 3), visualize = "scaled", 
                       #user.maxK = 10, #if you want to decide number of clusters yourself.
                       universe = background_list_sites)
dev.off()


clusters_df <-  data.frame()
for (i in 1:length(best$enrichList)) {
  clusters_df <- rbind(clusters_df, 
                       data.frame(best$enrichList[[i]], cluster = names(best[["enrichList"]])[i]))
}

clusters_df$size <-  as.numeric(clusters_df$size)
clusters_df$pvalue <-  as.numeric(clusters_df$pvalue)

theme_set(theme_cowplot())
colors <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF",
            "#8491B4FF", "#91D1C2FF","#DC0000FF", "#7E6148FF", "#B09C85FF")

for (i in 1:length(best$enrichList)) {
  ggdata <-  clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],]
  ggdata <- ggdata[order(ggdata$pvalue,decreasing = TRUE),]
  ordered_kinases <- factor(ggdata$kinase,levels = ggdata$kinase)
  
  g <- ggplot(ggdata, 
              aes(x = ordered_kinases, y = size))+
    geom_col(aes(fill = pvalue),width = 0.75)+
    geom_text(aes(label = kinase, y = 0.2), hjust = 0) +
    scale_fill_gradient(low = "#91D1C2FF",
                        high = "#3C5488FF")+
    labs(title = names(best[["enrichList"]])[i])+
    ylab("pathway")+
    theme(axis.text.y = element_blank())+
    coord_flip()
  
  tiff(filename = paste0("../analysis/Clusters/site_centric/DEG/site_centric_best_",supp, "_", names(best[["enrichList"]])[i], ".tiff"),
       width = 6 * 300, 
       height = 2 * nrow(clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],]) * 300,
       res = 300,
       compression = "lzw")
  print(g)
  dev.off()
  
  assign(paste0("g", i), g)
  
}

clustobjs <-  as.data.frame(best[["clustObj"]][["cluster"]])
colnames(clustobjs) <-  "cluster"
clustobjs$protein <- rownames(clustobjs)
rownames(clustobjs) <-  NULL
clustobjs <-  clustobjs[c(2,1)]

write.table(clusters_df, paste0("../analysis/Clusters/site_centric/DEG/site_centric_enrichlist_", supp, ".txt"), sep = "\t", row.names = FALSE)
write.table(clustobjs, paste0("../analysis/Clusters/site_centric/DEG/site_centric_clustobjs_", supp, ".txt"), sep = "\t", row.names = FALSE)

###individual visualizations----
i=4
tiff(filename = paste0("analysis/Clusters/site_centric_best_",names(best[["enrichList"]])[i], "_degs.tiff"),
     width = 8 * 300, 
     height = 6* 300,
     res = 300,
     compression = "lzw")
g4
dev.off()

### all in one figure ----
clusters_df$cluster <- clusters_df$cluster %>%
  gsub("cluster ", "", .) %>%
  as.numeric(.)
clusters_df <- clusters_df[order(clusters_df$cluster,decreasing = FALSE),]
clusters_df$cluster <- as.factor(clusters_df$cluster)

tiff(filename = paste0("../analysis/Clusters/site_centric/DEG/site_centric_best_",supp, "_ALL",".tiff"),
     width = 5 * 300, 
     height = 4 * 300,
     res = 300,
     compression = "lzw")
ggplot(clusters_df, aes(x = cluster, y = kinase)) +
  geom_point(aes(size = size,
                 fill = cluster, alpha = -log10(pvalue)),
             shape = 21) +
  scale_fill_manual(values =c("#c95640",
                              "#4baf90",
                              "#d54795",
                              "#71b249",
                              "#9a64ca",
                              "#ce9944",
                              "#6788cc",
                              "#7b7f39",
                              "#9e4769"), guide = "none")+
  scale_alpha(range = c(0.3, 0.9)) +
  #scale_fill_discrete(guide="none")+
  scale_size(range = c(3,8))+
  theme_bw()+
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 0),
        axis.title.x = element_blank(), axis.title.y = element_blank())
dev.off()






#PHOSR WAY - ANOVA####

#~~~~~~~~~~~~~~~~~~~~
## Gene-Centric#######
#~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### prepare pathways ----------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pathways = as.list(reactomePATHID2EXTID)
path_names = as.list(reactomePATHID2NAME)
name_id = match(names(pathways), names(path_names))
names(pathways) = unlist(path_names)[name_id]

pathways = pathways[which(grepl("Homo sapiens", names(pathways), ignore.case = TRUE))]

pathways = lapply(pathways, function(path) {
  gene_name = unname(getSYMBOL(path, data = "org.Hs.eg"))
  toupper(unique(gene_name))
})# perform ClueR to identify optimal number of clusters
names(pathways) <- gsub("_", " ", gsub("REACTOME_|Homo sapiens: ", "", names(pathways)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### phosr tutorial method for dif.phosphosites --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# select differentially phosphorylated sites 

#sites.p <- matANOVA(ppe_norm, grps) #pvals calculated by ANOVA
#ppm <- meanAbundance(ppe_norm, grps)

sites.p <- matANOVA(norm_intensity_filter, grps) #pvals calculated by ANOVA
ppm <- meanAbundance(norm_intensity_filter, grps)

sel <- which((sites.p < 0.05))
#& (rowSums(abs(ppm) > 1) != 0)) --> norm_intensity minus values
ppm_filtered <- ppm[sel,]


# summarise phosphosites information into gene level

ppm_gene_anova <- phosCollapse(ppm_filtered, 
                               sapply(strsplit(rownames(ppm_filtered), ";"), "[[", 2), 
                               stat = apply(abs(ppm_filtered), 1, max), by = "max")


##prepare for heatmaps
ppm_filteredM <- as.data.frame(ppm_filtered)
ppm_gene_anovaM <- as.data.frame(ppm_gene_anova)
full_rows <- apply(ppm_filteredM, 1, paste, collapse="_")
full_rows <- as.data.frame(full_rows)
full_rows$id <- rownames(full_rows)
colnames(full_rows) <- c("values", "id")
full_rows2 <- apply(ppm_gene_anovaM, 1, paste, collapse="_")
full_rows2 <- as.data.frame(full_rows2)
full_rows2$names <- rownames(full_rows2)
colnames(full_rows2) <- c("values", "names")
merged <- merge(full_rows2, full_rows, by = "values")
##




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### perform ClueR to identify optimal number of clusters ------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ppm_gene <-  ppm_gene_anova


#RNGkind("L'Ecuyer-CMRG")

set.seed(19) 
c1 <- runClue(ppm_gene, annotation=pathways, 
              kRange = seq(2,10), rep = 5, effectiveSize = c(5, 100), 
              pvalueCutoff = 0.05, alpha = 0.5,
              universe = background_list,
              standardise = T) # (x-mean)/sd

c1[["Tc"]] <- as.matrix(c1[["Tc"]]) #without this, sometimes you get empty plots in the best graph


# Visualise the evaluation results
data <- data.frame(Success=as.numeric(c1$evlMat), Freq=rep(seq(2,10), each=5))
myplot <- ggplot(data, aes(x=Freq, y=Success)) + 
  geom_boxplot(aes(x = factor(Freq), fill="gray")) +
  stat_smooth(method="loess", colour="red", linewidth=3, span = 0.5) +
  xlab("# of cluster") + 
  ylab("Enrichment score") + 
  #scale_x_discrete(labels= c(10:20))+
  theme_classic()
myplot

supp = "anova"
tiff(filename = paste0("../analysis/Clusters/gene_centric/Anova/gene_centric_boxplot_",
                       supp, ".tiff"),
     width = 5* 300, 
     height = 3* 300,
     res = 300,
     compression = "lzw")
myplot
dev.off()


tiff(filename = paste0("../analysis/Clusters/gene_centric/Anova/gene_centric_best_",
                       supp, ".tiff"),
     width =  15* 300, 
     height = 6* 300,
     res = 300,
     compression = "lzw")
#dev.new()
set.seed(19)
best <- clustOptimal(c1, rep=5, mfrow=c(2, 5), visualize = T, 
                     #user.maxK = 10, #if you want to decide number of clusters yourself.
                     universe = background_list)
dev.off()
#use custom function for scaling x-axis
tiff(filename = paste0("../analysis/Clusters/gene_centric/Anova/gene_centric_best_",
                       supp, "_scaled.tiff"),
     width =  15* 300, 
     height = 6* 300,
     res = 300,
     compression = "lzw")
set.seed(19)
best <- clustOptimalÖO(c1, rep=5, mfrow=c(2, 5), visualize = "scaled", 
                       #user.maxK = 10, #if you want to decide number of clusters yourself.
                       universe = background_list)
dev.off()


clusters_df <-  data.frame()
for (i in 1:length(best$enrichList)) {
  clusters_df <- rbind(clusters_df, 
                       data.frame(best$enrichList[[i]], 
                                  cluster = names(best[["enrichList"]])[i]))
}

clusters_df$size <-  as.numeric(clusters_df$size)
clusters_df$pvalue <-  as.numeric(clusters_df$pvalue)

theme_set(theme_cowplot())


for (i in 1:length(best$enrichList)) {
  ggdata <-  clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],]
  if (nrow(ggdata) < 10) {
    ggdata <- ggdata[order(ggdata$pvalue,decreasing = TRUE),]
    ordered_kinases <- factor(ggdata$kinase,levels = ggdata$kinase)
    g <- ggplot(ggdata, 
                aes(x = ordered_kinases, y = size))+
      geom_col(aes(fill = pvalue))+
      geom_text(aes(label = kinase, y = 0.2), hjust = 0) +
      scale_fill_gradient(low = "#91D1C2FF",
                          high = "#3C5488FF")+
      labs(title = names(best[["enrichList"]])[i])+
      ylab("pathway")+
      theme(axis.text.y = element_blank())+
      coord_flip()
    
    tiff(filename = paste0("../analysis/Clusters/gene_centric/Anova/gene_centric_best_",supp, "_", names(best[["enrichList"]])[i], ".tiff"),
         width = 6 * 300, 
         height = 0.5 * nrow(clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],]) * 300,
         res = 300,
         compression = "lzw")
    print(g)
    dev.off()
    
    assign(paste0("g", i), g)
  } else {
    ggdata <- ggdata[order(ggdata$pvalue,decreasing = F),]
    ggdata <- ggdata[1:10,]
    ggdata <- ggdata[order(ggdata$pvalue,decreasing = TRUE),]
    ordered_kinases <- factor(ggdata$kinase,levels = ggdata$kinase)
    g <- ggplot(ggdata, 
                aes(x = ordered_kinases, y = size))+
      geom_col(aes(fill = pvalue))+
      geom_text(aes(label = kinase, y = 0.2), hjust = 0) +
      scale_fill_gradient(low = "#91D1C2FF",
                          high = "#3C5488FF")+
      labs(title = names(best[["enrichList"]])[i])+
      ylab("pathway")+
      theme(axis.text.y = element_blank())+
      coord_flip()
    
    tiff(filename = paste0("../analysis/Clusters/gene_centric/Anova/gene_centric_best_",supp, "_", names(best[["enrichList"]])[i], ".tiff"),
         width = 6 * 300, 
         height = 4 * 300,
         res = 300,
         compression = "lzw")
    print(g)
    dev.off()
    
    assign(paste0("g", i), g)
  }
}

clustobjs <-  as.data.frame(best[["clustObj"]][["cluster"]])
colnames(clustobjs) <-  "cluster"
clustobjs$protein <- rownames(clustobjs)
rownames(clustobjs) <-  NULL
clustobjs <-  clustobjs[c(2,1)]


write.table(clusters_df, paste0("../analysis/Clusters/gene_centric/Anova/gene_centric_enrichlist_", supp, ".txt"), sep = "\t", row.names = FALSE)
write.table(clustobjs, paste0("../analysis/Clusters/gene_centric/Anova/gene_centric_clustobjects_", supp, ".txt"), sep = "\t", row.names = FALSE)


####individual visualizations----
names(best[["enrichList"]])

i=6
tiff(filename = paste0("../analysis/Clusters/gene_centric/Anova/gene_centric_best_",supp, "_", names(best[["enrichList"]])[i], ".tiff"),
     width = 6 * 300, 
     height = 1 * nrow(clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],]) * 300,
     res = 300,
     compression = "lzw")
ggdata <-  clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],]
ggdata <- ggdata[order(ggdata$pvalue,decreasing = TRUE),]
ordered_kinases <- factor(ggdata$kinase,levels = ggdata$kinase)

ggplot(ggdata, 
       aes(x = ordered_kinases, y = size))+
  geom_col(aes(fill = pvalue))+
  geom_text(aes(label = kinase, y = 0.2), hjust = 0) +
  scale_fill_gradient(low = "#91D1C2FF",
                      high = "#3C5488FF")+
  labs(title = names(best[["enrichList"]])[i])+
  ylab("pathway")+
  theme(axis.text.y = element_blank())+
  coord_flip()
dev.off()

### all in one figure ----
clusters_df$cluster <- clusters_df$cluster %>%
  gsub("cluster ", "", .) %>%
  as.numeric(.)
clusters_df <- clusters_df[order(clusters_df$cluster,decreasing = FALSE),]
clusters_df$cluster <- as.factor(clusters_df$cluster)

tiff(filename = paste0("../analysis/Clusters/gene_centric/Anova/gene_centric_best_",supp, "_ALL",".tiff"),
     width = 14 * 300, 
     height = 15 * 300,
     res = 300,
     compression = "lzw")
ggplot(clusters_df, aes(x = cluster, y = kinase)) +
  geom_point(aes(size = size,
                 fill = cluster, alpha = -log10(pvalue)),
             shape = 21) +
  # scale_fill_manual(values =c("#c95640",
  #                             "#4baf90",
  #                             "#d54795",
  #                             "#71b249",
  #                             "#9a64ca",
  #                             "#ce9944",
  #                             "#6788cc",
  #                             "#7b7f39",
  #                             "#9e4769"), guide = "none")+
  scale_alpha(range = c(0.3, 0.9)) +
  scale_fill_discrete(guide="none")+
  scale_size(range = c(3,8))+
  theme_bw()+
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 0),
        axis.title.x = element_blank(), axis.title.y = element_blank())
dev.off()


#### if too many, get the top5 from each ----
clusters_df_sub <-  data.frame()
#i=4
for (i in 1:length(best$enrichList)) {
  clusters_df_sub <- rbind(clusters_df_sub, 
                           data.frame(best$enrichList[[i]], cluster = names(best[["enrichList"]])[i])[1:5,])
}
clusters_df_sub$size <-  as.numeric(clusters_df_sub$size)
clusters_df_sub$pvalue <-  as.numeric(clusters_df_sub$pvalue)

##continue 
input <- clusters_df_sub
input <- na.omit(input)

input$cluster <- input$cluster %>%
  gsub("cluster ", "", .) %>%
  as.numeric(.)
input <- input[order(input$cluster,decreasing = FALSE),]
input$cluster <- as.factor(input$cluster)

tiff(filename = paste0("../analysis/Clusters/gene_centric/Anova/gene_centric_best_",supp, "_ALL_Sub",".tiff"),
     width = 10 * 300, 
     height = 10 * 300,
     res = 300,
     compression = "lzw")
ggplot(input, aes(x = cluster, y = kinase)) +
  geom_point(aes(size = size,
                 fill = cluster, alpha = -log10(pvalue)),
             shape = 21) +
  # scale_fill_manual(values =c("#c95640",
  #                     "#4baf90",
  #                     "#d54795",
  #                     "#71b249",
  #                     "#9a64ca",
  #                     "#ce9944",
  #                     "#6788cc",
  #                     "#7b7f39",
  #                     "#9e4769"), guide = "none")+
  scale_alpha(range = c(0.3, 0.9)) +
  scale_fill_discrete(guide="none")+
  scale_size(range = c(3,8))+
  theme_bw()+
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 0),
        axis.title.x = element_blank(), axis.title.y = element_blank())
dev.off()
#~~~~~~~~~~~~~~~~~~~~
## Site-Centric#######
#~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### phosr tutorial method for dif.phosphosites --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df_pl <-  as.data.frame(ppm_filtered)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### perform ClueR to identify optimal number of clusters ------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_pl$id <-  rownames(df_pl)

df_pl <- df_pl %>%
  separate(id, c("Protein", "Gene", "Residue.Both", "Peptide", "delete"), ";")

df_pl = df_pl[, c(9,10,1:8,11,13)]
colnames(df_pl)<- c("uniprot_id", "name", "X0000","X0010","X0030","X0060","X0300",
                    "X0600","X0900","X1800", "PSite","ID")
df_pl$ID <- rownames(df_pl)

## collpase same site by mean (sequence sometime slightly different)
df_pl2 <- plyr::ddply(df_pl, .(uniprot_id, PSite), dplyr::summarise,
                      name = max(name),
                      X0000 = mean(X0000),
                      X0010 = mean(X0010),
                      X0030 = mean(X0030),
                      X0060 = mean(X0060),
                      X0300 = mean(X0300),
                      X0600 = mean(X0600),
                      X0900 = mean(X0900),
                      X1800 = mean(X1800),
                      ID = max(ID)
)

rownames(df_pl2) <- paste(df_pl2$name, df_pl2$PSite, "", sep = ";")
df_pl2 <- df_pl2[c(4:11)]
df_pl2 <-  as.matrix(df_pl2)


##prepare for heatmaps
ppm_gene_anovaM <- as.data.frame(df_pl2)
full_rows2 <- apply(ppm_gene_anovaM, 1, paste, collapse="_")
full_rows2 <- as.data.frame(full_rows2)
full_rows2$names <- rownames(full_rows2)
colnames(full_rows2) <- c("values", "names")
merged2 <- merge(full_rows2, full_rows, by = "values")
##


# perform ClueR to identify optimal number of clusters
set.seed(321)
c3 <- runClue(df_pl2, annotation=PhosphoSite.human, 
              kRange = seq(2,10), rep = 5, effectiveSize = c(5, 100), 
              pvalueCutoff = 0.05, alpha = 0.5,
              universe = background_list_sites)


c3[["Tc"]] <- as.matrix(c3[["Tc"]]) #without this, sometimes you get empty plots in the best graph


# Visualise the evaluation results
data <- data.frame(Success=as.numeric(c3$evlMat), Freq=rep(seq(2,10), each=5))
myplot <- ggplot(data, aes(x=Freq, y=Success)) + 
  geom_boxplot(aes(x = factor(Freq), fill="gray")) +
  stat_smooth(method="loess", colour="red", linewidth=3, span = 0.5) +
  xlab("# of cluster") + 
  ylab("Enrichment score") + 
  #scale_x_discrete(labels= c(2:20))+
  theme_classic()
myplot

tiff(filename = paste0("../analysis/Clusters/site_centric/Anova/site_centric_boxplot_",supp, ".tiff"),    
     width = 5 * 300, 
     height = 3* 300,
     res = 300,
     compression = "lzw")
myplot
dev.off()

tiff(filename = paste0("../analysis/Clusters/site_centric/Anova/site_centric_best_",supp, ".tiff"),    
     width = 9 * 300, 
     height = 3* 300,
     res = 300,
     compression = "lzw")
set.seed(321)
best <- clustOptimal(c3, rep=1000, mfrow=c(2,6), visualize = T, user.maxK = 12,
                     universe = background_list_sites)
dev.off()

# let's try more clusters 
set.seed(321)
best <- clustOptimalÖO(c3, rep=15, mfrow=c(2, 6), 
                     user.maxK = 12,
                     visualize = "scaled", 
                     universe = background_list_sites)

tiff(filename = paste0("../analysis/Clusters/site_centric/Anova/site_centric_best_",supp, "_scaled.tiff"),    
     width = 3 * 300, 
     height = 9 * 300,
     res = 300,
     compression = "lzw")
set.seed(321)
best <- clustOptimalÖO(c3, rep=5, mfrow=c(3, 1), visualize = "scaled", 
                       #user.maxK = 10, #if you want to decide number of clusters yourself.
                       universe = background_list_sites)
dev.off()



clusters_df <-  data.frame()
for (i in 1:length(best$enrichList)) {
  clusters_df <- rbind(clusters_df, 
                       data.frame(best$enrichList[[i]], cluster = names(best[["enrichList"]])[i]))
}

clusters_df$size <-  as.numeric(clusters_df$size)
clusters_df$pvalue <-  as.numeric(clusters_df$pvalue)

theme_set(theme_cowplot())

for (i in 1:length(best$enrichList)) {
  ggdata <-  clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],]
  ggdata <- ggdata[order(ggdata$pvalue,decreasing = TRUE),]
  ordered_kinases <- factor(ggdata$kinase,levels = ggdata$kinase)
  
  g <- ggplot(ggdata, 
              aes(x = ordered_kinases, y = size))+
    geom_col(aes(fill = pvalue))+
    geom_text(aes(label = kinase, y = 0.2), hjust = 0) +
    scale_fill_gradient(low = "#91D1C2FF",
                        high = "#3C5488FF")+
    labs(title = names(best[["enrichList"]])[i])+
    ylab("pathway")+
    theme(axis.text.y = element_blank())+
    coord_flip()
  
  tiff(filename = paste0("../analysis/Clusters/site_centric/Anova/site_centric_best_",supp, "_", names(best[["enrichList"]])[i], ".tiff"),
       width = 4 * 300, 
       height = 1 * nrow(clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],]) * 300,
       res = 300,
       compression = "lzw")
  print(g)
  dev.off()
  
  assign(paste0("g", i), g)
  
}

clustobjs <-  as.data.frame(best[["clustObj"]][["cluster"]])
colnames(clustobjs) <-  "cluster"
clustobjs$protein <- rownames(clustobjs)
rownames(clustobjs) <-  NULL
clustobjs <-  clustobjs[c(2,1)]

write.table(clusters_df, paste0("../analysis/Clusters/site_centric/Anova/site_centric_enrichlist_", supp, ".txt"), sep = "\t", row.names = FALSE)
write.table(clustobjs, paste0("../analysis/Clusters/site_centric/Anova/site_centric_clustobjs_", supp, ".txt"), sep = "\t", row.names = FALSE)

####individual visualizations----


i=1
tiff(filename = paste0("../analysis/Clusters/site-centric/Anova/site_centric_best_",supp, "_", names(best[["enrichList"]])[i], ".tiff"),
     width = 4 * 300, 
     height = 2* 300,
     res = 300,
     compression = "lzw")
ggdata <-  clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],]
ggdata <- ggdata[order(ggdata$pvalue,decreasing = TRUE),]
ordered_kinases <- factor(ggdata$kinase,levels = ggdata$kinase)

g <- ggplot(ggdata, 
            aes(x = ordered_kinases, y = size))+
  geom_col(aes(fill = pvalue))+
  geom_text(aes(label = kinase, y = 0.2), hjust = 0) +
  scale_fill_gradient(low = "#91D1C2FF",
                      high = "#3C5488FF")+
  labs(title = names(best[["enrichList"]])[i])+
  ylab("pathway")+
  theme(axis.text.y = element_blank())+
  coord_flip()

print(g)
dev.off()

### all in one figure ----
clusters_df$cluster <- clusters_df$cluster %>%
  gsub("cluster ", "", .) %>%
  as.numeric(.)
clusters_df <- clusters_df[order(clusters_df$cluster,decreasing = FALSE),]
clusters_df$cluster <- as.factor(clusters_df$cluster)

tiff(filename = paste0("../analysis/Clusters/site_centric/Anova/site_centric_best_",supp, "_ALL",".tiff"),
     width = 5 * 300, 
     height = 4 * 300,
     res = 300,
     compression = "lzw")
ggplot(clusters_df, aes(x = cluster, y = kinase)) +
  geom_point(aes(size = size,
                 fill = cluster, alpha = -log10(pvalue)),
             shape = 21) +
  scale_fill_manual(values =c("#c95640",
                              "#4baf90",
                              "#d54795",
                              "#71b249",
                              "#9a64ca",
                              "#ce9944",
                              "#6788cc",
                              "#7b7f39",
                              "#9e4769"), guide = "none")+
  scale_alpha(range = c(0.3, 0.9)) +
  #scale_fill_discrete(guide="none")+
  scale_size(range = c(3,8))+
  theme_bw()+
  theme(axis.text = element_text(size = 14), axis.text.x = element_text(angle = 0),
        axis.title.x = element_blank(), axis.title.y = element_blank())
dev.off()
