
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

setwd("/Users/ozgeosmanoglu/Nextcloud/Phosphoproteom/phosphoproteom analysis 2/")


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
    print(as.list(test[[i]]))
    top.names <- append(top.names, rownames(input2[[i]][input2[[i]]$PValue<0.05 & abs(input2[[i]]$logFC)>0,]))
  }
}

norm_intensity2 = norm_intensity[!(row.names(norm_intensity) %in% not.uniq[,1]), ]
norm_intensity2 <- as.data.frame(norm_intensity2)

top.norm_intensity <- norm_intensity2[top.names,]

top.norm_intensity<- top.norm_intensity[!duplicated(top.norm_intensity), ]

top.norm_intensity <- as.matrix(top.norm_intensity)


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
                         sapply(strsplit(rownames(top.norm_intensity), ";"), "[[", 2), 
                         stat = apply(abs(top.norm_intensity), 1, max), by = "max")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### perform ClueR to identify optimal number of clusters ------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ppm_gene <- ppm_gene_degs

RNGkind("L'Ecuyer-CMRG")
set.seed(123) 

c1 <- runClue(ppm_gene, annotation=pathways, 
              kRange = seq(2,10), rep = 5, effectiveSize = c(5, 100), 
              pvalueCutoff = 0.05, alpha = 0.5)

# Visualise the evaluation results
data <- data.frame(Success=as.numeric(c1$evlMat), Freq=rep(seq(2,10), each=5))
myplot <- ggplot(data, aes(x=Freq, y=Success)) + 
  geom_boxplot(aes(x = factor(Freq), fill="gray")) +
  stat_smooth(method="loess", colour="red", linewidth=3, span = 0.5) +
  xlab("# of cluster") + 
  ylab("Enrichment score") + 
  theme_classic()
myplot


tiff(filename = paste0("analysis/Clusters/gene_centric_boxplot","_degs", ".tiff"),
    width = 5* 300, 
    height = 3* 300,
    res = 300,
    compression = "lzw")
myplot
dev.off()


tiff(filename = paste0("analysis/Clusters/gene_centric_best","_degs", ".tiff"),
     width = 8 * 300, 
     height = 6* 300,
     res = 300,
     compression = "lzw")
set.seed(123)
best <- clustOptimal(c1, rep=5, mfrow=c(2, 3), visualize = TRUE)
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
  g <- ggplot(clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],], 
         aes(x = kinase, y = size))+
    geom_col(aes(fill = pvalue),width = 0.75)+
    geom_text(aes(label = kinase, y = 0.2), hjust = 0) +
    scale_fill_gradient(low = "#91D1C2FF",
                        high = "#3C5488FF")+
    labs(title = names(best[["enrichList"]])[i])+
    ylab("pathway")+
    theme(axis.text.y = element_blank())+
    coord_flip()
  
  tiff(filename = paste0("analysis/Clusters/gene_centric_best_",names(best[["enrichList"]])[i], "_degs.tiff"),
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

write.table(clusters_df, "analysis/Clusters/gene_centric_degs_enrichlist.txt", sep = "\t", row.names = FALSE)
write.table(clustobjs, "analysis/Clusters/gene_centric_degs_clustobjects.txt", sep = "\t", row.names = FALSE)

###individual visualizations----
i=4
tiff(filename = paste0("analysis/Clusters/gene_centric_best_",names(best[["enrichList"]])[i], "_degs.tiff"),
     width = 8 * 300, 
     height = 0.6 * nrow(clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],]) * 300,
     res = 300,
     compression = "lzw")
g4
dev.off()




#~~~~~~~~~~~~~~~~~~~~
# Site-CENTRIC#######
#~~~~~~~~~~~~~~~~~~~~

RNGkind("L'Ecuyer-CMRG")

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
              rep = 5, effectiveSize = c(5, 100), pvalueCutoff = 0.05, alpha = 0.5)

# Visualise the evaluation results
data <- data.frame(Success=as.numeric(c3$evlMat), Freq=rep(2:10, each=5))

myplot <- ggplot(data, aes(x=Freq, y=Success)) + geom_boxplot(aes(x = factor(Freq), fill="gray"))+
  stat_smooth(method="loess", colour="red", linewidth=3, span = 0.5) + xlab("# of cluster")+ ylab("Enrichment score")+theme_classic()



tiff(filename = paste0("analysis/Clusters/site_centric_boxplot","_degs", ".tiff"),
     width = 5 * 300, 
     height = 3* 300,
     res = 300,
     compression = "lzw")
myplot
dev.off()


tiff(filename = paste0("analysis/Clusters/site_centric_best","_degs", ".tiff"),
     width = 8 * 300, 
     height = 6* 300,
     res = 300,
     compression = "lzw")
set.seed(321)
best <- clustOptimal(c3, rep=15, mfrow=c(2, 3), visualize = TRUE)
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
  g <- ggplot(clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],], 
              aes(x = kinase, y = size))+
    geom_col(aes(fill = pvalue),width = 0.75)+
    geom_text(aes(label = kinase, y = 0.2), hjust = 0) +
    scale_fill_gradient(low = "#91D1C2FF",
                        high = "#3C5488FF")+
    labs(title = names(best[["enrichList"]])[i])+
    ylab("pathway")+
    theme(axis.text.y = element_blank())+
    coord_flip()
  
  tiff(filename = paste0("analysis/Clusters/site_centric_best_",names(best[["enrichList"]])[i], "_degs.tiff"),
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

write.table(clusters_df, "analysis/Clusters/site_centric_degs_enrichlist.txt", sep = "\t", row.names = FALSE)
write.table(clustobjs, "analysis/Clusters/site_centric_degs_clustobjects.txt", sep = "\t", row.names = FALSE)

###individual visualizations----


i=4
tiff(filename = paste0("analysis/Clusters/site_centric_best_",names(best[["enrichList"]])[i], "_degs.tiff"),
     width = 8 * 300, 
     height = 6* 300,
     res = 300,
     compression = "lzw")
g4
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

ppe_norm2 <- ppe_norm[!(rownames(ppe_norm) %in% not.uniq[,1]),]
sites.p <- matANOVA(ppe_norm2, grps) #pvals calculated by ANOVA
ppm <- meanAbundance(ppe_norm2, grps)

sel <- which((sites.p < 0.05)
             & (rowSums(abs(ppm) > 1) != 0))
ppm_filtered <- ppm[sel,]

# summarise phosphosites information into gene level
ppm_gene_anova <- phosCollapse(ppm_filtered, 
                               sapply(strsplit(rownames(ppm_filtered), ";"), "[[", 2), 
                               stat = apply(abs(ppm_filtered), 1, max), by = "max")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### perform ClueR to identify optimal number of clusters ------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ppm_gene <-  ppm_gene_anova

RNGkind("L'Ecuyer-CMRG")
set.seed(123) 

c1 <- runClue(ppm_gene, annotation=pathways, 
              kRange = seq(2,10), rep = 5, effectiveSize = c(5, 100), 
              pvalueCutoff = 0.05, alpha = 0.5)

# Visualise the evaluation results
data <- data.frame(Success=as.numeric(c1$evlMat), Freq=rep(seq(2,10), each=5))
myplot <- ggplot(data, aes(x=Freq, y=Success)) + 
  geom_boxplot(aes(x = factor(Freq), fill="gray")) +
  stat_smooth(method="loess", colour="red", linewidth=3, span = 0.5) +
  xlab("# of cluster") + 
  ylab("Enrichment score") + 
  theme_classic()
myplot


tiff(filename = paste0("analysis/Clusters/gene_centric_boxplot","_anova", ".tiff"),
     width = 5* 300, 
     height = 3* 300,
     res = 300,
     compression = "lzw")
myplot
dev.off()


tiff(filename = paste0("analysis/Clusters/gene_centric_best","_anova", ".tiff"),
     width = 8 * 300, 
     height = 6* 300,
     res = 300,
     compression = "lzw")
set.seed(123)
best <- clustOptimal(c1, rep=5, mfrow=c(2, 3), visualize = TRUE)
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
  g <- ggplot(clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],], 
              aes(x = kinase, y = size))+
    geom_col(aes(fill = pvalue),width = 0.75)+
    geom_text(aes(label = kinase, y = 0.2), hjust = 0) +
    scale_fill_gradient(low = "#91D1C2FF",
                        high = "#3C5488FF")+
    labs(title = names(best[["enrichList"]])[i])+
    ylab("pathway")+
    theme(axis.text.y = element_blank())+
    coord_flip()
  
  tiff(filename = paste0("analysis/Clusters/gene_centric_best_",names(best[["enrichList"]])[i], "_anova.tiff"),
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

write.table(clusters_df, "analysis/Clusters/gene_centric_anova_enrichlist.txt", sep = "\t", row.names = FALSE)
write.table(clustobjs, "analysis/Clusters/gene_centric_anova_clustobjects.txt", sep = "\t", row.names = FALSE)


####individual visualizations----
i=4
tiff(filename = paste0("analysis/Clusters/gene_centric_best_",names(best[["enrichList"]])[i], "_anova.tiff"),
     width = 8 * 300, 
     height = 0.6 * nrow(clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],]) * 300,
     res = 300,
     compression = "lzw")
g4
dev.off()




#~~~~~~~~~~~~~~~~~~~~
## Site-Centric#######
#~~~~~~~~~~~~~~~~~~~~

RNGkind("L'Ecuyer-CMRG")

# PhosphoSite.human2 = mapply(function(kinase) {
#   gsub("(.*)(;[A-Z])([0-9]+;)", "\\1;\\3", kinase)
# }, PhosphoSite.human)

Depod.human <- read.csv("Depod_phosphatases.txt", sep = "\t")
Depod.human <- Depod.human[is.na(Depod.human$Dephosphosites) == FALSE,]
Depod.human <- Depod.human[c(1:917),]
Depod.human <- Depod.human[Depod.human$allto != "",]

mylist <-  list()

list("ACP1" = Depod.human[Depod.human$Phosphataseentrynames == "ACP1","substrate"])
for (i in Depod.human$Phosphataseentrynames){
  myinput <- list(Depod.human[Depod.human$kinasesubs == Depod.human[Depod.human$Phosphataseentrynames ==i, "kinasesubs"],"allto"])
  names(myinput) <- as.character(i)
  mylist <- append(mylist, myinput)
  
}

temp<-unique(mylist)

Depod.human$Phosphataseentrynames=Depod.human[Depod.human$Phosphataseentrynames,"substrate"])

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

# perform ClueR to identify optimal number of clusters
set.seed(321)
c3 <- runClue(df_pl2, annotation=PhosphoSite.human, kRange = 2:10, 
              rep = 5, effectiveSize = c(5, 100), pvalueCutoff = 0.05, alpha = 0.5)

# Visualise the evaluation results
data <- data.frame(Success=as.numeric(c3$evlMat), Freq=rep(2:10, each=5))

myplot <- ggplot(data, aes(x=Freq, y=Success)) + geom_boxplot(aes(x = factor(Freq), fill="gray"))+
  stat_smooth(method="loess", colour="red", linewidth=3, span = 0.5) + xlab("# of cluster")+ ylab("Enrichment score")+theme_classic()



tiff(filename = paste0("analysis/Clusters/site_centric_boxplot","_anova", ".tiff"),
     width = 5 * 300, 
     height = 3* 300,
     res = 300,
     compression = "lzw")
myplot
dev.off()


tiff(filename = paste0("analysis/Clusters/site_centric_best","_anova", ".tiff"),
     width = 8 * 300, 
     height = 6* 300,
     res = 300,
     compression = "lzw")
set.seed(321)
best <- clustOptimal(c3, rep=15, mfrow=c(2, 3), visualize = TRUE)
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
  g <- ggplot(clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],], 
              aes(x = kinase, y = size))+
    geom_col(aes(fill = pvalue),width = 0.75)+
    geom_text(aes(label = kinase, y = 0.2), hjust = 0) +
    scale_fill_gradient(low = "#91D1C2FF",
                        high = "#3C5488FF")+
    labs(title = names(best[["enrichList"]])[i])+
    ylab("pathway")+
    theme(axis.text.y = element_blank())+
    coord_flip()
  
  tiff(filename = paste0("analysis/Clusters/site_centric_best_",
                         names(best[["enrichList"]])[i], "_anova.tiff"),
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

write.table(clusters_df, "analysis/Clusters/site_centric_anova_enrichlist.txt", sep = "\t", row.names = FALSE)
write.table(clustobjs, "analysis/Clusters/site_centric_anova_clustobjects.txt", sep = "\t", row.names = FALSE)

####individual visualizations----


i=4
tiff(filename = paste0("analysis/Clusters/site_centric_best_",names(best[["enrichList"]])[i], "_anova.tiff"),
     width = 8 * 300, 
     height = 6* 300,
     res = 300,
     compression = "lzw")
g4
dev.off()

