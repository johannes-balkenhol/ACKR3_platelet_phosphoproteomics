# Load libraries ----
suppressPackageStartupMessages({
  library(annotate)
  library(ClueR)
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(ggsci)
  library(org.Hs.eg.db)
  library(parallel)
  library(pheatmap)
  library(PhosR)
  library(plyr)
  library(RColorBrewer)
  library(reactome.db)
  library(tidyr)
})

data(PhosphoSitePlus)

background_list <-unique(sapply(strsplit(rownames(norm_intensity_collapse), ";"), "[[", 2))
background_list_sites <-unique(paste0(sapply(strsplit(rownames(norm_intensity_filter), ";"), "[[", 2), 
                                      ";", sapply(strsplit(rownames(norm_intensity_filter), ";"), "[[", 3), ";"))


#setwd("/Users/ozgeosmanoglu/Nextcloud/Phosphoproteom/phosphoproteom validation 2/scripts")
#setwd("C:/Users/job37yv/Research/Phospho/phospho vaildation 2/scripts/")  ## JB Laptop RVZ


#OUR WAY - no test####

#~~~~~~~~~~~~~~~~~~~~
## Gene-Centric#######
#~~~~~~~~~~~~~~~~~~~~


#OUR WAY - DEG####

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

## 1. use all logFCs ####
supp = "all_FC_no_test" ##define which method
folder = "allFC"
logFC_table <- top.all
logFC_table$ID <-  rownames(logFC_table)
logFC_table<- logFC_table[!duplicated(logFC_table), ]
logFC_table <-  logFC_table[grep("logFC", colnames(logFC_table))]
top.norm_intensity.median <- logFC_table
ppm_gene_all <- phosCollapse(top.norm_intensity.median, 
                             sapply(strsplit(rownames(top.norm_intensity.median), ";"), "[[", 2),
                             stat = apply(abs(top.norm_intensity.median), 1, max), by = "max")
ppm_gene_all <- as.data.frame(ppm_gene_all)
ppm_gene_input <- ppm_gene_all #change according to method
logFC.00 <- 1
ppm_gene_input <- cbind(logFC.00, ppm_gene_input)

###end###



#####################################################################################

### 2. get significantly regulated phosphopeptides ####

input = list(top.filter.10, top.filter.600, top.filter.1800, 
             top.filter.10.dmso.vs.0s, top.filter.600.dmso.vs.0s, top.filter.1800.dmso.vs.0s,
             top.filter.10.cxcr7.vs.0s, top.filter.600.cxcr7.vs.0s, top.filter.1800.cxcr7.vs.0s)

names_input = c("10", "600", "1800", 
                "10.dmso.vs.0s", "600.dmso.vs.0s", "1800.dmso.vs.0s",
                "10.cxcr7.vs.0s", "600.cxcr7.vs.0s", "1800.cxcr7.vs.0s")


### select top differential sites
input2 <- input[c(1:3)]

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

top.names <- top.names[!duplicated(top.names)]


### 2.1. use FCs of DEGs ####
supp = "deg_FC" ##define which method
folder = "DEG/FC"
top.logfc <- top.all[top.names,]
top.logfc$ID <-  rownames(top.logfc)
top.logfc<- top.logfc[!duplicated(top.logfc), ]
top.logfc <-  top.logfc[grep("logFC", colnames(top.logfc))]
top.norm_intensity.median <- top.logfc
ppm_gene_FCs <- phosCollapse(top.norm_intensity.median, 
                             sapply(strsplit(rownames(top.norm_intensity.median), ";"), "[[", 2),
                             stat = apply(abs(top.norm_intensity.median), 1, max), by = "max")
ppm_gene_FCs <- as.data.frame(ppm_gene_FCs) #use this for logfcs

#take top names, use log2FCs
ppm_gene_input <- ppm_gene_FCs #change according to method
logFC.00 <- 1
ppm_gene_input <- cbind(logFC.00, ppm_gene_input)
###end###



### 2.2. use norm. abundances of DEGs ####

norm_intensity2 <- as.data.frame(norm_intensity_filter)
top.norm_intensity <- norm_intensity2[top.names,]
top.norm_intensity$ID <-  rownames(top.norm_intensity)
top.norm_intensity<- top.norm_intensity[!duplicated(top.norm_intensity), ]
top.norm_intensity <- as.matrix(top.norm_intensity[1:70])

########### aggregate by mean per group 
x_tt <- as.factor(grps)
t00 <- rowMedians(top.norm_intensity[,1:10])
t10 <- rowMedians(top.norm_intensity[,11:20])
t10wt <- rowMedians(top.norm_intensity[,21:30])
t600 <- rowMedians(top.norm_intensity[,31:40])
t600wt <- rowMedians(top.norm_intensity[,41:50])
t1800 <- rowMedians(top.norm_intensity[,51:60])
t1800wt <- rowMedians(top.norm_intensity[,61:70])

top.norm_intensity.median <- cbind(t00,t10,t10wt,t600,t600wt,t1800,t1800wt)
rownames(top.norm_intensity.median) <- rownames(top.norm_intensity)
colnames(top.norm_intensity.median) <- c("t00_ctrl","t10_cxcr7","t10_dmso","t600_cxcr7","t600_dmso","t1800_cxcr7","t1800_dmso")
top.norm_intensity.median <- as.data.frame(top.norm_intensity.median)

##collapse for names
ppm_gene <- phosCollapse(top.norm_intensity.median, 
                         sapply(strsplit(rownames(top.norm_intensity.median), ";"), "[[", 2),
                         stat = apply(abs(top.norm_intensity.median), 1, max), by = "max")
ppm_gene <- as.data.frame(ppm_gene)

# Select columns with names containing "cxcr7" or "ctrl"


#### 2.2.1 take top names, DMSO abundances----
ppm_gene_dmso <- dplyr::select(ppm_gene, contains("ctrl") | contains("dmso"))
ppm_gene_input <- ppm_gene_dmso #change according to method
supp = "deg_dmso" #describe method
folder = "DEG/DMSO"

#### 2.2.2. take top names, CXCR7 abundances----
ppm_gene_cxcr7 <- dplyr::select(ppm_gene, contains("ctrl") | contains("cxcr7"))
ppm_gene_input <- ppm_gene_cxcr7 #change according to method
supp = "deg_cxcr7" #describe method
folder = "DEG/CXCR7"




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 3. perform ClueR to identify optimal number of clusters ------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set.seed(123) 
c1 <- runClue(ppm_gene_input, annotation=pathways, 
              kRange = seq(10,20), rep = 5, effectiveSize = c(5, 100), 
              pvalueCutoff = 0.05, alpha = 0.5, 
              universe = background_list)

c1[["Tc"]] <- as.matrix(c1[["Tc"]]) #without this, sometimes you get empty plots in the best graph

# ##all annotations are already in c1, not overrepresented or clustered but just annotation of whatever sigficantly abundant
# all_hemostasis <- c1[["annotation"]][["Hemostasis"]]
# all_platelet_degranulation <- c1[["annotation"]][["Platelet degranulation "]]
# all_platelet_signaling <- c1[["annotation"]][["Platelet activation, signaling and aggregation"]]

# Visualise the evaluation results
data <- data.frame(Success=as.numeric(c1$evlMat), Freq=rep(seq(1,11), each=5))
myplot <- ggplot(data, aes(x=Freq, y=Success)) + 
  geom_boxplot(aes(x = factor(Freq), fill="gray")) +
  stat_smooth(method="loess", colour="red", linewidth=3, span = 0.5) +
  xlab("# of cluster") + 
  ylab("Enrichment score") + 
  scale_x_discrete(labels= c(10:20))+
  theme_classic()
myplot

#be careful with the path
tiff(filename = paste0("../analysis/Clusters/gene_centric/", folder, "/gene_centric_boxplot_",supp, ".tiff"),
     width = 5* 300, 
     height = 3* 300,
     res = 300,
     compression = "lzw")
myplot
dev.off()


#always arrange the size so that one plot is 3x3.
tiff(filename = paste0("../analysis/Clusters/gene_centric/", folder, "/gene_centric_best_",supp, ".tiff"),
     width = 12 * 300, 
     height = 9* 300,
     res = 300,
     compression = "lzw")

set.seed(123) 
best <- clustOptimal(c1, rep=5, mfrow=c(3,4), visualize = T,
                     universe = background_list) #default is the rownames of the dataset.

dev.off()



##
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
    geom_col(aes(fill = pvalue))+
    geom_text(aes(label = kinase, y = 0.2), hjust = 0) +
    scale_fill_gradient(low = "#91D1C2FF",
                        high = "#3C5488FF")+
    labs(title = names(best[["enrichList"]])[i])+
    ylab("pathway")+
    theme(axis.text.y = element_blank())+
    coord_flip()
  tiff(filename = paste0("../analysis/Clusters/gene_centric/", folder, "/gene_centric_best_",supp, "_", names(best[["enrichList"]])[i], ".tiff"),
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

write.table(clusters_df, paste0("../analysis/Clusters/gene_centric/", folder,  "/gene_centric_enrichlist_", supp, ".txt"), sep = "\t", row.names = FALSE)
write.table(clustobjs, paste0("../analysis/Clusters/gene_centric/", folder, "/gene_centric_clustobjects_", supp, ".txt"), sep = "\t", row.names = FALSE)

### 4. individual visualizations----
i=7
ggdata <-  clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],]
ggdata <- ggdata[order(ggdata$pvalue,decreasing = TRUE),]
ordered_kinases <- factor(ggdata$kinase,levels = ggdata$kinase)

tiff(filename = paste0("../analysis/Clusters/gene_centric/", folder, "/gene_centric_best_",supp, "_", names(best[["enrichList"]])[i], ".tiff"),
     width = 8 * 300, 
     height = 1.5 * nrow(clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],]) * 300,
     res = 300,
     compression = "lzw")

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


### 1. all logFC ####
supp = "all_FC_no_test" ##define which method
folder = "allFC"
df_pl <- logFC_table #all logFC
df_pl$id <-  rownames(df_pl)
df_pl <- df_pl %>%
  separate(id, c("Protein", "Gene", "Residue.Both", "Peptide", "delete"), ";")

df_pl = df_pl[, c(4,5,1:3,6,8)]
colnames(df_pl)<- c("uniprot_id", "name", "X0010", "X0600","X1800", "PSite","ID")
df_pl2 <- plyr::ddply(df_pl, .(uniprot_id, PSite), dplyr::summarise,
                      name = max(name),
                      X0010= mean(X0010),
                      X0600 = mean(X0600),
                      X1800 = mean(X1800),
                      ID = max(ID)
)
rownames(df_pl2) <- paste(df_pl2$name, df_pl2$PSite, "", sep = ";")
df_pl2 <- df_pl2[c(4:6)]
df_pl2_input <- df_pl2
logFC.00 <- 1
df_pl2_input <- cbind(logFC.00, df_pl2_input)
##end


### 2. DEG FCs ####
supp = "deg_FC" ##define which method
folder = "DEG/FC"
df_pl <- top.logfc 
df_pl$id <-  rownames(df_pl)
df_pl <- df_pl %>%
  separate(id, c("Protein", "Gene", "Residue.Both", "Peptide", "delete"), ";")

df_pl = df_pl[, c(4,5,1:3,6,8)]
colnames(df_pl)<- c("uniprot_id", "name", "X0010", "X0600","X1800", "PSite","ID")
df_pl2 <- plyr::ddply(df_pl, .(uniprot_id, PSite), dplyr::summarise,
                      name = max(name),
                      X0010= mean(X0010),
                      X0600 = mean(X0600),
                      X1800 = mean(X1800),
                      ID = max(ID)
)
rownames(df_pl2) <- paste(df_pl2$name, df_pl2$PSite, "", sep = ";")
df_pl2 <- df_pl2[c(4:6)]
df_pl2_input <- df_pl2
logFC.00 <- 1
df_pl2_input <- cbind(logFC.00, df_pl2_input)
##end


### 3. for norm_intensity df ####
df_pl <-as.data.frame(top.norm_intensity.median) ##different for logFc and norm abundances
df_pl$id <-  rownames(df_pl)

df_pl <- df_pl %>%
  separate(id, c("Protein", "Gene", "Residue.Both", "Peptide", "delete"), ";")

df_pl = df_pl[, c(8,9,1:7,10,12)]
colnames(df_pl)<- c("uniprot_id", "name", "X0000","X0010_CXCR7","X0010_DMSO","X0600_CXCR7","X0600_DMSO",
                    "X1800_CXCR7","X1800_DMSO", "PSite","ID")
df_pl$ID <- rownames(df_pl)
df_pl2 <- plyr::ddply(df_pl, .(uniprot_id, PSite), dplyr::summarise,
                      name = max(name),
                      X0000 = mean(X0000),
                      X0010_CXCR7 = mean(X0010_CXCR7),
                      X0010_DMSO = mean(X0010_DMSO),
                      X0600_CXCR7 = mean(X0600_CXCR7),
                      X0600_DMSO = mean(X0600_DMSO),
                      X1800_CXCR7 = mean(X1800_CXCR7),
                      X1800_DMSO = mean(X1800_DMSO),
                      ID = max(ID)
)

rownames(df_pl2) <- paste(df_pl2$name, df_pl2$PSite, "", sep = ";")
df_pl2 <- df_pl2[c(4:10)]

#### 3.1. CXCR7 ----
df_pl2_cxcr7 <- df_pl2[c(1,2,4,6)] 
df_pl2_cxcr7 <-  as.matrix(df_pl2_cxcr7)
df_pl2_input <- df_pl2_cxcr7
supp = "deg_cxcr7" #describe method
folder = "DEG/CXCR7"


#### 3.2. DMSO ----
df_pl2_dmso <- df_pl2[c(1,3,5,7)]
df_pl2_dmso <-as.matrix(df_pl2_dmso)
df_pl2_input <- df_pl2_dmso
supp = "deg_dmso" #describe method
folder = "DEG/DMSO"



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 4. perform ClueR to identify optimal number of clusters ------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set.seed(321)
c3 <- runClue(df_pl2_input, annotation=PhosphoSite.human, kRange = 2:10, 
              rep = 5, effectiveSize = c(3, 100), pvalueCutoff = 0.05, alpha = 0.5,
              universe = names(gene_list_pl))
c3[["Tc"]] <- as.matrix(c3[["Tc"]])

# Visualise the evaluation results
data <- data.frame(Success=as.numeric(c3$evlMat), Freq=rep(2:10, each=5))

myplot <- ggplot(data, aes(x=Freq, y=Success)) + geom_boxplot(aes(x = factor(Freq), fill="gray"))+
  stat_smooth(method="loess", colour="red", linewidth=3, span = 0.5) + xlab("# of cluster")+ ylab("Enrichment score")+theme_classic()

myplot
#be careful with the path
tiff(filename = paste0("../analysis/Clusters/site_centric/", folder, "/site_centric_boxplot_",supp, ".tiff"),
     width = 5 * 300, 
     height = 3* 300,
     res = 300,
     compression = "lzw")
myplot
dev.off()

tiff(filename = paste0("../analysis/Clusters/site_centric/", folder,"/site_centric_best_",supp, ".tiff"),
     width = 3 * 300, 
     height = 9* 300,
     res = 300,
     compression = "lzw")
set.seed(321)
best <- clustOptimal(c3, rep=15, mfrow=c(3,1), visualize = TRUE, universe = names(gene_list_pl))
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
  tiff(filename = paste0("../analysis/Clusters/site_centric/", folder, "/site_centric_best_",supp, "_", names(best[["enrichList"]])[i], ".tiff"),
       width = 6 * 300, 
       height = 1.5 * nrow(clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],]) * 300,
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

write.table(clusters_df, paste0("../analysis/Clusters/site_centric/", folder, "/site_centric_enrichlist_", supp, ".txt"), sep = "\t", row.names = FALSE)
write.table(clustobjs, paste0("../analysis/Clusters/site_centric/", folder, "/site_centric_clustobjects_", supp, ".txt"), sep = "\t", row.names = FALSE)

### individual visualizations----

i=4
tiff(filename = paste0("analysis/Clusters/site_centric_best_",names(best[["enrichList"]])[i], "_degs.tiff"),
     width = 8 * 300, 
     height = 6* 300,
     res = 300,
     compression = "lzw")
g4
dev.off()









####################################################################################################
####################################################################################################

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

### 1. select differentially phosphorylated sites ----
ppe_norm <-  as.data.frame(norm_intensity_filter)

#separate for cxcr7 and dmso
ppe_norm_cxcr7 <- dplyr::select(ppe_norm, contains("ctrl") | contains("cxcr7"))
ppe_norm_cxcr7 <- ppe_norm_cxcr7[c(1:10, 11:20, 31:40, 21:30)]
ppe_norm_dmso <- dplyr::select(ppe_norm, contains("ctrl") | contains("dmso"))
ppe_norm_dmso <- ppe_norm_dmso[c(1:10, 11:20, 31:40, 21:30)]
grps_dmso <- grps[c(1:10,21:30, 61:70, 41:50)]
grps_cxcr7 <- grps[c(1:10, 11:20, 51:60, 31:40)]

#### 1.1. for dmso ----
supp = "anova_dmso" ##define which method
folder = "Anova/DMSO"
sites.p <- matANOVA(ppe_norm_dmso, grps_dmso) #pvals calculated by ANOVA
ppm <- meanAbundance(ppe_norm_dmso, grps_dmso)
sel <- which((sites.p < 0.05))
#& (rowSums(abs(ppm) > 1) != 0))
ppm_filtered <- ppm[sel,]


#### 1.2. for cxcr7----
supp = "anova_cxcr7" ##define which method
folder = "Anova/CXCR7"
sites.p <- matANOVA(ppe_norm_cxcr7, grps_cxcr7) #pvals calculated by ANOVA
ppm <- meanAbundance(ppe_norm_cxcr7, grps_cxcr7)
sel <- which((sites.p < 0.05))
#& (rowSums(abs(ppm) > 1) != 0))
ppm_filtered <- ppm[sel,]


#### 1.3. alltogether - trying intensity difference norm_int(cxcr7)-norm_int(ctrl)----
supp = "anova_diff_all" ##define which method
folder = "Anova/Diff"

unique_donors <- unique(sample_info$donor)
unique_time <- unique(sample_info$time[11:70])
result_diffs <- data.frame(matrix(NA, nrow = nrow(ppe_norm), ncol = 0))

for(time in unique_time) {
  for (donor in unique_donors) {
    control_cols <- grep(paste0("x" , time, "sek_", "DMSO", "_", donor), colnames(ppe_norm))
    treatment_cols <- grep(paste0("x" , time, "sek_", "CXCR7", "_", donor), colnames(ppe_norm))
    diffs <- ppe_norm[, treatment_cols, drop = FALSE] - ppe_norm[, control_cols, drop = FALSE]
    colnames(diffs) <- paste0("x",time, "sek", "_DIFF_", donor)
    result_diffs <- cbind(result_diffs, diffs)
  }
}

ppe_norm_diff <- merge(ppe_norm[1:10], result_diffs, by = "row.names", all = TRUE)
row.names(ppe_norm_diff) <-  ppe_norm_diff$Row.names
ppe_norm_diff <- ppe_norm_diff %>% dplyr::select(-Row.names)
ppe_norm_diff[1:10] <- 0
#write.table(ppe_norm_diff, "../data/processed_data/ppe_norm_diff.txt", sep = "\t")
grps_diff <- gsub("_[0-9][0-9]", "", colnames(ppe_norm_diff))

sites.p <- matANOVA(ppe_norm_diff, grps_diff) #pvals calculated by ANOVA
ppm <- meanAbundance(ppe_norm_diff, grps_diff)
sel <- which((sites.p < 0.05))
# & (rowSums(abs(ppm) > 0) != 0))
ppm_filtered <- ppm[sel,]



##### 1.3.2 or no anova, all diffs----
folder = "Diff_all"
supp = "diff_all"
ppm_filtered <- ppm #comes from 1.3.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### 2. perform ClueR to identify optimal number of clusters ------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### summarise phosphosites information into gene level ----
ppm_gene_anova <- phosCollapse(ppm_filtered, 
                               sapply(strsplit(rownames(ppm_filtered), ";"), "[[", 2), 
                               stat = apply(abs(ppm_filtered), 1, max), by = "max")

ppm_gene <-  ppm_gene_anova


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


set.seed(123) 
c1 <- runClue(ppm_gene, annotation=pathways, 
              kRange = seq(2,20), rep = 5, effectiveSize = c(5, 100), 
              pvalueCutoff = 0.05, alpha = 0.5,
              universe = background_list,
              standardise = T)
c1[["Tc"]] <- as.matrix(c1[["Tc"]]) #without this, sometimes you get empty plots in the best graph


# Visualise the evaluation results
data <- data.frame(Success=as.numeric(c1$evlMat), Freq=rep(seq(1,19), each=5))
myplot <- ggplot(data, aes(x=Freq, y=Success)) + 
  geom_boxplot(aes(x = factor(Freq), fill="gray")) +
  stat_smooth(method="loess", colour="red", linewidth=3, span = 0.5) +
  xlab("# of cluster") + 
  ylab("Enrichment score") + 
  scale_x_discrete(labels= c(2:20))+
  theme_classic()
myplot

tiff(filename = paste0("../analysis/Clusters/gene_centric/", folder, "/gene_centric_boxplot_",supp, ".tiff"),
     width = 5* 300, 
     height = 3* 300,
     res = 300,
     compression = "lzw")
myplot
dev.off()

tiff(filename = paste0("../analysis/Clusters/gene_centric/", folder, "/gene_centric_best_",supp, ".tiff"),
     width = 12* 300, 
     height = 6* 300,
     res = 300,
     compression = "lzw")
set.seed(123)
best <- clustOptimal(c1, rep=5, mfrow=c(2, 4), visualize = T,
                     universe = background_list)
dev.off()

tiff(filename = paste0("../analysis/Clusters/gene_centric/", folder, "/gene_centric_best_",supp, "_scaled.tiff"),
     width = 12* 300, 
     height = 6* 300,
     res = 300,
     compression = "lzw")
set.seed(123)
best <- clustOptimalÖO(c1, rep=5, mfrow=c(2, 4), visualize = "scaled", 
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

i=1

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
    tiff(filename = paste0("../analysis/Clusters/gene_centric/", folder, "/gene_centric_best_",supp, "_", names(best[["enrichList"]])[i], ".tiff"),
         width = 6 * 300, 
         height = 0.75 * nrow(clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],]) * 300,
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
    assign(paste0("g", i), g)
    
    tiff(filename = paste0("../analysis/Clusters/gene_centric/", folder, "/gene_centric_best_",supp, "_", names(best[["enrichList"]])[i], ".tiff"),
         width = 6 * 300, 
         height = 6 * 300,
         res = 300,
         compression = "lzw")
    print(g)
    dev.off()
  }
}


clustobjs <-  as.data.frame(best[["clustObj"]][["cluster"]])
colnames(clustobjs) <-  "cluster"
clustobjs$protein <- rownames(clustobjs)
rownames(clustobjs) <-  NULL
clustobjs <-  clustobjs[c(2,1)]


write.table(clusters_df, paste0("../analysis/Clusters/gene_centric/", folder, "/gene_centric_enrichlist_", supp, ".txt"), sep = "\t", row.names = FALSE)
write.table(clustobjs, paste0("../analysis/Clusters/gene_centric/", folder, "/gene_centric_clustobjects_", supp, ".txt"), sep = "\t", row.names = FALSE)


### 3. individual visualizations----
names(best[["enrichList"]])
i = 1
tiff(filename = paste0("../analysis/Clusters/gene_centric/", folder, "/gene_centric_best_",supp, "_", names(best[["enrichList"]])[i], ".tiff"),
     width = 5 * 300, 
     height = 3 * 300,
     res = 300,
     compression = "lzw")
ggdata <-  clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],]
ggdata <- ggdata[order(ggdata$pvalue,decreasing = TRUE),]
ordered_kinases <- factor(ggdata$kinase,levels = ggdata$kinase)
#ggdata$kinase[1] <- "Diseases of signal transduction by growth\nfactor receptors and second messengers"
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
print(g)
dev.off()


### all in one figure ----

#if too many, get the top5 from each
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

#for the diffall
input$kinase[1] <- "Regulation of IGF transport and uptake by IGFBPs"

#[1] "Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)"

tiff(filename = paste0("../analysis/Clusters/gene_centric/", folder, "/gene_centric_best_",supp, "_Sub",".tiff"),
     width = 12 * 300, 
     height = 8 * 300,
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

### heatmaps of clusters ----
dataset <- norm_intensity_filter
dataset_df <- as.data.frame(dataset)
name <- sapply(strsplit(rownames(dataset_df), ";"), "[[", 2)
site <- sapply(strsplit(rownames(dataset_df), ";"), "[[", 3)
dataset_df$namesite <- paste0(name,";",site,";")

top.all.input <- top.all


protein_function <- function(description, proteins, cellH) {
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
  #targets_norm_abundance2 <-  targets_norm_abundance
  
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
           color = colorRampPalette(c("navy", "white", "red"))(100),
           #gaps_row=c(5,10,15,20,25,30,35,40,45,50, 55),
           scale="row",
           na_col = "grey",
           #breaks = mat_breaks,
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
           clustering_distance_rows="correlation",
           clustering_distance_cols="euclidean",
           clustering_method="complete",
           main = gsub("\\.", " ", description),
           display_numbers = test_labels,
           number_color = "green", 
           fontsize_number = 8,
           cellheight=cellH, cellwidth = 30)
  return(heatmap)
  
  
  
}

x <- "cluster 4"
prot_spec_anova_cluster <- as.data.frame(best$enrichList[[x]]) #CHANGE CLUSTER NAME
prot_spec_anova_cluster <- prot_spec_anova_cluster[7,]
cluster_prots <- unique(unlist(strsplit(prot_spec_anova_cluster$substrates, split = "\\|")))
description <- prot_spec_anova_cluster$kinase #just pathway name
#description <- paste(x, prot_spec_anova_cluster[rev(order(prot_spec_anova_cluster$size)),"kinase"][1])
proteinlist <- paste(cluster_prots, collapse = ";|") #to get all sites
proteins <- merged[merged$names %in% cluster_prots, "id"] #to get max sites

# to get the whole pathway
allprots <- pathways[[prot_spec_anova_cluster$kinase]]
allprots   <-  pathways[["Downstream TCR signaling"]]
description <- "cGMP effects"
proteinlist <- paste(allprots, collapse = ";|")

top.all.input <- top.all

protein_function(description, proteins, 10)

# Plot heatmap
tiff(filename = paste0("../analysis/Clusters/gene_centric/", folder, "/Heatmap_",gsub(" ", "_", description),".tiff"),
     width = 10 * 300, 
     height = 15 * 300,
     res = 300,
     compression = "lzw")
protein_function(description, proteins, 10)
dev.off()

# ## from GO ORA
# library(org.Hs.eg.db)
# library(GO.db)
# 
# union_sig_df <- data.frame(ids = union_sig, 
#                            protein = sapply(strsplit(union_sig, ";"), "[[", 2))

# # or you get the whole GO term of interest and we plot everything that belongs to it
# library(stringr)
# 
# # get the go term
# interesting_GOs <-  list("ERBB signaling" ="GO:0038127",
#                          "second-messenger-mediated signaling" = "GO:0019932",
#                          "cytoskeleton organization" = "GO:0007010",
#                          "negative regulation of apoptotic signaling pathway" = "GO:2001234")
# 
# results <- AnnotationDbi::select(org.Hs.eg.db, keys=c(interesting_GOs[[4]]), columns = c('SYMBOL'), keytype = "GOALL")
# symbols <- unique(results$SYMBOL)
# protein <- paste(symbols, collapse = ";|")
# 
# description <- names(interesting_GOs[4])


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

#foo <- data.frame(do.call('rbind', strsplit(as.character(df_pl$id),';',fixed=TRUE)))
#
df_pl <- cbind(df_pl, foo) %>% select(-id)
df_pl <- df_pl %>%
 separate(id, c("Protein", "Gene", "Residue.Both", "Peptide", "delete"), ";")

df_pl = df_pl[, c(5,6,1:4,7,9)]
colnames(df_pl)<- c("uniprot_id", "name", "X0000","X0010","X0600", "X1800","PSite","ID")
df_pl$ID <- rownames(df_pl)

#df_pl <- df_pl[df_pl$name %in% cluster8_genes$protein,]

## collpase same site by mean (sequence sometime slightly different)
df_pl2 <- plyr::ddply(df_pl, .(uniprot_id, PSite), dplyr::summarise,
                      name = max(name),
                      X0000 = mean(X0000),
                      X0010 = mean(X0010),
                      X0600 = mean(X0600),
                      X1800 = mean(X1800),
                      ID = max(ID)
)

rownames(df_pl2) <- paste(df_pl2$name, df_pl2$PSite, "", sep = ";")
df_pl2 <- df_pl2[c(4:7)]
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
c3 <- runClue(df_pl2, annotation=PhosphoSite.human, kRange = 2:20, 
              rep = 5, effectiveSize = c(3, 100), pvalueCutoff = 0.05, alpha = 0.5,
              universe = background_list_sites, 
              standardise = T)
#for diff effective size change to c(3,100)

# Visualise the evaluation results
data <- data.frame(Success=as.numeric(c3$evlMat), Freq=rep(2:20, each=5))

myplot <- ggplot(data, aes(x=Freq, y=Success)) + geom_boxplot(aes(x = factor(Freq), fill="gray"))+
  stat_smooth(method="loess", colour="red", linewidth=3, span = 0.5) + xlab("# of cluster")+ ylab("Enrichment score")+theme_classic()

myplot

tiff(filename = paste0("../analysis/Clusters/site_centric/", folder, "/site_centric_boxplot_",supp, ".tiff"),
     width = 5 * 300, 
     height = 3* 300,
     res = 300,
     compression = "lzw")
myplot
dev.off()


tiff(filename = paste0("../analysis/Clusters/site_centric/", folder, "/site_centric_best_",supp, ".tiff"),
     width = 12 * 300, 
     height = 3* 300,
     res = 300,
     compression = "lzw")
set.seed(321)
best <- clustOptimal(c3, rep=5, mfrow=c(1,4), visualize = T, 
                     universe = background_list_sites)
dev.off()

tiff(filename = paste0("../analysis/Clusters/site_centric/", folder, "/site_centric_best_",supp, "_scaled.tiff"),
     width = 12 * 300, 
     height = 3* 300,
     res = 300,
     compression = "lzw")
set.seed(321)
best <- clustOptimalÖO(c3, rep=5, mfrow=c(1,4), visualize = "scaled", 
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
    geom_col(aes(fill = pvalue),width = 0.75)+
    geom_text(aes(label = kinase, y = 0.2), hjust = 0) +
    scale_fill_gradient(low = "#91D1C2FF",
                        high = "#3C5488FF")+
    labs(title = names(best[["enrichList"]])[i])+
    ylab("pathway")+
    theme(axis.text.y = element_blank())+
    coord_flip()
  tiff(filename = paste0("../analysis/Clusters/site_centric/", folder, "/site_centric_best_",supp, "_", names(best[["enrichList"]])[i], ".tiff"),
       width = 4 * 300, 
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


write.table(clusters_df, paste0("../analysis/Clusters/site_centric/", folder, "/site_centric_enrichlist_", supp, ".txt"), sep = "\t", row.names = FALSE)
write.table(clustobjs, paste0("../analysis/Clusters/site_centric/", folder, "/site_centric_clustobjects_", supp, ".txt"), sep = "\t", row.names = FALSE)


####individual visualizations----
names(best[["enrichList"]])
i=1
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
  ylab("")+
  xlab("")+
  theme(axis.text.y = element_blank())+
  coord_flip()

tiff(filename = paste0("../analysis/Clusters/site_centric/", folder, "/site_centric_best_",supp, "_", names(best[["enrichList"]])[i], ".tiff"),
     width = 4 * 300, 
     height = 2 * nrow(clusters_df[clusters_df$cluster == names(best[["enrichList"]])[i],]) * 300,
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

clusters_df

tiff(filename = paste0("../analysis/Clusters/site_centric/", folder, "/site_centric_best_",supp, "_ALL",".tiff"),
     width = 5 * 300, 
     height = 4 * 300,
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

### heatmaps of clusters ----

x <- "cluster 2"
prot_spec_anova_cluster <- as.data.frame(best$enrichList$`cluster 2`)[1,] #CHANGE CLUSTER NAME
cluster_prots <- unique(unlist(strsplit(prot_spec_anova_cluster$substrates, split = "\\|")))
description <- x
description <- paste(x, prot_spec_anova_cluster[rev(order(prot_spec_anova_cluster$size)),"kinase"][1])
proteins <- merged2[merged2$names %in% cluster_prots, "id"] #merged2 for site-centric
top.all.input <- top.all

tiff(filename = paste0("../analysis/Clusters/gene_centric/", folder, "/Heatmap_",gsub(" ", "_", description),".tiff"),
     width = 10 * 300, 
     height = 15 * 300,
     res = 300,
     compression = "lzw")
protein_function(description, proteins)
dev.off()
