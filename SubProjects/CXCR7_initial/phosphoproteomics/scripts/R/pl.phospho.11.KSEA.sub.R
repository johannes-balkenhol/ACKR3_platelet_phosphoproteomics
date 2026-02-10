
############################################
#### script for annotation of phosphoproteom data 
## Kinase Substrate enrichment analysis
##



#install.packages("KSEAapp")
library(KSEAapp)
library(plyr)
library(dplyr)
library("tidyr")
library(stringr)
library(ggplot2)
library("ggsci")
library("cowplot")

#load phosphosite etc. data (downloaded from github KSEA) IMPORTANT: use fold change, not log2-transformed
KSData = read.csv("/Users/ozgeosmanoglu/Nextcloud/CXCR7_platelet_analysis/CXCR7_validation/phosphoproteomics/analysis/KSEA/PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv")
#KSData = read.csv("analysis/KSEA/PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv")

#setwd("/Users/ozgeosmanoglu/Nextcloud/Phosphoproteom/phosphoproteom analysis 2/scripts")


names_input = c("10", "30", "60", "300", "600", "900", "1800")
dfs_input = list(top.filter.10, top.filter.30, top.filter.60, top.filter.300, top.filter.600, top.filter.900, top.filter.1800)

for (i in 1:length(dfs_input)) {
  
  my_df1 <- dfs_input[[i]]
  my_df1$Peptide <- sapply(strsplit(rownames(my_df1), ";"), "[[", 4)
  my_df1$FC <- 2^my_df1$logFC
  my_df1 <- my_df1[c(1,2,7,6,5,8)]
  colnames(my_df1) <- c("Protein", "Gene", "Peptide", "Residue.Both", "p", "FC"  )
  PX<-my_df1
  PX$Residue.Both <- str_replace_all(PX$Residue.Both, fixed("|"), ";" )
  rownames(PX)<- NULL
  
  KSData.dataset <- KSEA.KS_table(KSData, PX, NetworKIN=FALSE)
  assign(paste0("KSData.dataset", names_input[[i]]), KSData.dataset)
  write.table(eval(as.name(paste0("KSData.dataset", names_input[[i]]))), file = paste0("../analysis/KSEA/", names_input[i], "_Kinase_Substrate_Links.csv"), sep = ",", quote = F, 
              row.names = F, col.names = T)
  
  Scores <- KSEA.Scores(KSData, PX, NetworKIN=FALSE)
  assign(paste0("Scores", names_input[[i]]), Scores)
  write.table(eval(as.name(paste0("Scores", names_input[[i]]))), file = paste0("../analysis/KSEA/", names_input[i], "_KSEA_Kinase_Scores.csv"), sep = ",", quote = F, 
              row.names = F, col.names = T)
  
  networkinput <- KSData.dataset[KSData.dataset$Kinase.Gene %in% Scores[Scores$p.value<0.05,"Kinase.Gene"],]
  assign(paste0("networkinput", names_input[[i]]), networkinput)
  networkinput$pValue <- "0.001"
  write.table(networkinput[c(3,5,6,1,2)], file = paste0("../analysis/KSEA/", names_input[i], "_networkinput.csv"), sep = ",",
              row.names = F, col.names = T)
  
  tiff(filename = paste0("../analysis/KSEA/", names_input[i] ,"_barplot2.tiff"),
       width = 5 * 300, 
       height = 5 * 300,
       res = 300,
       compression = "lzw")
  
  KSEA.Barplot(KSData, PX, NetworKIN=FALSE, NetworKIN.cutoff=5, m.cutoff=5, p.cutoff=0.05, export=FALSE)
  
  dev.off()
  
}


KSEA.Heatmap(list(Scores10, Scores30, Scores60, Scores300, Scores600, Scores900, Scores1800),
             sample.labels = c("10", "30", "60", "300", "600", "900", "1800"),
             m.cutoff = 5, stats="p.value",p.cutoff=0.05, sample.cluster=F)


KSEASignKinases <- unique(c(Scores10[Scores10$m >= 5 & Scores10$p.value < 0.05,"Kinase.Gene"],
                            Scores30[Scores30$m >= 5 & Scores30$p.value < 0.05,"Kinase.Gene"],
                            Scores60[Scores60$m >= 5 & Scores60$p.value < 0.05,"Kinase.Gene"],
                            Scores300[Scores300$m >= 5 & Scores300$p.value < 0.05,"Kinase.Gene"],
                            Scores600[Scores600$m >= 5 & Scores600$p.value < 0.05,"Kinase.Gene"],
                            Scores900[Scores900$m >= 5 & Scores900$p.value < 0.05,"Kinase.Gene"],
                            Scores1800[Scores1800$m >= 5 & Scores1800$p.value < 0.05,"Kinase.Gene"]))



###plot specific kinase activities over time

theme_set(theme_cowplot())


#merge all timepoints
all_kinase_scores <- data.frame()

for (i in 1:length(names_input)) {
  df <- eval(as.name(paste("Scores", names_input[i], sep="")))
  df$time <- names_input[i]
  all_kinase_scores <- rbind(all_kinase_scores, df)
}

#get scores for kinases of interest    

kinases <- list("PRKACA", "PRKG1", "PRKCA", "PRKAA1","AKT1", "MTOR", "CDK5", "CDK1","CDK2", "MAPKAPK2")
kinases <- kinases[c(1:5)] #for pks
kinases <- kinases[c(6:9)] #for ups

Kinase_scores <-  data.frame()
for (i in 1:length(kinases)) {
  Kinase_scores <- rbind(Kinase_scores, 
                         all_kinase_scores[all_kinase_scores$Kinase.Gene == kinases[[i]],])
}

Kinase_scores$time <-  factor(Kinase_scores$time)
Kinase_scores$time <- ordered(Kinase_scores$time, levels = c("10","30","60","300","600","900","1800"))

#library(wesanderson)


pks <- ggplot(Kinase_scores, aes(x = time, y = z.score, group = Kinase.Gene, color = Kinase.Gene))+
  geom_line(size = 1.1) +
  geom_point(size=4, alpha = 1) +
  geom_point(data = Kinase_scores[Kinase_scores$p.value < 0.05, ], 
             shape = "*", size=6, color = "white")+
  #labs(title = paste(kinase, "KSEA score"))+
  theme(plot.title = element_text(hjust = 0.5, size = 16, ))+
  #scale_colour_manual(values = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF"))+
  scale_color_manual(values = c("#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF"))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray")+
  scale_y_continuous(name ="z-score", 
                     limits = c(-3,3),
                     breaks=c(-3, -2, -1, 0, 1, 2, 2 ,3))

ups <- ggplot(Kinase_scores, aes(x = time, y = z.score, group = Kinase.Gene, color = Kinase.Gene))+
  geom_line(size = 1.1) +
  geom_point(size=4, alpha = 1) +
  geom_point(data = Kinase_scores[Kinase_scores$p.value < 0.05, ], 
             shape = "*", size=6, color = "white")+
  #labs(title = paste(kinase, "KSEA score"))+
  theme(plot.title = element_text(hjust = 0.5, size = 16, ))+
  scale_colour_manual(values = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF"))+
  #scale_color_manual(values = c("#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF"))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray")+
  scale_y_continuous(name ="z-score", 
                     limits = c(-3,3),
                     breaks=c(-3, -2, -1, 0, 1, 2, 2 ,3))

tiff(filename = paste0("../analysis/KSEA/pks", ".tiff"),
     width = 8 * 300, 
     height = 6* 300,
     res = 300,
     compression = "lzw")
pks
dev.off() 

tiff(filename = paste0("../analysis/KSEA/ups", ".tiff"),
     width = 8 * 300, 
     height = 6* 300,
     res = 300,
     compression = "lzw")
ups
dev.off()

#scale_color_manual(values=wes_palette(n=4, name="GrandBudapest2"))


# pal_npg("nrc")(10)
[1] "#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF" "#8491B4FF" "#91D1C2FF"
[8] "#DC0000FF" "#7E6148FF" "#B09C85FF"



