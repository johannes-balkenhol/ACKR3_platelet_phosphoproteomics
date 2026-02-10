
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


#load phosphosite etc. data (downloaded from github KSEA) IMPORTANT: use fold change, not log2-transformed
#KSData = read.csv("C:/Users/Ozge/Nextcloud/Phosphoproteom/phosphoproteom validation/analysis/KSEA/PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv") #ÖO
KSData = read.csv("../../analysis/KSEA/PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv") ##JB



dfs_input = list(top.filter.10, top.filter.600, top.filter.1800, 
                 top.filter.10.dmso.vs.0s, top.filter.600.dmso.vs.0s, top.filter.1800.dmso.vs.0s,
                 top.filter.10.cxcr7.vs.0s, top.filter.600.cxcr7.vs.0s, top.filter.1800.cxcr7.vs.0s)

names_input = c("10", "600", "1800", 
                "10.dmso.vs.0s", "600.dmso.vs.0s", "1800.dmso.vs.0s",
                "10.cxcr7.vs.0s", "600.cxcr7.vs.0s", "1800.cxcr7.vs.0s")


#i=2
#load input data

for (i in 1:length(dfs_input)) {
  
  my_df1 <- dfs_input[[i]]
  my_df1$Peptide <- sapply(strsplit(rownames(my_df1), ";"), "[[", 4)
  my_df1$FC <- 2^my_df1$logFC
  my_df1 <- my_df1[c(1,2,7,6,5,8)]
  colnames(my_df1) <- c("Protein", "Gene", "Peptide", "Residue.Both", "p", "FC"  )
  PX<-my_df1
  PX$Residue.Both <- str_replace_all(PX$Residue.Both, fixed("|"), ";" )
  rownames(PX)<- NULL
  
  KSData.dataset <- KSEA.KS_table(KSData, PX, NetworKIN=FALSE, NetworKIN.cutoff=5)
  assign(paste0("KSData.dataset", names_input[[i]]), KSData.dataset)
  write.table(eval(as.name(paste0("KSData.dataset", names_input[[i]]))), file = paste0("../../analysis/KSEA/", names_input[i], "_Kinase_Substrate_Links.csv"), sep = ",", quote = F, 
              row.names = F, col.names = T)
  
  
  Scores <- KSEA.Scores(KSData, PX, NetworKIN=FALSE, NetworKIN.cutoff=5)
  assign(paste0("Scores", names_input[[i]]), Scores)
  write.table(eval(as.name(paste0("Scores", names_input[[i]]))), file = paste0("../../analysis/KSEA/", names_input[i], "_KSEA_Kinase_Scores.csv"), sep = ",", quote = F, 
              row.names = F, col.names = T)
  
  
  #Barplot <- KSEA.Barplot(KSData, PX, NetworKIN=FALSE, NetworKIN.cutoff=5, m.cutoff=5, p.cutoff=0.05, export=FALSE)
  
  tiff(filename = paste0("../../analysis/KSEA/", names_input[i] ,"_barplot.tiff"),
       width = 5 * 300, 
       height = 5 * 300,
       res = 300,
       compression = "lzw")
  
  KSEA.Barplot(KSData, PX, NetworKIN=FALSE, NetworKIN.cutoff=5, m.cutoff=5, p.cutoff=0.05, export=FALSE)
  
  dev.off()
  
}

KSEA.Heatmap(list(Scores10, Scores600, Scores1800),
             sample.labels = c("10", "600", "1800"),
             m.cutoff = 5, stats="p.value",p.cutoff=0.05, sample.cluster=F)






