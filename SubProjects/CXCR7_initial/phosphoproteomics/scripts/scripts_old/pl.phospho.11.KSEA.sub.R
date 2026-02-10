
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
KSData = read.csv("C:/Users/Ozge/Nextcloud/Phosphoproteom/phosphoproteom validation/analysis/KSEA/PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv")
#KSData = read.csv("analysis/KSEA/PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv")

setwd("C:/Users/Ozge/Nextcloud/Phosphoproteom/phosphoproteom analysis 2/")


names_input = c("10", "30", "60", "300", "600", "900", "1800")
dfs_input = list(top.filter.10, top.filter.30, top.filter.60, top.filter.300, top.filter.600, top.filter.900, top.filter.1800)



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
write.table(eval(as.name(paste0("KSData.dataset", names_input[[i]]))), file = paste0("analysis/KSEA/", names_input[i], "_Kinase_Substrate_Links.csv"), sep = ",", quote = F, 
            row.names = F, col.names = T)


Scores <- KSEA.Scores(KSData, PX, NetworKIN=FALSE, NetworKIN.cutoff=5)
assign(paste0("Scores", names_input[[i]]), Scores)
write.table(eval(as.name(paste0("Scores", names_input[[i]]))), file = paste0("analysis/KSEA/", names_input[i], "_KSEA_Kinase_Scores.csv"), sep = ",", quote = F, 
            row.names = F, col.names = T)


#Barplot <- KSEA.Barplot(KSData, PX, NetworKIN=FALSE, NetworKIN.cutoff=5, m.cutoff=5, p.cutoff=0.05, export=FALSE)

tiff(filename = paste0("analysis/KSEA/", names_input[i] ,"_barplot.tiff"),
     width = 5 * 300, 
     height = 5 * 300,
     res = 300,
     compression = "lzw")

KSEA.Barplot(KSData, PX, NetworKIN=FALSE, NetworKIN.cutoff=5, m.cutoff=5, p.cutoff=0.05, export=FALSE)

dev.off()

}





