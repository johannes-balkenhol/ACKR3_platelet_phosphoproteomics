
############################################
#### script for annotation of phosphoproteom data 
## Kinase Substrate enrichment analysis
##



#install.packages("KSEAapp")

setwd("C:/Users/Ozge/Nextcloud/Phosphoproteom/phosphoproteom validation/analysis/KSEA/")


library(KSEAapp)
library(plyr)
library(dplyr)
library("tidyr")
library(stringr)


#load phosphosite etc. data (downloaded from github KSEA) IMPORTANT: use fold change, not log2-transformed
KSData = read.csv("C:/Users/Ozge/Nextcloud/Phosphoproteom/phosphoproteom validation/analysis/KSEA/PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv")
KSData = read.csv("../analysis/KSEA/PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv")

#load input data
desc <- "top.filter.10"
my_df1 <- subset(top.1800, select = c(id, logFC, adj.P.Val))
rownames(my_df1) <- NULL
my_df1 <- my_df1 %>%
  separate(id, c("Protein", "Gene", "Residue.Both", "Peptide", "delete"), ";")
my_df1$FC <- 2^my_df1$logFC
## get rid of site without annotation
my_df1 = my_df1[-which(my_df1$Residue.Both == ""),]

my_df1 <- my_df1[c(1,2,4,3,7,8)]

colnames(my_df1) <- c("Protein", "Gene", "Peptide", "Residue.Both", "p", "FC"  )

## collpase same site by mean (sequence sometime slightly different)
my_df1<- plyr::ddply(my_df1, .(Protein, Residue.Both), dplyr::summarise,
              Gene = max(Gene),
              Peptide = Peptide,
              #Average = mean(Average),
              p = mean(p),
              FC = mean(FC)
              )

my_df1 <- my_df1[!duplicated(my_df1), ]

PX<-my_df1

PX <- PX[c(1,3,4,2,5,6)]

PX$Residue.Both <- str_replace_all(PX$Residue.Both, fixed("|"), ";" )



#KSEA.Complete(KSData, PX, NetworKIN=FALSE, NetworKIN.cutoff=5, m.cutoff=5, p.cutoff=0.05) ##this gives exact result as the website

#or#

KSData.dataset <- KSEA.KS_table(KSData, PX, NetworKIN=FALSE, NetworKIN.cutoff=5)

write.table(KSData.dataset, file = paste0("../analysis/KSEA/", desc, "_Kinase_Substrate_Links.csv"), sep = ",", quote = F, 
            row.names = F, col.names = T)


Scores <- KSEA.Scores(KSData, PX, NetworKIN=FALSE, NetworKIN.cutoff=5)

write.table(Scores, file = paste0("../analysis/KSEA/", desc, "_KSEA_Kinase_Scores.csv"), sep = ",", quote = F, 
            row.names = F, col.names = T)




Barplot <- KSEA.Barplot(KSData, PX, NetworKIN=FALSE, NetworKIN.cutoff=5, m.cutoff=5, p.cutoff=0.05, export=FALSE)


tiff(filename = paste0("../analysis/KSEA/", desc,"_barplot.tiff"),
     width = 5 * 300, 
     height = 5 * 300,
     res = 300,
     compression = "lzw")

KSEA.Barplot(KSData, PX, NetworKIN=FALSE, NetworKIN.cutoff=5, m.cutoff=5, p.cutoff=0.05, export=FALSE)

dev.off()





