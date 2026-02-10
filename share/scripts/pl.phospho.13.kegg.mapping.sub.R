library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggnewscale)
library(ggridges)
library(europepmc)
library(plyr)
library(dplyr)
library(signatureSearch)
require(DOSE)
library("ReactomePA")
library(fgsea)
library("pathview")
library(PhosR)
library(RColorBrewer)


####extra plots and specific terms#######

################################################################
#### prepare input with Uniprot for GO
input = list(top.collapse.10, top.collapse.600, top.collapse.1800, 
             top.collapse.10.dmso.vs.0s, top.collapse.600.dmso.vs.0s, top.collapse.1800.dmso.vs.0s,
             top.collapse.10.cxcr7.vs.0s, top.collapse.600.cxcr7.vs.0s, top.collapse.1800.cxcr7.vs.0s)

## or

## for 3
# input = list(top.filter.10, top.filter.600, top.filter.1800, 
#              top.filter.10.dmso.vs.0s, top.filter.600.dmso.vs.0s, top.filter.1800.dmso.vs.0s,
#              top.filter.10.cxcr7.vs.0s, top.filter.600.cxcr7.vs.0s, top.filter.1800.cxcr7.vs.0s)

names_input = c("10", "600", "1800", 
                "10.dmso.vs.0s", "600.dmso.vs.0s", "1800.dmso.vs.0s",
                "10.cxcr7.vs.0s", "600.cxcr7.vs.0s", "1800.cxcr7.vs.0s")


### order inputs
#  input.2 <- top.collapse.600[ order(row.names(top.collapse.600)), ]

for (i in 1:length(input)) {
  input[[i]] <- input[[i]][input[[i]]$PValue<=0.05,]
  input[[i]] <- input[[i]][order(row.names(input[[i]])),]
  input[[i]]$rownames <- rownames(input[[i]])
  colnames(input[[i]])[2] <- paste0("logfc.", names_input[i])
  
}

###### Analysis 
## define dataframe with logFC's and meanlogFC's
Tcsign <- as.data.frame(merge(input[[1]][c(2,6)], input[[2]][c(2,6)],  by ="rownames", all = TRUE))
Tcsign <- as.data.frame(merge(Tcsign, input[[3]][c(2,6)],  by ="rownames", all = TRUE))
#Tcsign <- as.data.frame(merge(Tcsign, input[[4]][c(2,6)],  by ="rownames", all = TRUE))
#Tcsign <- as.data.frame(merge(Tcsign, input[[5]][c(2,6)],  by ="rownames", all = TRUE))
#Tcsign <- as.data.frame(merge(Tcsign, input[[6]][c(2,6)],  by ="rownames", all = TRUE))
#Tcsign <- as.data.frame(merge(Tcsign, input[[7]][c(2,6)],  by ="rownames", all = TRUE))



rownames(Tcsign) <- Tcsign$rownames
Tcsign$rownames <- NULL


#Tcsign <- as.data.frame(cbind(input[[1]]$meanlogFC, input[[2]]$meanlogFC, input[[3]]$meanlogFC))
#rownames(Tcsign) <- rownames(input[[1]])
#colnames(Tcsign) <- c("10s", "600s", "1800s")
ID = paste(paste(sapply(strsplit(rownames(Tcsign), ";"), "[[", 1)),paste(sapply(strsplit(rownames(Tcsign), ";"), "[[", 2)),sep=";")
Tcsign.gene <- phosCollapse(Tcsign, id=ID, 
                            stat=apply(abs(Tcsign), 1, max), by = "max")

Tcsign.gene[is.na(Tcsign.gene)] <- 0

logFC <- Tcsign.gene[,1:3]
uniprot <- sapply(strsplit(rownames(Tcsign.gene), ";"), "[[", 1)
name <- sapply(strsplit(rownames(Tcsign.gene), ";"), "[[", 2)

df <- data.frame(uniprot, name, logFC)
colnames(df)<-c("uniprot","name", "t10","t600","t1800")




original_gene_list <- df
row.names(original_gene_list) <- df$uniprot
gene_list_uniprot<-na.omit(original_gene_list)
colnames(gene_list_uniprot)<-c("uniprot","name","t10","t600","t1800")

## sort the list in decreasing order (required for clusterProfiler)
#gene_list_uniprot = sort(gene_list_uniprot, decreasing = TRUE)

################################################################
#### prepare input with Unirprot 2 Symbol for GO
original_gene_list2 <- logFC
row.names(original_gene_list2) <- name
gene_list_name<-na.omit(original_gene_list2)
## sort the list in decreasing order (required for clusterProfiler)
colnames(gene_list_name)<-c("t10","t600","t1800")


################################################################
#### prepare input with Uniprot 2 Entrez Gene ID for all timepoints together
ids<-bitr(row.names(original_gene_list), fromType = "UNIPROT", toType = "ENTREZID", OrgDb=organism)

dedup_ids = ids[!duplicated(ids[c("UNIPROT")]),]
#dedup_ids = ids[!duplicated(ids[c("UNIPROT")]),]
#dedup_ids = ids

df2 = df[df$uniprot %in% dedup_ids$UNIPROT,]
#df2$Y = dedup_ids$ENTREZID
colnames(dedup_ids) = c("uniprot", "ENTREZID")
df3 <- as.data.frame(merge(df2, dedup_ids, by="uniprot"))
gene_list_entrez <- df3[,3:5]
row.names(gene_list_entrez) <- df3$ENTREZID
#kegg_gene_list<-na.omit(kegg_gene_list)
#gene_list_entrez = sort(gene_list_entrez, decreasing = TRUE)


df3[df3 == 0] <- NA


setwd("../analysis/KEGG Mapping/")

#platelet activation
hsa04611 <- pathview(gene.data  = gene_list_entrez[,c(1:3)],
                     pathway.id = "hsa04611",
                     species    = "hsa",
                     limit      = list(gene=max(abs(gene_list_entrez))),
                     gene.idtype="entrez",
                     low = list(gene = "navy", cpd = "blue"),
                     mid = list(gene = "#E0E0E0", cpd = "gray"),
                     high = list(gene = "red", cpd = "yellow"))

#pi3k-akt

hsa04151 <- pathview(gene.data  = gene_list_entrez[,c(1:3)],
                     pathway.id = "hsa04151",
                     species    = "hsa",
                     limit      = list(gene=max(abs(gene_list_entrez))),
                     gene.idtype="entrez",
                     low = list(gene = "navy", cpd = "blue"),
                     mid = list(gene = "#E0E0E0", cpd = "gray"),
                     high = list(gene = "red", cpd = "yellow"))


#mtor
hsa04150 <- pathview(gene.data  = gene_list_entrez[,c(1:3)],
                     pathway.id = "hsa04150",
                     species    = "hsa",
                     limit      = list(gene=max(abs(gene_list_entrez))),
                     gene.idtype="entrez",
                     low = list(gene = "navy", cpd = "blue"),
                     mid = list(gene = "#E0E0E0", cpd = "gray"),
                     high = list(gene = "red", cpd = "yellow"))

#calcium signaling
hsa04020 <- pathview(gene.data  = gene_list_entrez[,c(1:3)],
                     pathway.id = "hsa04020",
                     species    = "hsa",
                     limit      = list(gene=max(abs(gene_list_entrez))),
                     gene.idtype="entrez",
                     low = list(gene = "navy", cpd = "blue"),
                     mid = list(gene = "#E0E0E0", cpd = "gray"),
                     high = list(gene = "red", cpd = "yellow"))


#Fluid shear stress and atherosclerosis
hsa05418 <- pathview(gene.data  = gene_list_entrez[,c(1:3)],
                     pathway.id = "hsa05418",
                     species    = "hsa",
                     limit      = list(gene=max(abs(gene_list_entrez))),
                     gene.idtype="entrez",
                     low = list(gene = "navy", cpd = "blue"),
                     mid = list(gene = "#E0E0E0", cpd = "gray"),
                     high = list(gene = "red", cpd = "yellow"))

#MAPK signaling pathway
hsa04010 <- pathview(gene.data  = gene_list_entrez[,c(1:3)],
                     pathway.id = "hsa04010",
                     species    = "hsa",
                     limit      = list(gene=max(abs(gene_list_entrez))),
                     gene.idtype="entrez",
                     low = list(gene = "navy", cpd = "blue"),
                     mid = list(gene = "#E0E0E0", cpd = "gray"),
                     high = list(gene = "red", cpd = "yellow"))

#Regulation of actin cytoskeleton
hsa04810 <- pathview(gene.data  = gene_list_entrez[,c(1:3)],
                     pathway.id = "hsa04810",
                     species    = "hsa",
                     limit      = list(gene=max(abs(gene_list_entrez))),
                     gene.idtype="entrez",
                     low = list(gene = "navy", cpd = "blue"),
                     mid = list(gene = "#E0E0E0", cpd = "gray"),
                     high = list(gene = "red", cpd = "yellow"))

#AMPK signaling pathway
hsa04152 <- pathview(gene.data  = gene_list_entrez[,c(1:3)],
                     pathway.id = "hsa04152",
                     species    = "hsa",
                     limit      = list(gene=max(abs(gene_list_entrez))),
                     gene.idtype="entrez",
                     low = list(gene = "navy", cpd = "blue"),
                     mid = list(gene = "#E0E0E0", cpd = "gray"),
                     high = list(gene = "red", cpd = "yellow"))

#Apoptosis
hsa04210 <- pathview(gene.data  = gene_list_entrez[,c(1:3)],
                     pathway.id = "hsa04210",
                     species    = "hsa",
                     limit      = list(gene=max(abs(gene_list_entrez))),
                     gene.idtype="entrez",
                     low = list(gene = "navy", cpd = "blue"),
                     mid = list(gene = "#E0E0E0", cpd = "gray"),
                     high = list(gene = "red", cpd = "yellow"))

#Complement and coagulation cascade
hsa04610 <- pathview(gene.data  = gene_list_entrez[,c(1:3)],
                     pathway.id = "hsa04610",
                     species    = "hsa",
                     limit      = list(gene=max(abs(gene_list_entrez))),
                     gene.idtype="entrez",
                     low = list(gene = "navy", cpd = "blue"),
                     mid = list(gene = "#E0E0E0", cpd = "gray"),
                     high = list(gene = "red", cpd = "yellow"))

# Arachidonic acid metabolism
hsa00590 <- pathview(gene.data  = gene_list_entrez[,c(1:3)],
                     pathway.id = "hsa00590",
                     species    = "hsa",
                     limit      = list(gene=max(abs(gene_list_entrez))),
                     gene.idtype="entrez",
                     low = list(gene = "navy", cpd = "blue"),
                     mid = list(gene = "#E0E0E0", cpd = "gray"),
                     high = list(gene = "red", cpd = "yellow"))

# cAMP signaling
hsa04024 <- pathview(gene.data  = gene_list_entrez[,c(1:3)],
                     pathway.id = "hsa04024",
                     species    = "hsa",
                     limit      = list(gene=max(abs(gene_list_entrez))),
                     gene.idtype="entrez",
                     low = list(gene = "navy", cpd = "blue"),
                     mid = list(gene = "#E0E0E0", cpd = "gray"),
                     high = list(gene = "red", cpd = "yellow"))


