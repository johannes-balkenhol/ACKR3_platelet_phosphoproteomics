################################################################
##### Export for GSEA


##### get intesity data
norm_intensity <- SummarizedExperiment::assay(ppe4, "normalised")
UniprotID<- sapply(strsplit(rownames(ppe), ";"), "[[", 1)
GeneSymbol <- sapply(strsplit(rownames(ppe), ";"), "[[", 2)



##### collpase with collpase funiton
norm_intensity.collapse <- phosCollapse(norm_intensity, id=sapply(strsplit(rownames(norm_intensity), ";"), "[[", 1), 
                        stat=apply(abs(norm_intensity), 1, max), by = "max")

x_tt <- as.factor(grps)

## prepare to tables for GSEA
## 1. intensity table
## 2. sample table



##### prepare table for GO enrichment 
df_pl = top.600
#df_pl = df_pl[, c(7,7,2,1,5)]
df_pl = df_pl[, c(7,7,2,1,5)]
#uniprot, gene_name, avrexpr, log2fold, pvlaue
#f_pl$uniprot.1 = rownames(df_pl)
colnames(df_pl)<- c("Uniprot_ID", "Name", "Average", "Log2FoldChange", "PValue")
df_pl$Name <- sapply(strsplit(rownames(df_pl), ";"), "[[", 2)
#colnames(df_pl)<- c("Uniprot_ID","Average", "Log2FoldChange", "PValue")



##### collpase phosphosites to the protein
top.collapse <- ddply(df_pl, .(Uniprot_ID, Name), summarise,
              Average = max(Average),
              minLog2FoldChange = min(Log2FoldChange),
			  maxLog2FoldChange = max(Log2FoldChange),
              AbsLog2FoldChange = max(abs(Log2FoldChange)),
              PValue = min(PValue))

top.collapse[abs(top.collapse[, "minLog2FoldChange"]) < abs(top.collapse[, "maxLog2FoldChange"]),]

top.collapse$Log2FoldChange <- 0
top.collapse$Log2FoldChange[abs(top.collapse[, "minLog2FoldChange"]) < abs(top.collapse[, "maxLog2FoldChange"])] = top.collapse$maxLog2FoldChange[abs(top.collapse[, "minLog2FoldChange"]) < abs(top.collapse[, "maxLog2FoldChange"])]

top.collapse$Log2FoldChange[abs(top.collapse[, "minLog2FoldChange"]) > abs(top.collapse[, "maxLog2FoldChange"])] = top.collapse$minLog2FoldChange[abs(top.collapse[, "minLog2FoldChange"]) > abs(top.collapse[, "maxLog2FoldChange"])]

top.collapse$Log2FoldChange[abs(top.collapse[, "minLog2FoldChange"]) == abs(top.collapse[, "maxLog2FoldChange"])] = top.collapse$minLog2FoldChange[abs(top.collapse[, "minLog2FoldChange"]) == abs(top.collapse[, "maxLog2FoldChange"])]



################################################################
##### GSEA in R

#### load packages
BiocManager::install("clusterProfiler", version = "3.8")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)


##### get top tables data and abundance data
top.rownames <- c(rownames(top.10[1:20,]),rownames(top.600[1:20,]),rownames(top.1800[1:20,]))
top.norm_intensity <- norm_intensity[top.rownames,]

rownames(top.norm_intensity) <- sapply(strsplit(rownames(top.norm_intensity), ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})






##### annotation 

organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

#### prepare input

# we want the log2 fold change 
original_gene_list <- top.collapse$Log2FoldChange
# name the vector
names(original_gene_list) <- top.collapse$Uniprot_ID
# omit any NA values 
gene_list<-na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)



#### Gene Set Enrichment
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "UNIPROT", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.1, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")
			 
			 
#### require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)


#### enrichment map 
emapplot(gse, showCategory = 10)







