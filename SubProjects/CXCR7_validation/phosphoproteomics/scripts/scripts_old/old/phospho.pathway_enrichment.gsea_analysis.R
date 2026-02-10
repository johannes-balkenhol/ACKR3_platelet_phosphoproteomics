################################################################
##### Export for GSEA
##### online guide https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
##### AND
##### https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html


################################################################
#### load packages
BiocManager::install("clusterProfiler", version = "3.8")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
BiocManager::install("ggnewscale")
BiocManager::install("ggridges")
BiocManager::install("europepmc")


library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggnewscale)
library(ggridges)
library(europepmc)



################################################################
##### prepare different table types for GSEA
##### 1. top tables with log2foldchanges of psites
##### 2. collapsed psites to uniprot ID's with collpased log2folchanges


##### 1. get top tables data and abundance data
top.rownames <- c(rownames(top.10[1:20,]),rownames(top.600[1:20,]),rownames(top.1800[1:20,]))
top.norm_intensity <- norm_intensity[top.rownames,]

rownames(top.norm_intensity) <- sapply(strsplit(rownames(top.norm_intensity), ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})




##### 2. prepare collpased psite tables for gsea
df_pl = top.1800
#df_pl = df_pl[, c(7,7,2,1,5)]
df_pl = df_pl[, c(7,7,2,1,5)]
#uniprot, gene_name, avrexpr, log2fold, pvlaue
#f_pl$uniprot.1 = rownames(df_pl)
colnames(df_pl)<- c("Uniprot_ID", "Name", "Average", "Log2FoldChange", "PValue")
df_pl$Name <- sapply(strsplit(rownames(df_pl), ";"), "[[", 2)
#colnames(df_pl)<- c("Uniprot_ID","Average", "Log2FoldChange", "PValue")


## collpase phosphosites to the protein
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



##### annotation 
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

#### prepare input
original_gene_list <- top.collapse$Log2FoldChange
#names(original_gene_list) <- top.collapse$Uniprot_ID
names(original_gene_list) <- top.collapse$Name
gene_list<-na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)



################################################################
#### Gene Set Enrichment
keytypes(org.Hs.eg.db)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.1, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")
			 
			 
#### require(DOSE)
## we use ggplot2 to add x axis labels (ex: ridgeplot)
require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)


#### enrichment map
gse_pw <- pairwise_termsim(gse)
emapplot(gse_pw, showCategory = 40, color = "p.adjust", layout = "nicely")
emapplot(gse_pw, color = "p.adjust", layout = "kk")


#compare_cluster_KEGG <- compareCluster(geneClusters = gene_list, 
#                                       fun = "enrichKEGG",
#                                       organism = "hsa",
#                                       pAdjustMethod = "BH",
#                                       universe = gene_list.top,
#                                       qvalueCutoff = 0.05)
#d <- GOSemSim::godata("org.Hs.eg.db", ont = "B")
#compare_cluster_GO_emap <- enrichplot::pairwise_termsim(compare_cluster_GO, semData = #d)
#emapplot(compare_cluster_GO_emap)


#### Category Netplot
#### Cnetplot
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 10)

#### Ridgeplot
ridgeplot(gse) + labs(x = "enrichment distribution")

#### GSEA plot
# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gse$Description[1:20]
gse$ID[1]
gse@geneSets[[gse$ID[1]]]


top.10[which(top.10$symbol %in% gse@geneSets[[gse$ID[1]]]),]
##https://www.uniprot.org/uniprotkb/O43768/entry ENSA S67 (specified) -| PP2A 
top.600[which(top.600$symbol %in% gse@geneSets[[gse$ID[1]]]),]
top.10[grep("^PPP*", top.10$symbol),]
top.10[grep("^2A5B*", top.10$symbol),]
top.10[grep("O43768", top.10$uniprot_id),]
raw_abundance2[grep("PP2A*", rownames(raw_abundance2)),]
raw_abundance2[grep("P62714*", rownames(raw_abundance2)),]
raw_abundance2[grep("2A5B*", rownames(raw_abundance2)),]

raw_abundance2[grep("^O43768", rownames(raw_abundance2)),]
raw_abundance2[grep("CXCR*", rownames(raw_abundance2)),]
test <- SummarizedExperiment::assay(ppe_filtered,"Quantification")
test[grep("CXCR*", rownames(test)),]


gseaplot(gse, by = "all", title = gse$Description[3], geneSetID = 3)

#### PubMed trend of enriched terms
terms <- gse$Description[1:3]
pmcplot(terms, 2010:2018, proportion=FALSE)


################################################################
#### KEGG Gene Set Enrichment Analysis


#### prepare input
#### prepare input
original_gene_list <- top.collapse$Log2FoldChange
#names(original_gene_list) <- top.collapse$Uniprot_ID
names(original_gene_list) <- top.collapse$Name
gene_list<-na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)




ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)

dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

df2 = top.collapse[top.collapse$Name %in% dedup_ids$SYMBOL,]
df2$Y = dedup_ids$ENTREZID
kegg_gene_list <- df2$Log2FoldChange
names(kegg_gene_list) <- df2$Y
kegg_gene_list<-na.omit(kegg_gene_list)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)


kegg_organism = "hsa"
## organism list http://www.genome.jp/kegg/catalog/org_list.html

options(clusterProfiler.download.method = "wininet")
## according to https://github.com/YuLab-SMU/clusterProfiler/issues/472
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = "hsa",
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.08,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

#### dotplot
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)


#### emapplot
kk_pw <- pairwise_termsim(kk2)
emapplot(kk_pw, showCategory = 40, color = "p.adjust", layout = "nicely")
emapplot(kk_pw, color = "p.adjust", layout = "kk")

#### Category Netplot
#### Cnetplot
cnetplot(kk2, categorySize="pvalue", foldChange=gene_list, showCategory = 10,node_label = "all",)

#### Ridgeplot
ridgeplot(kk2) + labs(x = "enrichment distribution")

#### GSEA plot
# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
kk2$Description[1:20]
kk2$ID[2]
kk2@geneSets[[kk2$ID[2]]]
top.600[which(top.600$symbol %in% kk2@geneSets[[kk2$ID[2]]]),]

gseaplot(kk2, by = "all", title = kk2$Description[5], geneSetID = 5)


##### Pathview
library(pathview)
kk2$ID
#  "hsa05120" "hsa05168" "hsa05161" "hsa04630" "hsa04936"
hsa <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04936", species = kegg_organism)

hsa <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04936", species = kegg_organism, kegg.native = F)


knitr::include_graphics("hsa05120.pathview.png")
