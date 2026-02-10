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
install.packages('plyr')
BiocManager::install("signatureSearch")
BiocManager::install("ReactomePA")
BiocManager::install("fgsea")




library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggnewscale)
library(ggridges)
library(europepmc)
library(plyr)
library(signatureSearch)
require(DOSE)
library("ReactomePA")
library(fgsea)



##### annotation 
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

################################################################
##### use pl_phospho.7.filter.collpase.psite.R
##### for enrichment choose between different collpasing stragtegies
##### 1. absolute value collpasing per uniprot_Id
##### 2. mean value per uniprot_id
##### 3. site sepcific collpasing (select psite wiht highest absolute logFC over timecourse per uniprot_id)

## for 1 and 2
input = list(top.collapse.10, top.collapse.600, top.collapse.1800, 
top.collapse.10.dmso.vs.0s, top.collapse.600.dmso.vs.0s, top.collapse.1800.dmso.vs.0s,
top.collapse.10.cxcr7.vs.0s, top.collapse.600.cxcr7.vs.0s, top.collapse.1800.cxcr7.vs.0s)

## or

## for 3
input = list(top.filter.10, top.filter.600, top.filter.1800, 
top.filter.10.dmso.vs.0s, top.filter.600.dmso.vs.0s, top.filter.1800.dmso.vs.0s,
top.filter.10.cxcr7.vs.0s, top.filter.600.cxcr7.vs.0s, top.filter.1800.cxcr7.vs.0s)


names_input = c("10", "600", "1800", 
"10.dmso.vs.0s", "600.dmso.vs.0s", "1800.dmso.vs.0s",
"10.cxcr7.vs.0s", "600.cxcr7.vs.0s", "1800.cxcr7.vs.0s")


### order inputs
#  input.2 <- top.collapse.600[ order(row.names(top.collapse.600)), ]

for (i in 1:length(input)) {
  input[[i]] <- input[[i]][order(row.names(input[[i]])),]
}

###### Analysis 
## define dataframe with logFC's and meanlogFC's
Tc <- as.data.frame(cbind(input[[1]]$logFC, input[[2]]$logFC, input[[3]]$logFC,
                          input[[4]]$logFC, input[[5]]$logFC, input[[6]]$logFC,
                          input[[7]]$logFC, input[[8]]$logFC, input[[9]]$logFC))
rownames(Tc) <- rownames(input[[1]])
colnames(Tc) <- names_input

#Tc <- as.data.frame(cbind(input[[1]]$meanlogFC, input[[2]]$meanlogFC, input[[3]]$meanlogFC))
#rownames(Tc) <- rownames(input[[1]])
#colnames(Tc) <- c("10s", "600s", "1800s")
ID = paste(paste(sapply(strsplit(rownames(Tc), ";"), "[[", 1)),paste(sapply(strsplit(rownames(Tc), ";"), "[[", 2)),sep=";")
Tc.gene <- phosCollapse(Tc, id=ID, 
                        stat=apply(abs(Tc), 1, max), by = "max")


#############################################################
## iteration
## select which dataset to choose
i = 1


################################################################
##### GSEA in R
#logFC <- top.collapse$Log2FoldChange
#uniprot <- top.collapse$Uniprot_ID
#name <-  top.collapse$Name

logFC <- Tc.gene[,i]
uniprot <- sapply(strsplit(rownames(Tc.gene), ";"), "[[", 1)
name <- sapply(strsplit(rownames(Tc.gene), ";"), "[[", 2)

df <- data.frame(uniprot, name, logFC)

################################################################
#### prepare input with Uniprot for GO
original_gene_list <- logFC
names(original_gene_list) <- uniprot
gene_list_uniprot<-na.omit(original_gene_list)
## sort the list in decreasing order (required for clusterProfiler)
gene_list_uniprot = sort(gene_list_uniprot, decreasing = TRUE)

################################################################
#### prepare input with Unirprot 2 Symbol for GO
original_gene_list2 <- logFC
names(original_gene_list2) <- name
gene_list_name<-na.omit(original_gene_list2)
## sort the list in decreasing order (required for clusterProfiler)
gene_list_name = sort(gene_list_name, decreasing = TRUE)

################################################################
#### prepare input with Uniprot 2 Entrez Gene ID for Reactome and GO
ids<-bitr(names(original_gene_list), fromType = "UNIPROT", toType = "ENTREZID", OrgDb=organism)

dedup_ids = ids[!duplicated(ids[c("UNIPROT")]),]
#dedup_ids = ids[!duplicated(ids[c("UNIPROT")]),]
#dedup_ids = ids

df2 = df[df$uniprot %in% dedup_ids$UNIPROT,]
#df2$Y = dedup_ids$ENTREZID
colnames(dedup_ids) = c("uniprot", "ENTREZID")
df3 <- as.data.frame(merge(df2, dedup_ids, by="uniprot"))
gene_list_entrez <- df3$logFC
names(gene_list_entrez) <- df3$ENTREZID
#kegg_gene_list<-na.omit(kegg_gene_list)
gene_list_entrez = sort(gene_list_entrez, decreasing = TRUE)



################################################################
################################################################
#### GSEA GO
#### input: gene_list_name or gene_list (unirprot)

#keytypes(org.Hs.eg.db)

gse.go <- gseGO(geneList=gene_list_uniprot, 
             ont ="BP", 
             keyType = "UNIPROT", 
             nPerm = 100000, 
             minGSSize = 10, 
             maxGSSize = 150, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")
			 
res.go <- summary(gse.go)
head(res.go, n=20)


res.go.10 <- res.go
gse.go.10 <- gse.go
res.go.600 <- res.go
gse.go.600 <- gse.go
res.go.1800 <- res.go
gse.go.1800 <- gse.go
res.go.10.dmso.vs.0s <- res.go
gse.go.10.dmso.vs.0s <- gse.go


res.go$symbol <- rep(0, dim(res.go)[1])
for (i in 1:dim(res.go)[1]) {
    core_enrich <- res.go[i,]$core_enrichment
    core_list <- strsplit(core_enrich,"/")
    id_pw<-bitr(core_list[[1]], fromType = "UNIPROT", toType = "SYMBOL", OrgDb=organism)[2]
    res.go[i,]$symbol <- paste(id_pw$SYMBOL, collapse = '/')
}














################################################################
################################################################
#### GSEA Reactome
#### https://rdrr.io/bioc/signatureSearch/man/gseReactome.html
#### https://bioconductor.riken.jp/packages/3.3/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html


gse.r <- gsePathway(
                   geneList=gene_list_entrez,
                   nPerm=100000,
                   organism = "human",
                   minGSSize=10,
                   maxGSSize =100,
                   pvalueCutoff=0.5,
                   pAdjustMethod="none",
                   verbose=TRUE
)


res.r <- summary(gse.r)
head(res.r, n=20)


res.r.10 <- res.r
gse.r.10 <- gse.r
res.r.600 <- res.r
gse.r.600 <- gse.r
res.r.1800 <- res.r
gse.r.1800 <- gse.r
res.r.10.dmso.vs.0s <- res.r
gse.r.10.dmso.vs.0s <- gse.r
res.r.600.dmso.vs.0s <- res.r
gse.r.600.dmso.vs.0s <- gse.r
res.r.1800.dmso.vs.0s <- res.r
gse.r.1800.dmso.vs.0s <- gse.r




## translate pathway core enrichemnt into uniprot IDs to evaluate (e.g. put directly to string)
res.r$symbol <- rep(0, dim(res.r)[1])
for (i in 1:dim(res.r)[1]) {
    core_enrich <- res.r[i,]$core_enrichment
    core_list <- strsplit(core_enrich,"/")
    id_pw<-bitr(core_list[[1]], fromType = "ENTREZID", toType = "SYMBOL", OrgDb=organism)[2]
    res.r[i,]$symbol <- paste(id_pw$SYMBOL, collapse = '/')
}


res.r$uniprot <- rep(0, dim(res.r)[1])
for (i in 1:dim(res.r)[1]) {
    core_enrich <- res.r[i,]$core_enrichment
    core_list <- strsplit(core_enrich,"/")
    id_pw<-bitr(core_list[[1]], fromType = "ENTREZID", toType = "UNIPROT", OrgDb=organism)[2]
    res.r[i,]$uniprot <- paste(id_pw$UNIPROT, collapse = '/')
}






################################################################
################################################################
#### analysis & visualization


## check some pathway
dev.new()
Sphingolipid de novo biosynthesis
gseaplot(gseR, geneSetID = "R-HSA-1660661")
dev.new()
Plasma lipoprotein assembly
gseaplot(gseR, geneSetID = "R-HSA-8963898")

top.600[grep("P02647", rownames(top.600)),]

gseR@geneList[["5007"]]

dev.new()
viewPathway("Plasma lipoprotein assembly", readable=TRUE, foldChange=kegg_gene_list)

## pathways of interest





#### require(DOSE)
## we use ggplot2 to add x axis labels (ex: ridgeplot)
dev.new()
dotplot(gseR, showCategory=30, split=".sign") + facet_grid(.~.sign)


#### enrichment map
gse_pw <- pairwise_termsim(gseR)
dev.new()
#png(file=paste0("../analysis/gsea_analysis/", "gsea_emaplot_cdc42",".png"), width = 10, height = 10, units = 'in', res = 600)
emapplot(gse_pw,
showCategory = 10,
color = "pvalue",
layout = "nicely",
node_scale = 0.2,
node_label_size = 0.8,
line_scale = 0.2,
min_edge = 0.2,
cex_label_category = 0.8,
cex_line = 0.2,
cex_category = 0.8)
dev.off()




#emapplot(gse_pw, color = "p.adjust", layout = "kk")


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
dev.new()
cnetplot(gseR, categorySize="pvalue", foldChange=gene_list, showCategory = 5)

#### Ridgeplot
dev.new()
#png(file=paste0("../analysis/gsea_analysis/", "gse.r.10_ridgeplot_",".png"), width = 7, height = 10, units = 'in', res = 600)
ridgeplot(ges.go.10.bp,  showCategory = 20,core_enrichment = TRUE,orderBy = "pvalue", label_format = 40, decreasing = TRUE)
+ labs(x = "enrichment distribution")
# + labs(x = "enrichment distribution") + theme(text=element_text(size=16,  family="Comic Sans MS"))
#label <- word_wrap(label, 45, linesep = NULL)
dev.off()


#res.go.1800.bp <- res.go
#ges.go.1800.bp <- gse.go