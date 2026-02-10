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
library(signatureSearch)
require(DOSE)
library("ReactomePA")
library(fgsea)



##### annotation 
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

################################################################
##### prepare different table types for GSEA
##### 1. top tables with log2foldchanges of psites
##### 2. collapsed psites to uniprot ID's with collpased log2folchanges


##### prepare collpased psite tables for gsea
## choose here the table of desire !!!!
#df_pl = top.600
df_pl = top.600.cxcr7.vs.0s
## get rid of site without annotation
df_pl = df_pl[-which(df_pl$psite == ""),]
df_pl = df_pl[, c(7,8,2,1,5,9)]

colnames(df_pl)<- c("Uniprot_ID", "Name", "Average", "Log2FoldChange", "PValue", "PSite")

## collpase same site by mean (sequence sometime slightly different)
df_pl<- ddply(df_pl, .(Uniprot_ID, PSite), summarise,
              Name = max(Name),
              Average = mean(Average),
              Log2FoldChange = mean(Log2FoldChange),
              PValue = mean(PValue))


## collpase phosphosites to the protein
top.collapse <- ddply(df_pl, .(Uniprot_ID, Name), summarise,
              Average = max(Average),
              minLog2FoldChange = min(Log2FoldChange),
              maxLog2FoldChange = max(Log2FoldChange),
              AbsLog2FoldChange = max(abs(Log2FoldChange)),
              MeanLog2FoldChange = mean(Log2FoldChange),
              PValue = min(PValue))


# get the optima (maxima extend; maxima and minima)
top.collapse$Log2FoldChange <- 0
top.collapse$Log2FoldChange[abs(top.collapse[, "minLog2FoldChange"]) < abs(top.collapse[, "maxLog2FoldChange"])] = top.collapse$maxLog2FoldChange[abs(top.collapse[, "minLog2FoldChange"]) < abs(top.collapse[, "maxLog2FoldChange"])]

top.collapse$Log2FoldChange[abs(top.collapse[, "minLog2FoldChange"]) > abs(top.collapse[, "maxLog2FoldChange"])] = top.collapse$minLog2FoldChange[abs(top.collapse[, "minLog2FoldChange"]) > abs(top.collapse[, "maxLog2FoldChange"])]

top.collapse$Log2FoldChange[abs(top.collapse[, "minLog2FoldChange"]) == abs(top.collapse[, "maxLog2FoldChange"])] = top.collapse$minLog2FoldChange[abs(top.collapse[, "minLog2FoldChange"]) == abs(top.collapse[, "maxLog2FoldChange"])]

# or instead get the mean of all psites per protein mean (performs good)
top.collapse$Log2FoldChange <- top.collapse[, "MeanLog2FoldChange"]

################################################################
##### GSEA in R

################################################################
#### prepare input with Uniprot for GO

original_gene_list <- top.collapse$Log2FoldChange
names(original_gene_list) <- top.collapse$Uniprot_ID
gene_list<-na.omit(original_gene_list)
## sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)



################################################################
#### prepare input with Unirprot 2 Symbol for GO

original_gene_list2 <- top.collapse$Log2FoldChange
names(original_gene_list2) <- top.collapse$Name
gene_list_name<-na.omit(original_gene_list2)
## sort the list in decreasing order (required for clusterProfiler)
gene_list_name = sort(gene_list_name, decreasing = TRUE)


################################################################
#### prepare input with Uniprot 2 Entrez Gene ID for Reactome and GO

ids<-bitr(names(original_gene_list), fromType = "UNIPROT", toType = "ENTREZID", OrgDb=organism)

dedup_ids = ids[!duplicated(ids[c("UNIPROT")]),]
#dedup_ids = ids[!duplicated(ids[c("UNIPROT")]),]
#dedup_ids = ids

df2 = top.collapse[top.collapse$Uniprot_ID %in% dedup_ids$UNIPROT,]
#df2$Y = dedup_ids$ENTREZID
colnames(dedup_ids) = c("Uniprot_ID", "ENTREZID")
df3 <- as.data.frame(merge(df2, dedup_ids, by="Uniprot_ID"))
entrez_gene_list <- df3$Log2FoldChange
names(entrez_gene_list) <- df3$ENTREZID
#kegg_gene_list<-na.omit(kegg_gene_list)
entrez_gene_list = sort(entrez_gene_list, decreasing = TRUE)


################################################################
################################################################
#### GSEA GO
#### input: gene_list_name or gene_list (unirprot)

#keytypes(org.Hs.eg.db)

gse.go <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "UNIPROT", 
             nPerm = 100000, 
             minGSSize = 5, 
             maxGSSize = 300, 
             pvalueCutoff = 0.1, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "BH")
			 
res.go <- summary(gse.go)
head(res.go, n=20)


res.go.10.dmso.vs.0s <- res.go
gse.go.10.dmso.vs.0s <- gse.go
res.go.10.cxcr7.vs.0s <- res.go
gse.go.10.cxcr7.vs.0s <- gse.go
res.go.600.dmso.vs.0s <- res.go
gse.go.600.dmso.vs.0s <- gse.go


################################################################
################################################################
#### GSEA Reactome
#### https://rdrr.io/bioc/signatureSearch/man/gseReactome.html
#### https://bioconductor.riken.jp/packages/3.3/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html

#gseR <- gseReactome(
#              geneList=kegg_gene_list,
#              organism = "mouse",
#              exponent = 1,
#              nPerm = 1000,
#              minGSSize = 10,
#              maxGSSize = 80,
#              pvalueCutoff = 0.05,
#              pAdjustMethod = "BH",
#              verbose = TRUE,
#              readable = FALSE
#)


gse.r <- gsePathway(
                   geneList=entrez_gene_list,
                   nPerm=10000,
                   organism = "human",
                   minGSSize=5,
                   maxGSSize =100,
                   pvalueCutoff=0.8,
                   pAdjustMethod="BH",
                   verbose=TRUE
)


res.r <- summary(gse.r)
head(res.r, n=20)

res.r.10.dmso.vs.0s <- res.r
gse.r.10.dmso.vs.0s <- gse.r
res.r.10.cxcr7.vs.0s <- res.r
gse.r.10.cxcr7.vs.0s <- gse.r
res.r.600.dmso.vs.0s <- res.r
gse.r.600.dmso.vs.0s <- gse.r




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
## t.10
                         ID                               Description setSize enrichmentScore       NES      pvalue    p.adjust    qvalue
R-HSA-2871837 R-HSA-2871837           FCERI mediated NF-kB activation      15      -0.6882059 -1.846436 0.002195895 0.002195895 0.7510697
R-HSA-111465   R-HSA-111465   Apoptotic cleavage of cellular proteins      16       0.6343063  1.838961 0.003855833 0.003855833 0.7510697
R-HSA-392154   R-HSA-392154 Nitric oxide stimulates guanylate cyclase       5      -0.8653425 -1.705969 0.004297967 0.004297967 0.7510697
R-HSA-418457   R-HSA-418457                              cGMP effects       5      -0.8653425 -1.705969 0.004297967 0.004297967 0.7510697
R-HSA-8863795 R-HSA-8863795         Downregulation of ERBB2 signaling       8      -0.7675403 -1.732064 0.006497954 0.006497954 0.7781515
R-HSA-75153     R-HSA-75153                 Apoptotic execution phase      18       0.5836617  1.748870 0.007181663 0.007181663 0.7781515
R-HSA-192105   R-HSA-192105    Synthesis of bile acids and bile salts       6      -0.8135226 -1.689770 0.007792648 0.007792648 0.7781515
R-HSA-1660661 R-HSA-1660661         Sphingolipid de novo biosynthesis       6      -0.8006537 -1.663040 0.010503134 0.010503134 0.8307331
R-HSA-9013407 R-HSA-9013407                         RHOH GTPase cycle      15      -0.6334185 -1.699443 0.010696134 0.010696134 0.8307331
R-HSA-381038   R-HSA-381038         XBP1(S) activates chaperone genes       9       0.6922236  1.686905 0.015932315 0.015932315 0.8697238
R-HSA-381070   R-HSA-381070            IRE1alpha activates chaperones       9       0.6922236  1.686905 0.015932315 0.015932315 0.8697238
R-HSA-211945   R-HSA-211945  Phase I - Functionalization of compounds       5       0.7991738  1.610206 0.017159239 0.017159239 0.8697238
R-HSA-381119   R-HSA-381119           Unfolded Protein Response (UPR)      10       0.6638867  1.671932 0.018349235 0.018349235 0.8697238
R-HSA-5357956 R-HSA-5357956  TNFR1-induced NFkappaB signaling pathway      10      -0.6828832 -1.640995 0.018676618 0.018676618 0.8697238
R-HSA-194068   R-HSA-194068        Bile acid and bile salt metabolism       7      -0.7432709 -1.610450 0.020912832 0.020912832 0.8697238
R-HSA-4086398 R-HSA-4086398                              Ca2+ pathway       9      -0.6953626 -1.620854 0.021835263 0.021835263 0.8697238
R-HSA-5607764 R-HSA-5607764               CLEC7A (Dectin-1) signaling      18      -0.5703245 -1.605438 0.022043246 0.022043246 0.8697238
R-HSA-5621481 R-HSA-5621481            C-type lectin receptors (CLRs)      25      -0.5177172 -1.581348 0.022396321 0.022396321 0.8697238
R-HSA-264876   R-HSA-264876                        Insulin processing       6       0.7383019  1.575034 0.029564846 0.029564846 0.9328747
R-HSA-9610379 R-HSA-9610379                          HCMV Late Events       5      -0.7863624 -1.550264 0.030371032 0.030371032 0.9328747

## t.600
R-HSA-453274   R-HSA-453274                                                                Mitotic G2-G2/M phases      29      -0.5830109
R-HSA-69275     R-HSA-69275                                                                       G2/M Transition      29      -0.5830109
R-HSA-2565942 R-HSA-2565942                                        Regulation of PLK1 Activity at G2/M Transition      19      -0.6284222
R-HSA-392518   R-HSA-392518                                                                  Signal amplification       5       0.8692722
R-HSA-69278     R-HSA-69278                                                                   Cell Cycle, Mitotic      68      -0.4309030
R-HSA-1660661 R-HSA-1660661                                                     Sphingolipid de novo biosynthesis       6       0.8296704
R-HSA-111933   R-HSA-111933                                                             Calmodulin induced events      10      -0.7120602
R-HSA-111997   R-HSA-111997                                                                           CaM pathway      10      -0.7120602
R-HSA-1169410 R-HSA-1169410                                           Antiviral mechanism by IFN-stimulated genes      11      -0.6881215
R-HSA-5610787 R-HSA-5610787                                                                  Hedgehog 'off' state      16      -0.6107687
R-HSA-168255   R-HSA-168255                                                                   Influenza Infection      10       0.6776430
R-HSA-8936459 R-HSA-8936459 RUNX1 regulates genes involved in megakaryocyte differentiation and platelet function       6       0.7857667
R-HSA-2559583 R-HSA-2559583                                                                   Cellular Senescence      16      -0.5963594
R-HSA-8963898 R-HSA-8963898                                                           Plasma lipoprotein assembly       5      -0.8413184
R-HSA-109581   R-HSA-109581                                                                             Apoptosis      39      -0.4767082
R-HSA-380259   R-HSA-380259                                                  Loss of Nlp from mitotic centrosomes      15      -0.6039361
R-HSA-380270   R-HSA-380270                              Recruitment of mitotic centrosome proteins and complexes      15      -0.6039361
R-HSA-380284   R-HSA-380284 Loss of proteins required for interphase microtubule organization from the centrosome      15      -0.6039361
R-HSA-380287   R-HSA-380287                                                                 Centrosome maturation      15      -0.6039361
R-HSA-8854518 R-HSA-8854518                                                              AURKA Activation by TPX2      15      -0.6039361


## translate pathway core enrichemnt into uniprot IDs to evaluate (e.g. put directly to string)
core_enrich <- res.r.1800["R-HSA-8856828",]$core_enrichment
core_list <- strsplit(core_enrich,"/")
id_pw<-bitr(core_list[[1]], fromType = "ENTREZID", toType = "UNIPROT", OrgDb=organism)

id_pw$UNIPROT


res.go.10.bp["GO:1901222",]$core_enrichment



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