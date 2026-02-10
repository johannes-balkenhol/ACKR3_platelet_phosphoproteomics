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

setwd("F:/Masterarbeit_BioWis/Proteomics/pipeline_online/scripts/R")
#setwd("D:/Eigene Datein/Proteomics/MK_PL_Rhoa_Cdc42_KO/scripts/R/")

################################################################
##### prepare different table types for GSEA
##### collapsed psites to uniprot ID's with collpased log2folchanges



##### annotation 
organism = "org.Mm.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)


##### prepare collpased psite tables for gsea
df_pl = top_mk_rhoa_norm_name
df_pl = df_pl[, c(1,8,3,2,6)]
colnames(df_pl)<- c("Uniprot_ID", "Name", "Average", "Log2FoldChange", "PValue")



################################################################
################################################################
#### GSEA GO
#### GO Gene Set Enrichment Analysis



#### prepare input
top.collapse <- df_pl

original_gene_list <- top.collapse$Log2FoldChange
#names(original_gene_list) <- top.collapse$Uniprot_ID
names(original_gene_list) <- top.collapse$Name
gene_list<-na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)





### enrichment
keytypes(org.Mm.eg.db)

gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 5, 
             maxGSSize = 80, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "BH")


#### require(DOSE)
## we use ggplot2 to add x axis labels (ex: ridgeplot)
dev.new()
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
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 5)

#### Ridgeplot
ridgeplot(gse) + labs(x = "enrichment distribution")

#### GSEA plot
# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gse$Description[1:40]
gse$ID[27]
gse@geneSets[[gse$ID[27]]]


top_mk_rhoa_norm_name[which(top_mk_rhoa_norm_name$external_gene_name %in% gse@geneSets[[gse$ID[1]]]),]
##https://www.uniprot.org/uniprotkb/O43768/entry ENSA S67 (specified) -| PP2A 
top.600[which(top.600$symbol %in% gse@geneSets[[gse$ID[1]]]),]
top.10[grep("^PPP*", top.10$symbol),]
top.10[grep("^2A5B*", top.10$symbol),]
top.10[grep("P62714*", top.10$symbol),]
raw_abundance2[grep("PP2A*", rownames(raw_abundance2)),]
raw_abundance2[grep("P62714*", rownames(raw_abundance2)),]
raw_abundance2[grep("2A5B*", rownames(raw_abundance2)),]

raw_abundance2[grep("^O43768", rownames(raw_abundance2)),]
raw_abundance2[grep("CXCR*", rownames(raw_abundance2)),]
test <- SummarizedExperiment::assay(ppe_filtered,"Quantification")
test[grep("CXCR*", rownames(test)),]


gseaplot(gse, by = "all", title = gse$Description[27], geneSetID = 27)

#### PubMed trend of enriched terms
terms <- gse$Description[1:3]
pmcplot(terms, 2010:2018, proportion=FALSE)


################################################################
################################################################
#### GSEA KEGG
#### KEGG Gene Set Enrichment Analysis


##### prepare collpased psite tables for gsea
df_pl = top_mk_rhoa_norm_name
df_pl = df_pl[, c(1,8,3,2,6)]
colnames(df_pl)<- c("Uniprot_ID", "Name", "Average", "Log2FoldChange", "PValue")



################################################################
#### prepare input
top.collapse <- df_pl
#### prepare input
#### prepare input
original_gene_list <- top.collapse$Log2FoldChange
names(original_gene_list) <- top.collapse$Uniprot_ID
#names(original_gene_list) <- top.collapse$Name
gene_list<-na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)




ids<-bitr(names(original_gene_list), fromType = "UNIPROT", toType = "ENTREZID", OrgDb=organism)

dedup_ids = ids[!duplicated(ids[c("UNIPROT")]),]
#dedup_ids = ids[!duplicated(ids[c("UNIPROT")]),]
#dedup_ids = ids

df2 = top.collapse[top.collapse$Uniprot_ID %in% dedup_ids$UNIPROT,]
#df2$Y = dedup_ids$ENTREZID
colnames(dedup_ids) = c("Uniprot_ID", "ENTREZID")
df3 <- as.data.frame(merge(df2, dedup_ids, by="Uniprot_ID"))
kegg_gene_list <- df3$Log2FoldChange
names(kegg_gene_list) <- df3$ENTREZID
kegg_gene_list<-na.omit(kegg_gene_list)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)


kegg_organism = "mmu"
## organism list http://www.genome.jp/kegg/catalog/org_list.html

options(clusterProfiler.download.method = "wininet")
## according to https://github.com/YuLab-SMU/clusterProfiler/issues/472
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = "mmu",
               nPerm        = 1000,
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
#  "hsa05120" "hsa05168" "hsa05161" "hsa04630" "mmu05171"
hsa <- pathview(gene.data=kegg_gene_list, pathway.id="mmu05171", species = kegg_organism)

hsa <- pathview(gene.data=kegg_gene_list, pathway.id="mmu05171", species = kegg_organism, kegg.native = F)


knitr::include_graphics("hsa05120.pathview.png")



################################################################
################################################################
#### GSEA Reactome
#### https://rdrr.io/bioc/signatureSearch/man/gseReactome.html
#### https://bioconductor.riken.jp/packages/3.3/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html


##### prepare collpased psite tables for gsea
df_pl = top_mk_cdc42_norm_name
df_pl = df_pl[, c(1,8,3,2,6)]
colnames(df_pl)<- c("Uniprot_ID", "Name", "Average", "Log2FoldChange", "PValue")



################################################################
#### prepare input
top.collapse <- df_pl
#### prepare input
#### prepare input
original_gene_list <- top.collapse$Log2FoldChange
names(original_gene_list) <- top.collapse$Uniprot_ID
#names(original_gene_list) <- top.collapse$Name
gene_list<-na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)



ids<-bitr(names(original_gene_list), fromType = "UNIPROT", toType = "ENTREZID", OrgDb=organism)

dedup_ids = ids[!duplicated(ids[c("UNIPROT")]),]
#dedup_ids = ids[!duplicated(ids[c("UNIPROT")]),]
#dedup_ids = ids

df2 = top.collapse[top.collapse$Uniprot_ID %in% dedup_ids$UNIPROT,]
#df2$Y = dedup_ids$ENTREZID
colnames(dedup_ids) = c("Uniprot_ID", "ENTREZID")
df3 <- as.data.frame(merge(df2, dedup_ids, by="Uniprot_ID"))
kegg_gene_list <- df3$Log2FoldChange
names(kegg_gene_list) <- df3$ENTREZID
kegg_gene_list<-na.omit(kegg_gene_list)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)




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

gseR <- gsePathway(
                   geneList=kegg_gene_list,
                   nPerm=50000,
                   organism = "mouse",
                   minGSSize=5,
                   maxGSSize = 100,
                   pvalueCutoff=0.2,
                   pAdjustMethod="BH",
                   verbose=TRUE
)

res <- summary(gseR)
head(res, n=20)



ID                                    Description setSize enrichmentScore       NES       pvalue   p.adjust    qvalues rank
R-MMU-114608   R-MMU-114608                        Platelet degranulation       46       0.5995068  2.173990 6.453278e-05 0.02340972 0.01999764  338
R-MMU-76005     R-MMU-76005   Response to elevated platelet cytosolic Ca2+      49       0.5885149  2.159829 6.585013e-05 0.02340972 0.01999764  338
R-MMU-418594   R-MMU-418594                  G alpha (i) signalling events      28       0.6328696  2.045862 1.190051e-04 0.02820421 0.02409332  372
R-MMU-437239   R-MMU-437239                        Recycling pathway of L1      23       0.7175308  2.211925 1.737821e-04 0.03088976 0.02638744  604
R-MMU-388396   R-MMU-388396                     GPCR downstream signalling      58       0.5009646  1.906426 3.433594e-04 0.04120306 0.03519752  476
R-MMU-372790   R-MMU-372790                              Signaling by GPCR      60       0.4908004  1.879005 3.477051e-04 0.04120306 0.03519752  476
R-MMU-157858   R-MMU-157858        Gap junction trafficking and regulation      15       0.7420719  2.052846 4.918570e-04 0.04995862 0.04267692  590
dev.new()
gseaplot(gseR, geneSetID = "R-MMU-76005")
dev.new()
gseaplot(gseR, geneSetID = "R-MMU-388396")

ID                                                                  Description setSize enrichmentScore       NES       pvalue   p.adjust
R-MMU-75153     R-MMU-75153                                                    Apoptotic execution phase      29       0.7830871  2.303737 4.598124e-05 0.01247018
R-MMU-5357801 R-MMU-5357801                                                        Programmed Cell Death      60       0.6185414  2.110986 4.616805e-05 0.01247018
R-MMU-1799339 R-MMU-1799339                  SRP-dependent cotranslational protein targeting to membrane      74      -0.5565153 -1.938716 7.011639e-05 0.01247018
R-MMU-975956   R-MMU-975956 Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC)      78      -0.5448040 -1.914976 7.015575e-05 0.01247018
R-MMU-109581   R-MMU-109581                                                                    Apoptosis      47       0.6499275  2.115191 9.184001e-05 0.01305965
R-MMU-111465   R-MMU-111465                                      Apoptotic cleavage of cellular proteins      18       0.7713928  2.037866 2.747001e-04 0.03255196
R-MMU-2559584 R-MMU-2559584               Formation of Senescence-Associated Heterochromatin Foci (SAHF)       7       0.8939629  1.875140 7.253933e-04 0.06621638
R-MMU-5368287 R-MMU-5368287                                                    Mitochondrial translation      24      -0.6747432 -1.873275 7.450507e-04 0.06621638

dev.new()
gseaplot(gseR, geneSetID = "R-MMU-75153")
dev.new()
gseaplot(gseR, geneSetID = "R-MMU-1799339")


dev.new()
viewPathway("G alpha (i) signalling events", readable=TRUE, foldChange=kegg_gene_list)


#### require(DOSE)
## we use ggplot2 to add x axis labels (ex: ridgeplot)
dev.new()
dotplot(gseR, showCategory=30, split=".sign") + facet_grid(.~.sign)


#### enrichment map
gse_pw <- pairwise_termsim(gseR)
dev.new()
png(file=paste0("../../analysis/gsea_analysis/", "gsea_emaplot_cdc42",".png"), width = 10, height = 10, units = 'in', res = 600)
emapplot(gse_pw,
showCategory = 50,
color = "p.adjust",
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
png(file=paste0("../../analysis/gsea_analysis/", "gsea_ridgeplot_cdc42",".png"), width = 7, height = 8, units = 'in', res = 600)
ridgeplot(gseR,  showCategory = 20,core_enrichment = TRUE,orderBy = "p.adjust", label_format = 30, decreasing = TRUE)
# + labs(x = "enrichment distribution") + theme(text=element_text(size=16,  family="Comic Sans MS"))
dev.off()
