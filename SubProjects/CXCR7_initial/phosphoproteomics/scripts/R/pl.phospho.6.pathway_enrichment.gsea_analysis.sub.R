################################################################
##### Export for GSEA
##### online guide https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
##### AND
##### https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html


################################################################
#### load packages
BiocManager::install("clusterProfiler", version = "3.16")
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
library("pathview")
library(PhosR)



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
input = list(top.collapse.10, top.collapse.30, top.collapse.60, top.collapse.300, top.collapse.600, top.collapse.900, 
             top.collapse.1800)


## or

## for 3
#input = list(top.filter.10, top.filter.30, top.filter.60, top.filter.300,top.filter.600, top.filter.900, top.filter.1800)



names_input = c("10", "30", "60", "300", "600", "900", "1800")

### order inputs
#  input.2 <- top.collapse.600[ order(row.names(top.collapse.600)), ]

for (i in 1:length(input)) {
  input[[i]] <- input[[i]][order(row.names(input[[i]])),]
  
}

###### Analysis 
## define dataframe with logFC's and meanlogFC's
Tc <- as.data.frame(cbind(input[[1]]$logFC, input[[2]]$logFC, input[[3]]$logFC,
                          input[[4]]$logFC, input[[5]]$logFC, input[[6]]$logFC,
                          input[[7]]$logFC))
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
#i = 1

################################################################
##### GSEA in R
#logFC <- top.collapse$Log2FoldChange
#uniprot <- top.collapse$Uniprot_ID
#name <-  top.collapse$Name

#for (i in 1:length(input)) {




##function does all.
new.function <- function(i) {
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
                  #nPerm = 100000, 
                  minGSSize = 10, 
                  maxGSSize = 150, 
                  pvalueCutoff = 1, 
                  verbose = TRUE, 
                  OrgDb = org.Hs.eg.db, 
                  pAdjustMethod = "BH")
  
  gse.go_genename <- setReadable(gse.go, OrgDb = org.Hs.eg.db, keyType="UNIPROT")
  
  assign(paste0("gse.go", names_input[[i]]), gse.go_genename)
  
  write.table(gse.go_genename, file = paste0("../analysis/GSEA/",names_input[[i]], "_GO_all.txt"), sep = "\t", quote = F, 
              row.names = F, col.names = T)
  
  
  #### require(DOSE)
  ## we use ggplot2 to add x axis labels (ex: ridgeplot)
  
  dot<- dotplot(gse.go_genename, showCategory=10, split=".sign") + facet_grid(.~.sign)
  
  tiff(filename = paste0("../analysis/GSEA/", names_input[[i]], "GOBP_dotplot.tiff"),
       width = 12 * 300, 
       height = 12 * 300,
       res = 300,
       compression = "lzw")
  print(dot)
  
  dev.off()
  
  #### enrichment map
  gse_pw <- pairwise_termsim(gse.go_genename)
  emap<- emapplot(gse_pw,
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
  tiff(filename = paste0("../analysis/GSEA/", names_input[[i]], "gosim.tiff"),
       width = 12 * 300, 
       height = 12 * 300,
       res = 300,
       compression = "lzw")
  
  print(emap)
  #png(file=paste0("../analysis/gsea_analysis/", "gsea_emaplot_cdc42",".png"), width = 10, height = 10, units = 'in', res = 600)
  
  dev.off()
  
  
  ##############################
  
  
  # res.go <- summary(gse.go)
  # head(res.go, n=20)
  # 
  # 
  # gse.go@result %>% 
  #   arrange(pvalue) %>% 
  #   head(8)
  # 
  # assign(paste0("res.go", names_input[[i]]), res.go)
  # 
  # res.go$symbol <- rep(0, dim(res.go)[1])
  # for (i in 1:dim(res.go)[1]) {
  #     core_enrich <- res.go[i,]$core_enrichment
  #     core_list <- strsplit(core_enrich,"/")
  #     id_pw<-bitr(core_list[[1]], fromType = "UNIPROT", toType = "SYMBOL", OrgDb=organism)[2]
  #     res.go[i,]$symbol <- paste(id_pw$SYMBOL, collapse = '/')
  # }
  
  
  
  
  
  
  
  
  
  
  ################################################################
  ################################################################
  #### GSEA Reactome
  #### https://rdrr.io/bioc/signatureSearch/man/gseReactome.html
  #### https://bioconductor.riken.jp/packages/3.3/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html
  
  
  
  
  gse.r <- gsePathway(
    geneList=gene_list_entrez,
    #nPerm=100000,
    organism = "human",
    minGSSize=10,
    maxGSSize =100,
    pvalueCutoff=1, #to get the table, filtering can be done after
    pAdjustMethod="BH",
    verbose=TRUE)
  
  gse.r_genename <- setReadable(gse.r, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  
  assign(paste0("gse.r", names_input[[i]]), gse.r_genename, envir = globalenv())
  
  write.table(gse.r_genename, file = paste0("../analysis/GSEA/",names_input[[i]], "_Reactome_all.txt"), sep = "\t", quote = F, 
              row.names = F, col.names = T)
  
  
  
  # res.r <- summary(gse.r)
  # head(res.r, n=20)
  # 
  # assign(paste0("res.r", names_input[[i]]), res.r)
  # assign(paste0("gse.r", names_input[[i]]), gse.r)
  
  
  
  #### require(DOSE)
  ## we use ggplot2 to add x axis labels (ex: ridgeplot)
  dot.r<-dotplot(gse.r_genename, showCategory=10, split=".sign") + facet_grid(.~.sign)
  
  
  tiff(filename = paste0("../analysis/GSEA/", names_input[[i]], "Reactome_dotplot.tiff"),
       width = 12 * 300, 
       height = 12 * 300,
       res = 300,
       compression = "lzw")
  print(dot.r)
  dev.off()
  
  #### enrichment map
  gse_pw.r <- pairwise_termsim(gse.r_genename)
  emap.r<-emapplot(gse_pw.r,
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
  tiff(filename = paste0("../analysis/GSEA/", names_input[[i]], "gosim_r.tiff"),
       width = 12 * 300, 
       height = 12 * 300,
       res = 300,
       compression = "lzw")
  print(emap.r)
  #png(file=paste0("../analysis/gsea_analysis/", "gsea_emaplot_cdc42",".png"), width = 10, height = 10, units = 'in', res = 600)
  
  dev.off()
  
  #### GSEA KEGG###
  gse.k <- gseKEGG(geneList     = gene_list_uniprot,
                   organism     = 'hsa',
                   keyType = "uniprot",
                   minGSSize    = 10,
                   maxGSSize = 150,
                   pvalueCutoff=1, #to get the table, filtering can be done after
                   pAdjustMethod="BH",
                   verbose=TRUE)
  gse.k_genename <- setReadable(gse.k, OrgDb = org.Hs.eg.db, keyType="UNIPROT")
  
  assign(paste0("gse.k", names_input[[i]]), gse.k_genename, envir = globalenv())
  
  write.table(gse.k_genename, file = paste0("../analysis/GSEA/",names_input[[i]], "_Kegg_all.txt"), sep = "\t", quote = F, 
              row.names = F, col.names = T)
  
  dot.k <- dotplot(gse.k_genename, showCategory=10, split=".sign") + facet_grid(.~.sign)
  
  tiff(filename = paste0("../analysis/GSEA/", names_input[[i]], "Kegg_dotplot.tiff"),
       width = 12 * 300, 
       height = 12 * 300,
       res = 300,
       compression = "lzw")
  print(dot.k)
  dev.off()
  
  #### enrichment map
  gse_pw.k <- pairwise_termsim(gse.r_genename)
  emap.k<-emapplot(gse_pw.k,
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
  
  tiff(filename = paste0("../analysis/GSEA/", names_input[[i]], "gosim_kegg.tiff"),
       width = 12 * 300, 
       height = 12 * 300,
       res = 300,
       compression = "lzw")
  print(emap.k)
}



#if you want to run for specific time point change heere.
for (j in 1:length(input)) {
  new.function(j)}


#new.function(6)

