################################################################
##### Export for GSEA
##### online guide https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
##### AND
##### https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html


# Load Packages ----
suppressPackageStartupMessages({
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
})

# Annotation ----
organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# Prepare input  ----
names_input <- list("top.10", "top.600", "top.1800", 
                    "top.10.cxcr7.vs.0s", "top.600.cxcr7.vs.0s", "top.1800.cxcr7.vs.0s",
                    "top.10.dmso.vs.0s", "top.600.dmso.vs.0s", "top.1800.dmso.vs.0s")

input <- list(top.10, top.600, top.1800, 
                  top.10.cxcr7.vs.0s, top.600.cxcr7.vs.0s, top.1800.cxcr7.vs.0s,
                  top.10.dmso.vs.0s, top.600.dmso.vs.0s, top.1800.dmso.vs.0s)

for (i in 1:length(input)) {
  input[[i]] <- input[[i]][order(row.names(input[[i]])),]
}

Tc <- as.matrix(cbind(input[[1]]$logFC, input[[2]]$logFC, input[[3]]$logFC,
                      input[[4]]$logFC, input[[5]]$logFC, input[[6]]$logFC,
                      input[[7]]$logFC, input[[8]]$logFC, input[[9]]$logFC))
rownames(Tc) <- rownames(input[[1]])
colnames(Tc) <- names_input

# GSEA function ----
GSEA.function <- function(i) {
  # Choose dataset ----
  logFC <- Tc[,i]
  uniprot <-  rownames(Tc)
  df <- data.frame(uniprot, logFC)
  
  # Prepare input with Uniprot ----
  original_gene_list <- logFC
  names(original_gene_list) <- uniprot
  gene_list_uniprot<-na.omit(original_gene_list)
  # sort the list in decreasing order ---- 
  gene_list_uniprot = sort(gene_list_uniprot, decreasing = TRUE) #(required for clusterProfiler)
  # Map ids from Uniprot to EntrezID ----
  ids<-bitr(names(original_gene_list), fromType = "UNIPROT", toType = "ENTREZID", OrgDb=organism)
  dedup_ids = ids[!duplicated(ids[c("UNIPROT")]),]
  df2 = df[df$uniprot %in% dedup_ids$UNIPROT,]
  colnames(dedup_ids) = c("uniprot", "ENTREZID")
  df3 <- as.data.frame(merge(df2, dedup_ids, by="uniprot"))
  gene_list_entrez <- df3$logFC
  names(gene_list_entrez) <- df3$ENTREZID
  gene_list_entrez = sort(gene_list_entrez, decreasing = TRUE)
  
  # GSEA GO ----
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
  
  ## Save Results ----
  assign(paste0("gse.go", names_input[[i]]), gse.go_genename)
  write.table(gse.go_genename, file = paste0("../analysis/GSEA/",names_input[[i]], "_GO_all.txt"), sep = "\t", quote = F, 
              row.names = F, col.names = T)
  ## Plot Results ----
  dot<- dotplot(gse.go_genename, showCategory=10, split=".sign") + facet_grid(.~.sign)
  
  tiff(filename = paste0("../analysis/GSEA/", names_input[[i]], "GOBP_dotplot.tiff"),
       width = 10 * 300, 
       height = 10 * 300,
       res = 300,
       compression = "lzw")
  print(dot)
  dev.off()
  
  # GSEA Reactome ----
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
  
  ## Save results ----
  assign(paste0("gse.r", names_input[[i]]), gse.r_genename, envir = globalenv())
  write.table(gse.r_genename, file = paste0("../analysis/GSEA/",names_input[[i]], "_Reactome_all.txt"), sep = "\t", quote = F, 
              row.names = F, col.names = T)
  
  ## Plot results ----
  dot.r<-dotplot(gse.r_genename, showCategory=10, split=".sign") + facet_grid(.~.sign)
  
  tiff(filename = paste0("../analysis/GSEA/", names_input[[i]], "Reactome_dotplot.tiff"),
       width = 10 * 300, 
       height = 10 * 300,
       res = 300,
       compression = "lzw")
  print(dot.r)
  dev.off()
  
  # GSEA KEGG ----
  gse.k <- gseKEGG(geneList     = gene_list_uniprot,
                   organism     = 'hsa',
                   keyType = "uniprot",
                   minGSSize    = 10,
                   maxGSSize = 150,
                   pvalueCutoff=1, #to get the table, filtering can be done after
                   pAdjustMethod="BH",
                   verbose=TRUE)
  
  gse.k_genename <- setReadable(gse.k, OrgDb = org.Hs.eg.db, keyType="UNIPROT")
  
  ## Save results ----
  assign(paste0("gse.k", names_input[[i]]), gse.k_genename, envir = globalenv())
  write.table(gse.k_genename, file = paste0("../analysis/GSEA/",names_input[[i]], "_Kegg_all.txt"), sep = "\t", quote = F, 
              row.names = F, col.names = T)
  
  ## Plot results ----
  dot.k <- dotplot(gse.k_genename, showCategory=10, split=".sign") + facet_grid(.~.sign)
  
  tiff(filename = paste0("../analysis/GSEA/", names_input[[i]], "Kegg_dotplot.tiff"),
       width = 10 * 300, 
       height = 10 * 300,
       res = 300,
       compression = "lzw")
  print(dot.k)
  dev.off()
}

# Run GSEA.funtion for all timepoints ----
for (j in 1:length(input)) {
  GSEA.function(j)
}

# Run GSEA.funtion for individual timepoints ----
GSEA.function(1)


