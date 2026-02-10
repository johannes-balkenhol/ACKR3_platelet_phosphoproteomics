




############################################
#### script for annotation of phosphoproteom data 
## GO and KEGG enrichtment
##


BiocManager::install("org.Hs.eg.db")
BiocManager::install("reactome.db")

#### load packges


suppressPackageStartupMessages({
  library(calibrate)
  library(limma)
  library(directPA)
  library(org.Rn.eg.db)
  library(reactome.db)
  library(annotate)
  library(PhosR)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(cowplot)
})


desc= "top.10"

df_pl = top.1800[c("uniprot", "symbol", "AveExpr", "logFC", "adj.P.Val")]
df_pl = top.collapse.1800[c("uniprot_id", "name", "Average", "logFC", "PValue")]

#uniprot, gene_name, avrexpr, log2fold, pvlaue
#df_pl$uniprot = rownames(df_pl)
colnames(df_pl)<- c("Uniprot_ID", "Name", "Average", "Log2FoldChange", "PValue")


##################################################################
### go enrichment

#GO OVERREPRESENTATION ANALYSIS

# Exctract significant results (padj < 0.05)
sig_genes_df_pl = subset(df_pl, PValue < 0.05)
sig_genes_df_pl_up = subset(sig_genes_df_pl, Log2FoldChange > 0)
sig_genes_df_pl_down = subset(sig_genes_df_pl, Log2FoldChange < 0)
# From significant results, we want to filter on log2fold change
genes_pl <- sig_genes_df_pl$Log2FoldChange
genes_pl_up <- sig_genes_df_pl_up$Log2FoldChange
genes_pl_down <- sig_genes_df_pl_down$Log2FoldChange
# Name the vector
names(genes_pl) <- sig_genes_df_pl$Uniprot_ID
names(genes_pl_up) <- sig_genes_df_pl_up$Uniprot_ID
names(genes_pl_down) <- sig_genes_df_pl_down$Uniprot_ID
# omit NA values
genes_pl <- na.omit(genes_pl)
genes_pl <- names(genes_pl)
genes_pl_up <- na.omit(genes_pl_up)
genes_pl_up <- names(genes_pl_up)
genes_pl_down <- na.omit(genes_pl_down)
genes_pl_down <- names(genes_pl_down)

go_enrich_pl <- enrichGO(gene = genes_pl,
                          universe = names(gene_list_pl),
                          OrgDb = org.Hs.eg.db, 
                          keyType = 'UNIPROT',
                          readable = T,
                          ont = "BP",
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.05)


#go overrep plots: save result to file####

pl.tab = go_enrich_pl@result


write.table(pl.tab, file = paste0(desc, "_GOBP_CUSTOM_all.txt"), sep = "\t", quote = F, 
            row.names = F, col.names = T)



#go overrep plots: custom background####

#dotplot(go_enrich_pl, title = paste(desc, "GO Overrepresentation", sep=" "))
#barplot(go_enrich_pl, 
 #       drop = TRUE, 
  #      x= "GeneRatio",
   #     showCategory = 10, 
    #    title = paste(desc, "GO Overrepresentation", sep=" "),
     #   font.size = 18,
#)


pl_dot <- dotplot(go_enrich_pl, title = paste(desc,"GO Overrepresentation", sep=" "), font.size=19)

pl_bar <- barplot(go_enrich_pl, 
                   drop = TRUE, 
                   x= "GeneRatio",
                   showCategory = 10, 
                   title = paste(desc, "GO Overrepresentation", sep=" "),
                   font.size = 18,
)


#save custom plots####


plot_list <- list(pl_bar, pl_dot)
for (i in 1:2) {
  file_name = paste(desc,"_GOBP_CUSTOM_", i , ".tiff", sep="")
  save_plot(file_name, plot_list[[i]], base_width = 12, base_height = 10)
}





#KEGG PATHWAY OVERREPRESENTATION####

#kegg overrep input prep: id conversion, significant DEGs, upregulated and downregulated genes#######

#mk
ids_pl <-bitr(original_gene_list_pl, fromType = "UNIPROT", 
              toType = "ENTREZID", OrgDb="org.Hs.eg.db") # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)

dedup_ids_pl = ids_pl[!duplicated(ids_pl[c("UNIPROT")]),]
colnames(dedup_ids_pl) = c("Uniprot_ID", "ENTREZID")
# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2_pl = df_pl[df_pl$Uniprot_ID %in% dedup_ids_pl$Uniprot_ID,]
# Create a new column in df2 with the corresponding ENTREZ IDs
df2_pl <- merge(dedup_ids_pl, df2_pl, by = 'Uniprot_ID')
#df2_pl$Y = dedup_ids_pl$ENTREZID
# Create a vector of the gene unuiverse
kegg_gene_list_pl <- df2_pl$Log2FoldChange
# Name vector with ENTREZ ids
names(kegg_gene_list_pl) <- df2_pl$ENTREZID
# omit any NA values 
kegg_gene_list_pl<-na.omit(kegg_gene_list_pl)
# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list_pl = sort(kegg_gene_list_pl, decreasing = TRUE)
# Exctract significant results from df2
kegg_sig_genes_df_pl = subset(df2_pl, df2_pl$PValue < 0.05)
kegg_sig_genes_df_pl_up = subset(kegg_sig_genes_df_pl, kegg_sig_genes_df_pl$Log2FoldChange > 0)
kegg_sig_genes_df_pl_down = subset(kegg_sig_genes_df_pl, kegg_sig_genes_df_pl$Log2FoldChange < 0)
# From significant results, we want to filter on log2fold change
kegg_genes_pl <- kegg_sig_genes_df_pl$Log2FoldChange
kegg_genes_pl_up <- kegg_sig_genes_df_pl_up$Log2FoldChange
kegg_genes_pl_down <- kegg_sig_genes_df_pl_down$Log2FoldChange
# Name the vector with the CONVERTED ID!
names(kegg_genes_pl) <- kegg_sig_genes_df_pl$ENTREZID
names(kegg_genes_pl_up) <- kegg_sig_genes_df_pl_up$ENTREZID
names(kegg_genes_pl_down) <- kegg_sig_genes_df_pl_down$ENTREZID
# omit NA values
kegg_genes_pl <- na.omit(kegg_genes_pl)
kegg_genes_pl_up <- na.omit(kegg_genes_pl_up)
kegg_genes_pl_down <- na.omit(kegg_genes_pl_down)
kegg_genes_pl <- names(kegg_genes_pl)
kegg_genes_pl_up <- names(kegg_genes_pl_up)
kegg_genes_pl_down <- names(kegg_genes_pl_down)

  
##kegg overrep analysis: custom background####NOTHING SIGNIFICANT COMES OUT####

kegg_organism = "hsa"



kk_pl <- enrichKEGG(gene=kegg_genes_pl, universe=names(kegg_gene_list_pl),
                    organism=kegg_organism, pvalueCutoff = 0.05, keyType = "ncbi-geneid")





pl.tab_kk = kk_pl@result

write.table(pl.tab_kk, file = paste0(desc, "_KEGG_CUSTOM_all.txt"), sep = "\t", quote = F, 
            row.names = F, col.names = T)

pl_kk_bar <- barplot(pl.tab_kk, 
                     showCategory = 10, 
                     title = paste(desc, "Enriched KEGG Pathways", sep = " "),
                     font.size = 18)

pl_kk_dot <- dotplot(kk_mk, 
                     showCategory = 10, 
                     title = paste(desc, "Enriched KEGG Pathways", sep = " "),
                     font.size = 18)
plot_list2 <- list(pl_kk_bar, pl_kk_dot)
for (i in 1:2) {
  file_name = paste(desc,"_KEGG_CUSTOM_", i , ".tiff", sep="")
  save_plot(file_name, plot_list2[[i]], base_width = 12, base_height = 10)
}



##################################################################
### go enrichment with log2fold change


setwd("F:/Masterarbeit_BioWis/Proteomics/pipeline/analysis/GO enrichment")
filename= "go_logfc.txt"
desc <- tools::file_path_sans_ext(filename)

df_go = read.csv(filename, header=TRUE, sep="\t")

colnames(df_go)<- c("GO_name", "GO_domain", "NP", "NSP", "PSP", "avg.logfc", "sum.logfc", "avg.apm_change", "sum.apm_change")

# example by https://www.biostars.org/p/481038/
df_test <- data.frame(gene = c("A", "B", "C","E", "F", "G"), fc = c(-1,-2,-3,3,2,1), pval = c(0.01, 0.05, 0.09, 0.01,0.1,0.08))

df_go = sort(df_go, decreasing = TRUE)

ggplot(data=df_go, aes(x=reorder(GO_name, -avg.logfc), y=avg.logfc, fill = sum.logfc, )) +
geom_bar(stat="identity") + 
coord_flip()




##################################################################
### Rectome enrichment
## accroding to https://pyanglab.github.io/PhosR/articles/PhosR.html

suppressPackageStartupMessages({
  library(org.Rn.eg.db)
  library(reactome.db)
  library(org.Hs.eg.db)
  library(plyr)
})


## define dataframe with logFC's
Tc <- as.data.frame(cbind(top.10$logFC, top.600$logFC, top.1800$logFC))
rownames(Tc) <- rownames(top.10)
Tc <-  Tc[!sapply(strsplit(rownames(Tc), ";"), "[[", 3) == "",]

## reduce phosphoste to unique sites
Tc <-  Tc[!sapply(strsplit(rownames(Tc), ";"), "[[", 3) == "",]

gene_name <- sapply(strsplit(rownames(Tc), ";"), "[[", 2)
site <- sapply(strsplit(rownames(Tc), ";"), "[[", 3)
#rownames(Tc) <- gsub("(.*)(;[A-Z])([0-9]+)(;)", "\\1;\\3;", rownames(Tc))
Tc$ID <- paste(gene_name, site, sep= ";")

Tc <- ddply(Tc, "ID", summarize, V1 = mean(V1),
                             V2 = mean(V2),
                             V3 = mean(V3))
							 
#gene_name <- sapply(strsplit(rownames(Tc), ";"), "[[", 2)
#site <- sapply(strsplit(rownames(Tc), ";"), "[[", 3)
#rownames(Tc) <- paste(gene_name, site, sep= ";")

rownames(Tc) <- Tc$ID

Tc <- Tc[,-1]
colnames(Tc) <- c("top.10", "top.600", "top.1800")


#rownames(Tc) <- gsub("(.*)(;[A-Z])([0-9]+)(;)", "\\1;\\3;", rownames(Tc))



#Tc.gene <- phosCollapse(Tc, id=gsub(";.+", "", rownames(Tc)), 
#                        stat=apply(abs(Tc), 1, max), by = "max")
#rownames(Tc) <- sapply(strsplit(rownames(Tc), ";"), "[[", 2)
Tc.gene <- phosCollapse(Tc, id=sapply(strsplit(rownames(Tc), ";"), "[[", 1), 
                        stat=apply(abs(Tc), 1, max), by = "max")
geneSet <- names(sort(Tc.gene[,3], 
                        decreasing = TRUE))[seq(round(nrow(Tc.gene) * 0.1))]

head(geneSet)


pathways = as.list(reactomePATHID2EXTID)
path_names = as.list(reactomePATHID2NAME)
name_id = match(names(pathways), names(path_names))
names(pathways) = unlist(path_names)[name_id]

pathways = pathways[which(grepl("Homo sapiens", names(pathways), ignore.case = TRUE))]

pathways = lapply(pathways, function(path) {
    gene_name = unname(getSYMBOL(path, data = "org.Hs.eg"))
    toupper(unique(gene_name))
})



path1 <- pathwayOverrepresent(geneSet, annotation=pathways, 
                                universe = rownames(Tc.gene), alter = "greater")
path2 <- pathwayRankBasedEnrichment(Tc.gene[,3], 
                                    annotation=pathways, 
                                    alter = "greater")


lp1 <- -log10(as.numeric(path2[names(pathways),1]))
lp2 <- -log10(as.numeric(path1[names(pathways),1]))
par(mfrow=c(1,1))
plot(lp1, lp2, ylab="Overrepresentation (-log10 pvalue)", xlab="Rank-based enrichment (-log10 pvalue)", main="Comparison of 1D pathway analyses", xlim = c(0, 10))

# select highly enriched pathways
sel <- which(lp1 > 0.2 & lp2 > 0.2)
textxy(lp1[sel], lp2[sel], gsub("_", " ", gsub("REACTOME_", "", names(pathways)))[sel])





##directional analysis, dont really get it

library("directPA")
data("PhosphoSitePlus")
data("PhosphoELM")

# 2D direction site-centric kinase activity analyses
par(mfrow=c(1,1))
dpa1 <- directPA(Tc[,c(1,3)], direction=0, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dpa2 <- directPA(Tc[,c(1,2)], direction=pi*7/4, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dpa3 <- directPA(Tc[,c(2,3)], direction=pi*7/4, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
				 
				 
# top activated kinases
dpa1$pathways[1:5,]

dpa2$pathways[1:5,]

dpa3$pathways[1:5,]


z1 <- perturbPlot2d(Tc=Tc[,c(2,3)], 
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}),
                    cex=1, xlim=c(-2, 4), ylim=c(-2, 4), 
                    main="Kinase perturbation analysis")


###Kinase analysis

suppressPackageStartupMessages({
  
  library(PhosR)
  library(dplyr)
  library(ggplot2)
  library(GGally)
  library(ggpubr)
  library(calibrate)
  library(network)
})

data("KinaseMotifs")
data("KinaseFamily")






norm_intensity_10 <- norm_intensity[,c(1:10, 21:30)]
norm_intensity_1800 <- norm_intensity[,c(1:10, 41:50)]
norm_intensity_600 <- norm_intensity[,c(1:10, 61:70)]


ppe_mean <- meanAbundance(norm_intensity_1800, grps=gsub("_[0-9][0-9]", "", colnames(norm_intensity_1800)))
aov <- matANOVA(mat=norm_intensity_1800, grps=gsub("_[0-9][0-9]", "", colnames(norm_intensity_1800)))
idx <- (aov < 0.05) & (rowSums(ppe_mean > 0.5) > 0)
ppe_reg <- norm_intensity_1800[idx, ,drop = FALSE]

ppe.phos.std <- standardise(ppe_reg)
rownames(ppe.phos.std) <- paste0(GeneSymbol(ppe), ";", Residue(ppe), Site(ppe), ";")[idx]
ppe.phos.seq <- Sequence(ppe)[idx]

ppe.matrices <- kinaseSubstrateScore(substrate.list = PhosphoSite.human, 
                                    mat = ppe.phos.std, seqs =ppe.phos.seq, 
                                    numMotif = 5, numSub = 1, verbose = TRUE)



set.seed(1)
ppe.predMat <- kinaseSubstratePred(ppe.matrices, top=30, verbose = TRUE) 



plotKinaseNetwork(KSR = ppe.matrices, predMatrix = ppe.predMat, threshold = 0.9, color = kinase_all_color)
