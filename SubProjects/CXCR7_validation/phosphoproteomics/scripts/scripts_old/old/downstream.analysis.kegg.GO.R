############################################
#### script for annotation of phosphoproteom data 
## GO and KEGG enrichtment
##




#### load packges

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(cowplot)


#### load data
#setwd("F:/Masterarbeit_BioWis/Proteomics/pipeline/analysis/GO enrichment")
#filename= "top_pl_strap_norm_name.txt"
#desc <- tools::file_path_sans_ext(filename)
#### reading in input from deseq2####
#df_pl = read.csv(filename, header=TRUE, sep="\t")


df_pl = top.600
df_pl = df_pl[, c(7,7,2,1,5)] #this is for mutant files and files with more than 5 columns
#uniprot, gene_name, avrexpr, log2fold, pvlaue
#df_pl = df_pl[, c(1,2,18,19,20)] #this is for local files and files with more than 5 columns
df_pl$uniprot.1 = rownames(df_pl)

colnames(df_pl)<- c("Uniprot_ID", "Name", "Average", "Log2FoldChange", "PValue")



##################################################################
### go enrichment

#GO OVERREPRESENTATION ANALYSIS


#GENE ONTOLOGY BIOLOGICAL PROCESS OVERREPRESENTATION ANALYSIS####

#go overrep input prep: significant DEGs, upregulated and downregulated genes#####

#mk
# we want the log2 fold change 
original_gene_list_pl <- df_pl$Log2FoldChange
# name the vector
names(original_gene_list_pl) <- df_pl$Uniprot_ID
# omit any NA values 
gene_list_pl<-na.omit(original_gene_list_pl)
# sort the list in decreasing order (required for clusterProfiler)
gene_list_pl = sort(gene_list_pl, decreasing = TRUE)
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



#go overrep analysis: custom background####

go_enrich_pl <- enrichGO(gene = genes_pl_up,
                          universe = names(gene_list_pl),
                          OrgDb = org.Mm.eg.db, 
                          keyType = 'UNIPROT',
                          readable = T,
                          ont = "BP",
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.05)  


#go overrep plots: save result to file####

pl.tab = go_enrich_pl@result


write.table(pl.tab, file = paste0(filename, "_GOBP_CUSTOM_all.txt"), sep = "\t", quote = F, 
            row.names = F, col.names = T)



#go overrep plots: custom background####

pl_dot <- dotplot(go_enrich_pl, title = paste(desc,"GO Overrepresentation", sep=" "))

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