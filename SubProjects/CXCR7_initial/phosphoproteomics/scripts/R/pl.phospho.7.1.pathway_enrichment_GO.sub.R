




############################################
#### script for annotation of phosphoproteom data 
## GO and KEGG enrichtment
##


BiocManager::install("org.Hs.eg.db")
BiocManager::install("reactome.db")
remotes::install_github("YuLab-SMU/createKEGGdb")


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
  library(createKEGGdb)
})

i=6
#setwd("C:/Users/Ozge/Nextcloud/Phosphoproteom/phosphoproteom analysis 2/scripts")
names_input = c("10", "30", "60", "300", "600", "900", "1800")
input = list(top.collapse.10, top.collapse.30, top.collapse.60, top.collapse.300, top.collapse.600, 
             top.collapse.900, top.collapse.1800)
#names_input = c("10")
for (i in 1:length(names_input)) {
  df_pl <- input[[i]]
  desc= names_input[[i]]
  #df_pl = top.10[c("uniprot", "symbol", "AveExpr", "logFC", "adj.P.Val")]
  df_pl = df_pl[c("uniprot_id", "name", "Average", "logFC", "PValue")]
  
  #uniprot, gene_name, avrexpr, log2fold, pvlaue
  #df_pl$uniprot = rownames(df_pl)
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
  #original_gene_list_pl  <- sapply(strsplit(raw_abundance$peptide_id, ";"), "[[", 1) ##ÖO I use everything that was measured, from the raw abundance data
  gene_list_pl<-na.omit(original_gene_list_pl)
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list_pl = sort(gene_list_pl, decreasing = TRUE)
  
  # Exctract significant results (padj < 0.05)
  sig_genes_df_pl = subset(df_pl, PValue < 0.05)
  sig_genes_df_pl = subset(sig_genes_df_pl, abs(Log2FoldChange) > 0.5)
  sig_genes_df_pl_up = subset(sig_genes_df_pl, Log2FoldChange > 0.5)
  sig_genes_df_pl_down = subset(sig_genes_df_pl, Log2FoldChange < -0.5)
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
  
  
  #genes_10sec_agonist_exclusive <- c("P51397","Q00577","Q9UEW8","Q9H3Z4","P07359","P05121","Q9UHY8","O75379","Q15049","Q9UDT6","P61978","Q96C24","Q9UNZ2","Q9Y2L6","Q9HBI1","Q9UIB8","O75122","Q8ND56","Q96Q42","Q53ET0","Q13541","Q05682","Q96SB3","Q9HCH0","Q27J81","Q6P1L5","Q96JH8","Q7L7X3","P31749","Q9Y490","P23528","O95810","Q96HC4","O75427","Q9BZQ8","Q86UU1","Q9ULH1","Q7L591","Q9C0C9","Q8N4C8","Q96AG3","O43561","O43639","P31751","Q06210","Q96IG2","Q8IZD0","Q05209","P50395","Q14847","P27338","Q5TCZ1","O15126","Q96Q42","Q9C0B5","Q13418","Q9H5N1","Q9NQG7","Q9BYI3","Q13233","Q9P289","Q86X10","Q9Y2L6","Q8TEW0","Q99719","Q5JSH3","Q9HBL0","O43399")
  #genes_600sec_agonist_exclusive <- c("Q9C0C9","Q9UDY2","Q9NWQ8","Q4G0F5","Q01449","Q9BZ72","Q5SW79","O60573","Q9Y6D5","Q8IXS8","Q5VZ89","Q9NZN5","Q96P48","Q8NEN9","P12931","P50851","Q99961","Q86UX7","O15027","Q8IZD0","P21291","O43741","Q9H4L5","P21333","Q13586","Q7Z422","Q5M8T2","O15013","Q6Q0C0","Q92619","Q6WCQ1","Q13283","Q14432","Q9Y4E1","Q8WXF7","P29692","Q8N699","Q99501","P06396","Q86VR2","Q96Q42","Q8NEN9","Q0JRZ9","Q9HBL0","O60841","Q9HC56","Q86VQ1","O75116","O43150","Q9NUY8","P10909","Q8TDB6","Q8WYL5","Q9Y2J2","Q2PPJ7","P16157","Q86YW5","P07359","P22059","Q9Y2Q0","Q9HBL0")
  #genes_1800sec_agonist_exclusive <- c("Q9BWH6","Q9H0B6","Q96KC8","Q9BYI3","Q9UKG1","Q6PJF5","P00488","Q9HBI1","P27816","Q96D71","Q9NUP1","Q7Z460","Q9P266","Q92835","Q13615","Q9NRF8","A0FGR8","Q9P035","Q9UDY2","Q8WVT3","Q9C0D7","Q53ET0","Q5T5U3","Q8IWW6","Q8N9U0","P35611","Q6VY07","Q9NZN3","Q9NYI0","O43765","O60229","Q9Y6D5","O95425","O14639","P20645","Q9UDY2","P49585","Q9H4M9","Q9H6S0","Q8WW12","P23528","Q66K74","Q6ZS17","Q99961","Q9C0C9","O15013","Q92539","Q14432","O43639","Q7Z460","Q96NA2","Q9UBS8","P11274","Q9Y2L6","P07996","P02775","Q7RTP6","P00747","Q6R327","Q14699","O75376","Q5VY43","Q9C0B0","O60678","P61026","P22314","Q6UN15","Q8WYL5","Q7L1W4","Q15036","Q5TBA9","Q9Y2X7","Q9BZ67","O95425","Q9UQL6","P50502","Q8NFH8","O43488","Q7L4E1","Q9Y608","Q674X7","P26639","O15427","P17252","Q96SB3","Q9Y3L3","O60763","O14745","P04792","P54105","O94929","Q15746","Q13586","P20645","Q92974","Q8WYR4","Q16799","Q8N4C8","Q9Y2L6","Q13442","Q5UE93","Q15691","Q99698","Q7Z401","Q99501")
  
  #go overrep analysis: custom background####
  
  #genes_pl <- genes_600sec_agonist_exclusive
  #desc= "top.600_agonist excusive"
  
  go_enrich_pl <- enrichGO(gene = genes_pl,
                           universe = names(gene_list_pl),
                           OrgDb = org.Hs.eg.db, 
                           keyType = 'UNIPROT',
                           readable = T,
                           ont = "ALL",
                           pvalueCutoff = 0.05, 
                           qvalueCutoff = 0.2,
                           pAdjustMethod = "BH")
  
  #assign(paste0("go_enrich_pl", names_input[[i]]), go_enrich_pl)
  #go overrep plots: save result to file####
  
  pl.tab = go_enrich_pl@result
  
  
  write.table(pl.tab, file = paste0("../analysis/Gene Ontology ORA/", desc, "_GOBP_CUSTOM_all.txt"), sep = "\t", quote = F, 
              row.names = F, col.names = T)
  
  
  
  #go overrep plots: custom background####
  
  
  tiff(filename = paste0("../analysis/Gene Ontology ORA/", desc, "dotplot.tiff"),
       width = 10 * 300, 
       height = 10 * 300,
       res = 300,
       compression = "lzw")
  
  dotplot(go_enrich_pl, title = paste(desc,"GO Overrepresentation", sep=" "), font.size=19)
  
  dev.off()
  
  
  tiff(filename = paste0("../analysis/Gene Ontology ORA/", desc, "barplot.tiff"),
       width = 10 * 300, 
       height = 10 * 300,
       res = 300,
       compression = "lzw")
  
  barplot(go_enrich_pl, 
          drop = TRUE, 
          x= "GeneRatio",
          showCategory = 10, 
          title = paste(desc, "GO Overrepresentation", sep=" "),
          font.size = 18,
  )
  
  dev.off()
  
  #save custom plots####
  
  
  #plot_list <- list(pl_bar, pl_dot)
  #for (i in 1:2) {
  # file_name = paste("analysis/Gene Ontology ORA/",desc,"_GOBP_CUSTOM_", i , ".tiff", sep="")
  #save_plot(file_name, plot_list[[i]], base_width = 12, base_height = 10)}
  
  
  
  
  #KEGG PATHWAY OVERREPRESENTATION####
  
  #kegg overrep input prep: id conversion, significant DEGs, upregulated and downregulated genes#######
  
  #mk
  
  
  # install the packages
  #remotes::install_github("YuLab-SMU/clusterProfiler") 
  #remotes::install_github("YuLab-SMU/createKEGGdb")
  # import the library and create a KEGG database locally 
  
  
  #species <-"hsa"
  #createKEGGdb::create_kegg_db(species)
  # You will get KEGG.db_1.0.tar.gz file in your working directory
  
  #keytypes(org.Hs.eg.db)
  
  
  
  
  ids_pl <-bitr(names(original_gene_list_pl), fromType = "UNIPROT",
                toType = "ENTREZID", OrgDb=org.Hs.eg.db) # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
  
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
  
  write.table(pl.tab_kk, file = paste0("../analysis/KEGG ORA/",desc, "_KEGG_CUSTOM_all.txt"), sep = "\t", quote = F, 
              row.names = F, col.names = T)
  
  
  tiff(filename = paste0("../analysis/KEGG ORA/", desc, "dotplot.tiff"),
       width = 10 * 300, 
       height = 10 * 300,
       res = 300,
       compression = "lzw")
  
  dotplot(kk_pl, title = paste(desc,"KEGG Enriched Pathways", sep=" "), font.size=19)
  
  dev.off()
  
  
  tiff(filename = paste0("../analysis/KEGG ORA/", desc, "barplot.tiff"),
       width = 10 * 300, 
       height = 10 * 300,
       res = 300,
       compression = "lzw")
  
  barplot(kk_pl, 
          drop = TRUE, 
          x= "GeneRatio",
          showCategory = 10, 
          title = paste(desc, "KEGG Enriched Pathways", sep=" "),
          font.size = 18,
  )
  
  dev.off()
  
}















# plot_list2 <- list(pl_kk_bar, pl_kk_dot)
# for (i in 1:2) {
#   file_name = paste(desc,"_KEGG_CUSTOM_", i , ".tiff", sep="")
#   save_plot(file_name, plot_list2[[i]], base_width = 12, base_height = 10)
# }



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




