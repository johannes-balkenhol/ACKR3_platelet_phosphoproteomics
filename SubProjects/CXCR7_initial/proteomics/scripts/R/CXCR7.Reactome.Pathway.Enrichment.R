############################################
#### script for annotation of phosphoproteom data 
## GO and KEGG enrichtment

BiocManager::install("org.Hs.eg.db")
BiocManager::install("reactome.db")
install.packages("sjmisc")
remotes::install_github("JosephCrispell/basicPlotteR")

#  Load packges ----
suppressPackageStartupMessages({
  library(calibrate)
  library(ggplot2)
  library(limma)
  library(directPA)
  library(reactome.db)
  library(annotate)
  library(PhosR)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(cowplot)
  require(dplyr)
  library(plyr)
  library('calibrate')
  library("basicPlotteR")
  library('sjmisc')
  library(stringr)
})



##################################################################
### Rectome enrichment
## accroding to https://pyanglab.github.io/PhosR/articles/PhosR.html
#input = list(top.collapse.10, top.collapse.30, top.collapse.60, top.collapse.300, top.collpase.600, 
# top.collapse.900, top.collapse.1800)


## define input  ----
input = list(top.10, top.30, top.60, top.300, top.600, top.900, top.1800)
names_input = c("10", "30", "60", "300", "600", "900", "1800")

for (i in 1:length(input)) {
  input[[i]] <- input[[i]][order(row.names(input[[i]])),]
}

Tc <- as.matrix(cbind(input[[1]]$logFC, input[[2]]$logFC, input[[3]]$logFC,
                      input[[4]]$logFC, input[[5]]$logFC, input[[6]]$logFC,
                      input[[7]]$logFC))

rownames(Tc) <- rownames(input[[1]])
colnames(Tc) <- names_input


rownames(Tc) <- mapIds(org.Hs.eg.db, keys=rownames(Tc), column="SYMBOL", keytype="UNIPROT", multiVals="first")



## select which dataset to choose ----
pathways = as.list(reactomePATHID2EXTID)
path_names = as.list(reactomePATHID2NAME)
name_id = match(names(pathways), names(path_names))
names(pathways) = unlist(path_names)[name_id]

pathways = pathways[which(grepl("Homo sapiens", names(pathways), ignore.case = TRUE))]

pathways = lapply(pathways, function(path) {
  gene_name = unname(getSYMBOL(path, data = "org.Hs.eg"))
  toupper(unique(gene_name))
})

## select up or downregulation ----
class = "UP/more"
order = "greater"
colors <- c("#FDAE61", "red")


class = "DOWN/less"
order = "less"
colors <- c("lightblue2", "blue")



#  ENRICHED REACTOME PATHWAYS; MAKE FIGURES; SAVE TO TIFF
for (i in 1:length(input)) {
  geneSet <- names(sort(Tc[,i],  #decreasing False for enrichment of downregulated
                        decreasing = T))[seq(round(nrow(Tc) * 0.05))]

  head(geneSet)
  
  
  path1 <- pathwayOverrepresent(geneSet, annotation=pathways, 
                                universe = rownames(Tc), alter = order)
  path2 <- pathwayRankBasedEnrichment(Tc[,i], 
                                      annotation=pathways, 
                                      alter = order)
  
  lp1 <- -log10(as.numeric(path1[names(pathways),1]))
  lp2 <- -log10(as.numeric(path2[names(pathways),1]))
  

  ######### Visulatization
  pathways2 <- as.data.frame(t(data.frame(lapply(pathways, function(x) length(x)))))
  path3 <- as.data.frame(path2)
  path4 <- merge(path3, pathways2, by = "row.names", all = FALSE)
  colnames(path4) <- c("pathway", "pvalue", "number.substrates", "substrates", "pw.size")
  path4 <- na.omit(path4)
  path4$number.substrates <- as.numeric(path4$number.substrates)
  path4$pw.size <- as.numeric(path4$pw.size)
  path4$pvalue <- as.numeric(path4$pvalue)
  path4$ratio <- path4$number.substrates/path4$pw.size
  path4$pathway <- gsub("_", " ", gsub("REACTOME_|Homo.sapiens..", "", path4$pathway))
  #path4 <- path4[order(path4$pvalue),]
  path4 <- path4[path4$pvalue<0.05,]
  path4 <- path4[order(path4$ratio, decreasing = T),]
  write.table(path4, file = paste0("../analysis/Reactome_enrichment/", class, names_input[[i]] ,"_pathways.txt"), 
              sep = '\t', row.names = FALSE)
  path4 <- path4[1:10,]
  path4$pvalue <- round(path4$pvalue, digits=4)
  path4 <- path4[order(path4$pvalue),]
  path4$pathway <- gsub("\\.", " ", path4$pathway)
  
  #dev.new()
  ggp <- ggplot(path4[10:1,], aes(x=1:10,y=ratio,fill=pvalue)) + 
    coord_flip() +
    geom_bar(stat = "identity") + 
    scale_fill_gradient(low = colors[1], high = colors[2], 
                        limits = c(min(path4[1:10,]$pvalue), max(path4[1:10,]$pvalue)), 
                        name = "pvalue", 
                        guide = guide_colorbar(barwidth = , 
                                               barheight = 10, 
                                               title.position = "top", 
                                               title.hjust = 0.5,
                                               reverse = TRUE))  +
    geom_text(aes(label = pathway, y = 0.005),
              vjust = 0, colour = "black",
              position = position_dodge(0.1), size =8, hjust = 'left') +
    labs(title = paste0(names_input[i]," sec"), size = 20, x = "pathways", y = "gene ratio") +
    theme_cowplot() +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          axis.text.x=element_text(size=20),
          axis.title.x = element_text(size =20),
          axis.title.y = element_text(size = 20)
    )
  tiff(filename = paste0("../analysis/Reactome_enrichment/", class, names_input[[i]] ,"_barplot.tiff"),
       width = 12 * 300, 
       height = 8 * 300,
       res = 300,
       compression = "lzw")
  
  print(ggp)
  dev.off()
  
  assign(paste0("path4.", names_input[[i]]), path4)
  
}






