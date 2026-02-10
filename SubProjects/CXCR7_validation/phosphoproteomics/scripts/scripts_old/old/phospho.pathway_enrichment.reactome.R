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
  require(dplyr)
  library(plyr)
  library(reactome.db)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library('calibrate')
  library("basicPlotteR")
  library('sjmisc')
})




setwd("F:/Masterarbeit_BioWis/Phosphoproteomics/Phosphoproteom validation/scripts/")




##################################################################
### Rectome enrichment
## accroding to https://pyanglab.github.io/PhosR/articles/PhosR.html






## define dataframe with logFC's
Tc <- as.data.frame(cbind(top.10$logFC, top.600$logFC, top.1800$logFC))
rownames(Tc) <- rownames(ppe)

## reduce phosphoste to unique sites
Tc <-  Tc[!sapply(strsplit(rownames(Tc), ";"), "[[", 3) == "",]

gene_name <- sapply(strsplit(rownames(Tc), ";"), "[[", 2)
site <- sapply(strsplit(rownames(Tc), ";"), "[[", 3)
#rownames(Tc) <- paste(gene_name, site, sep= ";")
#rownames(Tc) <- gsub("(.*)(;[A-Z])([0-9]+)(;)", "\\1;\\3;", rownames(Tc))
Tc$ID <- paste(gene_name, site, sep= ";")

Tc <- ddply(Tc, "ID", summarize, V1 = mean(V1),
                             V2 = mean(V2),
                             V3 = mean(V3))

rownames(Tc) <- Tc$ID

Tc <- Tc[,-1]
colnames(Tc) <- c("top.10", "top.600", "top.1800")


#rownames(Tc) <- gsub("(.*)(;[A-Z])([0-9]+)(;)", "\\1;\\3;", rownames(Tc))

########JB use other collpase function

#Tc.gene <- phosCollapse(Tc, id=gsub(";.+", "", rownames(Tc)), 
#                        stat=apply(abs(Tc), 1, max), by = "max")
#rownames(Tc) <- sapply(strsplit(rownames(Tc), ";"), "[[", 2)
Tc.gene <- phosCollapse(Tc, id=sapply(strsplit(rownames(Tc), ";"), "[[", 1), 
                        stat=apply(abs(Tc), 1, max), by = "max")
geneSet <- names(sort(Tc.gene[,1], 
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
path2 <- pathwayRankBasedEnrichment(Tc.gene[,1], 
                                    annotation=pathways, 
                                    alter = "greater")

dev.new()
#png(filename="../analysis/Reactome_enrichment/reactome_enrichment.png") 
lp1 <- -log10(as.numeric(path2[names(pathways),1]))
lp2 <- -log10(as.numeric(path1[names(pathways),1]))
par(mfrow=c(1,1))
plot(lp1, lp2, ylab="Overrepresentation (-log10 pvalue)", xlab="Rank-based enrichment (-log10 pvalue)", main="Comparison of 1D pathway analyses top.1800")

# select highly enriched pathways
sel <- which(lp1 > 0.3 & lp2 > 0.1)
textxy(lp1[sel], lp2[sel], gsub("_", " ", gsub("REACTOME_", "", names(pathways)))[sel], cex= 0.8, srt=-20)
#dev.off()



dev.new()
#png(file=paste0("../analysis/reactome_enrichment/", "reactome_rhoa_vs_cdc42",".png"), width = 5, height = 5, units = 'in', res = 1200)
#par(mfrow=c(1,1))
plot(lp1, lp2, ylab="Overrepresentation (-log10 pvalue)", xlab="Rank-based enrichment (-log10 pvalue)", main="", xlim = c(0, 5), ylim = c(0, 2.6), col = "grey", bg = "grey", )

# select highly enriched pathways
sel <- which(lp1 > 1.3 & lp2 > 0.9)
#textxy(lp1[sel], lp2[sel], gsub("_", " ", gsub("REACTOME_|Mus musculus: ", "", names(pathways)))[sel], cex= 0.8,srt=0)
label <- gsub("_", " ", gsub("REACTOME_|Mus musculus: ", "", names(pathways)))[sel]
label <- word_wrap(label, 45, linesep = NULL)
addTextLabels(lp1[sel], lp2[sel], label, cex.label=0.9, col.label="black", lty=2, col.line=rgb(0,0,0, 0.5),keepLabelsInside = TRUE)
dev.off()




rea1 <- pathways[sel][["Mus musculus: Microtubule-dependent trafficking of connexons from Golgi to the plasma membrane"]]
rea1 <- rea1[rea1 %in% geneSet]
write.table(rea1, "../../analysis/Reactome_enrichment/rea1.txt", sep="\t",  row.names=FALSE, quote=FALSE)


top.10[grep("RAC", rownames(top.10),]


# 2D direction site-centric kinase activity analyses
par(mfrow=c(1,3))
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


z1 <- perturbPlot2d(Tc=Tc[,c(2,3)], 
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}),
                    cex=1, xlim=c(-2, 4), ylim=c(-2, 4), 
                    main="Kinase perturbation analysis")