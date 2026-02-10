############################################
#### script for annotation of phosphoproteom data 
## GO and KEGG enrichtment
##


BiocManager::install("org.Hs.eg.db")
BiocManager::install("reactome.db")
install.packages("sjmisc")
remotes::install_github("JosephCrispell/basicPlotteR")

#### load packges
suppressPackageStartupMessages({
  library(calibrate)
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
})


#################################################
##### prepare collpased psite tables for gsea
##### prepare collpased psite tables for gsea
input = list(top.10, top.600, top.1800, 
top.10.dmso.vs.0s, top.600.dmso.vs.0s, top.1800.dmso.vs.0s,
top.10.cxcr7.vs.0s, top.600.cxcr7.vs.0s, top.1800.cxcr7.vs.0s)

names_input = c("10", "600", "1800", 
"10.dmso.vs.0s", "600.dmso.vs.0s", "1800.dmso.vs.0s",
"10.cxcr7.vs.0s", "600.cxcr7.vs.0s", "1800.cxcr7.vs.0s")

for (i in 1:length(input)) {
    #df_pl = top.10.cxcr7.vs.0s
    df_pl = input[[i]]
    ## get rid of site without annotation
    df_pl = df_pl[-which(df_pl$psite == ""),]
    df_pl = df_pl[, c(7,8,2,1,5,9)]

    colnames(df_pl)<- c("uniprot_id", "name", "Average", "logFC", "PValue", "PSite")

    assign(paste0("top.filter.", names_input[[i]]), df_pl)

    ## collpase phosphosites to the protein
    top.collapse <- ddply(df_pl, .(uniprot_id, name), summarise,
                  Average = max(Average),
                  minlogFC = min(logFC),
                  maxlogFC = max(logFC),
                  AbslogFC = max(abs(logFC)),
                  meanlogFC = mean(logFC),
                  PValue = min(PValue))

    top.collapse$logFC <- 0
    top.collapse$logFC[abs(top.collapse[, "minlogFC"]) < abs(top.collapse[, "maxlogFC"])] = 
                                  top.collapse$maxlogFC[abs(top.collapse[, "minlogFC"]) < abs(top.collapse[, "maxlogFC"])]

    top.collapse$logFC[abs(top.collapse[, "minlogFC"]) > abs(top.collapse[, "maxlogFC"])] = 
                                  top.collapse$minlogFC[abs(top.collapse[, "minlogFC"]) > abs(top.collapse[, "maxlogFC"])]

    top.collapse$logFC[abs(top.collapse[, "minlogFC"]) == abs(top.collapse[, "maxlogFC"])] = 
                                  top.collapse$minlogFC[abs(top.collapse[, "minlogFC"]) == abs(top.collapse[, "maxlogFC"])]

    ## optinal: set the mean
    top.collapse$logFC = top.collapse$meanlogFC

    ## repeat for each input
    #colnames(top.collapse) <- c("Average","logFC","p.value","uniprot_id","name")
    rownames(top.collapse) <- paste(top.collapse.10$uniprot_id, top.collapse.10$name,sep=";")
    assign(paste0("top.collapse.", names_input[[i]]), top.collapse[,c(3,9,8,1,2)])
}


##################################################################
### Rectome enrichment
## accroding to https://pyanglab.github.io/PhosR/articles/PhosR.html
input = list(top.collapse.10, top.600, top.collapse.1800, 
top.collapse.10.dmso.vs.0s, top.collapse.600.dmso.vs.0s, top.collapse.1800.dmso.vs.0s,
top.collapse.10.cxcr7.vs.0s, top.collapse.600.cxcr7.vs.0s, top.collapse.1800.cxcr7.vs.0s)

## or

input = list(top.filter.10, top.filter.600, top.filter.1800, 
top.filter.10.dmso.vs.0s, top.filter.600.dmso.vs.0s, top.filter.1800.dmso.vs.0s,
top.filter.10.cxcr7.vs.0s, top.filter.600.cxcr7.vs.0s, top.filter.1800.cxcr7.vs.0s)


### order inputs
#  input.2 <- top.collapse.600[ order(row.names(top.collapse.600)), ]

for (i in 1:length(input)) {
  input[[i]] <- input[[i]][order(row.names(input[[i]])),]
}

###### Analysis
## define dataframe with logFC's
Tc <- as.data.frame(cbind(input[[4]]$logFC, input[[5]]$logFC, input[[6]]$logFC))
rownames(Tc) <- rownames(input[[1]])
colnames(Tc) <- c("10s", "600s", "1800s")

TcP <- as.data.frame(cbind(input[[4]]$PValue, input[[5]]$PValue, input[[6]]$PValue))
rownames(TcP) <- rownames(input[[1]])
colnames(TcP) <- c("10s", "600s", "1800s")


Tc.gene <- phosCollapse(Tc,   id=paste(sapply(strsplit(rownames(Tc), ";"), "[[", 2)), 
                        stat=apply(abs(Tc), 1, max), by = "max")

TcP.gene <- phosCollapse(TcP, id=paste(sapply(strsplit(rownames(TcP), ";"), "[[", 2)),
                        stat=apply(TcP, 1, min), by = "min")

#geneSet <- names(sort(Tc.gene[,1], 
#                        decreasing = FALSE))[seq(round(nrow(Tc.gene) * 0.025))]
geneSet <- names(sort(TcP.gene[,1], 
                        decreasing = FALSE))[TcP.gene[,1]<0.05]
#geneSet <- input[[1]]$name[input[[4]]$PValue<0.05]

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
                                universe = rownames(TcP.gene), alter = "greater")
path2 <- pathwayRankBasedEnrichment(Tc.gene[,1], 
                                    annotation=pathways, 
                                    alter = "greater")

lp1 <- -log10(as.numeric(path2[names(pathways),1]))
lp2 <- -log10(as.numeric(path1[names(pathways),1]))

dev.new()
#png(file=paste0("../analysis/reactome_enrichment/", "reactome_top.10",".png"), width = 5, height = 5, units = 'in', res = 1200)
#par(mfrow=c(1,1))
plot(lp1, lp2, ylab="Overrepresentation (-log10 pvalue)", xlab="Rank-based enrichment (-log10 pvalue)", main="", xlim = c(0, 5), 
ylim = c(0, 2.6), col = "grey", bg = "grey", )

# select highly enriched pathways
sel <- which(lp1 > 1.3 & lp2 >= 0.0)
#textxy(lp1[sel], lp2[sel], gsub("_", " ", gsub("REACTOME_|Mus musculus: ", "", names(pathways)))[sel], cex= 0.8,srt=0)
label <- gsub("_", " ", gsub("REACTOME_|Homo sapiens: ", "", names(pathways)))[sel]
label <- word_wrap(label, 45, linesep = NULL)
addTextLabels(lp1[sel], lp2[sel], label, cex.label=0.9, col.label="black", lty=2, col.line=rgb(0,0,0, 0.5),keepLabelsInside = TRUE)
dev.off()



rea1 <- pathways[sel][["Mus musculus: Microtubule-dependent trafficking of connexons from Golgi to the plasma membrane"]]
rea1 <- rea1[rea1 %in% geneSet]
write.table(rea1, "../../analysis/Reactome_enrichment/rea1.txt", sep="\t",  row.names=FALSE, quote=FALSE)


top.10[grep("RAC", rownames(top.10),]

#############################################
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