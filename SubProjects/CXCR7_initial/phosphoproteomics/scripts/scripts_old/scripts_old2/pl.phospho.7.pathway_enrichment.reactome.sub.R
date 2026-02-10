############################################
#### script for annotation of phosphoproteom data 
## GO and KEGG enrichtment

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



##################################################################
### Rectome enrichment
## accroding to https://pyanglab.github.io/PhosR/articles/PhosR.html
#input = list(top.collapse.10, top.collapse.30, top.collapse.60, top.collapse.300, top.collpase.600, 
            # top.collapse.900, top.collapse.1800)


## or

input = list(top.filter.10, top.filter.30, top.filter.60, top.filter.300, top.filter.600,
             top.filter.900, top.filter.1800)


names_input = c("10", "30", "60", "300", "600", "900", "1800")


### order inputs
#  input.2 <- top.collapse.600[ order(row.names(top.collapse.600)), ]

for (i in 1:length(input)) {
  input[[i]] <- input[[i]][order(row.names(input[[i]])),]
}

############ Analysis
## define dataframe with logFC's
Tc <- as.data.frame(cbind(input[[1]]$logFC, input[[2]]$logFC, input[[3]]$logFC,
                          input[[4]]$logFC, input[[5]]$logFC, input[[6]]$logFC,
                          input[[7]]$logFC))
rownames(Tc) <- rownames(input[[4]])
colnames(Tc) <- names_input

TcP <- as.data.frame(cbind(input[[1]]$PValue, input[[2]]$PValue, input[[3]]$PValue,
                           input[[4]]$PValue, input[[5]]$PValue, input[[6]]$PValue,
                           input[[7]]$PValue))
rownames(TcP) <- rownames(input[[4]])
colnames(TcP) <- names_input


## select which dataset to choose
pathways = as.list(reactomePATHID2EXTID)
path_names = as.list(reactomePATHID2NAME)
name_id = match(names(pathways), names(path_names))
names(pathways) = unlist(path_names)[name_id]

pathways = pathways[which(grepl("Homo sapiens", names(pathways), ignore.case = TRUE))]

pathways = lapply(pathways, function(path) {
    gene_name = unname(getSYMBOL(path, data = "org.Hs.eg"))
    toupper(unique(gene_name))
})


#  ENRICHED REACTOME PATHWAYS; MAKE FIGURES; SAVE TO TIFF
for (i in 1:length(input)) {

  Tc.gene <- phosCollapse(Tc, id=paste(sapply(strsplit(rownames(Tc), ";"), "[[", 2)), 
                          stat=apply(abs(Tc), 1, max), by = "max")

  TcP.gene <- phosCollapse(TcP, id=paste(sapply(strsplit(rownames(TcP), ";"), "[[", 2)),
                          stat=apply(TcP, 1, min), by = "min")

  #geneSet <- names(sort(Tc.gene[,1], 
  #                        decreasing = FALSE))[seq(round(nrow(Tc.gene) * 0.025))]
  geneSet <- names(sort(TcP.gene[,i], 
                          decreasing = FALSE))[TcP.gene[,i]<0.05]
  #geneSet <- input[[1]]$name[input[[4]]$PValue<0.05]

  head(geneSet)


  path1 <- pathwayOverrepresent(geneSet, annotation=pathways, 
                                  universe = rownames(TcP.gene), alter = "greater")
  path2 <- pathwayRankBasedEnrichment(Tc.gene[,i], 
                                      annotation=pathways, 
                                      alter = "greater")

  lp1 <- -log10(as.numeric(path1[names(pathways),1]))
  lp2 <- -log10(as.numeric(path2[names(pathways),1]))

  

  ######### Visualiztaion1
  #dev.new()
  tiff(filename = paste0("analysis/Reactome_enrichment/", names_input[[i]] ,"_plot.tiff"),
    width = 5 * 300, 
    height = 5 * 300,
    res = 300,
    compression = "lzw")
  #png(file=paste0("../analysis/reactome_enrichment/", "reactome_top.10",".png"), width = 5, height = 5, units = 'in', res = 1200)
  #par(mfrow=c(1,1))
  plot(lp1, lp2, xlab="Overrepresentation (-log10 pvalue)", ylab="Rank-based enrichment (-log10 pvalue)", main="", xlim = c(0, 5), 
  ylim = c(0, 2.6), col = "grey", bg = "grey", )

  # select highly enriched pathways
  sel <- which(lp1 > 0.3 & lp2 >= 0.3)
  #textxy(lp1[sel], lp2[sel], gsub("_", " ", gsub("REACTOME_|Mus musculus: ", "", names(pathways)))[sel], cex= 0.8,srt=0)
  label <- gsub("_", " ", gsub("REACTOME_|Homo sapiens: ", "", names(pathways)))[sel]
  label <- word_wrap(label, 45, linesep = NULL)
  addTextLabels(lp1[sel], lp2[sel], label, cex.label=0.9, col.label="black", lty=2, col.line=rgb(0,0,0, 0.5),keepLabelsInside = TRUE)
  dev.off()

  ######### Visulatization2
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
  path4 <- path4[order(path4$pvalue),]
  path4$pvalue <- round(path4$pvalue, digits=4)

  #dev.new()

  ggp <- ggplot(path4[10:1,], aes(x=1:10,y=ratio,fill=pvalue)) + 
                    coord_flip() +
                    geom_bar(stat = "identity") + 
                    scale_fill_gradient(low = "red", high = "blue", 
                        limits = c(min(path4[1:10,]$pvalue), max(path4[1:10,]$pvalue)), 
                        name = "pvalue", 
                        guide = guide_colorbar(barwidth = , 
                                              barheight = 10, 
                                              title.position = "top", 
                                              title.hjust = 0.5,
                                              reverse = TRUE))  +
                    geom_text(aes(label = pathway, y = 0.005),
                    , vjust = 0, colour = "black",
                    position = position_dodge(0.1), size = 6, hjust = 'left') +
                    labs(title = "", size = 8, x = "pathways", y = "gene ratio") +
                    theme_cowplot() +
                    theme(axis.text.y = element_blank(), 
                          axis.ticks.y = element_blank(),
                          axis.text.x=element_text(size=10)
                    )
  tiff(filename = paste0("analysis/Reactome_enrichment/", names_input[[i]] ,"_barplot.tiff"),
    width = 12 * 300, 
    height = 8 * 300,
    res = 300,
    compression = "lzw")

  print(ggp)
  dev.off()
  
  assign(paste0("path4.", names_input[[i]]), path4)

}








#### extraxt certain pathways
rea1 <- pathways["Homo sapiens: Signaling by ROBO receptors"]
path2["Homo sapiens: Signaling by ROBO receptors"]
rea1 <- rea1[rea1 %in% geneSet]
write.table(rea1, "../../analysis/Reactome_enrichment/rea1.txt", sep="\t",  row.names=FALSE, quote=FALSE)

#### extract certain proteins
top.all[grep("ARHGA", rownames(top.all)),]



barplot(path3$number ~ path3$pvalue, horiz = TRUE)


pl_bar <- barplot(go_enrich_pl, 
                   drop = TRUE, 
                   x= "GeneRatio",
                   showCategory = 10, 
                   title = paste(desc, "GO Overrepresentation", sep=" "),
                   font.size = 18,
)





#############################################
#############################################
###### 2D direction site-centric kinase activity analyses
data("PhosphoSitePlus")
data("phospho_L6_ratio_pe")
data("SPSs")

Tc1 <- Tc
rownames(Tc1) <- paste(sapply(strsplit(rownames(Tc), ";"), "[[", 2), sapply(strsplit(rownames(Tc), ";"), "[[", 3), "", sep = ";")
rownames(Tc1) <- lapply(rownames(Tc1), function(x){gsub(";[STY]", ";", x)})




#dev.new()

tiff(filename = paste0("analysis/DPA/DPA.tiff"),
     width = 14 * 300, 
     height = 4 * 300,
     res = 300,
     compression = "lzw")

par(mfrow=c(1,3))
dpa1 <- directPA(Tc1[,c(1,5)], direction=0, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dpa2 <- directPA(Tc1[,c(1,7)], direction=0, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dpa3 <- directPA(Tc1[,c(5,7)], direction=0, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dev.off()

##### match psites in data and in phophositeplus
library(rlist)
test <- lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)})
for (i in 1:length(test)) {
  if (i == 1 & i <= 2) { 
    test2 = as.list(test[[i]])
  } else {
    print(as.list(test[[i]]))
    test2 <- append(test2, as.list(test[[i]]))
  }
}

which(rownames(Tc1) %in% test2)

sites.psite <- as.data.frame(t(as.data.frame(test2)))
sites.data <- as.data.frame(rownames(Tc1))
colnames(sites.psite) <- "psite"
colnames(sites.data) <- "psite"

which(sites.psite$psite %in% sites.data$psite)

sites.data$psite <- str_trim(sites.data$psite, side = c("both"))
sites.psite$psite <- str_trim(sites.psite$psite, side = c("both"))

########### top activated kinases
dpa1$pathways[1:5,]

dpa2$pathways[1:5,]

dpa3$pathways[1:5,]


tiff(filename = paste0("analysis/DPA/Kinase_pert.tiff"),
     width = 14 * 300, 
     height = 4 * 300,
     res = 300,
     compression = "lzw")

par(mfrow=c(1,3))


z1 <- perturbPlot2d(Tc=Tc1[,c(1,5)], 
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}),
                    cex=1, xlim=c(-2, 4), ylim=c(-2, 4), 
                    main="Kinase perturbation analysis")

z2 <- perturbPlot2d(Tc=Tc1[,c(1,7)], 
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}),
                    cex=1, xlim=c(-2, 4), ylim=c(-2, 4), 
                    main="Kinase perturbation analysis")

z3 <- perturbPlot2d(Tc=Tc1[,c(5,7)], 
                  annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}),
                    cex=1, xlim=c(-2, 4), ylim=c(-2, 4), 
                    main="Kinase perturbation analysis")
dev.off()

#############################################
#############################################
###### 2D direction site-centric pathway analyses
data("PhosphoSitePlus")
data("phospho_L6_ratio_pe")
data("SPSs")

# load the proteomics dataset
data(PM)

# load pathway annotations
data(Pathways)


# direction pathway analysis in 3-dimensional space. Implemnted as rotating by contrast
# (1) test combined effect of all 3 treatments (stimulation and inhibitions) vs control (basal) 
# on the original direction.
dPA <- directPA(Tc=Tc.gene, direction=c(1,1,1), annotation=Pathways.reactome)

dPA$gst[order(unlist(dPA$gst[,1])),][1:20,]
# rank substrates on the direciton of interest
sort(dPA$gene.pvalues)[1:20]

# (2) test combined effect of all 3 treatments vs controls on direction c(1,-1, 0)
# this rotates Ins by 0 degree, Wmn by 90 degree, and MK by 45 degree.
dPA <- directPA(Tc=Tc.gene, direction=c(1,-1,0), annotation=Pathways.reactome)
dPA$gst[order(unlist(dPA$gst[,1])),][1:20,]
# rank substrates on the direciton of interest
sort(dPA$gene.pvalues)[1:20]

# (3) test combined effect of all 3 perturbations vs controls on direction c(1,-1, 1)
# this rotates Ins by 0 degree, Wmn by 90 degree, and MK by 0 degree.
dPA <- directPA(Tc=Tc.gene, direction=c(1,-1,1), annotation=Pathways.reactome)
dPA$gst[order(unlist(dPA$gst[,1])),][1:20,]
# rank substrates on the direciton of interest
sort(dPA$gene.pvalues)[1:20]


#############################################
#############################################
# Kinase site cluster over time course