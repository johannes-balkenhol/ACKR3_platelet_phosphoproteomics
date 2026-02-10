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
input = list(top.filter.10, top.filter.30, top.filter.60, top.filter.300, top.filter.600,
             top.filter.900, top.filter.1800)


names_input = c("10", "30", "60", "300", "600", "900", "1800")

for (i in 1:length(input)) {
  input[[i]] <- input[[i]][order(row.names(input[[i]])),]
}

Tc <- as.matrix(cbind(input[[1]]$logFC, input[[2]]$logFC, input[[3]]$logFC,
                      input[[4]]$logFC, input[[5]]$logFC, input[[6]]$logFC,
                      input[[7]]$logFC))

rownames(Tc) <- rownames(input[[1]])
colnames(Tc) <- names_input

## define input2 ----
input = list(top.collapse.10, top.collapse.30, top.collapse.60, top.collapse.300, top.collapse.600,
             top.collapse.900, top.collapse.1800)


names_input = c("10", "30", "60", "300", "600", "900", "1800")
for (i in 1:length(input)) {
  input[[i]] <- input[[i]][order(row.names(input[[i]])),]
  
}
names(input) <- names_input
Tc.gene <- as.matrix(cbind(input[[1]]$logFC, input[[2]]$logFC, input[[3]]$logFC,
                           input[[4]]$logFC, input[[5]]$logFC, input[[6]]$logFC,
                           input[[7]]$logFC))

rownames(Tc.gene) <- paste(sapply(strsplit(rownames(input[[1]]), ";"), "[[", 2))
colnames(Tc.gene) <- names_input

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
colors <- c("#FFC7A3", "#FF223B")


class = "DOWN/less"
order = "less"
colors <- c("lightcyan", "royalblue1")

i = 7
#  ENRICHED REACTOME PATHWAYS; MAKE FIGURES; SAVE TO TIFF
for (i in 1:length(input)) {
  geneSet <- names(sort(Tc.gene[,i],  #decreasing False for enrichment of downregulated
                        decreasing = T))[seq(round(nrow(Tc.gene) * 0.05))]
  
  # geneSet <- names(sort(Tc.gene[,1], 
  #                       decreasing = TRUE))[TcP.gene[,1]<0.05]
  #geneSet <- input[[1]]$name[input[[4]]$PValue<0.05]
  
  head(geneSet)
  
  
  path1 <- pathwayOverrepresent(geneSet, annotation=pathways, 
                                universe = rownames(Tc.gene), alter = order)
  path2 <- pathwayRankBasedEnrichment(Tc.gene[,i], 
                                      annotation=pathways, 
                                      alter = order)
  
  lp1 <- -log10(as.numeric(path1[names(pathways),1]))
  lp2 <- -log10(as.numeric(path2[names(pathways),1]))
  
  
  
  ######### Visualiztaion1
  #dev.new()
  # tiff(filename = paste0("../analysis/Reactome_enrichment/", class,  names_input[[i]] , "_plot.tiff"),
  #      width = 5 * 300,
  #      height = 5 * 300,
  #      res = 300,
  #      compression = "lzw")
  # #png(file=paste0("../analysis/reactome_enrichment/", "reactome_top.10",".png"), width = 5, height = 5, units = 'in', res = 1200)
  # #par(mfrow=c(1,1))
  # plot(lp1, lp2, xlab="Overrepresentation (-log10 pvalue)", ylab="Rank-based enrichment (-log10 pvalue)", main="", xlim = c(0, 5), 
  #      ylim = c(0, 2.6), col = "grey", bg = "grey", )
  # 
  # # select highly enriched pathways
  # sel <- which(lp1 > 1.3 & lp2 > 1.3)
  # #textxy(lp1[sel], lp2[sel], gsub("_", " ", gsub("REACTOME_|Mus musculus: ", "", names(pathways)))[sel], cex= 0.8,srt=0)
  # label <- gsub("_", " ", gsub("REACTOME_|Homo sapiens: ", "", names(pathways)))[sel]
  # label <- word_wrap(label, 45, linesep = NULL)
  # addTextLabels(lp1[sel], lp2[sel], label, cex.label=0.9, col.label="black", lty=2, col.line=rgb(0,0,0, 0.5),keepLabelsInside = TRUE)
  # dev.off()
  
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
  #path4 <- path4[order(path4$pvalue),]
  path4 <- path4[path4$pvalue<0.05,]
  path4 <- path4[order(path4$ratio, decreasing = T),]
  #write.table(path4, file = paste0("../analysis/Reactome_enrichment/", class, names_input[[i]] ,"_pathways.txt"), 
  #            sep = '\t', row.names = FALSE)
  #path4 <- path4[1:10,] #AutoTop10
  path4 <- path4[c(1,3,5:8,10,15,17,35),] # ManualTop10 - UP10:c(1:3, 5:11), DOWN10:c(1:3, 5,7,9,15,19,23,24),  UP600:c(1:2, 5:11, 14), DOWN600:c(1,2,5:9,12,18,22),  UP1800:c(1,3,5:8,10,15,17,35), DOWN1800:c(1:5,8,12,14,17,18)
  #path4[10, "pathway"] <- "Regulation.of.IGF.transport.and.uptake.by.IGFBPs"
  path4$pvalue <- round(path4$pvalue, digits=4)
  path4 <- path4[order(path4$pvalue),]
  path4$pathway <- gsub("\\.", " ", path4$pathway)
  
  #dev.new()
  ggp <- ggplot(path4[10:1,], aes(x=1:10,y=ratio,fill=pvalue)) + 
    coord_flip() +
    geom_bar(stat = "identity", colour="black") + 
    scale_fill_gradient(low = colors[1], high = colors[2], 
                        limits = c(min(path4[1:10,]$pvalue), max(path4[1:10,]$pvalue)), 
                        name = "p-value", 
                        guide = guide_colorbar(barwidth = , 
                                               barheight = 10, 
                                               title.position = "top", 
                                               title.hjust = 0.5,
                                               reverse = TRUE))  +
    geom_text(aes(label = pathway, y = 0.005),
              vjust = 0.5, colour = "black",
              position = position_dodge(0.1), size =8, hjust = 'left') +
    labs(title = paste0(names_input[i]," sec"), size = 35, x = "pathways", y = "gene ratio") +
    theme_cowplot() +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 20),
          axis.text.x=element_text(size=30),
          title = element_text(size = 30),
          axis.title.x = element_text(size =35),
          axis.title.y = element_text(size = 35)
    )
  tiff(filename = paste0("../analysis/Reactome_enrichment/", class, names_input[[i]] ,"_barplot2.tiff"),
       width = 12 * 300, 
       height = 8 * 300,
       res = 300,
       compression = "lzw")
  
  print(ggp)
  dev.off()
  
  ggsave(file = paste0("../analysis/Reactome_enrichment/", class, names_input[[i]] ,"_barplot2.pdf"),
       plot = ggp,
       width = 12, 
       height = 8,
       dpi = 300,
       device = "pdf")
  
  assign(paste0("path4.", names_input[[i]]), path4)
  
}










#### extraxt certain pathways
rea1 <- pathways["Homo sapiens: Signaling by ROBO receptors"]
path2["Homo sapiens: Signaling by ROBO receptors"]
rea1 <- rea1[rea1 %in% geneSet]
write.table(rea1, "../analysis/Reactome_enrichment/rea1.txt", sep="\t",  row.names=FALSE, quote=FALSE)

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

tiff(filename = paste0("../analysis/DPA/DPA.tiff"),
     width = 14 * 300, 
     height = 8 * 300,
     res = 300,
     compression = "lzw")

par(mfrow=c(3,7))
dpa1 <- directPA(Tc1[,c(1,2)], direction=0, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dpa2 <- directPA(Tc1[,c(1,3)], direction=0, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dpa3 <- directPA(Tc1[,c(1,4)], direction=0, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dpa4 <- directPA(Tc1[,c(1,5)], direction=0, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dpa5 <- directPA(Tc1[,c(1,6)], direction=0, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dpa6 <- directPA(Tc1[,c(1,7)], direction=0, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dpa7 <- directPA(Tc1[,c(2,3)], direction=0, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dpa8 <- directPA(Tc1[,c(2,4)], direction=0, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dpa9 <- directPA(Tc1[,c(2,5)], direction=0, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dpa10 <- directPA(Tc1[,c(2,6)], direction=0, 
                  annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                  main="Direction pathway analysis")
dpa11 <- directPA(Tc1[,c(2,7)], direction=0, 
                  annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                  main="Direction pathway analysis")
dpa12 <- directPA(Tc1[,c(3,4)], direction=0, 
                  annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                  main="Direction pathway analysis")
dpa13 <- directPA(Tc1[,c(3,5)], direction=0, 
                  annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                  main="Direction pathway analysis")
dpa14 <- directPA(Tc1[,c(3,6)], direction=0, 
                  annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                  main="Direction pathway analysis")
dpa15 <- directPA(Tc1[,c(3,7)], direction=0, 
                  annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                  main="Direction pathway analysis")
dpa16 <- directPA(Tc1[,c(4,5)], direction=0, 
                  annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                  main="Direction pathway analysis")
dpa17 <- directPA(Tc1[,c(4,6)], direction=0, 
                  annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                  main="Direction pathway analysis")
dpa18 <- directPA(Tc1[,c(4,7)], direction=0, 
                  annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                  main="Direction pathway analysis")
dpa19 <- directPA(Tc1[,c(5,6)], direction=0, 
                  annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                  main="Direction pathway analysis")
dpa20 <- directPA(Tc1[,c(5,7)], direction=0, 
                  annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                  main="Direction pathway analysis")
dpa21 <- directPA(Tc1[,c(6,7)], direction=0, 
                  annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                  main="Direction pathway analysis")
dev.off()

tiff(filename = paste0("../analysis/DPA/DPA_short.tiff"),
     width = 15 * 300, 
     height = 5 * 300,
     res = 300,
     compression = "lzw")

par(mfrow=c(1,3))
dpa4 <- directPA(Tc1[,c(1,5)], direction=0, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dpa6 <- directPA(Tc1[,c(1,7)], direction=0, 
                 annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                 main="Direction pathway analysis")
dpa20 <- directPA(Tc1[,c(5,7)], direction=0, 
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


tiff(filename = paste0("../analysis/DPA/Kinase_pert.tiff"),
     width = 24 * 300, 
     height = 12 * 300,
     res = 300,
     compression = "lzw")

par(mfrow=c(3,7))

z1 <- perturbPlot2d(Tc1[,c(1,2)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                    main="Kinase perturbation analysis")
z2 <- perturbPlot2d(Tc1[,c(1,3)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                    main="Kinase perturbation analysis")
z3 <- perturbPlot2d(Tc1[,c(1,4)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                    main="Kinase perturbation analysis")
z4 <- perturbPlot2d(Tc1[,c(1,5)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                    main="Kinase perturbation analysis")
z5 <- perturbPlot2d(Tc1[,c(1,6)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                    main="Kinase perturbation analysis")
z6 <- perturbPlot2d(Tc1[,c(1,7)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                    main="Kinase perturbation analysis")
z7 <- perturbPlot2d(Tc1[,c(2,3)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                    main="Kinase perturbation analysis")
z8 <- perturbPlot2d(Tc1[,c(2,4)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                    main="Kinase perturbation analysis")
z9 <- perturbPlot2d(Tc1[,c(2,5)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                    main="Kinase perturbation analysis")
z10 <- perturbPlot2d(Tc1[,c(2,6)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
                     annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                     main="Kinase perturbation analysis")
z11 <- perturbPlot2d(Tc1[,c(2,7)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
                     annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                     main="Kinase perturbation analysis")
z12 <- perturbPlot2d(Tc1[,c(3,4)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
                     annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                     main="Kinase perturbation analysis")
z13 <- perturbPlot2d(Tc1[,c(3,5)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
                     annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                     main="Kinase perturbation analysis")
z14 <- perturbPlot2d(Tc1[,c(3,6)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
                     annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                     main="Kinase perturbation analysis")
z15 <- perturbPlot2d(Tc1[,c(3,7)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
                     annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                     main="Kinase perturbation analysis")
z16 <- perturbPlot2d(Tc1[,c(4,5)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
                     annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                     main="Kinase perturbation analysis")
z17 <- perturbPlot2d(Tc1[,c(4,6)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
                     annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                     main="Kinase perturbation analysis")
z18 <- perturbPlot2d(Tc1[,c(4,7)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
                     annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                     main="Kinase perturbation analysis")
z19 <- perturbPlot2d(Tc1[,c(5,6)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
                     annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                     main="Kinase perturbation analysis")
z20 <- perturbPlot2d(Tc1[,c(5,7)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
                     annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                     main="Kinase perturbation analysis")
z21 <- perturbPlot2d(Tc1[,c(6,7)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
                     annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                     main="Kinase perturbation analysis")

dev.off()

tiff(filename = paste0("../analysis/DPA/Kinase_pert_short.tiff"),
     width = 15 * 300, 
     height = 5 * 300,
     res = 300,
     compression = "lzw")

par(mfrow=c(1,3))
z4 <- perturbPlot2d(Tc1[,c(1,5)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                    main="Kinase perturbation analysis")
z6 <- perturbPlot2d(Tc1[,c(1,7)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
                    annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
                    main="Kinase perturbation analysis")
z20 <- perturbPlot2d(Tc1[,c(5,7)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
                     annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
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
dPA <- directPA(Tc=Tc1[,c(1,5,7)], direction=c(1,1,1), annotation=Pathways.reactome)
perturbPlot3d(Tc1[,c(1,5,7)], cex=1, xlim=c(-2, 4), ylim=c(-2, 4),  
              annotation=lapply(PhosphoSite.human, function(x){gsub(";[STY]", ";", x)}), 
              main="Kinase perturbation analysis")
bda <- directExplorer2d(Tc=Tc[c(1,7)], annotation=PhosphoSite.human)
bda$gene.tab[order(bda$gene.tab[,"*-"]),][1:10,]


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