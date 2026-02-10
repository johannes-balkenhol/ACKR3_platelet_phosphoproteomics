install.packages("ClueR")

suppressPackageStartupMessages({
  library(parallel)
  library(ggplot2)
  library(ClueR)
  library(reactome.db)
  library(org.Mm.eg.db)
  library(annotate)
  library(PhosR)
})



############################ 
### example from website phosr
# take grouping information
# data("phospho_liverInsTC_RUV_pe")
data("phospho.liver.Ins.TC.ratio.RUV.pe")
ppe <- phospho.liver.Ins.TC.ratio.RUV.pe
ppe



grps <- sapply(strsplit(colnames(ppe), "_"), 
            function(x)x[3])

# select differentially phosphorylated sites
sites.p <- matANOVA(SummarizedExperiment::assay(ppe, "Quantification"), 
                    grps)
ppm <- meanAbundance(SummarizedExperiment::assay(ppe, "Quantification"), grps)
sel <- which((sites.p < 0.05) & (rowSums(abs(ppm) > 1) != 0))
ppm_filtered <- ppm[sel,]

# summarise phosphosites information into gene level
ppm_gene <- phosCollapse(ppm_filtered, 
            gsub(";.+", "", rownames(ppm_filtered)), 
                stat = apply(abs(ppm_filtered), 1, max), by = "max")



###################################
### select all differential sites

input = list(top.filter.10, top.filter.600, top.filter.1800, 
top.filter.10.dmso.vs.0s, top.filter.600.dmso.vs.0s, top.filter.1800.dmso.vs.0s,
top.filter.10.cxcr7.vs.0s, top.filter.600.cxcr7.vs.0s, top.filter.1800.cxcr7.vs.0s)

names_input = c("10", "600", "1800", 
"10.dmso.vs.0s", "600.dmso.vs.0s", "1800.dmso.vs.0s",
"10.cxcr7.vs.0s", "600.cxcr7.vs.0s", "1800.cxcr7.vs.0s")

#input_tables_names <- list("top.10", "top.600", "top.1800")

#top.rownames <- c(rownames(top.filter.10),rownames(top.600[1:20,]),rownames(top.1800[1:20,]))
#top.norm_intensity <- norm_intensity[top.rownames,]
#rownames(top.norm_intensity) <- sapply(strsplit(rownames(top.norm_intensity), ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})

cr = 1

###################################
### select top differential sites

for (i in 1:3) {
    #top.names <- append(top.names, rownames(input[[i]][input[[i]]$PValue<0.05,]))
    #assign(paste0(".", names_input[[i]]), )
    if (i == 1 & i <= 2) { 
      top.names <- rownames(input[[i]][input[[i]]$PValue<0.05,])
  } else {
    print(as.list(test[[i]]))
      top.names <- append(top.names, rownames(input[[i]][input[[i]]$PValue<0.05,]))
  }
}


norm_intensity2 <- as.data.frame(norm_intensity)
#norm_intensity2 <- meanAbundance(norm_intensity2, grps)
top.norm_intensity <- norm_intensity2[top.names,]


ppm_gene <- phosCollapse(top.norm_intensity, 
            sapply(strsplit(rownames(top.norm_intensity), ";"), "[[", 2), 
                stat = apply(abs(top.norm_intensity), 1, max), by = "max")


#sapply(strsplit(rownames(Tc), ";"), "[[", 2)

################################
### pathway cluster enrichment

pathways = as.list(reactomePATHID2EXTID)
path_names = as.list(reactomePATHID2NAME)
name_id = match(names(pathways), names(path_names))
names(pathways) = unlist(path_names)[name_id]

pathways = pathways[which(grepl("Homo sapiens", names(pathways), ignore.case = TRUE))]

pathways = lapply(pathways, function(path) {
    gene_name = unname(getSYMBOL(path, data = "org.Hs.eg"))
    toupper(unique(gene_name))
})# perform ClueR to identify optimal number of clusters
 names(pathways) <- gsub("_", " ", gsub("REACTOME_|Homo sapiens: ", "", names(pathways)))

#pathways = as.list(reactomePATHID2EXTID)
#pathways = pathways[which(grepl("Homo sapiens", names(pathways), ignore.case = TRUE))]
#pathways = pathways[which(grepl("R-HSA", names(pathways), ignore.case = TRUE))]
#pathways = lapply(pathways, function(path) {
#    gene_name = unname(getSYMBOL(path, data = "org.Hs.eg"))
#    toupper(unique(gene_name))
#})

RNGkind("L'Ecuyer-CMRG")
set.seed(123)
c1 <- runClue(ppm_gene, annotation=pathways, 
            kRange = seq(2,10), rep = 5, effectiveSize = c(5, 100), 
            pvalueCutoff = 0.05, alpha = 0.5)

# Visualise the evaluation results
data <- data.frame(Success=as.numeric(c1$evlMat), Freq=rep(seq(2,10), each=5))
myplot <- ggplot(data, aes(x=Freq, y=Success)) + 
    geom_boxplot(aes(x = factor(Freq), fill="gray")) +
    stat_smooth(method="loess", colour="red", size=3, span = 0.5) +
    xlab("# of cluster") + 
    ylab("Enrichment score") + 
    theme_classic()

set.seed(123)
best <- clustOptimal(c1, rep=10, mfrow=c(2, 3), visualize = TRUE)

best$enrichList

set.seed(1)
best <- clustOptimal(c3, rep=10, mfrow=c(2, 3), visualize = TRUE)

best$enrichList


######################
### psite cluster enrichment
RNGkind("L'Ecuyer-CMRG")
set.seed(1)

PhosphoSite.human2 = mapply(function(kinase) {
  gsub("(.*)(;[A-Z])([0-9]+;)", "\\1;\\3", kinase)
}, PhosphoSite.human)

# perform ClueR to identify optimal number of clusters
ppm_filtered <- top.norm_intensity

ppm_filtered <- phosCollapse(ppm_filtered, 
                paste(sapply(strsplit(rownames(ppm_filtered), ";"), "[[", 2), 
                      sapply(strsplit(rownames(ppm_filtered), ";"), "[[", 3), "",
                sep = ";"),
                stat = apply(abs(ppm_filtered), 1, mean), by = "mean")

#rownames(ppm_filtered) <- paste(sapply(strsplit(rownames(ppm_filtered), ";"), "[[", 2), 
#                          sapply(strsplit(rownames(ppm_filtered), ";"), "[[", 3), sep = ";")
rownames(ppm_filtered) <- lapply(rownames(ppm_filtered), function(x){gsub(";[STY]", ";", x)})


c3 <- runClue(ppm_filtered, annotation=PhosphoSite.human2, 
kRange = 2:10, rep = 5, effectiveSize = c(5, 100), pvalueCutoff = 0.05, alpha = 0.5)

# Visualise the evaluation results
data <- data.frame(Success=as.numeric(c3$evlMat), Freq=rep(2:10, each=5))

myplot <- ggplot(data, aes(x=Freq, y=Success)) + geom_boxplot(aes(x = factor(Freq), fill="gray"))+
  stat_smooth(method="loess", colour="red", size=3, span = 0.5) + xlab("# of cluster")+ ylab("Enrichment score")+theme_classic()
myplot

set.seed(1)
best <- clustOptimal(c3, rep=10, mfrow=c(2, 3), visualize = TRUE)
