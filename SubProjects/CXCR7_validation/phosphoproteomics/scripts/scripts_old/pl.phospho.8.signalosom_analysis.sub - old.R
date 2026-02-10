##############################################################
#########Signalome Analysis
### use mysql database
### get norm abundance
### use parserNorm.pl
### get norm abundance with single phosphosite


### load data

############# plot distributions for all proteins
top.rownames <- c(rownames(top.10[top.10$adj.P.Val<0.05,]),rownames(top.600[top.600$adj.P.Val<0.05,]),rownames(top.1800[top.1800$adj.P.Val<0.05,]))
top.rownames <- c(rownames(top.10[1:60,]),rownames(top.600[1:60,]),rownames(top.1800[1:120,]))
top.norm_intensity <- norm_intensity[top.rownames,]


## problem is that the short sequences often have more often have more psites
## but downstream analyss are optimized for one site only
## therefore via parserNorm.pl the multiple pstes are devided into single psite
norm_abundance2 <- read.table("norm_abundance_ss2.txt", sep="\t", header=TRUE, dec=".")
##or
norm_abundance2 <- top.norm_intensity

norm_abundance2 <- norm_abundance2[!duplicated(norm_abundance2), ]
norm_abundance3 <- norm_abundance2
##just if imported from text
#norm_abundance3 <- norm_abundance2[,c(-1)]
#rownames(norm_abundance3) <- norm_abundance2[,1]


### here we filter the rownames of differential expressed peptides 
### before the phosphosite have to be split to match (use the perl parser for this)
#rownames(top.300[top.300[, "adj.P.Val"] <0.05,])


##############################################################
### determine kinase substrate interaction

mat.reg <- as.matrix(norm_abundance3)

### standardize matrix
mat.std <- PhosR::standardise(mat.reg)

rownames(mat.std) <- sapply(strsplit(rownames(mat.std), ";"), function(x) { 
  gsub(" ", "", paste(toupper(x[2]), x[3], "", sep=";"))
  })
  


### Run PhosR kinase-substrate prediction with the default parameters
data("KinaseMotifs")
#seqs2 <- sapply(strsplit(rownames(mat.reg), ";"), "[[", 4)
seqs <- sapply(strsplit(rownames(mat.reg), ";"), function(x) { gsub(" ", "", x[4]) })

kssMat <- kinaseSubstrateScore(substrate.list = PhosphoSite.human, 
                               mat = mat.std, seqs = seqs,
                               numMotif = 5, numSub = 1, verbose = FALSE)


##############################################################
### export tables 

### table of the scoring matrix

kss1 <- kssMat[["combinedScoreMatrix"]]
kks_act <- kssMat[["ksActivityMatrix"]]

 #   kss1[kks1>0.8,]

rownames(kss1) <- rownames(mat.reg)

write.table(kss1, "kss_late.txt", sep="\t", , row.names=TRUE)



### table of kinase substrates log2fold over time
logfc <- read.table("logfc_intermediate.txt", sep="\t", header=TRUE, dec=".")
logfc2 <- logfc[!duplicated(logfc), ]
logfc3 <- logfc2[,c(-1)]
rownames(logfc3) <- logfc2[,1]


rownames(kss1) <- rownames(mat.reg)

filter <- apply(kss1, 1, function(x) length(x[x>0.7])>=1)
kss1_filt <- kss1[filter,]
kss1_filt2 <- as.data.frame(kss1_filt[kss1_filt[,"ROCK1"]>0.7,"ROCK1"])
colnames(kss1_filt2) <- "ROCK1"


## get the target of the kinases
kss_logfc_time <- logfc3[rownames(logfc3) %in% rownames(kss1_filt2),]

write.table(kss_logfc_time, "ROCK1_target_intermediate_logfc.txt", sep="\t", , row.names=TRUE)

## get the kinases
library(dplyr)
logfc_kss <- logfc3 %>% slice(grep("PRK", row.names(.)))





### table of all kinase subunits
pka_targets_psplus <- PhosphoSite.human[["PRKACA"]]
pkg_targets_psplus <- PhosphoSite.human[["PRKCG"]]

write.table(pka_targets_psplus, "pka_targets_psiteplus.txt", sep="\t", , row.names=TRUE)
write.table(pkg_targets_psplus, "pkg_targets_psiteplus.txt", sep="\t", , row.names=TRUE)

##############################################################
### create signalosom

### PhosR uses the ‘kinaseSubstratePred’ function to synthesise the scores generated from ‘kinaseSubstrateScore’ to predict the kinase-substrate relationships using an adaptive sampling-based positive-unlabeled learning method (Yang et al., 2018)
set.seed(1)
predMat <- kinaseSubstratePred(kssMat, top = 50, cs = 0.8, inclusion = 10, iter = 10, verbose = FALSE, )

filter <- apply(kss1, 1, function(x) length(x[x>0.8])>=1)
kss1 <- kss1[filter,]



colnames(predMat)
write.table(predMat, "predMat.txt", sep="\t", , row.names=TRUE, col.names=TRUE)

### Constructing Signalling Networks (Signalomes)
kinaseOI = c("PRKCA")
signalomesRes <- Signalomes(KSR = kssMat, 
                            predMatrix = predMat, 
                            exprsMat = mat.std, 
                            module_res = 10,
                            KOI = kinaseOI)



### generate palette
my_color_palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))
kinase_all_color <- my_color_palette(ncol(kssMat$combinedScoreMatrix))
names(kinase_all_color) <- colnames(kssMat$combinedScoreMatrix)
kinase_signalome_color <- kinase_all_color[colnames(predMat)]

plotSignalomeMap(signalomes = signalomesRes,
                 color = kinase_signalome_color)


### plot the signalome network that illustrates the connectivity between kinase signalome networks.
plotKinaseNetwork(KSR = kssMat, 
                  predMatrix = predMat, 
                  threshold = 0.95, 
                  color = kinase_all_color)
