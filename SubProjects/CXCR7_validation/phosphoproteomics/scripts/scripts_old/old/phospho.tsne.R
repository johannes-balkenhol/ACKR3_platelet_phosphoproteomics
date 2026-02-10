################################################################
################################################################
############ tSNE

##JB no cool yet
### preapre abundance data AND groups 
quantification <- SummarizedExperiment::assay(ppe0,"Quantification")
scaled <- SummarizedExperiment::assay(ppe,"scaled")
normalised <- SummarizedExperiment::assay(ppe,"normalised")


x_tt <- as.factor(grps)

### version 2
BiocManager::install(version='devel')
BiocManager::install("M3C")

# example 
BiocManager::install("M3C")
library(M3C)
BiocManager::install("scde")
library(scde)


# mita data
tsne(normalised,labels=x_tt)
