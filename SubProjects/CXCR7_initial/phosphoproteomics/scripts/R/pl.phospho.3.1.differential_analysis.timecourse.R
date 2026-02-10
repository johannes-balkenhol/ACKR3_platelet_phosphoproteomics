## https://nellev.github.io/2019timecourse-rnaseq-pipeline/manuscript.pdf
## https://nellev.github.io/moanin/articles/documentation.html

library(BiocManager)
# Need to use the development version for now.


BiocManager::install(version="devel")
BiocManager::valid()
BiocManager::install("timecoursedata")
BiocManager::install("moanin")
The following additional packages are needed for this workflow:
# From Github
library(moanin)
library(timecoursedata)
# From CRAN
library(NMF)
library(ggfortify)
# From Bioconductor
library(topGO)
library(biomaRt)
library(KEGGprofile)
library(BiocWorkflowTools) #Needed for compiling the .Rmd script

BiocManager::install(
c("NMF", "ggfortify", "topGO", "biomaRt",
"KEGGprofile", "BiocWorkflowTools", "timecoursedata"))



