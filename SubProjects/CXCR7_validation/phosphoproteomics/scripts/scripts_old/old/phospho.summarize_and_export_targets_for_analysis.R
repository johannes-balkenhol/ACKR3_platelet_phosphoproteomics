################################################################
################################################################
############ filter Mitas Interesting targets, annotate
############  and export for mysql
############
############ use script: analyze_phophoproteom_full.r
############ use R dataframe top and top.pka
library("dplyr")
########## loaddata

### go to script folder
setwd("C:/Users/Ozge/Nextcloud/Phosphoproteom/phosphoproteom validation/scripts") ##ÖO folder structure should fit.
setwd("D:/Eigene Datein/Phosphoproteom/phosphoproteom validation/scripts/")
setwd("F:/Masterarbeit_BioWis/Phosphoproteomics/Phosphoproteom validation/scripts/")



#####  test
UniprotID<- sapply(strsplit(rownames(ppe), ";"), "[[", 1)
GeneSymbol <- sapply(strsplit(rownames(ppe), ";"), "[[", 2)



##### summarize all targets
##### annotate and export for mysql

top.all <- top

colnames(top.all) <- c("ID",
"uniprot_id_0000","logfc_0010","p.value_0010","adj.p.value_0010","symbol_0010",
"uniprot_id_10","logfc_0030","p.value_0030","adj.p.value_0030","symbol_0030",
"uniprot_id_30","logfc_0060","p.value_0060","adj.p.value_0060","symbol_0060",
"uniprot_id_300","logfc_0300","p.value_0300","adj.p.value_0300","symbol_0300",
"uniprot_id_600","logfc_0600","p.value_0600","adj.p.value_0600","symbol_0600",
"uniprot_id_900","logfc_0900","p.value_0900","adj.p.value_0900","symbol_0900",
"uniprot_id_1800","logfc_1800","p.value_1800","adj.p.value_1800","symbol_1800")

top.all <- top.all %>% dplyr::select(sort(names(top.all)))

##### set a column for specified phosphosite
##### or export to mysql


top.all$Residue <- gsub("[0-9]","", sapply(strsplit(rownames(top.all), ";"), "[[", 3))
top.all$Site <- gsub("[A-Z]","", sapply(strsplit(rownames(top.all), ";"), "[[", 3))
top.all$Sequence <- sapply(strsplit(rownames(top.all), ";"), "[[", 4)


## read tables locally or
## RS: connect to mysql server and load the table from the mysql database
write.table(top.all, "top.all.csv", sep="\t", , row.names=FALSE)


 
top.all.sub <- read.table("top.all.csv", sep="\t", header=TRUE, dec=",")

top.all.sub$peptide_id <- sapply(strsplit(top.all.sub$ID, ";"), "[[", 5)

write.table(top.all.sub, "top.all.sub.csv", sep="\t", , row.names=FALSE)


################################################################
##### export norm abundance
##### export raw abundance
raw_intensity <- SummarizedExperiment::assay(ppe0,"Quantification")
scaled_intensity <-SummarizedExperiment::assay(ppe,"scaled")
norm_intensity <-SummarizedExperiment::assay(ppe,"normalised")

## read tables locally or
## RS: connect to mysql server and load the table from the mysql database
write.table(norm_intensity, "norm_intensity.csv", sep="\t", , row.names=TRUE)
write.table(raw_intensity, "raw_intensity.csv", sep="\t", , row.names=TRUE)

