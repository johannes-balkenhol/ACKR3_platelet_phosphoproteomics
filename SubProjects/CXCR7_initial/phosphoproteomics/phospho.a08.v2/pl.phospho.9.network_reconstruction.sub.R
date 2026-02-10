################################################################
##### PSIQUIC network reconstruction 
BiocManager::install("PSICQUIC")
library("PSICQUIC")
library(stringr)


##### set workign directory
setwd("C:/Users/Ozge/Nextcloud/Phosphoproteom/phosphoproteom validation/scripts") ##ÖO folder structure should fit.
setwd("D:/Eigene Datein/Phosphoproteom/phosphoproteom validation/scripts/")
setwd("F:/Masterarbeit_BioWis/Phosphoproteomics/Phosphoproteom validation/scripts/")


################################################################
##### get intesity data

UniprotID <- sapply(strsplit(rownames(ppe_norm), ";"), "[[", 1)
GeneSymbol <- sapply(strsplit(rownames(ppe_norm), ";"), "[[", 2)

##### collpase with collpase funiton
norm_intensity.collapse <- phosCollapse(norm_intensity, id=sapply(strsplit(rownames(norm_intensity), ";"), "[[", 1), 
                        stat=apply(abs(norm_intensity), 1, max), by = "max")
						
x_tt <- as.factor(grps)

#################################################################
##### psiquic connect
psicquic <- PSICQUIC()
providers(psicquic)


#### example
interactions(object,
                                  id=NA,
                                  species=NA,
                                  speciesExclusive=TRUE,
                                  type=NA,
                                  provider=NA,
                                  detectionMethod=NA,
                                  publicationID=NA,
                                  quiet=TRUE)

tbl <- interactions(psicquic, id=c("TP53", "MYC"), species="9606")
tbl[, c("provider", "type", "detectionMethod")]
tbl[grep("affinity", tbl$detectionMethod),
+ c("type", "publicationID", "firstAuthor", "confidenceScore", "provider")]

tbl.myc <- interactions(psicquic, "MYC", species="9606", publicationID="21150319")

table(tbl.myc$provider)

table(tbl.myc$confidenceScore)

idMapper <- IDMapper("9606")

tbl.myc <- addGeneInfo(idMapper,tbl.myc)

print(head(tbl.myc$A.name))

print(head(tbl.myc$B.name))

tbl.3 <- interactions(psicquic, id=c("ALK", "JAK3", "SHC3"),
+ species="9606", quiet=TRUE)

tbl.3g <- addGeneInfo(idMapper, tbl.3)

tbl.3gd <- with(tbl.3g, as.data.frame(table(detectionMethod, type, A.name, B.name, provider)))
print(tbl.3gd <- subset(tbl.3gd, Freq > 0))


################################################################
##### psiquic: get interactions of proteins/psites in individual top tables
UniprotID <- sapply(strsplit(rownames(top.10.sign), ";"), "[[", 1)
interactions.top.10 <- interactions(psicquic, id=UniprotID, species="9606", provider="IntAct")

#### filter the interactions
## read detection method table in /anaylsis/network_reconstruction/psiquic_folder
## to pull the interesting detection method

table(interactions.top.10$provider)
table(interactions.top.10$detectionMethod)
table(interactions.top.10$type)

detectionmethod <- c("psi-mi:MI:0407(direct interaction)", "psi-mi:MI:0915(physical association)", 
"psi-mi:MI:0407(direct interaction)")
interactions.top.10.a <- interactions.top.10[which(interactions.top.10$type %in% detectionmethod),]
interactions.top.10.a <- interactions.top.10.a[grep("uniprotkb", interactions.top.10.a$A),]
interactions.top.10.a <- interactions.top.10.a[grep("uniprotkb", interactions.top.10.a$B),]
interactions.top.10.a$A <- gsub("uniprotkb:", "", interactions.top.10.a$A)
interactions.top.10.a$B <- gsub("uniprotkb:", "", interactions.top.10.a$B)
interactions.top.10.a$gene_name1 <- sub(":(\S+)(gene name)", "\\1", interactions.top.10.a$aliasA,perl=TRUE)

interactions.top.10.a$gene_name1 <- str_extract(interactions.top.10.a$aliasA, "\\|uniprotkb:\\S+(gene name)")
interactions.top.10.a$gene_name1 <- gsub("\\|uniprotkb:(.+)\\(gene name","\\1",interactions.top.10.a$gene_name1)

interactions.top.10.a$gene_name2 <- str_extract(interactions.top.10.a$aliasB, "\\|uniprotkb:\\S+(gene name)")
interactions.top.10.a$gene_name2 <- gsub("\\|uniprotkb:(.+)\\(gene name","\\1",interactions.top.10.a$gene_name2)

interactions.top.10.a$score <- gsub(".*intact-miscore:(\\d+)","\\1",interactions.top.10.a$confidenceScore,,perl=TRUE)

interactions.top.10.b <- interactions.top.10.a[,c(1,2,17,18,19,12,7,8,9,13)]

colnames(interactions.top.10.b)[1] = "uniprot_id1"
colnames(interactions.top.10.b)[2] = "uniprot_id2"

interactions.top.10.c <- interactions.top.10.b[interactions.top.10.b$score > 0.5,]

################################################################
##### psiquic: get interactions of proteins/psites in top tables
UniprotID <- unique(sapply(strsplit(top.rownames, ";"), "[[", 1))
interactions.top.all <- interactions(psicquic, id=UniprotID, species="9606", provider="IntAct")


table(interactions.top.all$provider)
table(interactions.top.all$detectionMethod)
table(interactions.top.all$type)

detectionmethod <- c("psi-mi:MI:0407(direct interaction)", "psi-mi:MI:0915(physical association)", 
"psi-mi:MI:0407(direct interaction)")
interactions.top.all.a <- interactions.top.all[which(interactions.top.all$type %in% detectionmethod),]
interactions.top.all.a <- interactions.top.all.a[grep("uniprotkb", interactions.top.all.a$A),]
interactions.top.all.a <- interactions.top.all.a[grep("uniprotkb", interactions.top.all.a$B),]
interactions.top.all.a$A <- gsub("uniprotkb:", "", interactions.top.all.a$A)
interactions.top.all.a$B <- gsub("uniprotkb:", "", interactions.top.all.a$B)

interactions.top.all.a$gene_name1 <- str_extract(interactions.top.all.a$aliasA, "\\|uniprotkb:\\S+(gene name)")
interactions.top.all.a$gene_name1 <- gsub("\\|uniprotkb:(.+)\\(gene name","\\1",interactions.top.all.a$gene_name1)

interactions.top.all.a$gene_name2 <- str_extract(interactions.top.all.a$aliasB, "\\|uniprotkb:\\S+(gene name)")
interactions.top.all.a$gene_name2 <- gsub("\\|uniprotkb:(.+)\\(gene name","\\1",interactions.top.all.a$gene_name2)

interactions.top.all.a$score <- gsub(".*intact-miscore:(\\d+)","\\1",interactions.top.all.a$confidenceScore,,perl=TRUE)

interactions.top.all.b <- interactions.top.all.a[,c(1,2,17,18,19,12,7,8,9,13)]

colnames(interactions.top.all.b)[1] = "uniprot_id1"
colnames(interactions.top.all.b)[2] = "uniprot_id2"

interactions.top.all.c <- interactions.top.all.b[interactions.top.all.b$score > 0.5,]


sif <- data.frame(interactions.top.all.c$gene_name1, "interacts_with", interactions.top.all.c$gene_name2)
colnames(sif) <- c("gene_name1", "type", "gene_name2")

write.table(sif, "../analysis/cytsocape/sif.txt", sep="\t", row.names=FALSE)


################################################################
##### psiquic: get interactions of all proteins/psites

##### load packages
library(readr)
library(SparkR)
library(stringr)


##### download and import the interactions
UniprotID <- unique(sapply(strsplit(rownames(norm_intensity), ";"), "[[", 1))

## query psicquic from R (but its lsow somehow
interactions.all <- interactions(psicquic, id=UniprotID, species="9606", provider="IntAct")

## or download directly from intact
## https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/
## and import into dataframe
#interactions.all <- read.table("D:/Eigene Datein/Interactions/head_intact.txt", sep="\t", header=TRUE, dec=",")
interactions.all <- read.delim("D:/Eigene Datein/Interactions/intact_negative.txt", sep="\t", header=TRUE, skip = 0)
## file to big

f <- function(x, pos) subset(x, pattern = "taxid:9606")
interactions.all.test <- read_delim_chunked("D:/Eigene Datein/Interactions/intact_negative.txt", delim="\t", DataFrameCallback$new(f), chunk_size = 1000000)

f <- function(x, pos) x[,c(1,2,5,6,7,8,9,10,11,12,13,15,21,22)]
interactions.all.test <- read_delim_chunked("D:/Eigene Datein/Interactions/intact_negative.txt", delim="\t", DataFrameCallback$new(f), chunk_size = 1000000)

f <- function(x, pos) subset(x, select = grepl("9606", x$`Host organism(s)`))
interactions.all.test <- read_delim_chunked("D:/Eigene Datein/Interactions/intact_negative.txt", delim="\t", DataFrameCallback$new(f), chunk_size = 1000000)

## after long term trying this seems to work
f <- function(x, pos) x[,c(1,2,5,6,7,8,9,10,11,12,13,15,21,22)]
interactions.all <- read_delim_chunked("D:/Eigene Datein/Interactions/intact.txt", delim="\t", DataFrameCallback$new(f), chunk_size = 100000)
## memory limit 
## https://www.rdocumentation.org/packages/utils/versions/3.4.1/topics/memory.size
## memory.size(max = FALSE)
## memory.limit(size = NA)
##  gc()
## read with sparkR
## ## install.packages("SparkR")

## https://docs.databricks.com/spark/latest/sparkr/overview.html
## install.packages("config")
## https://cran.r-project.org/web/packages/config/vignettes/introduction.html
##

##### parse and pre-filter the imported interactions
colnames(interactions.all)
colnames(interactions.all) <- c("A","B","aliasA","aliasB","detectionMethod","firstAuthor","publicationID","taxonA","taxonB","type","sourceDatabases","confidenceScore","typeA","typeB")

#table(interactions.all$provider)
table(interactions.all$detectionMethod)
table(interactions.all$type)

detectionmethod <- c("psi-mi:MI:0407(direct interaction)", "psi-mi:MI:0915(physical association)", 
"psi-mi:MI:0407(direct interaction)")
interactions.all.a <- interactions.all
interactions.all.a$type <- gsub("\"","",interactions.all.a$type)
interactions.all.a <- interactions.all.a[which(interactions.all.a$type %in% detectionmethod),]
interactions.all.a <- interactions.all.a[grep("uniprotkb", interactions.all.a$A),]
interactions.all.a <- interactions.all.a[grep("uniprotkb", interactions.all.a$B),]
interactions.all.a$A <- gsub("uniprotkb:", "", interactions.all.a$A)
interactions.all.a$B <- gsub("uniprotkb:", "", interactions.all.a$B)

interactions.all.a$gene_name1 <- str_extract(interactions.all.a$aliasA, "\\|uniprotkb:\\S+(gene name)")
interactions.all.a$gene_name1 <- gsub("\\|uniprotkb:(.+)\\(gene name","\\1",interactions.all.a$gene_name1)

interactions.all.a$gene_name2 <- str_extract(interactions.all.a$aliasB, "\\|uniprotkb:\\S+(gene name)")
interactions.all.a$gene_name2 <- gsub("\\|uniprotkb:(.+)\\(gene name","\\1",interactions.all.a$gene_name2)

interactions.all.a$score <- gsub(".*intact-miscore:(\\d+)","\\1",interactions.all.a$confidenceScore,,perl=TRUE)



> colnames(interactions.all.a)
 [1] "A"               "B"               "aliasA"          "aliasB"         
 [5] "detectionMethod" "firstAuthor"     "publicationID"   "taxonA"         
 [9] "taxonB"          "type"            "sourceDatabases" "confidenceScore"
[13] "typeA"           "typeB"           "gene_name1"      "gene_name2"     
[17] "score" 

#####
> colnames(interactions.top.all.b)
 [1] "uniprot_id1"     "uniprot_id2"     "gene_name1"      "gene_name2"     
 [5] "score"           "type"            "detectionMethod" "firstAuthor"    
 [9] "publicationID"   "sourceDatabases"

interactions.all.b <- interactions.all.a[,c(1,2,15,16,17,10,5,6,7,11)]
#interactions.all.b <- interactions.all.a[,c(interactions.all.a$A,interactions.all.a$B,interactions.all.a$gene_name1,
#interactions.all.a$gene_name2,interactions.all.a$score,interactions.all.a$type,interactions.all.a$detectionMethod,
#interactions.all.a$firstAuthor,interactions.all.a$publicationID,interactions.all.a$sourceDatabase)]

colnames(interactions.all.b)[1] = "uniprot_id1"
colnames(interactions.all.b)[2] = "uniprot_id2"

write.table(interactions.all.bb, "../analysis/network_reconstruction/interactions.all.b.txt", sep="\t", row.names=FALSE, , quote = FALSE, , eol = "\n")
##

interactions.all.b <- read.delim("../analysis/network_reconstruction/interactions.all.b.txt", sep="\t", header=TRUE, skip = 0)

interactions.all.c <- interactions.all.b[interactions.all.b$score > 0.3,]
dim(interactions.all.c)

##### now filter for the PL proteins
#"P25106" cxcr7
#"P61073" cxcr4
interactions.all.c[interactions.all.c$uniprot_id1 == "P25106",]

UniprotID <- c(unique(sapply(strsplit(rownames(norm_intensity), ";"), "[[", 1)),"P25106","P61073")
interactions.all.d <- unique(interactions.all.c[which(interactions.all.c$uniprot_id1 %in% UniprotID),])
interactions.all.d <- unique(interactions.all.d[which(interactions.all.d$uniprot_id2 %in% UniprotID),])
dim(interactions.all.d)
interactions.all.d[interactions.all.d$uniprot_id1 == "P25106",]
#P49407 ARRB1
#P32121 ARRB2
#P61073 CXCR4
top.all[top.all$uniprot == "P61073",]
top.all[top.all$uniprot == "P49407",]
top.all[top.all$uniprot == "P32121",]

##### export to cytoscpae
## network file
sif <- data.frame(interactions.all.d$gene_name1, "interacts_with", interactions.all.d$gene_name2)
colnames(sif) <- c("gene_name1", "type", "gene_name2")
write.table(sif, "../analysis/cytsocape/sif.txt", sep="\t", row.names=FALSE, quote = FALSE)

## edge infromatin file
edge <- data.frame(paste0(interactions.all.d$gene_name1, " (", "interacts_with", ") ", interactions.all.d$gene_name2),interactions.all.d$score)
write.table(edge, "../analysis/cytsocape/edge.txt", sep="\t", row.names=FALSE, quote = FALSE)

## node information file
##### 2. get top tables data and abundance data


#### collapse sites
### with phoscollpase with seem to be faily
top.all.logfc <- top.all[,c(3,6,9)]
library(PhosR)
top.all.logfc.collapse <- phosCollapse(top.all.logfc, id=sapply(strsplit(rownames(norm_intensity), ";"), "[[", 1), 
                        stat=apply(abs(top.all.logfc), 1, max), by = "max")

top.all.pval <- top.all[,c(5,8,11)]
library(PhosR)
top.all.pval.collapse <- phosCollapse(top.all.pval, id=sapply(strsplit(rownames(norm_intensity), ";"), "[[", 1), 
                        stat=apply(abs(top.all.pval), 1, min), by = "min")


### or collpase self made
top.all.logfc.collapse <- ddply(top.all, .(uniprot), summarise,
              symbol = max(symbol),
              min.logFC.10 = min(logFC.10),
              max.logFC.10 = max(logFC.10),
              Abs.logFC.10 = max(abs(logFC.10)),
              min.logFC.600 = min(logFC.600),
              max.logFC.600 = max(logFC.600),
              Abs.logFC.600 = max(abs(logFC.600)),
              min.logFC.1800 = min(logFC.1800),
              max.logFC.1800 = max(logFC.1800),
              Abs.logFC.1800 = max(abs(logFC.1800))
             )

#top.collapse[abs(top.collapse[, "min.logFC.10"]) < abs(top.collapse[, "max.logFC.10"]),]

## 10 sec
top.all.logfc.collapse$logFC.10 <- 0
top.all.logfc.collapse$logFC.10[abs(top.all.logfc.collapse[, "min.logFC.10"]) < abs(top.all.logfc.collapse[, "max.logFC.10"])] = top.all.logfc.collapse$max.logFC.10[abs(top.all.logfc.collapse[, "min.logFC.10"]) < abs(top.all.logfc.collapse[, "max.logFC.10"])]

top.all.logfc.collapse$logFC.10[abs(top.all.logfc.collapse[, "min.logFC.10"]) > abs(top.all.logfc.collapse[, "max.logFC.10"])] = top.all.logfc.collapse$min.logFC.10[abs(top.all.logfc.collapse[, "min.logFC.10"]) > abs(top.all.logfc.collapse[, "max.logFC.10"])]

top.all.logfc.collapse$logFC.10[abs(top.all.logfc.collapse[, "min.logFC.10"]) == abs(top.all.logfc.collapse[, "max.logFC.10"])] = top.all.logfc.collapse$min.logFC.10[abs(top.all.logfc.collapse[, "min.logFC.10"]) == abs(top.all.logfc.collapse[, "max.logFC.10"])]
## 600 sec
top.all.logfc.collapse$logFC.600 <- 0
top.all.logfc.collapse$logFC.600[abs(top.all.logfc.collapse[, "min.logFC.600"]) < abs(top.all.logfc.collapse[, "max.logFC.600"])] = top.all.logfc.collapse$max.logFC.600[abs(top.all.logfc.collapse[, "min.logFC.600"]) < abs(top.all.logfc.collapse[, "max.logFC.600"])]

top.all.logfc.collapse$logFC.600[abs(top.all.logfc.collapse[, "min.logFC.600"]) > abs(top.all.logfc.collapse[, "max.logFC.600"])] = top.all.logfc.collapse$min.logFC.600[abs(top.all.logfc.collapse[, "min.logFC.600"]) > abs(top.all.logfc.collapse[, "max.logFC.600"])]

top.all.logfc.collapse$logFC.600[abs(top.all.logfc.collapse[, "min.logFC.600"]) == abs(top.all.logfc.collapse[, "max.logFC.600"])] = top.all.logfc.collapse$min.logFC.600[abs(top.all.logfc.collapse[, "min.logFC.600"]) == abs(top.all.logfc.collapse[, "max.logFC.600"])]
## 1800 sec
top.all.logfc.collapse$logFC.1800 <- 0
top.all.logfc.collapse$logFC.1800[abs(top.all.logfc.collapse[, "min.logFC.1800"]) < abs(top.all.logfc.collapse[, "max.logFC.1800"])] = top.all.logfc.collapse$max.logFC.1800[abs(top.all.logfc.collapse[, "min.logFC.1800"]) < abs(top.all.logfc.collapse[, "max.logFC.1800"])]

top.all.logfc.collapse$logFC.1800[abs(top.all.logfc.collapse[, "min.logFC.1800"]) > abs(top.all.logfc.collapse[, "max.logFC.1800"])] = top.all.logfc.collapse$min.logFC.1800[abs(top.all.logfc.collapse[, "min.logFC.1800"]) > abs(top.all.logfc.collapse[, "max.logFC.1800"])]

top.all.logfc.collapse$logFC.1800[abs(top.all.logfc.collapse[, "min.logFC.1800"]) == abs(top.all.logfc.collapse[, "max.logFC.1800"])] = top.all.logfc.collapse$min.logFC.1800[abs(top.all.logfc.collapse[, "min.logFC.1800"]) == abs(top.all.logfc.collapse[, "max.logFC.1800"])]

top.all.logfc.collapse.b <- as.data.frame(top.all.logfc.collapse[,c("logFC.10", "logFC.600", "logFC.1800")])
rownames(top.all.logfc.collapse.b) <- top.all.logfc.collapse$uniprot
top.all.logfc.collapse.b$uniprot <- top.all.logfc.collapse$uniprot
top.all.logfc.collapse.b$symbol <- top.all.logfc.collapse$symbol


interactors <- unique(c(interactions.all.d$uniprot_id1,interactions.all.d$uniprot_id2))

node <- as.data.frame(unique(top.all.logfc.collapse.b[which(rownames(top.all.logfc.collapse.b) %in% interactors),]))
node$uniprot <- rownames(node)
node1 <- merge(node, top.all[,c(1,3,5)], by = c("uniprot", "logFC.10"))
node2 <- merge(node1, top.all[,c(1,6,8)], by = c("uniprot", "logFC.600"))
node3 <- merge(node2, top.all[,c(1,9,11)], by = c("uniprot", "logFC.1800"))
node3 <- node3[,c(1,5,4,3,2,6,7,8)]
rownames(node3) <- node3$uniprot

## to get gene_name on top of uniprot sybols
library(dplyr)
df1 <- as.data.frame(interactions.all.d[,c(1,3)])
df2 <- as.data.frame(interactions.all.d[,c(2,4)])
colnames(df1) <- c("uniprot","gene_name")
colnames(df2) <- c("uniprot","gene_name")
int_names <- unique(rbind(df1,df2))
rownames(int_names) <- int_names$uniprot

node4 <- merge(node3,int_names,by="row.names")
node4 <- node4[,-10]

write.table(node4, "../analysis/cytsocape/node4.txt", sep="\t", row.names=FALSE, quote = FALSE)

