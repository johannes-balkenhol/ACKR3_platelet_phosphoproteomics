################################################################
##### PSIQUIC network reconstruction 
BiocManager::install("PSICQUIC")
library("PSICQUIC")


################################################################
##### get intesity data

norm_intensity <- SummarizedExperiment::assay(ppe4, "normalised")
UniprotID <- sapply(strsplit(rownames(ppe), ";"), "[[", 1)
GeneSymbol <- sapply(strsplit(rownames(ppe), ";"), "[[", 2)

##### collpase with collpase funiton
norm_intensity.collapse <- phosCollapse(norm_intensity, id=sapply(strsplit(rownames(norm_intensity), ";"), "[[", 1), 
                        stat=apply(abs(norm_intensity), 1, max), by = "max")
						
x_tt <- as.factor(grps)

#################################################################
##### psiquic example
psicquic <- PSICQUIC()
providers(psicquic)


## S4 method for signature 'PSICQUIC'
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
##### psiquic top tables
UniprotID <- sapply(strsplit(rownames(top.10.sign), ";"), "[[", 1)

interactions.top.10 <- interactions(psicquic, id=UniprotID, species="9606", provider="IntAct")



##read detection method table in /anaylsis/network_reconstruction/psiquic_folder
## to pull the interesting detection method

table(interactions.top.10$provider)
table(interactions.top.10$detectionMethod)
table(interactions.top.10$type)

detectionmethod <- c("psi-mi:MI:0407(direct interaction)", "psi-mi:MI:0915(physical association)", 
"psi-mi:MI:0407(direct interaction)")
interactions.top.10.a <- interactions.top.10[which(interactions.top.10$type %in% detectionmethod),]
interactions.top.10.a <- interactions.top.10[which(interactions.top.10$type %in% detectionmethod),]
interactions.top.10.a <- interactions.top.10.a[grep("uniprotkb", interactions.top.10.a$A),]
interactions.top.10.a <- interactions.top.10.a[grep("uniprotkb", interactions.top.10.a$B),]
interactions.top.10.a$A <- gsub("uniprotkb:", "", interactions.top.10.a$A)
interactions.top.10.a$B <- gsub("uniprotkb:", "", interactions.top.10.a$B)
interactions.top.10.a$gene_name1 <- sub(":(\S+)(gene name)", "\\1", interactions.top.10.a$aliasA,perl=TRUE)


## R pattern match is shot, we should use perl
#sub(".*?([0-9]+).*", "\\1", "aaa12xx99",perl=TRUE)
#gsub('([[:alpha:]]+)([0-9]+)([[:alpha:]]+)', '\\2', "aaa12xxx")

#sub("uniprotkb:(.*)?(gene name)", "\\2", interactions.top.10.a$aliasA,perl=TRUE)
#gregexpr("uniprotkb:(.*)?(gene name)", interactions.top.10.a$aliasA,perl=TRUE)
#sub("uniprotkb:(.*)?(gene name)", "\\1", interactions.top.10.a$aliasA)

#gsub("uniprotkb:", "",sapply(strsplit(interactions.top.10.a$aliasA, "\\|"), "[[", 2))

#grep("(gene name)", sapply(strsplit(interactions.top.10.a$aliasA, "\\|")), value=TRUE)

#sapply(strsplit(interactions.top.10.a$aliasA, "\\|")) %like% "gene_name"

#  lapply(strsplit(interactions.top.10.a$aliasA, " \\|"), 
#         \(x) x[grepl("^[GD]", x)])
#str_extract(interactions.top.10.a$aliasA, "\\|uniprotkb:.+(gene name)")
##

interactions.top.10.a$gene_name1 <- str_extract(interactions.top.10.a$aliasA, "\\|uniprotkb:\\S+(gene name)")
interactions.top.10.a$gene_name1 <- gsub("\\|uniprotkb:(.+)\\(gene name","\\1",interactions.top.10.a$gene_name1)

interactions.top.10.a$gene_name2 <- str_extract(interactions.top.10.a$aliasB, "\\|uniprotkb:\\S+(gene name)")
interactions.top.10.a$gene_name2 <- gsub("\\|uniprotkb:(.+)\\(gene name","\\1",interactions.top.10.a$gene_name2)

interactions.top.10.a$score <- gsub(".*intact-miscore:(\\d+)","\\1",interactions.top.10.a$confidenceScore,,perl=TRUE)

interactions.top.10.b <- interactions.top.10.a[,c(1,2,17,18,19,12,7,8,9,13)]

colnames(interactions.top.10.b)[1] = "uniprot_id1"
colnames(interactions.top.10.b)[2] = "uniprot_id2"

interactions.top.10.c <- interactions.top.10.b[interactions.top.10.b$score > 0.5,]