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

## or
top.rownames <- unique(c(rownames(top.10.sign),rownames(top.600.sign),rownames(top.1800.sign)))
top.norm_intensity <- norm_intensity[top.rownames,]

rownames(top.norm_intensity) <- sapply(strsplit(rownames(top.norm_intensity), ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})



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
table(interactions.top.10$provider)
