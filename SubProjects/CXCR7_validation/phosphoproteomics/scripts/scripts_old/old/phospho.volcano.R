############# Voclano Plots of Experiments
######## the volvano plot helps to get an overview of the differential expressed proteins
######## and is useful to compare different normalization strategies as well as different testing strategies like LRT ebayes, t test, fisher exact etc. 
######## , as well as differetn data fitting with linear model fit, or glm fit

#### use: determine_empirical_control_genes.r
#### use: differential_regulation_raw_data.r
#### use: RUV_normalize_and_differential_regulation.r
#### use: annotate_gene_name.r

if (!requireNamespace('BiocManager', quietly = TRUE))
install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')

library(EnhancedVolcano)
library(foreach)
library(doParallel)

######################################
######################################
### read data if not in the R workspace
setwd("F:/Masterarbeit_BioWis/Phosphoproteomics/Phosphoproteom validation/scripts/")
# e.g.
#top.1800 <- read.table("F:\\Masterarbeit_BioWis\\people\\Mita\\R analysis\\top.1800.txt", sep="\t", header=TRUE, dec=",")

### define input 

input_tables_names <- list("top.10", "top.600", "top.1800")
cr = 1
input_tables <- list(top.10, top.600, top.1800)

#input_tables <- list(top.10)

######################################
######################################
### modify labels and plot volcano
#cl <- 4
#registerDoParallel(cl)
#foreach(top.x=input_tables) %dopar% 

for (i in 1:length(input_tables)) {
	#cr = 1
	top.x = as.data.frame(input_tables[cr])

	#UniprotID <- sapply(strsplit(rownames(ppe0), ";"), "[[", 1)
	#GeneSymbol <- sapply(strsplit(rownames(ppe0), ";"), "[[", 2)
	#Residue <- gsub("[0-9]","", sapply(strsplit(rownames(ppe0), ";"), "[[", 3))
	#Site <- gsub("[A-Z]","", sapply(strsplit(rownames(ppe0), ";"), "[[", 3))
	#Sequence <- sapply(strsplit(rownames(ppe0), ";"), "[[", 4)

	UniprotID<- sapply(strsplit(rownames(top.x), ";"), "[[", 1)
	GeneSymbol <- sapply(strsplit(rownames(top.x), ";"), "[[", 2)
	Site <- sapply(strsplit(rownames(top.x), ";"), "[[", 3)
	id_site <- which(Site != "")
	top.x.b <- top.x[id_site,]

	label <- sapply(strsplit(rownames(top.x.b), ";"),  function(x){paste(x[[2]], x[[3]], sep=".")})

	# define different kinase targets
	pka_target1 <- rownames(top.10[1:10,])

	# create custom key-value pairs for different cell-types
	# this can be achieved with nested ifelse statements
	keyvals.shape <- ifelse(
		rownames(top.x.b) %in% pka_target1, 17,
			3)
	keyvals.shape[is.na(keyvals.shape)] <- 3
	names(keyvals.shape)[keyvals.shape == 3] <- 'sites'
	names(keyvals.shape)[keyvals.shape == 17] <- 'PKA target'

	windows.options(width=30, height=30)
	#dev.new()
	png(file=paste0("../analysis/Volcano_plots/", input_tables_names[cr],".png"), width = 1200, height = 1000, bg = "white")
	
	print(EnhancedVolcano(top.x.b,
	lab = label,
	x = 'logFC',
	y = 'adj.P.Val',
	xlab = bquote(~Log[2]~ 'fold change'),
	pCutoff = 0.05,
	FCcutoff = 0.5,
	cutoffLineType = 'twodash',
	cutoffLineWidth = 0.8,
	pointSize = 4.0,
	labSize = 3.0,
	colAlpha = 0.3,
	legendLabels=c('Not sig.','Not sig. Log2FC','adj.p.value',
	  'adj.p.value & Log2FC'),
	legendPosition = 'right',
	legendLabSize = 16,
	legendIconSize = 4.0,
	drawConnectors = TRUE,
	widthConnectors = 1,
	#maxoverlapsConnectors = 50,
	boxedLabels = TRUE,
	#shapeCustom = keyvals.shape,
	xlim = c(-3.4, 3.4),
	ylim = c(0,7.5),
	#title = "hello"
	title = input_tables_names[cr]
	))
#	+ coord_flip()
	cr = cr + 1
	dev.off()
	
	#Sys.sleep(2)
}
