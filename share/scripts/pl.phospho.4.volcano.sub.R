############# Voclano Plots of Experiments
######## the volvano plot helps to get an overview of the differential expressed proteins
######## and is useful to compare different normalization strategies as well as different testing strategies like LRT ebayes, t test, fisher exact etc. 
######## , as well as differetn data fitting with linear model fit, or glm fit

#### use: determine_empirical_control_genes.r
#### use: differential_regulation_raw_data.r
#### use: RUV_normalize_and_differential_regulation.r
#### use: annotate_gene_name.r


library(EnhancedVolcano)
library(foreach)
library(doParallel)

######################################
######################################

### define input 
input = list(top.filter.10, top.filter.600, top.filter.1800, 
             top.filter.10.dmso.vs.0s, top.filter.600.dmso.vs.0s, top.filter.1800.dmso.vs.0s,
             top.filter.10.cxcr7.vs.0s, top.filter.600.cxcr7.vs.0s, top.filter.1800.cxcr7.vs.0s)

names_input = c("10", "600", "1800", 
                "10.dmso.vs.0s", "600.dmso.vs.0s", "1800.dmso.vs.0s",
                "10.cxcr7.vs.0s", "600.cxcr7.vs.0s", "1800.cxcr7.vs.0s")

#input_tables_names <- list("top.10", "top.600", "top.1800")


cr = 1



######################################
######################################
### modify labels and plot volcano
#cl <- 4
#registerDoParallel(cl)
#foreach(top.x=input_tables) %dopar% 

for (i in 1:length(input)) {
  #cr = 1
  top.x = as.data.frame(input[i])
  
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
  top.x.b <- top.x.b[order(top.x.b$PValue),]
  top.x.b$adj.P.Val <- top.x.b$PValue #by filter adj.p.value is named pvalue
  
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
  
  #windows.options(width=30, height=30)
  #dev.new()
  png(file=paste0("../analysis/Volcano_plots/", names_input[i],".png"), 
      width = 1200, height = 1000,
      bg = "white")
  print(EnhancedVolcano(top.x.b,
                        lab = label,
                        selectLab = label[1:20],
                        x = 'logFC',
                        y = 'adj.P.Val',
                        xlab = bquote(~Log[2]~ 'fold change'),
                        pCutoff = 0.05,
                        FCcutoff = 0.5,
                        cutoffLineType = 'twodash',
                        cutoffLineWidth = 0.8,
                        pointSize = 8.0,
                        labSize = 8.0,
                        colAlpha = 0.3,
                        axisLabSize = 22,
                        legendLabels=c('Not sig.','Not sig. Log2FC','adj.p.value',
                                       'adj.p.value & Log2FC'),
                        legendPosition = 'right',
                        legendLabSize = 16,
                        legendIconSize = 4.0,
                        drawConnectors = TRUE,
                        widthConnectors = 0.5,
                        maxoverlapsConnectors = 50,
                        boxedLabels = TRUE,
                        #shapeCustom = keyvals.shape,
                        xlim = c(-4, 4),
                        ylim = c(0, -log10(min(top.x.b$PValue))*1.1),
                        #title = "hello"
                        title = names_input[i]
  ))
  #	+ coord_flip()
  dev.off()
  
  #Sys.sleep(2)
}
