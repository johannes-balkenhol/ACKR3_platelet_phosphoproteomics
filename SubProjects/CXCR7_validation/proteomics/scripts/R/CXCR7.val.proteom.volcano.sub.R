############# Voclano Plots of Experiments
######## the volvano plot helps to get an overview of the differential expressed proteins
######## and is useful to compare different normalization strategies as well as different testing strategies like LRT ebayes, t test, fisher exact etc. 
######## , as well as differetn data fitting with linear model fit, or glm fit

#### use: determine_empirical_control_genes.r
#### use: differential_regulation_raw_data.r
#### use: RUV_normalize_and_differential_regulation.r
#### use: annotate_gene_name.r

# Load packages ----
library(EnhancedVolcano)
library(foreach)
library(doParallel)
library(dplyr)
library(UniprotR)


# Define input ----

names_input <- list("top.10", "top.600", "top.1800", 
                     "top.10.cxcr7.vs.0s", "top.600.cxcr7.vs.0s", "top.1800.cxcr7.vs.0s",
                     "top.10.dmso.vs.0s", "top.600.dmso.vs.0s", "top.1800.dmso.vs.0s")

dfs_input <- list(top.10, top.600, top.1800, 
                  top.10.cxcr7.vs.0s, top.600.cxcr7.vs.0s, top.1800.cxcr7.vs.0s,
                  top.10.dmso.vs.0s, top.600.dmso.vs.0s, top.1800.dmso.vs.0s)


cr = 1
# Modify labels and plot volcano ----

for (i in 1:length(dfs_input)) {
  
  top.x = as.data.frame(dfs_input[i])
  
  UniprotID<- rownames(top.x)
  GeneSymbol <- top.x$symbol
  top.x.b <- top.x[order(top.x$P.Value),]
  rownames(top.x.b) <- paste0(UniprotID, ";", GeneSymbol)
  label <- rownames(top.x.b)
  sign_labels <- rownames(subset(top.x.b, P.Value <0.05 & abs(logFC)>0.5))
  
  png(file=paste0("../analysis/Volcano_plots/", names_input[cr],"_2.png"), width = 1200, height = 1000, bg = "white")
  
  print(EnhancedVolcano(top.x.b,
                        lab = label,
                        selectLab = label[label %in% sign_labels][1:20],
                        x = 'logFC',
                        y = 'P.Value',
                        xlab = bquote(~Log[2]~ 'fold change'),
                        pCutoff = 0.05,
                        FCcutoff = 0.5,
                        cutoffLineType = 'twodash',
                        cutoffLineWidth = 0.8,
                        pointSize = 6.0,
                        titleLabSize =24,
                        labSize = 8,
                        colAlpha = 0.3,
                        legendLabels=c('Not sig.','Not sig. Log2FC','P.value',
                                       'P.value & Log2FC'),
                        legendPosition = 'right',
                        legendLabSize = 20,
                        legendIconSize = 6.0,
                        drawConnectors = TRUE,
                        #widthConnectors = 0.5
                        #maxoverlapsConnectors = 50
                        #boxedLabels = TRUE,
                        xlim = c(- 4, 4),
                        #ylim = c(0, -log10(min(top.x.b$adj.P.Val))*1.1),
                        ylim = c(0, -log10(min(top.x.b$P.Value))),
                        title = names_input[cr]
  ))
  cr = cr + 1
  dev.off()
}

