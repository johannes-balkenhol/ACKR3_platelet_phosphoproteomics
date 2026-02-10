
library(edgeR)

scaled_x <- data.frame(2^SummarizedExperiment::assay(ppe_imputed_scaled,"scaled"))


x <-  DGEList(counts = scaled_x,lib.size = colSums(scaled_x), 
              samples = sample_info, group = grps,)


x <- calcNormFactors(x, method = "TMM")

x$samples$norm.factors


lcpm <- cpm(x, log=TRUE)


col.group <- as.factor(x$samples$time)
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)

plotMDS(lcpm, labels=x$samples$sample, col=col.group)
title(main="A. Sample groups")
