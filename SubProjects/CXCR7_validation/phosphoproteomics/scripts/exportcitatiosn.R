
used_packages <-  c("BiocManager", "ClueR", "clusterProfiler", "DESeq2", "dplyr", 
                    "EDASeq", "EnhancedVolcano", "enrichplot", "ggplot2", "ggpubr", 
                    "KSEAapp", "limma", "matrixStats", "OmnipathR", "org.Hs.eg.db", 
                    "pheatmap", "PhosR", "RColorBrewer", "reactome.db", "ReactomePA", 
                    "RMySQL", "ruv", "RUVSeq", "stats", "SummarizedExperiment", "tibble", 
                    "tidyr", "directPA","annotate", "basicPlotteR", "calibrate", 
                    "cowplot", "createKEGGdb", "doParallel", "foreach", "GGally", 
                    "ggsci", "ggVennDiagram", "grDevices", "network", "parallel", 
                    "plyr", "psych", "remotes", "reshape2", "rlist", "sjmisc", "sqldf", 
                    "stringr", "tools", "viridis")




for (i in seq_along(used_packages)) {
  capture.output(utils:::print.bibentry(citation(used_packages[i]), 
                                        style = "Bibtex"),
                 file = paste0("/Users/ozgeosmanoglu/Desktop/rpackages/endnote_import", 
                               i, ".bib"))
}


