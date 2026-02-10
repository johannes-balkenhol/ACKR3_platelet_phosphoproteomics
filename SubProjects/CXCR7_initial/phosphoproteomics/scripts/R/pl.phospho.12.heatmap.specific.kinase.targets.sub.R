
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(reshape2)
library(ggplot2)
library(PhosR)
library(dplyr)
library(matrixStats)
library(stringr)
library(cowplot)


#setwd("/Users/ozgeosmanoglu/Nextcloud/Phosphoproteom/phosphoproteom analysis 2/")



data(PhosphoSitePlus)



pka_targets_psplus <- PhosphoSite.human[["PRKACA"]]
pkg_targets_psplus <- PhosphoSite.human[["PRKG1"]]


#prep data####
dataset <- norm_intensity_filter #this is the same as norm_intensity
dataset_df <- as.data.frame(dataset)

#get the list of significant ids in all

input = list(top.filter.10, top.filter.30, top.filter.60, top.filter.300, top.filter.600, top.filter.900, top.filter.1800)
#top.10.dmso.vs.0s, top.600.dmso.vs.0s, top.1800.dmso.vs.0s,
#top.10.cxcr7.vs.0s, top.600.cxcr7.vs.0s, top.1800.cxcr7.vs.0s)

names_input = c("10", "30", "60", "300", "600", "900", "1800")
#"10.dmso.vs.0s", "600.dmso.vs.0s", "1800.dmso.vs.0s",
#"10.cxcr7.vs.0s", "600.cxcr7.vs.0s", "1800.cxcr7.vs.0s")

for (i in 1:length(input)) {
  #df_pl = top.10.cxcr7.vs.0s
  df_pl2 = input[[i]]
  df_pl2 <- df_pl2[df_pl2[, "PValue"] <0.05,]
  assign(paste0("top.filter.", names_input[[i]], ".sign"), df_pl2)
}

# union_sig_collapse <- Reduce(union, list(rownames(top.collapse.10.sign), rownames(top.collapse.30.sign), 
#                                 rownames(top.collapse.60.sign), row.names(top.collapse.300.sign),
#                                 rownames(top.collapse.600.sign), rownames(top.collapse.900.sign), 
#                                 rownames(top.collapse.1800.sign)))
union_sig <- Reduce(union, list(rownames(top.filter.10.sign), rownames(top.filter.30.sign), 
                                rownames(top.filter.60.sign), row.names(top.filter.300.sign),
                                rownames(top.filter.600.sign), rownames(top.filter.900.sign), 
                                rownames(top.filter.1800.sign)))
dataset_df <- dataset_df[rownames(dataset_df) %in% union_sig, ]

#if you want to get all targets

union_all <- Reduce(union, list(rownames(top.filter.10), rownames(top.filter.30), 
                                rownames(top.filter.60), row.names(top.filter.300),
                                rownames(top.filter.600), rownames(top.filter.900), 
                                rownames(top.filter.1800)))
dataset_df <- dataset_df[rownames(dataset_df) %in% union_all, ]





#continue
name <- sapply(strsplit(rownames(dataset_df), ";"), "[[", 2)
site <- sapply(strsplit(rownames(dataset_df), ";"), "[[", 3)
dataset_df$namesite <- paste0(name,";",site,";")


##specific kinases####
kinases = c("PRKACA", "PRKG1", "PRKAA1", "CDK1", "MAPKAPK2", "AKT1", "MTOR", "SRC", "CDK2", "PRKCA", "CDK5", "CSNK2A1")

#x = kinases[[1]]
kinase_function <- function(x) {
  kinase <- PhosphoSite.human[[x]]
  targets_norm_abundance <- dataset_df[dataset_df$namesite %in% kinase,1:53]
  targets_norm_abundance <-targets_norm_abundance[order(row.names(targets_norm_abundance)), ]
  targets_norm_abundance <- as.matrix(targets_norm_abundance)
  
  targets_pvals <- top.all[rownames(top.all) %in% rownames(targets_norm_abundance), c(5,8,11,14,17,20,23)]
  targets_pvals<-targets_pvals[order(row.names(targets_pvals)), ]
  
  targets_pvals <- as.matrix(targets_pvals)
  rownames(targets_pvals) <- sapply(strsplit(rownames(targets_pvals), ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})
  targets_pvals <- as.data.frame(targets_pvals)
  targets_pvals[targets_pvals <= 0.05] <- "*"
  targets_pvals[targets_pvals> 0.05] <- ""
  targets_pvals$adj.P.Val.00 <- ""
  targets_pvals <- targets_pvals[c(8,1:7)]
  
  targets_logfcs <- top.all[rownames(top.all) %in% rownames(targets_norm_abundance), c(3,6,9,12,15,18,21)]
  targets_logfcs<-targets_logfcs[order(row.names(targets_logfcs)), ]
  targets_logfcs <- as.matrix(targets_logfcs)
  rownames(targets_logfcs) <- sapply(strsplit(rownames(targets_logfcs), ";"),  
                                     function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})
  targets_logfcs <- as.data.frame(targets_logfcs)
  
  targets_logfcs<-round(targets_logfcs, digits = 1)
  
  for (j in 1:(length(colnames(targets_pvals))-1)) {
    print(j)
    for (i in 1:length(rownames(targets_logfcs))) {
      if (targets_pvals[i, j+1] == "") {
        targets_logfcs[i, j] <- ""}}}
  
  targets_logfcs$logFC.00<- ""
  targets_logfcs <- targets_logfcs[c(8,1:7)]
  targets_logfcs<-targets_logfcs[order(row.names(targets_logfcs)), ]
  
  
  
  rownames(targets_norm_abundance) <- sapply(strsplit(rownames(targets_norm_abundance), ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})
  
  x_tt <- as.factor(grps)
  
  
  t00 <- rowMedians(targets_norm_abundance[,1:7])
  t10 <- rowMedians(targets_norm_abundance[,8:13])
  t30 <- rowMedians(targets_norm_abundance[,14:20])
  t60 <- rowMedians(targets_norm_abundance[,21:27])
  t300 <- rowMedians(targets_norm_abundance[,28:34])
  t600 <- rowMedians(targets_norm_abundance[,35:40])
  t900 <- rowMedians(targets_norm_abundance[,41:47])
  t1800 <- rowMedians(targets_norm_abundance[,48:53])
  
  
  targets_norm_abundance.collapse <- cbind(t00,t10,t30,t60,t300,t600,t900,t1800)
  rownames(targets_norm_abundance.collapse) <- rownames(targets_norm_abundance)
  colnames(targets_norm_abundance.collapse) <- c("t00","t10","t30","t60","t300","t600","t900", "t1800")
  targets_norm_abundance.collapse <- as.data.frame(targets_norm_abundance.collapse)
  
  targets_norm_abundance2 <- targets_norm_abundance.collapse
  
  
  #dd <- melt(as.matrix(targets_norm_abundance2), variable.name = "normalised")
  #ggplot(dd, aes(value, colour = Var2)) + geom_density(bw = "SJ")
  
  
  
  ###############################################################################
  ## use median, as more robust than average, for the Z values in heatmap
  cal_z_score <- function(x){
    (x - mean(x)) / sd(x)
  }
  
  ########### coloring according to sample names
  col_groups <- substr(colnames(targets_norm_abundance2), 1, 5)
  table(col_groups)
  mat_col <- data.frame(time = col_groups)
  rownames(mat_col) <- colnames(targets_norm_abundance2)
  mat_colors <- list(time = brewer.pal(8, "BuPu"))
  #mat_colors <- list(group = brewer.pal(length(table(col_groups)), "Set3"))
  names(mat_colors$time) <- unique(col_groups)
  
  
  ############ create a simple heatmap
  targets_data_subset_norm <- t(apply(targets_norm_abundance2, 1, cal_z_score))
  #rownames(data_subset_norm) <- c(sapply(strsplit(rownames(data_subset_norm), ";"), "[[", c(1,2,3)), ";", sapply(strsplit(rownames(data_subset_norm), ";"), "[[", 2), ";" , sapply(strsplit(rownames(data_subset_norm), ";"), "[[", 3))
  #rownames(data_subset_norm) <- sapply(strsplit(rownames(data_subset_norm), ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})
  
  
  ############ include breaks in the heatmap 
  ############ for better visualization in tailed data
  ############ as we use a color code cutoff regarding quantiles
  
  mat_breaks <- seq(min(targets_norm_abundance2, na.rm=TRUE), max(targets_norm_abundance2, na.rm=TRUE), length.out = 20)
  
  
  # define function
  quantile_breaks <- function(xs, n) {
    breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm=TRUE)
    breaks[!duplicated(breaks)]
  }
  
  mat_breaks <- quantile_breaks(targets_data_subset_norm, n = 40)
  test_labels <- as.matrix(targets_logfcs) 
  
  # Define the threshold for dynamic text colors - new coloring
  threshold <- mat_breaks[10]
  # Generate dynamic text colors based on cell values
  number_colors <- ifelse(targets_data_subset_norm < threshold, "white", "black")
  
  ph <- pheatmap(
    mat = targets_data_subset_norm,
    color = colorRampPalette(c("navy", "white", "red"))(length(mat_breaks) - 1),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    breaks = mat_breaks,
    silent = TRUE # Suppress immediate plotting to adjust data
  )
  
  # Reorder `number_colors` to match clustering
  reordered_number_colors <- number_colors[ph$tree_row$order, ]
  
  myheatmap <- pheatmap(mat = targets_data_subset_norm, 
                        color = colorRampPalette(c("navy", "white", "red"))(length(mat_breaks) - 1),
                        gaps_row=c(5,10,15,20,25,30,35,40,45,50, 55),
                        scale="row",
                        na_col = "grey",
                        breaks = mat_breaks,
                        border_color = "white", 
                        show_colnames = TRUE, 
                        show_rownames = TRUE, 
                        annotation_col = mat_col, 
                        annotation_colors = mat_colors, 
                        drop_levels = TRUE, 
                        fontsize = 10, 
                        cluster_cols = FALSE,
                        cluster_rows = TRUE,
                        cex=1,
                        clustering_distance_rows="euclidean",
                        clustering_distance_cols="euclidean",
                        clustering_method="complete",
                        main = paste(x, "targets"),
                        display_numbers = test_labels,
                        number_color = reordered_number_colors, 
                        fontsize_number = 10, #increased font number
                        cellheight=12, cellwidth = 24)
  
  #new folder - shortened
  tiff(filename = paste0("../analysis/Heatmap/Shortened/Heatmap_", x, "_targets_sign.tiff", sep =""), #for kinases
       #tiff(filename = paste0("analysis/Heatmap/Heatmap_", protein, "_sign.tiff", sep =""), #for proteins
       width = 10 * 300, 
       height = 10 * 300,
       res = 300,
       compression = "lzw")
  
  print(myheatmap)
  
  dev.off()
  
  
  input_balloon <- targets_data_subset_norm[order(myheatmap$tree_row$order),]
  balloon_dd <- reshape2::melt(input_balloon, variable.name = "sample")
  targets_logfcs2 <- targets_logfcs[order(myheatmap$tree_row$order),]
  colnames(targets_logfcs2) <-  gsub("logFC.", "t", colnames(targets_logfcs2))
  balloon_dd2 <- reshape2::melt(as.matrix(targets_logfcs2), variable.name = "sample")
  colnames(balloon_dd2) <- c("Var1", "Var2", "logfc")
  balloon_dd_input <- cbind(balloon_dd, balloon_dd2)
  balloon_dd_input <- balloon_dd_input[c(1,2,3,6)]
  
  balloon_plot <-  ggplot(balloon_dd_input, aes(x=factor(Var2), y=factor(Var1), 
                                                color=value, size=value)) +
    geom_point() +    # plot as points
    geom_text(aes(label=logfc, fontface = "bold"), alpha=1.0, size=3, color = "green") +   # display the value next to the "balloons"
    #scale_alpha_continuous(range=c(0.9, 1)) +
    scale_color_gradient2(low = "navy", mid= "white", high = "red")+
    scale_size_area(max_size = 10) +
    ggtitle(paste(x, "targets"))+
    theme_cowplot() +
    theme(axis.line = element_blank(),            # disable axis lines
          axis.title = element_blank(),           # disable axis titles
          panel.border = element_blank(),         # disable panel border
          panel.grid.major.x = element_blank(),   # disable lines in grid on X-axis
          panel.grid.minor.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 90))  # disable lines in grid on X-axis
  
  # tiff(filename = paste0("../analysis/Heatmap/Balloon_", x, "_targets.tiff", sep =""), #for kinases
  #      #tiff(filename = paste0("analysis/Heatmap/Heatmap_", protein, "_sign.tiff", sep =""), #for proteins
  #      width = 6 * 300, 
  #      height = 6 * 300,
  #      res = 300,
  #      compression = "lzw")
  # 
  # print(balloon_plot)
  # 
  # dev.off()
  
}
###
i=1
for (i in 1:length(kinases)) {
  kinasename <-  kinases[[i]]
  kinase_function(kinasename)
}




##specific proteins####

#PDE3A;Q14432 in raw data
#VASP; P50552 in raw data
#ENSA; O43768 in raw data #not in at the end
#AMPK subunits/isoforms
##AMPK; Q9Y478 in raw data 
##AMPK1; Q13131 in raw data #not in at the end
##PRKAB2: O43741 in raw data 
##AMPKG3; Q9UGI9 NOT in raw data
##PRKAG1; P54619 NOT in raw data
##PRKAG2; Q9UGJ0 NOT in raw data
##AMPK2; P54646 NOT in raw data
#ACACA; Q13085 in raw data #not in at the end
#ACACB; O00763 NOT in raw data

protein_list <- list("Q14432", "P50552", "O43768", "Q9Y478", "Q13131", "O43741", 
                     "Q9UGI9", "P54619", "Q9UGJ0", "P54646", "Q13085", "O00763")
protein <- paste(protein_list, collapse="|")


##clusters## used directly norm_intensity_filter
x <- "cluster8"

#for taking full cluster

prot_spec_anova_cluster <- as.data.frame(clustobjs[clustobjs$cluster == "1",]) #CHANGE CLUSTER NAME
cluster_prots <- unique(unlist(strsplit(prot_spec_anova_cluster$protein, split = "\\|")))

##for mapping to full heatmap data in the heatmap script
clustobjswithIDs <- data.frame()
for (i in 1:length(unique(clustobjs$cluster))) {
  prot_spec_anova_cluster <- as.data.frame(clustobjs[clustobjs$cluster == i,]) #CHANGE CLUSTER NAME
  cluster_prots <- unique(unlist(strsplit(prot_spec_anova_cluster$protein, split = "\\|")))
  proteins <- merged[merged$names %in% cluster_prots, "id"] #merged2 for site-centric
  clustobjsadd <- data.frame(ID  = proteins, cluster = i)
  clustobjswithIDs <- rbind(clustobjswithIDs, clustobjsadd)
}


#for taking prots in the pathways
prot_spec_anova_cluster <- as.data.frame(best$enrichList$`cluster 8`[1:10,]) #CHANGE CLUSTER NAME
cluster_prots <- unique(unlist(strsplit(prot_spec_anova_cluster$substrates, split = "\\|")))

# mydf <- prot_spec_anova_cluster3_platelet
# mydf <- mydf %>% 
#   mutate(substrates = strsplit(as.character(substrates), "\\|")) %>% 
#   unnest(substrates)
# write.table(mydf, "../analysis/Clusters/gene_centric_best_cluster3_platelet_annotation.txt", sep = "\t")


##or for each process, separate

for (i in 1:nrow(prot_spec_anova_cluster)) {
  substrates <- unlist(strsplit(as.character(prot_spec_anova_cluster$substrates[i]), "\\|"))
  kinase <- as.character(prot_spec_anova_cluster$kinase[i])
  
  # Create a vector dynamically based on the row number
  vector_name <- paste("cluster_prots", "_",prot_spec_anova_cluster$kinase[i] , sep = "")
  
  # Assign values to the dynamically named vector
  assign(vector_name, substrates)
  
}

#cluster3_prots_platelet <- unique(unlist(strsplit(prot_spec_anova_cluster3_platelet$substrates, split = "\\|")))
#cluster3_prots_pi3kakt <- unique(unlist(strsplit(prot_spec_anova_cluster3_pi3kakt$substrates, split = "\\|")))
#cluster3_prots_rest <- unique(unlist(strsplit(prot_spec_anova_cluster3_rest$substrates, split = "\\|")))


#protein <- paste(sort(cluster_prots), collapse = "|")

proteins <- merged[merged$names %in% cluster_prots, "id"] #merged2 for site-centric
#proteins <- merged[str_detect(merged$id, paste(protein, ";", sep = "")),"id"]



##KSEA results
KSDataInput <- KSData.dataset10

KSDataInput$substrate <- paste0(KSDataInput$Substrate.Gene, ";", KSDataInput$Substrate.Mod)
KSDataAnnot <- KSDataInput[KSDataInput$Kinase.Gene %in% KSEASignKinases, c(6,1)]

mapdf <- data.frame(namesite = str_sub(dataset_df$namesite, end = -2), ID = rownames(dataset_df))

kinase = KSDataAnnot[(str_detect(KSDataAnnot$substrate, paste(mapdf$namesite, collapse = "|"))), ]
colnames(kinase) <-  c("namesite", "Kinase")
mapdf_sub <- merge(mapdf, kinase, by = "namesite")

x = unique(mapdf_sub$Kinase)[2]
KSSignSubIDs <- mapdf_sub$ID


##annot for heatmap with all kinases
mapdf_sub <- mapdf_sub[order(mapdf_sub$Kinase),]
mapdf_sub2 <-  aggregate(Kinase ~ ID, mapdf_sub, FUN = paste, collapse=",")
rownames(mapdf_sub2) <- sapply(strsplit(mapdf_sub2$ID, ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})
mapdf_sub2 <- mapdf_sub2[order(mapdf_sub2$Kinase),]
mapdf_sub2 <-  mapdf_sub2 %>% select(-ID)

#separate ksea kinases

x = unique(mapdf_sub$Kinase)[5]
KSSignSubIDs <- mapdf_sub[mapdf_sub$Kinase == x, "ID"]


#protein <-  protein_list[[1]]
#protein <- "GAPDH"
#protein_function <- function(x) {
#uncomment for specific prot names:
#targets_norm_abundance <- dataset_df[str_detect(row.names(dataset_df), KSEA600subs), 1:53]
#for cluster phosphosites use following directly:
targets_norm_abundance <- dataset_df[rownames(dataset_df) %in% KSSignSubIDs, 1:53]
#continue from here
targets_norm_abundance <-targets_norm_abundance[order(row.names(targets_norm_abundance)), ]
targets_norm_abundance <- as.matrix(targets_norm_abundance)

targets_pvals <- top.all[rownames(top.all) %in% rownames(targets_norm_abundance), c(5,8,11,14,17,20,23)]
targets_pvals<-targets_pvals[order(row.names(targets_pvals)), ]

targets_pvals <- as.matrix(targets_pvals)
rownames(targets_pvals) <- sapply(strsplit(rownames(targets_pvals), ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})
targets_pvals <- as.data.frame(targets_pvals)
targets_pvals[targets_pvals <= 0.05] <- "*"
targets_pvals[targets_pvals> 0.05] <- ""
targets_pvals$adj.P.Val.00 <- ""
targets_pvals <- targets_pvals[c(8,1:7)]

targets_logfcs <- top.all[rownames(top.all) %in% rownames(targets_norm_abundance), c(3,6,9,12,15,18,21)]
targets_logfcs<-targets_logfcs[order(row.names(targets_logfcs)), ]
targets_logfcs <- as.matrix(targets_logfcs)
rownames(targets_logfcs) <- sapply(strsplit(rownames(targets_logfcs), ";"),  
                                   function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})
targets_logfcs <- as.data.frame(targets_logfcs)

targets_logfcs<-round(targets_logfcs, digits = 2)

for (j in 1:(length(colnames(targets_pvals))-1)) {
  print(j)
  for (i in 1:length(rownames(targets_logfcs))) {
    if (targets_pvals[i, j+1] == "") {
      targets_logfcs[i, j] <- ""}}}

targets_logfcs$logFC.00<- ""
targets_logfcs <- targets_logfcs[c(8,1:7)]
targets_logfcs<-targets_logfcs[order(row.names(targets_logfcs)), ]

rownames(targets_norm_abundance) <- sapply(strsplit(rownames(targets_norm_abundance), ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})
x_tt <- as.factor(grps)

t00 <- rowMedians(targets_norm_abundance[,1:7])
t10 <- rowMedians(targets_norm_abundance[,8:13])
t30 <- rowMedians(targets_norm_abundance[,14:20])
t60 <- rowMedians(targets_norm_abundance[,21:27])
t300 <- rowMedians(targets_norm_abundance[,28:34])
t600 <- rowMedians(targets_norm_abundance[,35:40])
t900 <- rowMedians(targets_norm_abundance[,41:47])
t1800 <- rowMedians(targets_norm_abundance[,48:53])

targets_norm_abundance.collapse <- cbind(t00,t10,t30,t60,t300,t600,t900,t1800)
rownames(targets_norm_abundance.collapse) <- rownames(targets_norm_abundance)
colnames(targets_norm_abundance.collapse) <- c("t00","t10","t30","t60","t300","t600","t900", "t1800")
targets_norm_abundance.collapse <- as.data.frame(targets_norm_abundance.collapse)

targets_norm_abundance2 <- targets_norm_abundance.collapse


#dd <- melt(as.matrix(targets_norm_abundance2), variable.name = "normalised")
#ggplot(dd, aes(value, colour = Var2)) + geom_density(bw = "SJ")



###############################################################################
## use median, as more robust than average, for the Z values in heatmap
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

########### coloring according to sample names
col_groups <- substr(colnames(targets_norm_abundance2), 1, 5)
table(col_groups)


mat_col <- data.frame(time = col_groups)
rownames(mat_col) <- colnames(targets_norm_abundance2)

mat_colors <- list(time = brewer.pal(8, "BuPu"))
#mat_colors <- list(group = brewer.pal(length(table(col_groups)), "Set3"))

names(mat_colors$time) <- unique(col_groups)

############ create a simple heatmap
targets_data_subset_norm <- t(apply(targets_norm_abundance2, 1, cal_z_score))
#rownames(data_subset_norm) <- c(sapply(strsplit(rownames(data_subset_norm), ";"), "[[", c(1,2,3)), ";", sapply(strsplit(rownames(data_subset_norm), ";"), "[[", 2), ";" , sapply(strsplit(rownames(data_subset_norm), ";"), "[[", 3))
#rownames(data_subset_norm) <- sapply(strsplit(rownames(data_subset_norm), ";"),  function(x){paste(x[[1]], x[[2]], x[[3]], sep="_")})


############ include breaks in the heatmap 
############ for better visualization in tailed data
############ as we use a color code cutoff regarding quantiles

mat_breaks <- seq(min(targets_norm_abundance2, na.rm=TRUE), max(targets_norm_abundance2, na.rm=TRUE), length.out = 20)


# define function
quantile_breaks <- function(xs, n) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm=TRUE)
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(targets_data_subset_norm, n = 40)


test_labels <- as.matrix(targets_logfcs) 

myheatmap <- pheatmap(mat = targets_data_subset_norm, 
                      color = colorRampPalette(c("navy", "white", "red"))(length(mat_breaks) - 1),
                      gaps_row=c(5,10,15,20,25,30,35,40,45,50, 55),
                      scale="row",
                      na_col = "grey",
                      breaks = mat_breaks,
                      border_color = "white", 
                      show_colnames = TRUE, 
                      show_rownames = TRUE, 
                      annotation_col = mat_col, 
                      annotation_colors = mat_colors, 
                      #annotation_row = mapdf_sub2,
                      drop_levels = TRUE, 
                      fontsize = 10, 
                      cluster_cols = FALSE,
                      cluster_rows = T,
                      cex=1,
                      clustering_distance_rows="euclidean",
                      clustering_distance_cols="euclidean",
                      clustering_method="complete",
                      main = paste(x, "sites"),
                      display_numbers = test_labels,
                      number_color = "green", 
                      fontsize_number = 7.5,
                      cellheight=10, cellwidth = 30
)

tiff(filename = paste0("../analysis/KSEA/Heatmap_", x, "_substrates.tiff"), #for proteins
     #tiff(filename = paste0("../analysis/Clusters/gene_centric/Anova/gene_centric_best_", x, "_anova.tiff", sep =""), #for proteins
     #tiff(filename = paste0("../analysis/Heatmap/Heatmap_", x, "_sign.tiff", sep =""), #for proteins
     width = 10 * 300, 
     height = 10 * 300,
     res = 300,
     compression = "lzw")

print(myheatmap)

dev.off()





##BALLOON PLOT

input_balloon <- targets_data_subset_norm[order(myheatmap$tree_row$order),]
balloon_dd <- reshape2::melt(input_balloon, variable.name = "sample")
targets_logfcs2 <- targets_logfcs[order(myheatmap$tree_row$order),]
colnames(targets_logfcs2) <-  gsub("logFC.", "t", colnames(targets_logfcs2))
balloon_dd2 <- reshape2::melt(as.matrix(targets_logfcs2), variable.name = "sample")
colnames(balloon_dd2) <- c("Var1", "Var2", "logfc")
balloon_dd_input <- cbind(balloon_dd, balloon_dd2)
balloon_dd_input <- balloon_dd_input[c(1,2,3,6)]

balloon_plot <-  ggplot(balloon_dd_input, aes(x=factor(Var2), y=factor(Var1), 
                                              color=value, size = value)) +
  geom_point()+
  #geom_point(shape = 21, color = "black") +    # plot as points
  geom_text(aes(label=logfc, fontface = "bold"), alpha=1.0, size=5, color = "green") +   # display the value next to the "balloons"
  scale_alpha_continuous(range=c(0.9, 1)) +
  scale_color_gradient2(low = "navy", mid= "white", high = "red")+
  scale_size_area(max_size = 15) +
  ggtitle(paste(x, "targets"))+
  theme_cowplot() +
  theme(axis.line = element_blank(),            # disable axis lines
        axis.title = element_blank(),           # disable axis titles
        panel.border = element_blank(),         # disable panel border
        panel.grid.major.x = element_blank(),   # disable lines in grid on X-axis
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90))  # disable lines in grid on X-axis

tiff(filename = paste0("../analysis/Clusters/gene_centric/Anova/site_centric_best_", x, "_anova_balloon.tiff", sep =""), #for proteins
     #tiff(filename = paste0("analysis/Heatmap/Heatmap_", protein, "_sign.tiff", sep =""), #for proteins
     width = 8 * 300, 
     height = 16 * 300,
     res = 300,
     compression = "lzw")

print(balloon_plot)

dev.off()
}
###
for (i in 1:length(protein_list)) {
  protein <-  protein_list[[i]]
  protein_function(protein)
}



#######


