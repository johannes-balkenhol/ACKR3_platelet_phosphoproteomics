

Tc = c1$Tc
clustObj = best$clustObj
min.mem = 0
mfrow = c(1, 1)
new.window = FALSE
llwd=3
	
#FF0030

  clusterindex <- clustObj$cluster
  memship <- clustObj$membership
  memship[memship < min.mem] <- -1
  colorindex <- integer(dim(Tc)[[1]])
  cols <- c("#FF8F00", "#FFA700", "#FFBF00", "#FFD700", 
              "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00", "#AFFF00", 
              "#97FF00", "#80FF00", "#68FF00", "#50FF00", "#38FF00", 
              "#20FF00", "#08FF00", "#00FF10", "#00FF28", "#00FF40", 
              "#00FF58", "#00FF70", "#00FF87", "#00FF9F", "#00FFB7", 
              "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", "#00CFFF", 
              "#00B7FF", "#009FFF", "#0087FF", "#0070FF", "#0058FF", 
              "#0040FF", "#0028FF", "#0010FF", "#0800FF", "#2000FF", 
              "#3800FF", "#5000FF", "#6800FF", "#8000FF", "#9700FF", 
              "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", "#FF00EF", 
              "#FF0030", "#FF0018")

  colorseq <- seq(0, 1, length = length(cols))

j = 1

tmp <- Tc[clusterindex == j, ]
tmpmem <- memship[clusterindex == j, j]
ymin <- min(tmp)
ymax <- max(tmp)
plot(x = NA, xlim = c(1 -0.5, dim(Tc)[[2]] +0.5), 
     ylim = c(ymin, ymax), xlab = "Time Course", ylab = "Standardized Profile", 
     main = paste("Cluster", j, "; size=", nrow(tmp)))


df = data.frame()
for (jj in 1:(length(colorseq) - 1)) {

  tmpcol <- (tmpmem >= colorseq[jj] & tmpmem <= colorseq[jj + 1])
  print(tmpcol)
  range <- paste0("[", colorseq[jj], ",", colorseq[jj + 1], "]")
  tmpind <- which(tmpcol)
  assign("tmpind", tmpind)
  prots <- paste(names(tmpind), collapse = ",")
  output <- c(range, cols[jj], prots)
  df <- rbind(df, output)
  assign("color_ranges" , df)
  for (k in 1:length(tmpind)) {
    print(paste(tmp[tmpind[k], ], cols[jj]))
    graphics::lines(tmp[tmpind[k], ], col = cols[jj],lwd=llwd)
}
}

plot(x = NA, xlim = c(0,3), ylim = c(-2,2),
     xlab = "Time Course", ylab = "Standardized Profile", 
     main = paste("Cluster", j, "; size=", nrow(tmp)))
graphics::lines(tmp[tmpind[45], ], col = "red",lwd=llwd)



#fuzzPlot(c1$Tc, clustObj = best$clustObj, mfrow = c(2,2))


##manual plotting 

library(viridis)
genes_in_cluster_1 <- names(best[["clustObj"]][["cluster"]][best[["clustObj"]][["cluster"]] == 1])
genes_in_cluster_2 <- names(best[["clustObj"]][["cluster"]][best[["clustObj"]][["cluster"]] == 2])

membership_values <- as.data.frame(best[["clustObj"]][["membership"]])
colnames(membership_values) <-  c("cluster1", "cluster2")


membership_values_1 <- membership_values[rownames(membership_values) %in% genes_in_cluster_1,]
for (i in 1:nrow(color_ranges)) { 
  membership_values_1$color <- membership_values_1[membership_values_1$cluster1 %in% ]
  
  membership_values_2 <- membership_values[rownames(membership_values) %in% genes_in_cluster_2,]
  
  cluster1_TC <- c1$Tc[rownames(c1$Tc) %in% genes_in_cluster_1,]
  cluster1_TC <- merge(cluster1_TC, membership_values_1, by = 'row.names')
  rownames(cluster1_TC) <- cluster1_TC$Row.names
  cluster1_TC <- cluster1_TC[-1]
  
  cluster2_TC <- c1$Tc[rownames(c1$Tc) %in% genes_in_cluster_2,]
  cluster2_TC <- merge(cluster2_TC, membership_values_2, by = 'row.names')
  rownames(cluster2_TC) <- cluster2_TC$Row.names
  cluster2_TC <- cluster2_TC[-1]
  
  transposed_data1 <- t(cluster1_TC[1:3])
  transposed_data2 <- t(cluster2_TC[1:3])
  
  par(mfrow =c(2,1)) 
  # Create a line plot
  confidence_scores1 <- cluster1_TC$cluster1
  confidence_scores2 <- cluster2_TC$cluster2
  
  # Create a color gradient using viridis
  color_gradient1 <- viridis(length(confidence_scores1))
  names(color_gradient1) <-  cluster1_TC$cluster1
  
  color_gradient2 <- viridis(length(confidence_scores2))
  names(color_gradient2) <-  cluster2_TC$cluster2
  
  cols <- c("#FF8F00", "#FFA700", "#FFBF00", "#FFD700", 
            "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00", "#AFFF00", 
            "#97FF00", "#80FF00", "#68FF00", "#50FF00", "#38FF00", 
            "#20FF00", "#08FF00", "#00FF10", "#00FF28", "#00FF40", 
            "#00FF58", "#00FF70", "#00FF87", "#00FF9F", "#00FFB7", 
            "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", "#00CFFF", 
            "#00B7FF", "#009FFF", "#0087FF", "#0070FF", "#0058FF", 
            "#0040FF", "#0028FF", "#0010FF", "#0800FF", "#2000FF", 
            "#3800FF", "#5000FF", "#6800FF", "#8000FF", "#9700FF", 
            "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", "#FF00EF", 
            "#FF0030", "#FF0018")
  
  colorseq <- seq(0, 1, length = length(cols))
  
  
  tiff(filename = paste0("../analysis/Clusters/gene_centric_best","_degs_foldchange_manual", ".tiff"),
       width = 4 * 300, 
       height = 8* 300,
       res = 300,
       compression = "lzw")
  
  par(mfrow =c(2,1)) 
  matplot(transposed_data1, type = "l", col = cols, lty = 1, 
          xlab = "Time Points", ylab = "Expression Value",
          main = paste0("cluster1 (n=", length(genes_in_cluster_1), ")"))
  
  matplot(transposed_data2, type = "l", col = cols, lty = 1, 
          xlab = "Time Points", ylab = "Expression Value",
          main = paste0("cluster2 (n=", length(genes_in_cluster_2), ")"))
  
  dev.off()
  
  
