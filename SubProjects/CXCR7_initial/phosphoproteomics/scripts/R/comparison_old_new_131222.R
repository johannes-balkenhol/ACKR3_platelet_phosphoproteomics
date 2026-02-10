####################################################################################
##### comparison with old data
library(VennDiagram)
library("ggpubr")



## check intersection
#top.all$uniprot_old <- sapply(strsplit(rownames(top.all), ";"), "[[", 1)
old_data<- read.table("reanalysis_141122/top.all_olddata_ÖO_121222.txt")
old_data$site_old <- sapply(strsplit(rownames(old_data), ";"), "[[", 3)
old_data$sequence_old <- sapply(strsplit(rownames(old_data), ";"), "[[", 4)
old_data$id_old <- paste(old_data$uniprot, old_data$site_old, old_data$sequence_old, sep = ";")
old_data$id2_old <- paste(old_data$uniprot, old_data$site_old, sep = ";")


old_data_raw <- raw_abundance2
old_data_raw$uniprot <- sapply(strsplit(rownames(old_data_raw), ";"), "[[", 1)
old_data_raw$site <- sapply(strsplit(rownames(old_data_raw), ";"), "[[", 3)
old_data_raw$sequence<- sapply(strsplit(rownames(old_data_raw), ";"), "[[", 4)
old_data_raw$id_oldraw <- paste(old_data_raw$uniprot, old_data_raw$site, old_data_raw$sequence, sep = ";")
old_data_raw$id2_oldraw <- paste(old_data_raw$uniprot, old_data_raw$site, sep = ";")


new_data<- read.table("../phosphoproteom validation/data/processed_data/top.all.txt", 
                      sep="\t", header=TRUE, dec=".")
new_data$uniprot_new <- sapply(strsplit(rownames(new_data), ";"), "[[", 1)
new_data$site_new <- sapply(strsplit(rownames(new_data), ";"), "[[", 3)
new_data$sequence_new <- sapply(strsplit(rownames(new_data), ";"), "[[", 4)
new_data$id_new <- paste(new_data$uniprot_new, new_data$site_new, new_data$sequence_new, sep = ";")
new_data$id2_new <- paste(new_data$uniprot_new, new_data$site_new, sep = ";")

new_data_toCtrl0 <- read.table("../phosphoproteom validation/data/processed_data/top.all_toCtrl0.txt", 
                               sep="\t", header=TRUE, dec=".")
new_data_toCtrl0$uniprot_new <- sapply(strsplit(rownames(new_data_toCtrl0), ";"), "[[", 1)
new_data_toCtrl0$site_new <- sapply(strsplit(rownames(new_data_toCtrl0), ";"), "[[", 3)
new_data_toCtrl0$sequence_new <- sapply(strsplit(rownames(new_data_toCtrl0), ";"), "[[", 4)
new_data_toCtrl0$id_new <- paste(new_data_toCtrl0$uniprot_new, new_data_toCtrl0$site_new, new_data_toCtrl0$sequence_new, sep = ";")
new_data_toCtrl0$id2_new <- paste(new_data_toCtrl0$uniprot_new, new_data_toCtrl0$site_new, sep = ";")

raw_abundance_new <- read.table("../phosphoproteom validation/data/raw_data/A08_val_phosphoR_v3.txt", sep="\t", header=TRUE, dec=",")
raw_abundance2_new <- raw_abundance_new[,c(-1)]
rownames(raw_abundance2_new) <- raw_abundance_new[,1]
raw_abundance2_new <- as.data.frame(raw_abundance2_new)
new_data_raw <- raw_abundance2_new
new_data_raw$uniprot <- sapply(strsplit(rownames(new_data_raw), ";"), "[[", 1)
new_data_raw$site <- sapply(strsplit(rownames(new_data_raw), ";"), "[[", 3)
new_data_raw$sequence<- sapply(strsplit(rownames(new_data_raw), ";"), "[[", 4)
new_data_raw$id_newraw <- paste(new_data_raw$uniprot, new_data_raw$site, new_data_raw$sequence, sep = ";")
new_data_raw$id2_newraw <- paste(new_data_raw$uniprot, new_data_raw$site, sep = ";")

length(intersect(old_data_raw$id_oldraw, new_data_raw$id_newraw))
length(intersect(old_data_raw$id2_oldraw, new_data_raw$id2_newraw))

length(intersect(old_data$id_old, new_data$id_new))
length(intersect(old_data$id2_old, new_data$id2_new))


common_ids<- intersect(old_data$id_old, new_data$id_new)
common_uniprotsites <-intersect(old_data$id2_old, new_data$id2_new)

common_ids_raw<- intersect(old_data_raw$id_oldraw, new_data_raw$id_newraw)
common_uniprotsites_raw <-intersect(old_data_raw$id2_oldraw, new_data_raw$id2_newraw)



venn.diagram(
  x = list(old = old_data$id2_old, new = new_data$id2_new),
  main = "Overlap after filtering",
  sub = "Uniprot;Site",
    category.names = c("old data" , "new data"),
  filename = 'reanalysis_141122/venn_total2.tiff',
  height = 3000, 
  width = 3300,
  output = FALSE,
  print.mode = c("raw", "percent"),
  resolution = 500, 
  imagetype="tiff" ,
  compression = "lzw",
  disable.logging = TRUE, 
  col=c("#440154ff", '#21908dff'),
  fill = c("#440154ff", "#21908dff"),
  main.cex = 1.6,
  sub.cex = 1.5,
  alpha = 0.2,
  cex = 1.2,
  cat.cex = 1.1,
  fontfamily = "sans",
  cat.fontfamily = "sans",
  main.fontfamily = "sans",
  sub.fontfamily	= "sans"
)


##CORRELATIONS###################################################################################
##filter top tables by the common ids

old_data_top_common_rows <- which(old_data$id_old %in% common_ids)
old_data_top_common <- old_data[old_data_top_common_rows, ]
old_data_top_common$id <- row.names(old_data_top_common)


new_data_top_common_rows <- which(new_data$id_new %in% common_ids)
new_data_top_common <- new_data[new_data_top_common_rows, ]
new_data_top_common$id <- row.names(new_data_top_common)

new_data_toCtrl0_top_common_rows <- which(new_data_toCtrl0$id_new %in% common_ids)
new_data_toCtrl0_top_common <- new_data_toCtrl0[new_data_toCtrl0_top_common_rows, ]
new_data_toCtrl0_top_common$id <- row.names(new_data_toCtrl0_top_common)

new_old_merged <- merge(old_data_top_common, new_data_top_common, by.x = "id_old", by.y = "id_new", all = TRUE)
new2_old_merged <- merge(old_data_top_common, new_data_toCtrl0_top_common, by.x = "id_old", by.y = "id_new", all = TRUE)


write.table(new_old_merged, file = "new_old_merged.txt", sep = "\t", row.names = F)
write.table(new_data_top_common, file = "new_data_top_common.txt", sep = "\t", row.names = T)
write.table(old_data_top_common, file = "old_data_top_common.txt", sep = "\t", row.names = T)
write.table(new_data_toCtrl0_top_common, file = "new_data_toCtrl0_top_common.txt", sep = "\t", row.names = T)

# > colnames(new_old_merged)
# [1] "id_old"           "uniprot.x"        "symbol.x"         "logFC.10.x"       "P.Value.10.x"    
# [6] "adj.P.Val.10.x"   "logFC.30"         "P.Value.30"       "adj.P.Val.30"     "logFC.60"        
# [11] "P.Value.60"       "adj.P.Val.60"     "logFC.300"        "P.Value.300"      "adj.P.Val.300"   
# [16] "logFC.600.x"      "P.Value.600.x"    "adj.P.Val.600.x"  "logFC.900"        "P.Value.900"     
# [21] "adj.P.Val.900"    "logFC.1800.x"     "P.Value.1800.x"   "adj.P.Val.1800.x" "site_old"        
# [26] "sequence_old"     "id2_old"          "id.x"             "uniprot.y"        "symbol.y"        
# [31] "logFC.10.y"       "P.Value.10.y"     "adj.P.Val.10.y"   "logFC.600.y"      "P.Value.600.y"   
# [36] "adj.P.Val.600.y"  "logFC.1800.y"     "P.Value.1800.y"   "adj.P.Val.1800.y" "uniprot_new"     
# [41] "site_new"         "sequence_new"     "id2_new"          "id.y"            

ggpubr::ggqqplot(new_old_merged$logFC.10.x)
shapiro.test(new_old_merged$logFC.10.x)
ggpubr::ggqqplot(new_old_merged$logFC.600.x)
shapiro.test(new_old_merged$logFC.600.x)
ggpubr::ggqqplot(new_old_merged$logFC.1800.x)
shapiro.test(new_old_merged$logFC.1800.x)


ggpubr::ggqqplot(new_old_merged$logFC.10.y)
shapiro.test(new_old_merged$logFC.10.y)
ggpubr::ggqqplot(new_old_merged$logFC.600.y)
shapiro.test(new_old_merged$logFC.600.y)
ggpubr::ggqqplot(new_old_merged$logFC.1800.y)
shapiro.test(new_old_merged$logFC.1800.y)


cor.test(new_old_merged$logFC.10.x, new_old_merged$logFC.10.y, method = "kendall")
cor.test(new_old_merged$logFC.600.x, new_old_merged$logFC.600.y, method = "kendall")
cor.test(new_old_merged$logFC.1800.x, new_old_merged$logFC.1800.y, method = "kendall")


a<-ggscatter(new_old_merged, x = "logFC.10.x", y = "logFC.10.y", 
          title = "t = 10 sec",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "kendall",
          xlab = "LogFC old", ylab = "logFC new")

b<-ggscatter(new_old_merged, x = "logFC.600.x", y = "logFC.600.y", 
          title = "t = 600 sec",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "kendall",
          xlab = "LogFC old", ylab = "logFC new")
c<-ggscatter(new_old_merged, x = "logFC.1800.x", y = "logFC.1800.y", 
          title = "t = 1800 sec",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "kendall",
          xlab = "LogFC old", ylab = "logFC new")

abc<- ggpubr::ggarrange(a, b, c, nrow = 3)

#compare old results to new to ctrl0

ggpubr::ggqqplot(new2_old_merged$logFC.10.x)
shapiro.test(new2_old_merged$logFC.10.x)
ggpubr::ggqqplot(new2_old_merged$logFC.600.x)
shapiro.test(new2_old_merged$logFC.600.x)
ggpubr::ggqqplot(new2_old_merged$logFC.1800.x)
shapiro.test(new2_old_merged$logFC.1800.x)


ggpubr::ggqqplot(new2_old_merged$logFC.10.y)
shapiro.test(new2_old_merged$logFC.10.y)
ggpubr::ggqqplot(new2_old_merged$logFC.600.y)
shapiro.test(new2_old_merged$logFC.600.y)
ggpubr::ggqqplot(new2_old_merged$logFC.1800.y)
shapiro.test(new2_old_merged$logFC.1800.y)


cor.test(new2_old_merged$logFC.10.x, new2_old_merged$logFC.10.y, method = "kendall")
cor.test(new2_old_merged$logFC.600.x, new2_old_merged$logFC.600.y, method = "kendall")
cor.test(new2_old_merged$logFC.1800.x, new2_old_merged$logFC.1800.y, method = "kendall")


d<-ggscatter(new2_old_merged, x = "logFC.10.x", y = "logFC.10.y", 
          title = "t = 10 sec",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "kendall",
          xlab = "LogFC old", ylab = "logFC new")

e<-ggscatter(new2_old_merged, x = "logFC.600.x", y = "logFC.600.y", 
          title = "t = 600 sec",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "kendall",
          xlab = "LogFC old", ylab = "logFC new")
f<-ggscatter(new2_old_merged, x = "logFC.1800.x", y = "logFC.1800.y", 
          title = "t = 1800 sec",
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "kendall",
          xlab = "LogFC old", ylab = "logFC new")

def<- ggpubr::ggarrange(d, e, f, nrow = 3)


ggpubr::ggarrange(abc, def, nrow = 1)


####significants

compare10s <- new_old_merged[, c(1:6, 31:33)]
compare10s.sign <- compare10s[compare10s[, "P.Value.10.x"] < 0.05 & compare10s[, "P.Value.10.y"] < 0.05,]

compare600s <- new_old_merged[, c(1:3,16:18, 34:36)]
compare600s.sign <- compare600s[compare600s[, "P.Value.600.x"] < 0.05 & compare600s[, "P.Value.600.y"] < 0.05,]

compare1800s <- new_old_merged[, c(1:3, 22:24, 37:39)]
compare1800s.sign <- compare1800s[compare1800s[, "P.Value.1800.x"] < 0.05 & compare1800s[, "P.Value.1800.y"] < 0.05,]



ggpubr::ggqqplot(compare10s.sign$logFC.10.x)
shapiro.test(compare10s.sign$logFC.10.x)
ggpubr::ggqqplot(compare600s.sign$logFC.600.x)
shapiro.test(compare600s.sign$logFC.600.x)
ggpubr::ggqqplot(compare1800s.sign$logFC.1800.x)
shapiro.test(compare1800s.sign$logFC.1800.x)


ggpubr::ggqqplot(compare10s.sign$logFC.10.y)
shapiro.test(compare10s.sign$logFC.10.y)
ggpubr::ggqqplot(compare600s.sign$logFC.600.y)
shapiro.test(compare600s.sign$logFC.600.y)
ggpubr::ggqqplot(compare1800s.sign$logFC.1800.y)
shapiro.test(compare1800s.sign$logFC.1800.y)


cor.test(compare10s.sign$logFC.10.x, compare10s.sign$logFC.10.y, method = "kendall")
cor.test(compare600s.sign$logFC.600.x, compare600s.sign$logFC.600.y, method = "kendall")
cor.test(compare1800s.sign$logFC.1800.x, compare1800s.sign$logFC.1800.y, method = "kendall")


g<-ggscatter(compare10s.sign, x = "logFC.10.x", y = "logFC.10.y", 
             title = "t = 10 sec",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = TRUE, cor.method = "kendall",
             xlab = "LogFC old", ylab = "logFC new")

h<-ggscatter(compare600s.sign, x = "logFC.600.x", y = "logFC.600.y", 
             title = "t = 600 sec",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = TRUE, cor.method = "kendall",
             xlab = "LogFC old", ylab = "logFC new")
i<-ggscatter(compare1800s.sign, x = "logFC.1800.x", y = "logFC.1800.y", 
             title = "t = 1800 sec",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = TRUE, cor.method = "kendall",
             xlab = "LogFC old", ylab = "logFC new")

ghi<- ggpubr::ggarrange(g, h, i, nrow = 3)








compare10s2 <- new2_old_merged[, c(1:6, 31:33)]
compare10s2.sign <- compare10s2[compare10s2[, "P.Value.10.x"] < 0.05 & compare10s2[, "P.Value.10.y"] < 0.05,]

compare600s2 <- new2_old_merged[, c(1:3,16:18, 34:36)]
compare600s2.sign <- compare600s2[compare600s2[, "P.Value.600.x"] < 0.05 & compare600s2[, "P.Value.600.y"] < 0.05,]

compare1800s2 <- new2_old_merged[, c(1:3, 22:24, 37:39)]
compare1800s2.sign <- compare1800s2[compare1800s2[, "P.Value.1800.x"] < 0.05 & compare1800s2[, "P.Value.1800.y"] < 0.05,]

ggpubr::ggqqplot(compare10s2.sign$logFC.10.x)
shapiro.test(compare10s2.sign$logFC.10.x)
ggpubr::ggqqplot(compare600s2.sign$logFC.600.x)
shapiro.test(compare600s2.sign$logFC.600.x)
ggpubr::ggqqplot(compare1800s2.sign$logFC.1800.x)
shapiro.test(compare1800s2.sign$logFC.1800.x)


ggpubr::ggqqplot(compare10s2.sign$logFC.10.y)
shapiro.test(compare10s2.sign$logFC.10.y)
ggpubr::ggqqplot(compare600s2.sign$logFC.600.y)
shapiro.test(compare600s2.sign$logFC.600.y)
ggpubr::ggqqplot(compare1800s2.sign$logFC.1800.y)
shapiro.test(compare1800s2.sign$logFC.1800.y)


cor.test(compare10s2.sign$logFC.10.x, compare10s2.sign$logFC.10.y, method = "kendall")
cor.test(compare600s2.sign$logFC.600.x, compare600s2.sign$logFC.600.y, method = "kendall")
cor.test(compare1800s2.sign$logFC.1800.x, compare1800s2.sign$logFC.1800.y, method = "kendall")



j<-ggscatter(compare10s2.sign, x = "logFC.10.x", y = "logFC.10.y", 
             title = "t = 10 sec",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = TRUE, cor.method = "kendall",
             xlab = "LogFC old", ylab = "logFC new")

k<-ggscatter(compare600s2.sign, x = "logFC.600.x", y = "logFC.600.y", 
             title = "t = 600 sec",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = TRUE, cor.method = "kendall",
             xlab = "LogFC old", ylab = "logFC new")
l<-ggscatter(compare1800s2.sign, x = "logFC.1800.x", y = "logFC.1800.y", 
             title = "t = 1800 sec",
             add = "reg.line", conf.int = TRUE, 
             cor.coef = TRUE, cor.method = "kendall",
             xlab = "LogFC old", ylab = "logFC new")

jkl<- ggpubr::ggarrange(j, k, l, nrow = 3)


ggpubr::ggarrange(ghi, jkl, nrow = 1)
