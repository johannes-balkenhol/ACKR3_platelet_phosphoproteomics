######### compare empirical

empirical_topall_old <- read.table("../phosphoproteom analysis 2/data/processed_data/empirical_topall.txt", 
                        sep="\t", header=TRUE, dec=",")


empirical.old <- paste(sapply(strsplit(empirical_topall_old$x, ";"), "[[", 2),
                    sapply(strsplit(empirical_topall_old$x, ";"), "[[", 3), sep = ";")

empirical.new <- paste(sapply(strsplit(empirical_topall, ";"), "[[", 2),
                    sapply(strsplit(empirical_topall, ";"), "[[", 3), sep = ";")




sel <- which(empirical.new %in% empirical.old)
empirical_match <- empirical_topall[sel]

sel <- which(empirical.old %in% empirical.new)
empirical_topall_old <- as.data.frame(empirical_topall_old)
empirical_match <- empirical_topall_old[sel,]

write.table(empirical_match, "data/processed_data/empirical_match.txt", sep="\t", , row.names=FALSE)



############ compare psites
top.all_old <- read.table("../phosphoproteom analysis 2/data/processed_data/top.all.txt", 
                        sep="\t", header=TRUE, dec=",")

top.all_old.sign <- top.all_old[top.all_old,]
top.all_old.sign <- top.all_old[top.all_old,]

top.all_old_names <- paste(sapply(strsplit(rownames(top.all_old), ";"), "[[", 2),
                    sapply(strsplit(rownames(top.all_old), ";"), "[[", 3), sep = ";")

top.all_names <- paste(sapply(strsplit(rownames(top.all), ";"), "[[", 2),
                    sapply(strsplit(rownames(top.all), ";"), "[[", 3), sep = ";")

which(top.all_names %in% top.all_old_names)
length(which(top.all_names %in% top.all_old_names))
length(top.all_names)
length(top.all_old_names)



#3787 length(top.all_names)
#4208 length(top.all_old_names)
#2456 overlap

############ compare significant psites
top.all_old.sign <- top.all_old[top.all_old$P.Value.600<0.05,]
top.all.sign <- top.all_old[top.all$P.Value.600<0.05,]

top.all_old_names <- paste(sapply(strsplit(rownames(top.all_old.sign), ";"), "[[", 2),
                    sapply(strsplit(rownames(top.all_old.sign), ";"), "[[", 3), sep = ";")

top.all_names <- paste(sapply(strsplit(rownames(top.all.sign), ";"), "[[", 2),
                    sapply(strsplit(rownames(top.all.sign), ";"), "[[", 3), sep = ";")

which(top.all_names %in% top.all_old_names)

length(which(top.all_names %in% top.all_old_names))
length(top.all_names)
length(top.all_old_names)