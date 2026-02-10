#################################################
##### prepare collpased psite tables for gsea
##### prepare collpased psite tables for gsea
input = list(top.10, top.600, top.1800, 
top.10.dmso.vs.0s, top.600.dmso.vs.0s, top.1800.dmso.vs.0s,
top.10.cxcr7.vs.0s, top.600.cxcr7.vs.0s, top.1800.cxcr7.vs.0s)

names_input = c("10", "600", "1800", 
"10.dmso.vs.0s", "600.dmso.vs.0s", "1800.dmso.vs.0s",
"10.cxcr7.vs.0s", "600.cxcr7.vs.0s", "1800.cxcr7.vs.0s")

for (i in 1:length(input)) {
    #df_pl = top.10.cxcr7.vs.0s
    df_pl = input[[i]]
    ## get rid of site without annotation
    #df_pl = df_pl[-which(df_pl$psite == ""),]
    df_pl = df_pl[, c(7,8,2,1,5,9,10)]

    colnames(df_pl)<- c("uniprot_id", "name", "Average", "logFC", "PValue", "PSite","ID")

    ## collpase same site by mean (sequence sometime slightly different)
    df_pl2 <- plyr::ddply(df_pl, .(uniprot_id, PSite), dplyr::summarise,
                name = max(name),
                Average = mean(Average),
                logFC = mean(logFC),
                PValue = mean(PValue),
                ID = max(ID)
                )

    rownames(df_pl2) <- df_pl2$ID
    df_pl2 <- df_pl2[c(1,3,4,5,6,2)]

    assign(paste0("top.filter.", names_input[[i]]), df_pl2)

    ## collpase phosphosites to the protein
    top.collapse <- ddply(df_pl, .(uniprot_id, name), summarise,
                  Average = max(Average),
                  minlogFC = min(logFC),
                  maxlogFC = max(logFC),
                  AbslogFC = max(abs(logFC)),
                  meanlogFC = mean(logFC),
                  PValue = min(PValue))

    top.collapse$logFC <- 0
    top.collapse$logFC[abs(top.collapse[, "minlogFC"]) < abs(top.collapse[, "maxlogFC"])] = 
                                  top.collapse$maxlogFC[abs(top.collapse[, "minlogFC"]) < abs(top.collapse[, "maxlogFC"])]

    top.collapse$logFC[abs(top.collapse[, "minlogFC"]) > abs(top.collapse[, "maxlogFC"])] = 
                                  top.collapse$minlogFC[abs(top.collapse[, "minlogFC"]) > abs(top.collapse[, "maxlogFC"])]

    top.collapse$logFC[abs(top.collapse[, "minlogFC"]) == abs(top.collapse[, "maxlogFC"])] = 
                                  top.collapse$minlogFC[abs(top.collapse[, "minlogFC"]) == abs(top.collapse[, "maxlogFC"])]

    ## optinal: set the mean
    #top.collapse$logFC = top.collapse$meanlogFC

    ## repeat for each input
    #colnames(top.collapse) <- c("Average","logFC","p.value","uniprot_id","name")
    rownames(top.collapse) <- paste(top.collapse$uniprot_id, top.collapse$name,sep=";")

    assign(paste0("top.collapse.", names_input[[i]]), top.collapse[,c(3,9,8,1,2)])
}


write.table(top.collapse, "../data/processed_data/top.collapse.txt", sep = "\t")

