KSEA.Barplot <- function (KSData, PX, NetworKIN, NetworKIN.cutoff, m.cutoff, 
          p.cutoff, export) 
{
  if (length(grep(";", PX$Residue.Both)) == 0) {
    new = PX
    colnames(new)[c(2, 4)] = c("SUB_GENE", "SUB_MOD_RSD")
    new$log2FC = log2(abs(as.numeric(as.character(new$FC))))
    new = new[complete.cases(new$log2FC), ]
  }
  else {
    double = PX[grep(";", PX$Residue.Both), ]
    residues = as.character(double$Residue.Both)
    residues = as.matrix(residues, ncol = 1)
    split = strsplit(residues, split = ";")
    x = sapply(split, length)
    single = data.frame(Protein = rep(double$Protein, x), 
                        Gene = rep(double$Gene, x), Peptide = rep(double$Peptide, 
                                                                  x), Residue.Both = unlist(split), p = rep(double$p, 
                                                                                                            x), FC = rep(double$FC, x))
    new = PX[-grep(";", PX$Residue.Both), ]
    new = rbind(new, single)
    colnames(new)[c(2, 4)] = c("SUB_GENE", "SUB_MOD_RSD")
    new$log2FC = log2(abs(as.numeric(as.character(new$FC))))
    new = new[complete.cases(new$log2FC), ]
  }
  if (NetworKIN == TRUE) {
    KSData.filtered = KSData[grep("[a-z]", KSData$Source), 
    ]
    KSData.filtered = KSData.filtered[(KSData.filtered$networkin_score >= 
                                         NetworKIN.cutoff), ]
  }
  else {
    KSData.filtered = KSData[grep("PhosphoSitePlus", KSData$Source), 
    ]
  }
  KSData.dataset = merge(KSData.filtered, new)
  KSData.dataset = KSData.dataset[order(KSData.dataset$GENE), 
  ]
  KSData.dataset$Uniprot.noIsoform = sapply(KSData.dataset$KIN_ACC_ID, 
                                            function(x) unlist(strsplit(as.character(x), split = "-"))[1])
  KSData.dataset.abbrev = KSData.dataset[, c(5, 1, 2, 16:19, 
                                             14)]
  colnames(KSData.dataset.abbrev) = c("Kinase.Gene", "Substrate.Gene", 
                                      "Substrate.Mod", "Peptide", "p", "FC", "log2FC", "Source")
  KSData.dataset.abbrev = KSData.dataset.abbrev[order(KSData.dataset.abbrev$Kinase.Gene, 
                                                      KSData.dataset.abbrev$Substrate.Gene, KSData.dataset.abbrev$Substrate.Mod, 
                                                      KSData.dataset.abbrev$p), ]
  KSData.dataset.abbrev = aggregate(log2FC ~ Kinase.Gene + 
                                      Substrate.Gene + Substrate.Mod + Source, data = KSData.dataset.abbrev, 
                                    FUN = mean)
  KSData.dataset.abbrev = KSData.dataset.abbrev[order(KSData.dataset.abbrev$Kinase.Gene), 
  ]
  kinase.list = as.vector(KSData.dataset.abbrev$Kinase.Gene)
  kinase.list = as.matrix(table(kinase.list))
  Mean.FC = aggregate(log2FC ~ Kinase.Gene, data = KSData.dataset.abbrev, 
                      FUN = mean)
  Mean.FC = Mean.FC[order(Mean.FC[, 1]), ]
  Mean.FC$mS = Mean.FC[, 2]
  Mean.FC$Enrichment = Mean.FC$mS/abs(mean(new$log2FC, na.rm = T))
  Mean.FC$m = kinase.list
  Mean.FC$z.score = ((Mean.FC$mS - mean(new$log2FC, na.rm = T)) * 
                       sqrt(Mean.FC$m))/sd(new$log2FC, na.rm = T)
  Mean.FC$p.value = pnorm(-abs(Mean.FC$z.score))
  Mean.FC$FDR = p.adjust(Mean.FC$p.value, method = "fdr")
  Mean.FC.filtered = Mean.FC[(Mean.FC$m >= m.cutoff), -2]
  Mean.FC.filtered = Mean.FC.filtered[order(Mean.FC.filtered$z.score), 
  ]
  plot.height = length(Mean.FC.filtered$z.score)^0.55
  Mean.FC.filtered$color = "black"
  Mean.FC.filtered[(Mean.FC.filtered$p.value < p.cutoff) & 
                     (Mean.FC.filtered$z.score < 0), ncol(Mean.FC.filtered)] = "blue"
  Mean.FC.filtered[(Mean.FC.filtered$p.value < p.cutoff) & 
                     (Mean.FC.filtered$z.score > 0), ncol(Mean.FC.filtered)] = "red"
  if (export == "TRUE") {
    tiff("KSEA Bar Plot.tiff", width = 6 * 300, height = 300 * 
           plot.height, res = 300, pointsize = 13)
    par(mai = c(1, 1, 0.4, 0.4))
    barplot(as.numeric(Mean.FC.filtered$z.score), col = Mean.FC.filtered$color, 
            border = NA, xpd = F, cex.names = 0.8, cex.axis = 0.8, 
            xlab = "Kinase enrichment (z-score)", names.arg = Mean.FC.filtered$Kinase.Gene, 
            horiz = T, las = 1)
    dev.off()
  }
  else {
    par(mgp = c(3, 0.25, 0))
    barplot(as.numeric(Mean.FC.filtered$z.score), col = Mean.FC.filtered$color, 
            border = "black", xpd = T, cex.names = 0.8, cex.axis = 1.2, xlim = c(-3,3),
            #xlab = "Kinase enrichment (z-score)", 
            names.arg = Mean.FC.filtered$Kinase.Gene,
            cex.lab = 1,
            horiz = T, las = 1, xaxt = "n")
    axis(1, at = seq(-3, 3, by = 1), tcl = -0.1, cex.axis = 1.1) # Shorten tick length
    title(xlab = "Kinase enrichment (z-score)", cex.lab = 1.4, line = 1.4)
    

  }
}


KSEA.Barplot(KSData, PX, NetworKIN=FALSE, NetworKIN.cutoff=5, m.cutoff=5, p.cutoff=0.05, export=FALSE)
