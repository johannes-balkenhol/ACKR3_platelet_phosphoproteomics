
library(ggplot2)
library(cowplot)
library(gridExtra)


# 1. Load Phosphodata 1 ----
phospho1 <- read.delim("/Users/ozgeosmanoglu/Nextcloud/CXCR7_platelet_analysis/CXCR7_initial/phosphoproteomics/data/processed_data/top.all.txt", dec = ",")
phospho1 <- data.frame(phospho1)

# change rownames
uniprot <- sapply(strsplit(rownames(phospho1), ";"), `[`, 1)
name <- sapply(strsplit(rownames(phospho1), ";"), `[`, 2)
site <- sapply(strsplit(rownames(phospho1), ";"), `[`, 3)
namesite <- paste0(name, ";", site)
rownames(phospho1) <- namesite



# 2. Load Phosphodata 2 ----
phospho2 <- read.delim("/Users/ozgeosmanoglu/Nextcloud/CXCR7_platelet_analysis/CXCR7_validation/phosphoproteomics/data/processed_data/top.all.cxcr7.vs.0s.txt")
phospho2 <- data.frame(phospho2)
uniprot <- sapply(strsplit(rownames(phospho2), ";"), `[`, 1)
name <- sapply(strsplit(rownames(phospho2), ";"), `[`, 2)
site <- sapply(strsplit(rownames(phospho2), ";"), `[`, 3)
namesite <- paste0(name, ";", site)
rownames(phospho2) <- namesite

# 3. Get Core Phosphopeptide Set ----

CoreSet <- intersect(rownames(phospho1), rownames(phospho2))
CoreDf1 <- phospho1[CoreSet,]
CoreDf2 <- phospho2[CoreSet,]

# 4. Get Counts
phospho1_counts <- data.frame(dataset = rep("initial", 6),
                              value = c(nrow(phospho1[phospho1$logFC.10 > 0.5 & phospho1$adj.P.Val.10 < 0.05,]),
                                        nrow(phospho1[phospho1$logFC.600 > 0.5 & phospho1$adj.P.Val.600 < 0.05,]),
                                        nrow(phospho1[phospho1$logFC.1800 > 0.5 & phospho1$adj.P.Val.1800 < 0.05,]),
                                        -1*nrow(phospho1[phospho1$logFC.10 < -0.5 & phospho1$adj.P.Val.10 < 0.05,]),
                                        -1*nrow(phospho1[phospho1$logFC.600 < -0.5 & phospho1$adj.P.Val.600 < 0.05,]),
                                        -1*nrow(phospho1[phospho1$logFC.1800 < -0.5 & phospho1$adj.P.Val.1800 < 0.05,])),
                              condition = rep(c("10s", "600s", "1800s"), 2),
                              variable = c(rep("Upregulated", 3), rep("Downregulated", 3)))

phospho2_counts <- data.frame(dataset = rep("initial", 6),
                              value = c(nrow(phospho2[phospho2$logFC.10 > 0.5 & phospho2$adj.P.Val.10 < 0.05,]),
                                        nrow(phospho2[phospho2$logFC.600 > 0.5 & phospho2$adj.P.Val.600 < 0.05,]),
                                        nrow(phospho2[phospho2$logFC.1800 > 0.5 & phospho2$adj.P.Val.1800 < 0.05,]),
                                        -1*nrow(phospho2[phospho2$logFC.10 < -0.5 & phospho2$adj.P.Val.10 < 0.05,]),
                                        -1*nrow(phospho2[phospho2$logFC.600 < -0.5 & phospho2$adj.P.Val.600 < 0.05,]),
                                        -1*nrow(phospho2[phospho2$logFC.1800 < -0.5 & phospho2$adj.P.Val.1800 < 0.05,])),
                              condition = rep(c("10s", "600s", "1800s"), 2),
                              variable = c(rep("Upregulated", 3), rep("Downregulated", 3)))

CoreDf1_counts <- data.frame(dataset = rep("core", 6),
                             value = c(nrow(CoreDf1[CoreDf1$logFC.10 > 0.5 & CoreDf1$adj.P.Val.10 < 0.05,]),
                                       nrow(CoreDf1[CoreDf1$logFC.600 > 0.5 & CoreDf1$adj.P.Val.600 < 0.05,]),
                                       nrow(CoreDf1[CoreDf1$logFC.1800 > 0.5 & CoreDf1$adj.P.Val.1800 < 0.05,]),
                                       -1*nrow(CoreDf1[CoreDf1$logFC.10 < -0.5 & CoreDf1$adj.P.Val.10 < 0.05,]),
                                       -1*nrow(CoreDf1[CoreDf1$logFC.600 < -0.5 & CoreDf1$adj.P.Val.600 < 0.05,]),
                                       -1*nrow(CoreDf1[CoreDf1$logFC.1800 < -0.5 & CoreDf1$adj.P.Val.1800 < 0.05,])),
                             condition = rep(c("10s", "600s", "1800s"), 2),
                             variable = c(rep("Upregulated", 3), rep("Downregulated", 3)))


CoreDf2_counts <- data.frame(dataset = rep("core", 6),
                             value = c(nrow(CoreDf2[CoreDf2$logFC.10 > 0.5 & CoreDf2$adj.P.Val.10 < 0.05,]),
                                       nrow(CoreDf2[CoreDf2$logFC.600 > 0.5 & CoreDf2$adj.P.Val.600 < 0.05,]),
                                       nrow(CoreDf2[CoreDf2$logFC.1800 > 0.5 & CoreDf2$adj.P.Val.1800 < 0.05,]),
                                       -1*nrow(CoreDf2[CoreDf2$logFC.10 < -0.5 & CoreDf2$adj.P.Val.10 < 0.05,]),
                                       -1*nrow(CoreDf2[CoreDf2$logFC.600 < -0.5 & CoreDf2$adj.P.Val.600 < 0.05,]),
                                       -1*nrow(CoreDf2[CoreDf2$logFC.1800 < -0.5 & CoreDf2$adj.P.Val.1800 < 0.05,])),
                             condition = rep(c("10s", "600s", "1800s"), 2),
                             variable = c(rep("Upregulated", 3), rep("Downregulated", 3)))

phospho1Combined_counts <- rbind(CoreDf1_counts, phospho1_counts)
phospho2Combined_counts <- rbind(phospho2_counts, CoreDf2_counts)

phospho1Combined_counts$condition <- factor(phospho1Combined_counts$condition, 
                                            levels = c("10s", "600s", "1800s"))  # Replace with your desired order
phospho2Combined_counts$condition <- factor(phospho2Combined_counts$condition, 
                                            levels = c("10s", "600s", "1800s"))  # Replace with your desired order


# Barplot for Initial Data ----
InitialPlot <- ggplot(phospho1Combined_counts, aes(x = condition, y = value)) +
  # Full-width bars for 'initial'
  geom_bar(
    data = subset(phospho1Combined_counts, dataset == "initial"),
    aes(fill = variable),
    color = NA,
    stat = "identity",
    alpha = 0.4,
    width = 0.85 # Full-width for 'initial'
  ) +
  # Narrower bars for 'core'
  geom_bar(
    data = subset(phospho1Combined_counts, dataset == "core"),
    aes(fill = variable),
    stat = "identity",
    color = NA,
    width = 0.7 # Narrower for 'core'
  ) +
  # Horizontal reference line
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.25) +
  # Use a light theme
  theme_cowplot() +
  ylim(-400, 400) +
  # Set axis labels
  labs(
    x = "Time",
    y = "Value",
    fill = "Variable",
    title = "Initial Dataset"
  ) +
  # Adjust fill colors
  scale_fill_manual(values = c("Upregulated" = "tomato", "Downregulated" = "blue"))+
  theme(legend.position="none",
        plot.title = element_text(size = 14))

# Barplot for Validation Data ----
ValidationPlot <- ggplot(phospho2Combined_counts, aes(x = condition, y = value)) +
  # Full-width bars for 'initial'
  geom_bar(
    data = subset(phospho2Combined_counts, dataset == "initial"),
    aes(fill = variable),
    color = NA,
    stat = "identity",
    alpha = 0.4,
    width = 0.85 # Full-width for 'initial'
  ) +
  # Narrower bars for 'core'
  geom_bar(
    data = subset(phospho2Combined_counts, dataset == "core"),
    aes(fill = variable),
    stat = "identity",
    color = NA,
    width = 0.7 # Narrower for 'core'
  ) +
  # Horizontal reference line
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.25) +
  # Use a light theme
  theme_cowplot() +
  ylim(-400, 400) +
  # Set axis labels
  labs(
    x = "Time",
    y = "Value",
    fill = "Variable",
    title = "Validation Dataset"
  ) +
  # Adjust fill colors
  scale_fill_manual(values = c("Upregulated" = "tomato", "Downregulated" = "blue"))+
  theme(legend.position="none",
        plot.title = element_text(size = 14))

combFig <- grid.arrange(InitialPlot, ValidationPlot, ncol =2)

ggsave(plot = InitialPlot,
       filename = "/Users/ozgeosmanoglu/Nextcloud/CXCR7_platelet_analysis/analysis/DEGbarplotInitial.tiff", 
       device = "tiff", width = 4, height = 6, dpi = 300, bg = "white")

ggsave(plot = ValidationPlot,
       filename = "/Users/ozgeosmanoglu/Nextcloud/CXCR7_platelet_analysis/analysis/DEGbarplotValidation.tiff", 
       device = "tiff", width = 4, height = 6, dpi = 300, bg = "white")

ggsave(plot = combFig,
       filename = "/Users/ozgeosmanoglu/Nextcloud/CXCR7_platelet_analysis/analysis/DEGbarplots.tiff", 
       device = "tiff", width = 8, height = 6, dpi = 300)
   



