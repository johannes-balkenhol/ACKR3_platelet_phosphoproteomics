
devtools::install_github("ByrumLab/proteoDA", 
                         dependencies = TRUE, 
                         build_vignettes = TRUE)

library(proteoDA)
library(ComplexHeatmap)

# Split input data into protein intensity data and annotation data
#intensity_data <- input_data[,5:21] # select columns 5 to 21 --> raw_abundance
intensity_data <- raw_abundance
annotation_data <-  peptide_info
annotation_data$uniprot_id <- annotation_data$Accession
#annotation_data <- input_data[,1:4] # select columns 1 to 4  --> peptide_info
# Match up row names of metadata with column names of data

## define header
donor_nr <- gsub("Donor", "",sapply(strsplit(colnames(raw_abundance), "_"), "[[", 3))
time_point <- sapply(strsplit(colnames(raw_abundance), "_"), "[[", 1)
timepoint_fac <- as.factor(time_point)
condition <- sapply(strsplit(colnames(raw_abundance), "_"), "[[", 2)

sample_metadata <- data.frame(colnames(raw_abundance), condition, donor_nr, paste0(condition, "_", time_point))
colnames(sample_metadata) <- c("data_column_name", "sample", "batch", "group")
rownames(sample_metadata) <- sample_metadata$data_column_name

# Assemble into DAList
raw <- DAList(data = intensity_data,
              annotation = annotation_data,
              metadata = sample_metadata)

# Filter out unneeded samples and proteins with too much missing data
filtered <- raw |>
  #filter_samples(group != "Pool") |>
  zero_to_missing() |>
  filter_proteins_by_proportion(min_prop = 0.66,
                                grouping_column = "group")
# Make the normalization report
write_norm_report(filtered,
                  grouping_column = "group", output_dir = "../analysis", filename = "normalization_report_on_raw.pdf")

# Normalize
normalized <- normalize_data(filtered, 
                             norm_method = "log2")

# Make the quality control report
write_qc_report(normalized,
                color_column = "group", output_dir = "../analysis", filename = "QC_report_raw(log2).pdf")

# Turn metadata column into a factor with desired levels
normalized$metadata$group <- factor(normalized$metadata$group, 
                                    levels = c("Ctrl_0000", "CXCR7_0010", "CXCR7_0600", 
                                               "CXCR7_1800", "DMSO_0010", "DMSO_0600", "DMSO_1800"))

# Add a statistical design, fit the model, and extract results
final <- normalized |>
  add_design(design_formula = ~ group) |>
  fit_limma_model() |>
  extract_DA_results()



# Export results
write_limma_tables(final, output_dir = "../analysis")
write_limma_plots(final,
                  grouping_column = "group", output_dir = "../analysis")




colorBatch = function(batch){
  batchCol =  unlist(ifelse(length(unique(batch)) == 1, rainbow(1, start = 0.5), list(rainbow(length(unique(batch))))))
  names(batchCol) = sort(unique(batch))
  return(batchCol)
}

colorGroup = function(groups){
  groupCol = viridis::viridis(length(unique(groups)))
  names(groupCol) = sort(unique(groups))
  return(groupCol)
}

heatmapCorr = function(data, groups, batch, sampleLabels){
  cor_mat <- cor(data, use = "pairwise.complete.obs")
  
  ColAnn <- HeatmapAnnotation(Sample = groups, 
                              col = list(Sample = colorGroup(groups)), 
                              annotation_legend_param = list(
                                Sample = list(
                                  title = "Sample",
                                  at = sort(unique(groups)),
                                  labels = paste("Group", sort(unique(groups)))
                                )
                              ))
  RowAnn <- rowAnnotation(Batch = batch, 
                          col = list(Batch = colorBatch(batch)), 
                          annotation_legend_param = list(
                            Batch = list(
                              title = "Batch",
                              at = sort(unique(batch)),
                              labels = paste("Batch", sort(unique(batch)))
                            )
                          ))
  hm_corr = Heatmap(cor_mat, 
                    col = circlize::colorRamp2(seq(min(cor_mat), 1, ((1 - min(cor_mat))/7)),RColorBrewer::brewer.pal(8, "Reds")), 
                    heatmap_legend_param = list(color_bar = "continuous", 
                                                legend_direction = "horizontal", 
                                                legend_width = unit(5, "cm"), 
                                                title_position = "topcenter"), 
                    name = "Pearson correlation", 
                    column_names_gp = gpar(fontsize = 12), 
                    row_names_gp = gpar(fontsize = 12), 
                    top_annotation = ColAnn, 
                    left_annotation = RowAnn,
                    column_labels = sampleLabels, row_labels = sampleLabels
  )
  
  draw(hm_corr, heatmap_legend_side = "top")
}

##prep data
data <- as.data.frame(raw_abundance)

#or
data <- as.data.frame(scaled_data)

##or
data <- as.data.frame(norm_intensity)


#continue
donor_nr <- sapply(strsplit(colnames(data), "_"), "[[", 3)
time_point <- sapply(strsplit(colnames(data), "_"), "[[", 1)
timepoint_fac <- as.factor(time_point)
condition <- sapply(strsplit(colnames(data), "_"), "[[", 2)
colnames(data) <- paste(time_point, condition, donor_nr, sep = "_")

#prep metadata
metadata <-  data.frame(donor_nr, timepoint_fac, condition, colnames(data))
colnames(metadata) <- c("donor", "time", "condition", "sample")
metadata$group <- paste0(metadata$time, "_", metadata$condition)
rownames(metadata) <- colnames(data)

#prep heatmap inputs
groups <- metadata$group
batch <- metadata$donor
sampleLabels <- metadata$sample

tiff("../analysis/PCA/correlation_scaled_style2.tiff", 
     width = 20, height = 12, units = 'in', res = 600, compression = "lzw")
heatmapCorr(data, groups, batch, sampleLabels)
dev.off()


