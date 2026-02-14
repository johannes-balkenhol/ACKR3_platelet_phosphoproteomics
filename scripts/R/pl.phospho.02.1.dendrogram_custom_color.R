# ============================================================
# Publication-Quality Dendrogram with Custom Colors
# Replace your plotQC(..., panel="dendrogram") calls with this
# ============================================================

library(dendextend)

# ============================================================
# FUNCTION: Create a colored dendrogram
# ============================================================
plot_colored_dendrogram <- function(data, ppe = NULL, title = "", filename = NULL) {
  
  # --- Clustering (same as typical plotQC) ---
  # data is test2 or test4: rows = phosphosites, columns = samples
  # We transpose so that samples are clustered
  d <- dist(t(data))
  hc <- hclust(d, method = "complete")  # change method if needed: "ward.D2", "average", etc.
  dend <- as.dendrogram(hc)
  
  # --- Get labels ---
  # These will be the column names of test2/test4 (e.g. "x0sek_Ctrl_01")
  labs <- labels(dend)
  
  # If ppe is provided, use its colnames for labels (like your original plotQC call)
  if (!is.null(ppe)) {
    # Match the dendrogram order to ppe columns
    # plotQC used: sapply(strsplit(colnames(ppe), "_"), "[[", 1)
    # This extracts the first element, e.g. "x0sek" from "x0sek_Ctrl_01"
    # But for coloring we still use the full column names
    labs <- labels(dend)
  }
  
  # --- Rename labels to be more intuitive ---
  # e.g. "x10sek_DMSO_04" -> "DMSO_10s_04"
  #       "x0sek_Ctrl_01"  -> "Control_0s_01"
  rename_label <- function(label) {
    # Extract timepoint
    tp <- sub("^x(\\d+)sek_.*", "\\1s", label)
    # Extract treatment
    tr <- ifelse(grepl("Ctrl", label), "Control",
          ifelse(grepl("CXCR7", label), "ACKR3",
          ifelse(grepl("DMSO", label), "DMSO", "Unknown")))
    # Extract replicate number
    rep <- sub(".*_(\\d+)$", "\\1", label)
    paste0(tr, "_", tp, "_", rep)
  }
  
  # Apply new labels to dendrogram
  new_labs <- sapply(labs, rename_label)
  labels(dend) <- new_labs
  
  # --- Extract timepoint and treatment from ORIGINAL column names ---
  get_timepoint <- function(label) {
    if (grepl("^x1800sek", label)) return("x1800sek")
    if (grepl("^x600sek", label))  return("x600sek")
    if (grepl("^x10sek", label))   return("x10sek")
    if (grepl("^x0sek", label))    return("x0sek")
    return(NA)
  }
  
  get_treatment <- function(label) {
    if (grepl("Ctrl", label))  return("Ctrl")
    if (grepl("CXCR7", label)) return("CXCR7")
    if (grepl("DMSO", label))  return("DMSO")
    return(NA)
  }
  
  # ==========================================================
  # COLOR PALETTE — EDIT THESE TO CHANGE COLORS
  # ==========================================================
  # Option 1: Color by TREATMENT (simple, 3 colors)
  # treatment_colors <- c(
  #   "Ctrl"  = "#4DAF4A",
  #   "DMSO"  = "#377EB8",
  #   "CXCR7" = "#E41A1C"
  # )
  # label_cols <- sapply(labs, function(l) treatment_colors[get_treatment(l)])
  # legend_colors <- treatment_colors
  # legend_title <- "Treatment"
  
  # Option 2: Color by TIMEPOINT x TREATMENT (detailed, 7 colors)
  combo_colors <- c(
    "x0sek_Ctrl"      = "#1B9E77",   # teal
    "x10sek_DMSO"     = "#D95F02",   # orange
    "x10sek_CXCR7"    = "#E6AB02",   # gold
    "x600sek_DMSO"    = "#7570B3",   # purple
    "x600sek_CXCR7"   = "#E41A1C",   # red
    "x1800sek_DMSO"   = "#377EB8",   # blue
    "x1800sek_CXCR7"  = "#66A61E"    # green
  )
  
  # Intuitive legend names (matching the new label format)
  combo_legend_names <- c(
    "x0sek_Ctrl"      = "Control  0s",
    "x10sek_DMSO"     = "DMSO  10s",
    "x10sek_CXCR7"    = "ACKR3  10s",
    "x600sek_DMSO"    = "DMSO  600s",
    "x600sek_CXCR7"   = "ACKR3  600s",
    "x1800sek_DMSO"   = "DMSO  1800s",
    "x1800sek_CXCR7"  = "ACKR3  1800s"
  )
  
  get_combo <- function(label) {
    tp <- get_timepoint(label)
    tr <- get_treatment(label)
    paste0(tp, "_", tr)
  }
  
  label_cols <- sapply(labs, function(l) {
    combo <- get_combo(l)
    col <- combo_colors[combo]
    if (is.na(col)) return("grey40")
    return(col)
  })
  
  legend_colors <- combo_colors
  legend_title <- "Group"
  # ==========================================================
  
  # --- Apply colors ---
  labels_colors(dend) <- label_cols
  dend <- set(dend, "labels_cex", 0.55)     # label size
  dend <- set(dend, "branches_lwd", 0.8)    # branch line width
  
  # --- Plot ---
  if (!is.null(filename)) {
    tiff(filename = filename,
         width = 14 * 300, height = 4.5 * 300,
         res = 300, compression = "lzw")
  }
  
  # Extra right margin for legend outside the plot
  par(mar = c(10, 4, 3, 10), xpd = TRUE)
  plot(dend, main = title, ylab = "Height")
  
  # --- Legend (outside plot, right side) ---
  # Only show groups that actually exist in the data
  present <- unique(sapply(labs, get_combo))
  present <- present[present %in% names(legend_colors)]
  # Sort legend logically: by timepoint then treatment
  present <- present[order(
    as.numeric(sub("x(\\d+)sek_.*", "\\1", present)),
    sub("x\\d+sek_(.*)", "\\1", present)
  )]
  legend("right",
         legend = combo_legend_names[present],
         col    = legend_colors[present],
         pch    = 15, cex = 0.65, bty = "n",
         title  = "Group",
         inset  = c(-0.12, 0))    # push legend outside plot area
  
  if (!is.null(filename)) dev.off()
}


# ============================================================
# USAGE — replace your plotQC calls with these:
# ============================================================

# Dendrogram after normalization
# Original: plotQC(test4, grps=grps, labels = sapply(strsplit(colnames(ppe), "_"), "[[",1), panel="dendrogram")
plot_colored_dendrogram(
  data     = test4,
  ppe      = ppe,
  title    = "Hierarchical Clustering - After Normalization",
  filename = "analysis/PCA/dendrogram_after_normalization.tiff"
)

# Dendrogram imputed & scaled
# Original: plotQC(test2, grps=grps, labels = sapply(strsplit(colnames(ppe), "_"), "[[",1), panel="dendrogram")
plot_colored_dendrogram(
  data     = test2,
  ppe      = ppe,
  title    = "Hierarchical Clustering - Imputed & Scaled",
  filename = "analysis/PCA/dendrogram_imputed_scaled.tiff"
)
