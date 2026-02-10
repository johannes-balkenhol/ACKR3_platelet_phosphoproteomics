# Johannes Balkenhol, MIT License

#### Script Overview ####
# This script processes phosphoproteomics and proteomics data, generates Venn diagrams to compare different conditions,
# and prepares for automation.

#### Inputs ####
# 1. Phosphoproteomics data files (e.g., top.10.dmso.vs.cxcr7.txt).
# 2. Proteomics data files (e.g., top.all.txt).
# 3. Folder paths for validation and initial experiment data.

#### Outputs ####
# 1. Processed and merged data frames (e.g., top.all.val.cxcr7.vs.0s).
# 2. Proteomics table (e.g. top.all.prot.val)
# 3. Venn diagrams comparing different conditions.

#### Required Packages ####
if (!require("VennDiagram")) install.packages("VennDiagram")
if (!require("dplyr")) install.packages("dplyr")

# Load necessary libraries
library(VennDiagram)
library(dplyr)

#### Data Preparation and Processing ####

# Function to read data, adjust structure, and set rownames based on the 'id' column
read_and_set_rowname <- function(folder_path, file_name, is_init = FALSE) {
  file_path <- file.path(folder_path, file_name)
  if (file.exists(file_path)) {
    data <- read.table(file_path, header = TRUE, sep = "\t")
    
    if (is_init && !("psite" %in% colnames(data))) {
      # Preprocess init files: add 'psite' column
      data$psite <- sapply(strsplit(data$id, ";"), function(x) x[3])
    }
    
    rownames(data) <- data$id
    return(data)
  } else {
    return(NULL)  # Return NULL if the file doesn't exist
  }
}

# Function to collapse rows by New ID
collapse_by_id <- function(df) {
  if (!is.null(df)) {
    df$NewID <- sapply(strsplit(rownames(df), ";"), function(x) paste(x[1:3], collapse = ";"))
    collapsed_df <- df %>%
      group_by(NewID) %>%
      slice_max(abs(logFC), with_ties = FALSE) %>%
      ungroup() %>%
      as.data.frame()
    rownames(collapsed_df) <- collapsed_df$NewID
    collapsed_df <- collapsed_df %>%
      dplyr::select(-NewID)
    return(collapsed_df)
  } else {
    return(NULL)
  }
}

# Function to rename columns with the time and comparison tag
rename_columns <- function(df, time, comparison) {
  if (!is.null(df)) {
    colnames(df)[colnames(df) %in% c("AveExpr", "logFC", "P.Value", "adj.P.Val", "t", "B")] <- 
      paste(colnames(df)[colnames(df) %in% c("AveExpr", "logFC", "P.Value", "adj.P.Val", "t", "B")], 
            time, comparison, sep = "_")
    return(df)
  } else {
    return(NULL)
  }
}

# Function to merge and reorder dataframes
merge_and_reorder <- function(dataframes, output_name) {
  # Filter out NULL dataframes
  dataframes <- Filter(Negate(is.null), dataframes)
  
  if (length(dataframes) < 2) {
    stop("At least two valid dataframes are required for merging.")
  }
  
  # Step 1: Merge the dataframes on the common columns
  merged_df <- dataframes[[1]]
  for (i in 2:length(dataframes)) {
    merged_df <- merged_df %>%
      inner_join(dataframes[[i]], by = c("uniprot", "symbol", "psite", "id"))
  }
  
  # Step 2: Set rownames based on uniprot, symbol, psite
  merged_df$NewRowNames <- paste(merged_df$uniprot, merged_df$symbol, merged_df$psite, sep = ";")
  rownames(merged_df) <- merged_df$NewRowNames
  merged_df$NewRowNames <- NULL
  
  # Step 3: Reorder columns to place uniprot, symbol, psite, id at the beginning
  all_columns <- colnames(merged_df)
  desired_order <- c("uniprot", "symbol", "psite", "id", setdiff(all_columns, c("uniprot", "symbol", "psite", "id")))
  merged_df <- merged_df[, desired_order]
  
  # Step 4: Assign the merged and reordered dataframe to the specified output name
  assign(output_name, merged_df, envir = .GlobalEnv)
}

# Main function to handle the entire process
process_data <- function(folder_path, file_names, experiment_type) {
  is_init <- experiment_type == "init"
  
  # Filter the file_names to only include existing files
  existing_files <- file_names[file.exists(file.path(folder_path, file_names))]
  
  # Extract comparisons and times from the existing files
  comparisons <- sapply(existing_files, function(fn) {
    strsplit(fn, split = "\\.txt")[[1]][1] %>%
      sub("top\\.[0-9]+\\.", "", .)
  })

  times <- sapply(existing_files, function(fn) {
    sub("top\\.", "", fn) %>%
      sub("\\..*", "", .) %>%
      gsub("[^0-9]", "", .)
  })
  
  collapsed_data <- list()
  for (i in 1:length(existing_files)) {
    data <- read_and_set_rowname(folder_path, existing_files[i], is_init = is_init)
    collapsed <- collapse_by_id(data)
    renamed <- rename_columns(collapsed, times[i], comparisons[i])
    collapsed_data[[paste(comparisons[i], times[i], sep = "_")]] <- renamed
  }
  
  unique_comparisons <- unique(comparisons)
  for (comp in unique_comparisons) {
    relevant_data <- collapsed_data[grep(comp, names(collapsed_data))]
    output_name <- paste0("top.all.", experiment_type, ".", comp)
    merge_and_reorder(relevant_data, output_name)
  }
}

# Example usage for validation (val) and initial (init) experiments
folder_path_phos_val <- "../../../phosphoproteomics/data/processed_data/"
folder_path_phos_init <- "../../../../CXCR7_initial/phosphoproteomics/data/processed_data/"

file_names <- c("top.10.dmso.vs.cxcr7.txt", "top.600.dmso.vs.cxcr7.txt", "top.1800.dmso.vs.cxcr7.txt",
                "top.10.cxcr7.vs.0s.txt", "top.600.cxcr7.vs.0s.txt", "top.1800.cxcr7.vs.0s.txt",
                "top.10.dmso.vs.0s.txt", "top.600.dmso.vs.0s.txt", "top.1800.dmso.vs.0s.txt")

# Process data for 'val' (validation) and 'init' (initial) experiments
process_data(folder_path_phos_val, file_names, "val")
process_data(folder_path_phos_init, file_names, "init")

#### Load Proteomics Data and Prepare Files ####

# Define folder paths for proteomics data
folder_path_prot_val <- "../../../proteomics/data/processed_data/"
folder_path_prot_init <- "../../../../CXCR7_initial/proteomics/data/processed_data/"

# Define the file name
file_name <- "top.all.txt"

# Load proteomics data for validation (val) and initial (init) experiments
top.all.prot.val <- read.table(file.path(folder_path_prot_val, file_name), header = TRUE, sep = "\t", row.names = 1)
top.all.prot.init <- read.table(file.path(folder_path_prot_init, file_name), header = TRUE, sep = "\t", row.names = 1)

# Display the first few rows of each dataframe to confirm successful loading
head(top.all.prot.val)
head(top.all.prot.init)

#### Analysis: Venn Diagrams ####

# Venn Diagram for Phosphoproteomics: Init vs Val CXCR7 vs 0s
venn_diagram_phospho <- function(init_df, val_df, title) {
  init_rownames <- rownames(init_df)
  val_rownames <- rownames(val_df)
  
  venn.plot <- venn.diagram(
    x = list("Init CXCR7 vs 0s" = init_rownames, "Val CXCR7 vs 0s" = val_rownames),
    category.names = c("Init CXCR7 vs 0s", "Val CXCR7 vs 0s"),
    filename = NULL,
    output = TRUE,
    main = title
  )
  
  dev.new()
  grid.draw(venn.plot)
}

venn_diagram_phospho(top.all.init.cxcr7.vs.0s, top.all.val.cxcr7.vs.0s, "Venn Diagram of Overlap: Init vs Val CXCR7 vs 0s")

# Venn Diagram for Proteomics: Validation vs Initial
venn_diagram_proteomics <- function(prot_val, prot_init, title) {
  prot_val_rownames <- rownames(prot_val)
  prot_init_rownames <- rownames(prot_init)
  
  venn.plot <- venn.diagram(
    x = list("Validation Proteomics" = prot_val_rownames, "Initial Proteomics" = prot_init_rownames),
    category.names = c("Validation Proteomics", "Initial Proteomics"),
    filename = NULL,
    output = TRUE,
    main = title
  )
  
  dev.new()
  grid.draw(venn.plot)
}

venn_diagram_proteomics(top.all.prot.val, top.all.prot.init, "Venn Diagram of Proteomics: Validation vs Initial")

# Venn Diagram for Proteomics and Phosphoproteomics: All comparisons
venn_diagram_all <- function(prot_val, prot_init, phos_init, phos_val, title) {
  prot_val_rownames <- rownames(prot_val)
  prot_init_rownames <- rownames(prot_init)
  init_phos_uniprot <- phos_init$uniprot
  val_phos_uniprot <- phos_val$uniprot
  
  venn.plot <- venn.diagram(
    x = list(
      "Validation Proteomics" = prot_val_rownames,
      "Initial Proteomics" = prot_init_rownames,
      "Init Phosphoproteomics CXCR7 vs 0s" = init_phos_uniprot,
      "Val Phosphoproteomics CXCR7 vs 0s" = val_phos_uniprot
    ),
    category.names = c("Validation Proteomics", "Initial Proteomics", "Init Phospho CXCR7 vs 0s", "Val Phospho CXCR7 vs 0s"),
    filename = NULL,
    output = TRUE,
    main = title
  )
  
  dev.new()
  grid.draw(venn.plot)
}

venn_diagram_all(top.all.prot.val, top.all.prot.init, top.all.init.cxcr7.vs.0s, top.all.val.cxcr7.vs.0s, "Venn Diagram of Uniprot ID Overlap")

# Venn Diagram for Phosphoproteomics: Init CXCR7 vs 0s, Val CXCR7 vs 0s, Val DMSO vs 0s
venn_diagram_phospho_comparison <- function(init_df, val_df, val_df2, title) {
  init_rownames <- rownames(init_df)
  val_rownames <- rownames(val_df)
  val2_rownames <- rownames(val_df2)
  
  venn.plot <- venn.diagram(
    x = list(
      "Init CXCR7 vs 0s" = init_rownames, 
      "Val CXCR7 vs 0s" = val_rownames, 
      "Val DMSO vs 0s" = val2_rownames
    ),
    category.names = c("Init CXCR7 vs 0s", "Val CXCR7 vs 0s", "Val DMSO vs 0s"),
    filename = NULL,
    output = TRUE,
    main = title
  )
  
  dev.new()
  grid.draw(venn.plot)
}

venn_diagram_phospho_comparison(top.all.init.cxcr7.vs.0s, top.all.val.cxcr7.vs.0s, top.all.val.dmso.vs.0s, "Venn Diagram of Overlap: (Init CXCR7 vs 0s) vs (Val CXCR7 vs 0s) vs (Val DMSO vs 0s)")
