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
      select(-NewID)
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
  # Determine if this is the "init" experiment for preprocessing
  is_init <- experiment_type == "init"
  
  # Filter the file_names to only include existing files
  existing_files <- file_names[file.exists(file.path(folder_path, file_names))]
  
  # Extract comparisons and times from the existing files
  comparisons <- sapply(existing_files, function(fn) {
    # Extract the comparison part (e.g., dmso.vs.cxcr7)
    strsplit(fn, split = "\\.txt")[[1]][1] %>%
      sub("top\\.[0-9]+\\.", "", .)  # Remove the "top." and time number part
  })

  times <- sapply(existing_files, function(fn) {
    # Extract the time part (e.g., 10, 600, 1800)
    sub("top\\.", "", fn) %>%
      sub("\\..*", "", .) %>%
      gsub("[^0-9]", "", .)
  })
  
  # Load and process each file
  collapsed_data <- list()
  for (i in 1:length(existing_files)) {
    data <- read_and_set_rowname(folder_path, existing_files[i], is_init = is_init)
    collapsed <- collapse_by_id(data)
    renamed <- rename_columns(collapsed, times[i], comparisons[i])
    # Store each dataframe with a unique key, e.g., "dmso.vs.cxcr7_10"
    collapsed_data[[paste(comparisons[i], times[i], sep = "_")]] <- renamed
  }
  
  # Merge and reorder for each comparison
  unique_comparisons <- unique(comparisons)
  print(unique_comparisons)
  for (comp in unique_comparisons) {
    relevant_data <- collapsed_data[grep(comp, names(collapsed_data))]

    output_name <- paste0("top.all.", experiment_type, ".", comp)
    merge_and_reorder(relevant_data, output_name)
  }
}

# Example usage for validation (val) and initial (init) experiments
folder_path_phos_val <- "../../../phosphoproteomics/data/processed_data/"
folder_path_phos_init <- "../../../../CXCR7_initial/phosphoproteomics/data/processed_data/"

file_names <- c("top.10.dmso.vs.cxcr7.txt", "top.600.dmso.vs.cxcr7.txt", "top.1800.dmso.vs.cxcr7.1800.txt",
                "top.10.cxcr7.vs.0s.txt", "top.600.cxcr7.vs.0s.txt", "top.1800.cxcr7.vs.0s.txt",
                "top.10.dmso.vs.0s.txt", "top.600.dmso.vs.0s.txt", "top.1800.dmso.vs.0s.txt")

# Process data for 'val' (validation) and 'init' (initial) experiments
process_data(folder_path_phos_val, file_names, "val")
process_data(folder_path_phos_init, file_names, "init")
