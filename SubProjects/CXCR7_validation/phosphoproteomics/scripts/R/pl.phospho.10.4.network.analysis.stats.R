# Define the relevant columns for comparison
adj_p_val_columns <- c("adj.P.Val_10_cxcr7.vs.0s", "adj.P.Val_600_cxcr7.vs.0s", "adj.P.Val_1800_cxcr7.vs.0s",
                       "adj.P.Val_10_dmso.vs.0s", "adj.P.Val_600_dmso.vs.0s", "adj.P.Val_1800_dmso.vs.0s")

# Function to count unique uniprot and psite values based on a condition
count_unique_values <- function(df, columns) {
  result <- data.frame(Column = character(), Dataframe = character(), Unique_Uniprot = integer(), Unique_Psite = integer())
  
  for (column in columns) {
    if (column %in% colnames(df)) {
      filtered_data <- df %>%
        filter(!!sym(column) < 0.05)
      
      count_unique_uniprot <- n_distinct(filtered_data$uniprot)
      count_unique_psite <- n_distinct(filtered_data$psite)
      
      result <- rbind(result, data.frame(Column = column, Dataframe = deparse(substitute(df)), 
                                         Unique_Uniprot = count_unique_uniprot, Unique_Psite = count_unique_psite))
    }
  }
  
  return(result)
}

# Apply the function to all relevant dataframes
results_init_cxcr7 <- count_unique_values(top.all.init.cxcr7.vs.0s, adj_p_val_columns)
results_val_cxcr7 <- count_unique_values(top.all.val.cxcr7.vs.0s, adj_p_val_columns)
results_val_dmso <- count_unique_values(top.all.val.dmso.vs.0s, adj_p_val_columns)

# Combine the results into a comparative table
comparative_results <- rbind(results_init_cxcr7, results_val_cxcr7, results_val_dmso)

# Display the comparative table
comparative_results
