### not working properly
### Özge did it manually

df <- raw_abundance2

df$uniprot_id <- sapply(strsplit(rownames(ppe0), ";"), "[[", 1)
df$sequence <- sapply(strsplit(rownames(ppe0), ";"), "[[", 4)
df$psite <- sapply(strsplit(rownames(ppe0), ";"), "[[", 3)


# Find duplicated combinations of psite and sequence
duplicated_combinations <- duplicated(df[, c("psites", "sequence")]) |
                           duplicated(df[, c("psites", "sequence")], fromLast = TRUE)

# Remove rows with duplicated combinations
cleaned_df <- df[!duplicated_combinations, ]


# Output the cleaned data frame
print(cleaned_df)