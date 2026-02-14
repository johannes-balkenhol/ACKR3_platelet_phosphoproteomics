#################################################
##### filter and collpase the norm_intensity table for downstream analysis

##### filter
df <- as.data.frame(norm_intensity)

names_input = c("norm_intensity")

df$uniprot_id <- sapply(strsplit(rownames(df), ";"), "[[", 1)
df$PSite <- sapply(strsplit(rownames(df), ";"), "[[", 3)
df$ID <- rownames(df)

value_cols <- colnames(df)[startsWith(colnames(df), "x")]

#summary_df <- ddply(df, .(uniprot_id, PSite), colwise(mean), .cols = value_cols)
summary_df <- df %>%
  group_by(uniprot_id, PSite) %>%
  summarise(across(value_cols, mean), ID = max(ID))

summary_df = as.data.frame(summary_df)
rownames(summary_df) <- summary_df$ID

norm_intensity_filter <- as.matrix(select(summary_df, -c(ID, PSite, uniprot_id)))

write.table(norm_intensity_filter, "../data/processed_data/norm_intensity_filter.txt", sep = "\t")


##### collapse
df <- as.data.frame(norm_intensity)

names_input = c("norm_intensity")

df$uniprot_id <- sapply(strsplit(rownames(df), ";"), "[[", 1)
df$symbol <- sapply(strsplit(rownames(df), ";"), "[[", 2)
df$PSite <- sapply(strsplit(rownames(df), ";"), "[[", 3)
df$ID <- rownames(df)

value_cols <- colnames(df)[startsWith(colnames(df), "x")]

#summary_df <- ddply(df, .(uniprot_id, PSite), colwise(mean), .cols = value_cols)
summary_df <- df %>%
  group_by(uniprot_id, symbol) %>%
  summarise(across(value_cols, mean), ID = max(ID))

summary_df = as.data.frame(summary_df)
rownames(summary_df) <- summary_df$ID

norm_intensity_collapse <- as.matrix(select(summary_df, -c(ID, symbol, uniprot_id)))
