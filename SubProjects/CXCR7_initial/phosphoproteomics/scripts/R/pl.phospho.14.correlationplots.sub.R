
library(ggplot2)

data <- norm_intensity
#data <- norm_intensity
# Calculating the correlation matrix
correlation_matrix <- cor(data, use = "pairwise.complete.obs")
# Transform the correlation matrix into a long format
long_corr <- reshape2::melt(correlation_matrix)
# Plot with correlation values
tiff("../analysis/PCA/correlation_normalized.tiff", 
     width = 10, height = 8, units = 'in', res = 600, compression = "lzw")
ggplot(long_corr, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), size = 2, color = "black") +
  scale_fill_gradient2(low = "white", high = "red4", mid = "mistyrose", 
                       midpoint = median(long_corr$value), limit = c(min(long_corr$value), max(long_corr$value)), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_blank()) +
  labs(title = "Correlation Matrix of Samples")
dev.off()
