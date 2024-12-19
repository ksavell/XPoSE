# Function to calculate Jaccard index
jaccard_index <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  return(intersection / union)
}

# Calculate Jaccard indices between all pairs of iterations
num_iterations <- length(all_indices)
jaccard_matrix <- matrix(NA, nrow = num_iterations, ncol = num_iterations)

for (i in 1:num_iterations) {
  for (j in i:num_iterations) {
    jaccard_matrix[i, j] <- jaccard_index(all_indices[[i]], all_indices[[j]])
    jaccard_matrix[j, i] <- jaccard_matrix[i, j]
  }
}

# Convert the matrix to a data frame for visualization
jaccard_df <- melt(jaccard_matrix)
colnames(jaccard_df) <- c("Iteration1", "Iteration2", "JaccardIndex")

# Plot a heatmap of the Jaccard indices
ggplot(jaccard_df, aes(x = Iteration1, y = Iteration2, fill = JaccardIndex)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(0, 1), space = "Lab", 
                       name="Jaccard Index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Jaccard Index Heatmap", x = "Iteration", y = "Iteration")