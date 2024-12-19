# delta variance test

# this is to compare the variance in the expression of each gene in pseudobulks and pseudo-replicates

# load the data

load("glut.RData")

library(Seurat)
library(Libra)

dv <- calculate_delta_variance(glut,
                               replicate_col = "ratID",
                               cell_type_col = "cluster_name",
                               label_col = "group")

# Rank the scores
scores <- dv[["ITL23"]] %>%
  arrange(desc(DV)) %>%
  mutate(rank = row_number())

# Plot using ggplot2
ggplot(scores, aes(x = rank, y = DV, label = gene)) +
  geom_point() +
  geom_text(nudge_y = 0.2) +
  labs(title = "Ranked Scores", x = "Rank", y = "Score") +
  theme_minimal()
