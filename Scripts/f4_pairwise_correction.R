# figure 4 correcting for multiple comparisons 4I-J

# import the pvalues from pairwise comparisons

f4 <- read.csv("~/Library/CloudStorage/Box-Box/RM_Projects/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/NeuroResource_November2024/DataForFigures/F4/F4_pairwise.csv", row.names = 1)

# Combine all p-values into a single vector for correction
all_p_values <- as.vector(as.matrix(f4))

# Apply multiple comparisons correction (e.g., Benjamini-Hochberg)
adjusted_p_values <- p.adjust(all_p_values, method = "BH")

# Reshape adjusted p-values back into the original data frame structure
adjusted_df <- data.frame(
  HCvsNA = adjusted_p_values[1:nrow(f4)],
  HCvsA = adjusted_p_values[(nrow(f4) + 1):(2 * nrow(f4))],
  AvsNA = adjusted_p_values[(2 * nrow(f4) + 1):(3 * nrow(f4))]
)
rownames(adjusted_df) <- rownames(f4)

# write to new csv
write.csv(adjusted_df, "F4_pairwise_BHcorrection.csv")
