
library(ComplexHeatmap)
library(tidyverse)

#Load in csv files of population deg comparisons per cluster

ITL23_A_HC <- read.csv("~/Library/CloudStorage/Box-Box/RM_Projects/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/NeuroResource_January2025/DataForFigures/F5/population_degs/Group_Active_HC_2371subset/ITL23_group_Active_Homecage.csv")
ITL23_A_NA <- read.csv("~/Library/CloudStorage/Box-Box/RM_Projects/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/NeuroResource_January2025/DataForFigures/F5/population_degs/Group_Active_Nonactive_2371subset/ITL23_group_Active_Non-active.csv")
ITL23_NA_HC <- read.csv("~/Library/CloudStorage/Box-Box/RM_Projects/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/NeuroResource_January2025/DataForFigures/F5/population_degs/Group_Nonactive_Homecage_2371subset/ITL23_group_Non-active_Homecage.csv")


A_HC <- ITL23_A_HC
A_NA <- ITL23_A_NA
NA_HC <- ITL23_NA_HC

# Extract columns with the specified suffix

A_HC <- A_HC[, c("gene", "log2FoldChange", "padj")]  # Select only the specific columns you need from df1
A_NA <- A_NA[, c("gene", "log2FoldChange", "padj")] 
NA_HC <- NA_HC[, c("gene", "log2FoldChange", "padj")] 


# A_HC vs A_NA ------------------------------------------------------------

# Merge columns from both dataframes into a single dataframe, filling missing values with zeros
merged_df <- merge(A_NA, A_HC, by = "gene", suffixes = c("_A_NA", "_A_HC"))
  
merged_df <- merged_df[complete.cases(merged_df), ]
merged_df <- merged_df %>%
  filter(padj_A_NA < 0.05 | padj_A_HC < 0.05)
  
# Create a significance category
merged_df <- merged_df %>%
  mutate(
    significance = case_when(
      padj_A_NA < 0.05 & padj_A_HC < 0.05 ~ "Both",
      padj_A_NA < 0.05 ~ "A_NA only",
      padj_A_HC < 0.05 ~ "A_HC only",
      TRUE ~ "Not Significant"
    )
  )

write.csv(merged_df, "AHC_ANA_degcomparison_ITL23.csv")

# Scatter plot with color coding
ggplot(merged_df, aes(x = log2FoldChange_A_NA, y = log2FoldChange_A_HC, color = significance)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Both" = "black", "A_NA only" = "blue", "A_HC only" = "red", "Not Significant" = "black")) +
  labs(
    x = "Log2 Fold Change (A vs NA)",
    y = "Log2 Fold Change (A vs HC)",
    title = "Comparison of Fold Changes with Significance Categories"
  ) +
  theme_classic()



# A_HC vs NA_HC -----------------------------------------------------------

# Merge columns from both dataframes into a single dataframe, filling missing values with zeros
merged_df <- merge(A_HC, NA_HC, by = "gene", suffixes = c("_A_HC", "_NA_HC"))

merged_df <- merged_df[complete.cases(merged_df), ]
merged_df <- merged_df %>%
  filter(padj_A_HC < 0.05 | padj_NA_HC < 0.05)

# Create a significance category
merged_df <- merged_df %>%
  mutate(
    significance = case_when(
      padj_A_HC < 0.05 & padj_NA_HC < 0.05 ~ "Both",
      padj_A_HC < 0.05 ~ "A_HC only",
      padj_NA_HC < 0.05 ~ "NA_HC only",
      TRUE ~ "Not Significant"
    )
  )

write.csv(merged_df, "AHC_NAHC_degcomparison_ITL23.csv")

ggplot(merged_df, aes(x = log2FoldChange_A_HC, y = log2FoldChange_NA_HC, color = significance)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Both" = "black", "A_NA only" = "blue", "A_HC only" = "red", "Not Significant" = "black")) +
  labs(
    x = "Log2 Fold Change (A vs HC)",
    y = "Log2 Fold Change (NA vs HC)",
    title = "Comparison of Fold Changes with Significance Categories"
  ) +
  theme_classic()
