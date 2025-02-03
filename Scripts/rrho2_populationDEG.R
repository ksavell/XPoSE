
library(colorRamp2)
library(circlize)
library(ComplexHeatmap)
library(RRHO2)

#Load in csv files of population deg comparisons per cluster

ITL23_A_HC <- read.csv("~/Library/CloudStorage/Box-Box/RM_Projects/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/NeuroResource_January2025/DataForFigures/F5/population_degs/ITL23_group_Active_Homecage.csv")
ITL23_A_NA <- read.csv("~/Library/CloudStorage/Box-Box/RM_Projects/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/NeuroResource_January2025/DataForFigures/F5/population_degs/ITL23_group_Active_Nonactive.csv")
ITL23_NA_HC <- read.csv("~/Library/CloudStorage/Box-Box/RM_Projects/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/NeuroResource_January2025/DataForFigures/F5/population_degs/ITL23_group_Nonactive_Homecage.csv")


A_HC <- ITL23_A_HC
A_NA <- ITL23_A_NA
NA_HC <- ITL23_NA_HC

# Extract columns with the specified suffix

A_HC <- A_HC[, c("gene", "log2FoldChange", "padj")]  # Select only the specific columns you need from df1
A_NA <- A_NA[, c("gene", "log2FoldChange", "padj")] 
NA_HC <- NA_HC[, c("gene", "log2FoldChange", "padj")] 

# Merge columns from both dataframes into a single dataframe, filling missing values with zeros
merged_df <- merge(A_NA, A_HC, by = "gene", suffixes = c("_A_NA", "_A_HC"))
  
merged_df <- merged_df[complete.cases(merged_df), ]
  
# DDE for df1
merged_df[[paste0("DDE_", "A_NA")]] <- NA
merged_df[[paste0("DDE_", "A_NA")]] <- ifelse(
  merged_df[["log2FoldChange_A_NA"]] > 0, 
  -log10(merged_df[["padj_A_NA"]]), 
  ifelse(
    merged_df[["log2FoldChange_A_NA"]] == 0 | is.na(merged_df[["log2FoldChange_A_NA"]]),
    0,
    log10(merged_df[["padj_A_NA"]])
  )
)

# DDE for df2
merged_df[[paste0("DDE_", "A_HC")]] <- NA
merged_df[[paste0("DDE_", "A_HC")]] <- ifelse(
  merged_df[["log2FoldChange_A_HC"]] > 0, 
  -log10(merged_df[["padj_A_HC"]]), 
  ifelse(
    merged_df[["log2FoldChange_A_HC"]] == 0 | is.na(merged_df[["log2FoldChange_A_HC"]]),
    0,
    log10(merged_df[["padj_A_HC"]])
  )
)


# Return the modified data frame
return(merged_df)

gene_names <- merged_df[["gene"]]

gene_list1 <- merged_df[, c("gene", "DDE_A_NA")]
gene_list2 <- merged_df[, c("gene", "DDE_A_HC")]

# Now initialize the RRHO2 object
RRHO_obj <- RRHO2_initialize(gene_list1, gene_list2, labels = c("A_NA", "A_HC"), log10.ind = TRUE, stepsize = 100)


rrho_matrix <- as.matrix(RRHO_obj[["hypermat"]])


custom_colors <- colorRamp2(
  c(0, 1.3, 1.305, 500, 1000, 2300),  # Corresponding -log10(pval) values (e.g., 0.05, 0.01, etc.)
  c("gray90", "gray90", "paleturquoise", "dodgerblue1","royalblue1", "royalblue4")  # Colors mapped to these values
)

rrho_matrix[is.na(rrho_matrix)] <- 0


RRHO2_heatmap(RRHO_obj)

Heatmap(rrho_matrix, col = custom_colors, name = "-log10(pval)", cluster_rows = FALSE, cluster_columns = FALSE)

heatmap <- Heatmap(rrho_matrix, col = custom_colors, name = "-log10(pval)", cluster_rows = FALSE, cluster_columns = FALSE)

pdf("heatmap_NAA_AHC.pdf", height = 20, width = 20)


pdf("heatmap_NAA_AHC.pdf", height = 20, width = 20)
draw(heatmap, heatmap_legend_side = "right")  # Ensure it's drawn
dev.off()


RRHO2_heatmap(RRHO_obj, colorGradient = )
dev.off()
