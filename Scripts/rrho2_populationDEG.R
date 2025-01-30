
library(devtools)
library(RRHO2)

#Load in csv files of population deg comparisons per cluster

ITL23_A_HC <- read.csv("/Users/holmesar/Library/CloudStorage/Box-Box/RM_Projects/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/NeuroResource_January2025/DataForFigures/F5/population_degs/ITL23_group_Active_Homecage.csv")
ITL23_A_NA <- read.csv("/Users/holmesar/Library/CloudStorage/Box-Box/RM_Projects/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/NeuroResource_January2025/DataForFigures/F5/population_degs/ITL23_group_Active_Nonactive.csv")
ITL23_NA_HC <- read.csv("/Users/holmesar/Library/CloudStorage/Box-Box/RM_Projects/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/NeuroResource_January2025/DataForFigures/F5/population_degs/ITL23_group_Nonactive_Homecage.csv")


A_HC <- ITL23_A_HC
A_NA <- ITL23_A_NA
NA_HC <- ITL23_NA_HC

# Extract columns with the specified suffix

A_HC <- A_HC[, c("gene", "log2FoldChange", "padj")]  # Select only the specific columns you need from df1
A_NA <- A_NA[, c("gene", "log2FoldChange", "padj")] 
NA_HC <- NA_HC[, c("gene", "log2FoldChange", "padj")] 

# Merge columns from both dataframes into a single dataframe, filling missing values with zeros
merged_df <- merge(A_NA, NA_HC, by = "gene", suffixes = c("_A_NA", "_NA_HC"))
  
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
merged_df[[paste0("DDE_", "NA_HC")]] <- NA
merged_df[[paste0("DDE_", "NA_HC")]] <- ifelse(
  merged_df[["log2FoldChange_NA_HC"]] > 0, 
  -log10(merged_df[["padj_NA_HC"]]), 
  ifelse(
    merged_df[["log2FoldChange_NA_HC"]] == 0 | is.na(merged_df[["log2FoldChange_NA_HC"]]),
    0,
    log10(merged_df[["padj_NA_HC"]])
  )
)


# Return the modified data frame
return(merged_df)

gene_names <- merged_df[["gene"]]

gene_list1 <- merged_df[, c("gene", "DDE_A_NA")]
gene_list2 <- merged_df[, c("gene", "DDE_NA_HC")]

# Now initialize the RRHO2 object
RRHO_obj <- RRHO2_initialize(gene_list1, gene_list2, labels = c("A_NA", "NA_HC"), log10.ind = TRUE, stepsize = 10)

# RRHO_obj <- RRHO2_initialize(
#   gene_list1,
#   gene_list2,
#   stepsize = defaultStepSize(list1, list2),
#   labels = NULL,
#   log10.ind = FALSE,
#   multipleTesting = "none",
#   boundary = 0.1,
#   method = "hyper"
# )
pdf("bin10_ANA_NAHC_rrho.pdf", height = 20, width = 20)
RRHO2_heatmap(RRHO_obj, colorGradient = )
dev.off()
