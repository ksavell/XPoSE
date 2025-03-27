# Figure S4
library(Seurat)
library(tidyverse)

# load the downsampled object
load("all_downsampledforde_03212025.RData") # loaded ident is group

# Panel A: a: summary DEGs for Wilcoxon (FindMarkers) ---------------------


clusters <- unique(all$cluster_name)

for (cl in clusters) {
  
  fm <- FindMarkers(all,
                    ident.1 =  "Active",
                    ident.2 = "Homecage")
  fm$score <- ifelse(fm$p_val_adj < 0.05, 1, 0)
  fm$score_updn <- ifelse(fm$p_val_adj < 0.05 & fm$avg_log2FC > 0, 1,  # Upregulated
                                  ifelse(fm$p_val_adj < 0.05 & fm$avg_log2FC < 0, 2,  # Downregulated
                                         0))
  fm_name <- paste0("Findmarkers_cluster_",cl,".csv")
  write.csv(fm, file = fm_name)
}


# ITL23 correspondence plot -----------------------------------------------

fm_itl23 <- read.csv("FM_downsamplebyratgroup/Findmarkers_cluster_ITL23.csv",row.names = 1)

deseq_itl23 <- read.csv("DESeq_downsamplebyratgroup/ITL23_group_Active_Homecage.csv",
                        row.names = 1)

common_rows <- intersect(rownames(fm_itl23), rownames(deseq_itl23))

# Subset both to only those common row names
fm_sub <- fm_itl23[common_rows, , drop = FALSE]
deseq_sub <- deseq_itl23[common_rows, , drop = FALSE]
deseq_sub <- deseq_sub[match(rownames(fm_sub), rownames(deseq_sub)), , drop = FALSE]


# Combine them side-by-side
merged_df <- cbind(fm_sub, deseq_sub)

merged_df <- merged_df %>%
  mutate(sig_category = case_when(
    p_val_adj < 0.05 & padj < 0.05 ~ "both",
    p_val_adj < 0.05 & (is.na(padj) | padj >= 0.05) ~ "fm_only",
    padj < 0.05 & (is.na(p_val_adj) | p_val_adj >= 0.05) ~ "deseq_only",
    TRUE ~ "not_significant"
  ))

write.csv(merged_df, "figureS4_itl23_fm_deseq_correspondence.csv")
