# Cluster object

# Info --------------------------------------------------------------------

# This script creates
#       * hc only object
#       * hc/nc combined object
#       * pre-clustered object with mapmycell annotation

# Loading -----------------------------------------------------------------
# loads all required packages
library(Seurat)
library(tidyverse)

source("Scripts/Functions/cluster_first.R")
source("Scripts/Functions/clstr_vln.R")
source("Scripts/Functions/subset_reclust.R")
source("Scripts/Functions/join_mapmycells.R")

date_stamp <- format(Sys.Date(), "_%m%d%Y")

# load in initial 'combined' object that is output of createobject_01.R

load("QC/all_filtered_06222026.RData")

mmc_files <- list(
  c1 = "input/C1_06222026_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1782177798832.csv",
  c2 = "input/C2_06222026_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1782177805732.csv"
)

combined_all <- join_mapmycells(combined_all, mmc_files)

save(combined_all, file = paste0("all_filtered_mmc",date_stamp,".RData"))

# Calculate subclass_confident from bootstrapping probability --------------
#   prob_threshold <- 0.5
# combined_all$subclass_confident <- ifelse(
#   is.na(combined_all$subclass_bootstrapping_probability) |
#     combined_all$subclass_bootstrapping_probability < prob_threshold,
#   "low_confidence", "confident"
# )

combined_all <- subset(combined_all, celltype == "neuron")
combined_all <- subset(combined_all, nFeature_RNA > 1000)
# classes to keep
class_keep <- c("01 IT-ET Glut", "02 NP-CT-L6b Glut", "06 CTX-CGE GABA", "07 CTX-MGE GABA","08 CNU-MGE GABA")
combined_all <- subset(combined_all, subset = class_name %in% class_keep)
# combined_all <- subset(combined_all, subclass_confident == "confident")
combined_all$subclass_confident <- "confident"

# Explore homecage clusters ----------------------------------------------

combined_hc <- subset(combined_all, subset = experience == 'HC')

combined_hc <- cluster_first(combined_hc)

pdf(paste0("hc_first_clustering", date_stamp, ".pdf"), height = 10, width = 10)
DimPlot(combined_hc, label = TRUE)
dev.off()

pdf(paste0("hc_qc_check_first_clustering", date_stamp, ".pdf"), height = 50, width = 50)
VlnPlot(combined_hc, c("nFeature_RNA", "nCount_RNA",
                      "Slc17a7", "Gad1",
                      "Snap25", "Mbp", "Gja1", "Col5a3", "Gpc5",
                      "Rfx3", "Rorb",
                      "Syt6", "Ctgf",
                      "Tshz2", "Kcnc2",
                      "Sst", "Chodl",
                      "Ppp1r1b", "Meis2",
                      "Vip", "Lamp5",
                      "Ndnf", "Cck", "Unc5b",
                      "Prox1","Drd1","Drd2","Drd3"), pt.size = 0)
dev.off()

write.csv(table(combined_hc$subclass_name, combined_hc$RNA_snn_res.2),
          file = paste0("hc_first_clustering_subclass_by_cluster", date_stamp, ".csv"))

write.csv(table(combined_hc$class_name, combined_hc$RNA_snn_res.2),
          file = paste0("hc_first_clustering_class_by_cluster", date_stamp, ".csv"))


# Second clustering -------------------------------------------------------
keep_1 <- as.character(c(0:20,22:23))  # 21 is Mbp+

all <- subset_reclust(combined_hc, clust_tokeep = keep_1,
                      neigh_dim = 1:50, umap_dim = 1:40, res = 5)

pdf(paste0("hc_second_clustering", date_stamp, ".pdf"), height = 10, width = 10)
DimPlot(all, label = TRUE)
dev.off()

pdf(paste0("hc_qc_check_second_clustering", date_stamp, ".pdf"), height = 50, width = 50)
VlnPlot(all, c("nFeature_RNA", "nCount_RNA",
               "Slc17a7", "Gad1",
               "Mbp", "Gja1", "Col5a3", "Gpc5",
               "Rfx3", "Rorb",
               "Syt6", "Ctgf",
               "Tshz2", "Kcnc2",
               "Sst", "Chodl",
               "Drd1", "Drd2", "Drd3",
               "Vip", "Lamp5",
               "Egfr", "Unc5b"), pt.size = 0)
dev.off()

write.csv(table(all$subclass, all$RNA_snn_res.5),
          file = paste0("hc_second_clustering_subclass_by_cluster", date_stamp, ".csv"))

write.csv(table(all$class, all$RNA_snn_res.5),
          file = paste0("hc_second_clustering_class_by_cluster", date_stamp, ".csv"))

# Third clustering -------------------------------------------------------
# Subclass composition: confident cells only, weighted by bootstrap probability
# Used for hierarchical clustering, dot plot, and auto-annotation ----------
confident_meta <- all@meta.data[all@meta.data$subclass_confident == "confident", ]
cluster_levels  <- levels(all)
subclass_levels <- sort(unique(confident_meta$subclass_name))

# Weighted table: sum bootstrap probabilities per subclass per cluster
subclass_wtab <- sapply(cluster_levels, function(cl) {
  cells <- confident_meta[confident_meta$RNA_snn_res.5 == cl, ]
  tapply(cells$subclass_bootstrapping_probability, cells$subclass_name,
         FUN = sum, default = 0)[subclass_levels]
})
subclass_wtab[is.na(subclass_wtab)] <- 0
rownames(subclass_wtab) <- subclass_levels

# Normalise to column proportions
subclass_prop <- sweep(subclass_wtab, 2, colSums(subclass_wtab), "/")

# Hierarchical clustering of clusters (columns) and subclasses (rows)
col_hc <- hclust(dist(t(subclass_prop)))
row_hc <- hclust(dist(subclass_prop))
subclass_ordered <- subclass_prop[row_hc$order, col_hc$order]

# --- Dot plot (supplemental figure) ---
# Size  = % of cluster cells (confident only) assigned to that subclass
# Color = average bootstrapping probability for those cells
dot_df <- do.call(rbind, lapply(cluster_levels, function(cl) {
  cells <- confident_meta[confident_meta$RNA_snn_res.5 == cl, ]
  do.call(rbind, lapply(subclass_levels, function(sb) {
    sub_cells <- cells[cells$subclass_name == sb, ]
    data.frame(
      Cluster  = cl,
      Subclass = sb,
      pct      = nrow(sub_cells) / nrow(cells) * 100,
      avg_prob = if (nrow(sub_cells) > 0) mean(sub_cells$subclass_bootstrapping_probability) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}))

dot_df$Cluster  <- factor(dot_df$Cluster,  levels = colnames(subclass_ordered))
dot_df$Subclass <- factor(dot_df$Subclass, levels = rownames(subclass_ordered))

p_dot <- ggplot(dot_df[dot_df$pct > 0, ],
                aes(x = Cluster, y = Subclass, size = pct, color = avg_prob)) +
  geom_point() +
  scale_size_continuous(name = "% of cluster", range = c(0.5, 8)) +
  scale_color_viridis_c(option = "magma", name = "Avg bootstrap\nprobability",
                        limits = c(0.4, 1)) +
  labs(title = "Subclass Composition by Cluster",
       subtitle = "Weighted by probability",
       x = "Cluster (RNA_snn_res.5)", y = "MapMyCells Subclass") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8))

pdf(paste0("hc_second_clustering_subclass_dotplot", date_stamp, ".pdf"),
    width = 22, height = 14)
print(p_dot)
dev.off()

message("Subclass dot plot saved to second_clustering_subclass_dotplot", date_stamp, ".pdf")

# Auto-annotate clusters by dominant subclass (probability-weighted argmax) ----
# Uses the same weighted proportion matrix as the dot plot.
# Clusters where the top two subclasses are within 5% of each other are flagged.
ids <- apply(subclass_prop, 2, function(x) {
  sorted <- sort(x, decreasing = TRUE)
  top1   <- sorted[1]
  top2   <- sorted[2]
  if ((top1 - top2) < 0.05) {
    paste0(names(top1), "?")
  } else {
    names(top1)
  }
})

ids_df <- data.frame(
  cluster    = names(ids),
  auto_label = ids,
  top_pct    = round(apply(subclass_prop, 2, max) * 100, 1),
  flagged    = grepl("\\?", ids),
  stringsAsFactors = FALSE
)

message("Auto-annotation summary (flagged = near-tie within 5%):")
print(ids_df)

write.csv(ids_df,
          file = paste0("hc_third_clustering_auto_annotation", date_stamp, ".csv"),
          row.names = FALSE)

# >> REVIEW auto-annotation CSV and edit ids before proceeding
# >> Override any assignments by name, e.g.: ids["12"] <- "Pvalb"
ids["23"] <- "006 L4/5 IT CTX Glut"
ids["27"] <- "022 L5 ET CTX Glut"


all <- RenameIdents(all, ids)
all$cluster_name_allen <- paste(all@active.ident)

save(all, file = paste0("HC_Allen_annotated", date_stamp, ".RData"))


# Dot plot: auto-annotated clusters vs subclass ---------------------------
confident_meta2       <- all@meta.data[all@meta.data$subclass_confident == "confident", ]
cluster_levels2       <- levels(all)
subclass_levels2      <- sort(unique(confident_meta2$subclass_name))

dot_df2 <- do.call(rbind, lapply(cluster_levels2, function(cl) {
  cells <- confident_meta2[confident_meta2$cluster_name == cl, ]
  do.call(rbind, lapply(subclass_levels2, function(sb) {
    sub_cells <- cells[cells$subclass_name == sb, ]
    data.frame(
      Cluster  = cl,
      Subclass = sb,
      pct      = nrow(sub_cells) / max(nrow(cells), 1) * 100,
      avg_prob = if (nrow(sub_cells) > 0) mean(sub_cells$subclass_bootstrapping_probability) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}))

dot_df2$Cluster  <- factor(dot_df2$Cluster,  levels = cluster_levels2)
dot_df2$Subclass <- factor(dot_df2$Subclass, levels = subclass_levels2)

p_dot2 <- ggplot(dot_df2[dot_df2$pct > 0, ],
                 aes(x = Cluster, y = Subclass, size = pct, color = avg_prob)) +
  geom_point() +
  scale_size_continuous(name = "% of cluster", range = c(0.5, 8)) +
  scale_color_viridis_c(option = "magma", name = "Avg bootstrap\nprobability",
                        limits = c(0.4, 1)) +
  labs(title = "Subclass Composition by Auto-Annotated Cluster",
       subtitle = "Confident cells only (bootstrap prob ≥ 0.5); weighted by probability",
       x = "Cluster", y = "MapMyCells Subclass") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8))

pdf(paste0("HC_annotated_cluster_subclass_dotplot", date_stamp, ".pdf"),
    width = 22, height = 14)
print(p_dot2)
dev.off()

message("Auto-annotated dot plot saved to annotated_cluster_subclass_dotplot", date_stamp, ".pdf")



# Collapsed clustering for main figures -----------------------------------
# subclass_keep <- c("006 L4/5 IT CTX Glut", "030 L6 CT CTX Glut", "052 Pvalb Gaba", 
#                 "022 L5 ET CTX Glut", "056 Sst Chodl Gaba", "053 Sst Gaba", 
#                 "007 L2/3 IT CTX Glut", "022 L5 ET CTX Glut", "010 IT AON-TT-DP Glut", 
#                 "049 Lamp5 Gaba", "051 Pvalb chandelier Gaba", "032 L5 NP CTX Glut",       
#                 "004 L6 IT CTX Glut", "029 L6b CTX Glut", "046 Vip Gaba", "047 Sncg Gaba")
# all <- subset(all, subset = cluster_name_allen %in% subclass_keep)

all <- RenameIdents(all,
                    "004 L6 IT CTX Glut"           = "ITL56",
                    # "005 L5 IT CTX Glut"           = "ITL56",
                    "006 L4/5 IT CTX Glut"         = "ITL56",
                    "007 L2/3 IT CTX Glut"        = "ITL23",
                    "010 IT AON-TT-DP Glut"        = "ITvm",
                    
                    "022 L5 ET CTX Glut"           = "ETL5",
                    "029 L6b CTX Glut"             = "CTL6",
                    "030 L6 CT CTX Glut"           = "CTL6",
                    "032 L5 NP CTX Glut"           = "NPL5",
                    
                    "046 Vip Gaba"                 = "CGE",
                    "047 Sncg Gaba"                = "CGE",
                    "049 Lamp5 Gaba"               = "CGE",
                    
                    "051 Pvalb chandelier Gaba"    = "Pvalb",
                    "052 Pvalb Gaba"               = "Pvalb",
                    "053 Sst Gaba"                 = "Sst",
                    "056 Sst Chodl Gaba"           = "Sst"
                    )
                   

all$cluster_name <- paste(all@active.ident)

save(all, file = paste0("HC_annotated", date_stamp, ".RData"))

# Dot plot: collapsed clusters vs subclass --------------------------------
confident_meta3  <- all@meta.data[all@meta.data$subclass_confident == "confident", ]
cluster_levels3  <- levels(all)
subclass_levels3 <- sort(unique(confident_meta3$subclass_name))

dot_df3 <- do.call(rbind, lapply(cluster_levels3, function(cl) {
  cells <- confident_meta3[confident_meta3$cluster_name == cl, ]
  do.call(rbind, lapply(subclass_levels3, function(sb) {
    sub_cells <- cells[cells$subclass_name == sb, ]
    data.frame(
      Cluster  = cl,
      Subclass = sb,
      pct      = nrow(sub_cells) / max(nrow(cells), 1) * 100,
      avg_prob = if (nrow(sub_cells) > 0) mean(sub_cells$subclass_bootstrapping_probability) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}))

dot_df3$Cluster  <- factor(dot_df3$Cluster,  levels = cluster_levels3)
dot_df3$Subclass <- factor(dot_df3$Subclass, levels = subclass_levels3)

p_dot3 <- ggplot(dot_df3[dot_df3$pct > 0, ],
                 aes(x = Cluster, y = Subclass, size = pct, color = avg_prob)) +
  geom_point() +
  scale_size_continuous(name = "% of cluster", range = c(0.5, 8)) +
  scale_color_viridis_c(option = "magma", name = "Avg bootstrap\nprobability",
                        limits = c(0.4, 1)) +
  labs(title = "Subclass Composition by Collapsed Cluster",
       subtitle = "Confident cells only (bootstrap prob ≥ 0.5); weighted by probability",
       x = "Collapsed Cluster", y = "MapMyCells Subclass") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8))

pdf(paste0("HC_collapsed_cluster_subclass_dotplot", date_stamp, ".pdf"),
    width = 18, height = 14)
print(p_dot3)
dev.off()


# combined object ---------------------------------------------------------

combined_all <- cluster_first(combined_all)

pdf(paste0("combined_first_clustering", date_stamp, ".pdf"), height = 10, width = 10)
DimPlot(combined_all, label = TRUE)
dev.off()

pdf(paste0("combined_qc_check_first_clustering", date_stamp, ".pdf"), height = 50, width = 50)
VlnPlot(combined_all, c("nFeature_RNA", "nCount_RNA",
                       "Slc17a7", "Gad1",
                       "Snap25", "Mbp", "Gja1", "Col5a3", "Gpc5",
                       "Rfx3", "Rorb",
                       "Syt6", "Ctgf",
                       "Tshz2", "Kcnc2",
                       "Sst", "Chodl",
                       "Ppp1r1b", "Meis2",
                       "Vip", "Lamp5",
                       "Ndnf", "Cck", "Unc5b",
                       "Prox1"), pt.size = 0)
dev.off()

write.csv(table(combined_all$subclass_name, combined_all$RNA_snn_res.2),
          file = paste0("combined_first_clustering_subclass_by_cluster", date_stamp, ".csv"))

write.csv(table(combined_all$class_name, combined_all$RNA_snn_res.2),
          file = paste0("combined_first_clustering_class_by_cluster", date_stamp, ".csv"))


# Second clustering -------------------------------------------------------
keep_1 <- as.character(c(0:23,26:27))  # 24 is mural, 25 is Mpb+

all <- subset_reclust(combined_all, clust_tokeep = keep_1,
                      neigh_dim = 1:50, umap_dim = 1:40, res = 5)

pdf(paste0("combined_second_clustering", date_stamp, ".pdf"), height = 10, width = 10)
DimPlot(all, label = TRUE)
dev.off()

pdf(paste0("combined_qc_check_second_clustering", date_stamp, ".pdf"), height = 50, width = 50)
VlnPlot(all, c("nFeature_RNA", "nCount_RNA",
               "Slc17a7", "Gad1",
               "Mbp", "Gja1", "Col5a3", "Gpc5",
               "Rfx3", "Rorb",
               "Syt6", "Ctgf",
               "Tshz2", "Kcnc2",
               "Sst", "Chodl",
               "Drd1", "Drd2", "Drd3",
               "Vip", "Lamp5",
               "Egfr", "Unc5b"), pt.size = 0)
dev.off()

write.csv(table(all$subclass, all$RNA_snn_res.5),
          file = paste0("combined_second_clustering_subclass_by_cluster", date_stamp, ".csv"))

write.csv(table(all$class, all$RNA_snn_res.5),
          file = paste0("combined_second_clustering_class_by_cluster", date_stamp, ".csv"))

# Subclass composition: confident cells only, weighted by bootstrap probability
# Used for hierarchical clustering, dot plot, and auto-annotation ----------
confident_meta <- all@meta.data[all@meta.data$subclass_confident == "confident", ]
cluster_levels  <- levels(all)
subclass_levels <- sort(unique(confident_meta$subclass_name))

# Weighted table: sum bootstrap probabilities per subclass per cluster
subclass_wtab <- sapply(cluster_levels, function(cl) {
  cells <- confident_meta[confident_meta$RNA_snn_res.5 == cl, ]
  tapply(cells$subclass_bootstrapping_probability, cells$subclass_name,
         FUN = sum, default = 0)[subclass_levels]
})
subclass_wtab[is.na(subclass_wtab)] <- 0
rownames(subclass_wtab) <- subclass_levels

# Normalise to column proportions
subclass_prop <- sweep(subclass_wtab, 2, colSums(subclass_wtab), "/")

# Hierarchical clustering of clusters (columns) and subclasses (rows)
col_hc <- hclust(dist(t(subclass_prop)))
row_hc <- hclust(dist(subclass_prop))
subclass_ordered <- subclass_prop[row_hc$order, col_hc$order]

# --- Dot plot (supplemental figure) ---
# Size  = % of cluster cells (confident only) assigned to that subclass
# Color = average bootstrapping probability for those cells
dot_df <- do.call(rbind, lapply(cluster_levels, function(cl) {
  cells <- confident_meta[confident_meta$RNA_snn_res.5 == cl, ]
  do.call(rbind, lapply(subclass_levels, function(sb) {
    sub_cells <- cells[cells$subclass_name == sb, ]
    data.frame(
      Cluster  = cl,
      Subclass = sb,
      pct      = nrow(sub_cells) / nrow(cells) * 100,
      avg_prob = if (nrow(sub_cells) > 0) mean(sub_cells$subclass_bootstrapping_probability) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}))

dot_df$Cluster  <- factor(dot_df$Cluster,  levels = colnames(subclass_ordered))
dot_df$Subclass <- factor(dot_df$Subclass, levels = rownames(subclass_ordered))

p_dot <- ggplot(dot_df[dot_df$pct > 0, ],
                aes(x = Cluster, y = Subclass, size = pct, color = avg_prob)) +
  geom_point() +
  scale_size_continuous(name = "% of cluster", range = c(0.5, 8)) +
  scale_color_viridis_c(option = "magma", name = "Avg bootstrap\nprobability",
                        limits = c(0.4, 1)) +
  labs(title = "Subclass Composition by Cluster",
       subtitle = "Confident cells only (bootstrap prob ≥ 0.5); weighted by probability",
       x = "Cluster (RNA_snn_res.5)", y = "MapMyCells Subclass") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8))

pdf(paste0("combined_second_clustering_subclass_dotplot", date_stamp, ".pdf"),
    width = 22, height = 14)
print(p_dot)
dev.off()

# Auto-annotate clusters by dominant subclass (probability-weighted argmax) ----
# Uses the same weighted proportion matrix as the dot plot.
# Clusters where the top two subclasses are within 5% of each other are flagged.
ids <- apply(subclass_prop, 2, function(x) {
  sorted <- sort(x, decreasing = TRUE)
  top1   <- sorted[1]
  top2   <- sorted[2]
  if ((top1 - top2) < 0.05) {
    paste0(names(top1), "?")
  } else {
    names(top1)
  }
})

ids_df <- data.frame(
  cluster    = names(ids),
  auto_label = ids,
  top_pct    = round(apply(subclass_prop, 2, max) * 100, 1),
  flagged    = grepl("\\?", ids),
  stringsAsFactors = FALSE
)

message("Auto-annotation summary (flagged = near-tie within 5%):")
print(ids_df)

write.csv(ids_df,
          file = paste0("combined_third_clustering_auto_annotation", date_stamp, ".csv"),
          row.names = FALSE)

# >> REVIEW auto-annotation CSV and edit ids before proceeding
# >> Override any assignments by name, e.g.: ids["12"] <- "Pvalb"

all <- RenameIdents(all, ids)
all$cluster_name_allen <- paste(all@active.ident)

save(all, file = paste0("combined_Allen_annotated", date_stamp, ".RData"))


# Dot plot: auto-annotated clusters vs subclass ---------------------------
confident_meta2       <- all@meta.data[all@meta.data$subclass_confident == "confident", ]
cluster_levels2       <- levels(all)
subclass_levels2      <- sort(unique(confident_meta2$subclass_name))

dot_df2 <- do.call(rbind, lapply(cluster_levels2, function(cl) {
  cells <- confident_meta2[confident_meta2$cluster_name == cl, ]
  do.call(rbind, lapply(subclass_levels2, function(sb) {
    sub_cells <- cells[cells$subclass_name == sb, ]
    data.frame(
      Cluster  = cl,
      Subclass = sb,
      pct      = nrow(sub_cells) / max(nrow(cells), 1) * 100,
      avg_prob = if (nrow(sub_cells) > 0) mean(sub_cells$subclass_bootstrapping_probability) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}))

dot_df2$Cluster  <- factor(dot_df2$Cluster,  levels = cluster_levels2)
dot_df2$Subclass <- factor(dot_df2$Subclass, levels = subclass_levels2)

p_dot2 <- ggplot(dot_df2[dot_df2$pct > 0, ],
                 aes(x = Cluster, y = Subclass, size = pct, color = avg_prob)) +
  geom_point() +
  scale_size_continuous(name = "% of cluster", range = c(0.5, 8)) +
  scale_color_viridis_c(option = "magma", name = "Avg bootstrap\nprobability",
                        limits = c(0.4, 1)) +
  labs(title = "Subclass Composition by Auto-Annotated Cluster",
       subtitle = "Confident cells only (bootstrap prob ≥ 0.5); weighted by probability",
       x = "Cluster", y = "MapMyCells Subclass") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8))

pdf(paste0("combined_annotated_cluster_subclass_dotplot", date_stamp, ".pdf"),
    width = 22, height = 14)
print(p_dot2)
dev.off()

# Collapsed clustering for main figures -----------------------------------

all <- RenameIdents(all,
                    "004 L6 IT CTX Glut"           = "ITL56",
                    # "005 L5 IT CTX Glut"           = "ITL56",
                    "006 L4/5 IT CTX Glut"         = "ITL56",
                    "007 L2/3 IT CTX Glut"        = "ITL23",
                    "010 IT AON-TT-DP Glut"        = "ITvm",
                    
                    "022 L5 ET CTX Glut"           = "ETL5",
                    "029 L6b CTX Glut"             = "CTL6",
                    "030 L6 CT CTX Glut"           = "CTL6",
                    "032 L5 NP CTX Glut"           = "NPL5",
                    
                    "046 Vip Gaba"                 = "CGE",
                    "047 Sncg Gaba"                = "CGE",
                    "049 Lamp5 Gaba"               = "CGE",
                    
                    "051 Pvalb chandelier Gaba"    = "Pvalb",
                    "052 Pvalb Gaba"               = "Pvalb",
                    "053 Sst Gaba"                 = "Sst",
                    "056 Sst Chodl Gaba"           = "Sst"
)

all$cluster_name <- paste(all@active.ident)

save(all, file = paste0("combined_annotated", date_stamp, ".RData"))

# Dot plot: collapsed clusters vs subclass --------------------------------
confident_meta3  <- all@meta.data[all@meta.data$subclass_confident == "confident", ]
cluster_levels3  <- levels(all)
subclass_levels3 <- sort(unique(confident_meta3$subclass_name))

dot_df3 <- do.call(rbind, lapply(cluster_levels3, function(cl) {
  cells <- confident_meta3[confident_meta3$cluster_name == cl, ]
  do.call(rbind, lapply(subclass_levels3, function(sb) {
    sub_cells <- cells[cells$subclass_name == sb, ]
    data.frame(
      Cluster  = cl,
      Subclass = sb,
      pct      = nrow(sub_cells) / max(nrow(cells), 1) * 100,
      avg_prob = if (nrow(sub_cells) > 0) mean(sub_cells$subclass_bootstrapping_probability) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}))

dot_df3$Cluster  <- factor(dot_df3$Cluster,  levels = cluster_levels3)
dot_df3$Subclass <- factor(dot_df3$Subclass, levels = subclass_levels3)

p_dot3 <- ggplot(dot_df3[dot_df3$pct > 0, ],
                 aes(x = Cluster, y = Subclass, size = pct, color = avg_prob)) +
  geom_point() +
  scale_size_continuous(name = "% of cluster", range = c(0.5, 8)) +
  scale_color_viridis_c(option = "magma", name = "Avg bootstrap\nprobability",
                        limits = c(0.4, 1)) +
  labs(title = "Subclass Composition by Collapsed Cluster",
       subtitle = "Confident cells only (bootstrap prob ≥ 0.5); weighted by probability",
       x = "Collapsed Cluster", y = "MapMyCells Subclass") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8))

pdf(paste0("combined_collapsed_cluster_subclass_dotplot", date_stamp, ".pdf"),
    width = 18, height = 14)
print(p_dot3)
dev.off()

# Generate marker tables for each object ----------------------------------

load(paste0("HC_annotated", date_stamp,".RData")) # join layers error
hc_markers <- FindAllMarkers(all)
hc_top10 <- hc_markers %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 10, with_ties = FALSE) %>%
  arrange(cluster, p_val_adj)
hc_top10$gene <- paste0("'", hc_top10$gene)
write.csv(hc_top10, "hc_top10_markers.csv", row.names = FALSE)

load("combined_annotated", date_stamp, ".RData")
all_markers <- FindAllMarkers(all)
all_top10 <- all_markers %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 10, with_ties = FALSE) %>%
  arrange(cluster, p_val_adj)
all_top10$gene <- paste0("'", all_top10$gene) # prevents date conversion in excel
write.csv(all_top10, "all_top10_markers.csv")
