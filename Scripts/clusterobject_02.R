# Cluster HC-only and HC/NC-combined objects
# Updated 07/2026:  
#           mapmycells join (2 cartridges) pulling subclass + supertype;
#           neuron/other assignment; iterative clustering with manual review;
#           all-cells auto-annotation weighted by bootstrap probability;
#           subclass + supertype ("expand") collapse for main figures

# Info --------------------------------------------------------------------
# This script creates:
#   * hc only object
#   * hc/nc combined object
#   * pre-clustered object with mapmycell annotation
#
# For each object it:
#   1. Loads cleaned filtered object and joins MapMyCells annotations inline
#   2. Assigns neuron/other based on class annotation
#   3. Subsets to neurons, nFeature_RNA > 1000, and classes to keep
#   4. Iterative clustering with manual cluster review steps
#   5. All-cells auto-annotation (subclass + supertype), weighted by prob
#   6. Subclass + supertype ("expand") collapse for main figures

# Load packages -----------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(ggalluvial)
})

source("Scripts/Functions/cluster_first.R")
source("Scripts/Functions/clstr_vln.R")
source("Scripts/Functions/subset_reclust.R")

# Date stamp for outputs --------------------------------------------------
date_stamp <- format(Sys.Date(), "_%m%d%Y")

# Helper: read and join MapMyCells CSVs to a single Seurat object ---------
join_mapmycells <- function(seur_obj, mapping_files) {
  mapping <- do.call(rbind, lapply(names(mapping_files), function(cart) {
    f <- mapping_files[[cart]]
    message("Reading ", f)
    df <- read.csv(f, comment.char = "#", stringsAsFactors = FALSE)
    df$cart <- cart
    df
  }))
  
  message("Cells in object: ", ncol(seur_obj))
  message("Cells in mapping: ", nrow(mapping))
  message("Cells with no mapping: ", sum(!colnames(seur_obj) %in% mapping$cell_id))
  
  meta         <- seur_obj@meta.data
  meta$cell_id <- rownames(meta)
  meta         <- left_join(meta,
                            mapping[, c("cell_id", "class_name", "subclass_name", "subclass_bootstrapping_probability",
                                        "supertype_name", "supertype_bootstrapping_probability")],
                            by = "cell_id")
  rownames(meta)     <- meta$cell_id
  meta$cell_id       <- NULL
  seur_obj@meta.data <- meta
  
  seur_obj$class    <- seur_obj@meta.data$class_name
  seur_obj$subclass <- seur_obj@meta.data$subclass_name
  
  message("NAs in class_name after join: ", sum(is.na(seur_obj@meta.data$class_name)))
  
  # Assign neuron vs other (Allen CCN classes >= 30 are non-neuronal)
  meta              <- seur_obj@meta.data
  meta$numeric_part <- as.numeric(sub(" .*", "", meta$class))
  meta$celltype     <- ifelse(meta$numeric_part >= 30, "other", "neuron")
  meta$numeric_part <- NULL
  seur_obj@meta.data <- meta
  
  return(seur_obj)
}

# =========================================================================
# Composition + auto-annotation helpers, generalized over a "level"
# (subclass / supertype). All cells used (no confidence filter); weighted
# by the level's bootstrap probability.
# =========================================================================

# Build probability-weighted proportion matrix (rows = level, cols = cluster)
weighted_prop <- function(obj, cluster_col, level_name_col, prob_col) {
  meta         <- obj@meta.data
  cluster_lvls <- levels(obj)
  level_lvls   <- sort(unique(meta[[level_name_col]]))
  
  wtab <- sapply(cluster_lvls, function(cl) {
    cells <- meta[meta[[cluster_col]] == cl, ]
    tapply(cells[[prob_col]], cells[[level_name_col]],
           FUN = sum, default = 0)[level_lvls]
  })
  wtab[is.na(wtab)] <- 0
  rownames(wtab) <- level_lvls
  
  prop <- sweep(wtab, 2, colSums(wtab), "/")
  prop[is.na(prop)] <- 0
  prop
}

# Long-format dot-plot data (size = % of cluster, color = avg bootstrap prob)
build_dot_df <- function(obj, cluster_col, level_name_col, prob_col) {
  meta         <- obj@meta.data
  cluster_lvls <- levels(factor(meta[[cluster_col]]))
  level_lvls   <- sort(unique(meta[[level_name_col]]))
  
  do.call(rbind, lapply(cluster_lvls, function(cl) {
    cells <- meta[meta[[cluster_col]] == cl, ]
    do.call(rbind, lapply(level_lvls, function(lv) {
      lv_cells <- cells[cells[[level_name_col]] == lv, ]
      data.frame(
        Cluster  = cl,
        Level    = lv,
        pct      = nrow(lv_cells) / max(nrow(cells), 1) * 100,
        avg_prob = if (nrow(lv_cells) > 0) mean(lv_cells[[prob_col]]) else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  }))
}

# Dot plot from long df, with cluster/level ordering supplied
make_dot_plot <- function(dot_df, cluster_order, level_order,
                          level_label, title, subtitle) {
  dot_df$Cluster <- factor(dot_df$Cluster, levels = cluster_order)
  dot_df$Level   <- factor(dot_df$Level,   levels = level_order)
  
  ggplot(dot_df[dot_df$pct > 0, ],
         aes(x = Cluster, y = Level, size = pct, color = avg_prob)) +
    geom_point() +
    scale_size_continuous(name = "% of cluster", range = c(0.5, 8)) +
    scale_color_viridis_c(option = "magma", name = "Avg bootstrap\nprobability",
                          limits = c(0.4, 1)) +
    labs(title = title, subtitle = subtitle,
         x = "Cluster", y = level_label) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 8))
}

# Auto-annotate: dominant level per cluster by weighted argmax (no ? flag).
# Records near-tie info (top1 - top2 < 0.05) in the CSV only.
auto_annotate <- function(prop) {
  labels <- apply(prop, 2, function(x) names(which.max(x)))
  margin <- apply(prop, 2, function(x) {
    s <- sort(x, decreasing = TRUE); s[1] - s[2]
  })
  data.frame(
    cluster    = colnames(prop),
    auto_label = labels,
    top_pct    = round(apply(prop, 2, max) * 100, 1),
    margin_pct = round(margin * 100, 1),
    near_tie   = margin < 0.05,
    stringsAsFactors = FALSE
  )
}

# Full pipeline for one level: prop matrix -> ordered dot plot (by res.5) ->
# annotation CSV -> named ident vector. Returns everything needed downstream.
# 'prefix' tags output files per object (e.g. "hc", "combined").
process_level <- function(obj, prefix, level_tag, level_name_col, prob_col,
                          cluster_col = "RNA_snn_res.5", level_label) {
  prop <- weighted_prop(obj, cluster_col, level_name_col, prob_col)
  
  col_hc  <- hclust(dist(t(prop)))
  row_hc  <- hclust(dist(prop))
  ordered <- prop[row_hc$order, col_hc$order]
  
  dot_df <- build_dot_df(obj, cluster_col, level_name_col, prob_col)
  p_dot  <- make_dot_plot(dot_df, colnames(ordered), rownames(ordered),
                          level_label,
                          paste0(level_label, " Composition by Cluster"),
                          "All cells; weighted by bootstrap probability")
  
  pdf(paste0(prefix, "_second_clustering_", level_tag, "_dotplot", date_stamp, ".pdf"),
      width = 22, height = 14)
  print(p_dot)
  dev.off()
  message("Dot plot saved to ", prefix, "_second_clustering_", level_tag, "_dotplot", date_stamp, ".pdf")
  
  ids_df <- auto_annotate(prop)
  write.csv(ids_df,
            file = paste0(prefix, "_second_clustering_", level_tag, "_auto_annotation", date_stamp, ".csv"),
            row.names = FALSE)
  message("Auto-annotation (", prefix, " ", level_tag, ") — near_tie flags margin < 5%:")
  print(ids_df)
  
  ids <- setNames(ids_df$auto_label, ids_df$cluster)
  list(prop = prop, ids = ids, ids_df = ids_df)
}

# Subclass -> short-name collapse map (from vmPFC/dmPFC script) -----------
subclass_collapse_map <- c(
  "007 L2/3 IT CTX Glut"        = "ITL23",
  "030 L6 CT CTX Glut"          = "CTL6",
  "010 IT AON-TT-DP Glut"       = "ITvm",
  "056 Sst Chodl Gaba"          = "SstChodl",
  "052 Pvalb Gaba"              = "Pvalb",
  "006 L4/5 IT CTX Glut"        = "ITL5",
  "004 L6 IT CTX Glut"          = "ITL6",
  "022 L5 ET CTX Glut"          = "ETL5",
  "053 Sst Gaba"                = "Sst",
  "029 L6b CTX Glut"            = "CTL6b",
  "032 L5 NP CTX Glut"          = "NPL5",
  "005 L5 IT CTX Glut"          = "ITL5",
  "046 Vip Gaba"                = "Vip",
  "049 Lamp5 Gaba"              = "Lamp5",
  "047 Sncg Gaba"               = "Sncg",
  "051 Pvalb chandelier Gaba"   = "PvalbChand",
  "050 Lamp5 Lhx6 Gaba"         = "Lamp5"
)

# Supertype -> short-name "expand" collapse map (from vmPFC/dmPFC script) -
supertype_collapse_map <- c(
  "0013 L6 IT CTX Glut_1"          = "ITL6_1",
  "0014 L6 IT CTX Glut_2"          = "ITL6_2",
  "0017 L6 IT CTX Glut_5"          = "ITL6_5",
  "0018 L5 IT CTX Glut_1"          = "ITL5_1",
  "0022 L5 IT CTX Glut_5"          = "ITL5_5",
  "0023 L4/5 IT CTX Glut_1"        = "ITL45_1",
  "0024 L4/5 IT CTX Glut_2"        = "ITL45_2",
  "0025 L4/5 IT CTX Glut_3"        = "ITL45_3",
  "0026 L4/5 IT CTX Glut_4"        = "ITL45_4",
  "0027 L4/5 IT CTX Glut_5"        = "ITL45_5",
  "0029 L2/3 IT CTX Glut_1"        = "ITL23_1",
  "0030 L2/3 IT CTX Glut_2"        = "ITL23_2",
  "0032 L2/3 IT CTX Glut_4"        = "ITL23_4",
  "0046 IT AON-TT-DP Glut_1"       = "ITLvm_1",
  "0047 IT AON-TT-DP Glut_2"       = "ITLvm_2",
  "0090 L5 ET CTX Glut_1"          = "ETL5_1",
  "0091 L5 ET CTX Glut_2"          = "ETL5_2",
  "0094 L5 ET CTX Glut_5"          = "ETL5_5",
  "0110 L6b CTX Glut_1"            = "CTL6b_1",
  "0114 L6 CT CTX Glut_1"          = "CTL6_1",
  "0115 L6 CT CTX Glut_2"          = "CTL6_2",
  "0116 L6 CT CTX Glut_3"          = "CTL6_3",
  "0118 L6 CT CTX Glut_5"          = "CTL6_5",
  "0123 L5 NP CTX Glut_2"          = "NPL5_2",
  "0173 Vip Gaba_1"                = "Vip_1",
  "0177 Vip Gaba_5"                = "Vip_5",
  "0187 Sncg Gaba_3"               = "Sncg_3",
  "0199 Lamp5 Gaba_1"              = "Lamp5_1",
  "0203 Lamp5 Lhx6 Gaba_1"         = "Lamp5_Lhx_1",
  "0204 Pvalb chandelier Gaba_1"   = "PvalbChand_1",
  "0205 Pvalb Gaba_1"              = "Pvalb_1",
  "0207 Pvalb Gaba_3"              = "Pvalb_3",
  "0208 Pvalb Gaba_4"              = "Pvalb_4",
  "0212 Pvalb Gaba_8"              = "Pvalb_8",
  "0214 Sst Gaba_1"                = "Sst_1",
  "0217 Sst Gaba_4"                = "Sst_4",
  "0218 Sst Gaba_5"                = "Sst_5",
  "0223 Sst Gaba_10"               = "Sst_10",
  "0225 Sst Gaba_12"               = "Sst_12",
  "0227 Sst Gaba_14"               = "Sst_14",
  "0228 Sst Gaba_15"               = "Sst_15",
  "0241 Sst Chodl Gaba_4"          = "SstChodl_4"
)

# Apply a collapse map to active idents, keeping only mappings whose source
# idents are actually present (RenameIdents errors on absent idents).
apply_collapse <- function(obj, collapse_map) {
  present <- collapse_map[names(collapse_map) %in% levels(obj)]
  if (length(present) == 0) {
    warning("No collapse-map idents present in object; returning unchanged.")
    return(obj)
  }
  do.call(RenameIdents, c(list(object = obj), as.list(present)))
}

# Annotate one clustered object: subclass + supertype auto-annotation,
# per-level dot plots, collapse to cluster_name (subclass) and
# cluster_name_expand (supertype), plus alluvials. Saves annotated RData.
annotate_object <- function(obj, prefix) {
  
  # --- Auto-annotate subclass and supertype -----------------------------
  subclass_res <- process_level(obj, prefix, "subclass", "subclass_name",
                                "subclass_bootstrapping_probability",
                                level_label = "MapMyCells Subclass")
  supertype_res <- process_level(obj, prefix, "supertype", "supertype_name",
                                 "supertype_bootstrapping_probability",
                                 level_label = "MapMyCells Supertype")
  
  # >> REVIEW auto-annotation CSVs before proceeding (no manual overrides)
  
  obj$cluster_name_subclass  <- unname(subclass_res$ids[as.character(obj$RNA_snn_res.5)])
  obj$cluster_name_supertype <- unname(supertype_res$ids[as.character(obj$RNA_snn_res.5)])
  
  save(obj, file = paste0(prefix, "_Allen_annotated", date_stamp, ".RData"))
  
  # --- Dot plots: auto-annotated clusters vs level ----------------------
  annot_subclass_df <- build_dot_df(obj, "cluster_name_subclass",
                                    "subclass_name",
                                    "subclass_bootstrapping_probability")
  p_annot_sub <- make_dot_plot(
    annot_subclass_df,
    cluster_order = sort(unique(obj$cluster_name_subclass)),
    level_order   = sort(unique(obj$subclass_name)),
    level_label   = "MapMyCells Subclass",
    title         = "Subclass Composition by Auto-Annotated Cluster",
    subtitle      = "All cells; weighted by bootstrap probability")
  
  pdf(paste0(prefix, "_annotated_cluster_subclass_dotplot", date_stamp, ".pdf"),
      width = 22, height = 14)
  print(p_annot_sub)
  dev.off()
  message("Saved ", prefix, "_annotated_cluster_subclass_dotplot", date_stamp, ".pdf")
  
  annot_supertype_df <- build_dot_df(obj, "cluster_name_supertype",
                                     "supertype_name",
                                     "supertype_bootstrapping_probability")
  p_annot_super <- make_dot_plot(
    annot_supertype_df,
    cluster_order = sort(unique(obj$cluster_name_supertype)),
    level_order   = sort(unique(obj$supertype_name)),
    level_label   = "MapMyCells Supertype",
    title         = "Supertype Composition by Auto-Annotated Cluster",
    subtitle      = "All cells; weighted by bootstrap probability")
  
  pdf(paste0(prefix, "_annotated_cluster_supertype_dotplot", date_stamp, ".pdf"),
      width = 22, height = 14)
  print(p_annot_super)
  dev.off()
  message("Saved ", prefix, "_annotated_cluster_supertype_dotplot", date_stamp, ".pdf")
  
  # --- Alluvial: subclass -> supertype ----------------------------------
  alluv_df <- as.data.frame(
    table(subclass  = obj$cluster_name_subclass,
          supertype = obj$cluster_name_supertype),
    stringsAsFactors = FALSE
  )
  alluv_df <- alluv_df[alluv_df$Freq > 0, ]
  
  p_alluv <- ggplot(alluv_df,
                    aes(axis1 = subclass, axis2 = supertype, y = Freq)) +
    geom_alluvium(aes(fill = subclass), alpha = 0.7, curve_type = "sigmoid") +
    geom_stratum(width = 0.35, fill = "grey90", color = "grey40") +
    geom_text(stat = "stratum",
              aes(label = after_stat(stratum)), size = 2.5) +
    scale_x_discrete(limits = c("Subclass", "Supertype"),
                     expand = c(0.1, 0.05)) +
    labs(title = "Cluster Identity: Subclass -> Supertype",
         subtitle = "Flow width = number of cells",
         y = "Cells") +
    theme_classic(base_size = 11) +
    theme(legend.position = "none")
  
  pdf(paste0(prefix, "_cluster_subclass_supertype_alluvial", date_stamp, ".pdf"),
      width = 14, height = 18)
  print(p_alluv)
  dev.off()
  message("Saved ", prefix, "_cluster_subclass_supertype_alluvial", date_stamp, ".pdf")
  
  # =======================================================================
  # Collapsed clustering for main figures
  # =======================================================================
  # Subclass-based collapse -> cluster_name
  Idents(obj) <- "cluster_name_subclass"
  obj <- apply_collapse(obj, subclass_collapse_map)
  obj$cluster_name <- paste(obj@active.ident)
  
  # Supertype-based collapse -> cluster_name_expand
  Idents(obj) <- "cluster_name_supertype"
  obj <- apply_collapse(obj, supertype_collapse_map)
  obj$cluster_name_expand <- paste(obj@active.ident)
  
  Idents(obj) <- "cluster_name"
  save(obj, file = paste0(prefix, "_annotated", date_stamp, ".RData"))
  
  # --- Dot plot: collapsed clusters vs subclass -------------------------
  collapsed_df <- build_dot_df(obj, "cluster_name", "subclass_name",
                               "subclass_bootstrapping_probability")
  p_collapsed <- make_dot_plot(
    collapsed_df,
    cluster_order = sort(unique(obj$cluster_name)),
    level_order   = sort(unique(obj$subclass_name)),
    level_label   = "MapMyCells Subclass",
    title         = "Subclass Composition by Collapsed Cluster",
    subtitle      = "All cells; weighted by bootstrap probability")
  
  pdf(paste0(prefix, "_collapsed_cluster_subclass_dotplot", date_stamp, ".pdf"),
      width = 18, height = 14)
  print(p_collapsed)
  dev.off()
  message("Saved ", prefix, "_collapsed_cluster_subclass_dotplot", date_stamp, ".pdf")
  
  # --- Alluvial: cluster_name -> expand -> subclass -> supertype --------
  alluv3_df <- as.data.frame(
    table(cluster   = obj$cluster_name,
          expand    = obj$cluster_name_expand,
          subclass  = obj$cluster_name_subclass,
          supertype = obj$cluster_name_supertype),
    stringsAsFactors = FALSE
  )
  alluv3_df <- alluv3_df[alluv3_df$Freq > 0, ]
  
  p_alluv3 <- ggplot(alluv3_df,
                     aes(axis1 = cluster, axis2 = subclass,
                         axis3 = expand, axis4 = supertype, y = Freq)) +
    geom_alluvium(aes(fill = cluster), alpha = 0.7, curve_type = "sigmoid") +
    geom_stratum(width = 0.35, fill = "grey90", color = "grey40") +
    geom_text(stat = "stratum",
              aes(label = after_stat(stratum)), size = 2.5) +
    scale_x_discrete(limits = c("Collapsed", "Subclass", "Expand", "Supertype"),
                     expand = c(0.1, 0.05)) +
    labs(title = "Cluster Identity: Collapsed -> Subclass -> Expand -> Supertype",
         subtitle = "Flow width = number of cells",
         y = "Cells") +
    theme_classic(base_size = 11) +
    theme(legend.position = "none")
  
  pdf(paste0(prefix, "_collapsed_subclass_supertype_expand_alluvial", date_stamp, ".pdf"),
      width = 18, height = 20)
  print(p_alluv3)
  dev.off()
  message("Saved ", prefix, "_collapsed_subclass_supertype_expand_alluvial", date_stamp, ".pdf")
  
  return(obj)
}

# QC marker set (shared) --------------------------------------------------
qc_markers <- c("nFeature_RNA", "nCount_RNA",
                "Slc17a7", "Gad1",
                "Snap25", "Mbp", "Gja1", "Col5a3", "Gpc5",
                "Rfx3", "Rorb",
                "Syt6", "Ctgf",
                "Tshz2", "Kcnc2",
                "Sst", "Chodl",
                "Vip", "Lamp5",
                "Ndnf", "Cck", "Unc5b",
                "Prox1")

# QC marker set for second clustering (post-QC subset) --------------------
qc_markers_second <- c("nFeature_RNA", "nCount_RNA",
                       "Slc17a7", "Gad1",
                       "Mbp", "Gja1", "Col5a3", "Gpc5",
                       "Rfx3", "Rorb",
                       "Syt6", "Ctgf",
                       "Tshz2", "Kcnc2",
                       "Sst", "Chodl",
                       "Vip", "Lamp5",
                       "Egfr", "Unc5b")

# =========================================================================
# Load, join MapMyCells, subset
# =========================================================================
# load in initial 'combined' object that is output of createobject_01.R
load("QC/all_filtered_06222026.RData")   # -> combined_all

mmc_files <- list(
  c1 = "input/C1_06222026_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1782177798832.csv",
  c2 = "input/C2_06222026_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1782177805732.csv"
)

message("--- Joining MapMyCells ---")
combined_all <- join_mapmycells(combined_all, mmc_files)

save(combined_all, file = paste0("all_filtered_mmc", date_stamp, ".RData"))

# Subset to neurons, feature count, classes to keep -----------------------
combined_all <- subset(combined_all, celltype == "neuron")
combined_all <- subset(combined_all, nFeature_RNA > 1000)
class_keep <- c("01 IT-ET Glut", "02 NP-CT-L6b Glut", "06 CTX-CGE GABA", "07 CTX-MGE GABA", "08 CNU-MGE GABA")
combined_all <- subset(combined_all, subset = class_name %in% class_keep)

# =========================================================================
# HC-only object
# =========================================================================
combined_hc <- subset(combined_all, subset = experience == 'HC')

# First clustering --------------------------------------------------------
combined_hc <- cluster_first(combined_hc)

pdf(paste0("hc_first_clustering", date_stamp, ".pdf"), height = 10, width = 10)
DimPlot(combined_hc, label = TRUE)
dev.off()

pdf(paste0("hc_qc_check_first_clustering", date_stamp, ".pdf"), height = 50, width = 50)
VlnPlot(combined_hc, qc_markers, pt.size = 0)
dev.off()

write.csv(table(combined_hc$subclass_name, combined_hc$RNA_snn_res.2),
          file = paste0("hc_first_clustering_subclass_by_cluster", date_stamp, ".csv"))
write.csv(table(combined_hc$class_name, combined_hc$RNA_snn_res.2),
          file = paste0("hc_first_clustering_class_by_cluster", date_stamp, ".csv"))

# >> REVIEW hc_first_clustering and hc_qc_check_first_clustering PDFs
# >> Fill in clusters to keep before proceeding

# Second clustering -------------------------------------------------------
keep_1 <- as.character(c(0:20, 22:23))  # 21 is Mbp+

all <- subset_reclust(combined_hc, clust_tokeep = keep_1,
                      neigh_dim = 1:50, umap_dim = 1:40, res = 5)

pdf(paste0("hc_second_clustering", date_stamp, ".pdf"), height = 10, width = 10)
DimPlot(all, label = TRUE)
dev.off()

pdf(paste0("hc_qc_check_second_clustering", date_stamp, ".pdf"), height = 50, width = 50)
VlnPlot(all, qc_markers_second, pt.size = 0)
dev.off()

write.csv(table(all$subclass, all$RNA_snn_res.5),
          file = paste0("hc_second_clustering_subclass_by_cluster", date_stamp, ".csv"))
write.csv(table(all$class, all$RNA_snn_res.5),
          file = paste0("hc_second_clustering_class_by_cluster", date_stamp, ".csv"))

# Annotate + collapse -----------------------------------------------------
all <- annotate_object(all, prefix = "hc")

# =========================================================================
# HC/NC combined object
# =========================================================================
combined_all <- cluster_first(combined_all)

pdf(paste0("combined_first_clustering", date_stamp, ".pdf"), height = 10, width = 10)
DimPlot(combined_all, label = TRUE)
dev.off()

pdf(paste0("combined_qc_check_first_clustering", date_stamp, ".pdf"), height = 50, width = 50)
VlnPlot(combined_all, qc_markers, pt.size = 0)
dev.off()

write.csv(table(combined_all$subclass_name, combined_all$RNA_snn_res.2),
          file = paste0("combined_first_clustering_subclass_by_cluster", date_stamp, ".csv"))
write.csv(table(combined_all$class_name, combined_all$RNA_snn_res.2),
          file = paste0("combined_first_clustering_class_by_cluster", date_stamp, ".csv"))

# >> REVIEW combined_first_clustering and combined_qc_check_first_clustering PDFs
# >> Fill in clusters to keep before proceeding

# Second clustering -------------------------------------------------------
keep_1 <- as.character(c(0:23, 26:27))  # 24 is mural, 25 is Mbp+

all <- subset_reclust(combined_all, clust_tokeep = keep_1,
                      neigh_dim = 1:50, umap_dim = 1:40, res = 5)

pdf(paste0("combined_second_clustering", date_stamp, ".pdf"), height = 10, width = 10)
DimPlot(all, label = TRUE)
dev.off()

pdf(paste0("combined_qc_check_second_clustering", date_stamp, ".pdf"), height = 50, width = 50)
VlnPlot(all, qc_markers_second, pt.size = 0)
dev.off()

write.csv(table(all$subclass, all$RNA_snn_res.5),
          file = paste0("combined_second_clustering_subclass_by_cluster", date_stamp, ".csv"))
write.csv(table(all$class, all$RNA_snn_res.5),
          file = paste0("combined_second_clustering_class_by_cluster", date_stamp, ".csv"))

# Annotate + collapse -----------------------------------------------------
all <- annotate_object(all, prefix = "combined")
