# Marker expression figures and tables for the POC dataset

# Loading -----------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(writexl)

# Paths -------------------------------------------------------------------------------
input_file <- 'output/01_metadata_clustering_qc/poc_hc_annotated.rds'
poc_hc <- readRDS(input_file)

output_dir <- 'output/02_population_characterization/poc/'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Excitatory/Inhibitory marker expression ---------------------------------------------
legend_colors <- c('#D1D1D1', '#2D00FF')

# Create the DimPlot for Slc17a7
FeaturePlot(poc_hc, features = c('Slc17a7'), 
            cols = legend_colors, min.cutoff = 'q1') +
  theme_void() +   # Removes the background grid
  theme(axis.title = element_blank(),  # Removes axis titles
        legend.position = 'right',      # Positions the legend to the right
        #legend.title = element_blank(), # Optional: Remove legend title
        plot.title = element_blank())

Slc17a7_umap <- FeaturePlot(poc_hc, features = c('Slc17a7'), 
                            cols = legend_colors, min.cutoff = 'q1') +
  theme_void() +   # Removes the background grid
  theme(axis.title = element_blank(),  # Removes axis titles
        legend.position = 'right',      # Positions the legend to the right
        #legend.title = element_blank(), # Optional: Remove legend title
        plot.title = element_blank())

ggsave(file.path(output_dir, 'Slc17a7.svg'), Slc17a7_umap, width = 10, height = 10)

# Create the DimPlot for Gad1
FeaturePlot(poc_hc, features = c('Gad1'), 
            cols = legend_colors, min.cutoff = 'q1') +
  theme_void() +   # Removes the background grid
  theme(axis.title = element_blank(),  # Removes axis titles
        legend.position = 'right',      # Positions the legend to the right
        #legend.title = element_blank(), # Optional: Remove legend title
        plot.title = element_blank())

Gad1_umap <- FeaturePlot(poc_hc, features = c('Gad1'), 
                         cols = legend_colors, min.cutoff = 'q1') +
  theme_void() +   # Removes the background grid
  theme(axis.title = element_blank(),  # Removes axis titles
        legend.position = 'right',      # Positions the legend to the right
        #legend.title = element_blank(), # Optional: Remove legend title
        plot.title = element_blank())

ggsave(file.path(output_dir, 'Gad1.svg'), Gad1_umap, width = 10, height = 10)

# Marker gene expression --------------------------------------------------------------
marker_genes <- c('Rfx3', 'Cux2', # 'ITL23'
                  'Rorb', 'Slc7a11', # 'ITL5'
                  'Col6a1', 'Col6a2', # 'ITL6'
                  'Ndst4', 'Nrp2',# 'ITvm'
                  'Syt6', 'Foxp2', # 'CTL6'
                  'Ctgf', 'Cplx3', # 'CTL6b'
                  'Gpc5', 'Fezf2', # 'ETL5'
                  'Tshz2', 'Htr4', # 'NPL5'
                  'F2r', 'Kcnc2', # 'Pvalb'
                  'Sst', 'Elfn1', # 'Sst'
                  'Slc6a1', 'Unc5b', # 'PvalbChand'
                  'Chodl', 'Nos1', # 'SstChodl'
                  'Vip', 'Prox1', # 'Vip'
                  'Lamp5', 'Egfr', # 'Lamp5'
                  'Htr3a', 'Frem1' # 'Sncg'
)

cluster_order <- c(
            'ITL23', 
            'ITL5', 
            'ITL6', 
            'ITvm', 
            'CTL6', 
            'CTL6b', 
            'ETL5', 
            'NPL5', 
            'Pvalb', 
            'Sst', 
            'PvalbChand', 
            'SstChodl', 
            'Vip', 
            'Lamp5', 
            'Sncg'
)

# Compute average expression per cluster
marker_genes <- marker_genes
avg_expr <- AverageExpression(poc_hc, features = marker_genes, group.by = 'cluster_name')$RNA
expr_matrix <- as.matrix(avg_expr)
expr_matrix <- expr_matrix[marker_genes, cluster_order, drop = FALSE]

# Create table to show z-scored expression values
scaled_expr_matrix <- t(apply(expr_matrix, 1, scale))
colnames(scaled_expr_matrix) <- colnames(expr_matrix)

write.csv(scaled_expr_matrix, file.path(output_dir, 'poc_scaled_expr_matrix.csv'), row.names = TRUE) 

# General marker table ----------------------------------------------------------------
seurat_obj <- poc_hc
print(levels(Idents(seurat_obj)))

# Find all marker genes
markers <- FindAllMarkers(
  object = seurat_obj,
  only.pos = FALSE,        # keep both positive and negative log2FC genes
  test.use = 'wilcox',
  logfc.threshold = 0,     # keep all tested genes
  min.pct = 0
)

# Choose only the top 20 hits
top20_markers <- markers %>%
  filter(
    !grepl('gm|rik', gene, ignore.case = TRUE),   # filter out pseudogenes
    avg_log2FC > 0
  ) %>%
  mutate(
    specificity = pct.1 - pct.2,
    marker_score = avg_log2FC * specificity
  ) %>%
  group_by(cluster) %>%
  slice_max(
    order_by = marker_score,
    n = 20,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(cluster, desc(marker_score))

# Format table for Excel
marker_table <- top20_markers %>%
  select(
    p_val,
    avg_log2FC,
    pct.1,
    pct.2,
    p_val_adj,
    cluster,
    gene
  ) %>%
  mutate(
    gene = paste0("'", gene),
    cluster = as.character(cluster)
  ) %>%
  arrange(
    factor(cluster, levels = levels(Idents(seurat_obj))),
    p_val_adj,
    desc(abs(avg_log2FC))
  )

# Save table
write_xlsx(marker_table, file.path(output_dir, 'Seurat_marker_table.xlsx'))

