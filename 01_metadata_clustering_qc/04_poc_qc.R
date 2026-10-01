# QC for POC dataset

# Loading -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(cluster)
  library(ggplot2)
  library(patchwork)
  library(RANN)
  library(vegan)
})

# Paths -------------------------------------------------------------------------------
input_file <- 'output/01_metadata_clustering_qc/poc_combined_annotated.rds'
input_file2 <- 'output/01_metadata_clustering_qc/poc_hc_annotated.rds'

poc_combined <- readRDS(input_file)
poc_hc <- readRDS(input_file2)

output_dir <- 'output/01_metadata_clustering_qc/poc_qc'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Settings ----------------------------------------------------------------------------
batch_col      <- 'capture'        # capture: C1 / C2
celltype_col   <- 'cluster_name'   # 15 fine types
ratid_col      <- 'ratID'          # animal
reduction      <- 'pca'
n_dims         <- 40               # matches UMAP/clustering embedding
no_struct_band <- 0.05             # silhouette 'no-structure' shaded zone (not a cutoff)
max_pairs      <- 20000            # sampled pairs per (cell type x comparison)
min_cells      <- 20               # skip a type if a capture has fewer cells than this
seed           <- 42
date_tag       <- format(Sys.Date(), '_%m%d%Y')
pal_group <- c('capture' = '#999999', 'cell type' = '#2d8cb8', 'ratID' = '#4d4d4d')

set.seed(seed)
emb  <- Embeddings(poc_hc, reduction)[, 1:n_dims]
meta <- poc_hc@meta.data
stopifnot(identical(rownames(emb), rownames(meta)))

# Shared distance matrix
d <- dist(emb)   # euclidean on PCA dims (~780MB for ~9.9k cells)

cluster_colors <- c(
  'ITL23'      = '#2EBF5E',
  'ITL5'       = '#50B2AD',
  'ITL6'       = '#58D2CF', 
  'ITvm'       = '#B1DE7D',
  'CTL6'       = '#2D8CB8',
  'CTL6b'      = '#7044AA',
  'ETL5'       = '#0D5A8B',
  'NPL5'       = '#3E9E64',
  'Pvalb'      = '#B9342C',
  'Sst'        = '#FF9900',
  'PvalbChand' = '#FF2D4E',
  'SstChodl'   = '#B1B10C',
  'Vip'        = '#B864CC',
  'Lamp5'      = '#DA808C',
  'Sncg'       = '#D3408D'
)

cluster_order <- c('ITL23', 
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

poc_hc$cluster_name <- factor(
  as.character(poc_hc$cluster_name),
  levels = rev(cluster_order)
)

# QC plot by cluster ------------------------------------------------------------------
# Genes expressed
p <- VlnPlot(
  poc_combined,
  features = 'nFeature_RNA',
  group.by = 'cluster_name',
  cols = cluster_colors,
  sort = FALSE,
  pt.size = 0,
  combine = FALSE)[[1]] +
  scale_x_discrete(
    labels = function(x) {sprintf("<span style='color:%s'>%s</span>",
                                  cluster_colors[x],
                                  x)}) +
  scale_y_continuous(
    breaks = seq(1000, 5000, by = 1000),
    labels = function(x) x / 1000,
    expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = NULL,
    y = 'Genes expressed (\u00D710\u00B3)') +
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 7,
                                   family = 'Arial',
                                   angle = 0),
        axis.text.y = ggtext::element_markdown(size = 7,
                                               family = 'Arial'),
        axis.title.x = element_text(size = 8,
                                    family = 'Arial'),
        axis.title.y = element_blank(),
        plot.title = element_blank(),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.5)) +
  coord_flip()

# Set violin outlines to 0.25
violin_layers <- vapply(
  p$layers,
  function(layer) inherits(layer$geom, 'GeomViolin'),
  logical(1)
)

for (i in which(violin_layers)) {p$layers[[i]]$aes_params$linewidth <- 0.25}

quartz(
  type = 'pdf',
  file = file.path(output_dir, 'FS3C_genesexpressed.pdf'),
  width = 1.5,
  height = 2.5,
  family = 'Arial'
)

print(p)
dev.off()

# Transcripts expressed
p2 <- VlnPlot(
  poc_combined,
  features = 'nCount_RNA',
  group.by = 'cluster_name',
  cols = cluster_colors,
  sort = FALSE,
  pt.size = 0,
  combine = FALSE)[[1]] +
  scale_x_discrete(labels = function(x) {
    sprintf(
      "<span style='color:%s'>%s</span>",
      cluster_colors[x], x)}) +
  scale_y_continuous(
    breaks = seq(0, 15000, by = 5000),
    labels = function(x) x / 1000,
    expand = expansion(mult = c(0, 0))) +
  labs(
    x = NULL,
    y = 'Transcripts expressed (\u00D710\u00B3)') +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(size = 7,
                               family = 'Arial',
                               angle = 0),
    axis.text.y = ggtext::element_markdown(size = 7,
                                           family = 'Arial'),
    axis.title.x = element_text(size = 8,
                                family = 'Arial'),
    axis.title.y = element_blank(),
    plot.title = element_blank(),
    
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5)) +
  coord_flip(
    ylim = c(0, 15000))

violin_layers <- vapply(
  p2$layers,
  function(layer) inherits(layer$geom, 'GeomViolin'),
  logical(1))

for (i in which(violin_layers)) {p2$layers[[i]]$aes_params$linewidth <- 0.25}

quartz(
  type = 'pdf',
  file = file.path(output_dir, 'FS3C_transcriptsexpressed.pdf'),
  width = 1.5,
  height = 2.5,
  family = 'Arial'
)

print(p2)
dev.off()

# QC plot by sort order ---------------------------------------------------------------
poc_combined$sort_day <- dplyr::case_when(
  poc_combined$orig.ident %in% c('C1', 'C2') ~ 'Sort_day_1',
  TRUE ~ NA_character_
)

# Store sort_day in the intended order
poc_combined$sort_day <- factor(
  poc_combined$sort_day,
  levels = paste0('Sort_day_', 1)
)

sortdays <- unique(poc_combined$sort_day)

# Use the same limits for every sort-day plot
gene_max <- 5000  
transcript_max <- 15000

for (day in sortdays) {
  
  # Subset the Seurat object
  poc_combined <- subset(poc_combined, subset = sort_day == day)
  
  # Create the xpose_tag-to-ratID lookup
  xpose_map <- unique(
    data.frame(
      xpose_tag = as.character(poc_combined$xpose_tag),
      ratID = as.character(poc_combined$ratID)
    )
  )
  
  # Confirm that each xpose_tag maps to only one ratID
  if (anyDuplicated(xpose_map$xpose_tag)) {
    stop('At least one xpose_tag is associated with multiple ratID values.')
  }
  
  # Confirm xpose_tag format
  valid_tags <- grepl(
    "^xpose_tag_(0[1-9]|1[0-2])_mm$",
    xpose_map$xpose_tag
  )
  
  if (any(!valid_tags)) {
    stop(
      'Unexpected xpose_tag values: ',
      paste(
        unique(xpose_map$xpose_tag[!valid_tags]),
        collapse = ', '
      )
    )
  }
  
  # Extract the number from xposeTag##_mm
  xpose_map$xpose_number <- as.integer(
    sub(
      "^xpose_tag_(0[1-9]|1[0-2])_mm$",
      '\\1',
      xpose_map$xpose_tag
    )
  )
  
  # Sort xposeTags numerically
  xpose_map <- xpose_map[
    order(xpose_map$xpose_number),
  ]
  
  xpose_order <- xpose_map$xpose_tag
  
  # Named label vector: xpose_tag -> ratID
  rat_labels <- setNames(
    xpose_map$ratID,
    xpose_map$xpose_tag
  )
  
  # Reverse factor levels so the lowest xposeTag appears at the top
  poc_combined$xpose_tag <- factor(
    as.character(poc_combined$xpose_tag),
    levels = rev(xpose_order)
  )
  
  # Same violin color for every xpose tag
  xpose_cols <- setNames(
    rep('#808080', length(xpose_order)),
    xpose_order
  )
  
  # Genes expressed
  p_genes <- VlnPlot(
    poc_combined,
    features = 'nFeature_RNA',
    group.by = 'xpose_tag',
    cols = xpose_cols,
    sort = FALSE,
    pt.size = 0,
    combine = FALSE
  )[[1]] +
    scale_x_discrete(
      labels = rat_labels
    ) +
    scale_y_continuous(
      breaks = seq(1000, 5000, by = 1000),
      labels = function(x) x / 1000,
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      x = NULL,
      y = 'Genes expressed (\u00D710\u00B3)'
    ) +
    theme(
      legend.position = 'none',
      axis.text.x = element_text(
        size = 7,
        family = 'Arial',
        angle = 0,
        hjust = 0.5
      ),
      axis.text.y = element_text(
        size = 7,
        family = 'Arial',
        color = 'black'
      ),
      axis.title.x = element_text(
        size = 8,
        family = 'Arial'
      ),
      axis.title.y = element_blank(),
      plot.title = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks = element_line(linewidth = 0.5)
    ) +
    coord_flip(
      ylim = c(1000, gene_max)
    )
  
  violin_layers <- vapply(
    p_genes$layers,
    function(layer) inherits(layer$geom, 'GeomViolin'),
    logical(1)
  )
  
  for (i in which(violin_layers)) {
    p_genes$layers[[i]]$aes_params$linewidth <- 0.25
  }
  
  quartz(
    type = 'pdf',
    file = file.path(output_dir, (paste0(day, '_xpose_genesexpressed.pdf'))),
    width = 1.5,
    height = 2.5,
    family = 'Arial'
  )
  
  print(p_genes)
  dev.off()
  
  # Transcripts expressed
  p_counts <- VlnPlot(
    poc_combined,
    features = 'nCount_RNA',
    group.by = 'xpose_tag',
    cols = xpose_cols,
    sort = FALSE,
    pt.size = 0,
    combine = FALSE
  )[[1]] +
    scale_x_discrete(
      labels = rat_labels
    ) +
    scale_y_continuous(
      breaks = seq(0, 15000, by = 5000),
      labels = function(x) x / 1000,
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x = NULL,
      y = 'Transcripts expressed (\u00D710\u00B3)'
    ) +
    theme(
      legend.position = 'none',
      axis.text.x = element_text(
        size = 7,
        family = 'Arial',
        angle = 0,
        hjust = 0.5
      ),
      axis.text.y = element_text(
        size = 7,
        family = 'Arial',
        color = 'black'
      ),
      axis.title.x = element_text(
        size = 8,
        family = 'Arial'
      ),
      axis.title.y = element_blank(),
      plot.title = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      axis.ticks = element_line(linewidth = 0.5)
    ) +
    coord_flip(
      ylim = c(0, transcript_max)
    )
  
  violin_layers <- vapply(
    p_counts$layers,
    function(layer) inherits(layer$geom, 'GeomViolin'),
    logical(1)
  )
  
  for (i in which(violin_layers)) {
    p_counts$layers[[i]]$aes_params$linewidth <- 0.25
  }
  
  quartz(
    type = 'pdf',
    file = file.path(output_dir, (paste0(day, '_xpose_transcriptsexpressed.pdf'))),
    width = 1.5,
    height = 2.5,
    family = 'Arial'
  )
  
  print(p_counts)
  dev.off()
}

# Variance within/between captures per cluster ----------------------------------------
batch <- as.character(meta[[batch_col]])
ct    <- as.character(meta[[celltype_col]])
caps  <- sort(unique(batch)); stopifnot(length(caps) == 2)

pair_dist <- function(ia, ib)
  sqrt(rowSums((emb[ia, , drop = FALSE] - emb[ib, , drop = FALSE])^2))

sample_pairs <- function(pool_a, pool_b = NULL, n, same_pool) {
  if (same_pool) {
    if (length(pool_a) < 2) return(NULL)
    ia <- sample(pool_a, n, TRUE); ib <- sample(pool_a, n, TRUE)
    keep <- ia != ib; list(a = ia[keep], b = ib[keep])
  } else {
    if (length(pool_a) < 1 || length(pool_b) < 1) return(NULL)
    list(a = sample(pool_a, n, TRUE), b = sample(pool_b, n, TRUE))
  }
}

res_list <- list(); raw_list <- list()
for (type in sort(unique(ct))) {
  in_type <- which(ct == type)
  idx1 <- in_type[batch[in_type] == caps[1]]
  idx2 <- in_type[batch[in_type] == caps[2]]
  if (length(idx1) < min_cells || length(idx2) < min_cells) next
  
  w1 <- sample_pairs(idx1, n = max_pairs %/% 2, same_pool = TRUE)
  w2 <- sample_pairs(idx2, n = max_pairs %/% 2, same_pool = TRUE)
  bt <- sample_pairs(idx1, idx2, n = max_pairs, same_pool = FALSE)
  
  d_within  <- pair_dist(c(w1$a, w2$a), c(w1$b, w2$b))
  d_between <- pair_dist(bt$a, bt$b)
  
  res_list[[type]] <- data.frame(
    cell_type = type, n_c1 = length(idx1), n_c2 = length(idx2),
    within_mean = mean(d_within), between_mean = mean(d_between),
    ratio = mean(d_between) / mean(d_within),
    pct_diff = 100 * (mean(d_between) - mean(d_within)) / mean(d_within))
  raw_list[[type]] <- rbind(
    data.frame(cell_type = type, comparison = 'within capture',  dist = d_within),
    data.frame(cell_type = type, comparison = 'between capture', dist = d_between))
}

per_type <- do.call(rbind, res_list); rownames(per_type) <- NULL
raw_df   <- do.call(rbind, raw_list)
raw_df$comparison <- factor(raw_df$comparison, levels = c('within capture', 'between capture'))

pD <- ggplot(raw_df, aes(x = dist, color = comparison)) +
  geom_density(linewidth = 0.8, adjust = 1) +
  facet_wrap(~ cell_type, scales = 'free', ncol = 7) +
  scale_color_manual(values = c('within capture' = 'gray65', 'between capture' = 'black'), 
                     labels = c('within', 'between')) +
  labs(title = 'Variance within/between captures per cluster',
       x = 'Euclidean distance (PCA)', y = 'Density', color = NULL) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text (size = 15, hjust = 0.5, margin = margin(b = 8)),
        axis.title.x = element_text(size = 13, margin = margin(t = 8)),
        axis.title.y = element_text(size = 13, margin = margin(r = 8)),
        axis.text = element_text(size = 10, color = 'black'),
        axis.ticks = element_line(linewidth = 0.4, color = 'black'),
        axis.line = element_line(linewidth = 0.5, color = 'black'),
        strip.background = element_rect(fill = 'white', color = 'black', linewidth = 0.6),
        strip.text = ggtext::element_markdown(size = 11, margin = margin(t = 4, b = 4)),
        panel.spacing.x = unit(0.45, 'cm'),
        panel.spacing.y = unit(0.45, 'cm'),
        legend.position = c(0.93, 0.2),
        legend.text = element_text(size = 10),
        legend.key.width = unit(0.8, 'cm'),
        plot.margin = margin(t = 5, r = 10, b = 5, l = 5)
  )

ggsave(file.path(output_dir, (paste0('distance_distributions')), date_tag, '.pdf'), pD, width = 11, height = 6)
ggsave(file.path(output_dir, (paste0('distance_distributions')), date_tag, '.png'), pD, width = 11, height = 6, dpi = 300)

# Per-PC variance explained by grouping -----------------------------------------------
pc_poc_hc <- Embeddings(poc_hc, reduction)[, 1:n_dims]
r2_for <- function(col) {
  f <- factor(meta[[col]])
  apply(pc_poc_hc, 2, function(pc) summary(lm(pc ~ f))$r.squared)
}
r2_df <- data.frame(
  PC          = factor(paste0('PC', 1:n_dims), levels = paste0('PC', 1:n_dims)),
  capture     = r2_for(batch_col),
  `cell type` = r2_for(celltype_col),
  ratID       = r2_for(ratid_col),
  check.names = FALSE)
write.csv(r2_df, file.path(output_dir, paste0('PC_variance', date_tag, '.csv')), row.names = FALSE)

r2_long <- rbind(
  data.frame(PC = r2_df$PC, source = 'capture',   r2 = r2_df$capture),
  data.frame(PC = r2_df$PC, source = 'cell type', r2 = r2_df$`cell type`),
  data.frame(PC = r2_df$PC, source = 'ratID',     r2 = r2_df$ratID))
r2_long$source <- factor(r2_long$source, levels = c('capture', 'cell type', 'ratID'))

pC <- ggplot(r2_long, aes(PC, r2, fill = source)) +
  geom_col(position = 'dodge') +
  scale_fill_manual(values = pal_group) +
  scale_x_discrete() +
  scale_y_continuous(labels = function(x) paste0(x * 100, '%')) +
  labs(title = 'Variance explained per PC',
       subtitle = paste0('top ', n_dims, ' PCs (note: ratID has 4 levels vs captures)'),
       x = NULL, y = 'variance explained', fill = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
        legend.position = 'top')
ggsave(file.path(output_dir, (paste0('PCvariance', date_tag, '.pdf'))), pC, width = 9, height = 4.5)
ggsave(file.path(output_dir, (paste0('PCvariance', date_tag, '.png'))), pC, width = 9, height = 4.5, dpi = 300)


