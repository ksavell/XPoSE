# QC for Main dataset

# Loading -----------------------------------------------------------------------------
library(ggplot2)
library(Seurat)
library(dplyr)

# Paths -------------------------------------------------------------------------------
input_file <- 'output/01_metadata_clustering_qc/main_annotated.rds'
main <- readRDS(input_file)

output_dir <- 'output/01_metadata_clustering_qc/main_qc'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Settings ----------------------------------------------------------------------------
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

cluster_labels <- c(
  'ITL23'      = 'IT L2/3',
  'ITL5'       = 'IT L5',
  'ITL6'       = 'IT L6',
  'ITvm'       = 'IT vm',
  'CTL6'       = 'CT L6',
  'CTL6b'      = 'CT L6b',
  'ETL5'       = 'ET L5',
  'NPL5'       = 'NP L5',
  'Pvalb'      = 'Pvalb',
  'Sst'        = 'Sst',
  'PvalbChand' = 'Pvalb Chand',
  'SstChodl'   = 'Sst Chodl',
  'Vip'        = 'Vip',
  'Lamp5'      = 'Lamp5',
  'Sncg'       = 'Sncg'
)

# QC plot by cluster ------------------------------------------------------------------
# Genes expressed
p <- VlnPlot(
  main,
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
    breaks = seq(2000, 8000, by = 2000),
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
  main,
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
    breaks = seq(0, 60000, by = 15000),
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
    ylim = c(0, 60000))

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

main$sort_day <- dplyr::case_when(
  main$capture %in% c('dmPFC1', 'dmPFC2') ~ 'Sort_day_1',
  main$capture %in% c('dmPFC3', 'dmPFC4') ~ 'Sort_day_2',
  main$capture %in% c('vmPFC1', 'vmPFC2') ~ 'Sort_day_3',
  main$capture %in% c('vmPFC3', 'vmPFC4') ~ 'Sort_day_4',
  TRUE ~ NA_character_
)

# Store sort_day in the intended order
main$sort_day <- factor(
  main$sort_day,
  levels = paste0('Sort_day_', 1:4)
)

sortdays <- unique(main$sort_day)

# Use the same limits for every sort-day plot
gene_max <- 8000  # ceiling(max(main$nFeature_RNA, na.rm = TRUE) / 1000) * 1000
transcript_max <- 60000

for (day in sortdays) {
  
  # Subset the Seurat object
  main <- subset(main, subset = sort_day == day)
  
  # Create the xpose_tag-to-ratID lookup
  xpose_map <- unique(
    data.frame(
      xpose_tag = as.character(main$xpose_tag),
      ratID = as.character(main$ratID)
    )
  )
  
  # Confirm that each xpose_tag maps to only one ratID
  if (anyDuplicated(xpose_map$xpose_tag)) {
    stop('At least one xpose_tag is associated with multiple ratID values.')
  }
  
  # Confirm xpose_tag format
  valid_tags <- grepl(
    '^xpose_tag_(0[1-9]|1[0-2])_mm$',
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
      '^xpose_tag_(0[1-9]|1[0-2])_mm$',
      '\\1',
      xpose_map$xpose_tag
    )
  )
  
  # Sort xpose tags numerically
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
  main$xpose_tag <- factor(
    as.character(main$xpose_tag),
    levels = rev(xpose_order)
  )
  
  # Same violin color for every xpose tag
  xpose_cols <- setNames(
    rep('#808080', length(xpose_order)),
    xpose_order
  )
  
# Genes expressed
  p_genes <- VlnPlot(
    main,
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
      breaks = seq(2000, 8000, by = 2000),
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
      ylim = c(2000, gene_max)
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
}
  
# Transcripts expressed
  p_counts <- VlnPlot(
    main,
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
      breaks = seq(0, 60000, by = 15000),
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
