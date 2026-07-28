# Figure S2
library(Seurat)
library(tidyverse)
library(ggtext)
library(systemfonts)
library(scales)

cluster_cols = c('ITL23' = '#2EBF5E',
                 'ITL5' =	'#50B2AD',
                 'ITL6' =	'#58D2CF',
                 'ITvm' =	'#B1DE7D',
                 'CTL6' =	'#2D8CB8',
                 'CTL6b' =	'#7044AA',
                 'ETL5' =	'#0D5A8B',
                 'NPL5' =	'#3E9E64',
                 'Pvalb' =	'#B9342C',
                 'PvalbChand' =	'#FF2D4E',
                 'Sst' =	'#FF9900',
                 'SstChodl' =	'#B1B10C',
                 'Sncg' =	'#D3408D',
                 'Vip' =	'#B864CC',
                 'Lamp5' =	'#DA808C')

cl_order <- c("ITL23", "ITL5", "ITL6", 'ITvm', "CTL6", "CTL6b", "ETL5", "NPL5", 
              "Pvalb", "Sst","PvalbChand","SstChodl","Vip","Lamp5","Sncg") 
all$cluster_name <- factor(
  as.character(all$cluster_name),
  levels = rev(cl_order)
)


# FS2A --------------------------------------------------------------------

p <- VlnPlot(
  all,
  features = "nFeature_RNA",
  group.by = "cluster_name",
  cols = cluster_cols,
  sort = FALSE,
  pt.size = 0,
  combine = FALSE)[[1]] +
  scale_x_discrete(
    labels = function(x) {sprintf("<span style='color:%s'>%s</span>",
                                  cluster_cols[x],
                                  x)}) +
  scale_y_continuous(
    breaks = function(x) {seq(0,
                              ceiling(max(x) / 1000) * 1000,
                              by = 1000)},
    labels = function(x) x / 1000,
    expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = NULL,
    y = "Genes expressed (\u00D710\u00B3)") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 7,
                                   family = "Arial",
                                   angle = 0),
    axis.text.y = ggtext::element_markdown(size = 7,
                                           family = "Arial"),
    axis.title.x = element_text(size = 8,
                                family = "Arial"),
    axis.title.y = element_blank(),
    plot.title = element_blank(),
    
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5)) +
  coord_flip()

# Set violin outlines to 0.25
violin_layers <- vapply(
  p$layers,
  function(layer) inherits(layer$geom, "GeomViolin"),
  logical(1)
)

for (i in which(violin_layers)) {p$layers[[i]]$aes_params$linewidth <- 0.25}

quartz(
  type = "pdf",
  file = "FS2A_genesexpressed.pdf",
  width = 1.5,
  height = 2.5,
  family = "Arial"
)

print(p)
dev.off()


# FS2B --------------------------------------------------------------------

p2 <- VlnPlot(
  all,
  features = "nCount_RNA",
  group.by = "cluster_name",
  cols = cluster_cols,
  sort = FALSE,
  pt.size = 0,
  combine = FALSE)[[1]] +
  scale_x_discrete(labels = function(x) {
      sprintf(
        "<span style='color:%s'>%s</span>",
        cluster_cols[x], x)}) +
  scale_y_continuous(
    breaks = seq(0, 15000, by = 5000),
    labels = function(x) x / 1000,
    expand = expansion(mult = c(0, 0))) +
  labs(
    x = NULL,
    y = "Transcripts expressed (\u00D710\u00B3)") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 7,
                               family = "Arial",
                               angle = 0),
    axis.text.y = ggtext::element_markdown(size = 7,
                                           family = "Arial"),
    axis.title.x = element_text(size = 8,
                                family = "Arial"),
    axis.title.y = element_blank(),
    plot.title = element_blank(),
    
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5)) +
  coord_flip(
    ylim = c(0, 15000))

violin_layers <- vapply(
  p2$layers,
  function(layer) inherits(layer$geom, "GeomViolin"),
  logical(1))

for (i in which(violin_layers)) {p2$layers[[i]]$aes_params$linewidth <- 0.25}

quartz(
  type = "pdf",
  file = "FS2B_transcriptsexpressed.pdf",
  width = 1.5,
  height = 2.5,
  family = "Arial"
)

print(p2)
dev.off()


# added 07/08/2026

## Stacked bar plot + Prism-ready tables of Sample_tag contribution per cluster

load("~/Projects/XPoSE/hc_annotated_07202026.RData")
all <- obj

date_tag <- format(Sys.Date(), "_%m%d%Y")
n_tags   <- length(unique(as.character(all$Sample_tag)))
expected <- 1 / n_tags

## counts + proportions per cluster
counts <- all@meta.data %>%
  dplyr::select(cluster_name, Sample_tag) %>%
  dplyr::mutate(dplyr::across(everything(), as.character)) %>%
  dplyr::count(cluster_name, Sample_tag, name = "n") %>%
  tidyr::complete(cluster_name, Sample_tag, fill = list(n = 0)) %>%
  dplyr::group_by(cluster_name) %>%
  dplyr::mutate(cluster_total = sum(n),
                prop = n / cluster_total,
                percent = prop * 100) %>%
  dplyr::ungroup()

## ---- Prism-ready wide tables (clusters x tags) ----
write_prism_wide <- function(value_col, fname) {
  wide <- counts %>%
    dplyr::select(cluster_name, Sample_tag, dplyr::all_of(value_col)) %>%
    tidyr::pivot_wider(names_from = Sample_tag,
                       values_from = dplyr::all_of(value_col)) %>%
    dplyr::arrange(cluster_name)
  write.csv(wide, file.path(paste0(fname, date_tag, ".csv")),
            row.names = FALSE)
}

write_prism_wide("percent", "prism_stackedbar_percent")
# write_prism_wide("prop",    "prism_stackedbar_proportion")
# write_prism_wide("n",       "prism_stackedbar_counts")


# FS2 sort order ----------------------------------------------------------
load("~/Projects/XPoSE/combined_annotated_07202026.RData")

# Create the Sample_tag-to-ratID lookup
sample_map <- unique(
  data.frame(
    Sample_tag = as.character(obj$Sample_tag),
    ratID = as.character(obj$ratID)
  )
)

# Confirm that each Sample_tag maps to only one ratID
if (anyDuplicated(sample_map$Sample_tag)) {
  stop("At least one Sample_tag is associated with multiple ratID values.")
}

# Extract the number from SampleTag##_mm
sample_map$sample_number <- as.integer(
  sub(
    "^SampleTag([0-9]{2})_mm$",
    "\\1",
    sample_map$Sample_tag
  )
)

if (anyNA(sample_map$sample_number)) {
  stop("Some Sample_tag values do not match the expected SampleTag##_mm format.")
}

# Sort SampleTag02_mm through SampleTag09_mm
sample_map <- sample_map[
  order(sample_map$sample_number),
]

sample_order <- sample_map$Sample_tag

# Named label vector: Sample_tag -> ratID
rat_labels <- setNames(
  sample_map$ratID,
  sample_map$Sample_tag
)

# Reverse factor levels so 02 appears at the top after coord_flip()
obj$Sample_tag <- factor(
  as.character(obj$Sample_tag),
  levels = rev(sample_order)
)

# Same violin color for every sample
sample_cols <- setNames(
  rep("#808080", length(sample_order)),
  sample_order
)

p_counts <- VlnPlot(
  obj,
  features = "nCount_RNA",
  group.by = "Sample_tag",
  cols = sample_cols,
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
    y = "Transcripts expressed (\u00D710\u00B3)"
  ) +
  theme(
    legend.position = "none",
    
    axis.text.x = element_text(
      size = 7,
      family = "Arial",
      angle = 0,
      hjust = 0.5
    ),
    
    # These are now ratID labels
    axis.text.y = element_text(
      size = 7,
      family = "Arial",
      color = "black"
    ),
    
    axis.title.x = element_text(
      size = 8,
      family = "Arial"
    ),
    
    axis.title.y = element_blank(),
    plot.title = element_blank(),
    
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5)
  ) +
  coord_flip(
    ylim = c(0, 15000)
  )

violin_layers <- vapply(
  p_counts$layers,
  function(layer) inherits(layer$geom, "GeomViolin"),
  logical(1)
)

for (i in which(violin_layers)) {
  p_counts$layers[[i]]$aes_params$linewidth <- 0.25
}

quartz(
  type = "pdf",
  file = "sample_transcriptsexpressed.pdf",
  width = 1.2,
  height = 2.5,
  family = "Arial"
)

print(p_counts)
dev.off()

p_genes <- VlnPlot(
  obj,
  features = "nFeature_RNA",
  group.by = "Sample_tag",
  cols = sample_cols,
  sort = FALSE,
  pt.size = 0,
  combine = FALSE
)[[1]] +
  scale_x_discrete(
    labels = rat_labels
  ) +
  scale_y_continuous(
    breaks = seq(0, 5000, by = 1000),
    labels = function(x) x / 1000,
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = "Genes expressed (\u00D710\u00B3)"
  ) +
  theme(
    legend.position = "none",
    
    axis.text.x = element_text(
      size = 7,
      family = "Arial",
      angle = 0,
      hjust = 0.5
    ),
    
    # These are now ratID labels
    axis.text.y = element_text(
      size = 7,
      family = "Arial",
      color = "black"
    ),
    
    axis.title.x = element_text(
      size = 8,
      family = "Arial"
    ),
    
    axis.title.y = element_blank(),
    plot.title = element_blank(),
    
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5)
  ) +
  coord_flip(
    ylim = c(1000, 5000)
  )

violin_layers <- vapply(
  p_genes$layers,
  function(layer) inherits(layer$geom, "GeomViolin"),
  logical(1)
)

for (i in which(violin_layers)) {
  p_genes$layers[[i]]$aes_params$linewidth <- 0.25
}

quartz(
  type = "pdf",
  file = "sample_genesexpressed.pdf",
  width = 1.2,
  height = 2.5,
  family = "Arial"
)

print(p_genes)
dev.off()
 