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

load("input/hc_annotated_07202026.RData")
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
# Sample_tag-to-ratID lookup
sample_map <- unique(data.frame(
  Sample_tag = as.character(all$Sample_tag),
  ratID = as.character(all$ratID)
))

if (anyDuplicated(sample_map$Sample_tag))
  stop("At least one Sample_tag maps to multiple ratID values.")

sample_map$sample_number <- as.integer(sub("^SampleTag([0-9]{2})_mm$", "\\1",
                                           sample_map$Sample_tag))
sample_map <- sample_map[order(sample_map$sample_number), ]

sample_order <- sample_map$Sample_tag
rat_labels <- setNames(sample_map$ratID, sample_map$Sample_tag)

# Restrained sample palette; change hex codes if needed
sample_cols <- setNames(
  c("black", "gray20", "gray40", "gray60",
    "#F28E2B", "#E15759", "#B07AA1", "#9C755F"),
  sample_order
)

counts$cluster_name <- factor(counts$cluster_name, levels = rev(cl_order))
counts$Sample_tag <- factor(counts$Sample_tag, levels = rev(sample_order))

p_stack <- ggplot(counts, aes(x = cluster_name, y = percent, fill = Sample_tag)) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.25,
           position = position_stack(reverse = TRUE)) +
  scale_x_discrete(labels = function(x)
    sprintf("<span style='color:%s'>%s</span>", cluster_cols[x], x)) +
  scale_y_continuous(
    breaks = seq(0, 100, 25),
    expand = expansion(mult = c(0, 0)),
    labels = function(x) paste0(x, "%")) +
  scale_fill_manual(values = sample_cols, breaks = sample_order,
                    labels = rat_labels, drop = FALSE) +
  labs(x = NULL, y = "Sample composition", fill = NULL) +
  coord_flip() +
  theme_classic(base_family = "Arial") +
  theme(
    axis.text.x = element_text(size = 7, colour = "black", angle = 0),
    axis.text.y = ggtext::element_markdown(size = 7),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_blank(),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

quartz(type = "pdf", file = paste0("sample_composition_stackedbar", date_tag, ".pdf"),
       width = 2, height = 1.75, family = "Arial")
print(p_stack)
dev.off()

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
 

# Cluster matching to subclass and supertype  ---------------------------------
load("input/hc_annotated_07202026.RData")

supp_table <- obj@meta.data %>%
  count(cluster_name, cluster_name_subclass, cluster_name_supertype,
        name = "n_cells") %>%
  mutate(
    cluster_name = factor(cluster_name, levels = cl_order),
    percent_cells = round(100 * n_cells / sum(n_cells), 2)
  ) %>%
  arrange(cluster_name, cluster_name_subclass, cluster_name_supertype) %>%
  mutate(cluster_name = as.character(cluster_name)) %>%
  select(
    Cluster = cluster_name,
    Subclass = cluster_name_subclass,
    Supertype = cluster_name_supertype,
    `Percent of cells` = percent_cells
  )

write.csv(supp_table, "Supplemental_cluster_hierarchy.csv",
          row.names = FALSE, na = "")



# no integration needed ---------------------------------------------------

library(ggh4x)   # per-facet strip colors

## ---------------- settings ----------------

batch_col    <- "orig.ident"
celltype_col <- "cluster_name"
ratid_col    <- "ratID"
reduction    <- "pca"
n_dims       <- 40
max_pairs    <- 20000
min_cells    <- 20
seed         <- 42
date_tag     <- format(Sys.Date(), "_%m%d%Y")

## ---- fixed cluster palette + biological order ----
cluster_cols <- c('ITL23'='#2EBF5E','ITL5'='#50B2AD','ITL6'='#58D2CF','ITvm'='#B1DE7D',
                  'CTL6'='#2D8CB8','CTL6b'='#7044AA','ETL5'='#0D5A8B','NPL5'='#3E9E64',
                  'Pvalb'='#B9342C','PvalbChand'='#FF2D4E','Sst'='#FF9900','SstChodl'='#B1B10C',
                  'Sncg'='#D3408D','Vip'='#B864CC','Lamp5'='#DA808C')
cl_order <- c("ITL23","ITL5","ITL6","ITvm","CTL6","CTL6b","ETL5","NPL5",
              "Pvalb","Sst","PvalbChand","SstChodl","Vip","Lamp5","Sncg")

## ---- within/between colors (meaningful contrast = between gets accent) ----
wb_cols <- c("within" = "grey60", "between" = "black")

set.seed(seed)
emb  <- Embeddings(obj, reduction)[, 1:n_dims]
meta <- obj@meta.data
stopifnot(all(rownames(emb) == rownames(meta)))

## ============================================================
## Shared minimalist theme (Arial, thin lines, no grid/bg)
## ============================================================
theme_pub <- function() {
  theme_classic(base_size = 7, base_family = "Arial") +
    theme(
      axis.line       = element_line(linewidth = 0.5 / .pt),
      axis.ticks      = element_line(linewidth = 0.5 / .pt),
      axis.text       = element_text(size = 7, color = "black"),
      axis.title      = element_text(size = 8, color = "black"),
      plot.title      = element_blank(),
      plot.subtitle   = element_blank(),
      legend.title    = element_text(size = 7),
      legend.text     = element_text(size = 7),
      legend.key.size = unit(6, "pt"),
      strip.background = element_blank(),
      strip.text      = element_text(size = 7, color = "black"),
      panel.grid      = element_blank(),
      panel.background = element_blank(),
      plot.background  = element_blank()
    )
}
outline_w <- 1 / .pt   # 0.25pt data outlines

## ============================================================
## PLOT 1 — within vs between-capture distance distributions
##          faceted by cell type, biological order, on-palette
## ============================================================
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

raw_list <- list()
for (type in cl_order) {
  in_type <- which(ct == type)
  idx1 <- in_type[batch[in_type] == caps[1]]
  idx2 <- in_type[batch[in_type] == caps[2]]
  if (length(idx1) < min_cells || length(idx2) < min_cells) next
  w1 <- sample_pairs(idx1, n = max_pairs %/% 2, same_pool = TRUE)
  w2 <- sample_pairs(idx2, n = max_pairs %/% 2, same_pool = TRUE)
  bt <- sample_pairs(idx1, idx2, n = max_pairs, same_pool = FALSE)
  raw_list[[type]] <- rbind(
    data.frame(cell_type = type, comparison = "within",  dist = pair_dist(c(w1$a, w2$a), c(w1$b, w2$b))),
    data.frame(cell_type = type, comparison = "between", dist = pair_dist(bt$a, bt$b)))
}
raw_df <- do.call(rbind, raw_list)
raw_df$cell_type  <- factor(raw_df$cell_type, levels = cl_order)
raw_df$comparison <- factor(raw_df$comparison, levels = c("within", "between"))
write.csv(raw_df, paste0("distance_distributions_data", date_tag, ".csv"), row.names = FALSE)

# per-facet strip colors matching each cluster's palette color (biological order)
present   <- cl_order[cl_order %in% levels(droplevels(raw_df$cell_type))]
strip_fg  <- cluster_cols[present]
strip_spec <- strip_themed(
  text_x = elem_list_text(colour = strip_fg, face = "bold"),
  background_x = elem_list_rect(fill = NA, colour = NA))

# color by within/between; consistent across all facets
p1 <- ggplot(raw_df, aes(dist, color = comparison)) +
  geom_density(fill = NA, linewidth = outline_w) +
  facet_wrap2(~ cell_type, scales = "free", ncol = 7, strip = strip_spec) +
  scale_color_manual(values = wb_cols, name = NULL) +
  scale_x_continuous() +
  labs(x = "Euclidean distance (PCA)", y = "Density") +
  theme_pub() +
  theme(legend.position = c(0.9, 0.12))

quartz(type = "pdf", file = paste0("distance_distributions", date_tag, ".pdf"),
       width = 6.5, height = 3, family = "Arial")
print(p1); dev.off()

## ============================================================
## PLOT 2 — per-PC variance explained by grouping
##          cell type colored, capture/ratID neutral gray
## ============================================================
pc_all <- Embeddings(obj, reduction)[, 1:n_dims]
r2_for <- function(col) {
  f <- factor(meta[[col]])
  apply(pc_all, 2, function(pc) summary(lm(pc ~ f))$r.squared)
}
r2_long <- rbind(
  data.frame(PC = 1:n_dims, source = "capture",   r2 = r2_for(batch_col)),
  data.frame(PC = 1:n_dims, source = "cell type", r2 = r2_for(celltype_col)),
  data.frame(PC = 1:n_dims, source = "ratID",     r2 = r2_for(ratid_col)))
r2_long$source <- factor(r2_long$source, levels = c("cell type", "ratID", "capture"))
write.csv(r2_long, paste0("PCvariance_data", date_tag, ".csv"), row.names = FALSE)

# cell type colored, capture/ratID neutral gray
grp_cols <- c("capture" = "grey60", "cell type" = "#2D8CB8", "ratID" = "grey30")

p2 <- ggplot(r2_long, aes(PC, r2 * 100, fill = source)) +
  geom_col(position = "dodge", width = 0.8, linewidth = outline_w, color = NA) +
  scale_fill_manual(values = grp_cols, name = NULL) +
  scale_x_continuous(breaks = seq(0, n_dims, 10), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Principal component", y = "Variance explained (%)") +
  theme_pub() +
  theme(legend.position = c(0.82, 0.85))

quartz(type = "pdf", file = paste0("PCvariance", date_tag, ".pdf"),
       width = 2, height = 1.5, family = "Arial")
print(p2); dev.off()

cat("Saved (", date_tag, "):\n  distance_distributions", date_tag, ".pdf\n  PCvariance",
    date_tag, ".pdf\n  + matching _data CSVs\n", sep = "")

