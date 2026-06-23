# Assign metadata and create combined initial seurat object
# Authors:  Katherine Savell and Padmashri Saravanan

# Info --------------------------------------------------------------------

# This script takes 3 input files per cartridge from the alignment pipeline
# in SB
# count table RSEC_MolsPerCell.csv
# sample tag reads Sample_Tag_ReadsPerCell.csv
# sample tag calls Sample_Tag_Calls.csv

# Load packages -----------------------------------------------------------
library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)

# Date stamp for outputs --------------------------------------------------
date_stamp <- format(Sys.Date(), "_%m%d%Y")

# Load metadata lookup table ----------------------------------------------
# CSV columns: cart, sample_tag, ratID, sex, experience, group
# cart values: dmPFC1, dmPFC2, dmPFC3, dmPFC4
meta_lookup <- read.csv("input/metadata_lookup.csv", stringsAsFactors = FALSE)

# Trim any whitespace from character columns
meta_lookup$cart       <- trimws(meta_lookup$cart)
meta_lookup$sample_tag <- trimws(as.character(meta_lookup$sample_tag))

# Convert integer sample tag to BD format (e.g. 1 -> "SampleTag01_mm")
meta_lookup$sample_tag <- ifelse(
  nchar(meta_lookup$sample_tag) == 1,
  paste0("SampleTag0", meta_lookup$sample_tag, "_mm"),
  paste0("SampleTag",  meta_lookup$sample_tag, "_mm")
)

# Read gene count tables
counts_c1 <- read.csv("input/Combined_cart1-good-only_RSEC_MolsPerCell.csv", 
                      skip = 7, row.names = 1)
counts_c2 <- read.csv("input/Combined_cart2-good-only_RSEC_MolsPerCell.csv", 
                      skip = 7, row.names = 1)

# Read sample tag calls
tags_c1 <- read.csv("input/cart1-good-only_Sample_Tag_Calls.csv", 
                    skip = 7, row.names = 1)
tags_c2 <- read.csv("input/cart2-good-only_Sample_Tag_Calls.csv", 
                    skip = 7, row.names = 1)

#Read Sample tag reads
str_c1 <- read.csv("input/cart1-good-only_Sample_Tag_ReadsPerCell.csv", 
                       skip = 7, row.names = 1)
str_c2 <- read.csv("input/cart2-good-only_Sample_Tag_ReadsPerCell.csv", 
                       skip = 7, row.names = 1)

# Reorder all input files -------------------------------------------------
objects <- c("counts_c1","counts_c2",
             "tags_c1","tags_c2",
             "str_c1","str_c2")

invisible(lapply(objects, function(name) {
  obj <- get(name, envir = .GlobalEnv)
  if (is.data.frame(obj)) {
    assign(name, obj[order(rownames(obj)), ], envir = .GlobalEnv)
  }
}))

# Build named input lists -------------------------------------------------
counts  <- list(c1 = counts_c1, c2 = counts_c2)
tags    <- list(c1 = subset(tags_c1, select = -c(2)),
                c2 = subset(tags_c2, select = -c(2)))
streads <- list(c1 = str_c1, c2 = str_c2)

source("Scripts/Functions/create_seur.R")
source("Scripts/Functions/assign_combine_lookup.R")  # updated lookup-table version

# Create Seurat objects ---------------------------------------------------
# cart_label maps internal key to dmPFC label for orig.ident
cart_map <- c(c1 = "C1", c2 = "C2")

for (i in c("c1","c2")) {
  cart_label <- cart_map[[i]]
  
  assign(paste0("data_", i),
         create_seur(count      = counts[[i]],
                     cart_label = cart_label,
                     tag        = tags[[i]],
                     addSTreads = TRUE,
                     streads    = streads[[i]],
                     STrange    = c(2:9)),
         envir = .GlobalEnv)
}

# Assign metadata from lookup table ---------------------------------------
seur_list <- list(data_c1 = data_c1, data_c2 = data_c2)

for (l in seq_along(seur_list)) {
  seur_list[[l]] <- assign_combine(seur_list[[l]], meta_lookup)
}

# Within-sample doublet detection with scDblFinder ------------------------
# Run per cartridge, stratified by Sample_tag (one rat per tag).
# Automatic doublet rate estimation; residual within-sample rate expected
# to be low since cross-sample doublets already removed via tag calls.
run_scDblFinder <- function(seur_obj, cart_label) {
  sce <- as.SingleCellExperiment(seur_obj)
  
  # Identify and temporarily reassign singleton Sample_tags to "Undetermined"
  tag_counts <- table(seur_obj$Sample_tag)
  singleton_tags <- names(tag_counts[tag_counts == 1])
  
  temp_tags <- seur_obj$Sample_tag
  temp_tags[temp_tags %in% singleton_tags] <- "Undetermined"
  sce$temp_sample_tag <- temp_tags
  
  set.seed(42)
  sce <- scDblFinder(sce, samples = "temp_sample_tag")
  
  seur_obj$doublet       <- sce$scDblFinder.class
  seur_obj$doublet_score <- sce$scDblFinder.score
  
  message(cart_label, " — doublet summary by Sample_tag:")
  print(table(seur_obj$doublet, seur_obj$Sample_tag))
  
  return(seur_obj)
}

cart_labels <- c("c1", "c2")
for (l in seq_along(seur_list)) {
  seur_list[[l]] <- run_scDblFinder(seur_list[[l]], cart_labels[l])
}

# Plots -------------------------------------------------------------------

plot_doublets <- function(seur_obj, cart_label) {
  
  seur_obj <- NormalizeData(seur_obj, verbose = FALSE)
  seur_obj <- FindVariableFeatures(seur_obj, verbose = FALSE)
  seur_obj <- ScaleData(seur_obj, verbose = FALSE)
  seur_obj <- RunPCA(seur_obj, verbose = FALSE)
  seur_obj <- RunUMAP(seur_obj, dims = 1:20, verbose = FALSE)
  
  p1 <- FeaturePlot(seur_obj, features = "doublet_score", pt.size = 0.3) +
    scale_color_viridis_c(option = "magma") +
    ggtitle(paste0(cart_label, " — Doublet Score")) +
    theme_classic(base_size = 11)
  
  p2 <- DimPlot(seur_obj, group.by = "doublet", pt.size = 0.3,
                cols = c(singlet = "grey80", doublet = "#E63946")) +
    ggtitle(paste0(cart_label, " — Doublet Calls")) +
    theme_classic(base_size = 11)
  
  meta   <- seur_obj@meta.data
  bar_df <- as.data.frame(table(Sample_tag = meta$Sample_tag,
                                doublet    = meta$doublet))
  
  p3 <- ggplot(bar_df, aes(x = Sample_tag, y = Freq, fill = doublet)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = c(singlet = "grey80", doublet = "#E63946")) +
    scale_y_continuous(labels = scales::percent) +
    labs(title = paste0(cart_label, " — Doublet Rate per Sample Tag"),
         x = "Sample Tag (Individual)", y = "Proportion", fill = "Call") +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(list(score_umap = p1, call_umap = p2, bar = p3))
}

all_plots <- mapply(plot_doublets, seur_list, cart_labels, SIMPLIFY = FALSE)

# Cross-cartridge reproducibility panel -----------------------------------
meta_all <- do.call(rbind, lapply(seur_list, function(s) s@meta.data))

cross_df <- as.data.frame(table(Sample_tag = meta_all$Sample_tag,
                                doublet    = meta_all$doublet,
                                cart       = meta_all$orig.ident))

p_cross <- ggplot(cross_df, aes(x = Sample_tag, y = Freq, fill = doublet)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c(singlet = "grey80", doublet = "#E63946")) +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~cart, ncol = 2) +
  labs(title = "Within-Sample Doublet Rate by Cartridge (scDblFinder)",
       x = "Sample Tag (Individual)", y = "Proportion", fill = "Call") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank())

# Save plots --------------------------------------------------------------
dir.create("Plots", showWarnings = FALSE)
pdf(paste0("Plots/scDblFinder_QC", date_stamp, ".pdf"), width = 14, height = 10)

for (l in seq_along(all_plots)) {
  p <- all_plots[[l]]
  print((p$score_umap | p$call_umap) / p$bar +
          plot_annotation(title = paste0(cart_labels[l], " — Within-Sample Doublet Detection")))
}
print(p_cross)

dev.off()

message("Plots saved to Plots/scDblFinder_QC", date_stamp, ".pdf")

# Doublet summary table ---------------------------------------------------
summary_table <- do.call(rbind, lapply(seq_along(seur_list), function(l) {
  meta <- seur_list[[l]]@meta.data
  cart <- cart_labels[l]
  do.call(rbind, lapply(unique(meta$Sample_tag), function(st) {
    cells     <- meta[meta$Sample_tag == st, ]
    n_nuclei  <- nrow(cells)
    n_doublet <- sum(cells$doublet == "doublet")
    data.frame(
      cart          = cart,
      sample_tag    = st,
      ratID         = unique(cells$ratID),
      n_nuclei      = n_nuclei,
      n_doublets    = n_doublet,
      pct_doublets  = round(n_doublet / n_nuclei * 100, 2),
      stringsAsFactors = FALSE
    )
  }))
}))

write.csv(summary_table,
          file      = paste0("Plots/scDblFinder_summary", date_stamp, ".csv"),
          row.names = FALSE)

message("Summary table saved to Plots/scDblFinder_summary", date_stamp, ".csv")

# Merge all cartridges and save -------------------------------------------
combined_all <- merge(seur_list[[1]],
                      y        = seur_list[2:length(seur_list)],
                      merge.dr = FALSE)

# Save full unfiltered object
save(combined_all, file = paste0("all", date_stamp, ".RData"))

# Remove cross-sample multiplets/undetermined (sample tag based)
combined_all <- subset(combined_all, !(Sample_tag == "Multiplet" |
                                         Sample_tag == "Undetermined" |
                                         Sample_tag == "SampleTag10_mm" |
                                         Sample_tag == "SampleTag12_mm"))

# Remove within-sample doublets flagged by scDblFinder
combined_all <- subset(combined_all, doublet == "singlet")

save(combined_all, file = paste0("all_filtered", date_stamp, ".RData"))

message("Final cell count after all filtering: ", ncol(combined_all))
message("Experience breakdown:")
print(table(combined_all$experience))
