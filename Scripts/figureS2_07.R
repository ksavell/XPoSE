# Figure S2
library(Seurat)
library(tidyverse)


cluster_name = c('CTL6' = '#2C8CB9',
                   'PTL5' = '#0A5B8C',
                   'ITL23' = '#41B75F',
                   'ITL5' = '#5DBFC1',  
                   'ITL6' = '#3A8F87',
                   'NPL56' = '#3C9E64',
                   'CTL6b' = '#6F499D',
                   'Pvalb' = '#E66027',
                   'Sst' = '#F8991D',
                   'Meis2' = '#C52126',
                   'Vip' = '#A669AB',
                   'Lamp5' = '#DB808C',
                   'SstChodl' = '#B0B235',
                   'PvalbChand' = '#AD6C49')

cl_order <- c("ITL23", "ITL5", "ITL6", "CTL6", "CTL6b", "PTL5", "NPL56", 
              "Pvalb", "Sst","PvalbChand","SstChodl","Vip","Lamp5","Meis2") 
all$cluster_name <- factor(all$cluster_name, levels = cl_order)


# FS2A --------------------------------------------------------------------

pdf("FS2A_genesexpressed.pdf",
    width = 3,
    height = 3)
VlnPlot(
  all, 
  features = "nFeature_RNA", 
  pt.size = 0, 
  cols = cluster_name
) + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        title.position = "none") +
  coord_flip()
dev.off()


# FS2B --------------------------------------------------------------------

pdf("FS2B_transcriptsexpressed.pdf",
    width = 3,
    height = 3)
VlnPlot(
  all, 
  features = "nCount_RNA", 
  pt.size = 0, 
  cols = cluster_name
) + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        title.position = "none") +
  coord_flip()
dev.off()


# added 07/08/2026

## Stacked bar plot + Prism-ready tables of Sample_tag contribution per cluster

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

load("input/HC_annotated_07012026.RData")

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
