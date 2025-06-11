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
