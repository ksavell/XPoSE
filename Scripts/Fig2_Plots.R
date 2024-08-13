library(ggplot2)
library(dplyr)

load("~/Library/CloudStorage/Box-Box/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/DataForFigures/Robjects/glut.RData")
load("~/Library/CloudStorage/Box-Box/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/DataForFigures/Robjects/gaba.RData")

source("/Users/holmesar/Documents/XPoSE/Scripts/Functions/calc_prop.R")

glut_group <- calc_prop(seur_obj = glut, fact1 = 'ratID',
                        fact2 = 'cluster_name',
                        fact3 = 'group',
                        file_n = "glut_group.csv")

gaba_group <- calc_prop(seur_obj = gaba, fact1 = 'ratID',
                        fact2 = 'cluster_name',
                        fact3 = 'group',
                        file_n = "gaba_group.csv")

excitatory_clusters <- c("CTL6", "CTL6b", "ITL23", "ITL56", "NPL56", "PTL5")
inhibitory_clusters <- c("Lamp5", "Meis2", "Pvalb", "Sst", "SstChodl", "Vip")

glut <- glut_group %>%
  mutate(Identity = case_when(
    cluster_name %in% excitatory_clusters ~ "Excitatory",
    cluster_name %in% inhibitory_clusters ~ "Inhibitory",
    TRUE ~ "unknown"
  ))

gaba <- gaba_group %>%
  mutate(Identity = case_when(
    cluster_name %in% excitatory_clusters ~ "Excitatory",
    cluster_name %in% inhibitory_clusters ~ "Inhibitory",
    TRUE ~ "unknown"
  ))

glut_gaba_combined <- bind_rows(glut, gaba)

glut_gaba_combined$cluster_name <- factor(glut_gaba_combined$cluster_name, levels = c("PTL5", "CTL6", "ITL56", "NPL56", "ITL23", "CTL6b", "Meis2", "Lamp5", "Pvalb", "Sst", "SstChodl", "Vip"))

custom_colors <- c("PTL5" = "#15608f",
                   "CTL6" = "#2c8cb9",
                   "ITL56" = "#64c7c8",
                   "NPL56" = "#3c9e64",
                   "ITL23" = "#41b75f",
                   "CTL6b" = "#704a9e",
                   "Meis2" = "#c52126",
                   "Lamp5" = "#db7f8c",
                   "Pvalb" = "#e66027",
                   "Sst" = "#f8991c",
                   "SstChodl" = "#b0b336",
                   "Vip" = "#a669ab") 

ggplot(glut_gaba_combined) +
  geom_bar(aes(x = Identity, y = count, fill = cluster_name),
           position = "stack",
           stat = "identity") +
  facet_grid(~ ratID, switch = "x") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = custom_colors) +
  theme_classic() +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm"))
