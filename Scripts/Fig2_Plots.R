library(ggplot2)
library(dplyr)

data <- read.csv("/Users/holmesar/Desktop/Excitatory_Inhibitory_Combined.csv", header = TRUE)

data$Cluster <- factor(data$Cluster, levels = c("PTL5", "CTL6", "ITL56", "NPL56", "ITL23", "CTL6b", "Meis2", "Lamp5", "Pvalb", "Sst", "SstChodl", "Vip"))

data <- data %>%
  mutate(Context = case_when(
    Rat %in% c("R1", "R2", "R3", "R4") ~ "Homecage",
    Rat %in% c("Rat5", "Rat6", "Rat7", "Rat8") ~ "Novel context",
    TRUE ~ "Other"  # In case there are levels not covered by the conditions
  ))


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

ggplot(data) +
  geom_bar(aes(x = Identity, y = Count, fill = Cluster),
           position = "stack",
           stat = "identity") +
  facet_grid(~ Rat, switch = "x") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  scale_fill_manual(values = custom_colors) +
  theme_classic() +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm"))

