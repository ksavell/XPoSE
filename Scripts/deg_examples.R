# example gene heatmaps

source("~/XPoSE/Scripts/Functions/find_DEcounts.R")

highlights <- c("Vgf", "Ptprn", "Homer1", "Lingo1", "Fosb", "Nptx2", "Scg2",
                "Bdnf","Arc","Pcdh15","Lingo2")

find_DEcounts(ITL56_dds, coexp_AN, "ITL56", highlights)
find_DEcounts(ITL23_dds, coexp_AN, "ITL23", highlights)
find_DEcounts(CTL6_dds, coexp_AN, "CTL6", highlights)
find_DEcounts(PTL5_dds, coexp_AN, "PTL5", highlights)
find_DEcounts(Pvalb_dds, coexp_AN, "Pvalb", highlights)
find_DEcounts(Sst_dds, coexp_AN, "Sst", highlights)

