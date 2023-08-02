# Figure 1 plots

# This script generates:
# Panel G: sample tag read averages

# Load packages -----------------------------------------------------------
library(Seurat)
library(plotly)
library(tidyverse)

# Load data ---------------------------------------------------------------

# Combine seurat objects

combined <- merge(glut, gaba)

# Make df of sample tag reads and relevant metadata

make_stdf <- function(seur_obj) {
        IDs <- data.frame(seur_obj@assays[["RNA"]]@data@Dimnames[[2]],
                          cart = seur_obj@meta.data[["orig.ident"]])
        
        IDs$cart <- seur_obj@meta.data[["orig.ident"]]
        IDs$rat <- seur_obj@meta.data[["ratID"]]
        IDs$st <- seur_obj@meta.data[["Sample_tag"]]
        IDs$st2 <- seur_obj@meta.data[["SampleTag02_reads"]]
        IDs$st3 <- seur_obj@meta.data[["SampleTag03_reads"]]
        IDs$st4 <- seur_obj@meta.data[["SampleTag04_reads"]]
        IDs$st5 <- seur_obj@meta.data[["SampleTag05_reads"]]
        IDs$st6 <- seur_obj@meta.data[["SampleTag06_reads"]]
        IDs$st7 <- seur_obj@meta.data[["SampleTag07_reads"]]
        IDs$st8 <- seur_obj@meta.data[["SampleTag08_reads"]]
        IDs$st9 <- seur_obj@meta.data[["SampleTag09_reads"]]
        
        rownames(IDs) <- IDs$seur_obj.assays...RNA....data.Dimnames..2..
        #IDs <- IDs %>% select(cart, rat, st)
        
        return(IDs)
}

df <- make_stdf(combined)

#Check sample tag reads for nuclei in rest of figures
IDslist <- split(df, df$st)

mean <- sapply(IDslist, function(x) {
        # Find the names of the numeric columns
        numeric_cols <- sapply(x, is.numeric)
        
        # Ignore character columns
        numeric_cols <- numeric_cols & !sapply(x, is.character)
        
        # Find the means of the numeric columns
        means <- apply(x[, numeric_cols], 2, mean)
        
        # Return a named vector of means
        names(means) <- names(x)[numeric_cols]
        return(means)
})

write.csv(mean, '~/Output/stReadsMean.csv')
