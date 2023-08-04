# Figure 1 plots

# This script generates:
# Panel G: sample tag read averages*

# * denotes that the final plot was made in Prism with output generated in R

# Load packages -----------------------------------------------------------

library(Seurat)
library(plotly)
library(tidyverse)

# Load data ---------------------------------------------------------------

setwd("~/Output")

load("~/MakeObject/combined08032023all.RData")

# Make df of sample tag reads and relevant metadata

source("~/XPoSE/Scripts/Functions/make_stdf.R")

df <- make_stdf(combined)

# split by sample tag and find mean
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

# save the table
write.csv(mean, '~/Output/stReadsMean.csv')
