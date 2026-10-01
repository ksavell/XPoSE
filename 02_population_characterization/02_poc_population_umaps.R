# Population UMAP generation for POC dataset

# Loading -----------------------------------------------------------------------------
library(Seurat)
library(tidyverse)

source('02_population_characterization/functions/save_dimplot.R')

# Paths -------------------------------------------------------------------------------
input_file <- 'output/01_metadata_clustering_qc/poc_hc_annotated.rds'
input_file2 <- 'output/01_metadata_clustering_qc/poc_combined_annotated.rds'

poc_hc <- readRDS(input_file)
poc_combined <- readRDS(input_file2)

output_dir <- 'output/02_population_characterization/poc'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Settings ----------------------------------------------------------------------------
# Define colors
hex_list <- list(
  'cluster_name' = c('CTL6' = '#2D8CB8',
                     'CTL6b' = '#7044AA',
                     'ETL5' = '#0D5A8B',
                     'ITL23' = '#2EBF5E',
                     'ITL5' = '#50B2AD',  
                     'ITL6' = '#58D2CF',
                     'ITvm' = '#B1DE7D',
                     'NPL5' = '#3E9E64',
                     'Pvalb' = '#B9342C',
                     'PvalbChand' = '#FF2D4E',
                     'Sst' = '#FF9900',
                     'SstChodl' = '#B1B10C',
                     'Sncg' = '#D3408D',
                     'Vip' = '#B864CC',
                     'Lamp5' = '#DA808C'),
  'capture' = c('C1' = '#5A9BC7',
                   'C2' = '#E08A2D'),
  'experience' =  c('HC' = '#C0C0C0',
                    'NC' = '#3B958E'),
  'population' = c('non-active' = '#75C3BC',
              'active' =  '#00675F')
)

# UMAP generation ---------------------------------------------------------------------
save_dimplot(poc_hc, 
             groupby = 'cluster_name',
             file_n = 'poc',
             hex_list = hex_list)

# UMAP by capture per subject ---------------------------------------------------------
save_dimplot(poc_hc, 
             groupby = 'capture',
             splitby = 'ratID',
             file_n = 'poc',
             hex_list = hex_list)

# UMAP by experience ------------------------------------------------------------------
save_dimplot(poc_combined, 
             groupby = 'experience',
             file_n = 'poc',
             hex_list = hex_list)

# UMAP by population ------------------------------------------------------------------
NC <- subset(poc_combined, subset = experience == 'NC')

save_dimplot(NC, 
             groupby = 'population',
             file_n = 'NC', 
             hex_list = hex_list)

# Split by individual
save_dimplot(NC, 
             groupby = 'population',
             splitby = 'ratID',
             file_n = 'NC',
             hex_list = hex_list)

