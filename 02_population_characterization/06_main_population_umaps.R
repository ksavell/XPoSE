# Population UMAP generation for Main dataset

# Loading -----------------------------------------------------------------------------
library(Seurat)
library(tidyverse)

source('02_population_characterization/functions/calc_prop.R')
source('02_population_characterization/functions/make_stdf.R')

# Paths -------------------------------------------------------------------------------
main <- 'output/02_population_characterization/main_annotated.rds'

output_dir <- 'output/02_population_characterization'
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
  'experience' =  c('NT' = '#cfe0af',
                    'RT' = '#B24971',
                    'NC' = '#3B958E',
                    'N' = '#808080'),
  'population' = c('non-active' = '#e37a9e', # RT: Non-active
              'active' =  '#801743'), # RT: Active
  # 'population' = c('non-active' = '#75C3BC', # NC: Non-active
  #             'active' =  '#00675F'), # NC: Active
  'region' = c('dmPFC' = '#7B3294',
               'vmPFC' = '#E66101')
)

# UMAP generation ---------------------------------------------------------------------
save_dimplot(main, 
             groupby = 'cluster_name',
             file_n = 'main',
             hex_list = hex_list)

# Split by experience
save_dimplot(main, 
             groupby = 'cluster_name',
             splitby = 'experience',
             file_n = 'main',
             hex_list = hex_list)

# UMAP by experience and region -------------------------------------------------------
# Subset by experience comparison of choice
experience_set <- subset(
  x = main, subset = experience == 'N' | experience == 'NT')

save_dimplot(
  experience_set,
  groupby = 'experience',
  splitby = 'region',
  file_n = 'N_vs_NT',
  hex_list = hex_list
)

# UMAP by population and region ------------------------------------------------------------
# Subset by population comparison of choice
population_set <- subset(
  x = main, subset = experience == 'NC')

save_dimplot(
  population_set,
  groupby = 'population',
  splitby = 'region',
  file_n = 'NC_active_vs_NC_nonactive',
  hex_list = hex_list
)

