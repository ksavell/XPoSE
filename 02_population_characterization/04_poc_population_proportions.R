# Population proportion calculations and statistics for POC dataset

# Loading -----------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(tidyr)
library(purrr)

source('02_population_characterization/functions/calc_prop.R')

# Paths -------------------------------------------------------------------------------
input_file <- 'output/01_metadata_clustering_qc/poc_hc_annotated.rds'
input_file2 <- 'output/01_metadata_clustering_qc/poc_combined_annotated.rds'

poc_hc <- readRDS(input_file)
poc_combined <- readRDS(input_file2)

output_dir <- 'output/02_population_characterization/poc/'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Cluster proportions by capture ------------------------------------------------------
clust_prop_cart <- calc_prop(poc_hc, 
                             fact1 = 'ratID',
                             fact2 = 'cluster_name',
                             fact3 = 'capture') 

write.csv(clust_prop_cart, file.path(output_dir, 'poc_clust_prop_capture.csv'))

# CAPTURE COMPARISONS =================================================================
# Build per-rat/capture proportions ---------------------------------------------------
prop_df <- poc_hc@meta.data %>%          # or just `poc_combined` if it's already a data.frame
  count(ratID, capture, cluster_name, name = 'n') %>%
  group_by(ratID, capture) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  # make sure zero-count clusters become prop = 0, not missing
  complete(nesting(ratID, capture), cluster_name, fill = list(n = 0, prop = 0))

# Test runners ------------------------------------------------------------------------
# PAIRED: active vs non-active within the SAME NC rat.
run_paired <- function(df, g1, g2) {
  wide <- df %>%
    filter(capture %in% c(g1, g2)) %>%
    select(ratID, capture, cluster_name, prop) %>%
    pivot_wider(names_from = capture, values_from = prop)
  
  wide %>%
    group_by(cluster_name) %>%
    summarise(
      t = list(tryCatch(
        t.test(.data[[g1]], .data[[g2]], paired = TRUE),
        error = function(e) NULL)),
      mean_g1 = mean(.data[[g1]], na.rm = TRUE),
      mean_g2 = mean(.data[[g2]], na.rm = TRUE),
      .groups = 'drop') %>%
    mutate(
      p    = map_dbl(t, ~ if (is.null(.x)) NA_real_ else .x$p.value),
      stat = map_dbl(t, ~ if (is.null(.x)) NA_real_ else unname(.x$statistic))) %>%
    select(-t) %>%
    mutate(comparison = paste(g1, 'vs', g2, '(paired)'))
}

# Run the comparison ------------------------------------------------------------------
res1 <- run_paired(prop_df, 'C1', 'C2')

# FDR correction ----------------------------------------------------------------------
# Correcting across clusters WITHIN each comparison (15 tests each).
add_fdr <- function(x) mutate(x, padj = p.adjust(p, method = 'BH'))

res1 <- add_fdr(res1)

results <- bind_rows(res1) %>%
  arrange(comparison, padj)

print(results, n = Inf)
write.csv(results, file.path(output_dir, 'xpose_population_prop_stats_summary.csv'))

# Cluster proportions by population ----------------------------------------------------
clust_prop_pop <- calc_prop(poc_combined, 
                            fact1 = 'ratID',
                            fact2 = 'cluster_name',
                            fact3 = 'population') 

write.csv(clust_prop_pop, file.path(output_dir, 'poc_clust_prop_population.csv'))

# POPULATION COMPARISONS ===============================================================                     
# Build per-rat/population proportions -------------------------------------------------
prop_df <- poc_combined@meta.data %>%          # or just `poc_combined` if it's already a data.frame
  count(ratID, population, cluster_name, name = 'n') %>%
  group_by(ratID, population) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  # make sure zero-count clusters become prop = 0, not missing
  complete(nesting(ratID, population), cluster_name, fill = list(n = 0, prop = 0)) %>%
  # drop the rat x population cells that don't exist (e.g. HC + active)
  filter(!(population == 'all'        & !grepl('^HC', ratID)),
         !(population %in% c('active','non-active') & !grepl('^NC', ratID)))

# Test runners ------------------------------------------------------------------------
# UNPAIRED: all (HC rats) vs one NC population. Independent rats, so Welch t-test.
run_unpaired <- function(df, g1, g2) {
  df %>%
    filter(population %in% c(g1, g2)) %>%
    group_by(cluster_name) %>%
    summarise(
      t = list(tryCatch(
        t.test(prop ~ population, var.equal = FALSE),   # Welch
        error = function(e) NULL)),
      mean_g1 = mean(prop[population == g1]),
      mean_g2 = mean(prop[population == g2]),
      .groups = 'drop') %>%
    mutate(
      p    = map_dbl(t, ~ if (is.null(.x)) NA_real_ else .x$p.value),
      stat = map_dbl(t, ~ if (is.null(.x)) NA_real_ else unname(.x$statistic))) %>%
    select(-t) %>%
    mutate(comparison = paste(g1, 'vs', g2))
}

# PAIRED: active vs non-active within the SAME NC rat.
run_paired <- function(df, g1, g2) {
  wide <- df %>%
    filter(population %in% c(g1, g2)) %>%
    select(ratID, population, cluster_name, prop) %>%
    pivot_wider(names_from = population, values_from = prop)
  
  wide %>%
    group_by(cluster_name) %>%
    summarise(
      t = list(tryCatch(
        t.test(.data[[g1]], .data[[g2]], paired = TRUE),
        error = function(e) NULL)),
      mean_g1 = mean(.data[[g1]], na.rm = TRUE),
      mean_g2 = mean(.data[[g2]], na.rm = TRUE),
      .groups = 'drop') %>%
    mutate(
      p    = map_dbl(t, ~ if (is.null(.x)) NA_real_ else .x$p.value),
      stat = map_dbl(t, ~ if (is.null(.x)) NA_real_ else unname(.x$statistic))) %>%
    select(-t) %>%
    mutate(comparison = paste(g1, 'vs', g2, '(paired)'))
}

# Run the three comparisons -----------------------------------------------------------
res1 <- run_unpaired(prop_df, 'all', 'non-active')
res2 <- run_unpaired(prop_df, 'all', 'active')
res3 <- run_paired(prop_df, 'active', 'non-active')

# FDR correction ----------------------------------------------------------------------
# Correcting across clusters WITHIN each comparison (15 tests each).
add_fdr <- function(x) mutate(x, padj = p.adjust(p, method = 'BH'))

res1 <- add_fdr(res1)
res2 <- add_fdr(res2)
res3 <- add_fdr(res3)

results <- bind_rows(res1, res2, res3) %>%
  arrange(comparison, padj)

print(results, n = Inf)
write.csv(results, file.path(output_dir, 'population_prop_stats_summary.csv'))
