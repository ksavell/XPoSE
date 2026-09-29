# Population proportion calculations and statistics for Main dataset

# Loading -----------------------------------------------------------------------------

library(Seurat)
library(dplyr)
library(tidyr)
library(purrr)

source('functions/calc_prop.R')

# Load in clustered Main object that is output of createobject_01.R
load('dmvmpfc_annotated_07162026.RData')

# Cluster proportions by experience ---------------------------------------------------

# Subset by experience and region
region_set <- subset(x = all, subset = (experience == 'N'| experience == 'NT') | region == 'dmPFC')

clust_prop_region <- calc_prop(region_set, 
                                     fact1 = 'ratID',
                                     fact2 = 'cluster_name',
                                     fact3 = 'experience')

write.csv(clust_prop_region, file = 'output/main_experience_dmPFC.csv')

# Cluster proportions by region -------------------------------------------------------

# Subset by experience 
region_set <- subset(x = all, subset = experience == 'N')

clust_prop_exp_region <- calc_prop(region_set, 
                                     fact1 = 'ratID',
                                     fact2 = 'cluster_name',
                                     fact3 = 'region') 

write.csv(clust_prop_exp_region, file = 'output/main_experience_by_region.csv')

# EXPERIENCE COMPARISONS ==============================================================
# Build per-rat/experience proportions ------------------------------------------------

experience_comparison <- subset(all, subset = experience == 'N' | experience == 'NT')

experience_prop_df <- experience_comparison@meta.data %>%          
  count(ratID, experience, cluster_name, name = 'n') %>%
  group_by(ratID, experience) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  # make sure zero-count clusters become prop = 0, not missing
  complete(nesting(ratID, experience), cluster_name, fill = list(n = 0, prop = 0))

# Test runners ------------------------------------------------------------------------

# UNPAIRED: experience vs experience. Independent rats, so Welch t-test.
run_unpaired <- function(df, g1, g2) {
  df %>%
    filter(experience %in% c(g1, g2)) %>%
    group_by(cluster_name) %>%
    summarise(
      t = list(tryCatch(
        t.test(prop ~ experience, var.equal = FALSE),   # Welch
        error = function(e) NULL)),
      mean_g1 = mean(prop[experience == g1]),
      mean_g2 = mean(prop[experience == g2]),
      .groups = 'drop') %>%
    mutate(
      p    = map_dbl(t, ~ if (is.null(.x)) NA_real_ else .x$p.value),
      stat = map_dbl(t, ~ if (is.null(.x)) NA_real_ else unname(.x$statistic))) %>%
    select(-t) %>%
    mutate(comparison = paste(g1, 'vs', g2))
}

# Run the comparisons -----------------------------------------------------------

res1 <- run_unpaired(experience_prop_df, 'N', 'NT')

# FDR correction
# Correcting across clusters WITHIN each comparison (15 tests each).
add_fdr <- function(x) mutate(x, padj = p.adjust(p, method = 'BH'))
res1 <- add_fdr(res1)
results <- bind_rows(res1) %>%
  arrange(comparison, padj)

print(results, n = Inf)
write.csv(results, 'output/N_vs_NT_proportion_stats_summary.csv')

# WITHIN-EXPERIENCE POPULATION COMPARISONS ============================================
# Build per-rat/group proportions -----------------------------------------------------

population_comparison <- subset(all, subset = experience == 'RT')

population_prop_df <- population_comparison@meta.data %>%          
  count(ratID, population, cluster_name, name = 'n') %>%
  group_by(ratID, population) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  # make sure zero-count clusters become prop = 0, not missing
  complete(nesting(ratID, population), cluster_name, fill = list(n = 0, prop = 0))

# Test runners ------------------------------------------------------------------------

# PAIRED: active vs non-active within the SAME rat.
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

# Run the comparisons -----------------------------------------------------------

res1 <- run_paired(population_prop_df, 'active', 'non-active')

# FDR correction
# Correcting across clusters WITHIN each comparison (15 tests each).
add_fdr <- function(x) mutate(x, padj = p.adjust(p, method = 'BH'))
res1 <- add_fdr(res1)
results <- bind_rows(res1) %>%
  arrange(comparison, padj)

print(results, n = Inf)
write.csv(results, 'output/RT_A_vs_RT_NA_proportion_stats_summary.csv')

# REGION COMPARISONS ===================================================================
# Build per-rat/region proportions -----------------------------------------------------

region_comparison <- subset(all, subset = experience == 'N')

region_prop_df <- region_comparison@meta.data %>%          
  count(ratID, region, cluster_name, name = 'n') %>%
  group_by(ratID, region) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  # make sure zero-count clusters become prop = 0, not missing
  complete(nesting(ratID, region), cluster_name, fill = list(n = 0, prop = 0))

# Test runners ------------------------------------------------------------------------

# PAIRED: dmPFC vs vmPFC within the SAME rat.
run_paired <- function(df, g1, g2) {
  wide <- df %>%
    filter(region %in% c(g1, g2)) %>%
    select(ratID, region, cluster_name, prop) %>%
    pivot_wider(names_from = region, values_from = prop)
  
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

# Run the comparisons -----------------------------------------------------------

res1 <- run_paired(region_prop_df, 'dmPFC', 'vmPFC')

# FDR correction
# Correcting across clusters WITHIN each comparison (15 tests each).
add_fdr <- function(x) mutate(x, padj = p.adjust(p, method = 'BH'))
res1 <- add_fdr(res1)
results <- bind_rows(res1) %>%
  arrange(comparison, padj)

print(results, n = Inf)
write.csv(results, 'output/N_dmPFC_vs_vmPFC_proportion_stats_summary.csv')

