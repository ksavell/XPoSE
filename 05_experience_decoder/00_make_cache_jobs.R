#!/usr/bin/env Rscript

# ==============================================================================
# XPoSE-seq experience decoder: create nested-cache job manifest
# ==============================================================================
# Run from the repository root. Set DEG_CACHE_ROOT before running on Biowulf if
# you want a cache location other than /data/$USER/hierarchical_paired_deg_cache.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

pb_dir <- "output/03_differential_expression/03_main_de/pseudobulk"
script_dir <- "05_experience_decoder"
jobs_file <- file.path(script_dir, "cache_jobs.tsv")

user_name <- Sys.getenv("USER", unset = "")
default_cache_root <- if (nzchar(user_name)) {
  file.path("/data", user_name, "hierarchical_paired_deg_cache")
} else {
  "output/05_experience_decoder/nested_deg_cache"
}
cache_root <- Sys.getenv("DEG_CACHE_ROOT", unset = default_cache_root)

cluster_sets <- list(
  vmPFC = c("ITL23", "ITL5", "ITL6", "ITvm", "CTL6", "ETL5", "NPL5", "Pvalb", "Sst", "Vip"),
  dmPFC = c("ITL23", "ITL5", "ITL6", "CTL6", "ETL5", "NPL5", "Pvalb", "Sst")
)

jobs <- bind_rows(lapply(names(cluster_sets), function(region) {
  tibble(
    region = region,
    cluster = cluster_sets[[region]],
    pb_rds = file.path(pb_dir, paste0("main_pb_", region, ".rds")),
    meta_rds = file.path(pb_dir, paste0("main_meta_", region, ".rds")),
    cache_root = cache_root
  )
}))

if (nrow(jobs) != 18L) stop("Expected 18 cache jobs; generated ", nrow(jobs))
if (any(jobs$region == "dmPFC" & jobs$cluster == "ITvm")) stop("ITvm cannot be included for dmPFC")

dir.create(dirname(jobs_file), recursive = TRUE, showWarnings = FALSE)
write_tsv(jobs, jobs_file)

message("Wrote ", nrow(jobs), " cache jobs: ", jobs_file)
print(jobs %>% select(region, cluster, cache_root))
