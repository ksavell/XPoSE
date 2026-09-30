#!/usr/bin/env Rscript

# ==============================================================================
# XPoSE-seq experience decoder: validate nested-cache inputs
# ==============================================================================
# Validates the 18 region x cell-type jobs before expensive DESeq2 cache jobs are
# submitted. No DESeq2 models are fit here.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
jobs_file <- if (length(args) >= 1L) args[[1]] else "05_experience_decoder/cache_jobs.tsv"
if (!file.exists(jobs_file)) stop("Jobs file not found: ", jobs_file)

jobs <- read_tsv(jobs_file, show_col_types = FALSE)
required_job_cols <- c("region", "cluster", "pb_rds", "meta_rds", "cache_root")
missing_job_cols <- setdiff(required_job_cols, names(jobs))
if (length(missing_job_cols) > 0) {
  stop("Jobs file is missing columns: ", paste(missing_job_cols, collapse = ", "))
}
if (nrow(jobs) != 18L) stop("Expected exactly 18 jobs; found ", nrow(jobs))
if (any(jobs$region == "dmPFC" & jobs$cluster == "ITvm")) stop("ITvm cannot be included for dmPFC")
if (sum(jobs$region == "vmPFC") != 10L || sum(jobs$region == "dmPFC") != 8L) {
  stop("Expected 10 vmPFC jobs and 8 dmPFC jobs")
}
if (anyDuplicated(jobs[, c("region", "cluster")])) stop("Duplicate region x cluster rows in jobs file")

validate_one <- function(job) {
  region <- as.character(job$region)
  cluster <- as.character(job$cluster)
  pb_file <- as.character(job$pb_rds)
  meta_file <- as.character(job$meta_rds)

  if (!file.exists(pb_file)) stop(region, " / ", cluster, ": missing pseudobulk file: ", pb_file)
  if (!file.exists(meta_file)) stop(region, " / ", cluster, ": missing metadata file: ", meta_file)

  pb <- readRDS(pb_file)
  meta <- readRDS(meta_file)

  if (!is.list(pb) || !"RNA" %in% names(pb)) stop(region, " / ", cluster, ": pb_rds lacks $RNA")
  counts <- pb$RNA
  if (!is.data.frame(meta)) meta <- as.data.frame(meta)

  required_meta <- c("sample_id", "region", "population", "experience", "sex", "cluster_name", "ratID")
  missing_meta <- setdiff(required_meta, names(meta))
  if (length(missing_meta) > 0) {
    stop(region, " / ", cluster, ": metadata missing: ", paste(missing_meta, collapse = ", "))
  }

  meta[] <- lapply(meta, function(x) if (is.factor(x)) as.character(x) else x)
  sub <- meta %>%
    filter(
      .data$region == region,
      .data$cluster_name == cluster,
      .data$experience %in% c("RT", "NC"),
      .data$population %in% c("active", "non-active")
    )

  if (anyDuplicated(sub$sample_id)) stop(region, " / ", cluster, ": duplicate sample IDs")
  if (!all(sub$sample_id %in% colnames(counts))) stop(region, " / ", cluster, ": metadata/count columns are not aligned")
  if (nrow(counts) < 100L) stop(region, " / ", cluster, ": fewer than 100 genes in count matrix")

  rat_meta <- sub %>% distinct(ratID, experience, sex)
  if (nrow(rat_meta) != 12L) stop(region, " / ", cluster, ": expected 12 rats; found ", nrow(rat_meta))
  if (any(table(rat_meta$sex) != 6L)) stop(region, " / ", cluster, ": expected 6 rats per sex")

  sex_exp <- table(rat_meta$sex, rat_meta$experience)
  if (!all(c("RT", "NC") %in% colnames(sex_exp)) || any(sex_exp[, c("RT", "NC"), drop = FALSE] != 3L)) {
    stop(region, " / ", cluster, ": expected 3 RT and 3 NC rats within each sex")
  }

  completeness <- sub %>%
    count(ratID, population, name = "n") %>%
    complete(ratID = rat_meta$ratID, population = c("active", "non-active"), fill = list(n = 0L))
  if (any(completeness$n != 1L)) {
    stop(region, " / ", cluster, ": every rat must have exactly one Active and one Non-active sample")
  }

  sexes <- sort(unique(rat_meta$sex))
  if (length(sexes) != 2L) stop(region, " / ", cluster, ": expected exactly two sex strata")
  rats_1 <- rat_meta$ratID[rat_meta$sex == sexes[1]]
  rats_2 <- rat_meta$ratID[rat_meta$sex == sexes[2]]
  n_five_rat_subsets <- choose(length(rats_1), 2) * choose(length(rats_2), 3) +
    choose(length(rats_1), 3) * choose(length(rats_2), 2)
  n_same_sex_pairs <- choose(length(rats_1), 2) + choose(length(rats_2), 2)

  if (n_five_rat_subsets != 600L) stop(region, " / ", cluster, ": expected 600 five-rat subsets")
  if (n_same_sex_pairs != 30L) stop(region, " / ", cluster, ": expected 30 same-sex held-out pairs")

  nuclei_text <- if ("n_nuclei" %in% names(sub)) {
    paste0(" | nuclei range=", min(sub$n_nuclei), "-", max(sub$n_nuclei))
  } else {
    ""
  }

  message(
    "OK: ", region, " / ", cluster,
    " | rats=12 | samples=", nrow(sub),
    " | five-rat subsets=600 | held-out pairs=30",
    nuclei_text
  )

  invisible(TRUE)
}

for (i in seq_len(nrow(jobs))) validate_one(jobs[i, , drop = FALSE])
message("VALIDATION COMPLETE: all 18 cache jobs passed")
