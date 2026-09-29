#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
})

option_list <- list(
  make_option("--input_root", type = "character", default = NULL),
  make_option("--out_dir", type = "character", default = NULL)
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input_root) || is.null(opt$out_dir)) {
  stop("--input_root and --out_dir are required")
}

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

all_files <- list.files(
  opt$input_root,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

summary_files <- all_files[basename(all_files) == "decoder_region_summary.csv"]
exact_files <- all_files[basename(all_files) == "all_exact_pair_concordance.csv"]

if (length(summary_files) != 2L) {
  stop(
    "Expected exactly two decoder_region_summary.csv files (dmPFC and vmPFC); found ",
    length(summary_files)
  )
}

summary_table <- bind_rows(lapply(
  summary_files,
  function(file) read_csv(file, show_col_types = FALSE)
)) %>%
  arrange(match(region, c("dmPFC", "vmPFC")))

if (!setequal(summary_table$region, c("dmPFC", "vmPFC"))) {
  stop("Collected summaries must contain exactly dmPFC and vmPFC")
}
if (anyDuplicated(summary_table$region)) {
  stop("More than one summary row was found for a region")
}

write_csv(
  summary_table,
  file.path(opt$out_dir, "decoder_percent_correct.csv")
)

if (length(exact_files) == 2L) {
  exact_table <- bind_rows(lapply(
    exact_files,
    function(file) read_csv(file, show_col_types = FALSE)
  )) %>%
    arrange(region, perm_id)

  write_csv(
    exact_table,
    file.path(opt$out_dir, "all_exact_pair_concordance_by_region.csv")
  )
}

file.create(file.path(opt$out_dir, "_SUCCESS"))
