#!/usr/bin/env Rscript
# XPoSE-seq: POC active-fraction downsampling and differential expression
#
# One Biowulf job runs one cell type x target active percentage across all
# iterations, compares NC versus HC with pseudobulk DESeq2, and records the
# number of significant genes per iteration.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(ggplot2)
  library(readr)
})

# Command-line options ----------------------------------------------------
opt <- parse_args(OptionParser(option_list = list(
  make_option("--cluster", type = "character"),
  make_option("--percentage", type = "double",
              help = "Target active fraction of the NC pool, from 0 to 100"),
  make_option("--data_root", type = "character"),
  make_option("--out_root", type = "character"),
  make_option("--de_script", type = "character"),
  make_option("--run_id", type = "character"),
  make_option("--iterations", type = "integer", default = 100L),
  make_option("--alpha", type = "double", default = 0.05),
  make_option("--min_cell", type = "integer", default = 1L),
  make_option("--seed", type = "integer", default = 22L),
  make_option("--save_raw", action = "store_true", default = FALSE)
)))

required_args <- c("cluster", "percentage", "data_root", "out_root", "de_script", "run_id")
for (arg in required_args) {
  if (is.null(opt[[arg]]) || (is.character(opt[[arg]]) && !nzchar(opt[[arg]]))) {
    stop("Missing required option --", arg)
  }
}

source(opt$de_script)

# Ensure neurorestore/Libra is installed rather than the unrelated CRAN package.
suppressPackageStartupMessages(library(Libra))
if (!exists("to_pseudobulk", where = asNamespace("Libra"), inherits = FALSE)) {
  stop(
    "Wrong Libra package loaded: to_pseudobulk() is unavailable. ",
    "Install neurorestore/Libra."
  )
}

set.seed(opt$seed)

cl <- opt$cluster
pct <- opt$percentage

out_dir <- file.path(opt$out_root, paste0("run_", opt$run_id))
tally_dir <- file.path(out_dir, "tally")
summary_dir <- file.path(out_dir, "summary")
plot_dir <- file.path(out_dir, "plots")
raw_dir <- file.path(out_dir, "raw", cl)

for (d in c(tally_dir, summary_dir, plot_dir, if (opt$save_raw) raw_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# Downsample one cell type ------------------------------------------------
# For each NC rat, subsample active and non-active nuclei to the requested
# active fraction. The HC pool is then downsampled to the same total number of
# nuclei as the resulting NC pool.
group_downsample_percentage <- function(cl_obj, target_percentage, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  p <- target_percentage / 100
  md <- cl_obj@meta.data
  md$cell <- rownames(md)

  nc_rats <- unique(md$ratID[md$experience == "NC"])
  nc_cells <- character(0)

  for (rat in nc_rats) {
    active <- md$cell[md$ratID == rat & md$group == "active"]
    non_active <- md$cell[md$ratID == rat & md$group == "non-active"]

    n_active <- length(active)
    n_non_active <- length(non_active)

    if (p == 0) {
      keep_active <- 0L
      keep_non_active <- n_non_active
    } else if (p == 1) {
      keep_active <- n_active
      keep_non_active <- 0L
    } else {
      max_total <- floor(min(n_active / p, n_non_active / (1 - p)))
      keep_active <- round(p * max_total)
      keep_non_active <- max_total - keep_active
    }

    if (keep_active > 0) {
      nc_cells <- c(nc_cells, sample(active, keep_active))
    }
    if (keep_non_active > 0) {
      nc_cells <- c(nc_cells, sample(non_active, keep_non_active))
    }
  }

  hc_pool <- md$cell[md$experience == "HC"]
  n_hc <- min(length(nc_cells), length(hc_pool))
  hc_cells <- if (n_hc > 0) sample(hc_pool, n_hc) else character(0)

  c(nc_cells, hc_cells)
}

# Run one DE iteration ----------------------------------------------------
run_one_iteration <- function(obj, cluster_name) {
  suppressWarnings(suppressMessages(tryCatch({
    if (inherits(obj[["RNA"]], "Assay5")) {
      obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
    }

    invisible(capture.output(
      de_result <- single_factor_DESeq(
        obj,
        comp_vect = c("experience", "NC", "HC"),
        cluster = cluster_name,
        min_cell = opt$min_cell
      )
    ))

    list(result = de_result, error = NA_character_)
  }, error = function(e) {
    list(result = NULL, error = conditionMessage(e))
  })))
}

# Load one cell-type object ----------------------------------------------
cluster_file <- file.path(opt$data_root, paste0(cl, ".rds"))
if (!file.exists(cluster_file)) {
  stop("Cluster file not found: ", cluster_file)
}

cl_obj <- readRDS(cluster_file)

# Iterative downsampling --------------------------------------------------
n_iterations <- opt$iterations
all_indices <- vector("list", n_iterations)
all_seeds <- integer(n_iterations)
n_sig <- rep(NA_integer_, n_iterations)
error_message <- rep(NA_character_, n_iterations)
results <- vector("list", n_iterations)

for (iteration in seq_len(n_iterations)) {
  iteration_seed <- sample.int(10000, 1)
  all_seeds[iteration] <- iteration_seed

  chosen_cells <- group_downsample_percentage(
    cl_obj,
    pct,
    seed = iteration_seed
  )
  all_indices[[iteration]] <- chosen_cells

  sampled_obj <- subset(cl_obj, cells = chosen_cells)
  de_run <- run_one_iteration(sampled_obj, cl)

  results[[iteration]] <- de_run$result
  error_message[iteration] <- de_run$error

  if (!is.null(de_run$result)) {
    n_sig[iteration] <- sum(
      de_run$result$results$padj < opt$alpha,
      na.rm = TRUE
    )
  }
}

# Per-job outputs ---------------------------------------------------------
tag <- paste0(cl, "_", pct)

write_csv(
  data.frame(
    cluster = cl,
    percentage = pct,
    iteration = seq_len(n_iterations),
    n_sig = n_sig,
    error = error_message
  ),
  file.path(tally_dir, paste0("tally_", tag, ".csv"))
)

n_fail <- sum(is.na(n_sig))

write_csv(
  data.frame(
    cluster = cl,
    percentage = pct,
    mean_DE = mean(n_sig, na.rm = TRUE),
    sd_DE = if (sum(!is.na(n_sig)) > 1) sd(n_sig, na.rm = TRUE) else NA_real_,
    n_ok = sum(!is.na(n_sig)),
    n_fail = n_fail
  ),
  file.path(summary_dir, paste0("summary_", tag, ".csv"))
)

if (opt$save_raw) {
  saveRDS(
    list(indices = all_indices, seeds = all_seeds),
    file.path(raw_dir, paste0("indices_seeds_", pct, ".rds"))
  )
  saveRDS(
    results,
    file.path(raw_dir, paste0("results_", pct, ".rds"))
  )
}

if (sum(!is.na(n_sig)) > 0) {
  p <- ggplot(data.frame(n_sig = n_sig[!is.na(n_sig)]), aes(n_sig)) +
    geom_histogram(bins = 30, fill = "grey75", color = "white") +
    labs(
      title = paste0(cl, " | ", pct, "% active"),
      x = "DE genes per iteration (adjusted P < alpha)",
      y = "Count"
    ) +
    theme_classic(base_size = 12)

  ggsave(
    file.path(plot_dir, paste0("hist_", tag, ".png")),
    p,
    width = 6,
    height = 4,
    dpi = 300
  )
}
