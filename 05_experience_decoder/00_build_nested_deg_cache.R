#!/usr/bin/env Rscript

# ==============================================================================
# XPoSE-seq experience decoder: build nested DEG/expression cache
# ==============================================================================
# One job builds one region x cell-type cache used by the final experience
# decoder. The cache preserves the original leakage-free nested design:
#   - every possible sex-balanced five-rat training subset is enumerated (600);
#   - Active vs Non-active DESeq2 is rebuilt from raw pseudobulk counts for each
#     subset with the paired design ~ decoder_cond + decoder_rat;
#   - the 30 possible same-sex held-out pairs receive training-only expression
#     preprocessing for the Activity-difference representation.
#
# Run from the repository root. The input pseudobulk/meta objects are generated
# by 03_differential_expression/03_main_de/01_create_pseudobulk_meta.R.
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(matrixStats)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(parallel)
})

option_list <- list(
  make_option("--region", type = "character"),
  make_option("--cluster", type = "character"),
  make_option("--pb_rds", type = "character"),
  make_option("--meta_rds", type = "character"),
  make_option("--cache_root", type = "character"),
  make_option("--chunk_size", type = "integer", default = 20L),
  make_option("--n_cores", type = "integer", default = NA_integer_),
  make_option("--rebuild", action = "store_true", default = FALSE)
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("region", "cluster", "pb_rds", "meta_rds", "cache_root")
missing <- required[vapply(required, function(x) {
  is.null(opt[[x]]) || !nzchar(trimws(opt[[x]]))
}, logical(1))]
if (length(missing) > 0) {
  stop("Missing required option(s): ", paste(paste0("--", missing), collapse = ", "))
}
if (!opt$region %in% c("dmPFC", "vmPFC")) stop("--region must be dmPFC or vmPFC")
if (opt$region == "dmPFC" && opt$cluster == "ITvm") stop("ITvm is vmPFC-specific")
if (!file.exists(opt$pb_rds)) stop("Pseudobulk file not found: ", opt$pb_rds)
if (!file.exists(opt$meta_rds)) stop("Metadata file not found: ", opt$meta_rds)
if (opt$chunk_size < 1) stop("--chunk_size must be >= 1")

# Final-analysis constants ------------------------------------------------------
run_id <- "paired_deg_cache_5pct"
top_frac <- 0.05
active_label <- "active"
nonactive_label <- "non-active"
case_label <- "RT"
control_label <- "NC"
deg_padj_thresh <- 0.05
deg_lfc_thresh <- 0
deg_basemean_min <- 0
expression_basemean_min <- 5
detect_frac <- 0.50
detect_min_count <- 0
adjust_sex <- TRUE

n_cores <- opt$n_cores
if (is.na(n_cores)) {
  n_cores <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
}
if (is.na(n_cores) || n_cores < 1) n_cores <- 1L
if (.Platform$OS.type == "windows") n_cores <- 1L

sanitize_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
subset_key <- function(rats) paste(sort(as.character(rats)), collapse = "|")
pair_id <- function(rats) paste0("PAIR__", paste(sort(as.character(rats)), collapse = "__"))

cache_dir <- file.path(
  opt$cache_root,
  run_id,
  opt$region,
  sanitize_name(opt$cluster),
  "top_5pct"
)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cache_dir, "deg_chunks"), recursive = TRUE, showWarnings = FALSE)

final_cache_file <- file.path(cache_dir, "nested_deg_and_expression_cache.rds")
if (file.exists(final_cache_file) && !isTRUE(opt$rebuild)) {
  message("Final cache already exists; skipping: ", final_cache_file)
  quit(save = "no", status = 0)
}

# Load raw pseudobulk counts + metadata ----------------------------------------
pb <- readRDS(opt$pb_rds)
meta <- readRDS(opt$meta_rds)

if (!is.list(pb) || !"RNA" %in% names(pb)) stop("pb_rds must contain list element $RNA")
counts_all_input <- pb$RNA
if (!is.matrix(counts_all_input) && !inherits(counts_all_input, "Matrix")) {
  stop("pb$RNA must be a matrix-like genes x pseudobulk-samples object")
}
if (!is.data.frame(meta)) meta <- as.data.frame(meta)

required_meta <- c("sample_id", "region", "population", "experience", "sex", "cluster_name", "ratID")
missing_meta <- setdiff(required_meta, names(meta))
if (length(missing_meta) > 0) {
  stop("Metadata is missing required columns: ", paste(missing_meta, collapse = ", "))
}

meta[] <- lapply(meta, function(x) if (is.factor(x)) as.character(x) else x)
meta_use <- meta %>%
  filter(
    region == opt$region,
    cluster_name == opt$cluster,
    experience %in% c(case_label, control_label),
    population %in% c(active_label, nonactive_label)
  ) %>%
  mutate(
    decoder_experience = as.character(experience),
    decoder_population = as.character(population),
    decoder_rat = as.character(ratID),
    decoder_sex = as.character(sex)
  )

if (nrow(meta_use) == 0) stop("No RT/NC Active/Non-active samples found for ", opt$region, " / ", opt$cluster)
if (anyDuplicated(meta_use$sample_id)) stop("Duplicate sample_id values in filtered metadata")
if (!all(meta_use$sample_id %in% colnames(counts_all_input))) {
  stop("Some filtered metadata sample IDs are absent from the pseudobulk count matrix")
}

# Preserve metadata order and align counts exactly.
counts_all <- counts_all_input[, meta_use$sample_id, drop = FALSE]
meta_all <- meta_use[match(colnames(counts_all), meta_use$sample_id), , drop = FALSE]
stopifnot(identical(colnames(counts_all), meta_all$sample_id))

# Cohort validation -------------------------------------------------------------
rat_meta <- meta_all %>%
  distinct(decoder_rat, decoder_experience, decoder_sex) %>%
  arrange(decoder_sex, decoder_experience, decoder_rat)

if (nrow(rat_meta) != 12L) stop("Expected 12 RT/NC rats; found ", nrow(rat_meta))
if (any(table(rat_meta$decoder_sex) != 6L)) {
  stop("Expected exactly 6 rats per sex; observed: ", paste(capture.output(print(table(rat_meta$decoder_sex))), collapse = " "))
}
sex_exp <- table(rat_meta$decoder_sex, rat_meta$decoder_experience)
if (!all(c(case_label, control_label) %in% colnames(sex_exp)) || any(sex_exp[, c(case_label, control_label), drop = FALSE] != 3L)) {
  stop("Expected 3 RT and 3 NC rats within each sex; observed: ", paste(capture.output(print(sex_exp)), collapse = " | "))
}

completeness <- meta_all %>%
  count(decoder_rat, decoder_population, name = "n") %>%
  complete(
    decoder_rat = rat_meta$decoder_rat,
    decoder_population = c(active_label, nonactive_label),
    fill = list(n = 0L)
  )
if (any(completeness$n != 1L)) {
  stop("Each rat must have exactly one Active and one Non-active pseudobulk sample")
}
if (nrow(counts_all) < 100L) stop("Fewer than 100 genes are available")

# Every five-rat class subset possible under the balanced paired design ---------
sex_levels <- sort(unique(rat_meta$decoder_sex))
if (length(sex_levels) != 2L) stop("Expected exactly two sex strata")

rats_s1 <- sort(rat_meta$decoder_rat[rat_meta$decoder_sex == sex_levels[1]])
rats_s2 <- sort(rat_meta$decoder_rat[rat_meta$decoder_sex == sex_levels[2]])

subsets_2_3 <- unlist(
  lapply(combn(rats_s1, 2, simplify = FALSE), function(a) {
    lapply(combn(rats_s2, 3, simplify = FALSE), function(b) c(a, b))
  }),
  recursive = FALSE
)
subsets_3_2 <- unlist(
  lapply(combn(rats_s1, 3, simplify = FALSE), function(a) {
    lapply(combn(rats_s2, 2, simplify = FALSE), function(b) c(a, b))
  }),
  recursive = FALSE
)
training_subsets <- c(subsets_2_3, subsets_3_2)
keys <- vapply(training_subsets, subset_key, character(1))
training_subsets <- training_subsets[!duplicated(keys)]
names(training_subsets) <- vapply(training_subsets, subset_key, character(1))
if (length(training_subsets) != 600L) {
  stop("Expected 600 unique five-rat subsets; found ", length(training_subsets))
}
message("Five-rat DESeq2 subsets: ", length(training_subsets))

# Nested Active vs Non-active DESeq2 -------------------------------------------
run_subset_deseq <- function(rats) {
  sub_meta <- meta_all %>%
    filter(
      decoder_rat %in% rats,
      decoder_population %in% c(active_label, nonactive_label)
    ) %>%
    arrange(
      match(decoder_rat, sort(rats)),
      match(decoder_population, c(nonactive_label, active_label))
    )

  if (nrow(sub_meta) != 10L) {
    stop("Expected 10 samples for subset ", subset_key(rats), "; found ", nrow(sub_meta))
  }

  sub_counts <- counts_all[, sub_meta$sample_id, drop = FALSE]
  coldata <- data.frame(
    decoder_cond = factor(sub_meta$decoder_population, levels = c(nonactive_label, active_label)),
    decoder_rat = factor(sub_meta$decoder_rat),
    row.names = sub_meta$sample_id
  )

  dds <- DESeqDataSetFromMatrix(
    countData = round(sub_counts),
    colData = coldata,
    design = ~ decoder_cond + decoder_rat
  )

  fit <- NULL
  fit_used <- NA_character_
  errors <- character()
  for (fit_type in c("parametric", "local", "mean")) {
    attempt <- tryCatch(
      suppressWarnings(DESeq(dds, fitType = fit_type, quiet = TRUE)),
      error = function(e) e
    )
    if (!inherits(attempt, "error")) {
      fit <- attempt
      fit_used <- fit_type
      break
    }
    errors <- c(errors, paste0(fit_type, ": ", conditionMessage(attempt)))
  }
  if (is.null(fit)) {
    stop("DESeq2 failed for subset ", subset_key(rats), " | ", paste(errors, collapse = " || "))
  }

  res <- results(
    fit,
    contrast = c("decoder_cond", active_label, nonactive_label),
    alpha = deg_padj_thresh,
    cooksCutoff = FALSE,
    independentFiltering = FALSE
  )
  res_df <- as.data.frame(res)

  deg_genes <- rownames(res_df)[
    !is.na(res_df$padj) &
      !is.na(res_df$log2FoldChange) &
      !is.na(res_df$baseMean) &
      res_df$baseMean >= deg_basemean_min &
      res_df$padj < deg_padj_thresh &
      abs(res_df$log2FoldChange) > deg_lfc_thresh
  ]

  list(
    subset_key = subset_key(rats),
    rats = sort(rats),
    deg_genes = sort(unique(deg_genes)),
    n_deg = length(unique(deg_genes)),
    n_tested = nrow(res_df),
    n_padj_nonmissing = sum(!is.na(res_df$padj)),
    fit_type = fit_used
  )
}

chunk_indices <- split(
  seq_along(training_subsets),
  ceiling(seq_along(training_subsets) / opt$chunk_size)
)
chunk_files <- file.path(
  cache_dir,
  "deg_chunks",
  sprintf("deg_chunk_%03d.rds", seq_along(chunk_indices))
)

for (chunk_i in seq_along(chunk_indices)) {
  chunk_file <- chunk_files[chunk_i]
  idx <- chunk_indices[[chunk_i]]

  if (file.exists(chunk_file) && !isTRUE(opt$rebuild)) {
    message("Reusing DEG chunk ", chunk_i, "/", length(chunk_indices))
    next
  }

  message(
    "Running DEG chunk ", chunk_i, "/", length(chunk_indices),
    " | subsets ", min(idx), "-", max(idx)
  )
  chunk_subsets <- training_subsets[idx]

  if (n_cores > 1L && .Platform$OS.type != "windows") {
    answer <- mclapply(
      chunk_subsets,
      run_subset_deseq,
      mc.cores = min(n_cores, length(chunk_subsets)),
      mc.preschedule = FALSE
    )
  } else {
    answer <- lapply(chunk_subsets, run_subset_deseq)
  }

  names(answer) <- names(chunk_subsets)
  saveRDS(answer, chunk_file, compress = TRUE)
}

deg_by_subset <- do.call(c, lapply(chunk_files, readRDS))
deg_by_subset <- deg_by_subset[!duplicated(names(deg_by_subset))]
if (length(deg_by_subset) != 600L) {
  stop("Combined DEG cache expected 600 subsets; found ", length(deg_by_subset))
}

write_csv(
  bind_rows(lapply(deg_by_subset, function(x) {
    tibble(
      subset_key = x$subset_key,
      rats = paste(x$rats, collapse = ","),
      n_deg = x$n_deg,
      n_tested = x$n_tested,
      n_padj_nonmissing = x$n_padj_nonmissing,
      fit_type = x$fit_type
    )
  })),
  file.path(cache_dir, "deg_subset_summary.csv")
)

# Training-only Activity-difference expression cache ----------------------------
estimate_training_reference_size_factors <- function(train_counts, all_counts) {
  n_train <- ncol(train_counts)
  strict_valid <- rowSums(train_counts > 0) == n_train

  if (sum(strict_valid) >= 100L) {
    geo_means <- exp(rowMeans(log(train_counts[strict_valid, , drop = FALSE])))
    valid_genes <- rownames(train_counts)[strict_valid]
  } else {
    log_train <- log(train_counts)
    log_train[train_counts <= 0] <- NA_real_
    geo_all <- exp(rowMeans(log_train, na.rm = TRUE))
    positive_n <- rowSums(train_counts > 0)
    valid <- is.finite(geo_all) & geo_all > 0 & positive_n >= max(2L, ceiling(n_train / 2))
    geo_means <- geo_all[valid]
    valid_genes <- rownames(train_counts)[valid]
  }

  train_lib <- colSums(train_counts)
  all_lib <- colSums(all_counts)
  library_reference <- exp(mean(log(train_lib[train_lib > 0])))
  library_sf <- all_lib / library_reference

  if (length(valid_genes) < 20L) {
    sf <- library_sf
    sf[!is.finite(sf) | sf <= 0] <- 1
    names(sf) <- colnames(all_counts)
    return(sf)
  }

  ratio_counts <- all_counts[valid_genes, , drop = FALSE]
  sf <- vapply(seq_len(ncol(ratio_counts)), function(j) {
    x <- ratio_counts[, j]
    valid <- is.finite(x) & x > 0 & is.finite(geo_means) & geo_means > 0
    if (sum(valid) < 20L) return(NA_real_)
    exp(median(log(x[valid]) - log(geo_means[valid]), na.rm = TRUE))
  }, numeric(1))

  train_positions <- match(colnames(train_counts), colnames(all_counts))
  train_sf <- sf[train_positions]

  if (any(!is.finite(train_sf) | train_sf <= 0)) {
    sf <- library_sf
  } else {
    sf <- sf / exp(mean(log(train_sf)))
    bad <- !is.finite(sf) | sf <= 0
    if (any(bad)) sf[bad] <- library_sf[bad]
  }

  sf[!is.finite(sf) | sf <= 0] <- 1
  names(sf) <- colnames(all_counts)
  sf
}

residualize_training_sex <- function(train_mat, test_mat, sex_train, sex_test) {
  if (!isTRUE(adjust_sex)) {
    return(list(train = train_mat, test = test_mat, adjusted = FALSE))
  }

  sex_train <- as.character(sex_train)
  sex_test <- as.character(sex_test)
  sex_levels_train <- sort(unique(sex_train))

  if (length(sex_levels_train) < 2L) {
    return(list(train = train_mat, test = test_mat, adjusted = FALSE))
  }
  if (any(!sex_test %in% sex_levels_train)) stop("Held-out sex absent from training fold")

  design_train <- model.matrix(
    ~ sex,
    data = data.frame(sex = factor(sex_train, levels = sex_levels_train))
  )
  design_test <- model.matrix(
    ~ sex,
    data = data.frame(sex = factor(sex_test, levels = sex_levels_train))
  )
  design_test <- design_test[, colnames(design_train), drop = FALSE]

  coefs <- qr.coef(qr(design_train), train_mat)
  coefs[!is.finite(coefs)] <- 0
  nuisance <- setdiff(seq_len(ncol(design_train)), 1L)

  if (!length(nuisance)) {
    return(list(train = train_mat, test = test_mat, adjusted = FALSE))
  }

  list(
    train = train_mat - design_train[, nuisance, drop = FALSE] %*% coefs[nuisance, , drop = FALSE],
    test = test_mat - design_test[, nuisance, drop = FALSE] %*% coefs[nuisance, , drop = FALSE],
    adjusted = TRUE
  )
}

finalize_matrix <- function(train_mat, test_mat, adjusted_flag) {
  vars <- matrixStats::colVars(train_mat)
  names(vars) <- colnames(train_mat)
  ordered_genes <- names(vars)[is.finite(vars) & vars > 0]
  ordered_genes <- ordered_genes[order(vars[ordered_genes], decreasing = TRUE)]

  if (length(ordered_genes) < 1L) stop("No positive-variance genes")

  list(
    x_train_all = train_mat[, ordered_genes, drop = FALSE],
    x_test_all = test_mat[, ordered_genes, drop = FALSE],
    ordered_genes = ordered_genes,
    n_positive_variance = length(ordered_genes),
    sex_adjusted = adjusted_flag,
    normalization = "training-reference median-ratio log2 expression; Active minus Non-active within rat"
  )
}

sample_key <- meta_all %>%
  select(sample_id, decoder_rat, decoder_experience, decoder_sex, decoder_population)

prepare_activity_difference <- function(train_rats, test_rats) {
  ordered_rats <- c(train_rats, test_rats)

  samples <- sample_key %>%
    filter(
      decoder_rat %in% ordered_rats,
      decoder_population %in% c(active_label, nonactive_label)
    ) %>%
    mutate(
      rat_order = match(decoder_rat, ordered_rats),
      population_order = match(decoder_population, c(active_label, nonactive_label))
    ) %>%
    arrange(rat_order, population_order)

  all_counts <- counts_all[, samples$sample_id, drop = FALSE]
  train_mask <- samples$decoder_rat %in% train_rats
  train_counts <- all_counts[, train_mask, drop = FALSE]

  size_factors <- estimate_training_reference_size_factors(train_counts, all_counts)
  normalized <- sweep(all_counts, 2, size_factors, "/")

  keep <- rowMeans(normalized[, train_mask, drop = FALSE]) >= expression_basemean_min &
    rowMeans(train_counts > detect_min_count) > detect_frac
  if (sum(keep) < 10L) stop("Fewer than 10 genes pass Activity-difference expression filtering")

  expr <- t(log2(normalized[keep, , drop = FALSE] + 1))
  rownames(expr) <- samples$sample_id

  make_delta <- function(rat) {
    active_id <- samples$sample_id[samples$decoder_rat == rat & samples$decoder_population == active_label]
    nonactive_id <- samples$sample_id[samples$decoder_rat == rat & samples$decoder_population == nonactive_label]
    if (length(active_id) != 1L || length(nonactive_id) != 1L) {
      stop("Incomplete Active/Non-active pair for ", rat)
    }
    expr[active_id, , drop = FALSE] - expr[nonactive_id, , drop = FALSE]
  }

  delta <- do.call(rbind, lapply(ordered_rats, make_delta))
  rownames(delta) <- ordered_rats

  adjusted <- residualize_training_sex(
    delta[train_rats, , drop = FALSE],
    delta[test_rats, , drop = FALSE],
    rat_meta$decoder_sex[match(train_rats, rat_meta$decoder_rat)],
    rat_meta$decoder_sex[match(test_rats, rat_meta$decoder_rat)]
  )

  finalize_matrix(adjusted$train, adjusted$test, adjusted$adjusted)
}

candidate_pairs <- list()
for (sex_now in sex_levels) {
  rats_now <- sort(rat_meta$decoder_rat[rat_meta$decoder_sex == sex_now])
  candidate_pairs <- c(candidate_pairs, combn(rats_now, 2, simplify = FALSE))
}
if (length(candidate_pairs) != 30L) stop("Expected 30 same-sex held-out pairs")

expression_pair_cache <- setNames(
  vector("list", length(candidate_pairs)),
  vapply(candidate_pairs, pair_id, character(1))
)

for (i in seq_along(candidate_pairs)) {
  test_rats <- candidate_pairs[[i]]
  train_rats <- setdiff(rat_meta$decoder_rat, test_rats)
  fold_id <- pair_id(test_rats)

  message("Expression fold ", i, "/", length(candidate_pairs), " | ", fold_id)

  expression_pair_cache[[fold_id]] <- list(
    fold_id = fold_id,
    train_rats = train_rats,
    test_rats = test_rats,
    test_sexes = setNames(
      rat_meta$decoder_sex[match(test_rats, rat_meta$decoder_rat)],
      test_rats
    ),
    decoders = list(
      activity_difference = prepare_activity_difference(train_rats, test_rats)
    )
  )
}

# Final cache ------------------------------------------------------------------
cache_object <- list(
  version = "nested_deg_cache_main_pb_v1",
  region = opt$region,
  cluster = opt$cluster,
  run_id = run_id,
  input_pb_rds = normalizePath(opt$pb_rds),
  input_meta_rds = normalizePath(opt$meta_rds),
  rat_meta = rat_meta,
  deg_config = list(
    active_label = active_label,
    nonactive_label = nonactive_label,
    padj_thresh = deg_padj_thresh,
    lfc_thresh = deg_lfc_thresh,
    basemean_min = deg_basemean_min,
    design = "~ decoder_cond + decoder_rat",
    cooksCutoff = FALSE,
    independentFiltering = FALSE
  ),
  decoder_preprocess_config = list(
    cache_tag = "top_5pct",
    top_frac_legacy = top_frac,
    basemean_min = expression_basemean_min,
    detect_frac = detect_frac,
    detect_min_count = detect_min_count,
    adjust_sex = adjust_sex,
    representation = "activity_difference"
  ),
  deg_by_subset = deg_by_subset,
  expression_pair_cache = expression_pair_cache
)

saveRDS(cache_object, final_cache_file, compress = TRUE)

write_csv(
  tibble(
    region = opt$region,
    cluster = opt$cluster,
    n_rats = nrow(rat_meta),
    n_genes = nrow(counts_all),
    n_training_subsets = length(deg_by_subset),
    n_pair_candidates = length(expression_pair_cache),
    deg_padj_thresh = deg_padj_thresh,
    deg_lfc_thresh = deg_lfc_thresh,
    deg_design = "~ decoder_cond + decoder_rat",
    expression_basemean_min = expression_basemean_min,
    detect_frac = detect_frac,
    sex_adjusted = adjust_sex
  ),
  file.path(cache_dir, "cache_manifest.csv")
)

message("CACHE COMPLETE: ", final_cache_file)
