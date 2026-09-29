#!/usr/bin/env Rscript

# ==============================================================================
# XPoSE-seq experience decoder
# Nested regional RT + NC DEG-union activity-difference decoder
# ==============================================================================
# Final paper-facing analysis only. Feature selection remains nested inside each
# outer fold and all 400 exact within-sex RT/NC assignments.
# Primary metric: held-out RT/NC pair concordance.
# percent_correct = 100 * concordance.
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(glmnet)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(parallel)
})

option_list <- list(
  make_option("--region", type = "character", default = NULL),
  make_option("--clusters", type = "character", default = NULL),
  make_option("--cache_root", type = "character", default = NULL),
  make_option("--cache_run_id", type = "character", default = "paired_deg_cache_5pct"),
  make_option("--out_root", type = "character", default = NULL),
  make_option("--run_id", type = "character", default = "experience_decoder"),
  make_option("--top_frac", type = "double", default = 0.05),
  make_option("--case_label", type = "character", default = "RT"),
  make_option("--control_label", type = "character", default = "NC"),
  make_option("--lambda_choice", type = "character", default = "lambda.min"),
  make_option("--inner_max_folds", type = "integer", default = 5L),
  make_option("--min_train_per_class", type = "integer", default = 4L),
  make_option("--min_meta_clusters", type = "integer", default = 1L),
  make_option("--perm_chunk_size", type = "integer", default = 10L),
  make_option("--n_cores", type = "integer", default = NA_integer_),
  make_option("--resume", action = "store_true", default = TRUE),
  make_option("--no_resume", action = "store_false", dest = "resume")
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("region", "clusters", "cache_root", "out_root")
missing <- required[vapply(required, function(x) is.null(opt[[x]]) || !nzchar(trimws(opt[[x]])), logical(1))]
if (length(missing)) stop("Missing required option(s): ", paste(paste0("--", missing), collapse = ", "))
if (!opt$region %in% c("dmPFC", "vmPFC")) stop("--region must be dmPFC or vmPFC")
if (!opt$lambda_choice %in% c("lambda.min", "lambda.1se")) stop("--lambda_choice must be lambda.min or lambda.1se")
if (!is.finite(opt$top_frac) || opt$top_frac <= 0 || opt$top_frac >= 1) stop("--top_frac must be > 0 and < 1")
if (!is.finite(opt$perm_chunk_size) || opt$perm_chunk_size < 1) stop("--perm_chunk_size must be >= 1")

clusters <- unique(trimws(strsplit(opt$clusters, ",", fixed = TRUE)[[1]]))
clusters <- clusters[nzchar(clusters)]
if (length(clusters) < 2) stop("Hierarchical decoder requires at least 2 clusters")
if (opt$region == "dmPFC" && "ITvm" %in% clusters) stop("dmPFC must not include ITvm")

n_cores <- opt$n_cores
if (is.na(n_cores)) n_cores <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))
if (is.na(n_cores) || n_cores < 1) n_cores <- 1L
if (.Platform$OS.type == "windows") n_cores <- 1L

sanitize_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
frac_tag <- function(x) {
  pct <- x * 100
  txt <- if (abs(pct - round(pct)) < 1e-10) sprintf("%d", round(pct)) else sub("0+$", "", sub("\\.$", "", sprintf("%.3f", pct)))
  paste0("top_", txt, "pct")
}
subset_key <- function(rats) paste(sort(as.character(rats)), collapse = "|")

cfg <- list(
  region = opt$region,
  clusters = clusters,
  decoder = "activity_difference",
  case_label = opt$case_label,
  control_label = opt$control_label,
  lambda_choice = opt$lambda_choice,
  inner_max_folds = opt$inner_max_folds,
  min_train_per_class = opt$min_train_per_class,
  min_meta_clusters = opt$min_meta_clusters,
  top_frac = opt$top_frac,
  n_cores = n_cores,
  perm_chunk_size = opt$perm_chunk_size
)

out_dir <- file.path(opt$out_root, paste(sanitize_name(opt$run_id), cfg$region, frac_tag(cfg$top_frac), sep = "_"))
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "checkpoints"), recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Load nested DEG/expression caches
# ==============================================================================

cache_path <- function(cluster) file.path(
  opt$cache_root, sanitize_name(opt$cache_run_id), cfg$region,
  sanitize_name(cluster), frac_tag(cfg$top_frac),
  "nested_deg_and_expression_cache.rds"
)

cluster_cache <- setNames(lapply(cfg$clusters, function(cluster) {
  path <- cache_path(cluster)
  if (!file.exists(path)) stop("Missing cache for ", cluster, ": ", path)
  x <- readRDS(path)
  if (!identical(as.character(x$region), cfg$region)) stop("Cache region mismatch for ", cluster)
  if (!identical(as.character(x$cluster), cluster)) stop("Cache cluster mismatch for ", cluster)
  if (!length(x$expression_pair_cache)) stop("Expression pair cache is empty for ", cluster)
  if (!cfg$decoder %in% names(x$expression_pair_cache[[1]]$decoders)) stop("activity_difference representation absent for ", cluster)
  if (length(x$deg_by_subset) != 600L) stop("Expected 600 nested DEG subsets for ", cluster)
  x
}), cfg$clusters)

cohort_key <- function(x) sort(unique(paste(x$decoder_rat, x$decoder_experience, x$decoder_sex, sep = "|")))
rat_meta <- cluster_cache[[1]]$rat_meta %>% arrange(decoder_sex, decoder_experience, decoder_rat)
ref_key <- cohort_key(rat_meta)
for (cluster in cfg$clusters[-1]) {
  key <- cohort_key(cluster_cache[[cluster]]$rat_meta)
  if (length(key) != length(ref_key) || !setequal(key, ref_key)) stop("Rat cohort differs across cache clusters: ", cluster)
}
if (nrow(rat_meta) != 12L) stop("Expected 12 rats; found ", nrow(rat_meta))
if (any(table(rat_meta$decoder_sex) != 6L)) stop("Expected 6 rats per sex")

observed_labels <- setNames(as.character(rat_meta$decoder_experience), rat_meta$decoder_rat)

pair_ids <- names(cluster_cache[[1]]$expression_pair_cache)
if (length(pair_ids) != 30L) stop("Expected 30 same-sex candidate pairs; found ", length(pair_ids))
for (cluster in cfg$clusters[-1]) {
  if (!setequal(pair_ids, names(cluster_cache[[cluster]]$expression_pair_cache))) stop("Pair-cache IDs differ across clusters: ", cluster)
}

pair_defs <- setNames(lapply(pair_ids, function(fold_id) {
  x <- cluster_cache[[1]]$expression_pair_cache[[fold_id]]
  list(train_rats = x$train_rats, test_rats = x$test_rats, test_sexes = x$test_sexes)
}), pair_ids)

eligible_fold_ids <- function(labels) {
  keep <- vapply(pair_defs, function(x) {
    lab <- labels[x$test_rats]
    length(lab) == 2L && all(!is.na(lab)) && length(unique(lab)) == 2L
  }, logical(1))
  names(pair_defs)[keep]
}
if (length(eligible_fold_ids(observed_labels)) != 18L) stop("Observed assignment must yield 18 eligible folds")

# ==============================================================================
# Complete 400-assignment within-sex label space
# ==============================================================================

assignment_key <- function(labels, rat_order) paste(as.character(labels[rat_order]), collapse = "|")
build_assignments <- function() {
  rat_order <- rat_meta$decoder_rat
  strata <- split(rat_order, rat_meta$decoder_sex)
  case_counts <- vapply(strata, function(r) sum(observed_labels[r] == cfg$case_label), integer(1))
  combos <- lapply(seq_along(strata), function(i) combn(strata[[i]], case_counts[[i]], simplify = FALSE))
  grid <- expand.grid(lapply(combos, seq_along), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

  assignments <- lapply(seq_len(nrow(grid)), function(i) {
    lab <- setNames(rep(cfg$control_label, length(rat_order)), rat_order)
    for (j in seq_along(combos)) lab[combos[[j]][[grid[[j]][i]]]] <- cfg$case_label
    lab
  })

  keys <- vapply(assignments, assignment_key, character(1), rat_order = rat_order)
  assignments <- assignments[!duplicated(keys)]
  keys <- vapply(assignments, assignment_key, character(1), rat_order = rat_order)
  obs_idx <- match(assignment_key(observed_labels, rat_order), keys)
  if (is.na(obs_idx)) stop("Observed assignment not found in exact assignment space")
  if (obs_idx != 1L) assignments <- c(list(assignments[[obs_idx]]), assignments[-obs_idx])
  names(assignments) <- sprintf("P%04d", seq_along(assignments))
  if (length(assignments) != 400L) stop("Expected 400 exact assignments; found ", length(assignments))
  list(assignments = assignments, observed_id = names(assignments)[1])
}
perm_info <- build_assignments()

# ==============================================================================
# Training-fold regional RT + NC DEG union
# ==============================================================================

lookup_deg_genes <- function(cluster, rats) {
  item <- cluster_cache[[cluster]]$deg_by_subset[[subset_key(rats)]]
  if (is.null(item)) stop("Missing DEG cache entry for ", cluster, " subset ", subset_key(rats))
  sort(unique(as.character(item$deg_genes)))
}

build_region_union <- function(train_rats, labels, fold_id) {
  rt_rats <- train_rats[labels[train_rats] == cfg$case_label]
  nc_rats <- train_rats[labels[train_rats] == cfg$control_label]
  if (length(rt_rats) != 5L || length(nc_rats) != 5L) stop("Outer training fold must contain 5 RT and 5 NC rats")

  rt_by_cluster <- lapply(cfg$clusters, function(cl) lookup_deg_genes(cl, rt_rats))
  nc_by_cluster <- lapply(cfg$clusters, function(cl) lookup_deg_genes(cl, nc_rats))
  region_rt <- sort(unique(unlist(rt_by_cluster, use.names = FALSE)))
  region_nc <- sort(unique(unlist(nc_by_cluster, use.names = FALSE)))
  region_union <- union(region_rt, region_nc)

  available <- lapply(cfg$clusters, function(cl) {
    sort(unique(as.character(cluster_cache[[cl]]$expression_pair_cache[[fold_id]]$decoders[[cfg$decoder]]$ordered_genes)))
  })
  common_analyzable <- sort(Reduce(intersect, available))
  if (!length(common_analyzable)) stop("No common analyzable genes across clusters for fold ", fold_id)

  genes <- sort(intersect(region_union, common_analyzable))
  list(
    genes = genes,
    counts = tibble(
      fold_id = fold_id,
      n_rt_union_genes = length(region_rt),
      n_nc_union_genes = length(region_nc),
      n_raw_union_genes = length(region_union),
      n_common_analyzable_genes = length(common_analyzable),
      n_final_genes = length(genes),
      rt_training_rats = paste(sort(rt_rats), collapse = ","),
      nc_training_rats = paste(sort(nc_rats), collapse = ",")
    )
  )
}

# ==============================================================================
# Nested ridge classifiers
# ==============================================================================

make_inner_foldid <- function(labels, rat_ids) {
  y <- factor(labels, levels = c(cfg$control_label, cfg$case_label))
  tab <- table(y)
  if (length(tab) < 2 || min(tab) < cfg$min_train_per_class) return(NULL)
  k <- min(cfg$inner_max_folds, min(tab))
  if (k < 2) return(NULL)

  sex_map <- setNames(rat_meta$decoder_sex, rat_meta$decoder_rat)
  foldid <- integer(length(y))
  for (cl in levels(y)) {
    idx <- which(y == cl)
    idx <- idx[order(sex_map[rat_ids[idx]], rat_ids[idx])]
    foldid[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  foldid
}

balanced_weights <- function(y) {
  y <- droplevels(factor(y))
  tab <- table(y)
  length(y) / (2 * as.numeric(tab[as.character(y)]))
}

binomial_deviance <- function(y, prob, weights) {
  p <- pmin(pmax(prob, 1e-8), 1 - 1e-8)
  yy <- matrix(y, nrow = length(y), ncol = ncol(p))
  ww <- matrix(weights, nrow = length(weights), ncol = ncol(p))
  colSums(ww * (-2 * (yy * log(p) + (1 - yy) * log(1 - p)))) / colSums(ww)
}

fit_ridge <- function(x_train, y_train, train_rats, x_test, need_oof = FALSE) {
  x_train <- as.matrix(x_train)
  x_test <- as.matrix(x_test)
  if (ncol(x_train) < 1 || nrow(x_test) < 1) return(list(ok = FALSE))
  if (any(!is.finite(x_train)) || any(!is.finite(x_test))) return(list(ok = FALSE))

  y_factor <- factor(y_train, levels = c(cfg$control_label, cfg$case_label))
  if (length(table(y_factor)) < 2 || min(table(y_factor)) < cfg$min_train_per_class) return(list(ok = FALSE))
  y_numeric <- as.integer(y_factor == cfg$case_label)
  foldid <- make_inner_foldid(y_factor, train_rats)
  if (is.null(foldid)) return(list(ok = FALSE))

  path_fit <- tryCatch(suppressWarnings(glmnet(
    x_train, y_numeric, family = "binomial", alpha = 0,
    weights = balanced_weights(y_factor), standardize = TRUE, nlambda = 100
  )), error = function(e) NULL)
  if (is.null(path_fit) || length(path_fit$lambda) < 2) return(list(ok = FALSE))
  lambda_grid <- path_fit$lambda

  folds <- sort(unique(foldid))
  fold_dev <- matrix(NA_real_, nrow = length(folds), ncol = length(lambda_grid))
  oof <- matrix(NA_real_, nrow = nrow(x_train), ncol = length(lambda_grid))

  for (i in seq_along(folds)) {
    val <- which(foldid == folds[i])
    tr <- which(foldid != folds[i])
    ytr <- droplevels(y_factor[tr])
    yval <- droplevels(y_factor[val])
    if (length(unique(ytr)) < 2 || length(unique(yval)) < 2) return(list(ok = FALSE))

    fit_i <- tryCatch(suppressWarnings(glmnet(
      x_train[tr, , drop = FALSE], y_numeric[tr], family = "binomial", alpha = 0,
      lambda = lambda_grid, weights = balanced_weights(ytr), standardize = TRUE
    )), error = function(e) NULL)
    if (is.null(fit_i)) return(list(ok = FALSE))

    p_i <- tryCatch(as.matrix(predict(
      fit_i, newx = x_train[val, , drop = FALSE], s = lambda_grid, type = "response"
    )), error = function(e) NULL)
    if (is.null(p_i)) return(list(ok = FALSE))

    oof[val, ] <- p_i
    fold_dev[i, ] <- binomial_deviance(y_numeric[val], p_i, balanced_weights(yval))
  }

  mean_dev <- colMeans(fold_dev, na.rm = TRUE)
  se_dev <- apply(fold_dev, 2, sd, na.rm = TRUE) / sqrt(nrow(fold_dev))
  min_idx <- which.min(mean_dev)
  chosen <- if (cfg$lambda_choice == "lambda.min") {
    min_idx
  } else {
    candidates <- which(mean_dev <= mean_dev[min_idx] + se_dev[min_idx])
    candidates[which.max(lambda_grid[candidates])]
  }
  lambda_value <- lambda_grid[chosen]

  final_fit <- tryCatch(suppressWarnings(glmnet(
    x_train, y_numeric, family = "binomial", alpha = 0,
    lambda = lambda_grid, weights = balanced_weights(y_factor), standardize = TRUE
  )), error = function(e) NULL)
  if (is.null(final_fit)) return(list(ok = FALSE))

  test_prob <- tryCatch(as.numeric(predict(final_fit, newx = x_test, s = lambda_value, type = "response")), error = function(e) NULL)
  if (is.null(test_prob) || length(test_prob) != nrow(x_test) || any(!is.finite(test_prob))) return(list(ok = FALSE))

  oof_prob <- NULL
  if (need_oof) {
    oof_prob <- oof[, chosen]
    if (any(!is.finite(oof_prob))) return(list(ok = FALSE))
    names(oof_prob) <- train_rats
  }

  list(ok = TRUE, oof = oof_prob, test = setNames(test_prob, rownames(x_test)))
}

# ==============================================================================
# Held-out pair scoring
# ==============================================================================

pair_results <- function(predictions) {
  predictions %>%
    group_by(fold_id) %>%
    group_modify(~ {
      rt <- .x %>% filter(true_experience == cfg$case_label)
      nc <- .x %>% filter(true_experience == cfg$control_label)
      if (nrow(rt) != 1L || nrow(nc) != 1L) return(tibble(concordant = NA_real_))
      margin <- rt$probability_RT[[1]] - nc$probability_RT[[1]]
      tibble(
        rt_rat = rt$decoder_rat[[1]],
        nc_rat = nc$decoder_rat[[1]],
        sex = rt$decoder_sex[[1]],
        rt_probability = rt$probability_RT[[1]],
        nc_probability = nc$probability_RT[[1]],
        pair_margin = margin,
        concordant = ifelse(margin > 0, 1, ifelse(margin < 0, 0, 0.5))
      )
    }) %>%
    ungroup()
}

score_predictions <- function(predictions, perm_id) {
  pairs <- pair_results(predictions)
  if (nrow(pairs) != 18L || any(is.na(pairs$concordant))) stop("Expected 18 complete held-out pairs for ", perm_id)
  tibble(
    perm_id = perm_id,
    concordance = mean(pairs$concordant),
    percent_correct = 100 * mean(pairs$concordant),
    n_pairs = nrow(pairs),
    n_concordant = sum(pairs$concordant == 1),
    n_discordant = sum(pairs$concordant == 0),
    n_ties = sum(pairs$concordant == 0.5)
  )
}

# ==============================================================================
# Evaluate one label assignment
# ==============================================================================

evaluate_assignment <- function(labels, perm_id, detail = FALSE) {
  labels <- labels[rat_meta$decoder_rat]
  fold_ids <- eligible_fold_ids(labels)
  if (length(fold_ids) != 18L) stop("Assignment ", perm_id, " did not yield 18 eligible folds")

  pred_rows <- list()
  feature_rows <- list()

  for (fold_id in fold_ids) {
    fold <- pair_defs[[fold_id]]
    train_rats <- fold$train_rats
    test_rats <- fold$test_rats
    y_train <- unname(labels[train_rats])
    if (sum(y_train == cfg$case_label) != 5L || sum(y_train == cfg$control_label) != 5L) stop("Training fold is not 5 RT / 5 NC")

    nested <- build_region_union(train_rats, labels, fold_id)
    genes <- nested$genes
    feature_rows[[length(feature_rows) + 1L]] <- nested$counts %>% mutate(perm_id = perm_id, .before = 1)

    train_scores <- list()
    test_scores <- list()
    successful <- character()

    if (length(genes)) {
      for (cluster in cfg$clusters) {
        fm <- cluster_cache[[cluster]]$expression_pair_cache[[fold_id]]$decoders[[cfg$decoder]]
        if (!all(genes %in% colnames(fm$x_train_all)) || !all(genes %in% colnames(fm$x_test_all))) stop("Regional gene-panel invariant violated")

        fit1 <- fit_ridge(
          fm$x_train_all[, genes, drop = FALSE], y_train, train_rats,
          fm$x_test_all[, genes, drop = FALSE], need_oof = TRUE
        )
        if (!isTRUE(fit1$ok)) next
        successful <- c(successful, cluster)
        train_scores[[cluster]] <- fit1$oof[train_rats]
        test_scores[[cluster]] <- fit1$test[test_rats]
      }
    }

    successful <- cfg$clusters[cfg$clusters %in% unique(successful)]
    probabilities <- setNames(rep(0.5, length(test_rats)), test_rats)

    if (length(successful) >= cfg$min_meta_clusters) {
      xtr2 <- do.call(cbind, train_scores[successful])
      xte2 <- do.call(cbind, test_scores[successful])
      colnames(xtr2) <- successful
      colnames(xte2) <- successful
      rownames(xtr2) <- train_rats
      rownames(xte2) <- test_rats

      fit2 <- fit_ridge(xtr2, y_train, train_rats, xte2, need_oof = FALSE)
      if (isTRUE(fit2$ok)) probabilities <- fit2$test[test_rats]
    }

    pred_rows[[length(pred_rows) + 1L]] <- tibble(
      perm_id = perm_id,
      fold_id = fold_id,
      decoder_rat = test_rats,
      decoder_sex = unname(fold$test_sexes[test_rats]),
      true_experience = unname(labels[test_rats]),
      probability_RT = unname(probabilities[test_rats]),
      n_stage1_clusters = length(successful),
      n_gene_features = length(genes)
    )
  }

  predictions <- bind_rows(pred_rows)
  out <- list(score = score_predictions(predictions, perm_id), feature_counts = bind_rows(feature_rows))
  if (detail) {
    out$predictions <- predictions
    out$pairs <- pair_results(predictions)
  }
  out
}

# ==============================================================================
# Observed assignment
# ==============================================================================

observed_checkpoint <- file.path(out_dir, "checkpoints", "observed_assignment.rds")
if (isTRUE(opt$resume) && file.exists(observed_checkpoint)) {
  observed <- readRDS(observed_checkpoint)
} else {
  observed <- evaluate_assignment(perm_info$assignments[[1]], perm_info$observed_id, detail = TRUE)
  saveRDS(observed, observed_checkpoint)
}

write_csv(observed$predictions %>% mutate(region = cfg$region, .before = 1), file.path(out_dir, "tables", "observed_outer_predictions.csv"))
write_csv(observed$pairs %>% mutate(region = cfg$region, .before = 1), file.path(out_dir, "tables", "observed_pair_results.csv"))
write_csv(observed$feature_counts %>% mutate(region = cfg$region, .before = 1), file.path(out_dir, "tables", "observed_fold_feature_counts.csv"))

# ==============================================================================
# Exact 400-assignment inference
# ==============================================================================

indices <- seq_along(perm_info$assignments)
chunks <- split(indices, ceiling(indices / cfg$perm_chunk_size))
chunk_files <- file.path(out_dir, "checkpoints", sprintf("exact_chunk_%03d.rds", seq_along(chunks)))

for (ci in seq_along(chunks)) {
  idx <- chunks[[ci]]
  if (isTRUE(opt$resume) && file.exists(chunk_files[[ci]])) next

  worker <- function(i) tryCatch({
    result <- if (i == 1L) observed else evaluate_assignment(perm_info$assignments[[i]], names(perm_info$assignments)[i], detail = FALSE)
    list(ok = TRUE, score = result$score, error = NA_character_)
  }, error = function(e) list(ok = FALSE, score = NULL, error = conditionMessage(e)))

  ans <- if (n_cores > 1 && .Platform$OS.type != "windows") {
    mclapply(idx, worker, mc.cores = min(n_cores, length(idx)), mc.preschedule = FALSE)
  } else {
    lapply(idx, worker)
  }

  bad <- which(!vapply(ans, function(x) isTRUE(x$ok), logical(1)))
  if (length(bad)) {
    write_csv(
      tibble(chunk = ci, assignment_index = idx[bad], perm_id = names(perm_info$assignments)[idx[bad]], error = vapply(ans[bad], function(x) x$error, character(1))),
      file.path(out_dir, "tables", sprintf("exact_failures_chunk_%03d.csv", ci))
    )
    stop("Exact-permutation failure in chunk ", ci)
  }
  saveRDS(lapply(ans, `[[`, "score"), chunk_files[[ci]])
}

all_exact <- bind_rows(unlist(lapply(chunk_files, readRDS), recursive = FALSE)) %>%
  mutate(region = cfg$region, is_observed = perm_id == perm_info$observed_id, .before = 1)
if (nrow(all_exact) != 400L) stop("Expected 400 exact assignment results; found ", nrow(all_exact))
write_csv(all_exact, file.path(out_dir, "tables", "all_exact_pair_concordance.csv"))

obs <- all_exact %>% filter(is_observed)
if (nrow(obs) != 1L) stop("Expected exactly one observed assignment")
exact_p_upper <- mean(all_exact$concordance >= obs$concordance[[1]], na.rm = TRUE)

summary <- obs %>% transmute(
  region,
  decoder = cfg$decoder,
  feature_panel = "RT_NC_DEG_union",
  n_pairs,
  n_concordant,
  n_discordant,
  n_ties,
  percent_correct,
  concordance,
  exact_p_upper = exact_p_upper,
  n_exact_assignments = nrow(all_exact)
)
write_csv(summary, file.path(out_dir, "tables", "decoder_region_summary.csv"))
file.create(file.path(out_dir, "_SUCCESS"))
