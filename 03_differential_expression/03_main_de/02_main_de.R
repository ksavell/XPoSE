# XPoSE-seq: differential expression for the main dataset
#
# Run from the base XPoSE repository directory.
# This script runs the configured pseudobulk DESeq2 contrasts across cell types
# and saves full per-cell-type results tables plus DEG-count summaries.
# Gene Ontology analysis is handled separately in 04_relapse_transcriptional_characterization.

rm(list = ls())

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
})

# Paths -------------------------------------------------------------------
config_file <- "03_differential_expression/04_main_de/deseq2_run_config.csv"
output_root <- "output/03_differential_expression/04_main_de"

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

# Analysis settings -------------------------------------------------------
padj_threshold <- 0.05
cluster_col <- "cluster_name"
id_cols <- c("group", "experience", "sex", "region", cluster_col, "ratID")

# Helpers -----------------------------------------------------------------
prepare_deseq_data <- function(meta_df, counts_mat, cluster,
                               case_pair, ctrl_pair,
                               fact1, fact2, id_cols, cluster_col) {

  if (!cluster_col %in% names(meta_df)) {
    stop("Cluster column '", cluster_col, "' was not found in metadata.")
  }

  # Use a pre-existing sample_id when available. Otherwise reconstruct the
  # pseudobulk sample identifier from the same metadata fields used upstream.
  if (!"sample_id" %in% names(meta_df)) {
    missing_id_cols <- setdiff(id_cols, names(meta_df))
    if (length(missing_id_cols) > 0) {
      stop(
        "Cannot construct sample_id; missing metadata columns: ",
        paste(missing_id_cols, collapse = ", ")
      )
    }

    meta_df$sample_id <- do.call(
      paste,
      c(meta_df[id_cols], sep = "_")
    )
  }

  case_label <- paste(case_pair, collapse = "_")
  ctrl_label <- paste(ctrl_pair, collapse = "_")

  meta_use <- meta_df %>%
    filter(.data[[cluster_col]] == cluster) %>%
    filter(
      (.data[[fact1]] == case_pair[1] & .data[[fact2]] == case_pair[2]) |
        (.data[[fact1]] == ctrl_pair[1] & .data[[fact2]] == ctrl_pair[2])
    ) %>%
    mutate(
      cond = paste(.data[[fact1]], .data[[fact2]], sep = "_"),
      cond = factor(cond, levels = c(ctrl_label, case_label))
    ) %>%
    droplevels() %>%
    as.data.frame()

  rownames(meta_use) <- meta_use$sample_id

  common_ids <- intersect(colnames(counts_mat), rownames(meta_use))

  meta_use <- meta_use[common_ids, , drop = FALSE]
  counts_use <- counts_mat[, common_ids, drop = FALSE]

  if (!identical(colnames(counts_use), rownames(meta_use))) {
    stop("Pseudobulk count columns and metadata rows are not aligned.")
  }

  sample_counts <- table(meta_use$cond)

  if (length(sample_counts) != 2 || any(sample_counts < 2)) {
    stop(
      "Fewer than two pseudobulk samples remained in one or both conditions: ",
      paste(names(sample_counts), as.integer(sample_counts), sep = "=", collapse = ", ")
    )
  }

  list(
    meta = meta_use,
    counts = counts_use,
    case_label = case_label,
    ctrl_label = ctrl_label,
    n_case = unname(sample_counts[case_label]),
    n_ctrl = unname(sample_counts[ctrl_label])
  )
}

run_cluster_de <- function(dat, design_formula, padj_threshold) {

  if (!"cond" %in% all.vars(design_formula)) {
    stop("DESeq2 design must include 'cond'.")
  }

  dds <- DESeqDataSetFromMatrix(
    countData = dat$counts,
    colData = dat$meta,
    design = design_formula
  )

  dds <- DESeq(dds)

  res <- results(
    dds,
    contrast = c("cond", dat$case_label, dat$ctrl_label),
    alpha = padj_threshold,
    cooksCutoff = FALSE,
    independentFiltering = FALSE
  ) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    mutate(
      significant = !is.na(padj) & padj < padj_threshold,
      direction = case_when(
        significant & log2FoldChange > 0 ~ "up",
        significant & log2FoldChange < 0 ~ "down",
        TRUE ~ "not_significant"
      )
    )

  res
}

# Load run configuration --------------------------------------------------
run_config <- read.csv(
  config_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

names(run_config) <- trimws(names(run_config))
run_config <- run_config[
  !is.na(run_config$analysis_name) & nzchar(trimws(run_config$analysis_name)),
  ,
  drop = FALSE
]

required_cols <- c(
  "analysis_name",
  "case_fact1", "case_fact2",
  "ctrl_fact1", "ctrl_fact2",
  "fact1", "fact2",
  "design",
  "meta_rds", "pb_rds"
)

missing_cols <- setdiff(required_cols, names(run_config))
if (length(missing_cols) > 0) {
  stop(
    "Run config is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# Run configured contrasts -----------------------------------------------
all_summary <- list()
batch_log <- list()

for (i in seq_len(nrow(run_config))) {

  run <- run_config[i, ]
  analysis_name <- run$analysis_name

  message("Running ", analysis_name)

  analysis_dir <- file.path(output_root, analysis_name)
  results_dir <- file.path(analysis_dir, "results")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  run_status <- "success"
  run_error <- NA_character_

  tryCatch({

    meta_raw <- readRDS(run$meta_rds)
    pb_counts <- readRDS(run$pb_rds)
    counts_mat <- pb_counts[["RNA"]]

    if (is.null(counts_mat)) {
      stop("Pseudobulk RDS does not contain an 'RNA' count matrix.")
    }

    fact1 <- run$fact1
    fact2 <- run$fact2
    case_pair <- c(run$case_fact1, run$case_fact2)
    ctrl_pair <- c(run$ctrl_fact1, run$ctrl_fact2)
    design_formula <- as.formula(run$design)

    cluster_names <- sort(unique(meta_raw[[cluster_col]]))
    summary_rows <- vector("list", length(cluster_names))
    names(summary_rows) <- cluster_names

    for (cluster in cluster_names) {

      cluster_status <- "success"
      cluster_error <- NA_character_
      n_case <- NA_integer_
      n_ctrl <- NA_integer_
      n_tested <- NA_integer_
      n_sig <- NA_integer_
      n_up <- NA_integer_
      n_down <- NA_integer_

      tryCatch({

        dat <- prepare_deseq_data(
          meta_df = meta_raw,
          counts_mat = counts_mat,
          cluster = cluster,
          case_pair = case_pair,
          ctrl_pair = ctrl_pair,
          fact1 = fact1,
          fact2 = fact2,
          id_cols = id_cols,
          cluster_col = cluster_col
        )

        n_case <- dat$n_case
        n_ctrl <- dat$n_ctrl

        res <- run_cluster_de(
          dat = dat,
          design_formula = design_formula,
          padj_threshold = padj_threshold
        )

        n_tested <- nrow(res)
        n_up <- sum(res$direction == "up", na.rm = TRUE)
        n_down <- sum(res$direction == "down", na.rm = TRUE)
        n_sig <- n_up + n_down

        write.csv(
          res,
          file = file.path(
            results_dir,
            paste0(cluster, "_DESeq2_results.csv")
          ),
          row.names = FALSE
        )

      }, error = function(e) {
        cluster_status <<- "skipped"
        cluster_error <<- conditionMessage(e)
      })

      summary_rows[[cluster]] <- data.frame(
        analysis_name = analysis_name,
        cluster = cluster,
        case = paste(case_pair, collapse = "_"),
        control = paste(ctrl_pair, collapse = "_"),
        n_case = n_case,
        n_control = n_ctrl,
        n_tested = n_tested,
        n_sig = n_sig,
        n_up = n_up,
        n_down = n_down,
        status = cluster_status,
        reason = cluster_error,
        stringsAsFactors = FALSE
      )
    }

    deg_summary <- bind_rows(summary_rows)

    write.csv(
      deg_summary,
      file = file.path(analysis_dir, "DEG_summary_bycluster.csv"),
      row.names = FALSE
    )

    all_summary[[analysis_name]] <- deg_summary

  }, error = function(e) {
    run_status <<- "failed"
    run_error <<- conditionMessage(e)
  })

  batch_log[[analysis_name]] <- data.frame(
    analysis_name = analysis_name,
    status = run_status,
    reason = run_error,
    stringsAsFactors = FALSE
  )
}

# Combined summaries ------------------------------------------------------
if (length(all_summary) > 0) {
  write.csv(
    bind_rows(all_summary),
    file = file.path(output_root, "DEG_summary_all_contrasts.csv"),
    row.names = FALSE
  )
}

write.csv(
  bind_rows(batch_log),
  file = file.path(output_root, "batch_run_log.csv"),
  row.names = FALSE
)
