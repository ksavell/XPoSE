# XPoSE-seq: differential expression for the main dataset
#
# Run from the base XPoSE repository directory.
# This script runs configured pseudobulk DESeq2 contrasts across cell types,
# applies minimum nuclei/replicate requirements, and saves full DE results,
# DEG-count summaries, and per-rat nuclei-count audit tables.

rm(list = ls())

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
})

# Paths -------------------------------------------------------------------
config_file <- "03_differential_expression/03_main_de/deseq2_run_config.csv"
output_root <- "output/03_differential_expression/03_main_de"

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

# Analysis settings -------------------------------------------------------
padj_threshold <- 0.05
cluster_col <- "cluster_name"
min_nuclei_per_sample <- 10
min_rats_per_group <- 3

# Helpers -----------------------------------------------------------------
make_condition_label <- function(pair) {
  paste(pair, collapse = "_")
}

get_condition_rats <- function(meta_df, fact1, fact2, pair) {
  meta_df %>%
    filter(
      .data[[fact1]] == pair[1],
      .data[[fact2]] == pair[2]
    ) %>%
    distinct(ratID) %>%
    pull(ratID) %>%
    as.character()
}

build_nuclei_audit <- function(meta_df, cluster, case_pair, ctrl_pair,
                               fact1, fact2, design_formula,
                               min_nuclei, min_rats) {
  
  case_label <- make_condition_label(case_pair)
  ctrl_label <- make_condition_label(ctrl_pair)
  paired_design <- "ratID" %in% all.vars(design_formula)
  
  case_rats <- get_condition_rats(meta_df, fact1, fact2, case_pair)
  ctrl_rats <- get_condition_rats(meta_df, fact1, fact2, ctrl_pair)
  
  if (paired_design) {
    pair_rats <- sort(unique(c(case_rats, ctrl_rats)))
    expected <- bind_rows(
      data.frame(
        ratID = pair_rats,
        comparison_group = "case",
        cond = case_label,
        stringsAsFactors = FALSE
      ),
      data.frame(
        ratID = pair_rats,
        comparison_group = "control",
        cond = ctrl_label,
        stringsAsFactors = FALSE
      )
    )
  } else {
    expected <- bind_rows(
      data.frame(
        ratID = case_rats,
        comparison_group = "case",
        cond = case_label,
        stringsAsFactors = FALSE
      ),
      data.frame(
        ratID = ctrl_rats,
        comparison_group = "control",
        cond = ctrl_label,
        stringsAsFactors = FALSE
      )
    )
  }
  
  cluster_meta <- meta_df %>%
    filter(.data[[cluster_col]] == cluster) %>%
    filter(
      (.data[[fact1]] == case_pair[1] & .data[[fact2]] == case_pair[2]) |
        (.data[[fact1]] == ctrl_pair[1] & .data[[fact2]] == ctrl_pair[2])
    ) %>%
    mutate(
      cond = paste(.data[[fact1]], .data[[fact2]], sep = "_"),
      comparison_group = if_else(cond == case_label, "case", "control")
    ) %>%
    select(
      ratID,
      comparison_group,
      cond,
      sample_id,
      n_nuclei,
      any_of(c("sex", "experience", "group", "region"))
    ) %>%
    mutate(
      ratID = as.character(ratID),
      n_nuclei = as.integer(n_nuclei)
    )
  
  duplicate_keys <- cluster_meta %>%
    count(ratID, cond) %>%
    filter(n > 1)
  
  if (nrow(duplicate_keys) > 0) {
    stop(
      "Multiple pseudobulk samples were found for the same rat and condition in cluster ",
      cluster, "."
    )
  }
  
  audit <- expected %>%
    left_join(
      cluster_meta,
      by = c("ratID", "comparison_group", "cond")
    ) %>%
    mutate(
      n_nuclei = if_else(is.na(n_nuclei), 0L, as.integer(n_nuclei)),
      meets_min_nuclei = n_nuclei >= min_nuclei
    )
  
  if (paired_design) {
    complete_rats <- audit %>%
      group_by(ratID) %>%
      summarise(
        n_conditions = n_distinct(cond),
        all_conditions_meet_min =
          n_conditions == 2 & all(meets_min_nuclei),
        .groups = "drop"
      ) %>%
      filter(all_conditions_meet_min) %>%
      pull(ratID)
    
    audit <- audit %>%
      mutate(
        complete_pair = ratID %in% complete_rats,
        included_in_de = meets_min_nuclei & complete_pair
      )
  } else {
    audit <- audit %>%
      mutate(
        complete_pair = NA,
        included_in_de = meets_min_nuclei
      )
  }
  
  included_counts <- audit %>%
    filter(included_in_de) %>%
    count(cond, name = "n_rats")
  
  n_case <- included_counts$n_rats[match(case_label, included_counts$cond)]
  n_ctrl <- included_counts$n_rats[match(ctrl_label, included_counts$cond)]
  
  if (length(n_case) == 0 || is.na(n_case)) n_case <- 0L
  if (length(n_ctrl) == 0 || is.na(n_ctrl)) n_ctrl <- 0L
  
  eligible <- n_case >= min_rats && n_ctrl >= min_rats
  
  reason <- if (eligible) {
    NA_character_
  } else {
    paste0(
      "Requires >=", min_rats, " rats per group with >=", min_nuclei,
      " nuclei per pseudobulk sample; retained case=", n_case,
      ", control=", n_ctrl,
      if (paired_design) " after requiring complete pairs." else "."
    )
  }
  
  audit %>%
    mutate(
      cluster = cluster,
      case = case_label,
      control = ctrl_label,
      min_nuclei_required = min_nuclei,
      min_rats_required = min_rats,
      paired_design = paired_design
    ) %>%
    select(
      cluster,
      comparison_group,
      cond,
      ratID,
      any_of(c("sex", "experience", "group", "region")),
      sample_id,
      n_nuclei,
      meets_min_nuclei,
      complete_pair,
      included_in_de,
      case,
      control,
      paired_design,
      min_nuclei_required,
      min_rats_required
    ) %>%
    arrange(cluster, comparison_group, ratID) %>%
    structure(
      n_case = as.integer(n_case),
      n_ctrl = as.integer(n_ctrl),
      eligible = eligible,
      reason = reason
    )
}

prepare_deseq_data <- function(meta_df, counts_mat, cluster,
                               case_pair, ctrl_pair,
                               fact1, fact2, design_formula,
                               min_nuclei, min_rats) {
  
  if (!cluster_col %in% names(meta_df)) {
    stop("Cluster column '", cluster_col, "' was not found in metadata.")
  }
  
  required_meta <- c(
    "sample_id", "ratID", "n_nuclei", cluster_col,
    fact1, fact2
  )
  missing_meta <- setdiff(required_meta, names(meta_df))
  if (length(missing_meta) > 0) {
    stop(
      "Pseudobulk metadata are missing required columns: ",
      paste(missing_meta, collapse = ", ")
    )
  }
  
  case_label <- make_condition_label(case_pair)
  ctrl_label <- make_condition_label(ctrl_pair)
  
  audit <- build_nuclei_audit(
    meta_df = meta_df,
    cluster = cluster,
    case_pair = case_pair,
    ctrl_pair = ctrl_pair,
    fact1 = fact1,
    fact2 = fact2,
    design_formula = design_formula,
    min_nuclei = min_nuclei,
    min_rats = min_rats
  )
  
  n_case <- attr(audit, "n_case")
  n_ctrl <- attr(audit, "n_ctrl")
  eligible <- isTRUE(attr(audit, "eligible"))
  reason <- attr(audit, "reason")
  
  if (!eligible) {
    return(list(
      eligible = FALSE,
      reason = reason,
      audit = audit,
      n_case = n_case,
      n_ctrl = n_ctrl,
      case_label = case_label,
      ctrl_label = ctrl_label
    ))
  }
  
  keep_sample_ids <- audit %>%
    filter(included_in_de, !is.na(sample_id)) %>%
    pull(sample_id)
  
  meta_use <- meta_df %>%
    filter(.data[[cluster_col]] == cluster) %>%
    filter(sample_id %in% keep_sample_ids) %>%
    mutate(
      cond = paste(.data[[fact1]], .data[[fact2]], sep = "_"),
      cond = factor(cond, levels = c(ctrl_label, case_label)),
      ratID = factor(ratID)
    ) %>%
    droplevels() %>%
    as.data.frame()
  
  rownames(meta_use) <- meta_use$sample_id
  
  common_ids <- colnames(counts_mat)[
    colnames(counts_mat) %in% rownames(meta_use)
  ]
  
  meta_use <- meta_use[common_ids, , drop = FALSE]
  counts_use <- counts_mat[, common_ids, drop = FALSE]
  
  if (!identical(colnames(counts_use), rownames(meta_use))) {
    stop("Pseudobulk count columns and metadata rows are not aligned.")
  }
  
  if (ncol(counts_use) == 0) {
    stop("No eligible pseudobulk samples remained after nuclei-count filtering.")
  }
  
  list(
    eligible = TRUE,
    reason = NA_character_,
    audit = audit,
    meta = meta_use,
    counts = counts_use,
    case_label = case_label,
    ctrl_label = ctrl_label,
    n_case = n_case,
    n_ctrl = n_ctrl
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
  
  results(
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
    
    cluster_names <- sort(unique(as.character(meta_raw[[cluster_col]])))
    cluster_names <- cluster_names[!is.na(cluster_names) & nzchar(cluster_names)]
    
    summary_rows <- vector("list", length(cluster_names))
    names(summary_rows) <- cluster_names
    audit_rows <- vector("list", length(cluster_names))
    names(audit_rows) <- cluster_names
    
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
          design_formula = design_formula,
          min_nuclei = min_nuclei_per_sample,
          min_rats = min_rats_per_group
        )
        
        n_case <- dat$n_case
        n_ctrl <- dat$n_ctrl
        
        audit_rows[[cluster]] <- dat$audit %>%
          mutate(analysis_name = analysis_name, .before = 1)
        
        if (!dat$eligible) {
          cluster_status <- "skipped"
          cluster_error <- dat$reason
        } else {
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
        }
        
      }, error = function(e) {
        cluster_status <<- "skipped"
        cluster_error <<- conditionMessage(e)
      })
      
      summary_rows[[cluster]] <- data.frame(
        analysis_name = analysis_name,
        cluster = cluster,
        case = make_condition_label(case_pair),
        control = make_condition_label(ctrl_pair),
        min_nuclei_per_sample = min_nuclei_per_sample,
        min_rats_per_group = min_rats_per_group,
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
    
    if (length(audit_rows) > 0) {
      nuclei_audit <- bind_rows(audit_rows)
      
      write.csv(
        nuclei_audit,
        file = file.path(
          analysis_dir,
          "nuclei_counts_by_celltype_rat_group.csv"
        ),
        row.names = FALSE
      )
    }
    
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
