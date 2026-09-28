# XPoSE-seq: POC supplemental differential-expression analyses
#
# Supplemental analyses corresponding to Figure S4:
#   1. Compare pseudobulk DESeq2 and single-nucleus Wilcoxon DE results.
#   2. Assess DESeq2 robustness with leave-one-rat-out analysis.
#
# Run from the base XPoSE repository directory.

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
})

source("03_differential_expression/functions/poc_de_functions.R")

# Paths -------------------------------------------------------------------
input_file <- "output/01_metadata_clustering_qc/poc_combined_annotated.rds"
output_dir <- "output/03_differential_expression/02_poc_de_supp"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load POC object ---------------------------------------------------------
all <- readRDS(input_file)
DefaultAssay(all) <- "RNA"

if (inherits(all[["RNA"]], "Assay5")) {
  all <- JoinLayers(all, assay = "RNA")
}

# Standardize the groups used for POC DE without modifying the source object.
all$de_group <- case_when(
  all$experience == "HC" ~ "Homecage",
  all$experience == "NC" & tolower(all$group) == "active" ~ "Active",
  TRUE ~ NA_character_
)

# Retain populations meeting the primary POC DE inclusion requirement:
# >= 10 nuclei per rat and >= 3 rats in both Active and Homecage.
eligibility <- all@meta.data %>%
  filter(de_group %in% c("Active", "Homecage")) %>%
  count(cluster_name, de_group, ratID, name = "n_nuclei") %>%
  mutate(included = n_nuclei >= 10)

write.csv(
  eligibility,
  file.path(output_dir, "population_rat_nuclei_counts.csv"),
  row.names = FALSE
)

eligible_clusters <- eligibility %>%
  filter(included) %>%
  group_by(cluster_name, de_group) %>%
  summarise(n_rats = n_distinct(ratID), .groups = "drop") %>%
  filter(n_rats >= 3) %>%
  count(cluster_name, name = "n_groups") %>%
  filter(n_groups == 2) %>%
  pull(cluster_name)

all_rats <- sort(unique(all$ratID[all$de_group %in% c("Active", "Homecage")]))

# The shared single_factor_DESeq() uses counts > min_cell; therefore
# min_cell = 9 corresponds to the required minimum of 10 nuclei per rat.
min_cell_threshold <- 9

# Logs --------------------------------------------------------------------
skip_log <- data.frame(
  cluster = character(),
  excluded_rat = character(),
  reason = character(),
  stringsAsFactors = FALSE
)

error_log <- data.frame(
  cluster = character(),
  excluded_rat = character(),
  reason = character(),
  stringsAsFactors = FALSE
)

# Run full-data and leave-one-rat-out analyses ----------------------------
set.seed(22)

for (cl in eligible_clusters) {
  
  cluster_dir <- file.path(output_dir, paste0("Cluster_", cl))
  dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (excluded_rat in c("none", all_rats)) {
    
    exclusion_label <- if (excluded_rat == "none") {
      "Excluded_none"
    } else {
      paste0("Excluded_", excluded_rat)
    }
    
    exclusion_dir <- file.path(cluster_dir, exclusion_label)
    dir.create(exclusion_dir, recursive = TRUE, showWarnings = FALSE)
    
    subset_data <- if (excluded_rat == "none") {
      all
    } else {
      subset(all, subset = ratID != excluded_rat)
    }
    
    # The full-data analysis requires >=3 rats/group. After leaving one rat
    # out, >=2 rats/group is expected and is required for the robustness run.
    min_rats <- if (excluded_rat == "none") 3 else 2
    
    tryCatch({
      
      # Pseudobulk DESeq2 --------------------------------------------------
      de_results <- single_factor_DESeq(
        object = subset_data,
        comp_vect = c("group", "Active", "Homecage"),
        cluster = cl,
        min_cell = 10,
        min_rat = 2,
        keep_dds = TRUE
      )$results
      
      de_results <- as.data.frame(de_results)
      
      if (!"gene" %in% names(de_results)) {
        stop("DESeq2 results do not contain the expected 'gene' column.")
      }
      
      de_results$score <- ifelse(
        !is.na(de_results$padj) & de_results$padj < 0.05,
        1L,
        0L
      )
      
      de_results$score_updn <- ifelse(
        !is.na(de_results$padj) & de_results$padj < 0.05 &
          de_results$log2FoldChange > 0,
        1L,
        ifelse(
          !is.na(de_results$padj) & de_results$padj < 0.05 &
            de_results$log2FoldChange < 0,
          2L,
          0L
        )
      )
      
      write.csv(
        de_results,
        file.path(exclusion_dir, "deseq2_results.csv"),
        row.names = FALSE
      )
      
      # Single-nucleus Wilcoxon -------------------------------------------
      cluster_subset <- subset(
        subset_data,
        subset = cluster_name == cl & de_group %in% c("Active", "Homecage")
      )
      Idents(cluster_subset) <- "de_group"
      
      wilcoxon_results <- FindMarkers(
        cluster_subset,
        ident.1 = "Active",
        ident.2 = "Homecage",
        test.use = "wilcox"
      ) %>%
        rownames_to_column("gene")
      
      wilcoxon_results$score <- ifelse(
        !is.na(wilcoxon_results$p_val_adj) &
          wilcoxon_results$p_val_adj < 0.05,
        1L,
        0L
      )
      
      wilcoxon_results$score_updn <- ifelse(
        !is.na(wilcoxon_results$p_val_adj) &
          wilcoxon_results$p_val_adj < 0.05 &
          wilcoxon_results$avg_log2FC > 0,
        1L,
        ifelse(
          !is.na(wilcoxon_results$p_val_adj) &
            wilcoxon_results$p_val_adj < 0.05 &
            wilcoxon_results$avg_log2FC < 0,
          2L,
          0L
        )
      )
      
      write.csv(
        wilcoxon_results,
        file.path(exclusion_dir, "wilcoxon_results.csv"),
        row.names = FALSE
      )
      
    }, error = function(e) {
      if (grepl("^SKIP:", e$message)) {
        skip_log[nrow(skip_log) + 1, ] <<- list(cl, excluded_rat, e$message)
      } else {
        error_log[nrow(error_log) + 1, ] <<- list(cl, excluded_rat, e$message)
      }
    })
  }
}

if (nrow(skip_log) > 0) {
  write.csv(skip_log, file.path(output_dir, "skipped_analyses.csv"), row.names = FALSE)
}

if (nrow(error_log) > 0) {
  write.csv(error_log, file.path(output_dir, "analysis_errors.csv"), row.names = FALSE)
}

# Build gene x exclusion significance matrices ---------------------------
build_consistency_matrix <- function(cluster_dir, result_file) {
  
  exclusion_dirs <- list.dirs(cluster_dir, recursive = FALSE)
  exclusion_dirs <- exclusion_dirs[grepl("^Excluded_", basename(exclusion_dirs))]
  
  score_list <- lapply(exclusion_dirs, function(d) {
    file <- file.path(d, result_file)
    if (!file.exists(file)) return(NULL)
    
    results <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(c("gene", "score_updn") %in% names(results))) return(NULL)
    
    out <- results[, c("gene", "score_updn"), drop = FALSE]
    names(out)[2] <- basename(d)
    out
  })
  
  score_list <- Filter(Negate(is.null), score_list)
  if (length(score_list) == 0) return(NULL)
  
  combined <- Reduce(
    function(x, y) full_join(x, y, by = "gene"),
    score_list
  )
  
  combined[is.na(combined)] <- 0
  combined
}

# Classify genes as consistent or non-consistent across exclusions --------
classify_consistency <- function(matrix_df) {
  
  if (is.null(matrix_df) || nrow(matrix_df) == 0) return(NULL)
  
  genes <- matrix_df$gene
  scores <- as.matrix(matrix_df[, setdiff(names(matrix_df), "gene"), drop = FALSE])
  storage.mode(scores) <- "numeric"
  
  baseline_col <- match("Excluded_none", colnames(scores))
  if (is.na(baseline_col)) {
    stop("Consistency matrix is missing Excluded_none.")
  }
  
  baseline <- scores[, baseline_col]
  
  consistent_up <- genes[apply(scores, 1, function(x) all(x == 1))]
  consistent_down <- genes[apply(scores, 1, function(x) all(x == 2))]
  
  nonconsistent_up <- genes[
    baseline == 1 & apply(scores, 1, function(x) any(x == 1) && !all(x == 1))
  ]
  
  nonconsistent_down <- genes[
    baseline == 2 & apply(scores, 1, function(x) any(x == 2) && !all(x == 2))
  ]
  
  list(
    consistent_up = consistent_up,
    consistent_down = consistent_down,
    nonconsistent_up = nonconsistent_up,
    nonconsistent_down = nonconsistent_down
  )
}

save_gene_subset <- function(matrix_df, genes, file) {
  out <- matrix_df[matrix_df$gene %in% genes, , drop = FALSE]
  write.csv(out, file, row.names = FALSE)
}

# Summaries for supplemental figure ---------------------------------------
method_summary <- data.frame()
consistency_summary <- data.frame()

for (cl in eligible_clusters) {
  
  cluster_dir <- file.path(output_dir, paste0("Cluster_", cl))
  
  deseq_matrix <- build_consistency_matrix(cluster_dir, "deseq2_results.csv")
  wilcoxon_matrix <- build_consistency_matrix(cluster_dir, "wilcoxon_results.csv")
  
  if (!is.null(deseq_matrix)) {
    write.csv(
      deseq_matrix,
      file.path(cluster_dir, "deseq2_leave_one_out_matrix.csv"),
      row.names = FALSE
    )
  }
  
  if (!is.null(wilcoxon_matrix)) {
    write.csv(
      wilcoxon_matrix,
      file.path(cluster_dir, "wilcoxon_leave_one_out_matrix.csv"),
      row.names = FALSE
    )
  }
  
  # Figure S4A: full-data DESeq2 vs Wilcoxon DEG counts ------------------
  deseq_baseline <- file.path(cluster_dir, "Excluded_none", "deseq2_results.csv")
  wilcoxon_baseline <- file.path(cluster_dir, "Excluded_none", "wilcoxon_results.csv")
  
  count_up_down <- function(file) {
    if (!file.exists(file)) return(c(up = 0L, down = 0L))
    x <- read.csv(file, stringsAsFactors = FALSE)
    if (!"score_updn" %in% names(x)) return(c(up = 0L, down = 0L))
    c(
      up = sum(x$score_updn == 1, na.rm = TRUE),
      down = sum(x$score_updn == 2, na.rm = TRUE)
    )
  }
  
  d <- count_up_down(deseq_baseline)
  w <- count_up_down(wilcoxon_baseline)
  
  method_summary <- bind_rows(
    method_summary,
    data.frame(
      cluster = cl,
      DESeq2_up = unname(d["up"]),
      DESeq2_down = unname(d["down"]),
      Wilcoxon_up = unname(w["up"]),
      Wilcoxon_down = unname(w["down"])
    )
  )
  
  # Figure S4B: leave-one-rat-out DESeq2 consistency ----------------------
  if (!is.null(deseq_matrix)) {
    consistency <- classify_consistency(deseq_matrix)
    
    save_gene_subset(
      deseq_matrix,
      consistency$consistent_up,
      file.path(cluster_dir, "deseq2_consistent_up.csv")
    )
    save_gene_subset(
      deseq_matrix,
      consistency$consistent_down,
      file.path(cluster_dir, "deseq2_consistent_down.csv")
    )
    save_gene_subset(
      deseq_matrix,
      consistency$nonconsistent_up,
      file.path(cluster_dir, "deseq2_nonconsistent_up.csv")
    )
    save_gene_subset(
      deseq_matrix,
      consistency$nonconsistent_down,
      file.path(cluster_dir, "deseq2_nonconsistent_down.csv")
    )
    
    consistency_summary <- bind_rows(
      consistency_summary,
      data.frame(
        cluster = cl,
        consistent_up = length(consistency$consistent_up),
        consistent_down = length(consistency$consistent_down),
        nonconsistent_up = length(consistency$nonconsistent_up),
        nonconsistent_down = length(consistency$nonconsistent_down)
      )
    )
  }
}

write.csv(
  method_summary,
  file.path(output_dir, "summary_deseq2_vs_wilcoxon.csv"),
  row.names = FALSE
)

write.csv(
  consistency_summary,
  file.path(output_dir, "summary_leave_one_out_deseq2.csv"),
  row.names = FALSE
)
