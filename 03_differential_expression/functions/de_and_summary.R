#' Run single-factor pseudobulk DE across cell populations and summarize DEGs
#'
#' @param seur_obj Seurat object containing cluster_name, ratID, and comparison metadata.
#' @param pair Character vector: factor, experimental level, control level.
#' @param clusters Cell populations to test. Defaults to all cluster_name values.
#' @param min_cell Minimum nuclei required per rat/population.
#' @param min_rat Minimum rats required per comparison group after filtering.
#' @param output_dir Directory for per-population DE results and summary output.
#' @param save_dds Logical; save DESeq2 objects for successful comparisons.
#'
#' @return Invisibly returns the signed DEG-count summary table.
de_and_summary <- function(seur_obj,
                           pair,
                           clusters = sort(unique(seur_obj$cluster_name)),
                           min_cell = 10,
                           min_rat = 3,
                           output_dir = ".",
                           save_dds = TRUE) {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  comparison_tag <- paste(pair, collapse = "_")
  results_list <- list()
  result_i <- 1
  
  for (cl in clusters) {
    
    deseq2_results <- tryCatch(
      single_factor_DESeq(
        object = seur_obj,
        comp_vect = pair,
        cluster = cl,
        min_cell = min_cell,
        min_rat = min_rat
      ),
      error = function(e) {
        message(
          "Skipping ", cl, " for ", comparison_tag, ": ",
          conditionMessage(e)
        )
        NULL
      }
    )
    
    if (is.null(deseq2_results)) {
      next
    }
    
    de_tbl <- deseq2_results$results
    dds <- deseq2_results$dds
    
    if (!all(c("padj", "log2FoldChange") %in% colnames(de_tbl))) {
      warning("padj/log2FoldChange missing for ", cl, "; skipping.")
      next
    }
    
    de_tbl$padj <- as.numeric(de_tbl$padj)
    de_tbl$log2FoldChange <- as.numeric(de_tbl$log2FoldChange)
    
    n_up <- sum(
      !is.na(de_tbl$padj) &
        de_tbl$padj < 0.05 &
        de_tbl$log2FoldChange > 0
    )
    
    n_down <- sum(
      !is.na(de_tbl$padj) &
        de_tbl$padj < 0.05 &
        de_tbl$log2FoldChange < 0
    )
    
    # Preserve the original F5A convention: upregulated DEG counts are
    # positive and downregulated DEG counts are negative.
    results_list[[result_i]] <- data.frame(
      Category = comparison_tag,
      Observation = cl,
      Value = n_up,
      stringsAsFactors = FALSE
    )
    
    results_list[[result_i + 1]] <- data.frame(
      Category = comparison_tag,
      Observation = cl,
      Value = -n_down,
      stringsAsFactors = FALSE
    )
    
    result_i <- result_i + 2
    
    file_tag <- paste0(cl, "_", comparison_tag)
    
    write.csv(
      de_tbl,
      file.path(output_dir, paste0(file_tag, ".csv")),
      row.names = FALSE
    )
    
    if (save_dds) {
      saveRDS(
        dds,
        file.path(output_dir, paste0(file_tag, "_dds.rds"))
      )
    }
  }
  
  if (length(results_list) == 0) {
    warning("No successful DE comparisons for ", comparison_tag)
    return(invisible(NULL))
  }
  
  summary_df <- do.call(rbind, results_list)
  
  write.csv(
    summary_df,
    file.path(output_dir, paste0(comparison_tag, "_summary.csv")),
    row.names = FALSE
  )
  
  invisible(summary_df)
}
