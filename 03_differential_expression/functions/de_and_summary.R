#' run de and summarize
#'
#' @param seur_obj seurat object
#' @param pair can be a character vector or a list element. first element is
#' the factor, second element is the experimental group, and third element is
#' the control group
#'
#' @returns
#' @export
#'
#' @examples
de_and_summary <- function(seur_obj, pair) {
  results_list <- list()  # Reset results list for each comparison
  i <- 1  # Index for list tracking
  
  for (cl in clusters) {
    # Run DESeq2 analysis
    deseq2_results <- tryCatch({
      single_factor_DESeq(object = all, comp_vect = pair, cluster = cl, min_cell = 0)
    }, error = function(e) {
      message(paste("DESeq2 failed for cluster:", cl, "with pair:", paste(pair, collapse = "_")))
      NULL
    })
    
    
    deseq2_results_tbl <- deseq2_results[["results"]]
    dds <- deseq2_results[["dds"]]
    
    if (!("padj" %in% colnames(deseq2_results_tbl)) || !("log2FoldChange" %in% colnames(deseq2_results_tbl))) {
      warning(paste("Columns 'padj' or 'log2FoldChange' missing for cluster:", cl))
      next
    }
    
    deseq2_results_tbl$padj <- as.numeric(deseq2_results_tbl$padj)
    deseq2_results_tbl$log2FoldChange <- as.numeric(deseq2_results_tbl$log2FoldChange)
    
    deseq2_results_tbl$up_score <- ifelse(!is.na(deseq2_results_tbl$padj) & deseq2_results_tbl$padj < 0.05 & deseq2_results_tbl$log2FoldChange > 0, 1, 0)
    deseq2_results_tbl$dn_score <- ifelse(!is.na(deseq2_results_tbl$padj) & deseq2_results_tbl$padj < 0.05 & deseq2_results_tbl$log2FoldChange < 0, -1, 0)
    
    results_list[[i]] <- data.frame(Category = paste(pair, collapse = "_"), Observation = cl, Value = sum(deseq2_results_tbl$up_score, na.rm = TRUE))
    results_list[[i + 1]] <- data.frame(Category = paste(pair, collapse = "_"), Observation = cl, Value = sum(deseq2_results_tbl$dn_score, na.rm = TRUE))
    i <- i + 2
    
    fn <- paste0(cl, "_", paste(pair, collapse = "_"))
    if (nrow(deseq2_results_tbl) > 0) {
      write.csv(deseq2_results_tbl, paste0(fn, ".csv"), row.names = FALSE)
      saveRDS(dds, paste0(fn,".RDS"))
    } else {
      message(paste("Skipping file save: No DESeq2 results for", fn))
    }
  }
  
  # Combine results for this pair and save
  if (length(results_list) > 0) {
    results <- do.call(rbind, results_list)
    file_name <- paste0(pair[1], "_", pair[2], "_", pair[3], ".csv")
    write.csv(results, file_name, row.names = FALSE)
  }
}