# F5

source("Scripts/Functions/single_factor_DESeq.R")
load("~/Projects/XPoSE/all_10312024.RData")


clusters <- unique(all$cluster_name)
all$experience <- ifelse(all$group == "Homecage", "HC", "NC")

pairs_list <- list(
  c("experience","NC","HC"),
  c("group", "Non-active", "Homecage"),
  c("group", "Active", "Homecage"),
  c("group", "Active", "Non-active")
)

for (pair in pairs_list) {
  results_list <- list()  # Reset results list for each comparison
  i <- 1  # Index for list tracking
  
  for (cl in clusters) {
    # Run DESeq2 analysis
    deseq2_results <- tryCatch({
      single_factor_DESeq(object = all, comp_vect = pair, cluster = cl, min_cell = 1)
    }, error = function(e) NULL)
    
    if (is.null(deseq2_results) || !("deseq_results" %in% names(deseq2_results))) {
      next
    }
    
    deseq2_results_tbl <- deseq2_results[["deseq_results"]]
    
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
  }
  
  # Combine results for this pair and save
  if (length(results_list) > 0) {
    results <- do.call(rbind, results_list)
    file_name <- paste0(pair[1], "_", pair[2], "_", pair[3], ".csv")
    write.csv(results, file_name, row.names = FALSE)
  }
}


