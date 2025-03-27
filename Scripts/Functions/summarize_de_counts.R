#' Title
#'
#' @param directory directory with the deg output files
#' @param clusters cell type clusters
#' @param gene_list character vector of genes
#'
#' @returns
#' @export
#'
#' @examples
summarize_DEcounts <- function(directory, clusters, gene_list) {
  
  # Convert gene_list to a character vector
  gene_list <- unlist(strsplit(gene_list, ","))
  
  # Initialize empty lists to store results
  log2FC_HC_list <- list()
  log2FC_Active_list <- list()
  adjpval_list <- list()
  
  # Loop through clusters and process each file
  for (cluster in clusters) {
    file_path <- file.path(directory, paste0("Log2FC_Values_", cluster, "_group_Active_Homecage.csv"))
    
    if (!file.exists(file_path)) {
      warning("File not found: ", file_path)
      next
    }
    
    data <- read.csv(file_path, row.names = 1)
    
    # Select relevant columns and order them
    log2FC_HC <- data[, grep("log2FC_HC.[1-4].Homecage", colnames(data)), drop = FALSE]
    log2FC_HC <- log2FC_HC[, order(as.numeric(gsub(".*HC.([1-4]).Homecage", "\\1", colnames(log2FC_HC))))]
    
    log2FC_Active <- data[, grep("log2FC_NC.[1-4].Active", colnames(data)), drop = FALSE]
    log2FC_Active <- log2FC_Active[, order(as.numeric(gsub(".*NC.([1-4]).Active", "\\1", colnames(log2FC_Active))))]
    
    adjpval <- data[, "adjpval", drop = FALSE]
    
    # Check for missing gene names
    missing_genes <- setdiff(gene_list, rownames(data))
    if (length(missing_genes) > 0) {
      warning("The following genes are not found in the data for cluster ", cluster, ": ", paste(missing_genes, collapse = ", "))
    }
    
    # Reorder data frames based on gene_list
    log2FC_HC <- log2FC_HC[gene_list, , drop = FALSE]
    log2FC_Active <- log2FC_Active[gene_list, , drop = FALSE]
    adjpval <- adjpval[gene_list, , drop = FALSE]
    
    # Store results
    log2FC_HC_list[[cluster]] <- log2FC_HC
    log2FC_Active_list[[cluster]] <- log2FC_Active
    adjpval_list[[cluster]] <- adjpval
  }
  
  # Combine results across clusters
  log2FC_HC_summary <- do.call(cbind, log2FC_HC_list)
  log2FC_Active_summary <- do.call(cbind, log2FC_Active_list)
  adjpval_summary <- do.call(cbind, adjpval_list)
  
  # Save results
  write.csv(log2FC_HC_summary, file = file.path(directory, "Log2FC_summary_HC.csv"), row.names = TRUE)
  write.csv(log2FC_Active_summary, file = file.path(directory, "Log2FC_summary_Active.csv"), row.names = TRUE)
  write.csv(adjpval_summary, file = file.path(directory, "adjpval_summary.csv"), row.names = TRUE)
  
  message("Summary files saved in ", directory)
}