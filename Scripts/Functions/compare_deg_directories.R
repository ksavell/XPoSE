#' Compare deg directories
#'
#' @param dir1 directory of deg comparisons
#' @param dir2 another directory of deg comparisons
#' @param prefixes usually cluster names
#' @param output_dir define where to save
#' @param dir1_name name to use in plots for dir 1
#' @param dir2_name name to use in plots for dir 2
#'
#' @returns
#' @export
#'
#' @examples
compare_deg_directories <- function(dir1, dir2, prefixes = NULL, output_dir = NULL, 
                                    dir1_name = "Directory 1", dir2_name = "Directory 2") {
  all_results <- list()
  skipped <- c()
  
  for (prefix in prefixes) {
    # Match files that start with "prefix_"
    files1 <- list.files(dir1, pattern = paste0("^", prefix, "_"), full.names = TRUE)
    files2 <- list.files(dir2, pattern = paste0("^", prefix, "_"), full.names = TRUE)
    
    message("\nChecking prefix: ", prefix)
    message("  Found in dir1: ", length(files1), " file(s)")
    message("  Found in dir2: ", length(files2), " file(s)")
    
    if (length(files1) == 0 || length(files2) == 0) {
      warning("Missing file(s) for prefix: ", prefix)
      skipped <- c(skipped, prefix)
      all_results[[prefix]] <- list(data = NULL, plot = NULL, status = "Skipped - Missing files")
      next
    }
    
    # Use first matching file per dir
    file1 <- files1[1]
    file2 <- files2[1]
    
    # Read and filter columns
    df1 <- read.csv(file1)[, c("gene", "log2FoldChange", "padj")]
    df2 <- read.csv(file2)[, c("gene", "log2FoldChange", "padj")]
    
    # Merge and filter for significance
    merged <- merge(df1, df2, by = "gene", suffixes = c("_dir1", "_dir2")) %>%
      filter(padj_dir1 < 0.05 | padj_dir2 < 0.05) %>%
      mutate(significance = case_when(
        padj_dir1 < 0.05 & padj_dir2 < 0.05 ~ "Both Significant",
        padj_dir1 < 0.05 ~ paste(dir1_name, "Significant"),
        padj_dir2 < 0.05 ~ paste(dir2_name, "Significant"),
        TRUE ~ "Not Significant"
      ))
    
    # Split into up/downregulated
    upregulated <- merged %>%
      filter((padj_dir1 < 0.05 & log2FoldChange_dir1 > 0) |
               (padj_dir2 < 0.05 & log2FoldChange_dir2 > 0))
    
    downregulated <- merged %>%
      filter((padj_dir1 < 0.05 & log2FoldChange_dir1 < 0) |
               (padj_dir2 < 0.05 & log2FoldChange_dir2 < 0))
    
    # Define custom colors
    sig_colors <- setNames(
      c("black", "blue", "red", "gray"),
      c("Both Significant", 
        paste(dir1_name, "Significant"),
        paste(dir2_name, "Significant"),
        "Not Significant")
    )
    
    # Plot comparison
    p <- ggplot(merged, aes(x = log2FoldChange_dir1, y = log2FoldChange_dir2, color = significance)) +
      geom_point(alpha = 0.7) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
      scale_color_manual(values = sig_colors) +
      labs(
        title = paste("Comparison for prefix:", prefix),
        x = paste("Log2 Fold Change (", dir1_name, ")"),
        y = paste("Log2 Fold Change (", dir2_name, ")")
      ) +
      theme_classic()
    
    all_results[[prefix]] <- list(
      data = merged,
      upregulated = upregulated,
      downregulated = downregulated,
      plot = p,
      file1 = file1,
      file2 = file2
    )
    
    # Save outputs
    if (!is.null(output_dir)) {
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
      
      write.csv(merged, file.path(output_dir, paste0("comparison_", prefix, ".csv")), row.names = FALSE)
      write.csv(upregulated, file.path(output_dir, paste0("comparison_", prefix, "_upregulated.csv")), row.names = FALSE)
      write.csv(downregulated, file.path(output_dir, paste0("comparison_", prefix, "_downregulated.csv")), row.names = FALSE)
      
      ggsave(file.path(output_dir, paste0("plot_", prefix, ".png")), plot = p, width = 8, height = 6)
    }
  }
  
  if (length(skipped) > 0) {
    message("\nSkipped prefixes: ", paste(skipped, collapse = ", "))
  }
  
  return(list(results = all_results, skipped = skipped))
}
