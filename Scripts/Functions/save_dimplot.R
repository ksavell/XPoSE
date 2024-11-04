#' save dimplots for figures with specified colors
#'
#' @param seur_obj seurat object
#' @param groupby 
#' @param splitby 
#' @param file_n file name, if applying split or group it will auto incorportate
#' @param hex_list list of hex codes that are named by a 'groupby' value
#'
#' @return
#' @export
#'
#' @examples
save_dimplot <- function(seur_obj, 
                         groupby = NULL,  # Explicit groupby argument
                         splitby = NULL, 
                         file_n = NULL, 
                         hex_list = NULL) {
  
  # Ensure groupby and hex_list are provided
  if (is.null(groupby)) {
    stop("You must specify a groupby argument.")
  }
  if (is.null(hex_list) || !groupby %in% names(hex_list)) {
    stop("You must provide a valid hex list for the specified groupby.")
  }
  
  # Get the hex color mapping for the specified groupby
  hex_colors <- hex_list[[groupby]]
  
  # Get the unique values from the splitby column if it exists
  split_values <- if (!is.null(splitby)) unique(seur_obj@meta.data[[splitby]]) else NULL
  
  # Handle cases where splitby is NULL (no splitting) or has a valid value
  if (!is.null(split_values)) {
    # Loop through each value in the split.by field
    for (split_val in split_values) {
      # Subset the Seurat object for each split value
      subset_obj <- seur_obj[, seur_obj@meta.data[[splitby]] == split_val]
      
      # Generate the DimPlot for the specific groupby and split
      current_plot <- DimPlot(subset_obj, 
                              group.by = groupby, 
                              cols = hex_colors[names(hex_colors) %in% unique(subset_obj@meta.data[[groupby]])],
                              pt.size = 0.5,
                              shuffle = T) + 
        ggtitle(paste(groupby, "-", split_val)) +
        theme_void() +
        theme(legend.position = "none")
      
      # Define the output file name with the current split value and groupby
      output_file <- paste0(file_n, "_", split_val, "_", groupby, ".pdf")
      
      # Open a PDF file for each split and save the plot
      pdf(file = output_file, width = 10, height = 10)  # Set a fixed size for the PDFs
      
      # Print the specific plot to the PDF
      print(current_plot)
      
      # Close the PDF file
      dev.off()
    }
  } else {
    # If no splitting, generate a single plot for the groupby column
    current_plot <- DimPlot(seur_obj, 
                            group.by = groupby, 
                            cols = hex_colors[names(hex_colors) %in% unique(seur_obj@meta.data[[groupby]])],
                            pt.size = 0.5,
                            shuffle = T) + 
      ggtitle(groupby) +
      theme_void() +
      theme(legend.position = "none")
    
    # Define the output file name without splitting
    output_file <- paste0(file_n, "_", groupby, ".pdf")
    
    # Open a PDF file and save the specific plot
    pdf(file = output_file, width = 10, height = 10)  # Set a fixed size for the PDFs
    
    # Print the specific plot to the PDF
    print(current_plot)
    
    # Close the PDF file
    dev.off()
  }
}