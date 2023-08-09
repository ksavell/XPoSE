#' Plot DimPlot with no formatting for Figures
#'
#' @param seur_obj seurat object
#' @param splitby factor to split umap plot
#' @param splitl number of splits so the width is calculated correctly
#' @param file_n name the file, must include .PDF
#' @param glutcol logical, T for glut cluster colors
#' @param gabacol logical, T for gaba cluster colors
#' @param groupcol logical, T for group cluster colors
#' @param groupby factor to group the umap plot
#'
#' @return saves PDF with specific colors and split/group settings
#' @export
#'
#' @examples
save_dimplot <- function(seur_obj, groupby = NULL, splitby = NULL, splitl = 1, 
                         file_n = NULL, glutcol = FALSE,
                         gabacol = FALSE, groupcol = FALSE,
                         sexcol = FALSE){
  
  glut_hex <- c('#0A5B8C','#64C7C8',
                '#2C8CB9','#41B75F',
                '#6F499D','#3C9E64')
  
  gaba_hex <- c('#F8991D','#E66027',
                '#C03C82','#C52126',
                '#B0B235','#DB808C',
                '#A669AB')
  
  grp_hex <- c('#ae1e5b','#e6e6e6','gray60')
  
  sex_hex <- c("#2C5F2D","#97BC62FF")
  
  # Open the PDF file
  pdf(file = file_n,
      width = 10*splitl,
      height = 10)
  
  plot_list <- list()  # Create a list to store modified plots
  
  if (glutcol == TRUE) {
    plot <- DimPlot(seur_obj, 
                    group.by = groupby, 
                    split.by = splitby, 
                    cols = glut_hex)
    plot_list[["glut"]] <- plot
  }
  
  if (gabacol == TRUE) {
    plot <- DimPlot(seur_obj, 
                    group.by = groupby, 
                    split.by = splitby, 
                    cols = gaba_hex)
    plot_list[["gaba"]] <- plot
  }
  
  if (groupcol == TRUE) {
    plot <- DimPlot(seur_obj, 
                    group.by = groupby, 
                    split.by = splitby, 
                    cols = grp_hex)
    plot_list[["group"]] <- plot
  }
  
  if (sexcol == TRUE) {
    plot <- DimPlot(seur_obj, 
                    group.by = groupby, 
                    split.by = splitby, 
                    cols = sex_hex)
    plot_list[["sex"]] <- plot
  }
  
  # Customize and print the modified plots
  for (i in names(plot_list)) {
    modified_plot <- plot_list[[i]] +
      theme_void() + 
      theme(legend.position = "none") + 
      ggtitle(NULL)
    print(modified_plot)
  }
  
  # Close the PDF file
  dev.off()
}



