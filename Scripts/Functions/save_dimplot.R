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
      
        glut_hex <- c("#64C7C8",'#41B75F','#2C8CB9','#0A5B8C',
                       '#3C9E64','#6F499D')
        
        gaba_hex <- c('#E66027','#F8991D',
                       '#C03C82','#A669AB',
                        '#C52126','#DB808C',
                         '#B0B235')
        
        grp_hex <- c('#e6e6e6','gray60','#ae1e5b')
        
        sex_hex <- c("#2C5F2D","#97BC62FF")
        
        if (glutcol == TRUE) {
                pdf(file = file_n,
                    width = 10*splitl,
                    height = 10)
                print(DimPlot(seur_obj, 
                              group.by = groupby, 
                              split.by = splitby, 
                              cols = glut_hex) + 
                        theme_void() + 
                        theme(legend.position = "none") + 
                        ggtitle(NULL))
                dev.off()
        }
        
        if (gabacol == TRUE) {
                pdf(file = file_n,
                    width = 10*splitl,
                    height = 10)
                print(DimPlot(seur_obj, 
                              group.by = groupby, 
                              split.by = splitby, 
                              cols = gaba_hex) + 
                        theme_void() + 
                        theme(legend.position = "none") + 
                        ggtitle(NULL))
                dev.off()
        }
        if (groupcol == TRUE) {
                pdf(file = file_n,
                    width = 10*splitl,
                    height = 10)
                print(DimPlot(seur_obj, 
                              group.by = groupby, 
                              split.by = splitby, 
                              cols = grp_hex) + 
                        theme_void() + 
                        theme(legend.position = "none") + 
                        ggtitle(NULL))
                dev.off()
        }
        if (sexcol == TRUE) {
          pdf(file = file_n,
              width = 10*splitl,
              height = 10)
          print(DimPlot(seur_obj, 
                        group.by = groupby, 
                        split.by = splitby, 
                        cols = sex_hex) + 
                  theme_void() + 
                  theme(legend.position = "none") + 
                  ggtitle(NULL))
          dev.off()
        }
}
