#' Plot VlnPlot with correct formatting for Figures
#'
#' @param seur_obj seurat object
#' @param splitby factor to split plot
#' @param splitl number of splits so the width is calculated correctly #do you need this
#' @param file_n name the file, must include .PDF
#' @param glutcol logical, T for glut cluster colors
#' @param gabacol logical, T for gaba cluster colors
#' @param groupcol logical, T for group cluster colors
#' @param groupby factor to group the plot
#'
#' @return saves PDF with specific colors and split/group settings
#' @export
#'
#' @examples
save_vlnplot <- function(seur_obj, groupby = NULL, splitby = NULL, splitl = 1, 
                         file_n = NULL, glutcol = FALSE,
                         gabacol = FALSE, groupcol = FALSE,
                         sexcol = FALSE, feature = "nCount_RNA"){
      
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
                    width = 2*splitl,
                    height = 2)
                print(VlnPlot(seur_obj, 
                              group.by = groupby, 
                              split.by = splitby, 
                              cols = glut_hex,
                              features = feature,
                              pt.size = 0, 
                              ncol = 1,
                              flip = T) + 
                        coord_flip() + 
                        theme(legend.position = "none") + 
                        theme(axis.title = element_blank()) +
                        ggtitle(NULL))
                dev.off()
        }
        
        if (gabacol == TRUE) {
                pdf(file = file_n,
                    width = 2*splitl,
                    height = 2)
                print(VlnPlot(seur_obj, 
                              group.by = groupby, 
                              split.by = splitby, 
                              cols = gaba_hex,
                              features = feature) + 
                        theme_void() + 
                        theme(legend.position = "none") + 
                        ggtitle(NULL))
                dev.off()
        }
        if (groupcol == TRUE) {
                pdf(file = file_n,
                    width = 2*splitl,
                    height = 2)
                print(VlnPlot(seur_obj, 
                              group.by = groupby, 
                              split.by = splitby, 
                              cols = group_hex,
                              features = feature) + 
                        theme_void() + 
                        theme(legend.position = "none") + 
                        ggtitle(NULL))
                dev.off()
        }
        if (sexcol == TRUE) {
          pdf(file = file_n,
              width = 2*splitl,
              height = 2)
          print(VlnPlot(seur_obj, 
                        group.by = groupby, 
                        split.by = splitby, 
                        cols = sex_hex,
                        features = feature) + 
                  theme_void() + 
                  theme(legend.position = "none") + 
                  ggtitle(NULL))
          dev.off()
        }
}
