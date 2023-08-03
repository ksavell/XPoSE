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
#' @param plotw plot width
#' @param ploth plot height
#'
#' @return saves PDF with specific colors and split/group settings
#' @export
#'
#' @examples
save_vlnplot <- function(seur_obj, groupby = NULL, splitby = NULL, splitl = 1, 
                         file_n = NULL, glutcol = FALSE,
                         gabacol = FALSE, groupcol = FALSE,
                         sexcol = FALSE, 
                         feature = "nCount_RNA", vln_max = NULL,
                         plotw = 2, ploth = 2){
      
        glut_hex <- c("#64C7C8",'#41B75F','#2C8CB9','#0A5B8C',
                       '#3C9E64','#6F499D')
        
        gaba_hex <- c('#E66027','#F8991D',
                       '#C03C82','#A669AB',
                        '#C52126','#DB808C',
                         '#B0B235')
        
        grp_hex <- c('#e6e6e6','gray60','#ae1e5b')
        
        sex_hex <- c("#2C5F2D","#97BC62FF")
        
        if (glutcol == TRUE) {
          data <- VlnPlot(seur_obj, 
                          group.by = groupby, 
                          split.by = splitby, 
                          cols = glut_hex,
                          features = feature,
                          pt.size = 0, 
                          y.max = vln_max) +
            coord_flip() + 
            #theme_void() +
            theme(legend.position = "none",
                  axis.title = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank()) +
            ggtitle(NULL)
        }
        
        if (gabacol == TRUE) {
          data <- VlnPlot(seur_obj, 
                          group.by = groupby, 
                          split.by = splitby, 
                          cols = gaba_hex,
                          features = feature,
                          pt.size = 0, 
                          y.max = vln_max) +
            coord_flip() + 
            #theme_void() +
            theme(legend.position = "none",
                  axis.title = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank()) +
            ggtitle(NULL)
        }
        
        if (groupcol == TRUE) {
          data <- VlnPlot(seur_obj, 
                          group.by = groupby, 
                          split.by = splitby, 
                          cols = grp_hex,
                          features = feature,
                          pt.size = 0, 
                          y.max = vln_max) +
            coord_flip() + 
            #theme_void() +
            theme(legend.position = "none",
                  axis.title = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank()) +
            ggtitle(NULL)
        }
        
        if (sexcol == TRUE) {
          data <- VlnPlot(seur_obj, 
                          group.by = groupby, 
                          split.by = splitby, 
                          cols = sex_hex,
                          features = feature,
                          pt.size = 0, 
                          y.max = vln_max) +
            coord_flip() + 
            #theme_void() +
            theme(legend.position = "none",
                  axis.title = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank()) +
            ggtitle(NULL)
        }
        
          # now save the pdf
          pdf(file = file_n,
              width = plotw,
              height = ploth)
          print(data)
          dev.off()
        
}