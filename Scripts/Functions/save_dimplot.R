#' Plot DimPlot with no formatting for Figures
#'
#' @param seur_obj seurat object
#' @param groupby factor to group the seurat object
#' @param file_name name to save the pdf
#'
#' @return 
#' @export
#'
#' @examples
save_dimplot <- function(seur_obj, groupby = NULL, file_n = NULL, glutcol = FALSE,
                         gabacol = FALSE, groupcol = FALSE){
      
        glut_hex <- c("#64C7C8",'#41B75F','#2C8CB9','#0A5B8C',
                       '#3C9E64','#6F499D')
        
        gaba_hex <- c('#E66027','#F8991D',
                       '#C03C82','#A669AB',
                        '#C52126','#DB808C',
                         '#B0B235')
        
        grp_hex <- c('#C7C7C7','gray2','#D0397C')
        
        if (glutcol == TRUE) {
                pdf(file = file_n,
                    width = 10,
                    height = 10)
                print(DimPlot(seur_obj, group.by = groupby, cols = glut_hex) + 
                              theme_void() + 
                              theme(legend.position = "none") + 
                              ggtitle(NULL))
                dev.off()
        }
        
        if (gabacol == TRUE) {
                pdf(file = file_n,
                    width = 10,
                    height = 10)
                print(DimPlot(seur_obj, group.by = groupby, cols = gaba_hex) + 
                              theme_void() + 
                              theme(legend.position = "none") + 
                              ggtitle(NULL))
                dev.off()
        }
        if (groupcol == TRUE) {
                pdf(file = file_n,
                    width = 10,
                    height = 10)
                print(DimPlot(seur_obj, group.by = groupby, cols = grp_hex,
                              shuffle = T) + 
                              theme_void() + 
                              theme(legend.position = "none") + 
                              ggtitle(NULL))
                dev.off()
        }
        
}
