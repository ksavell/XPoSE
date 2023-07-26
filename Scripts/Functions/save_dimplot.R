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
save_dimplot <- function(seur_obj, groupby = NULL, file_n = NULL){
        pdf(file = file_n,
            width = 10,
            height = 10)
        print(DimPlot(seur_obj, group.by = groupby) + 
                theme_void() + 
                theme(legend.position = "none") + 
                ggtitle(NULL))
        dev.off()   
}
