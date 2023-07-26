#' Plot DimPlot with no formatting for F1
#'
#' @param seur_obj seurat object
#' @param groupby factor to group the seurat object
#' @param savename name to save the pdf
#'
#' @return 
#' @export
#'
#' @examples
save_dimplot <- function(seur_obj, groupby = NULL, savename = date()){
        pdf(file = savename,
            width = 10,
            height = 10)
        DimPlot(seur_obj, group.by = groupby) + 
                theme_void() + 
                theme(legend.position = "none") + 
                ggtitle(NULL)
        dev.off()   
}
