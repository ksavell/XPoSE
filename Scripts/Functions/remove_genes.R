#' Remove genes from seurat object
#'
#' @param seur_obj seurat object
#'
#' @return seurat object with genes removed
#' @export
#'
#' @examples
remove_genes <- function(seur_obj) {
        genes <- seur_obj@assays$RNA@counts@Dimnames[[1]]
        
        # search known pseudogene, microRNA, model gene, or other non-coding patterns
        indices_remove <- c(grep(pattern = "BC0", x = genes),
                            grep(pattern = ".ps", x = genes),
                            grep(pattern = "Rik", x = genes),
                            grep(pattern = "rik", x = genes),
                            grep(pattern = "Gm[0-9]", x = genes),
                            grep(pattern = "Mir[0-9]", x = genes),
                            grep(pattern = "LOC[0-9]", x = genes))
        
        # remove genes
        genes[-indices_remove] -> genes
        
        #subset here again to save the day :') 
        seur_obj <- subset(seur_obj, features = genes)
        return(seur_obj)
}