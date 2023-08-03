#' Make Sample Tag read dataframe
#'
#' @param seur_obj seurat object with sample tag reads as metadata
#'
#' @return dataframe of sample tag reads and sample tag call
#' @export
#'
#' @examples
make_stdf <- function(seur_obj) {
        IDs <- data.frame(seur_obj@assays[["RNA"]]@data@Dimnames[[2]],
                          cart = seur_obj@meta.data[["orig.ident"]])
        
        IDs$cart <- seur_obj@meta.data[["orig.ident"]]
        IDs$rat <- seur_obj@meta.data[["ratID"]]
        IDs$st <- seur_obj@meta.data[["Sample_tag"]]
        IDs$st2 <- seur_obj@meta.data[["SampleTag02_reads"]]
        IDs$st3 <- seur_obj@meta.data[["SampleTag03_reads"]]
        IDs$st4 <- seur_obj@meta.data[["SampleTag04_reads"]]
        IDs$st5 <- seur_obj@meta.data[["SampleTag05_reads"]]
        IDs$st6 <- seur_obj@meta.data[["SampleTag06_reads"]]
        IDs$st7 <- seur_obj@meta.data[["SampleTag07_reads"]]
        IDs$st8 <- seur_obj@meta.data[["SampleTag08_reads"]]
        IDs$st9 <- seur_obj@meta.data[["SampleTag09_reads"]]
        
        rownames(IDs) <- IDs$seur_obj.assays...RNA....data.Dimnames..2..
        #IDs <- IDs %>% select(cart, rat, st)
        
        return(IDs)
}