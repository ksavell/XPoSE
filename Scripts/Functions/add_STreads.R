#' Adds sample tag reads to metadata of Seurat object
#'
#' @param seurat_obj input Seurat object
#' @param reads_obj input reads Table
#' @param STrange input range of Sample tags
#'
#' @return Seurat object with assigned reads
#' @export
add_STreads <- function(seurat_obj, reads_obj, STrange = c(1:12)) {
  
  for (i in STrange) {
    colName <- paste0("stag", i, "_reads")
    seurat_obj[[colName]] <- cbind(reads_obj[[paste0("SampleTag", sprintf("%02d", i), "_mm.stAbO")]])
  }
  return(seurat_obj)
}