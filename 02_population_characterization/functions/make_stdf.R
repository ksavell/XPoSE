#' Make Sample Tag read dataframe
#'
#' @param seur_obj seurat object with sample tag reads as metadata
#'
#' @return dataframe of sample tag reads and sample tag call
#' @export
#'
#' @examples
make_stdf <- function(seur_obj) {
  IDs <- seur_obj@meta.data[, c(
    "capture",
    "ratID",
    "xpose_tag",
    "xpose_tag_02_reads",
    # "SampleTag03_reads",
    "xpose_tag_04_reads",
    # "SampleTag05_reads",
    "xpose_tag_06_reads",
    # "SampleTag07_reads",
    "xpose_tag_08_reads"
    # "SampleTag09_reads"
  )]
  
  # colnames(IDs) <- c(
  #   "cart", "rat", "st",
  #   "st2", "st3", "st4", "st5", "st6", "st7", "st8", "st9"
  # )
  
  colnames(IDs) <- c(
    "cart", "rat", "st",
    "st2", "st4", "st6", "st8"
  )
  
  return(IDs)
}
