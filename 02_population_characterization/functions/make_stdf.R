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
    "orig.ident",
    "ratID",
    "Sample_tag",
    "SampleTag02_reads",
    # "SampleTag03_reads",
    "SampleTag04_reads",
    # "SampleTag05_reads",
    "SampleTag06_reads",
    # "SampleTag07_reads",
    "SampleTag08_reads"
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