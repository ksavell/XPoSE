#' Create Seurat Object and add count and sample tag tables
#'
#' @param count inputs count table for cartridge
#' @param cart_no inputs cartridge number
#' @param tag inputs sample tag table for cartridge
#'
#' @return a seurat object containing count and sample tag tables
#' @export
#'
#' @examples
#' create_seur(counts_c1, 1, tags_c1)
create_seur <- function(count, cart_no, tag,
                        addSTreads = FALSE, streads = NULL, STrange = c(2:9)){
  
  # Create Seurat objects for each cartridge while transposing since BD output is opposite of Seurat input
  seur_obj <- CreateSeuratObject(counts =  t(as.data.frame(counts[[cart_no]])), project = paste("C",cart_no,sep=""))
  
  # Add Sample_tags as metadata to the Seurat object
  seur_obj[["Sample_tag"]] <- cbind(tag[1])
  # seur_obj$stag <- cbind(tag$Sample_tag)
  seur_obj[["percent_mt"]] <- PercentageFeatureSet(seur_obj, pattern = "^mt")
  if (addSTreads == TRUE) {
    reads_obj <- streads
    for (i in STrange) {
      colName <- paste("SampleTag0", i, "_reads", sep="")
      seur_obj[[colName]] <- NA 
      seur_obj[[colName]] <- cbind(reads_obj[[paste("SampleTag", sprintf("%02d", i), "_mm.stAbO", sep="")]])
    }
    return(seur_obj)
  }
  return(seur_obj)
}
