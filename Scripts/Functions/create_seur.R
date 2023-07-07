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
#' create_seur(counts.C1, 1, tags.C1)
create_seur <- function(count, cart_no, tag){
  
  # Create Seurat objects for each cartridge while transposing since BD output is opposite of Seurat input
  seur_obj <- CreateSeuratObject(counts =  t(as.data.frame(counts[cart_no])), project = paste("C",cart_no,sep=""))
  
  # Add Sample_tags as metadata to the Seurat object
  seur_obj[["Sample_tag"]] <- cbind(tag)
  seur_obj[["percent_mt"]] <- PercentageFeatureSet(seur_obj, pattern = "^mt.")
  return(seur_obj)
}

# creates separate seurat object for each count table present in environment
for(i in 1:length(counts)){
  assign(paste("data",".C",i,sep=""),
         create_seur(counts[i], i, tags[i]),
         envir = .GlobalEnv)
  
}