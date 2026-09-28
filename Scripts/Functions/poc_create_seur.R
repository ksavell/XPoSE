#' Create Seurat Object and add count and sample tag tables
#'
#' @param count       count table for this cartridge (data frame, cells x genes)
#' @param cart_label  meaningful cartridge label from the metadata lookup table;
#'                    used as orig.ident (e.g. "Cart3_IL")
#' @param tag         sample tag calls table for this cartridge
#' @param addSTreads  logical; whether to add per-tag read counts as metadata
#' @param streads     sample tag reads table (required if addSTreads = TRUE)
#' @param STrange     integer vector of sample tag indices to add (default 1:12)
#'
#' @return a Seurat object with count data, Sample_tag, percent_mt,
#'         and optionally per-tag read count columns in metadata
#' @export
#'
#' @examples
#' create_seur(count = counts_c1, cart_label = "Cart3_IL",
#'             tag = tags_c1, addSTreads = TRUE,
#'             streads = streads_c1, STrange = 1:12)

create_seur <- function(count, cart_label, tag,
                        addSTreads = FALSE, streads = NULL, STrange = c(1:12)) {
  
  # BD output is cells x genes; Seurat expects genes x cells — transpose
  seur_obj <- CreateSeuratObject(
    counts  = t(as.data.frame(count)),
    project = cart_label          # orig.ident set from lookup table cart column
  )
  
  # Add sample tag calls and mitochondrial percentage
  seur_obj[["Sample_tag"]]  <- tag[, 1, drop = FALSE]
  seur_obj[["percent_mt"]]  <- PercentageFeatureSet(seur_obj, pattern = "^mt")
  
  # Optionally add per-sample-tag read counts as metadata columns
  if (addSTreads) {
    for (i in STrange) {
      col_name  <- ifelse(i < 10,
                          paste0("SampleTag0", i, "_reads"),
                          paste0("SampleTag",  i, "_reads"))
      reads_col <- paste0("SampleTag", sprintf("%02d", i), "_mm.stAbO")
      seur_obj[[col_name]] <- streads[[reads_col]]
    }
  }
  
  return(seur_obj)
}