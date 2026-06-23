#' Assign experimental metadata to a Seurat object via a lookup table
#'
#' Replaces the previous hardcoded index-based assignment. Metadata is joined
#' from a CSV with columns: cart, sample_tag, ratID, sex, experience, group.
#' Matching is done on orig.ident (cart) and Sample_tag simultaneously,
#' so the same sample tag label used across cartridges is handled correctly.
#'
#' @param seur_obj    a Seurat object produced by create_seur()
#' @param meta_lookup data frame read from the metadata CSV
#'
#' @return the Seurat object with ratID, sex, experience, and group added
#'         to metadata; cells with no match (Multiplet, Undetermined, etc.)
#'         will have NA in these columns
#' @export

assign_combine <- function(seur_obj, meta_lookup) {
  
  meta <- seur_obj@meta.data
  
  # Join on both cart (orig.ident) and sample_tag so the same ST label
  # used across cartridges maps to the correct individual
  meta <- merge(
    meta,
    meta_lookup[, c("cart", "sample_tag", "ratID", "sex", "experience", "group")],
    by.x = c("orig.ident", "Sample_tag"),
    by.y = c("cart",        "sample_tag"),
    all.x = TRUE,    # keep all cells; non-matching (Multiplet etc.) get NA
    sort  = FALSE
  )
  
  # merge() reorders rows — restore original cell order
  rownames(meta) <- meta$row_names_col  # placeholder; see below
  
  # Safer: preserve rownames through merge via explicit column
  meta_orig           <- seur_obj@meta.data
  meta_orig$cell_id   <- rownames(meta_orig)
  
  meta_joined <- merge(
    meta_orig,
    meta_lookup[, c("cart", "sample_tag", "ratID", "sex", "experience", "group")],
    by.x  = c("orig.ident", "Sample_tag"),
    by.y  = c("cart",        "sample_tag"),
    all.x = TRUE,
    sort  = FALSE
  )
  
  # Restore original cell order
  meta_joined         <- meta_joined[match(meta_orig$cell_id, meta_joined$cell_id), ]
  rownames(meta_joined) <- meta_joined$cell_id
  meta_joined$cell_id <- NULL
  
  seur_obj@meta.data <- meta_joined
  
  return(seur_obj)
}