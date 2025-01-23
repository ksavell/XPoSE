#' Downsample 2 groups in 1:1 ratio by cluster
#'
#' @param seur_obj seurat object
#' @param group_to_subset group one, usually "Active"
#' @param group_to_blend group 2, usually "Non-active"
#' @param cluster cluster, usually used with a loop
#' @param number number of nuclei to pull for each group
#' @param seed reproducibility seed
#'
#' @returns
#' @export
#'
#' @examples
#' 
group_downsample_counts <- function(seur_obj, group_to_subset, group_to_blend, cluster, number = 1, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Extract metadata
  metadata <- seur_obj@meta.data
  
  # Check if the 'group' column exists
  if (!"group" %in% colnames(metadata) || !"cluster_name" %in% colnames(metadata)) {
    stop("The metadata must contain 'group' and 'cluster_name' columns.")
  }
  
  # Filter metadata for each group
  metadata_subset <- metadata %>%
    dplyr::filter(group == group_to_subset, cluster_name == cluster)
  metadata_blend <- metadata %>%
    dplyr::filter(group == group_to_blend, cluster_name == cluster)

  # Number of cells in each group
  num_blend <- as.numeric(nrow(metadata_blend))
  num_subset <- as.numeric(nrow(metadata_subset))
 
  # Downsample `group_to_subset` if it has more cells than desired
  if (num_subset >= number) {
    subset_cells <- metadata_subset %>%
      dplyr::slice_sample(n = number) %>%
      rownames()
  } else {
    subset_cells <- NULL
  }
  
  # Downsample `group_to_blend` if it has more cells than desired
  if (num_blend >= number) {
    blend_cells <- metadata_blend %>%
      dplyr::slice_sample(n = number) %>%
      rownames()
  } else {
    blend_cells <- NULL
  }

  # Combine all desired cells
  all_cells <- c(subset_cells, blend_cells)
  
  # Return the indices of the selected cells
  return(all_cells)
  }


