# subset specific group population in seurat object

#' Title
#'
#' @param seur_obj seurat object
#' @param group_to_subset group that you want to downsample
#' @param frac the proportion that you want to keep of the group to subset
#'
#' @return
#' @export
#'
#' @examples
library(Seurat)
library(dplyr)

# Define the down sampling function
group_downsample <- function(seur_obj, group_to_subset, frac = 1, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Extract metadata
  metadata <- seur_obj@meta.data
  
  # Filter metadata for the specific group to subset
  metadata_group <- metadata %>%
    filter(group == group_to_subset)
  
  # Downsample the cells in the specified group
  cells_subset <- metadata_group %>%
    slice_sample(prop = frac) %>%
    rownames()
  
  # Get cells from other groups
  other_groups_cells <- metadata %>%
    filter(group != group_to_subset) %>%
    rownames()
  
  # Combine all desired cells
  all_cells <- c(cells_subset, other_groups_cells)
  
  # Return the indices of the cells chosen
  return(all_cells)
}

