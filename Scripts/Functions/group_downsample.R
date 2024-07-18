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
group_downsample <- function(seur_obj, group_to_subset, group_to_blend, frac = 1, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Extract metadata
  metadata <- seur_obj@meta.data
  
  # Check if the group column exists
  if (!"group" %in% colnames(metadata)) {
    stop("The 'group' column does not exist in the metadata.")
  }
  
  # Filter metadata for the specific group to subset
  metadata_group <- metadata %>%
    filter(group == group_to_subset)
  
  # Get cells from blend groups
  other_groups_cells <- metadata %>%
    filter(group == group_to_blend) %>%
    rownames()
  
  # fraction adjusted based on blend group
  num_group <- nrow(metadata_group)
  num_blend <- nrow(metadata %>%
                      filter(group == group_to_blend))
  
  # Calculate adjusted number of group cells to downsample
  num_group_adj <- frac * num_blend / (1 - frac)
  frac_adj <- min(1, num_group_adj / num_group)  # Ensure the fraction does not exceed 1
  
  # Downsample the cells in the specified group
  cells_subset <- metadata_group %>%
    slice_sample(prop = frac_adj) %>%
    rownames()
  
  # grab homecage cells
  cells_homecage <- metadata %>%
    filter(experience == "HC") %>%
    rownames()
  
  # Combine all desired cells
  all_cells <- c(cells_subset, other_groups_cells, cells_homecage)
  
  # Return the indices of the cells chosen
  return(all_cells)
}


