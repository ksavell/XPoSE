# subset specific group population in seurat object

#' Title
#'
#' @param seur_obj seurat object
#' @param group_to_subset group that you want to downsample
#' @param frac the proportion that you want to keep of the group to subset
#' @param bio_rep metadata associated with biological replicate identity
#'
#' @return
#' @export
#'
#' @examples
group_downsample <- function(seur_obj, group_to_subset, frac = 1, bio_rep = "ratID") {
  
  # Extract metadata
  metadata <- seur_obj@meta.data
  
  # Filter metadata for the specific group to subset
  metadata_group <- metadata %>%
    filter(group == group_to_subset)
  
  # Get unique rats within the group to subset
  ids <- unique(metadata_group[[bio_rep]])
  
  # Initialize list to store subsets
  subset_list <- list()
  
  # Loop through each rat
  for (rep_id in ids) {
    # Filter cells for the current rat within the specific group
    cells_subset <- metadata_group %>%
      filter(.data[[bio_rep]] == rep_id) %>%
      slice_sample(prop = frac) %>%
      rownames()
    
    # Add subset to list
    subset_list[[rep_id]] <- cells_subset
  }
  
  # Combine all subsets for the specific group
  cells_to_subset <- unlist(subset_list)
  
  # Get cells from other groups
  other_groups_cells <- metadata %>%
    filter(group != group_to_subset) %>%
    rownames()
  
  # Combine all desired cells
  all_cells <- c(cells_to_subset, other_groups_cells)
  
  # Subset the Seurat object
  seurat_subset <- subset(seur_obj, cells = all_cells)
  
  return(seurat_subset)
}

