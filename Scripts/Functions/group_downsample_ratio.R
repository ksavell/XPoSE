
# Define the down sampling function with a dynamic approach for downsampling
#' Title
#'
#' @param seur_obj 
#' @param group_to_subset 
#' @param group_to_blend 
#' @param ratio 
#' @param seed 
#'
#' @return
#' @export
#'
#' @examples
group_downsample_ratio <- function(seur_obj, group_to_subset, group_to_blend, ratio = 1, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Extract metadata
  metadata <- seur_obj@meta.data
  
  # Check if the 'group' column exists
  if (!"group" %in% colnames(metadata)) {
    stop("The 'group' column does not exist in the metadata.")
  }
  
  # Filter metadata for each group
  metadata_subset <- metadata %>% filter(group == group_to_subset)
  metadata_blend <- metadata %>% filter(group == group_to_blend)
  
  # Calculate the desired number of cells in each group based on the ratio
  num_blend <- nrow(metadata_blend)
  num_subset <- nrow(metadata_subset)
  desired_num_subset <- round(num_blend * ratio)
  desired_num_blend <- round(num_subset / ratio)
  
  # Downsample `group_to_subset` if it has more cells than desired
  if (num_subset > desired_num_subset) {
    subset_cells <- metadata_subset %>%
      slice_sample(n = desired_num_subset) %>%
      rownames()
  } else {
    subset_cells <- rownames(metadata_subset)
  }
  
  # Downsample `group_to_blend` if it has more cells than desired
  if (num_blend > desired_num_blend) {
    blend_cells <- metadata_blend %>%
      slice_sample(n = desired_num_blend) %>%
      rownames()
  } else {
    blend_cells <- rownames(metadata_blend)
  }
  
  # Select homecage cells if needed
  cells_homecage <- metadata %>%
    filter(experience == "HC") %>%
    rownames()
  
  # Combine all desired cells
  all_cells <- c(subset_cells, blend_cells, cells_homecage)
  
  # Return the indices of the selected cells
  return(all_cells)
}
