
# This is an edited version of the function "percentactive_de.R"
# Fixed the formula to calculate the number of actives and nonactives 
# The formula used to solve for A and NA is A/(NA + A) = desired percentage of active cells 
# The formula is then A/(NA + A) = p/100, where p is the desired %
# Therefore, solving for A gives us --- A = (NA(P))/(100-p)
# Then, solving for NA gives us NA = (A(100-p))/(p)
# p is already known!


group_downsample <- function(seur_obj, group_to_subset, group_to_blend, percent = 1, seed = NULL) {
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
  metadata_subset <- metadata %>% dplyr::filter(group == group_to_subset)
  metadata_blend <- metadata %>% dplyr::filter(group == group_to_blend)
  metadata_hc <- metadata %>% dplyr::filter(group == "Homecage")
  
  # Number of cells in each group
  num_blend <- as.numeric(nrow(metadata_blend))
  num_subset <- as.numeric(nrow(metadata_subset))
  num_hc <-  as.numeric(nrow(metadata_hc))
  
  # Calculate the desired number of Active cells (A) based on the formula
  if (percent == 0) {
    desired_num_subset <- 0  # If 0% active, use no active cells
  } else if (percent == 100) {
    desired_num_subset <- num_subset  # If 0% active, use no active cells
  } else {
    desired_num_subset <- round((percent / 100) * num_subset) ## not working?
  }
  
  # Calculate the desired number of Non-active cells (NA)
  if (percent == 100) {
    desired_num_blend <- 0  # If 100% active, no non-active cells needed
  } else if (percent == 0) {
      desired_num_blend <- num_blend   
      } else {
    desired_num_blend <- round(((desired_num_subset) * (100 - percent)) / (percent)) ## not working?
  }
  
  # Calculate the number of HC cells to keep (match the number of NC)
  desired_num_hc <- desired_num_blend + desired_num_subset
  
  # Downsample `group_to_subset` if it has more cells than desired
  if (num_subset >= desired_num_subset) {
    subset_cells <- metadata_subset %>%
      dplyr::slice_sample(n = desired_num_subset) %>%
      rownames()
  } else {
    subset_cells <- rownames(metadata_subset)
  }
  
  # Downsample `group_to_blend` if it has more cells than desired
  if (num_blend >= desired_num_blend) {
    blend_cells <- metadata_blend %>%
      dplyr::slice_sample(n = desired_num_blend) %>%
      rownames()
  } else {
    blend_cells <- rownames(metadata_blend)
  }
  
  # Select homecage cells
  if (num_hc >= desired_num_hc) {
    hc_cells <- metadata_hc %>%
      dplyr::slice_sample(n = desired_num_hc) %>%
      rownames()
  } else {
    hc_cells <- rownames(metadata_hc)
  }

  # Combine all desired cells
  all_cells <- c(subset_cells, blend_cells, hc_cells)
  
  # Return the indices of the selected cells
  return(all_cells)
  }


