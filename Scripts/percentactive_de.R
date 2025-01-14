## percent active in NC

# Load necessary libraries
library(Seurat)
library(dplyr)
library(readr)

# Define key variables
iterations <- 100
seur_obj <- all  # Full Seurat object
seur_obj$experience <- ifelse(seur_obj$group == "Homecage", "HC", "NC")
clusters <- unique(seur_obj$cluster_name)
target_percentages <- c(20, 40, 60, 80)  # Adjust percentages as needed

# Loop through clusters
for (cl in clusters) {
  
  # Create a folder for the cluster
  dir_path <- paste0(cl, "/")
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  # Loop through different target percentages
  for (target_percentage in target_percentages) {
    
    all_indices <- list()
    all_seeds <- numeric(iterations)
    
    for (it in 1:iterations) {
      seed <- sample(1:10000, 1)  # Use a new seed for each iteration
      all_seeds[it] <- seed
      
      # Get downsampled cells to match target percentage
      chosen_cells <- group_downsample_percentage(seur_obj,
                                                  group_active = "Active",
                                                  group_nonactive = "Non-active",
                                                  target_percentage = target_percentage,
                                                  seed = seed)
      all_indices[[it]] <- chosen_cells
    }
    
    # Save indices and seeds for this percentage
    save(all_indices, all_seeds, file = paste0(dir_path, "downsampled_indices_and_seeds_", target_percentage, ".RData"))
    
    # Run DESeq2 for each downsampled iteration
    results <- list()
    
    for (j in 1:iterations) {
      chosen_cells <- all_indices[[j]]
      seurat_subset <- subset(seur_obj, cells = chosen_cells) 
      
      deseq2_results <- single_factor_DESeq(object = seurat_subset,
                                            comp_vect = c("experience", "NC", "HC"),
                                            cluster = cl,
                                            min_cell = 1)
      results[[j]] <- deseq2_results
    }
    
    # Save DESeq2 results for this percentage
    save(results, file = paste0(dir_path, "results_", target_percentage, ".RData"))
    
    # Optional: Call a tally function if needed
    # tally_iterations(results, cl, prop_frac = target_percentage)
  }
}

# need to tally results

for (cl in clusters) {
  clstr <- cl 
  
  files <- list.files(clstr, pattern = "^results_.*\\.RData$", full.names = TRUE)
  data_list <- lapply(files, function(file) {
    env <- new.env()
    load(file, envir = env)
    as.list(env)
  })
  
  names(data_list) <- basename(files)
  
  # remove the file extension from the names
  names(data_list) <- sub("\\.RData$", "", names(data_list))
  
  # now tally the results and save that csv 
  for (j in names(data_list)) {
    merged_results_df <- tally_iterations(data_list, j)
    
    # Check the result
    print(paste("Processing:", j))
    print(head(merged_results_df))
    
    file_path <- paste0(clstr, "/", j, "_iteration_tally.csv")
    print(paste("Saving to:", file_path))
    
    write.csv(merged_results_df, file = file_path, row.names = TRUE)
  }
  
  # I am not sure if you even need this part
  
  all_data_lists <- list()
  
  for (cl in clusters) {
    
    # Folder path is the same as cluster name in this case
    folder_id <- cl
    
    # Get a list of all CSV files in the folder
    csv_files <- list.files(path = cl, pattern = "results_.*_iteration_tally\\.csv", full.names = TRUE)
    
    # Initialize a list to store data frames for this specific folder
    data_list <- list()
    
    # Loop through each CSV file in the folder
    for (file in csv_files) {
      # Extract the unique number from the filename
      unique_number <- sub(".*results_(.*?)_iteration_tally\\.csv", "\\1", basename(file))
      
      # Create a unique name by combining the folder name and unique number
      data_name <- paste0(folder_id, "_", unique_number)
      
      # Read the CSV file into a data frame and store it in data_list
      data_list[[data_name]] <- read.csv(file)
    }
    
    # Store the data_list in all_data_lists with the cluster ID as the key
    all_data_lists[[folder_id]] <- data_list
  }
  
  # Check the structure to confirm
  print(names(all_data_lists))
  
  
  
  ## this part is important
  # Define clusters and numbers
  clusters <- names(all_data_lists)
  numbers <- c("0.035", "0.07", "0.105", "0.175", "0.245")
  
  # Initialize data frames to store the results: one for averages, one for SD, and one for nods
  average_sums_df <- data.frame(matrix(nrow = length(clusters), ncol = length(numbers)))
  sd_sums_df <- data.frame(matrix(nrow = length(clusters), ncol = length(numbers)))
  
  rownames(average_sums_df) <- clusters
  colnames(average_sums_df) <- numbers
  rownames(sd_sums_df) <- clusters
  colnames(sd_sums_df) <- numbers
  
  # Populate the data frames with mean, SD values, and nods data
  for (cluster in clusters) {
    for (num in numbers) {
      # Construct the name of the DataFrame in all_data_lists
      df_name <- paste0(cluster, "_", num)
      
      # Check if the DataFrame exists and compute the average and SD if it does
      if (df_name %in% names(all_data_lists[[cluster]])) {
        df <- all_data_lists[[cluster]][[df_name]]
        numeric_df <- df[, sapply(df, function(col) is.numeric(col) || is.integer(col))]
        
        # Calculate column sums, then compute mean and SD
        col_sums <- colSums(numeric_df, na.rm = TRUE)
        average_sums_df[cluster, num] <- mean(col_sums)
        sd_sums_df[cluster, num] <- sd(col_sums)  # Calculate SD instead of SEM
      } else {
        # Assign NA if the DataFrame is missing
        average_sums_df[cluster, num] <- NA
        sd_sums_df[cluster, num] <- NA
      }
    }
  }
  
  # Check the resulting data frames
  print(average_sums_df)
  print(sd_sums_df)
  
  # save as csv
  write.csv(average_sums_df, "iteration_plot_sum_average.csv")
  write.csv(sd_sums_df, "iteration_plot_stdev.csv")