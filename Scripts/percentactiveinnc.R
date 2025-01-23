# Load necessary libraries
library(Seurat)
library(dplyr)
library(readr)


# Load the Seurat object from an RData file
load("/Users/holmesar/Library/CloudStorage/Box-Box/RM_Projects/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/NeuroResource_January2025/DataForFigures/Robjects/all_10312024.RData")
source("Scripts/Functions/single_factor_DESeq.R")
source("Scripts/Functions/group_downsample.R")

# Define key variables
iterations <- 100
seur_obj <- all
seur_obj$experience <- ifelse(seur_obj$group == "Homecage", "HC", "NC")
clusters <- unique(seur_obj$cluster_name)
target_percentages <- c(0, 5, 10, 25, 50, 75, 100)  # Adjust percentages as needed


# Loop 1: define percent subsetting for 100 iterations --------------------


# First Loop: Calculate nuclei to keep and output RData file


# Loop through different target percentages
for (target_percentage in target_percentages) {
  
  # Initialize storage for this percentage
  all_indices <- list()
  all_seeds <- numeric(iterations)
  
  for (it in 1:iterations) {
    seed <- sample(1:10000, 1)  # Use a new seed for each iteration
    all_seeds[it] <- seed
    
    # Call group_downsample_edited
    chosen_cells <- group_downsample(
      seur_obj = seur_obj,
      group_to_subset = "Active",
      group_to_blend = "Non-active",
      percent = target_percentage,
      seed = seed
    )
    
    # Store chosen cells in the list
    all_indices[[it]] <- chosen_cells
  }
  
  # Save indices and seeds for this percentage in the working directory
  save(all_indices, all_seeds, file = paste0("downsampled_indices_and_seeds_", target_percentage, ".RData"))

}



# Loop 2: Run DESeq2 for each percent/iteration combo ---------------------


# Second Loop: Load downsampled indices and seeds, process results for each percentage, cluster, and iteration


# Loop through target percentages
for (target_percentage in target_percentages) {
  
  # Load the downsampled indices and seeds for the target percentage
  load(paste0("downsampled_indices_and_seeds_", target_percentage, ".RData"))
  
  # Loop through clusters
  for (cl in clusters) {
    
    # Create a directory for the cluster and percentage
    cluster_percentage_dir <- paste0(cl, "/")
    if (!dir.exists(cluster_percentage_dir)) {
      dir.create(cluster_percentage_dir, recursive = TRUE)
    }
    
    # Initialize a list to store results for each iteration
    results <- list()
    
    # Loop through iterations and run DESeq2 for each
    for (j in 1:iterations) {
     
      # Extract chosen cells for the current iteration
      chosen_cells <- all_indices[[j]]
      
      result <- tryCatch({
        seurat_subset <- subset(seur_obj, cells = chosen_cells)
        deseq2_results <- single_factor_DESeq(object = seurat_subset,
                                              comp_vect = c("experience", "NC", "HC"),
                                              cluster = cl,
                                              min_cell = 1)
        deseq2_results$deseq_results  # Extract result tibble
      }, error = function(e) {
        message(paste("Error in cluster", cl, "at iteration", j, "for", target_percentage, "%: ", e$message))
        NULL  # Return NULL for failed iterations
      })
      
      # # Subset the Seurat object
      # seurat_subset <- subset(seur_obj, cells = chosen_cells)
      # 
      # # Run DESeq2
      # deseq2_results <- single_factor_DESeq(object = seurat_subset,
      #                                       comp_vect = c("experience", "NC", "HC"),
      #                                       cluster = cl,
      #                                       min_cell = 1)
     
       # pull out only result tibble from the deseq2_result list
      results[[j]] <- deseq2_results$deseq_results
    }
    
    # Save the DESeq2 results for the current cluster and percentage
    save(results, file = paste0(cl, "/PercentActiveInNC_experience_NC_HC_results_", target_percentage, "_percent.RData"))
  }
}

# Tally clusters per percent and gene
for (cl in clusters) {
  clstr <- cl
  
  # List all result files matching the pattern
  files <- list.files(clstr, pattern = "^PercentActiveInNC_experience_NC_HC_results_\\d+_percent\\.RData$", full.names = TRUE)
  #files <- list.files(clstr, pattern = "^results_.*\\.RData$", full.names = TRUE)
  
  # Load the results from each file
  results <- lapply(files, function(file) {
    env <- new.env()
    load(file, envir = env)
    as.list(env)
  })
  
  # Give the list of results meaningful names
  names(results) <- basename(files)
  
  # Remove the file extension from the names
  names(results) <- sub("\\.RData$", "", names(results))
  
  # Now tally the results and save that CSV
  all_tallies <- list()  # Collect all tallies for further processing
  for (j in names(results)) {
    merged_results_df <- tally_iterations(results, j)
    
    # Check the result
    print(paste("Processing:", j))
    print(head(merged_results_df))
    
    # Save the results to a CSV file
    file_path <- paste0(clstr, "/", j, "_iteration_tally.csv")
    print(paste("Saving to:", file_path))
    
    write.csv(merged_results_df, file = file_path, row.names = TRUE)
    
  }
}

# Calculate averages and stdevs
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


# Define clusters and numbers
clusters <- names(all_data_lists)
numbers <- c("0", "5", "10", "25", "50", "75", "100")

# Initialize data frames to store the results: one for averages, one for SD, and one for nods
average_sums_df <- data.frame(matrix(nrow = length(clusters), ncol = length(numbers)))
sd_sums_df <- data.frame(matrix(nrow = length(clusters), ncol = length(numbers)))

rownames(average_sums_df) <- clusters
colnames(average_sums_df) <- numbers
rownames(sd_sums_df) <- clusters
colnames(sd_sums_df) <- numbers

# Populate the data frames with mean, SD values, and nods data
for (cl in clusters) {
  for (num in numbers) {
    # Construct the name of the DataFrame in all_data_lists
    df_name <- paste0(cl, "_", num,'_percent')
    
    
    # Check if the DataFrame exists and compute the average and SD if it does
    if (df_name %in% names(all_data_lists[[cl]])) {
      
      df <- all_data_lists[[cl]][[df_name]]
      numeric_df <- df[, sapply(df, function(col) is.numeric(col) || is.integer(col))]
      
      # Calculate column sums, then compute mean and SD
      col_sums <- colSums(numeric_df, na.rm = TRUE)
      average_sums_df[cl, num] <- mean(col_sums)
      sd_sums_df[cl, num] <- sd(col_sums)  # Calculate SD instead of SEM
    } 
    else {
      # Assign NA if the DataFrame is missing
      average_sums_df[cl, num] <- NA
      sd_sums_df[cl, num] <- NA
    }
  }
}

# Check the resulting data frames
print(average_sums_df)
print(sd_sums_df)

# save as csv
write.csv(average_sums_df, "iteration_plot_sum_average.csv")
write.csv(sd_sums_df, "iteration_plot_stdev.csv")

