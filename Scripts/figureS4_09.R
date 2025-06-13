# Figure S4
library(Seurat)
library(tidyverse)

load("~/Projects/XPoSE/all_10312024.RData")

# FS4A --------------------------------------------------------------------

# downsample NC:A, this takes a long time to run

# Define key variables
iterations <- 100
seur_obj <- all  # Full Seurat object
seur_obj$experience <- ifelse(seur_obj$group == "Homecage", "HC", "NC")
clusters <- unique(seur_obj$cluster_name)
target_percentages <- c(0, 5, 25, 50, 75, 100)  # Adjust percentages as needed

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
}

# Check the resulting data frames
print(average_sums_df)
print(sd_sums_df)

# save as csv
write.csv(average_sums_df, "iteration_plot_sum_average.csv")
write.csv(sd_sums_df, "iteration_plot_stdev.csv")

# FS4B --------------------------------------------------------------------

# wilcoxon data is saved under no_exclusion folder in FS4C output

# FS4C --------------------------------------------------------------------

# leave one out analysis

set.seed(22)

# List of all unique rat IDs
all_rats <- unique(all$ratID)
clusters <- unique(all$cluster_name)

all_results <- list()

for (cl in clusters) {
  
  # Create a directory for the cluster
  cluster_dir <- paste0("Cluster_", cl)
  if (!dir.exists(cluster_dir)) {
    dir.create(cluster_dir)
  }
  
  for (excluded_rat in c("none", all_rats)) {
    
    # Define subdirectory for the excluded rat
    rat_dir <- if (excluded_rat == "none") {
      file.path(cluster_dir, "Excluded_none")
    } else {
      file.path(cluster_dir, paste0("Excluded_", excluded_rat))
    }
    
    if (!dir.exists(rat_dir)) {
      dir.create(rat_dir)
    }
    
    # Step 1: Exclude the rat first
    subset_data <- if (excluded_rat == "none") {
      all  # Use full dataset if no rat is excluded
    } else {
      subset(all, ratID != excluded_rat)
    }
    
    # Step 2: Set correct identity class
    Idents(subset_data) <- "group"  # Ensure "Active" and "Homecage" are the current identities
    
    tryCatch({
      # DESeq2 Analysis (Active vs. Homecage)
      results <- single_factor_DESeq(
        object = subset_data,
        comp_vect = c("group", "Active", "Homecage"),
        cluster = cl,
        min_cell = 0,
        min_rat = 2,
        keep_dds = TRUE
      )
      
      # Extract results and dds
      deseq_results <- results$results
      dds <- results$dds
      
      # Ensure valid row names for deseq_results
      if (is.null(rownames(deseq_results))) {
        rownames(deseq_results) <- paste0("Gene_", seq_len(nrow(deseq_results)))
      } else if (any(is.na(rownames(deseq_results)))) {
        rownames(deseq_results)[is.na(rownames(deseq_results))] <- paste0("Gene_", which(is.na(rownames(deseq_results))))
      }
      
      # Add score column
      deseq_results$score <- ifelse(deseq_results$padj < 0.05, 1, 0)
      deseq_results$score_updn <- ifelse(deseq_results$padj < 0.05 & deseq_results$log2FoldChange > 0, 1,  # Upregulated
                                         ifelse(deseq_results$padj < 0.05 & deseq_results$log2FoldChange < 0, 2,  # Downregulated
                                                0))
      
      # Save DESeq2 results
      result_file <- file.path(rat_dir, paste0("DESeq_Results_Group_Active_Homecage_Excluded_", excluded_rat, ".csv"))
      write.csv(deseq_results, result_file, row.names = TRUE)
      
      result_dds <- file.path(rat_dir, paste0("DESeq2_Group_Active_Homecage_Excluded_", excluded_rat, ".rds"))
      saveRDS(dds, file = result_dds, compress = FALSE)
      
       # FindMarkers Analysis (for Seurat object)
       cluster_subset <- subset(subset_data, cluster_name == cl)  # Subset by cluster
       Idents(cluster_subset) <- "group"  # Set identity to group again
       fm_results <- FindMarkers(cluster_subset, ident.1 = "Active", ident.2 = "Homecage")  # Compare Active vs. Homecage
       fm_results$score <- ifelse(fm_results$p_val_adj < 0.05, 1, 0)
       fm_results$score_updn <- ifelse(fm_results$p_val_adj < 0.05 & fm_results$avg_log2FC > 0, 1,  # Upregulated
                                       ifelse(fm_results$p_val_adj < 0.05 & fm_results$avg_log2FC < 0, 2,  # Downregulated
                                              0))

      # # Save FindMarkers results
       fm_file <- file.path(rat_dir, paste0("FindMarkers_Results_Group_Active_Homecage_Excluded_", excluded_rat, ".csv"))
       write.csv(fm_results, fm_file, row.names = TRUE)

    }, error = function(e) {
      cat("Error in analysis for cluster:", cl, "and excluded rat:", excluded_rat, "\n")
      cat("Error message:", e$message, "\n")
    })
  }
}

for (cluster_dir in cluster_dirs) {
  cluster_name <- basename(cluster_dir)
  
  tryCatch({
    # Define filenames
    fm_updnscore_file <- file.path(cluster_dir, paste0("FindMarkers_combined_results_Cluster_", cluster_name, ".csv"))
    deseq_updnscore_file <- file.path(cluster_dir, paste0("DESeq_combined_results_Cluster_", cluster_name, ".csv"))
    
    # Read in score tables
    fm_updnscore <- read.csv(fm_updnscore_file, row.names = 1, check.names = FALSE)
    deseq_updnscore <- read.csv(deseq_updnscore_file, row.names = 1, check.names = FALSE)
    
    # Consistent (all 1 or all 2)
    fm_consistent_up <- rownames(fm_updnscore)[apply(fm_updnscore, 1, \(x) all(x == 1))]
    fm_consistent_down <- rownames(fm_updnscore)[apply(fm_updnscore, 1, \(x) all(x == 2))]
    
    deseq_consistent_up <- rownames(deseq_updnscore)[apply(deseq_updnscore, 1, \(x) all(x == 1))]
    deseq_consistent_down <- rownames(deseq_updnscore)[apply(deseq_updnscore, 1, \(x) all(x == 2))]
    
    # Non-consistent (includes mixed or missing agreement)
    nonconsist_logic <- function(x, val) {
      any(x == val) & !all(x == val) & (!"Excluded_none" %in% names(x) || x["Excluded_none"] != 0)
    }
    
    fm_nonconsistent_up <- rownames(fm_updnscore)[apply(fm_updnscore, 1, nonconsist_logic, val = 1)]
    fm_nonconsistent_down <- rownames(fm_updnscore)[apply(fm_updnscore, 1, nonconsist_logic, val = 2)]
    
    deseq_nonconsistent_up <- rownames(deseq_updnscore)[apply(deseq_updnscore, 1, nonconsist_logic, val = 1)]
    deseq_nonconsistent_down <- rownames(deseq_updnscore)[apply(deseq_updnscore, 1, nonconsist_logic, val = 2)]
    
    # Save filtered tables
    save_filtered <- function(df, genes, filename) {
      out <- df[rownames(df) %in% genes, , drop = FALSE]
      if (nrow(out) > 0) {
        write.csv(out, filename, row.names = TRUE)
      } else {
        cat("Saving empty file:", basename(filename), "\n")
        write.csv(data.frame(Gene = character(), Score = integer()), filename, row.names = FALSE)
      }
    }
    
    # # Define filenames and save (FM)
    # save_filtered(fm_updnscore, fm_consistent_up,
    #               file.path(cluster_dir, paste0("FindMarkers_ConsistentUp_Cluster_", cluster_name, ".csv")))
    # 
    # save_filtered(fm_updnscore, fm_consistent_down,
    #               file.path(cluster_dir, paste0("FindMarkers_ConsistentDown_Cluster_", cluster_name, ".csv")))
    # 
    # save_filtered(fm_updnscore, fm_nonconsistent_up,
    #               file.path(cluster_dir, paste0("FindMarkers_NonConsistentUp_Cluster_", cluster_name, ".csv")))
    # 
    # save_filtered(fm_updnscore, fm_nonconsistent_down,
    #               file.path(cluster_dir, paste0("FindMarkers_NonConsistentDown_Cluster_", cluster_name, ".csv")))
    # 
    # Define filenames and save (DESeq)
    save_filtered(deseq_updnscore, deseq_consistent_up,
                  file.path(cluster_dir, paste0("DESeq_ConsistentUp_Cluster_", cluster_name, ".csv")))
    
    save_filtered(deseq_updnscore, deseq_consistent_down,
                  file.path(cluster_dir, paste0("DESeq_ConsistentDown_Cluster_", cluster_name, ".csv")))
    
    save_filtered(deseq_updnscore, deseq_nonconsistent_up,
                  file.path(cluster_dir, paste0("DESeq_NonConsistentUp_Cluster_", cluster_name, ".csv")))
    
    save_filtered(deseq_updnscore, deseq_nonconsistent_down,
                  file.path(cluster_dir, paste0("DESeq_NonConsistentDown_Cluster_", cluster_name, ".csv")))
    
    cat("✅ Processed up/down consistent and non-consistent genes for cluster:", cluster_name, "\n")
    
  }, error = function(e) {
    cat("Error processing cluster:", cluster_name, "\n")
    cat("Message:", e$message, "\n")
  })
}

count_if_exists_rows <- function(filepath, filter_val, consistent = TRUE) {
  if (!file.exists(filepath)) return(0)
  
  df <- read.csv(filepath, row.names = 1, check.names = FALSE)
  
  if (consistent) {
    # All values must be the same
    sum(apply(df, 1, function(x) all(x == filter_val)))
  } else {
    # Any value matches filter_val
    sum(apply(df, 1, function(x) any(x == filter_val)))
  }
}


# Initialize detailed summary table
consistency_summary <- data.frame(
  Cluster = character(),
  consistent_DESeq_up = integer(),
  consistent_DESeq_down = integer(),
  nonconsistent_DESeq_up = integer(),
  nonconsistent_DESeq_down = integer(),
  consistent_FM_up = integer(),
  consistent_FM_down = integer(),
  nonconsistent_FM_up = integer(),
  nonconsistent_FM_down = integer(),
  stringsAsFactors = FALSE
)

# Loop to populate the new summary
for (cluster_dir in cluster_dirs) {
  cluster_name <- basename(cluster_dir)
  
  # File paths
  deseq_consistent_file <- file.path(cluster_dir, paste0("DESeq2_Results_updnConsistent_Group_Active_Homecage_Cluster_", cluster_name, ".csv"))
  deseq_nonconsistent_file <- file.path(cluster_dir, paste0("DESeq2_Results_updnNonConsistent_Group_Active_Homecage_Cluster_", cluster_name, ".csv"))
  fm_consistent_file <- file.path(cluster_dir, paste0("FindMarkers_Results_updnConsistent_Group_Active_Homecage_Cluster_", cluster_name, ".csv"))
  fm_nonconsistent_file <- file.path(cluster_dir, paste0("FindMarkers_Results_updnNonConsistent_Group_Active_Homecage_Cluster_", cluster_name, ".csv"))
  
  summary_row <- data.frame(
    Cluster = cluster_name,
    consistent_DESeq_up = count_if_exists_rows(deseq_consistent_file, 1, consistent = TRUE),
    consistent_DESeq_down = count_if_exists_rows(deseq_consistent_file, 2, consistent = TRUE),
    nonconsistent_DESeq_up = count_if_exists_rows(deseq_nonconsistent_file, 1, consistent = FALSE),
    nonconsistent_DESeq_down = count_if_exists_rows(deseq_nonconsistent_file, 2, consistent = FALSE),
    consistent_FM_up = count_if_exists_rows(fm_consistent_file, 1, consistent = TRUE),
    consistent_FM_down = count_if_exists_rows(fm_consistent_file, 2, consistent = TRUE),
    nonconsistent_FM_up = count_if_exists_rows(fm_nonconsistent_file, 1, consistent = FALSE),
    nonconsistent_FM_down = count_if_exists_rows(fm_nonconsistent_file, 2, consistent = FALSE),
    stringsAsFactors = FALSE
  )
  
  
  consistency_summary <- rbind(consistency_summary, summary_row)
}

# Save updated summary
summary_output_file <- file.path(base_dir, "Cluster_Consistency_Summary_Directional.csv")
write.csv(consistency_summary, summary_output_file, row.names = FALSE)
cat("Updated cluster-level summary with up/down direction saved to:", summary_output_file, "\n")