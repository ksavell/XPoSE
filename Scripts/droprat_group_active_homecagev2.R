# Droprat for group, Active vs Homecage

library(dplyr)
library(tidyr)
library(Seurat)

source("Scripts/Functions/single_factor_DESeq.R") 
load("all_10312024.RData")


# Run DESeq and FM for dropped iterations and no drop ---------------------

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
        keep_dds = FALSE
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
      
      # Save FindMarkers results
      fm_file <- file.path(rat_dir, paste0("FindMarkers_Results_Group_Active_Homecage_Excluded_", excluded_rat, ".csv"))
      write.csv(fm_results, fm_file, row.names = TRUE)
      
    }, error = function(e) {
      cat("Error in analysis for cluster:", cl, "and excluded rat:", excluded_rat, "\n")
      cat("Error message:", e$message, "\n")
    })
  }
}


# Calculate gene overlap --------------------------------------------------

# Function to process each cluster and merge DESeq and FindMarkers results
# Function to process each cluster and merge DESeq and FindMarkers results
process_cluster <- function(cluster_dir) {
  # Extract the cluster name from the directory
  cluster_name <- basename(cluster_dir)
  
  # Get a list of Excluded directories within the current cluster directory
  excluded_dirs <- list.dirs(cluster_dir, full.names = TRUE, recursive = FALSE)
  
  # Initialize empty data frames to store results for DESeq and FindMarkers
  deseq_results <- data.frame()
  findmarkers_results <- data.frame()
  
  # Loop through each excluded directory
  for (excluded_dir in excluded_dirs) {
    # Get the list of DESeq and FindMarkers files within the Excluded directory
    deseq_file <- list.files(excluded_dir, pattern = "DESeq_Results_Group_Active_Homecage_Excluded_.*\\.csv$", full.names = TRUE)
    findmarkers_file <- list.files(excluded_dir, pattern = "FindMarkers_Results_Group_Active_Homecage_Excluded_.*\\.csv$", full.names = TRUE)
    
    # Extract the sample name from the directory
    sample_name <- basename(excluded_dir)
    
    # Process DESeq file
    if (length(deseq_file) == 1) {
      deseq_data <- read.csv(deseq_file, stringsAsFactors = FALSE)
      
      # Extract the gene and score_updn columns
      deseq_data <- deseq_data %>% select(gene, score_updn)
      
      # Rename the column dynamically (fixing the `:=` error)
      colnames(deseq_data)[colnames(deseq_data) == "score_updn"] <- sample_name
      
      # Merge with the main results data frame
      if (nrow(deseq_results) == 0) {
        deseq_results <- deseq_data
      } else {
        deseq_results <- merge(deseq_results, deseq_data, by = "gene", all = TRUE)
      }
    }
    
    # Process FindMarkers file
    if (length(findmarkers_file) == 1) {
      findmarkers_data <- read.csv(findmarkers_file, row.names = 1, stringsAsFactors = FALSE)
      
      # Ensure row names are preserved
      findmarkers_data <- findmarkers_data %>% mutate(gene = rownames(findmarkers_data))
      
      # Extract only the score_updn column
      findmarkers_data <- findmarkers_data %>% select(gene, score_updn)
      
      # Rename the column dynamically (fixing the `:=` error)
      colnames(findmarkers_data)[colnames(findmarkers_data) == "score_updn"] <- sample_name
      
      # Merge with the main results data frame
      if (nrow(findmarkers_results) == 0) {
        findmarkers_results <- findmarkers_data
      } else {
        findmarkers_results <- merge(findmarkers_results, findmarkers_data, by = "gene", all = TRUE)
      }
    }
  }
  
  # Save the results to the cluster directory with cluster name in filenames
  if (nrow(deseq_results) > 0) {
    write.csv(deseq_results, file.path(cluster_dir, paste0("DESeq_combined_results_Cluster_", cluster_name, ".csv")), row.names = FALSE)
  }
  if (nrow(findmarkers_results) > 0) {
    write.csv(findmarkers_results, file.path(cluster_dir, paste0("FindMarkers_combined_results_Cluster_", cluster_name, ".csv")), row.names = FALSE)
  }
}


# Define the main directory containing "Cluster_[cl]" directories
base_dir <- "~/Projects/XPoSE/droprat"

# Get the list of "Cluster_[cl]" directories
cluster_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

# Loop through each "Cluster_[cl]" directory and process
for (cluster_dir in cluster_dirs) {
  process_cluster(cluster_dir)
}

# Find Consistent and Non-consistent Genes ---------------------------------
for (cluster_dir in cluster_dirs) {
  cluster_name <- basename(cluster_dir)
  
  tryCatch({
    # Define filenames with cluster name included
    fm_updnscore_file <- file.path(cluster_dir, paste0("FindMarkers_combined_results_Cluster_", cluster_name, ".csv"))
    deseq_updnscore_file <- file.path(cluster_dir, paste0("DESeq_combined_results_Cluster_", cluster_name, ".csv"))
    
    # Read in score files
    fm_updnscore <- read.csv(fm_updnscore_file, row.names = 1, check.names = FALSE)
    deseq_updnscore <- read.csv(deseq_updnscore_file, row.names = 1, check.names = FALSE)
    
    # Identify consistent genes (all 1 or 2)
    fm_consistent <- rownames(fm_updnscore)[apply(fm_updnscore, 1, function(row) all(row == 1 | row == 2))]
    deseq_consistent <- rownames(deseq_updnscore)[apply(deseq_updnscore, 1, function(row) all(row == 1 | row == 2))]
    
    # Identify non-consistent genes
    if ("Excluded_none" %in% colnames(fm_updnscore)) {
      fm_nonconsistent <- rownames(fm_updnscore)[apply(fm_updnscore, 1, function(row) 
        any(row == 1 | row == 2) & !all(row == 1 | row == 2) & row["Excluded_none"] != 0)]
    } else {
      fm_nonconsistent <- rownames(fm_updnscore)[apply(fm_updnscore, 1, function(row) 
        any(row == 1 | row == 2) & !all(row == 1 | row == 2))]
    }
    
    if ("Excluded_none" %in% colnames(deseq_updnscore)) {
      deseq_nonconsistent <- rownames(deseq_updnscore)[apply(deseq_updnscore, 1, function(row) 
        any(row == 1 | row == 2) & !all(row == 1 | row == 2) & row["Excluded_none"] != 0)]
    } else {
      deseq_nonconsistent <- rownames(deseq_updnscore)[apply(deseq_updnscore, 1, function(row) 
        any(row == 1 | row == 2) & !all(row == 1 | row == 2))]
    }
    
    # Debugging: Check if non-consistent gene lists are empty
    cat("Cluster:", cluster_name, " - FM Non-Consistent Genes:", length(fm_nonconsistent), "\n")
    cat("Cluster:", cluster_name, " - DESeq Non-Consistent Genes:", length(deseq_nonconsistent), "\n")
    
    # Filter updnscore files by consistent and non-consistent genes
    fm_updn_consistent <- fm_updnscore[rownames(fm_updnscore) %in% fm_consistent, , drop = FALSE]
    fm_updn_nonconsistent <- fm_updnscore[rownames(fm_updnscore) %in% fm_nonconsistent, , drop = FALSE]
    
    deseq_updn_consistent <- deseq_updnscore[rownames(deseq_updnscore) %in% deseq_consistent, , drop = FALSE]
    deseq_updn_nonconsistent <- deseq_updnscore[rownames(deseq_updnscore) %in% deseq_nonconsistent, , drop = FALSE]
    
    # Ensure non-consistent files are saved even if empty
    fm_nonconsistent_file <- file.path(cluster_dir, paste0("FindMarkers_Results_updnNonConsistent_Group_Active_Homecage_Cluster_", cluster_name, ".csv"))
    deseq_nonconsistent_file <- file.path(cluster_dir, paste0("DESeq2_Results_updnNonConsistent_Group_Active_Homecage_Cluster_", cluster_name, ".csv"))
    
    write.csv(fm_updn_consistent, file.path(cluster_dir, paste0("FindMarkers_Results_updnConsistent_Group_Active_Homecage_Cluster_", cluster_name, ".csv")), row.names = TRUE)
    write.csv(deseq_updn_consistent, file.path(cluster_dir, paste0("DESeq2_Results_updnConsistent_Group_Active_Homecage_Cluster_", cluster_name, ".csv")), row.names = TRUE)
    
    # **Ensuring non-consistent files save even if empty**
    if (nrow(fm_updn_nonconsistent) > 0) {
      write.csv(fm_updn_nonconsistent, fm_nonconsistent_file, row.names = TRUE)
    } else {
      cat("Cluster:", cluster_name, " - No Non-Consistent FindMarkers Genes. Saving empty file.\n")
      write.csv(data.frame(Gene = character(), Score = integer()), fm_nonconsistent_file, row.names = FALSE)
    }
    
    if (nrow(deseq_updn_nonconsistent) > 0) {
      write.csv(deseq_updn_nonconsistent, deseq_nonconsistent_file, row.names = TRUE)
    } else {
      cat("Cluster:", cluster_name, " - No Non-Consistent DESeq Genes. Saving empty file.\n")
      write.csv(data.frame(Gene = character(), Score = integer()), deseq_nonconsistent_file, row.names = FALSE)
    }
    
    cat("Processed consistent and non-consistent updn scores for cluster:", cluster_name, "\n")
    
  }, error = function(e) {
    cat("Error processing cluster:", cluster_name, "\n")
    cat("Error message:", e$message, "\n")
  })
}


# Save Summary Results ----------------------------------------------------

# Initialize summary table
consistency_summary <- data.frame(
  Cluster = character(),
  consistent_DESeq = integer(),
  nonconsistent_DESeq = integer(),
  consistent_FM = integer(),
  nonconsistent_FM = integer(),
  stringsAsFactors = FALSE
)

# Loop again to count consistent and non-consistent genes
for (cluster_dir in cluster_dirs) {
  cluster_name <- basename(cluster_dir)
  
  # File paths (must match saved filenames)
  deseq_consistent_file <- file.path(cluster_dir, paste0("DESeq2_Results_updnConsistent_Group_Active_Homecage_Cluster_", cluster_name, ".csv"))
  deseq_nonconsistent_file <- file.path(cluster_dir, paste0("DESeq2_Results_updnNonConsistent_Group_Active_Homecage_Cluster_", cluster_name, ".csv"))
  fm_consistent_file <- file.path(cluster_dir, paste0("FindMarkers_Results_updnConsistent_Group_Active_Homecage_Cluster_", cluster_name, ".csv"))
  fm_nonconsistent_file <- file.path(cluster_dir, paste0("FindMarkers_Results_updnNonConsistent_Group_Active_Homecage_Cluster_", cluster_name, ".csv"))
  
  # Read counts if files exist, else assume 0
  count_if_exists <- function(filepath) {
    if (file.exists(filepath)) {
      nrow(read.csv(filepath))
    } else {
      0
    }
  }
  
  summary_row <- data.frame(
    Cluster = cluster_name,
    consistent_DESeq = count_if_exists(deseq_consistent_file),
    nonconsistent_DESeq = count_if_exists(deseq_nonconsistent_file),
    consistent_FM = count_if_exists(fm_consistent_file),
    nonconsistent_FM = count_if_exists(fm_nonconsistent_file),
    stringsAsFactors = FALSE
  )
  
  consistency_summary <- rbind(consistency_summary, summary_row)
}

# Save the summary table
summary_output_file <- file.path(base_dir, "Cluster_Consistency_Summary.csv")
write.csv(consistency_summary, summary_output_file, row.names = FALSE)

cat("Cluster-level consistency summary saved to:", summary_output_file, "\n")


