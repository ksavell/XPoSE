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
    
    # Define filenames and save (FM)
    save_filtered(fm_updnscore, fm_consistent_up,
                  file.path(cluster_dir, paste0("FindMarkers_ConsistentUp_Cluster_", cluster_name, ".csv")))
    
    save_filtered(fm_updnscore, fm_consistent_down,
                  file.path(cluster_dir, paste0("FindMarkers_ConsistentDown_Cluster_", cluster_name, ".csv")))
    
    save_filtered(fm_updnscore, fm_nonconsistent_up,
                  file.path(cluster_dir, paste0("FindMarkers_NonConsistentUp_Cluster_", cluster_name, ".csv")))
    
    save_filtered(fm_updnscore, fm_nonconsistent_down,
                  file.path(cluster_dir, paste0("FindMarkers_NonConsistentDown_Cluster_", cluster_name, ".csv")))
    
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
    cat("❌ Error processing cluster:", cluster_name, "\n")
    cat("Message:", e$message, "\n")
  })
}


# Save Summary Results ----------------------------------------------------

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


