# drop one rat analysis

library(dplyr)
library(tidyr)
library(UpSetR)


source("Scripts/Functions/single_factor_DESeq.R") 
load("all_10312024.RData")

all$experience <- ifelse(all$group == "Homecage", "HC",
                         "NC")

# List of all unique rat IDs
all_rats <- unique(all$ratID)
clusters <- unique(all$cluster_name)

# Placeholder for results
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
    
    # Subset data based on exclusion
    subset_data <- if (excluded_rat == "none") {
      all  # Use full dataset for "none"
    } else {
      subset(all, ratID != excluded_rat)
    }
    
    tryCatch({
      # DESeq2 Analysis
      results <- single_factor_DESeq(
        object = subset_data,
        comp_vect = c("experience", "NC", "HC"),
        cluster = cl,
        min_cell = 1,
        min_rat = 2,
        keep_dds = FALSE
      )
      
      # Extract results and dds
      deseq_results <- results$deseq_results
      dds <- results$clust_tbl
      
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
      result_file <- file.path(rat_dir, paste0("DESeq_Results_Experience_NC_HC_Excluded_", excluded_rat, ".csv"))
      write.csv(deseq_results, result_file, row.names = TRUE)
      
      result_dds <- file.path(rat_dir, paste0("DESeq2_Experience_NC_HC_Excluded_", excluded_rat, ".rds"))
      saveRDS(dds, file = result_dds, compress = FALSE)
      
      # FindMarkers Analysis (for Seurat object)
      cluster_subset <- subset(subset_data, cluster_name == cl)  # Subset by cluster
      Idents(cluster_subset) <- "experience"
      fm_results <- FindMarkers(cluster_subset, ident.1 = "NC", ident.2 = "HC")
      fm_results$score <- ifelse(fm_results$p_val_adj < 0.05, 1, 0)
      fm_results$score_updn <- ifelse(fm_results$p_val_adj < 0.05 & fm_results$avg_log2FC > 0, 1,  # Upregulated
                                    ifelse(fm_results$p_val_adj < 0.05 & fm_results$avg_log2FC < 0, 2,  # Downregulated
                                           0))
      
      # Save FindMarkers results
      fm_file <- file.path(rat_dir, paste0("FindMarkers_Results_Experience_NC_HC_Excluded_", excluded_rat, ".csv"))
      write.csv(fm_results, fm_file, row.names = TRUE)
      
    }, error = function(e) {
      cat("Error in analysis for cluster:", cl, "and excluded rat:", excluded_rat, "\n")
      cat("Error message:", e$message, "\n")
    })
  }
}

# gene overlap ----------------------------------------------------

excluded_dirs <- c("Excluded_none", "Excluded_HC-1", "Excluded_HC-2", "Excluded_HC-3", "Excluded_HC-4", 
                   "Excluded_NC-1", "Excluded_NC-2", "Excluded_NC-3", "Excluded_NC-4")

# Define the clusters and directories
clusters <- list.dirs(path = getwd(), recursive = FALSE, full.names = FALSE)
clusters <- clusters[grepl("^Cluster_", clusters)] # Keep only cluster directories

for (cluster in clusters) {
  tryCatch({
    # Set the cluster directory
    cluster_dir <- file.path(getwd(), cluster)
    excluded_dirs <- list.dirs(cluster_dir, recursive = FALSE, full.names = TRUE)
    
    # Initialize lists to store data for each file type
    fm_binary_list <- list()
    fm_updn_list <- list()
    deseq_binary_list <- list()
    deseq_updn_list <- list()
    
    # Process each excluded directory
    for (excluded_dir in excluded_dirs) {
      exclusion_name <- basename(excluded_dir)  # Get the name of exclusion (e.g., Excluded_HC-1)
      
      # FindMarkers processing
      fm_file <- file.path(excluded_dir, paste0("FindMarkers_Results_Experience_NC_HC_", exclusion_name, ".csv"))
      if (file.exists(fm_file)) {
        fm_data <- read.csv(fm_file, row.names = 1)  # Read row names as gene
        
        # Ensure required columns exist
        if (all(c("score", "score_updn") %in% colnames(fm_data))) {
          fm_data$gene <- rownames(fm_data)
          fm_binary_list[[exclusion_name]] <- fm_data[, c("gene", "score")]
          fm_binary_list[[exclusion_name]] <- setNames(fm_binary_list[[exclusion_name]], c("gene", exclusion_name))
          
          fm_updn_list[[exclusion_name]] <- fm_data[, c("gene", "score_updn")]
          fm_updn_list[[exclusion_name]] <- setNames(fm_updn_list[[exclusion_name]], c("gene", exclusion_name))
        } else {
          cat("FindMarkers file missing required columns:", fm_file, "\n")
          next
        }
      } else {
        cat("FindMarkers file not found in directory:", excluded_dir, "\n")
      }
      
      # DESeq2 processing
      deseq_file <- file.path(excluded_dir, paste0("DESeq_Results_Experience_NC_HC_", exclusion_name, ".csv"))
      if (file.exists(deseq_file)) {
        deseq_data <- read.csv(deseq_file)  # Genes are already a column
        
        # Ensure required columns exist
        if (all(c("gene", "score", "score_updn") %in% colnames(deseq_data))) {
          deseq_binary_list[[exclusion_name]] <- deseq_data[, c("gene", "score")]
          deseq_binary_list[[exclusion_name]] <- setNames(deseq_binary_list[[exclusion_name]], c("gene", exclusion_name))
          
          deseq_updn_list[[exclusion_name]] <- deseq_data[, c("gene", "score_updn")]
          deseq_updn_list[[exclusion_name]] <- setNames(deseq_updn_list[[exclusion_name]], c("gene", exclusion_name))
        } else {
          cat("DESeq2 file missing required columns:", deseq_file, "\n")
          next
        }
      } else {
        cat("DESeq2 file not found in directory:", excluded_dir, "\n")
      }
    }
    
    # Function to merge dataframes by gene
    combine_tables <- function(data_list) {
      combined <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), data_list)
      combined[is.na(combined)] <- 0  # Replace NAs with 0
      return(combined)
    }
    
    # Combine tables for each method/score type
    if (length(fm_binary_list) > 0) {
      fm_binary_data <- combine_tables(fm_binary_list)
      fm_updn_data <- combine_tables(fm_updn_list)
      
      # Save combined tables
      write.csv(fm_binary_data, file.path(cluster_dir, paste0("FindMarkers_Results_binaryTally_Experience_NC_HC_", cluster, ".csv")), row.names = FALSE)
      write.csv(fm_updn_data, file.path(cluster_dir, paste0("FindMarkers_Results_updnTally_Experience_NC_HC_", cluster, ".csv")), row.names = FALSE)
    }
    
    if (length(deseq_binary_list) > 0) {
      deseq_binary_data <- combine_tables(deseq_binary_list)
      deseq_updn_data <- combine_tables(deseq_updn_list)
      
      # Save combined tables
      write.csv(deseq_binary_data, file.path(cluster_dir, paste0("DESeq2_Results_binaryTally_Experience_NC_HC_", cluster, ".csv")), row.names = FALSE)
      write.csv(deseq_updn_data, file.path(cluster_dir, paste0("DESeq2_Results_updnTally_Experience_NC_HC_", cluster, ".csv")), row.names = FALSE)
    }
    
    cat("Processing completed for cluster:", cluster, "\n")
    
  }, error = function(e) {
    cat("Error processing cluster:", cluster, "\n")
    cat("Error message:", e$message, "\n")
  })
}


# Find consistent and nonconsistent DEGs ----------------------------------

for (cluster in clusters) {
  tryCatch({
    # Define the cluster directory
    cluster_dir <- file.path(getwd(), cluster)
    if (!dir.exists(cluster_dir)) {
      dir.create(cluster_dir)
    }
    
    # Load score files
    fm_score_file <- file.path(cluster_dir, paste0("FindMarkers_Results_binaryTally_Experience_NC_HC_", cluster, ".csv"))
    fm_updnscore_file <- file.path(cluster_dir, paste0("FindMarkers_Results_updnTally_Experience_NC_HC_", cluster, ".csv"))
    deseq_score_file <- file.path(cluster_dir, paste0("DESeq2_Results_binaryTally_Experience_NC_HC_", cluster, ".csv"))
    deseq_updnscore_file <- file.path(cluster_dir, paste0("DESeq2_Results_updnTally_Experience_NC_HC_", cluster, ".csv"))
    
    # Check if required files exist
    if (!all(file.exists(c(fm_score_file, fm_updnscore_file, deseq_score_file, deseq_updnscore_file)))) {
      cat("Score files missing for cluster:", cluster, "\n")
      next
    }
    
    # Read in score and updnscore files
    fm_score <- read.csv(fm_score_file, row.names = 1, check.names = FALSE)
    fm_updnscore <- read.csv(fm_updnscore_file, row.names = 1, check.names = FALSE)
    deseq_score <- read.csv(deseq_score_file, row.names = 1, check.names = FALSE)
    deseq_updnscore <- read.csv(deseq_updnscore_file, row.names = 1, check.names = FALSE)
    
    # Identify consistent (all 1's) and nonconsistent (at least 1 1, not all 0's) genes
    fm_consistent <- rownames(fm_score)[apply(fm_score, 1, function(row) all(row == 1))]
    fm_nonconsistent <- rownames(fm_score)[apply(fm_score, 1, function(row) any(row == 1) & !all(row == 1))]
    
    deseq_consistent <- rownames(deseq_score)[apply(deseq_score, 1, function(row) all(row == 1))]
    deseq_nonconsistent <- rownames(deseq_score)[apply(deseq_score, 1, function(row) any(row == 1) & !all(row == 1))]
    
    # Filter updnscore files by consistent and nonconsistent genes
    fm_updn_consistent <- fm_updnscore[rownames(fm_updnscore) %in% fm_consistent, , drop = FALSE]
    fm_updn_nonconsistent <- fm_updnscore[rownames(fm_updnscore) %in% fm_nonconsistent, , drop = FALSE]
    
    deseq_updn_consistent <- deseq_updnscore[rownames(deseq_updnscore) %in% deseq_consistent, , drop = FALSE]
    deseq_updn_nonconsistent <- deseq_updnscore[rownames(deseq_updnscore) %in% deseq_nonconsistent, , drop = FALSE]
    
    # Save filtered updnscore files
    fm_consistent_file <- file.path(cluster_dir, paste0("FindMarkers_Results_updnConsistent_Experience_NC_HC_", cluster, ".csv"))
    fm_nonconsistent_file <- file.path(cluster_dir, paste0("FindMarkers_Results_updnNonConsistent_Experience_NC_HC_", cluster, ".csv"))
    write.csv(fm_updn_consistent, fm_consistent_file, row.names = TRUE)
    write.csv(fm_updn_nonconsistent, fm_nonconsistent_file, row.names = TRUE)
    
    deseq_consistent_file <- file.path(cluster_dir, paste0("DESeq2_Results_updnConsistent_Experience_NC_HC_", cluster, ".csv"))
    deseq_nonconsistent_file <- file.path(cluster_dir, paste0("DESeq2_Results_updnNonConsistent_Experience_NC_HC_", cluster, ".csv"))
    write.csv(deseq_updn_consistent, deseq_consistent_file, row.names = TRUE)
    write.csv(deseq_updn_nonconsistent, deseq_nonconsistent_file, row.names = TRUE)
    
    cat("Processed consistent and nonconsistent updn scores for cluster:", cluster, "\n")
    
  }, error = function(e) {
    cat("Error processing cluster:", cluster, "\n")
    cat("Error message:", e$message, "\n")
  })
}


# now calculate comparison table ------------------------------------------

# Define clusters
clusters <- list.dirs(path = getwd(), recursive = FALSE, full.names = FALSE)
clusters <- clusters[grepl("^Cluster_", clusters)]  # Keep only cluster directories

# Define overlap categories
overlap_categories <- c(
  "Both_Downregulated",
  "Neutral_FM_Down_DES2",
  "Up_FM_Down_DES2",
  "Down_FM_Neutral_DES2",
  "Both_Neutral",
  "Up_FM_Neutral_DES2",
  "Down_FM_Up_DES2",
  "Neutral_FM_Up_DES2",
  "Both_Upregulated"
)

# Initialize results list
consistent_results <- list()
nonconsistent_results <- list()

# Define the overlap function
define_overlap <- function(fm_data, deseq_data) {
  all_genes <- union(rownames(fm_data), rownames(deseq_data))
  fm_subset <- fm_data[all_genes, , drop = FALSE]
  deseq_subset <- deseq_data[all_genes, , drop = FALSE]
  
  # Replace NAs with 0 (Neutral)
  fm_subset[is.na(fm_subset)] <- 0
  deseq_subset[is.na(deseq_subset)] <- 0
  
  # Calculate overlap categories
  overlap <- ifelse(
    fm_subset == 2 & deseq_subset == 2, "Both_Downregulated",
    ifelse(
      fm_subset == 0 & deseq_subset == 2, "Neutral_FM_Down_DES2",
      ifelse(
        fm_subset == 1 & deseq_subset == 2, "Up_FM_Down_DES2",
        ifelse(
          fm_subset == 2 & deseq_subset == 0, "Down_FM_Neutral_DES2",
          ifelse(
            fm_subset == 0 & deseq_subset == 0, "Both_Neutral",
            ifelse(
              fm_subset == 1 & deseq_subset == 0, "Up_FM_Neutral_DES2",
              ifelse(
                fm_subset == 2 & deseq_subset == 1, "Down_FM_Up_DES2",
                ifelse(
                  fm_subset == 0 & deseq_subset == 1, "Neutral_FM_Up_DES2",
                  ifelse(
                    fm_subset == 1 & deseq_subset == 1, "Both_Upregulated",
                    NA
                  )
                )
              )
            )
          )
        )
      )
    )
  )
  
  return(overlap)
}

for (cluster in clusters) {
  tryCatch({
    # Define cluster directory
    #cluster_dir <- file.path(getwd(), cluster)
    cat("Processing cluster:", cluster, "\n")
    
    # Load consistent and non-consistent data
    fm_consistent_file <- paste0(cluster, "/FindMarkers_Results_updnConsistent_Experience_NC_HC_", cluster, ".csv")
    fm_nonconsistent_file <- paste0(cluster, "/FindMarkers_Results_updnNonConsistent_Experience_NC_HC_", cluster, ".csv")
    deseq_consistent_file <- paste0(cluster, "/DESeq2_Results_updnConsistent_Experience_NC_HC_", cluster, ".csv")
    deseq_nonconsistent_file <- paste0(cluster, "/DESeq2_Results_updnNonConsistent_Experience_NC_HC_", cluster, ".csv")
    
    fm_consistent <- read.csv(fm_consistent_file, row.names = 1, check.names = FALSE)
    fm_nonconsistent <- read.csv(fm_nonconsistent_file, row.names = 1, check.names = FALSE)
    deseq_consistent <- read.csv(deseq_consistent_file, row.names = 1, check.names = FALSE)
    deseq_nonconsistent <- read.csv(deseq_nonconsistent_file, row.names = 1, check.names = FALSE)
    
    # Calculate overlaps
    consistent_overlap <- define_overlap(fm_consistent, deseq_consistent)
    nonconsistent_overlap <- define_overlap(fm_nonconsistent, deseq_nonconsistent)
    
    # Count overlap categories
    consistent_counts <- table(factor(consistent_overlap, levels = overlap_categories))
    nonconsistent_counts <- table(factor(nonconsistent_overlap, levels = overlap_categories))
    
    # Create dataframes for counts, filling missing categories with 0
    consistent_df <- data.frame(Cluster = cluster, t(as.data.frame(consistent_counts)))
    nonconsistent_df <- data.frame(Cluster = cluster, t(as.data.frame(nonconsistent_counts)))
    
    # Append to results lists
    consistent_results[[cluster]] <- consistent_df
    nonconsistent_results[[cluster]] <- nonconsistent_df
    
  }, error = function(e) {
    cat("Error processing cluster:", cluster, "\n")
    cat("Error message:", e$message, "\n")
  })
}

clean_up_results <- function(results_list, overlap_categories) {
  # Initialize cleaned results with zero counts
  cleaned_results <- data.frame(matrix(0, nrow = length(results_list), ncol = length(overlap_categories)))
  colnames(cleaned_results) <- overlap_categories
  rownames(cleaned_results) <- names(results_list)
  
  # Process each cluster
  for (cluster in names(results_list)) {
    if (!is.null(results_list[[cluster]]) && nrow(results_list[[cluster]]) > 0) {
      cluster_data <- results_list[[cluster]]
      
      # Debugging: Print cluster data
      cat("\nProcessing cluster:", cluster, "\n")
      print(cluster_data)
      
      # Extract categories and counts
      categories <- as.character(cluster_data["Var1", ])
      counts <- as.numeric(cluster_data["Freq", ])
      
      # Map counts to the appropriate categories
      cluster_counts <- setNames(rep(0, length(overlap_categories)), overlap_categories)
      cluster_counts[categories] <- counts
      
      # Assign counts to the cleaned results
      cleaned_results[cluster, ] <- cluster_counts
    } else {
      cat("No valid data for cluster:", cluster, "\n")
    }
  }
  
  return(cleaned_results)
}

# Clean up consistent and non-consistent results
consistent_cleaned <- clean_up_results(consistent_results, overlap_categories)
nonconsistent_cleaned <- clean_up_results(nonconsistent_results, overlap_categories)

# Save cleaned summary results
write.csv(consistent_cleaned, file = "Consistent_MethodComp_Experience_NC_HC_Summary.csv", row.names = TRUE)
write.csv(nonconsistent_cleaned, file = "NonConsistent_MethodComp_Experience_NC_HC_Summary.csv", row.names = TRUE)
