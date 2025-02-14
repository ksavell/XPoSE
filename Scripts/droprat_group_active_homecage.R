# Droprat for group, Active vs Homecage

library(dplyr)
library(tidyr)
library(UpSetR)
library(Seurat)

source("Scripts/Functions/single_factor_DESeq.R") 
load("all_10312024.RData")


# Run DESeq and FM for dropped iterations and no drop ---------------------

set.seed(22)

# Define target size for downsampling
target_size <- 2371  

# Function to downsample a group
downsample_group <- function(seurat_obj, group, target_size) {
  cells <- WhichCells(seurat_obj, ident = group)
  if (length(cells) > target_size) {
    return(sample(cells, target_size))  # Downsample if more than target
  } else {
    return(cells)  # Keep as is if already equal or smaller
  }
}

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
    
    # Step 2: Perform downsampling after exclusion
    groups <- levels(subset_data)
    selected_cells <- unlist(lapply(groups, function(g) downsample_group(subset_data, g, target_size)))
    
    # Step 3: Subset the downsampled Seurat object
    subset_data <- subset(subset_data, cells = selected_cells)
    
    tryCatch({
      # DESeq2 Analysis (Active vs. Homecage)
      results <- single_factor_DESeq(
        object = subset_data,
        comp_vect = c("group", "Active", "Homecage"),
        cluster = cl,
        min_cell = 1,
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
      Idents(cluster_subset) <- "group"
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

excluded_dirs <- c("Excluded_none", "Excluded_HC-1", "Excluded_HC-2", "Excluded_HC-3", "Excluded_HC-4", 
                   "Excluded_Active-1", "Excluded_Active-2", "Excluded_Active-3", "Excluded_Active-4")

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
      fm_file <- file.path(excluded_dir, paste0("FindMarkers_Results_Group_Active_Homecage_", exclusion_name, ".csv"))
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
      deseq_file <- file.path(excluded_dir, paste0("DESeq_Results_Group_Active_Homecage_", exclusion_name, ".csv"))
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
      write.csv(fm_binary_data, file.path(cluster_dir, paste0("FindMarkers_Results_binaryTally_Group_Active_Homecage_", cluster, ".csv")), row.names = FALSE)
      write.csv(fm_updn_data, file.path(cluster_dir, paste0("FindMarkers_Results_updnTally_Group_Active_Homecage_", cluster, ".csv")), row.names = FALSE)
    }
    
    if (length(deseq_binary_list) > 0) {
      deseq_binary_data <- combine_tables(deseq_binary_list)
      deseq_updn_data <- combine_tables(deseq_updn_list)
      
      # Save combined tables
      write.csv(deseq_binary_data, file.path(cluster_dir, paste0("DESeq2_Results_binaryTally_Group_Active_Homecage_", cluster, ".csv")), row.names = FALSE)
      write.csv(deseq_updn_data, file.path(cluster_dir, paste0("DESeq2_Results_updnTally_Group_Active_Homecage_", cluster, ".csv")), row.names = FALSE)
    }
    
    cat("Processing completed for cluster:", cluster, "\n")
    
  }, error = function(e) {
    cat("Error processing cluster:", cluster, "\n")
    cat("Error message:", e$message, "\n")
  })
}


# Find Consistent and nonconsistent genes ---------------------------------

for (cluster in clusters) {
  tryCatch({
    # Define the cluster directory
    cluster_dir <- file.path(getwd(), cluster)
    if (!dir.exists(cluster_dir)) {
      dir.create(cluster_dir)
    }
    
    # Load score files
    fm_score_file <- file.path(cluster_dir, paste0("FindMarkers_Results_binaryTally_Group_Active_Homecage_", cluster, ".csv"))
    fm_updnscore_file <- file.path(cluster_dir, paste0("FindMarkers_Results_updnTally_Group_Active_Homecage_", cluster, ".csv"))
    deseq_score_file <- file.path(cluster_dir, paste0("DESeq2_Results_binaryTally_Group_Active_Homecage_", cluster, ".csv"))
    deseq_updnscore_file <- file.path(cluster_dir, paste0("DESeq2_Results_updnTally_Group_Active_Homecage_", cluster, ".csv"))
    
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
    
    # Ensure Excluded_none column exists
    if (!"Excluded_none" %in% colnames(fm_score) || !"Excluded_none" %in% colnames(deseq_score)) {
      stop("Excluded_none column not found in score files!")
    }
    
    # Identify consistent (all 1's)
    fm_consistent <- rownames(fm_score)[apply(fm_score, 1, function(row) all(row == 1))]
    deseq_consistent <- rownames(deseq_score)[apply(deseq_score, 1, function(row) all(row == 1))]
    
    # Identify non-consistent genes:
    # 1. Must have at least one `1`
    # 2. Must **not** be in the consistent list
    # 3. Must **not have Excluded_none == 0**
    fm_nonconsistent <- rownames(fm_score)[apply(fm_score, 1, function(row) any(row == 1) & !all(row == 1) & row["Excluded_none"] != 0)]
    deseq_nonconsistent <- rownames(deseq_score)[apply(deseq_score, 1, function(row) any(row == 1) & !all(row == 1) & row["Excluded_none"] != 0)]
    
    # Filter updnscore files by consistent and non-consistent genes
    fm_updn_consistent <- fm_updnscore[rownames(fm_updnscore) %in% fm_consistent, , drop = FALSE]
    fm_updn_nonconsistent <- fm_updnscore[rownames(fm_updnscore) %in% fm_nonconsistent, , drop = FALSE]
    
    deseq_updn_consistent <- deseq_updnscore[rownames(deseq_updnscore) %in% deseq_consistent, , drop = FALSE]
    deseq_updn_nonconsistent <- deseq_updnscore[rownames(deseq_updnscore) %in% deseq_nonconsistent, , drop = FALSE]
    
    # Save filtered updnscore files
    fm_consistent_file <- file.path(cluster_dir, paste0("FindMarkers_Results_updnConsistent_Group_Active_Homecage_", cluster, ".csv"))
    fm_nonconsistent_file <- file.path(cluster_dir, paste0("FindMarkers_Results_updnNonConsistent_Group_Active_Homecage_", cluster, ".csv"))
    write.csv(fm_updn_consistent, fm_consistent_file, row.names = TRUE)
    write.csv(fm_updn_nonconsistent, fm_nonconsistent_file, row.names = TRUE)
    
    deseq_consistent_file <- file.path(cluster_dir, paste0("DESeq2_Results_updnConsistent_Group_Active_Homecage_", cluster, ".csv"))
    deseq_nonconsistent_file <- file.path(cluster_dir, paste0("DESeq2_Results_updnNonConsistent_Group_Active_Homecage_", cluster, ".csv"))
    write.csv(deseq_updn_consistent, deseq_consistent_file, row.names = TRUE)
    write.csv(deseq_updn_nonconsistent, deseq_nonconsistent_file, row.names = TRUE)
    
    cat("Processed consistent and non-consistent updn scores for cluster:", cluster, "\n")
    
  }, error = function(e) {
    cat("Error processing cluster:", cluster, "\n")
    cat("Error message:", e$message, "\n")
  })
}


# Tally and save ----------------------------------------------------------

# Define the working directory (main directory containing cluster folders)
base_dir <- "Droprat"

# Get the list of cluster directories
clusters <- list.dirs(path = base_dir, recursive = FALSE, full.names = FALSE)
clusters <- clusters[grepl("^Cluster_", clusters)]  # Keep only cluster directories

# Initialize results data frame
summary_results <- data.frame(
  Cluster = clusters,
  consistent_DESeq = numeric(length(clusters)),
  consistent_FM = numeric(length(clusters)),
  nonconsistent_DESeq = numeric(length(clusters)),
  nonconsistent_FM = numeric(length(clusters)),
  stringsAsFactors = FALSE
)

proportion_results <- data.frame(
  Cluster = clusters,
  consistent_DESeq = numeric(length(clusters)),
  consistent_FM = numeric(length(clusters)),
  nonconsistent_DESeq = numeric(length(clusters)),
  nonconsistent_FM = numeric(length(clusters)),
  stringsAsFactors = FALSE
)

# Loop through each cluster
for (i in seq_along(clusters)) {
  cluster <- clusters[i]
  cat("Processing cluster:", cluster, "\n")
  
  # Define file paths
  fm_consistent_file <- file.path(base_dir, cluster, paste0("FindMarkers_Results_updnConsistent_Group_Active_Homecage_", cluster, ".csv"))
  fm_nonconsistent_file <- file.path(base_dir, cluster, paste0("FindMarkers_Results_updnNonConsistent_Group_Active_Homecage_", cluster, ".csv"))
  deseq_consistent_file <- file.path(base_dir, cluster, paste0("DESeq2_Results_updnConsistent_Group_Active_Homecage_", cluster, ".csv"))
  deseq_nonconsistent_file <- file.path(base_dir, cluster, paste0("DESeq2_Results_updnNonConsistent_Group_Active_Homecage_", cluster, ".csv"))
  
  # Read and count rows for each file (if it exists)
  consistent_FM_count <- if (file.exists(fm_consistent_file)) {
    nrow(read.csv(fm_consistent_file, header = TRUE, stringsAsFactors = FALSE))
  } else {
    cat("File not found:", fm_consistent_file, "\n")
    NA
  }
  
  nonconsistent_FM_count <- if (file.exists(fm_nonconsistent_file)) {
    nrow(read.csv(fm_nonconsistent_file, header = TRUE, stringsAsFactors = FALSE))
  } else {
    cat("File not found:", fm_nonconsistent_file, "\n")
    NA
  }
  
  consistent_DESeq_count <- if (file.exists(deseq_consistent_file)) {
    nrow(read.csv(deseq_consistent_file, header = TRUE, stringsAsFactors = FALSE))
  } else {
    cat("File not found:", deseq_consistent_file, "\n")
    NA
  }
  
  nonconsistent_DESeq_count <- if (file.exists(deseq_nonconsistent_file)) {
    nrow(read.csv(deseq_nonconsistent_file, header = TRUE, stringsAsFactors = FALSE))
  } else {
    cat("File not found:", deseq_nonconsistent_file, "\n")
    NA
  }
  
  # Save counts to summary_results
  summary_results$consistent_FM[i] <- consistent_FM_count
  summary_results$nonconsistent_FM[i] <- nonconsistent_FM_count
  summary_results$consistent_DESeq[i] <- consistent_DESeq_count
  summary_results$nonconsistent_DESeq[i] <- nonconsistent_DESeq_count
  
  # Calculate proportions
  total_DESeq <- consistent_DESeq_count + nonconsistent_DESeq_count
  total_FM <- consistent_FM_count + nonconsistent_FM_count
  
  proportion_results$proportion_consistent_DESeq[i] <- if (!is.na(total_DESeq) && total_DESeq > 0) {
    consistent_DESeq_count / total_DESeq
  } else {
    NA
  }
  
  proportion_results$proportion_nonconsistent_DESeq[i] <- if (!is.na(total_DESeq) && total_DESeq > 0) {
    nonconsistent_DESeq_count / total_DESeq
  } else {
    NA
  }
  
  proportion_results$proportion_consistent_FM[i] <- if (!is.na(total_FM) && total_FM > 0) {
    consistent_FM_count / total_FM
  } else {
    NA
  }
  
  proportion_results$proportion_nonconsistent_FM[i] <- if (!is.na(total_FM) && total_FM > 0) {
    nonconsistent_FM_count / total_FM
  } else {
    NA
  }
}

# Save the summary results and proportions to separate CSV files
summary_file <- file.path(base_dir, "MethodComp_Group_Active_Homecage_Summary.csv")
proportion_file <- file.path(base_dir, "MethodComp_Group_Active_Homecage_Proportions.csv")

write.csv(summary_results, file = summary_file, row.names = FALSE)
write.csv(proportion_results, file = proportion_file, row.names = FALSE)
