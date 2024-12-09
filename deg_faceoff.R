load("~/rstudio_projects/XPoSE/all_10312024.RData")
source("Scripts/Functions/single_factor_DESeq.R")

all$experience <- ifelse(all$group == "Homecage", "HC",
                         "NC")

deltavar <- calculate_delta_variance(
  all,
  replicate_col = "ratID",
  cell_type_col = "cluster_name",
  label_col = "experience",
  min_cells = 3,
  min_reps = 2,
  min_features = 0
)


# FindMarkers DEG ---------------------------------------------------------

# Create an empty list to store DEG results
findmarkers <- list()

# Unique clusters
Idents(all) <- "cluster_name"
clusters <- unique(Idents(all))

# Loop through each cluster and perform DEG analysis
for (cluster in clusters) {
  # Subset the Seurat object for each cluster
  sub_all <- subset(all, idents = cluster)
  Idents(sub_all) <- "experience"
  
  # Perform differential expression analysis
  findmarkers[[paste(cluster)]] <- FindMarkers(sub_all, 
                                                          ident.1 = "NC", 
                                                          ident.2 = "HC", 
                                                          min.pct = 0.1, 
                                                          logfc.threshold = 0.1)
}

# Name the list elements
names(findmarkers) <- paste(clusters, sep = " ")

for (name in names(findmarkers)) {
  write.csv(findmarkers[[name]], paste0(name, "_findmarkers.csv"))
}


# Now bring in pseudobulk and compare -------------------------------------

# Define a list to hold data frames for each cluster
pseudobulk <- list()

# Loop over each cluster to import 'gene' and 'padj' columns from experience_tally.csv
for (cl in clusters) {
  # Define the file path for experience_tally.csv in the current cluster folder
  experience_file <- file.path(cl, "experience_tally.csv")
  
  # Check if the experience_tally.csv file exists
  if (file.exists(experience_file)) {
    # Read the entire CSV file
    experience_data <- read.csv(experience_file)
    
    # Select only the 'gene' and 'padj' columns, if they exist in the data
    if(all(c("gene", "padj") %in% colnames(experience_data))) {
      selected_data <- experience_data[, c("gene", "padj", "log2FoldChange")]
      # Add the data frame to the list with the cluster name as the key
      pseudobulk[[cl]] <- selected_data
    } else {
      warning(paste("Required columns 'gene' and 'padj' not found in", experience_file))
    }
  } else {
    warning(paste("experience_tally.csv not found for cluster:", cl))
  }
}

# Merge them and compare --------------------------------------------------

comparison_results <- list()  # Initialize a list to store comparison results

for (cluster_name in names(findmarkers)) {
  if (!is.null(pseudobulk[[cluster_name]]) && !is.null(findmarkers[[cluster_name]])) {
    # Convert row names to a column for findmarkers
    deg_data <- data.frame(gene = rownames(findmarkers[[cluster_name]]), 
                           p_val_adj = as.numeric(findmarkers[[cluster_name]]$p_val_adj),
                           avg_log2FC = as.numeric(findmarkers[[cluster_name]]$avg_log2FC),
                           stringsAsFactors = FALSE)
    
    # Merge the two data frames on 'gene' using full outer join
    merged_data <- merge(pseudobulk[[cluster_name]], deg_data, by = "gene", all = TRUE)
    
        # Store the merged data in the list
    comparison_results[[cluster_name]] <- merged_data
  } else {
    warning(paste("No data available for cluster:", cluster_name, "in either pseudobulk or findmarkers lists"))
  }
}


# combine the delta variance value ----------------------------------------

# Loop through each cluster in comparison_results
for (cluster_name in names(comparison_results)) {
  # Check if the cluster exists in both deltavar and comparison_results
  if (!is.null(deltavar[[cluster_name]]) && !is.null(comparison_results[[cluster_name]])) {
    
    # Merge the DV column from deltavar into comparison_results for each cluster
    comparison_results[[cluster_name]] <- merge(comparison_results[[cluster_name]], 
                                                deltavar[[cluster_name]][c("gene", "DV")], 
                                                by = "gene", 
                                                all.x = TRUE, 
                                                all.y = FALSE)
  } else {
    warning(paste("Data for cluster", cluster_name, "is missing in either deltavar or comparison_results."))
  }
}


# plot it -----------------------------------------------------------------

library(dplyr)
library(ggplot2)

# Loop through each cluster's data frame in the list
comparison_results_transformed <- lapply(comparison_results, function(df) {
  df %>%
    mutate(p_val_adj = ifelse(is.na(p_val_adj), 1, p_val_adj),
           padj = ifelse(is.na(padj), 1, padj),
           DV = ifelse(is.na(DV), 0, DV),# Replace NA with 1
           findmarkers = -log10(p_val_adj),
           pseudobulk = -log10(padj))  # Apply -log10 transformation
})

# Handle p-values that were 0 which now are -Inf after log transformation
comparison_results_transformed <- lapply(comparison_results_transformed, function(df) {
  df$pseudobulk[is.infinite(df$pseudobulk)] <- max(df$pseudobulk, na.rm = TRUE) + 1
  return(df)
})
comparison_results_transformed <- lapply(comparison_results_transformed, function(df) {
  df$findmarkers[is.infinite(df$findmarkers)] <- max(df$findmarkers, na.rm = TRUE) + 1
  return(df)
})

# Combine all clusters into one data frame for visualization
all_data <- do.call(rbind, lapply(names(comparison_results_transformed), function(cluster_name) {
  cbind(cluster = cluster_name, comparison_results_transformed[[cluster_name]])
}))

# plots not working
# Create a scatter plot
pseudo_plot <- ggplot(all_data, aes(x = pseudobulk, y = DV, color = cluster)) +
  geom_point(alpha = 0.5) + 
  labs(title = "Scatter plot of DV vs Negative Log-transformed p-value from pseudobulk",
       x = "Negative Log10 Adjusted p-value (pseudobulk)", y = "DV Value") +
  theme_minimal() +
  facet_wrap(~ cluster)

fm_plot <- ggplot(all_data, aes(x = findmarkers, y = DV, color = cluster)) +
  geom_point(alpha = 0.5) + 
  labs(title = "Scatter plot of DV vs Negative Log-transformed p-value from findmarkers",
       x = "Negative Log10 Adjusted p-value (Findmarkers)", y = "DV Value") +
  theme_minimal() +
  facet_wrap(~ cluster)

pseudo_plot + fm_plot

ggplot(all_data, aes(x = findmarkers, y = pseudobulk, color = cluster)) +
  geom_point(alpha = 0.5) +
  labs(x = "Negative Log10 Adjusted p-value (Findmarkers)", 
       y = "Negative Log10 Adjusted p-value (Pseudobulk)") +
  theme_minimal() +
  facet_wrap(~ cluster)


# try rrho  ---------------------------------------------------------------

source("Scripts/Functions/prep_rrho.R")

library(RRHO2)

comparison_rrho <- prep_rrho_list(comparison_results)

pseudobulk_list <- list()
findmarkers_list <- list()
rrho_obj_list <- list()

## THIS IS SKIPPING SOME CLUSTERS, INVESTIGATE IF YOU WANT TO KEEP THIS ANALYSIS

for (clst in clusters) {
  # Check if the data frames are empty
  if (nrow(comparison_rrho[[clst]]) == 0) {
    warning(paste("Data frame for cluster", clst, "is empty. Skipping..."))
    next
  }
  
  # Continue processing if data is valid
  pseudobulk_list[[clst]] <- data.frame(Genes = comparison_rrho[[clst]]$gene,
                                        DDE = comparison_rrho[[clst]]$DDE_pseudobulk, stringsAsFactors = FALSE)
  findmarkers_list[[clst]] <- data.frame(Genes = comparison_rrho[[clst]]$gene,
                                         DDE = comparison_rrho[[clst]]$DDE_findmarkers, stringsAsFactors = FALSE)
  
  if (nrow(pseudobulk_list[[clst]]) > 0 && nrow(findmarkers_list[[clst]]) > 0) {
    rrho_obj_list[[clst]] <- RRHO2_initialize(pseudobulk_list[[clst]], findmarkers_list[[clst]], stepsize = 10,
                                              labels = c("pseudobulk", "findmarkers"), log10.ind = TRUE)
  } else {
    warning(paste("One of the data frames for cluster", clst, "is empty or invalid. Skipping..."))
  }
}

# Visualize the heatmap, here 

RRHO2_heatmap(rrho_obj_list[["ITL23"]])

# Visualize the Venn Diagram

RRHO2_vennDiagram(rrho_obj_list[["ITL5"]], type = "uu") # type is pattern like up up 
RRHO2_vennDiagram(rrho_obj_list[["ITL5"]], type = "dd")


# coefficient of variation deseq2 ------------------------------------------------
#you generated the deseq2 objects in the enrichment.deg.r file by adding keep_dds = T to the no downsample run

# Define a list to hold the paths for DESeq2 RDS files
rds_files <- paste0("DESeq2_", clusters, ".rds")

# Define lists to hold CV data frames
cv_data_frames <- list() 

# Loop to create cv_data_frames and handle NAs correctly
for (i in seq_along(rds_files)) {
  file_path <- rds_files[i]
  cluster_name <- clusters[i]
  
  # Check if the RDS file exists
  if (file.exists(file_path)) {
    dds <- readRDS(file_path)
    
    if ("DESeqDataSet" %in% class(dds)) {
      norm_counts <- counts(dds, normalized = TRUE)
      gene_means <- rowMeans(norm_counts)
      gene_sds <- apply(norm_counts, 1, sd)
      gene_cv <- gene_sds / gene_means
      
      # Extract results including log2FoldChange and padj
      deseq_results <- results(dds)
      
      # Create data frame
      cv_df <- data.frame(
        Genes = rownames(deseq_results),
        log2FoldChange = deseq_results$log2FoldChange,
        padj = deseq_results$padj,
        CoefficientOfVariation = gene_cv[rownames(deseq_results)],
        stringsAsFactors = FALSE
      )
      
      # Remove rows with NA in log2FoldChange, padj, or CoefficientOfVariation
      cv_df <- na.omit(cv_df)
      
      cv_data_frames[[cluster_name]] <- cv_df
    } else {
      warning(paste("DESeq2 analysis may not be complete for", cluster_name))
    }
  } else {
    warning(paste("RDS file not found for", cluster_name))
  }
}


# now plot it
plot_list <- list()  # To store ggplot objects for each cluster

for (cluster in names(cv_data_frames)) {
  # Retrieve the data frame from the list
  df <- cv_data_frames[[cluster]]
  
  # Generate a plot for this cluster
  p <- ggplot(df, aes(x = log2FoldChange, y = CoefficientOfVariation, color = padj < 0.05)) +
    geom_point(alpha = 0.6) +  # Use point plots with slight transparency
    scale_color_manual(values = c("gray", "red"), labels = c("Not Significant", "Significant")) +
    labs(title = paste("CV vs. log2FoldChange for Cluster", cluster),
         x = "log2FoldChange",
         y = "Coefficient of Variation (CV)",
         color = "DEG Significance") +
    theme_classic() +
    theme(legend.position = "right")+
    xlim(-5,5)+
    ylim(0,3)
  
  # Add plot to the list
  plot_list[[cluster]] <- p
}

# Optional: Combine all plots into a single grid layout
combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 2))

# Coefficient of Variation FindMarkers ------------------------------------

# Loop through each cluster's DEG results
library(Seurat)

for (cluster in names(findmarkers)) {
  # Get the DEG results for the current cluster
  degs_results <- findmarkers[[cluster]]
  
  # Get gene names from DEGs
  deg_genes <- rownames(degs_results)
  
  # Check if there are DEGs to process
  if (length(deg_genes) > 0) {
    # Subset the Seurat object for cells in the current cluster
    sub_all <- subset(all, idents = cluster)
    
    # Fetch normalized expression data for these genes from cells within the cluster
    expression_data <- FetchData(sub_all, vars = deg_genes)
    
    # Calculate mean and SD for each gene across cells in the cluster
    gene_means <- colMeans(expression_data)
    gene_sds <- apply(expression_data, 2, sd)
    
    # Calculate coefficient of variation, safely handle division by zero
    gene_cv <- sapply(seq_along(gene_means), function(i) {
      if (gene_means[i] != 0) {
        gene_sds[i] / gene_means[i]
      } else {
        NA  # Assign NA where mean is zero to avoid division by zero
      }
    })
    
    # Add CV to the DEGs results data frame
    degs_results$CV <- gene_cv[match(rownames(degs_results), names(gene_cv))]
    
    # Optionally, save updated DEG results with CV back to findmarkers list
    findmarkers[[cluster]] <- degs_results
    
    # Save updated DEG results to CSV
    write.csv(degs_results, paste0(cluster, "_findmarkers_with_CV.csv"), row.names = TRUE)
  } else {
    cat("No DEGs for cluster:", cluster, "\n")
  }
}


#visualize
library(gridExtra)  # For arranging multiple plots

# Assuming 'findmarkers' is the list containing DEG results for each cluster
plot_list <- list()  # To store ggplot objects for each cluster

for (cluster in names(findmarkers)) {
  # Retrieve the data frame from the list
  df <- findmarkers[[cluster]]
  
  # Generate a plot for this cluster
  p <- ggplot(df, aes(x = avg_log2FC, y = CV, color = p_val_adj < 0.05)) +
    geom_point(alpha = 0.6) +  # Use point plots with slight transparency
    scale_color_manual(values = c("gray", "red"), labels = c("Not Significant", "Significant")) +
    labs(title = paste("CV vs. log2FoldChange for Cluster", cluster),
         x = "log2FoldChange",
         y = "Coefficient of Variation (CV)",
         color = "DEG Significance") +
    theme_classic() +
    theme(legend.position = "right")+
    xlim(-7,7) +
    ylim(0,5)
  
  # Add plot to the list
  plot_list[[cluster]] <- p
}

# Optional: Combine all plots into a single grid layout
combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 2))



# Scatter plot to compare CV between methods ------------------------------

cv_method1 <- data.frame(Genes = cv_data_frames[["ITL23"]]$Genes,
                         CV = cv_data_frames[["ITL23"]]$CoefficientOfVariation,
                         padj = cv_data_frames[["ITL23"]]$padj)
cv_method2 <- data.frame(Genes = rownames(findmarkers[["ITL23"]]),
                         CV = findmarkers[["ITL23"]]$CV,
                         padj = findmarkers[["ITL23"]]$p_val_adj)

# Assume cv_method1 and cv_method2 are prepared and have columns `Genes` and `CV`
merged_cv <- merge(cv_method1, cv_method2, by = "Genes", suffixes = c("_pseudobulk", "_findmarkers"))

# Enhancing the filtering and categorization for significance
significant_merged_cv$Significance_Category <- with(significant_merged_cv, ifelse(
  padj_pseudobulk < 0.05 & padj_findmarkers < 0.05, "Significant in Both",
  ifelse(padj_pseudobulk < 0.05, "Unique to Method 1",
         ifelse(padj_findmarkers < 0.05, "Unique to Method 2", "Not Significant")
  )))

# You might want to remove the "Not Significant" category if it's not needed
significant_merged_cv <- subset(significant_merged_cv, Significance_Category != "Not Significant")


# Scatter plot to compare CVs
ggplot(significant_merged_cv, aes(x = CV_pseudobulk, y = CV_findmarkers, color = Significance_Category)) +
  geom_point(alpha = 0.6) +  # Use point plots with slight transparency
  scale_color_manual(values = c("red", "blue", "green"),
                     labels = c("Significant in Both", "Unique to Method 1", "Unique to Method 2")) +
  labs(title = "CV Comparison by Significance Category",
       x = "Coefficient of Variation (Pseudobulk)",
       y = "Coefficient of Variation (FindMarkers)",
       color = "Category") +
  theme_minimal() +
  theme(legend.position = "right") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")

