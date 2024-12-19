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
    degs_results$CV_FM <- gene_cv[match(rownames(degs_results), names(gene_cv))]
    
    # Optionally, save updated DEG results with CV back to findmarkers list
    findmarkers[[cluster]] <- degs_results
    
    # Save updated DEG results to CSV
    write.csv(degs_results, paste0(cluster, "_findmarkers_with_CV.csv"), row.names = TRUE)
  } else {
    cat("No DEGs for cluster:", cluster, "\n")
  }
}

# Now bring in pseudobulk and compare -------------------------------------

for (cl in clusters) {
  tryCatch({
    # Run DESeq2 for Non-active vs Homecage
    deseq2_results <- single_factor_DESeq(object = all,
                                          comp_vect = c("experience", "NC", "HC"),
                                          cluster = cl,
                                          min_cell = 1,
                                          keep_dds = TRUE)
    deseq2_results$score_column <- ifelse(deseq2_results$padj < 0.05, 1, 0)
    
    # Save the results in the cluster-specific directory
    deg_file <- file.path(paste0(cl, "_experience_NC_HC.csv"))
    write.csv(deseq2_results, file = deg_file)
  }, error = function(e) {
    message(paste("Error for cluster:", cl, "-", e$message))
  })
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


# Tally up and down for each method  ------------------------------------------------------------

# Initialize an empty data frame to store the tally results
comparison_tally <- data.frame(
  Cluster = character(),
  Method1_Upregulated = integer(),
  Method1_Downregulated = integer(),
  Method1_Neutral = integer(),
  Method2_Upregulated = integer(),
  Method2_Downregulated = integer(),
  Method2_Neutral = integer(),
  stringsAsFactors = FALSE
)

# Loop through each element in the list
for (cluster_name in names(comparison_results)) {
  # Extract the data frame for the current cluster
  df <- comparison_results[[cluster_name]]
  
  # Tally for Method 1
  method1_upregulated <- sum(df$log2FoldChange > 0 & df$padj < 0.05, na.rm = TRUE)
  method1_downregulated <- sum(df$log2FoldChange < 0 & df$padj < 0.05, na.rm = TRUE)
  method1_neutral <- sum(df$padj >= 0.05, na.rm = TRUE)
  
  # Tally for Method 2
  method2_upregulated <- sum(df$avg_log2FC > 0 & df$p_val_adj < 0.05, na.rm = TRUE)
  method2_downregulated <- sum(df$avg_log2FC < 0 & df$p_val_adj < 0.05, na.rm = TRUE)
  method2_neutral <- sum(df$p_val_adj >= 0.05, na.rm = TRUE)
  
  # Append the results to the tally data frame
  comparison_tally <- rbind(comparison_tally, data.frame(
    Cluster = cluster_name,
    Method1_Upregulated = method1_upregulated,
    Method1_Downregulated = method1_downregulated,
    Method1_Neutral = method1_neutral,
    Method2_Upregulated = method2_upregulated,
    Method2_Downregulated = method2_downregulated,
    Method2_Neutral = method2_neutral,
    stringsAsFactors = FALSE
  ))
}

# Save the results as a CSV
output_file <- "Method_Comparison_Tally.csv"
write.csv(comparison_tally, file = output_file, row.names = FALSE)

# Print the tally results
print(comparison_tally)


# Method head to head matrix, both 0.05 threshold -------------------------

# Initialize an empty data frame to store the tally results
comparison_tally <- data.frame(
  Cluster = character(),
  Both_Downregulated = integer(),
  Neutral_1_Down_2 = integer(),
  Up_1_Down_2 = integer(),
  Down_1_Neutral_2 = integer(),
  Both_Neutral = integer(),
  Up_1_Neutral_2 = integer(),
  Down_1_Up_2 = integer(),
  Neutral_1_Up_2 = integer(),
  Both_Upregulated = integer(),
  stringsAsFactors = FALSE
)

# Loop through each element in the list
for (cluster_name in names(comparison_results)) {
  # Extract the data frame for the current cluster
  df <- comparison_results[[cluster_name]]
  
  # Exclude rows where both methods have NA values
  valid_rows <- !(is.na(df$log2FoldChange) & is.na(df$avg_log2FC) & is.na(df$padj) & is.na(df$p_val_adj))
  df <- df[valid_rows, ]
  
  # Replace NA with neutral for Method 1
  method1_down <- ifelse(is.na(df$log2FoldChange) | is.na(df$padj), FALSE, df$log2FoldChange < 0 & df$padj < 0.05)
  method1_neutral <- ifelse(is.na(df$padj), TRUE, df$padj >= 0.05)
  method1_up <- ifelse(is.na(df$log2FoldChange) | is.na(df$padj), FALSE, df$log2FoldChange > 0 & df$padj < 0.05)
  
  # Replace NA with neutral for Method 2
  method2_down <- ifelse(is.na(df$avg_log2FC) | is.na(df$p_val_adj), FALSE, df$avg_log2FC < 0 & df$p_val_adj < 0.05)
  method2_neutral <- ifelse(is.na(df$p_val_adj), TRUE, df$p_val_adj >= 0.05)
  method2_up <- ifelse(is.na(df$avg_log2FC) | is.na(df$p_val_adj), FALSE, df$avg_log2FC > 0 & df$p_val_adj < 0.05)
  
  # Tally the combinations
  both_downregulated <- sum(method1_down & method2_down, na.rm = TRUE)
  neutral_1_down_2 <- sum(method1_neutral & method2_down, na.rm = TRUE)
  up_1_down_2 <- sum(method1_up & method2_down, na.rm = TRUE)
  down_1_neutral_2 <- sum(method1_down & method2_neutral, na.rm = TRUE)
  both_neutral <- sum(method1_neutral & method2_neutral, na.rm = TRUE)
  up_1_neutral_2 <- sum(method1_up & method2_neutral, na.rm = TRUE)
  down_1_up_2 <- sum(method1_down & method2_up, na.rm = TRUE)
  neutral_1_up_2 <- sum(method1_neutral & method2_up, na.rm = TRUE)
  both_upregulated <- sum(method1_up & method2_up, na.rm = TRUE)
  
  # Append the results to the tally data frame
  comparison_tally <- rbind(comparison_tally, data.frame(
    Cluster = cluster_name,
    Both_Downregulated = both_downregulated,
    Neutral_1_Down_2 = neutral_1_down_2,
    Up_1_Down_2 = up_1_down_2,
    Down_1_Neutral_2 = down_1_neutral_2,
    Both_Neutral = both_neutral,
    Up_1_Neutral_2 = up_1_neutral_2,
    Down_1_Up_2 = down_1_up_2,
    Neutral_1_Up_2 = neutral_1_up_2,
    Both_Upregulated = both_upregulated,
    stringsAsFactors = FALSE
  ))
}

# Save the results as a CSV
output_file <- "Method_Comparison_Matrix.csv"
write.csv(comparison_tally, file = output_file, row.names = FALSE)

# Print the tally results
print(comparison_tally)


# Get gene names for certain combos ---------------------------------------

# Initialize a list to store gene names for each combination
gene_lists <- list(
  Up_1_Neutral_2 = list(),
  Neutral_1_Down_2 = list(),
  Neutral_1_Up_2 = list()
)

# Loop through each element in the list
for (cluster_name in names(comparison_results)) {
  # Extract the data frame for the current cluster
  df <- comparison_results[[cluster_name]]
  
  # Exclude rows where both methods have NA values
  valid_rows <- !(is.na(df$log2FoldChange) & is.na(df$avg_log2FC) & is.na(df$padj) & is.na(df$p_val_adj))
  df <- df[valid_rows, ]
  
  # Replace NA with neutral for Method 1
  method1_down <- ifelse(is.na(df$log2FoldChange) | is.na(df$padj), FALSE, df$log2FoldChange < 0 & df$padj < 0.05)
  method1_neutral <- ifelse(is.na(df$padj), TRUE, df$padj >= 0.05)
  method1_up <- ifelse(is.na(df$log2FoldChange) | is.na(df$padj), FALSE, df$log2FoldChange > 0 & df$padj < 0.05)
  
  # Replace NA with neutral for Method 2
  method2_down <- ifelse(is.na(df$avg_log2FC) | is.na(df$p_val_adj), FALSE, df$avg_log2FC < 0 & df$p_val_adj < 0.05)
  method2_neutral <- ifelse(is.na(df$p_val_adj), TRUE, df$p_val_adj >= 0.05)
  method2_up <- ifelse(is.na(df$avg_log2FC) | is.na(df$p_val_adj), FALSE, df$avg_log2FC > 0 & df$p_val_adj < 0.05)
  
  # Extract gene names for each combination
  gene_lists$Up_1_Neutral_2[[cluster_name]] <- df$gene[method1_up & method2_neutral]
  gene_lists$Neutral_1_Down_2[[cluster_name]] <- df$gene[method1_neutral & method2_down]
  gene_lists$Neutral_1_Up_2[[cluster_name]] <- df$gene[method1_neutral & method2_up]
}

# Save gene lists as CSV files for each combination
for (combination in names(gene_lists)) {
  gene_data <- data.frame(
    Cluster = rep(names(gene_lists[[combination]]), 
                  times = sapply(gene_lists[[combination]], length)),
    Gene = unlist(gene_lists[[combination]], use.names = FALSE),
    stringsAsFactors = FALSE
  )
  write.csv(gene_data, file = paste0(combination, "_Gene_List.csv"), row.names = FALSE)
}

# Print the first few genes for each combination
lapply(gene_lists, function(x) lapply(x, head))


# calculate coeff variation for FindMarkers -------------------------------

# Initialize a list to store the results
gene_data_extracted <- list()

# Loop through each condition in gene_lists
for (condition in names(gene_lists)) {
  # Initialize a data frame to store results for this condition
  condition_results <- data.frame(
    Cluster = character(),
    Gene = character(),
    p_val_adj = numeric(),
    avg_log2FC = numeric(),
    CV_FM = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each cluster in the condition
  for (cluster_name in names(gene_lists[[condition]])) {
    # Extract the gene list for this cluster and condition
    gene_list <- gene_lists[[condition]][[cluster_name]]
    
    # Extract the corresponding cluster data from findmarkers
    cluster_data <- findmarkers[[cluster_name]]
    
    # Check if row names (genes) and required columns exist
    if (!is.null(cluster_data) && all(c("p_val_adj", "avg_log2FC", "CV_FM") %in% colnames(cluster_data))) {
      # Add row names as a "Gene" column
      cluster_data$Gene <- rownames(cluster_data)
      
      # Filter the cluster_data for genes in the current gene_list
      filtered_data <- cluster_data[cluster_data$Gene %in% gene_list, c("Gene", "p_val_adj", "avg_log2FC", "CV_FM")]
      
      # Add cluster name to the filtered data if it is non-empty
      if (nrow(filtered_data) > 0) {
        filtered_data <- cbind(Cluster = cluster_name, filtered_data)
        # Append to condition results
        condition_results <- rbind(condition_results, filtered_data)
      }
    } else {
      message(paste("Skipping cluster:", cluster_name, "due to missing required columns or data."))
    }
  }
  
  # Store the results for this condition
  gene_data_extracted[[condition]] <- condition_results
}

# Save the results to CSV files for each condition
for (condition in names(gene_data_extracted)) {
  write.csv(gene_data_extracted[[condition]], file = paste0(condition, "_Extracted_Gene_Data.csv"), row.names = FALSE)
}

# Print a summary of the results
lapply(gene_data_extracted, head)


# Plotting dispersion vs mean normalized counts in DESeq2 for a category -----------------


# plotting dispersion vs mean normalized counts
# Extract dispersion estimates from the DESeq2 object
dispersion_data <- data.frame(
  gene = rownames(DESeq2_ITL23),
  baseMean = mcols(DESeq2_ITL23)$baseMean,
  dispersion = mcols(DESeq2_ITL23)$dispGeneEst
)

# Extract genes for the specific cluster (ITL23) across all conditions
ITL23_Neutral_1_Up_2 <- dispersion_data[dispersion_data$gene %in% gene_lists[["Neutral_1_Up_2"]][["ITL23"]], ]
ITL23_Neutral_1_Down_2 <- dispersion_data[dispersion_data$gene %in% gene_lists[["Neutral_1_Down_2"]][["ITL23"]], ]


# Base plot: dispersion_data
base_plot <- ggplot(dispersion_data, aes(x = baseMean, y = dispersion)) +
  geom_point(color = "gray", alpha = 0.5) +  # Background points in gray
  scale_x_log10() +  # Log scale for baseMean
  scale_y_log10() +  # Log scale for dispersion
  labs(
    title = "Base Mean vs. Dispersion with ITL23 Overlay",
    x = "Base Mean (log scale)",
    y = "Dispersion (log scale)"
  ) +
  theme_minimal()

# Add ITL23_Neutral_1_Up_2 points
overlay_plot <- base_plot +
  geom_point(data = ITL23_Neutral_1_Up_2, aes(x = baseMean, y = dispersion), 
             color = "blue", alpha = 0.7, size = 2) +
  geom_text(data = ITL23_Neutral_1_Up_2, aes(x = baseMean, y = dispersion, label = gene),
            hjust = 0.5, vjust = -0.5, size = 3, color = "blue") +  # Add blue labels
  
  # Add ITL23_Neutral_1_Down_2 points
  geom_point(data = ITL23_Neutral_1_Down_2, aes(x = baseMean, y = dispersion), 
             color = "red", alpha = 0.7, size = 2) +
  geom_text(data = ITL23_Neutral_1_Down_2, aes(x = baseMean, y = dispersion, label = gene),
            hjust = 0.5, vjust = -0.5, size = 3, color = "red")  # Add red labels

# Display the plot
print(overlay_plot)


