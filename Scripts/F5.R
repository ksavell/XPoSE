# F5

source("Scripts/Functions/single_factor_DESeq.R")
load("~/Projects/XPoSE/all_10312024.RData")


clusters <- unique(all$cluster_name)
all$experience <- ifelse(all$group == "Homecage", "HC", "NC")

pairs_list <- list(
  c("experience","NC","HC"),
  c("group", "Non-active", "Homecage"),
  c("group", "Active", "Homecage"),
  c("group", "Active", "Non-active")
)


# flow function -----------------------------------------------------------


de_and_summary <- function(seur_obj, pair) {
  results_list <- list()  # Reset results list for each comparison
  i <- 1  # Index for list tracking
  
  for (cl in clusters) {
    # Run DESeq2 analysis
    deseq2_results <- tryCatch({
      single_factor_DESeq(object = all, comp_vect = pair, cluster = cl, min_cell = 0)
    }, error = function(e) {
      message(paste("DESeq2 failed for cluster:", cl, "with pair:", paste(pair, collapse = "_")))
      NULL
    })
    
    
    deseq2_results_tbl <- deseq2_results[["results"]]
    dds <- deseq2_results[["dds"]]
    
    if (!("padj" %in% colnames(deseq2_results_tbl)) || !("log2FoldChange" %in% colnames(deseq2_results_tbl))) {
      warning(paste("Columns 'padj' or 'log2FoldChange' missing for cluster:", cl))
      next
    }
    
    deseq2_results_tbl$padj <- as.numeric(deseq2_results_tbl$padj)
    deseq2_results_tbl$log2FoldChange <- as.numeric(deseq2_results_tbl$log2FoldChange)
    
    deseq2_results_tbl$up_score <- ifelse(!is.na(deseq2_results_tbl$padj) & deseq2_results_tbl$padj < 0.05 & deseq2_results_tbl$log2FoldChange > 0, 1, 0)
    deseq2_results_tbl$dn_score <- ifelse(!is.na(deseq2_results_tbl$padj) & deseq2_results_tbl$padj < 0.05 & deseq2_results_tbl$log2FoldChange < 0, -1, 0)
    
    results_list[[i]] <- data.frame(Category = paste(pair, collapse = "_"), Observation = cl, Value = sum(deseq2_results_tbl$up_score, na.rm = TRUE))
    results_list[[i + 1]] <- data.frame(Category = paste(pair, collapse = "_"), Observation = cl, Value = sum(deseq2_results_tbl$dn_score, na.rm = TRUE))
    i <- i + 2
    
    fn <- paste0(cl, "_", paste(pair, collapse = "_"))
    if (nrow(deseq2_results_tbl) > 0) {
      write.csv(deseq2_results_tbl, paste0(fn, ".csv"), row.names = FALSE)
      saveRDS(dds, paste0(fn,".RDS"))
    } else {
      message(paste("Skipping file save: No DESeq2 results for", fn))
    }
  }
  
  # Combine results for this pair and save
  if (length(results_list) > 0) {
    results <- do.call(rbind, results_list)
    file_name <- paste0(pair[1], "_", pair[2], "_", pair[3], ".csv")
    write.csv(results, file_name, row.names = FALSE)
  }
}



# NC vs HC ----------------------------------------------------------------
# number of cells
# Define target sizes
target_size <- 2371  # For two groups
other_target_size <- 4742  # For the third group (change as needed)

# Specify which two groups should be the same and which one is different
same_size_groups <- c("Active", "Non-active")  # Replace with actual group names
different_size_group <- "Homecage"  # Replace with actual group name

# Function to downsample a group with a specified target size
downsample_group <- function(seurat_obj, group, target_size) {
  cells <- WhichCells(seurat_obj, ident = group)
  if (length(cells) > target_size) {
    return(sample(cells, target_size))  # Downsample if more than target
  } else {
    return(cells)  # Keep as is if already equal or smaller
  }
}

# Get all identities
groups <- levels(all)

# Downsample groups with different targets
selected_cells <- unlist(lapply(groups, function(g) {
  if (g %in% same_size_groups) {
    downsample_group(all, g, target_size)
  } else if (g == different_size_group) {
    downsample_group(all, g, other_target_size)
  } else {
    WhichCells(all, ident = g)  # Keep other groups untouched
  }
}))

# Subset the Seurat object
all <- subset(all, cells = selected_cells)

# Check the new cell counts
table(all$group)

de_and_summary(all,
               pair = pairs_list[[1]])

# A vs HC and A vs NA -----------------------------------------

# Set seed 
set.seed(22)

# number of cells
target_size <- 2371

Idents(all) <- "group"

# Function to downsample a group
downsample_group <- function(seurat_obj, group, target_size) {
  cells <- WhichCells(seurat_obj, ident = group)
  if (length(cells) > target_size) {
    return(sample(cells, target_size))  # Downsample if more than target
  } else {
    return(cells)  # Keep as is if already equal or smaller
  }
}

# Get all identities
groups <- levels(all)

# Downsample each group
selected_cells <- unlist(lapply(groups, function(g) downsample_group(all, g, target_size)))

# Subset the Seurat object
all <- subset(all, cells = selected_cells)

# Check the new cell counts
table(Idents(all))

de_and_summary(all,
               pair = pairs_list[[2]])

de_and_summary(all,
               pair = pairs_list[[3]])

de_and_summary(all,
               pair = pairs_list[[4]])


#save downsampled all for Figure S4

save(all, file = "all_downsampledforde_03212025.RData")
