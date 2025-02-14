# example gene heatmaps

source("~/XPoSE/Scripts/Functions/find_DEcounts.R")

all_glut <- c("Actn1",
              "Egr3",
              "Fosb",
              "Homer1",
              "Inhba",
              "Nptx2",
              "Ppm1h",
              "Slc6a17",
              "Slc9a5",
              "Tet3")

glut_gaba <- c("Vgf","Prprn")

specific <- c("Olfm3","Ecel1", # Sst
              "Dnah2","Crh", # Pvalb
              "Efcab3","Nptx1", # CTL6
              "Pim1","Sik1", # ITL23
              "Prkcdbp", # ITL5
              "Gfra1","Lbh", # ITL6
              "Pdyn","Ctrl","Dusp1","Oacyl" # PTL5
              )
              
others <- c("Bdnf","Homer1","Scg2","Arc","Fosb")

together <- c("Actn1",
              "Egr3",
              "Fosb",
              "Homer1",
              "Inhba",
              "Nptx2",
              "Ppm1h",
              "Slc6a17",
              "Slc9a5",
              "Tet3",
              "Olfm3","Ecel1", # Sst
              "Dnah2","Crh", # Pvalb
              "Efcab3","Nptx1", # CTL6
              "Pim1","Sik1", # ITL23
              "Prkcdbp", # ITL5
              "Gfra1","Lbh", # ITL6
              "Pdyn","Ctrl","Dusp1","Oacyl", # PTL5
              "Bdnf","Homer1","Scg2","Arc","Fosb",
              "Vgf","Ptprn")

clusters <- c("ITL23","ITL5","ITL6","CTL6","PTL5","Pvalb","Sst")
# test
for (cl in clusters) {
find_DEcounts(directory = "~/Library/CloudStorage/Box-Box/RM_Projects/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/NeuroResource_January2025/DataForFigures/F5/population_degs/Group_Active_HC_2371subset", 
              cluster = cl, 
              de_path = "_group_Active_Homecage", 
              feature_list = together, 
              control_suffix = "HC")
}


# summarize ---------------------------------------------------------------

summarize_DEcounts <- function(directory, clusters) {
  
  # Initialize empty lists to store results
  log2FC_HC_list <- list()
  log2FC_Active_list <- list()
  adjpval_list <- list()
  
  # Loop through clusters and process each file
  for (cluster in clusters) {
    file_path <- file.path(directory, paste0("Log2FC_Values_", cluster, "_group_Active_Homecage.csv"))
    
    if (!file.exists(file_path)) {
      warning("File not found: ", file_path)
      next
    }
    
    data <- read.csv(file_path, row.names = 1)
    
    # Select relevant columns and order them
    log2FC_HC <- data[, grep("log2FC_HC.[1-4].Homecage", colnames(data)), drop = FALSE]
    log2FC_HC <- log2FC_HC[, order(as.numeric(gsub(".*HC.([1-4]).Homecage", "\\1", colnames(log2FC_HC))))]
    
    log2FC_Active <- data[, grep("log2FC_NC.[1-4].Active", colnames(data)), drop = FALSE]
    log2FC_Active <- log2FC_Active[, order(as.numeric(gsub(".*NC.([1-4]).Active", "\\1", colnames(log2FC_Active))))]
    
    adjpval <- data[, "adjpval", drop = FALSE]
    
    # Store results
    log2FC_HC_list[[cluster]] <- log2FC_HC
    log2FC_Active_list[[cluster]] <- log2FC_Active
    adjpval_list[[cluster]] <- adjpval
  }
  
  # Combine results across clusters
  log2FC_HC_summary <- do.call(cbind, log2FC_HC_list)
  log2FC_Active_summary <- do.call(cbind, log2FC_Active_list)
  adjpval_summary <- do.call(cbind, adjpval_list)
  
  # Save results
  write.csv(log2FC_HC_summary, file = file.path(directory, "Log2FC_summary_HC.csv"), row.names = TRUE)
  write.csv(log2FC_Active_summary, file = file.path(directory, "Log2FC_summary_Active.csv"), row.names = TRUE)
  write.csv(adjpval_summary, file = file.path(directory, "adjpval_summary.csv"), row.names = TRUE)
  
  message("Summary files saved in ", directory)
}

summarize_DEcounts("~/Library/CloudStorage/Box-Box/RM_Projects/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/NeuroResource_January2025/DataForFigures/F5/population_degs/Group_Active_HC_2371subset",
                  clusters = clusters)
