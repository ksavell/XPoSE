# Figure S4
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)


# FS4A --------------------------------------------------------------------

# wilcoxon data is saved under no_exclusion folder in FS4C output

# FS4B --------------------------------------------------------------------

# leave one out analysis
load("~/Projects/XPoSE/input/combined_annotated_07202026.RData")
source("Scripts/Functions/single_factor_DESeq.R")
date_stamp <- format(Sys.Date(), "_%m%d%Y")   

set.seed(22)
all <- obj

# List of all unique rat IDs
all_rats <- unique(all$ratID)
clusters <- unique(all$cluster_name)

if (inherits(all[["RNA"]], "Assay5")){
  all[["RNA"]] <- JoinLayers(all[["RNA"]])
}
DefaultAssay(all) <- "RNA"


# List of all unique rat IDs
all_rats <- unique(all$ratID)
clusters <- unique(all$cluster_name)

# Collect cluster directories so loop 2 can find them
cluster_dirs <- character(0)

# Running skip log (insufficient-replicate skips, distinct from real errors)
skip_log   <- data.frame(cluster = character(), excluded_rat = character(),
                         reason = character(), stringsAsFactors = FALSE)
# Running error log (genuine failures)
error_log  <- data.frame(cluster = character(), excluded_rat = character(),
                         reason = character(), stringsAsFactors = FALSE)

all_results <- list()

for (cl in clusters) {
  
  # Create a directory for the cluster
  cluster_dir <- paste0("Cluster_", cl)
  if (!dir.exists(cluster_dir)) {
    dir.create(cluster_dir)
  }
  cluster_dirs <- union(cluster_dirs, cluster_dir)
  
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
    Idents(subset_data) <- "group"
    
    tryCatch({
      # DESeq2 Analysis (Active vs. Homecage)
      results <- single_factor_DESeq(
        object    = subset_data,
        comp_vect = c("group", "active", "all"),
        cluster   = cl,
        min_cell  = 9,
        min_rat   = 2,
        keep_dds  = TRUE
      )
      
      deseq_results <- results$results
      dds <- results$dds
      
      # Ensure valid row names for deseq_results
      if (is.null(rownames(deseq_results))) {
        rownames(deseq_results) <- paste0("Gene_", seq_len(nrow(deseq_results)))
      } else if (any(is.na(rownames(deseq_results)))) {
        rownames(deseq_results)[is.na(rownames(deseq_results))] <-
          paste0("Gene_", which(is.na(rownames(deseq_results))))
      }
      
      # Add score columns
      deseq_results$score <- ifelse(deseq_results$padj < 0.05, 1, 0)
      deseq_results$score_updn <- ifelse(
        deseq_results$padj < 0.05 & deseq_results$log2FoldChange > 0, 1,
        ifelse(deseq_results$padj < 0.05 & deseq_results$log2FoldChange < 0, 2, 0))
      
      # Save DESeq2 results
      result_file <- file.path(rat_dir, paste0(
        "DESeq_Results_Group_Active_Homecage_Excluded_",
        excluded_rat, date_stamp, ".csv"))
      write.csv(deseq_results, result_file, row.names = TRUE)
      
      result_dds <- file.path(rat_dir, paste0(
        "DESeq2_Group_Active_Homecage_Excluded_",
        excluded_rat, date_stamp, ".rds"))
      saveRDS(dds, file = result_dds, compress = FALSE)
      
      # FindMarkers Analysis (Seurat object)
      # `all` is already joined + downcast to v4, so subset_data inherits
      # a clean v4 assay; a plain subset here is safe.
      cluster_subset <- subset(subset_data, cluster_name == cl)
      Idents(cluster_subset) <- "group"
      fm_results <- FindMarkers(cluster_subset,
                                ident.1 = "active",
                                ident.2 = "all")
      fm_results$score <- ifelse(fm_results$p_val_adj < 0.05, 1, 0)
      fm_results$score_updn <- ifelse(
        fm_results$p_val_adj < 0.05 & fm_results$avg_log2FC > 0, 1,
        ifelse(fm_results$p_val_adj < 0.05 & fm_results$avg_log2FC < 0, 2, 0))
      
      fm_file <- file.path(rat_dir, paste0(
        "FindMarkers_Results_Group_Active_Homecage_Excluded_",
        excluded_rat, date_stamp, ".csv"))
      write.csv(fm_results, fm_file, row.names = TRUE)
      
    }, error = function(e) {
      # Distinguish an insufficient-replicate SKIP from a real error via
      # the "SKIP:" prefix set inside single_factor_DESeq().
      if (grepl("^SKIP:", e$message)) {
        cat("SKIPPED cluster:", cl, "| excluded rat:", excluded_rat,
            "|", e$message, "\n")
        skip_log[nrow(skip_log) + 1, ] <<-
          list(cl, excluded_rat, e$message)
      } else {
        cat("ERROR in cluster:", cl, "| excluded rat:", excluded_rat,
            "\n  ", e$message, "\n")
        error_log[nrow(error_log) + 1, ] <<-
          list(cl, excluded_rat, e$message)
      }
    })
  }
}

# Write skip / error logs
if (nrow(skip_log) > 0) {
  write.csv(skip_log,
            paste0("Skipped_log", date_stamp, ".csv"), row.names = FALSE)
}
if (nrow(error_log) > 0) {
  write.csv(error_log,
            paste0("Error_log", date_stamp, ".csv"), row.names = FALSE)
}

build_combined <- function(cluster_dir, prefix) {
  # prefix: "DESeq" or "FindMarkers"
  # Find all per-exclusion result files for this method in the cluster
  pattern <- paste0("^", prefix,
                    "_Results_Group_Active_Homecage_Excluded_.*\\.csv$")
  excl_dirs <- list.dirs(cluster_dir, recursive = FALSE)
  excl_dirs <- excl_dirs[grepl("Excluded_", basename(excl_dirs))]
  
  score_list <- list()
  for (d in excl_dirs) {
    f <- list.files(d, pattern = pattern, full.names = TRUE)
    if (length(f) == 0) next
    f <- f[1]
    
    # Column name = the exclusion label (e.g. "Excluded_none", "Excluded_HC-1")
    excl_label <- basename(d)
    
    df <- read.csv(f, row.names = 1, check.names = FALSE)
    if (!("score_updn" %in% colnames(df))) next
    
    s <- data.frame(gene = rownames(df),
                    score = df$score_updn,
                    stringsAsFactors = FALSE)
    colnames(s)[2] <- excl_label
    score_list[[excl_label]] <- s
  }
  
  if (length(score_list) == 0) return(NULL)
  
  # Full outer join across all exclusions on gene
  combined <- Reduce(function(x, y) full_join(x, y, by = "gene"), score_list)
  combined <- combined %>% column_to_rownames("gene")
  
  # Genes absent from an exclusion (filtered out there) -> 0 (not significant)
  combined[is.na(combined)] <- 0
  combined
}

for (cluster_dir in cluster_dirs) {
  cluster_name <- basename(cluster_dir)
  
  deseq_combined <- build_combined(cluster_dir, "DESeq")
  if (!is.null(deseq_combined)) {
    write.csv(deseq_combined,
              file.path(cluster_dir, paste0(
                "DESeq_combined_results_", cluster_name,
                date_stamp, ".csv")),
              row.names = TRUE)
  }
  
  fm_combined <- build_combined(cluster_dir, "FindMarkers")
  if (!is.null(fm_combined)) {
    write.csv(fm_combined,
              file.path(cluster_dir, paste0(
                "FindMarkers_combined_results_", cluster_name,
                date_stamp, ".csv")),
              row.names = TRUE)
  }
}

for (cluster_dir in cluster_dirs) {
  cluster_name <- basename(cluster_dir)
  
  tryCatch({
    # Locate the combined tables written above (date-stamped)
    fm_updnscore_file <- list.files(
      cluster_dir,
      pattern = paste0("^FindMarkers_combined_results_",
                       cluster_name, ".*\\.csv$"),
      full.names = TRUE)
    deseq_updnscore_file <- list.files(
      cluster_dir,
      pattern = paste0("^DESeq_combined_results_",
                       cluster_name, ".*\\.csv$"),
      full.names = TRUE)
    
    have_fm    <- length(fm_updnscore_file) > 0
    have_deseq <- length(deseq_updnscore_file) > 0
    
    if (have_deseq) {
      deseq_updnscore <- read.csv(deseq_updnscore_file[1],
                                  row.names = 1, check.names = FALSE)
    }
    if (have_fm) {
      fm_updnscore <- read.csv(fm_updnscore_file[1],
                               row.names = 1, check.names = FALSE)
    }
    
    # Consistent (all 1 or all 2 across every exclusion)
    # Non-consistent: appears as val in at least one exclusion but not all,
    # and is not zero in the Excluded_none baseline.
    nonconsist_logic <- function(x, val) {
      any(x == val) & !all(x == val) &
        (!"Excluded_none" %in% names(x) || x["Excluded_none"] != 0)
    }
    
    save_filtered <- function(df, genes, filename) {
      out <- df[rownames(df) %in% genes, , drop = FALSE]
      if (nrow(out) > 0) {
        write.csv(out, filename, row.names = TRUE)
      } else {
        cat("Saving empty file:", basename(filename), "\n")
        write.csv(data.frame(Gene = character(), Score = integer()),
                  filename, row.names = FALSE)
      }
    }
    
    # ---- DESeq consistency ----
    if (have_deseq) {
      deseq_consistent_up   <- rownames(deseq_updnscore)[
        apply(deseq_updnscore, 1, \(x) all(x == 1))]
      deseq_consistent_down <- rownames(deseq_updnscore)[
        apply(deseq_updnscore, 1, \(x) all(x == 2))]
      deseq_nonconsistent_up   <- rownames(deseq_updnscore)[
        apply(deseq_updnscore, 1, nonconsist_logic, val = 1)]
      deseq_nonconsistent_down <- rownames(deseq_updnscore)[
        apply(deseq_updnscore, 1, nonconsist_logic, val = 2)]
      
      save_filtered(deseq_updnscore, deseq_consistent_up,
                    file.path(cluster_dir, paste0(
                      "DESeq_ConsistentUp_Cluster_", cluster_name,
                      date_stamp, ".csv")))
      save_filtered(deseq_updnscore, deseq_consistent_down,
                    file.path(cluster_dir, paste0(
                      "DESeq_ConsistentDown_Cluster_", cluster_name,
                      date_stamp, ".csv")))
      save_filtered(deseq_updnscore, deseq_nonconsistent_up,
                    file.path(cluster_dir, paste0(
                      "DESeq_NonConsistentUp_Cluster_", cluster_name,
                      date_stamp, ".csv")))
      save_filtered(deseq_updnscore, deseq_nonconsistent_down,
                    file.path(cluster_dir, paste0(
                      "DESeq_NonConsistentDown_Cluster_", cluster_name,
                      date_stamp, ".csv")))
    }
    
    # ---- FindMarkers consistency ----
    if (have_fm) {
      fm_consistent_up   <- rownames(fm_updnscore)[
        apply(fm_updnscore, 1, \(x) all(x == 1))]
      fm_consistent_down <- rownames(fm_updnscore)[
        apply(fm_updnscore, 1, \(x) all(x == 2))]
      fm_nonconsistent_up   <- rownames(fm_updnscore)[
        apply(fm_updnscore, 1, nonconsist_logic, val = 1)]
      fm_nonconsistent_down <- rownames(fm_updnscore)[
        apply(fm_updnscore, 1, nonconsist_logic, val = 2)]
      
      save_filtered(fm_updnscore, fm_consistent_up,
                    file.path(cluster_dir, paste0(
                      "FindMarkers_ConsistentUp_Cluster_", cluster_name,
                      date_stamp, ".csv")))
      save_filtered(fm_updnscore, fm_consistent_down,
                    file.path(cluster_dir, paste0(
                      "FindMarkers_ConsistentDown_Cluster_", cluster_name,
                      date_stamp, ".csv")))
      save_filtered(fm_updnscore, fm_nonconsistent_up,
                    file.path(cluster_dir, paste0(
                      "FindMarkers_NonConsistentUp_Cluster_", cluster_name,
                      date_stamp, ".csv")))
      save_filtered(fm_updnscore, fm_nonconsistent_down,
                    file.path(cluster_dir, paste0(
                      "FindMarkers_NonConsistentDown_Cluster_", cluster_name,
                      date_stamp, ".csv")))
    }
    
    cat("Processed consistency filtering for cluster:", cluster_name, "\n")
    
  }, error = function(e) {
    cat("Error processing cluster:", cluster_name, "\n")
    cat("Message:", e$message, "\n")
  })
}

count_if_exists_rows <- function(filepath, filter_val, consistent = TRUE) {
  if (!file.exists(filepath)) return(0)
  df <- read.csv(filepath, row.names = 1, check.names = FALSE)
  if (nrow(df) == 0) return(0)
  if (consistent) {
    sum(apply(df, 1, function(x) all(x == filter_val)))
  } else {
    sum(apply(df, 1, function(x) any(x == filter_val)))
  }
}


# --- Summary 1: exclude=none up/down tally, DESeq2 + FindMarkers, per cluster
# Reads the Excluded_none result files and counts score_updn (1 = up, 2 = down).
summary_updown <- data.frame(
  cluster    = character(),
  DESeq_Up   = integer(), DESeq_Down = integer(),
  FM_Up      = integer(), FM_Down    = integer(),
  stringsAsFactors = FALSE
)

tally_updn <- function(filepath) {
  # returns c(up, down); zeros if file missing / empty
  if (length(filepath) == 0 || !file.exists(filepath[1])) return(c(0L, 0L))
  df <- read.csv(filepath[1], row.names = 1, check.names = FALSE)
  if (!("score_updn" %in% colnames(df)) || nrow(df) == 0) return(c(0L, 0L))
  c(sum(df$score_updn == 1, na.rm = TRUE),
    sum(df$score_updn == 2, na.rm = TRUE))
}

for (cluster_dir in cluster_dirs) {
  cluster_name <- basename(cluster_dir)
  none_dir <- file.path(cluster_dir, "Excluded_none")
  
  deseq_file <- list.files(none_dir,
                           pattern = "^DESeq_Results_Group_Active_Homecage_Excluded_none.*\\.csv$",
                           full.names = TRUE)
  fm_file <- list.files(none_dir,
                        pattern = "^FindMarkers_Results_Group_Active_Homecage_Excluded_none.*\\.csv$",
                        full.names = TRUE)
  
  d <- tally_updn(deseq_file)
  f <- tally_updn(fm_file)
  
  summary_updown[nrow(summary_updown) + 1, ] <-
    list(cluster_name, d[1], d[2], f[1], f[2])
}

write.csv(summary_updown,
          paste0("Summary_UpDown_ExcludeNone", date_stamp, ".csv"),
          row.names = FALSE)

# --- Summary 2: DESeq2 consistent / non-consistent tally, per cluster
# Reads the consistency files written by loop 2.
summary_consistency <- data.frame(
  cluster            = character(),
  ConsistentUp       = integer(), ConsistentDown    = integer(),
  NonConsistentUp    = integer(), NonConsistentDown = integer(),
  stringsAsFactors = FALSE
)

count_rows <- function(filepath) {
  if (length(filepath) == 0 || !file.exists(filepath[1])) return(0L)
  df <- read.csv(filepath[1], row.names = 1, check.names = FALSE)
  nrow(df)
}

for (cluster_dir in cluster_dirs) {
  cluster_name <- basename(cluster_dir)
  
  cu <- list.files(cluster_dir,
                   pattern = paste0("^DESeq_ConsistentUp_Cluster_", cluster_name, ".*\\.csv$"),
                   full.names = TRUE)
  cd <- list.files(cluster_dir,
                   pattern = paste0("^DESeq_ConsistentDown_Cluster_", cluster_name, ".*\\.csv$"),
                   full.names = TRUE)
  ncu <- list.files(cluster_dir,
                    pattern = paste0("^DESeq_NonConsistentUp_Cluster_", cluster_name, ".*\\.csv$"),
                    full.names = TRUE)
  ncd <- list.files(cluster_dir,
                    pattern = paste0("^DESeq_NonConsistentDown_Cluster_", cluster_name, ".*\\.csv$"),
                    full.names = TRUE)
  
  summary_consistency[nrow(summary_consistency) + 1, ] <-
    list(cluster_name, count_rows(cu), count_rows(cd),
         count_rows(ncu), count_rows(ncd))
}

write.csv(summary_consistency,
          paste0("Summary_Consistency_DESeq", date_stamp, ".csv"),
          row.names = FALSE)

