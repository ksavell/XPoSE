# Run single-factor pseudobulk DE across cell populations and summarize DEGs
single_factor_DESeq <- function(object, comp_vect, cluster, min_cell = 10,
                                min_rat = 3, keep_dds = FALSE) {
  library(Seurat)
  library(Libra)
  library(dplyr)
  library(tibble)
  library(DESeq2)
  library(rlang)
  library(rlist)
  
  # Check that cluster exists
  if (!(cluster %in% unique(object$cluster_name))) {
    stop(
      "Cluster '", cluster, "' not found within object."
    )
  }
  
  # Require at least 2 biological replicates per population.
  # Default is 3 rats/population for primary DE analyses, but min_rat = 2 can be
  # explicitly supplied for leave-one-out sensitivity analyses.
  if (min_rat < 2) {
    min_rat <- 2
  }
  
  # Require at least 10 nuclei per population per rat
  if (min_cell < 10) {
    min_cell <- 10
  }
  
  # Check comparison factor
  if (!(comp_vect[1] %in% colnames(object@meta.data))) {
    stop(
      "Factor ", comp_vect[1],
      " not found within object metadata."
    )
  }
  
  # Check comparison levels
  if (FALSE %in% (comp_vect[2:3] %in% object[[comp_vect[1]]][, 1])) {
    stop(
      "Comparison contains levels not found in factor ",
      comp_vect[1], "."
    )
  }
  
  # Subset to cluster and comparison populations
  Idents(object) <- "cluster_name"
  sub_obj <- subset(object, idents = cluster)
  
  Idents(sub_obj) <- comp_vect[1]
  sub_obj <- subset(sub_obj, idents = comp_vect[2:3])
  
  # Build table of nuclei counts per rat and comparison population
  t_tbl <- data.frame(matrix(nrow = 1, ncol = 4))
  
  rat_g_c <- table(
    sub_obj[[comp_vect[1]]][, 1],
    sub_obj$cluster_name,
    sub_obj$ratID
  )
  
  colnames(t_tbl) <- c(
    comp_vect[1],
    "ratID",
    "counts",
    "Included"
  )
  
  curr_row <- 1
  
  for (rat in unique(sub_obj$ratID)) {
    for (level in comp_vect[2:3]) {
      
      if (rat_g_c[level, , rat] != 0) {
        
        t_tbl[curr_row, ] <- list(
          level,
          rat,
          rat_g_c[level, , rat],
          rat_g_c[level, , rat] >= min_cell
        )
        
        rownames(t_tbl)[curr_row] <- paste(
          rat,
          level,
          sep = ":"
        )
        
        curr_row <- curr_row + 1
      }
    }
  }
  
  # Summarize included biological replicates
  incl_tbl <- data.frame(matrix(0, nrow = 1, ncol = 2))
  
  rownames(incl_tbl) <- comp_vect[1]
  colnames(incl_tbl) <- comp_vect[2:3]
  
  rat_ttls <- c()
  
  for (i in seq_len(nrow(t_tbl))) {
    
    if (t_tbl[i, "Included"]) {
      incl_tbl[, t_tbl[i, comp_vect[1]]] <-
        incl_tbl[, t_tbl[i, comp_vect[1]]] + 1
    }
    
    if (!(t_tbl[i, comp_vect[1]] %in% names(rat_ttls))) {
      rat_ttls[t_tbl[i, comp_vect[1]]] <- 1
    } else {
      rat_ttls[t_tbl[i, comp_vect[1]]] <-
        rat_ttls[t_tbl[i, comp_vect[1]]] + 1
    }
  }
  
  # Check sample availability
  low_samp_flag <- FALSE
  
  for (level in comp_vect[2:3]) {
    
    if (incl_tbl[, level] < min_rat) {
      low_samp_flag <- TRUE
    }
    
    if (incl_tbl[, level] == 0) {
      stop(
        "Level ", level,
        " within factor ", comp_vect[1],
        " was fully excluded with min_cell = ", min_cell, "."
      )
    }
  }
  
  # Add included / total counts for display
  for (i in seq_along(rat_ttls)) {
    incl_tbl[names(rat_ttls)[i]] <- paste(
      incl_tbl[names(rat_ttls)[i]],
      "/",
      rat_ttls[i],
      sep = ""
    )
  }
  
  cat("Sample count: Include/Total\n")
  print(incl_tbl)
  
  if (low_samp_flag) {
    warning(
      "One or more populations contain fewer than ",
      min_rat,
      " included rats.\n",
      immediate. = TRUE
    )
  }
  
  # Format individual sample nuclei-count table
  prin_tbl <- as.data.frame(
    t_tbl %>%
      arrange(
        grepl(
          comp_vect[2],
          !!as.symbol(comp_vect[1])
        )
      )
  )
  
  for (i in seq_len(nrow(prin_tbl))) {
    rownames(prin_tbl)[i] <- paste(
      prin_tbl[i, "ratID"],
      prin_tbl[i, comp_vect[1]],
      sep = ":"
    )
  }
  
  cat("\n")
  cat("Individual Sample Nuclei Counts\n")
  print(prin_tbl)
  cat("\n")
  
  # Identify samples excluded for insufficient nuclei
  exclusion <- rownames(
    t_tbl[which(!t_tbl[, "Included"]), ]
  )
  
  if (length(exclusion) != 0) {
    
    if (length(exclusion) > 1) {
      cat(
        paste(exclusion, collapse = ", "),
        " excluded because counts < ",
        min_cell,
        "\n\n",
        sep = ""
      )
    } else {
      cat(
        exclusion,
        " excluded because counts < ",
        min_cell,
        "\n\n",
        sep = ""
      )
    }
    
  } else {
    cat("No rats excluded\n\n")
  }
  
  # Retain samples meeting nuclei-count requirement
  temp <- t_tbl[
    which(t_tbl[, "Included"]),
    ,
    drop = FALSE
  ]
  
  # Check that both comparison populations remain
  for (cmp in comp_vect[2:3]) {
    
    if (!(cmp %in% temp[, comp_vect[1]])) {
      stop(
        "Nuclei-count filtering excluded all rats from ",
        cmp,
        " for factor ",
        comp_vect[1],
        "."
      )
    }
  }
  
  # Require requested number of biological replicates after filtering
  grp_counts <- table(temp[, comp_vect[1]])
  
  if (any(grp_counts < min_rat)) {
    
    low <- names(grp_counts)[grp_counts < min_rat]
    
    stop(
      "SKIP: population(s) ",
      paste(low, collapse = ", "),
      " have < ",
      min_rat,
      " rats after filtering (",
      paste(
        paste0(
          names(grp_counts),
          "=",
          as.integer(grp_counts)
        ),
        collapse = ", "
      ),
      ")."
    )
  }
  
  # Remove excluded rats
  Idents(sub_obj) <- "ratID"
  sub_obj <- subset(
    sub_obj,
    idents = unique(temp$ratID)
  )
  
  # Pseudobulk counts by rat
  counts <- to_pseudobulk(
    sub_obj,
    replicate_col = "ratID",
    cell_type_col = "cluster_name",
    label_col = comp_vect[1]
  )[[1]][, rownames(temp), drop = FALSE]
  
  # DESeq2 design
  form_fact <- as.formula(
    paste("~", comp_vect[1])
  )
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = temp,
    design = form_fact
  )
  
  dds <- DESeq(dds)
  
  # Differential expression
  deseq_results <- results(
    dds,
    contrast = comp_vect,
    alpha = 0.05,
    cooksCutoff = FALSE,
    independentFiltering = FALSE
  )
  
  deseq_results <- deseq_results %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble()
  
  # Return DESeq object only when requested
  if (keep_dds) {
    return(
      list(
        dds = dds,
        results = deseq_results
      )
    )
  }
  
  return(
    list(
      results = deseq_results
    )
  )
}

de_and_summary <- function(seur_obj,
                           pair,
                           clusters = sort(unique(seur_obj$cluster_name)),
                           min_cell = 10,
                           min_rat = 3,
                           output_dir = ".",
                           save_dds = TRUE) {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  comparison_tag <- paste(pair, collapse = "_")
  results_list <- list()
  result_i <- 1
  
  for (cl in clusters) {
    
    deseq2_results <- tryCatch(
      single_factor_DESeq(
        object = seur_obj,
        comp_vect = pair,
        cluster = cl,
        min_cell = min_cell,
        min_rat = min_rat
      ),
      error = function(e) {
        message(
          "Skipping ", cl, " for ", comparison_tag, ": ",
          conditionMessage(e)
        )
        NULL
      }
    )
    
    if (is.null(deseq2_results)) {
      next
    }
    
    de_tbl <- deseq2_results$results
    dds <- deseq2_results$dds
    
    if (!all(c("padj", "log2FoldChange") %in% colnames(de_tbl))) {
      warning("padj/log2FoldChange missing for ", cl, "; skipping.")
      next
    }
    
    de_tbl$padj <- as.numeric(de_tbl$padj)
    de_tbl$log2FoldChange <- as.numeric(de_tbl$log2FoldChange)
    
    n_up <- sum(
      !is.na(de_tbl$padj) &
        de_tbl$padj < 0.05 &
        de_tbl$log2FoldChange > 0
    )
    
    n_down <- sum(
      !is.na(de_tbl$padj) &
        de_tbl$padj < 0.05 &
        de_tbl$log2FoldChange < 0
    )
    
    # Preserve the original F5A convention: upregulated DEG counts are
    # positive and downregulated DEG counts are negative.
    results_list[[result_i]] <- data.frame(
      Category = comparison_tag,
      Observation = cl,
      Value = n_up,
      stringsAsFactors = FALSE
    )
    
    results_list[[result_i + 1]] <- data.frame(
      Category = comparison_tag,
      Observation = cl,
      Value = -n_down,
      stringsAsFactors = FALSE
    )
    
    result_i <- result_i + 2
    
    file_tag <- paste0(cl, "_", comparison_tag)
    
    write.csv(
      de_tbl,
      file.path(output_dir, paste0(file_tag, ".csv")),
      row.names = FALSE
    )
    
    if (save_dds) {
      saveRDS(
        dds,
        file.path(output_dir, paste0(file_tag, "_dds.rds"))
      )
    }
  }
  
  if (length(results_list) == 0) {
    warning("No successful DE comparisons for ", comparison_tag)
    return(invisible(NULL))
  }
  
  summary_df <- do.call(rbind, results_list)
  
  write.csv(
    summary_df,
    file.path(output_dir, paste0(comparison_tag, "_summary.csv")),
    row.names = FALSE
  )
  
  invisible(summary_df)
}
