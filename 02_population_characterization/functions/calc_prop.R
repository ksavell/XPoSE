calc_prop <- function(seur_obj,
                      fact1 = "ratID",
                      fact2 = "celltype",
                      fact3 = NULL) {
  meta <- seur_obj@meta.data
  
  # ---- Input checks ----
  vars <- c(fact1, fact2, if (!is.null(fact3)) fact3)
  missing_vars <- setdiff(vars, colnames(meta))
  if (length(missing_vars) > 0) {
    stop("These metadata columns are missing: ",
         paste(missing_vars, collapse = ", "))
  }
  
  # ---- Build count table ----
  if (is.null(fact3)) {
    count_df <- as.data.frame(
      table(meta[[fact1]], meta[[fact2]]),
      stringsAsFactors = FALSE
    )
    colnames(count_df) <- c(fact1, fact2, "count")
  } else {
    count_df <- as.data.frame(
      table(meta[[fact1]], meta[[fact2]], meta[[fact3]]),
      stringsAsFactors = FALSE
    )
    colnames(count_df) <- c(fact1, fact2, fact3, "count")
  }
  
  count_df$count <- as.numeric(count_df$count)
  
  # ---- Define grouping vars for proportions ----
  # Same as your original:
  # - if no fact3: proportion within fact1 (e.g. ratID)
  # - if fact3:    proportion within (fact1, fact3) combo
  if (is.null(fact3)) {
    group_vars <- fact1
  } else {
    group_vars <- c(fact1, fact3)
  }
  
  # ---- Compute proportions within each group ----
  # 1) Build a grouping factor from group_vars
  group_fac <- interaction(count_df[, group_vars, drop = FALSE], drop = TRUE)
  
  # 2) Get per-group totals and divide
  totals <- ave(count_df$count, group_fac, FUN = sum)
  count_df$percent <- count_df$count / totals
  
  rownames(count_df) <- NULL
  count_df
}
