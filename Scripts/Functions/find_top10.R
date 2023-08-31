#' Title
#'
#' @param dds 
#' @param coexp 
#' @param cluster 
#'
#' @return
#' @export
#'
#' @examples
find_top10 <- function(dds, coexp, cluster){

        source("~/XPoSE/Scripts/Functions/sort_df.R")
        
        # extract all normalized counts for each sample
        count <- counts(dds, normalized=TRUE)
        
        # define the cluster and genes to pull, here only sig. upregulated
        FCvar <- paste0("log2FoldChange_", cluster)
        PVvar <- paste0("padj_", cluster)
        cond <- c(0, 0.05)
        AN_Ordered <- coexp[order(coexp[[FCvar]]),]
        test <- AN_Ordered %>%
                filter(
                        .data[[FCvar]] > cond[[1]],
                        .data[[PVvar]] < cond[[2]]
                )
        
        test <- sort_df(test, col_sort = paste0("padj_", cluster))
        test2 <- head(test, n = 50)
        
        
        # filter the 
        filtered_count <- as.data.frame(count[rownames(count) %in% rownames(test2),])
        filtered_count$average <- rowMeans(filtered_count[, c("NC-1:Non-active", 
                                                              "NC-2:Non-active", 
                                                              "NC-3:Non-active", 
                                                              "NC-4:Non-active")])
        
        for (i in 1:8) {
                new_column_name <- paste0("FC_", i)
                filtered_count[new_column_name] <- log2(filtered_count[i] / filtered_count$average)
        }
        
        # copy over the adj_pval from coexp table
        idx <- match(rownames(filtered_count), rownames(test2))
        filtered_count$adjpval <- test2[[PVvar]][ idx ]
        filtered_count$adjpvalT <- -log10(filtered_count$adjpval)
        filtered_count <- sort_df(filtered_count, col_sort = "adjpval")
        
        # save as csv
        write.csv(filtered_count, file = paste0("~/top10", cluster, ".csv"))
}