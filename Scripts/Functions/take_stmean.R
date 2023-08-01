#' Title
#'
#' @param stag_df dataframe with sample tag reads
#'
#' @return means of each sample tag reads for each sample tag call
#' @export
#'
#' @examples
take_stmean <- sapply(IDslist, function(stag_df) {
        IDs <- data.frame(combined@assays[["RNA"]]@data@Dimnames[[2]])
        IDs$cart <- combined@meta.data[["orig.ident"]]
        IDs$rat <- combined@meta.data[["ratID"]]
        IDs$st <- combined@meta.data[["Sample_tag"]]
        IDs$st2 <- combined@meta.data[["ST2_reads"]]
        IDs$st3 <- combined@meta.data[["ST3_reads"]]
        IDs$st4 <- combined@meta.data[["ST4_reads"]]
        IDs$st5 <- combined@meta.data[["ST5_reads"]]
        IDs$st6 <- combined@meta.data[["ST6_reads"]]
        IDs$st7 <- combined@meta.data[["ST7_reads"]]
        IDs$st8 <- combined@meta.data[["ST8_reads"]]
        IDs$st9 <- combined@meta.data[["ST9_reads"]]
        rownames(IDs) <- IDs$combined.assays...RNA....data.Dimnames..2..
        IDslist <- split(IDs, IDs$st)
        
        # Find the names of the numeric columns
        numeric_cols <- sapply(x, is.numeric)
        
        # Ignore character columns
        numeric_cols <- numeric_cols & !sapply(x, is.character)
        
        # Find the means of the numeric columns
        means <- apply(x[, numeric_cols], 2, mean)
        
        # Return a named vector of medians
        names(means) <- names(x)[numeric_cols]
        return(means)
})