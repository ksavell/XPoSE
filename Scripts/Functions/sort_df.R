#' Sort dataframe by rows or a column 
#'
#' @param df dataframe to be sorted
#' @param col_sort optional, the column to sort by
#'
#' @return a sorted dataframe by row or column 
#' @export
#'
#' @examples
sort_df <- function(df, col_sort = NULL) {
        if (is.null(col_sort)) {
                # Sort based on row names
                sorted_df <- df[order(row.names(df)), ]
        } else {
                if (col_sort %in% colnames(df)) {
                        # Sort based on the specified column
                        sorted_df <- df[order(df[, col_sort]), ]
                } else {
                        stop("Column name not found!")
                }
        }
        
        return(sorted_df)
}