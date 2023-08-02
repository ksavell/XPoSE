#' Prep Upset
#'
#' @param coexp table of results data
#' @param p_thres Value used to filter for significance
#'
#' @return list of combination matrices
#' @export
#'
#' @examples
prep_upset <- function(coexp, p_thres = 0.05){  
    #gets names of all fold and p-val cols
    p_cols <- names(coexp)[str_detect(names(coexp), "padj")]
    fold_cols <- names(coexp)[str_detect(names(coexp), "log2FoldChange")]
    
    #makes our base_df
    base_df <- data.frame(matrix(nrow = nrow(coexp), ncol = length(p_cols)))
    rownames(base_df) <- rownames(coexp)
    
    #returned at functions end, holds matrices
    mat_lst <- list()
    
    #makes df for up and down
    for (dir in c("up", "dwn")){
        
        #creates df from base
        df <- base_df
        for (i in 1:length(p_cols)){
            #gets name of current cluster
            clust <- str_split_1(p_cols[i], "_")[2]
            
            #fills significant cols with 1's
            if (dir == "up"){#sect <- paste()
                names(df)[i] <- paste(dir, "_", clust, sep = "")
                
                df[, paste(dir, "_", clust, sep = "")][
                  ((coexp[[p_cols[i]]] < p_thres) & (coexp[[fold_cols[i]]] > 0))
              ] <- 1
              
            } else {
                names(df)[i] <- paste(dir, "_", clust, sep = "")
                
                df[[paste(dir, "_", clust, sep = "")]][
                  ((coexp[[p_cols[i]]] < p_thres) & (coexp[[fold_cols[i]]] < 0))
                ] <- 1
            }
        }
        #fills evrything else w/ 0's
        df[is.na(df)] <- 0
        
        #makes matrices from frames
        distMat <- make_comb_mat(df, mode = c("distinct"))
        interMat <- make_comb_mat(df, mode = c("intersect"))
        
        #adds matrices to list and appends
        mat_lst <- list.append(mat_lst, dist = distMat, inter = interMat)
        names(mat_lst)[match("dist", names(mat_lst))] <- paste("distinct",
                                                               dir, sep = "_")
        names(mat_lst)[match("inter", names(mat_lst))] <- paste("intersect", 
                                                                dir, sep = "_")
    }
    
    #returns list
    return(mat_lst)
}
