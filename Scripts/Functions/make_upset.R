#' Make Upset
#'
#' @param matrc_lst list of combination matrices returned by prep_upset
#' @param order_vect order you want the upset to list clusters
#' @param export logical, marks whether you want to save as pdf (TRUE) or to 
#'               environment (FALSE). TRUE by default.
#' @param dest target destination for files, default is current directory. 
#'
#' @return nothing or a list of upsets
#' @export
#'
#' @examples
make_upset <- function(matrc_lst, order_vect, export = TRUE, dest = getwd()){
    #will return this list at end if export false
    if (!export){
        retrn_lst <- list()
    }
    
    clusts <- order_vect
    
    #loops through the list
    for (i in 1:length(matrc_lst)){
        
        #ensures clusters in order_vect are tied to matrix counterpart
        for (j in 1:length(clusts)){
            order_vect[j] <- attr(matrc_lst[[i]],
                                "dimnames")[[1]][str_detect(attr(matrc_lst[[i]],
                                                            "dimnames")[[1]],
                                                             clusts[j])]
        }
        
        #makes the upset!
        new_upS <- UpSet(matrc_lst[[i]],
                      top_annotation = upset_top_annotation(matrc_lst[[i]], 
                                                            add_numbers = TRUE),
                      right_annotation = upset_right_annotation(matrc_lst[[i]],
                                                            add_numbers = TRUE),
                      set_order = order_vect)
        
        #exports to file
        if (export){
            pdf(file = paste(dest, "/", names(matrc_lst)[i], ".pdf", sep = ""),
                width = 4,
                height = 4)
            UpSet(t(matrc_lst[[i]]), set_order = order_vect)
            print(new_upS)
            dev.off()
            
        #sends to list
        } else {
            retrn_lst <- list.append(retrn_lst, new_upS)
            names(retrn_lst)[i] <- names(matrc_lst)[i]
        }
    }
    
    #returns 
    if (!export){
       return(retrn_lst)
    }
}
