

make_volcano <- function(res_tbl, object, thres = 0.05, hover = FALSE, up_color = "red",
                         down_color = "blue", export = "FALSE"){
    library(dplyr)
    library(tibble)
    
    #makes new metadata
    obj <- object
    obj$cluster_group <- paste(obj$seurat_clusters, obj$group)
    
    #
    #tibble being used 
    if ((length(class(res_tbl)) > 1 & class(res_tbl)[1] == "tbl_df") | (length(class(res_tbl)) == 1 & "data.frame" %in% class(res_tbl))){
        #tibble
        #if str_dtect > 1
        if (length(class(res_tbl)) > 1 & class(res_tbl)[1] == "tbl_df"){
          print("tibbs")
        }
        if (length(class(res_tbl)) == 1 & "data.frame" %in% class(res_tbl)){
          print("data")
        }
      
        clusters <- c()
        
        #check for multiples using str_detect tech
        if (length(names(res_tbl)[str_detect(names(res_tbl), "padj")]) > 1){
            
            #makes our option list
            for (name in names(res_tbl)[str_detect(names(res_tbl), "padj")]){
                clusters <- append(clusters, str_split_1(name, "padj_")[2])
            }
            
            clst <- ""
            
            #Pick a cluster
            while (clst == "" | !(clst %in% clusters)) {
                if (clst != ""){
                  cat("Invalid input. Please enter the name of an element",
                      "in the list.\n")
                }
                cat("Which section would you like to use. Options below.\n")
                sapply(clusters, print, quote = FALSE)
                clst = readline()
            }
            res <- res_tbl[, names(res_tbl)[str_detect(names(res_tbl), clst)]]
            #tibby
            if (!("gene" %in% colnames(res)) & ((length(class(res_tbl)) > 1 & class(res_tbl)[1] == "tbl_df") | ("gene" %in% colnames(res_tbl) & length(class(res_tbl)) == 1 & "data.frame" %in% class(res_tbl)))){
              res <- add_column(res, gene = res_tbl[, "gene"], .before = 1)
            
            #frame
            }
            # } else if (!("gene" %in% colnames(res)) & "gene" %in% colnames(res_tbl)){
            #   res <- add_column(res, gene = res_tbl[, "gene"], .before = 1)
            # }
        } else {
            res <- res_tbl
        }
        #print(res)
          #filter
        #creates new column
        res$keep <- FALSE
        
        #marks significant rows as TRUE
        for (name in names(res)[str_detect(names(res), "padj")]){
          res$keep[res[, name] < thres] <- TRUE
        }
        
        #removes unmarked rows
        res <- res[res$keep, ]
        
        #small check
        cat("There should be a TRUE under this line\n",
            !(FALSE %in% res$keep), "\n", sep = "")
        
        #removes extra col
        res <- res[, colnames(res) != "keep"] 
        print(res)
        
        res$diffE[res$log2FoldChange > 0 & res$padj < 0.05] <- "UP"
        res$diffE[res$log2FoldChange < 0 & res$padj < 0.05] <- "DOWN"
        res$delabel<- NA
        
        print(res)
      
    } else {
        #list being used
        if (class(res_tbl) == "list"){
            #list
            if (length(res_tbl) > 1){
                elem <- ""
                #tell user to pick elem in the list 
                while (elem == "" | !(elem %in% names(res_tbl))){
                    if (elem != ""){
                        cat("Invalid input. Please enter the name of an element",
                            "in the list.\n")
                    }else {
                        cat("Which element of your list do you want to use?\n",
                            "Elems include:\n", sep = "")
                    }
                    sapply(names(res_tbl), print, quote = FALSE)
                    elem = readline()
                }
                
                #runs function on chosen table
                return(make_volcano(res_tbl[[elem]], object))
              
            #single elem list   
            } else {
                #tbl <- res_tbl[[1]]
                return(make_volcano(res_tbl[[1]], object))
            }
        } else {    
        # #data.frame being used
        # } else if (class(res_tbl) == "data.frame"){
        #     #data frame
        #     print("data")
        #   
        #     clusters <- c()
        #   
        #     #check for multiples using str_detect tech
        #     if (length(names(res_tbl)[str_detect(names(res_tbl), "padj")]) > 1){
        #       
        #         #makes our option list
        #         for (name in names(res_tbl)[str_detect(names(res_tbl), "padj")]){
        #             clusters <- append(clusters, str_split_1(name, "padj_")[2])
        #         }
        #       
        #         clst <- ""
        #       
        #         #Pick a cluster
        #         while (clst == "" | !(clst %in% clusters)) {
        #             if (clst != ""){
        #               cat("Invalid input. Please enter the name of an element",
        #                   "in the list.\n")
        #             }
        #             cat("Which section would you like to use. Options below.\n")
        #             sapply(clusters, print, quote = FALSE)
        #             clst = readline()
        #         }
        #         #res <- clst
        #         res <- res_tbl[, names(res_tbl)[str_detect(names(res_tbl), clst)]]
        #         if (!("gene" %in% colnames(res)) & "gene" %in% colnames(res_tbl)){
        #             res <- add_column(res, gene = res_tbl[, "gene"], .before = 1)
        #         }
        #     } else {
        #       res <- res_tbl
        #     }
        #     #print(res)
        #     #filter
        #     #creates new column
        #     res$keep <- FALSE
        #     
        #     #marks significant rows as TRUE
        #     for (name in names(res)[str_detect(names(res), "padj")]){
        #       res$keep[res[, name] < thres] <- TRUE
        #     }
        #     
        #     #removes unmarked rows
        #     res <- res[res$keep, ]
        #     
        #     #small check
        #     cat("There should be a TRUE under this line\n",
        #         !(FALSE %in% res$keep), "\n", sep = "")
        #     
        #     #removes extra col
        #     res <- res[, colnames(res) != "keep"] 
        #     print(res)
        #     
        #     #labels for changes
        #     res$diffE[res$log2FoldChange > 0 & res$padj < 0.05] <- "UP"
        #     res$diffE[res$log2FoldChange < 0 & res$padj < 0.05] <- "DOWN"
        #     res$delabel<- NA
        #     
        #     
        # 
        # #not supported class 
        # } else {
            #error message
            stop("res_tbl, '", deparse(substitute(res_tbl)), "', is not a supported ",
                 "class. Please use an object of class 'tibble', 'list', or ", 
                 "'data.frame'")
        }
    }
}


