make_coexp <- function(seur_obj1, seur_obj2, threshold, factor, comp_vect, 
                       filter, p_thres = 0.05, incl_all = FALSE, 
                       genes_sorted = FALSE){
    #Necessary scripts for function
    source("~/XPoSE/Scripts/Functions/run_pseudobulk.R")
    source("~/XPoSE/Scripts/Functions/execute_DESeq.R")
    source("~/XPoSE/Scripts/Functions/prep_merge_fast.R")
    source("~/XPoSE/Scripts/Functions/prep_merge.R")
    source("~/XPoSE/Scripts/Functions/merge_results.R")
  
    #for process time checking
    timer <- proc.time()
    
    #a series of error messages
    #code last
    
    #list of objects for abstraction purposes
    obj_list <- list(obj1 = seur_obj1, obj2 = seur_obj2)
    names(obj_list) <- c(deparse(substitute(seur_obj1)),
                         deparse(substitute(seur_obj2)))
    
    #In the case that the gene lists are not the same, this ensures the merge
    #doesn't cause any errors
    gene_vect <- c()
    gene_lngth_s <- TRUE
    
    #checks for different lengths or different genes within
    if ((length(seur_obj1@assays[["RNA"]]@counts@Dimnames[[1]]) !=
         length(seur_obj2@assays[["RNA"]]@counts@Dimnames[[1]])) |
        (FALSE %in% (seur_obj1@assays[["RNA"]]@counts@Dimnames[[1]] ==
                     seur_obj2@assays[["RNA"]]@counts@Dimnames[[1]]))){
       
        #gene vectors are not same length
        gene_lngth_s <- FALSE
        
        #make the unique of c(obj1, obj2) and sorts
        gene_vect <- sort(unique(c(seur_obj1@assays[["RNA"]]@counts@Dimnames[[1]],
                                   seur_obj2@assays[["RNA"]]@counts@Dimnames[[1]])),
                          method = "quick")
    }
    
    #list of merged objects to be merged into co-expression table
    merg_lst <- list()
    
    ind <- 1
    #pseudobulks, runs DESeq and merges tables
    for (obj in obj_list){
        #differing lengths or genes within
        if (gene_lngth_s){
            gene_vect <- obj@assays[["RNA"]]@counts@Dimnames[[1]]
        }
        
        #gets results tables
        obj_rt <- execute_DESeq(run_pseudobulk(obj, threshold, factor, comp_vect, incl_all),
                                obj, factor, comp_vect)
        
        #checks if data sorted if genes_sorted is TRUE
        if (genes_sorted &
            !(FALSE %in% (obj@assays[["RNA"]]@counts@Dimnames[[1]] ==
                          sort(obj@assays[["RNA"]]@counts@Dimnames[[1]],
                               method = "quick")))){
            #runs version w/ binary search and merges
            obj_rt <- prep_merge_fast(obj_rt, obj, gene_vect)
          
          #List unsorted or user did not specify
        }else {
            #run unsorted version and merges
            obj_rt <- prep_merge(obj_rt, obj, gene_vect)
        }
        
        #append to merg_lst and name
        merg_lst <- list.append(merg_lst, obj_rt)
        names(merg_lst)[ind] <- paste("thing", ind, sep = "")
        
        #incs index
        ind <- ind + 1
    }
    
    #combines lists
    merg_lst <- append(merg_lst[[1]], merg_lst[[2]])
    
    #tells user how long function took
    elpd <- (proc.time() - timer)[[3]]
    print(paste("Process took", round(elpd, digits = 3), "seconds"),
          quote = FALSE)
    
    #merge the objs
    return(merge_results(merg_lst, thres = p_thres, filter = filter))
}
