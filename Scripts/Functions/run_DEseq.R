#' Get Result Tibbles 2
#' 
#' This version combines prep_DEseq and get_result_tbl. Additionally, this has 
#' default parameters to make your life easier (half the time). The meta_tbl 
#' param can be left alone regardless of object unless you want to use different
#' metadata table. 
#'
#' @param data_tab Aggregate table for the given object, returned by 'get_cnt_tbls'
#' @param object SE
#' @param dataname name of object, used a lot for once
#' @param meta_tbl The table of metadata
#'
#' @return a result table
#' @export
#'
#' @examples
#' #makes result table for a cluster in the glut object
#' my_results <- run_DEseq()
#' #makes result table for a cluster in the GABA object
#' my_results <- run_DEseq(gaba_tabs, gaba, "Gaba")
#' #result table for a cluster in glut but w/ a different metadata table
#' my_results <- run_DEseq(meta_tbl = my_meta_tbl)
run_DEseq <- function(data_tab = glut_tabs, object = glut, 
                      dataname = "glut", meta_tbl = metaData){
  #Gets count table
  cnt_tbl <- prep_DEseq(data_tab, object, dataname)
  
  #gets and returns result table
  return(get_result_tbl(meta_tbl, cnt_tbl, object))
}