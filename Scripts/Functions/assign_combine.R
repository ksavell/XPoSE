#' Assigns Metadata for each experiment group, ratID and sex
#' 
#' Note: The assignment of specific sample tags to each group is manually
#' assigned here. One would need to go in and make changes to sample tags
#' used respectively.
#'
#' @param seur_list list with Seurat objects (use create_seur.R to produce this)
#' @param cart_no cartridge number to assign
#' @param ratIDs IDs of rats used
#' @param st_list list of Sample Tags
#' @param rat_sex sex of rats used
#'
#' @return list with assigned Seurat metadata for cartridge specified
#' @export
assign_combine <- function(seur_list, cart_no, ratIDs, st_list, rat_sex){
  
  temp <- seur_list[[cart_no]]
  temp$group <- rep(NA,length(temp$Sample_tag))
  temp$sex <- temp$group
  temp$ratID <- temp$group
  
  for (s in 2:length(st_list)) {
    if (s %in% c(2,4,6,8)) {
      temp$group[temp$Sample_tag==st_list[s]] <- "homecage"
    }
    else if(s %in% c(3,7)){
      if (cart_no == 1) {
        temp$group[temp$Sample_tag==st_list[s]] <- "negative"
      }
      else {
        temp$group[temp$Sample_tag==st_list[s]] <- "positive"
      }
    }
    else if (s %in% c(5,9)) {
      if (cart_no == 1) {
        temp$group[temp$Sample_tag==st_list[s]] <- "positive"
      }
      else {
        temp$group[temp$Sample_tag==st_list[s]] <- "negative"
      }
    }
    
    temp$sex[temp$Sample_tag==st_list[s]] <- rat_sex[s]
    temp$ratID[temp$Sample_tag==st_list[s]] <- ratIDs[s]
  }
  seur_list[[cart_no]] <- temp
  return(seur_list)
}