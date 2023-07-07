# Plots of how many genes detected etc ----

#' Plots how many genes detected 
#'
#' @param feature "genes" or "transcripts" depending on which feature user wants to input
#' @param cart_no cartridge number = 1 or 2
#' @param seur_objs a list of seurat objects to plot
#'
#' @return a violin plot of genes detected by feature of choice
#' @export
#'
#' @examples
#' plot_detection("genes", 1, seur_list)
plot_detection <- function(feature, cart_no, seur_objs = seur_list){
  
  if (feature == "genes"){
    VlnPlot(seur_objs[[cart_no]], features = "nFeature_RNA", group.by = "Sample_tag") +
      ggtitle(paste("Cart",cart_no," genes", sep=""))
  }
  else {
    VlnPlot(seur_list[[cart_no]], features = "nCount_RNA", group.by = "Sample_tag") +
      ggtitle(paste("Cart",cart_no," transcripts", sep=""))
  }
}