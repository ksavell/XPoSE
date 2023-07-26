#' Explore markers during clustering
#'
#' @param seur_obj the seurat object you are clustering
#' @param qc plots nFeatures and nCount
#' @param all plots markers for neurons and glia contamination
#' @param glut plots glut neuron markers
#' @param gaba plots gaba neuron markers
#'
#' @return vln plot of qc plots or markers for all, glut, or gaba neurons
#' @export
#'
#' @examples
clstr_vln <- function(seur_obj, qc = FALSE, all = FALSE, 
                             glut = FALSE, gaba = FALSE) {
        plot_list <- list()
        
        if (qc == TRUE) {
                plot_list[["qc"]] <- VlnPlot(seur_obj, 
                                             features = c("nFeature_RNA", "nCount_RNA"))
        }
        if (all == TRUE) {
                plot_list[["all"]] <- VlnPlot(seur_obj, 
                                              features = c("Slc17a7", "Gad1",
                                                           "Snap25", "Mbp",
                                                           "Gja1", "Col5a3"))
        }
        if (glut == TRUE) {
                plot_list[["glut"]] <- VlnPlot(seur_obj, 
                                               features = c("Slc17a7", "Gpc5",
                                                            "Rfx3", "Rorb",
                                                            "Syt6", "Ctgf",
                                                            "Tshz2"))
        }
        if (gaba == TRUE) {
                plot_list[["gaba"]] <- VlnPlot(seur_obj, 
                                               features = c("Gad1","Kcnc2",
                                                            "Sst","Chodl",
                                                            "Ppp1r1b", "Meis2",
                                                            "Vip", "Lamp5", 
                                                            "Ndnf","Cck"))
        }
        
        print(plot_list)
}
