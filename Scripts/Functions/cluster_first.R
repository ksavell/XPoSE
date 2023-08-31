#' Cluster initial object
#'
#' @param seur_obj seurat object, unclustered
#' @param pca_dim number of PCs to calculate
#' @param neigh_dim number of PCs to use in FindNeighbor
#' @param umap_dim number of PCs to use in UMAP
#' @param res resolution for FindClusters
#'
#' @return
#' @export
#'
#' @examples
cluster_first <- function(seur_obj, pca_dim = 1:50, neigh_dim = 1:50, umap_dim = 1:50, res = 2) {
       
        mitogenes <- c("ATP6","COX1","COX2","CYTB","ND1","ND2","ND4","ND5","ND6")
        seur_obj$percMito <- PercentageFeatureSet(seur_obj,
                                                  features = mitogenes)
        seur_obj <- subset(seur_obj, subset = percMito < 10)
        seur_obj$percMito <- NULL
        
        seur_obj <- seur_obj %>%
                NormalizeData() %>%
                FindVariableFeatures(selection.method = "vst", 
                                     nfeatures = 2000) %>%
                ScaleData(verbose = FALSE) %>%
                RunPCA(ndims.print = 1:50, 
                       nfeatures.print = 50) %>%
                FindNeighbors(dims = neigh_dim) %>%
                RunUMAP(reduction = "pca", 
                        dims = umap_dim) %>%
                FindClusters(resolution = res)
}