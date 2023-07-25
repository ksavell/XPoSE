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