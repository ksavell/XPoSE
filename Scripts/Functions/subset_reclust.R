#' Subset and recluster
#'
#' @param seur_obj seurat object containing clusters to subset
#' @param clust_tokeep idents to keep in subset
#' @param pca_dim number of PCs to calculate
#' @param neigh_dim number of PCs to use in FindNeighbor
#' @param umap_dim number of PCs to use in UMAP
#' @param res resolution for FindClusters
#'
#' @return
#' @export
#'
#' @examples
subset_reclust <- function(seur_obj, clust_tokeep, pca_dim = 1:50, 
                           neigh_dim = 1:50, umap_dim = 1:50, res = 2) {
        seur_obj <- seur_obj %>%
                subset(idents = c(clust_tokeep)) %>%
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