#' Creates Heatmap and writes the CSV file of the top markers
#'
#' @param seur_obj Seurat object of interest
#' @param csvfilen file name for the CSV
#' @param pdffilen file name for the Heatmap PDF
#' @param groupcol a vector of colors to use for the heatmap according to cluster
#'
#' @return csv file + pdf of heatmap
#' @export
make_heatmap <- function(seur_obj, csvfilen = "topgenemarkers", pdffilen = "heatmap", groupcol){
    
    data.markers <- FindAllMarkers(seur_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    data.markers %>% 
        group_by(cluster) %>% 
        top_n(n = 10, wt = avg_log2FC) -> top10.data
    
    top10.data.df <- as.data.frame(top10.data)
    
    write.csv(top10.data.df, 
              file = paste0(csvfilen,".csv"))
    
    data.genes <- data.frame(table(top10.data$gene))
    colnames(data.genes)[1] = "Gene"
    
    pdf(file = paste0(pdffilen,".pdf"),
        width = 3.5,
        height = 4.5)
    data.heatmap <- DoHeatmap(object = seur_obj, 
                              features = top10.data$gene,
                              label = (FALSE),
                              group.colors = groupcol,
                              raster = TRUE) +
        theme(axis.title = element_blank(),
              axis.text.y = element_text(face = "italic", size = 5),
              legend.position = "none")
    
    print(data.heatmap)
    dev.off()
}
