#From Kareem

glut.markers <- FindAllMarkers(glut.subset2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
glut.markers %>% 
        group_by(cluster) %>% 
        top_n(n = 10, wt = avg_log2FC) -> top10.glut

top10.glut.df <- as.data.frame(top10.glut)

write.csv(top10.glut.df, file = "/Users/woodskad/Library/CloudStorage/Box-Box/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/S3_panels/top10.glut.csv")

glut.genes <- data.frame(table(top10.glut$gene))
colnames(glut.genes)[1] = "Gene"
df2 <- data.frame(Gene = "Slc17a7", Freq = 1)
glut.genes2 = rbind(glut.genes, df2)

pdf(file = "/Users/woodskad/Library/CloudStorage/Box-Box/mRFP-snSeq/Project1_XPoSEseq/XPoSEseq_manuscript/S3_panels/Heatmaps/glut_top10genes_heatmap052423.pdf",
    width = 3.5,
    height = 4.5)
glut.heatmap <- DoHeatmap(object = glut.subset2, 
                          features = top10.glut$gene,
                          group.colors =  c("#64C7C8", "#41B75F", "#2C8CB9", "#0A5B8C", "#3C9E64", "#6F499D"),
                          label = (FALSE),
                          raster = TRUE)+
        
        theme(axis.title = element_blank(), axis.text.y = element_text(face = "italic", size = 5), legend.position = "none")+
        
        glut.heatmap
dev.off()