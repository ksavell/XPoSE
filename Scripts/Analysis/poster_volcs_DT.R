I56_m <- merge_results(res[1])

I56_m$diffE[I56_m$log2FoldChange > 0 & I56_m$padj < 0.05] <- "UP"
I56_m$diffE[I56_m$log2FoldChange < 0 & I56_m$padj < 0.05] <- "DOWN"
I56_m$delabel<- NA

volc_5 <- ggplot(data=(I56_m),
       aes(x=log2FoldChange,
           y=-log10(padj),
           col=diffE,
           label=delabel))+
  geom_point(size= 3.0 )+
  theme_classic()+
  geom_text()+
  scale_color_manual(values=c('7A28CB','blue', 'transparent'))+
  theme(text=element_text(size=16))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0, size = 1.5)+
  xlab("Log2 (fold change)")+
  ylab("Adjusted Pvalue")
volc_5 + xlim(-5,5)

I23_m <- merge_results(res[2])
I23_m$diffE[I23_m$log2FoldChange > 0 & I23_m$padj < 0.05] <- "UP"
I23_m$diffE[I23_m$log2FoldChange < 0 & I23_m$padj < 0.05] <- "DOWN"
I23_m$delabel<- NA

volc_2 <- ggplot(data=(I23_m),
       aes(x=log2FoldChange,
           y=-log10(padj),
           col=diffE,
           label=delabel))+
  geom_point(size= 3.0 )+
  theme_classic()+
  geom_text()+
  scale_color_manual(values=c('7A28CB','FC6666', 'transparent'))+
  theme(text=element_text(size=16))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_vline(xintercept = 0, size = 1.5)+
  xlab("Log2 (fold change)")+
  ylab("Adjusted Pvalue")
volc_2 + xlim(-5,5)
