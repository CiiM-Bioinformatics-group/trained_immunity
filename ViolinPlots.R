#tcols <- colorcode
#mono <- seurat.obj
#genes <- markers to show

Idents(mono)<- "new.id"
violin.plot.data <- get.violin.data1(mono, gene)

get.violin.data1 <- function(seurat, genes) {
  output = data.frame(gene = character(0), value= numeric(0), ident = character(0), tech=character(0))
  #####################################################################################################
  for (gene in genes) {
    if(any(gene == seurat@assays$RNA@data@Dimnames[[1]])){
      data.use = data.frame(FetchData(seurat,gene))
      data.use = t(data.use)
      data.melt=data.frame(rep(gene, length(seurat@active.ident)))
      colnames(data.melt)[1]="gene"
      data.melt$value=as.numeric(data.use[1,1:length(seurat@active.ident)])
      data.melt$id=names(data.use)[1:length(seurat@active.ident)]
      data.melt$ident=seurat@active.ident

      data.melt$value=as.numeric(as.character(data.melt$value))
      output = rbind(output, data.melt)
    } else {
      data.melt=data.frame(
        gene = rep(gene, seurat@assays$RNA@data@Dim[2]),
        value = rep(0, seurat@assays$RNA@data@Dim[2]), 
        ident = seurat@active.ident
      )
      output = rbind(output, data.melt)
    }
  }
  return(output)
}

p <- ggplot(violin.plot.data, aes(factor(ident),value)) +
  geom_violin(scale="width",adjust=1,trim=TRUE,aes(fill=factor(ident)),show.legend = F,draw_quantiles = c(0.25, 0.5, 0.75), color = "grey30") + 
  ylab("") + xlab("") + theme_bw()+
  coord_flip() +facet_wrap(~ gene,scales = "free_x", ncol = length(levels(violin.plot.data$gene))) + theme_bw()+
  theme(strip.text.x = element_text(size=12, angle=-90),
        axis.text.y = element_text(size=12),
        strip.background = element_blank(),
        panel.spacing.x = unit(c(0.4), "lines"),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=-90, vjust=0.1),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0.5,0.4,0.5,1.3), "cm")) +
  scale_x_discrete(limits = rev(levels(violin.plot.data$ident)), position = "top")+
  scale_fill_manual(values= tcols)




