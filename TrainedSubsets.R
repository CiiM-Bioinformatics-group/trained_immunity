### assign trained subsets based on the 6 TI markers
library(Seurat)
library(dplyr)
library(Matrix)
library(Matrix.utils)
library(ggplot2)
library(reshape2)
library(pryr)
library(RColorBrewer)
library(tidyr)
library(cowplot)
library(pheatmap)

mono.pbmc.t2<- readRDS("mono.pbmc.t2.rds")
features <- c("TNF","IL1B","IL6","CXCL9","CXCL10","CXCL11")
Idents(mono.pbmc.t2) <- "celltype.short"
DefaultAssay(mono.pbmc.t2) <- "RNA"

hla.s <- mono.pbmc.t2[features,]
exp.hla <- as.data.frame(hla.s@assays$RNA@data)
anno.hla <- as.data.frame(cbind(as.character(hla.s$condition),as.character(hla.s$celltype),as.character(hla.s$time)))
row.names(anno.hla) <- colnames(hla.s)
names(anno.hla) <- c("condition","celltype","time")
Idents(hla.s) <- "time.stim"
rpmi.hla <- subset(hla.s, idents="T2.S1")
mono.mm <- rowMeans(rpmi.hla)



t2 <- (exp.hla - mono.mm)
t2t <- t(t2)
t2t<-as.data.frame(t2t)
get<-t2t[t2t$IL1B>0 | t2t$TNF>0 | t2t$IL6>0 | t2t$CXCL9 >0| t2t$CXCL10 >0| t2t$CXCL11 >0 ,]

dim(get)
RPMI.index <- grep(pattern = "s1", x = rownames(get), value = FALSE)
get <- get[-RPMI.index, ]
get.t <- t(get)

anno <- as.data.frame(cbind(as.character(hla.s$condition),as.character(hla.s$celltype),as.character(hla.s$cells)))
row.names(anno) <- colnames(hla.s)
names(anno) <- c("Stimuli","celltype","Tissue")
anno$Cell.IDs = "other"
anno$Cell.IDs[anno$celltype.sort %in% "Macrophages-1"] = "Mac1"
anno$Cell.IDs[anno$celltype.sort %in% "Macrophages-2"] = "Mac2"
anno <- anno[anno$Stimuli != "RPMI",]
anno <- anno[,-2]
anno$Cell.IDs <- as.factor(anno$Cell.IDs)
ann_colors = list(
Cell.IDs = c(`M1 Mpg`="#F8766D", `M2 Mpg`="#AA9E00", other="#09BC5C"),
    Stimuli = c(BG= "#4DAF4A", MDP= "#377EB8", UA="#E41A1C", oxLDL="#FFFF33"),
    Tissue = c(MONO="black",PBMC="lightgrey")
)

hm <- pheatmap(get.t, annotation_col = anno,breaks=seq(-5,5,0.1),show_colnames=F, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),annotation_colors = ann_colors)


plot(hm$tree_col)
rect.hclust(hm$tree_col, k=4)
cut <- cutree(hm$tree_col, k=4)

t1 <- names(cut[cut==1])
t2 <- names(cut[cut==2])
t3 <- names(cut[cut==3])
t4 <- names(cut[cut==4])


mono.pbmc.t2$newTI <- "NT"
mono.pbmc.t2$newTI[mono.pbmc.t2$stim %in% "RPMI.T2"] <- "RPMI"
mono.pbmc.t2$newTI[mono.pbmc.t2$cellnames %in% t4] <- "MC"
mono.pbmc.t2$newTI[mono.pbmc.t2$cellnames %in% t1] <- "MCI"
mono.pbmc.t2$newTI[mono.pbmc.t2$cellnames %in% t3] <- "MCI"

mono.pbmc.t2$newTI4 <- "NT"
mono.pbmc.t2$newTI4[mono.pbmc.t2$stim %in% "RPMI.T2"] <- "RPMI"
mono.pbmc.t2$newTI4[mono.pbmc.t2$cellnames %in% t4] <- "MTC"
mono.pbmc.t2$newTI4[mono.pbmc.t2$cellnames %in% t1] <- "MPC.weak"
mono.pbmc.t2$newTI4[mono.pbmc.t2$cellnames %in% t3] <- "MPC.strong"

NTvRPMI <- FindMarkers(mono.pbmc.t2, ident.1 = "NT", ident.2 = "RPMI", group.by="newTI")
MTCvRPMI <- FindMarkers(mono.pbmc.t2, ident.1 = "MC", ident.2 = "RPMI", group.by="newTI")
MPCvRPMI <- FindMarkers(mono.pbmc.t2, ident.1 = "MCI", ident.2 = "RPMI", group.by="newTI")

mks <- FindAllMarkers(object = mono.pbmc.t2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- mks %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

DoHeatmap(object = mono.pbmc.t2, features = top10$gene, size=3) + NoLegend()

mono.pbmc.t2$newTI.cond <- paste(mono.pbmc$newTI, mono.pbmc$condition, sep=".")

stimulus <- c("BG","UA","oxLDL","MDP")

for(index1 in stimulus){
 index1
 a <- paste("MCI",index1,sep=".")
 b <- paste("MC",index1,sep=".")
 c <- paste("NT",index1,sep=".")

 a.markers <- FindMarkers(object = mono.pbmc.t2, ident.1 = a, ident.2="RPMI.RPMI", group.by = "newTI.cond", verbose = FALSE, logfc.threshold = 0.1)
 b.markers <- FindMarkers(object = mono.pbmc.t2, ident.1 = b, ident.2="RPMI.RPMI", group.by = "newTI.cond", verbose = FALSE, logfc.threshold = 0.1)
 c.markers <- FindMarkers(object = mono.pbmc.t2, ident.1 = c, ident.2="RPMI.RPMI", group.by = "newTI.cond", verbose = FALSE, logfc.threshold = 0.1)
 
 a.markers <- a.markers[a.markers$p_val<0.05,]
 b.markers <- b.markers[b.markers$p_val<0.05,]
 c.markers <- c.markers[c.markers$p_val<0.05,]
 
 name1 <- paste("./gMCI.",index1,".csv",sep="")
 name2 <- paste("./gMC.",index1,".csv",sep="")
 name3 <- paste("./gNT.",index1,".csv",sep="")
 write.csv(a.markers,file=name1, quote=F)
 write.csv(b.markers,file=name2, quote=F)
 write.csv(c.markers,file=name3, quote=F)
}




