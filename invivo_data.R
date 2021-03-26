
setwd("/path/to/invivo/data")
options(future.globals.maxSize = 10737418240)  # 10240*1024^2 = 10737418240  for 10GB

library(Seurat)
library(dplyr)
library(Matrix)
library(Matrix.utils)
library(ggplot2)
library(reshape2)
library(pryr)
library(RColorBrewer)

library(cowplot)

ReadBowen <- function(bdir,bid,btype,bmeta=NULL){
	bdata <- Read10X(data.dir = bdir)
	bsample <- CreateSeuratObject(counts = bdata, min.cells = 3, min.features = 100, project = "10X",meta.data=bmeta)
	#bsample$individual <- individual
	bsample$sampid <- bid
	bsample$samptype <- btype
	return(bsample)
}



pool1.info <- read.table(file="./de_conv_genotype_based/compare/pool1_convert.txt", header=T, sep="\t", stringsAsFactors=F,row.names=1)
pool2.info <- read.table(file="./de_conv_genotype_based/compare/pool2_convert.txt", header=T, sep="\t", stringsAsFactors=F,row.names=1)
pool3.info <- read.table(file="./de_conv_genotype_based/compare/pool3_convert.txt", header=T, sep="\t", stringsAsFactors=F,row.names=1)
pool4.info <- read.table(file="./de_conv_genotype_based/compare/pool4_convert.txt", header=T, sep="\t", stringsAsFactors=F,row.names=1)
pool5.info <- read.table(file="./de_conv_genotype_based/compare/pool5_convert.txt", header=T, sep="\t", stringsAsFactors=F,row.names=1)
pool6.info <- read.table(file="./de_conv_genotype_based/compare/pool6_convert.txt", header=T, sep="\t", stringsAsFactors=F,row.names=1)

######### pool 1 ##########
pbmc.pool <- ReadBowen(bdir="./pool_1",bid="pool1",btype="pool1",bmeta=pool1.info)
Idents(pbmc.pool)<-"Sample_assign"
pbmc.pool[["percent.mt"]] <- PercentageFeatureSet(object = pbmc.pool, pattern = "^MT-")
pdf("raw_stats_pool1.pdf",height=10,width=8)
VlnPlot(pbmc.pool,features = c("nFeature_RNA","nCount_RNA","percent.mt"),pt.size=0.1,ncol=1)
dev.off()

pbmc.pool <- subset(x = pbmc.pool, subset = nFeature_RNA > 100 & nFeature_RNA < 2500 & nCount_RNA < 10000 & percent.mt < 15)
pbmc.pool <- subset(pbmc.pool, idents=c("300BCG064", "300BCG096", "300BCG197"))

pbmc.pool$time="before_vaccination"
pbmc.pool$stim="S.aureus"
pbmc.pool$stim[pbmc.pool$Sample_assign=="300BCG197"]="RPMI"

pbmc.pool$pool = "1"
pbmc.pool1 <- pbmc.pool

######### pool 2 ##########
pbmc.pool <- ReadBowen(bdir="./pool_2",bid="pool2",btype="pool2",bmeta=pool2.info)
Idents(pbmc.pool)<-"Sample_assign"
pbmc.pool[["percent.mt"]] <- PercentageFeatureSet(object = pbmc.pool, pattern = "^MT-")
pdf("raw_stats_pool2.pdf",height=10,width=8)
VlnPlot(pbmc.pool,features = c("nFeature_RNA","nCount_RNA","percent.mt"),pt.size=0.1,ncol=1)
dev.off()

pbmc.pool <- subset(x = pbmc.pool, subset = nFeature_RNA > 100 & nFeature_RNA < 2500 & nCount_RNA < 10000 & percent.mt < 15)
pbmc.pool <- subset(pbmc.pool, idents=c("300BCG064", "300BCG096", "300BCG197"))

pbmc.pool$time="3m_after_vaccination"
pbmc.pool$stim="S.aureus"
pbmc.pool$stim[pbmc.pool$Sample_assign=="300BCG197"]="RPMI"

pbmc.pool$pool = "2"
pbmc.pool2 <- pbmc.pool


######### pool 3 ##########
pbmc.pool <- ReadBowen(bdir="./pool_3",bid="pool3",btype="pool3",bmeta=pool3.info)
Idents(pbmc.pool)<-"Sample_assign"
pbmc.pool[["percent.mt"]] <- PercentageFeatureSet(object = pbmc.pool, pattern = "^MT-")
#p0 <- VlnPlot(pbmc.pool,features = c("nFeature_RNA","nCount_RNA","percent.mt"),pt.size=0)
pdf("raw_stats_pool3.pdf",height=10,width=8)
VlnPlot(pbmc.pool,features = c("nFeature_RNA","nCount_RNA","percent.mt"),pt.size=0.1,ncol=1)
dev.off()

pbmc.pool <- subset(x = pbmc.pool, subset = nFeature_RNA > 100 & nFeature_RNA < 2500 & nCount_RNA < 10000 & percent.mt < 15)
pbmc.pool <- subset(pbmc.pool, idents=c("300BCG064", "300BCG096", "300BCG197"))

pbmc.pool$time="before_vaccination"
#pbmc.pool$time[pbmc.pool1$Sample_assign=="300BCG096"]=""
#pbmc.pool$time[pbmc.pool1$Sample_assign=="300BCG197"]=""

pbmc.pool$stim="RPMI"
pbmc.pool$stim[pbmc.pool$Sample_assign=="300BCG096"]="LPS"
pbmc.pool$stim[pbmc.pool$Sample_assign=="300BCG197"]="S.aureus"

pbmc.pool$pool = "3"
pbmc.pool3 <- pbmc.pool

######### pool 4 ##########
pbmc.pool <- ReadBowen(bdir="./pool_4",bid="pool4",btype="pool4",bmeta=pool4.info)
Idents(pbmc.pool)<-"Sample_assign"
pbmc.pool[["percent.mt"]] <- PercentageFeatureSet(object = pbmc.pool, pattern = "^MT-")
#p0 <- VlnPlot(pbmc.pool,features = c("nFeature_RNA","nCount_RNA","percent.mt"),pt.size=0)
pdf("raw_stats_pool4.pdf",height=10,width=8)
VlnPlot(pbmc.pool,features = c("nFeature_RNA","nCount_RNA","percent.mt"),pt.size=0.1,ncol=1)
dev.off()

pbmc.pool <- subset(x = pbmc.pool, subset = nFeature_RNA > 100 & nFeature_RNA < 2500 & nCount_RNA < 10000 & percent.mt < 15)
pbmc.pool <- subset(pbmc.pool, idents=c("300BCG064", "300BCG096", "300BCG197"))

pbmc.pool$time="3m_after_vaccination"
#pbmc.pool$time[pbmc.pool1$Sample_assign=="300BCG096"]=""
#pbmc.pool$time[pbmc.pool1$Sample_assign=="300BCG197"]=""

pbmc.pool$stim="RPMI"
pbmc.pool$stim[pbmc.pool$Sample_assign=="300BCG096"]="LPS"
pbmc.pool$stim[pbmc.pool$Sample_assign=="300BCG197"]="S.aureus"

pbmc.pool$pool = "4"
pbmc.pool4 <- pbmc.pool


######### pool 5 ##########
pbmc.pool <- ReadBowen(bdir="./pool_5",bid="pool5",btype="pool5",bmeta=pool5.info)
Idents(pbmc.pool)<-"Sample_assign"
pbmc.pool[["percent.mt"]] <- PercentageFeatureSet(object = pbmc.pool, pattern = "^MT-")
#p0 <- VlnPlot(pbmc.pool,features = c("nFeature_RNA","nCount_RNA","percent.mt"),pt.size=0)
pdf("raw_stats_pool5.pdf",height=10,width=8)
VlnPlot(pbmc.pool,features = c("nFeature_RNA","nCount_RNA","percent.mt"),pt.size=0.1,ncol=1)
dev.off()

pbmc.pool <- subset(x = pbmc.pool, subset = nFeature_RNA > 100 & nFeature_RNA < 2500 & nCount_RNA < 10000 & percent.mt < 15)
pbmc.pool <- subset(pbmc.pool, idents=c("300BCG064", "300BCG096", "300BCG197"))

pbmc.pool$time="before_vaccination"
#pbmc.pool$time[pbmc.pool1$Sample_assign=="300BCG096"]=""
#pbmc.pool$time[pbmc.pool1$Sample_assign=="300BCG197"]=""

pbmc.pool$stim="LPS"
pbmc.pool$stim[pbmc.pool$Sample_assign=="300BCG096"]="RPMI"
#pbmc.pool$stim[pbmc.pool$Sample_assign=="300BCG197"]="S.aureus"

pbmc.pool$pool = "5"
pbmc.pool5 <- pbmc.pool



######### pool 6 ##########
pbmc.pool <- ReadBowen(bdir="./pool_6",bid="pool6",btype="pool6",bmeta=pool6.info)
Idents(pbmc.pool)<-"Sample_assign"
pbmc.pool[["percent.mt"]] <- PercentageFeatureSet(object = pbmc.pool, pattern = "^MT-")
#p0 <- VlnPlot(pbmc.pool,features = c("nFeature_RNA","nCount_RNA","percent.mt"),pt.size=0)
pdf("raw_stats_pool5.pdf",height=10,width=8)
VlnPlot(pbmc.pool,features = c("nFeature_RNA","nCount_RNA","percent.mt"),pt.size=0.1,ncol=1)
dev.off()

pbmc.pool <- subset(x = pbmc.pool, subset = nFeature_RNA > 100 & nFeature_RNA < 2500 & nCount_RNA < 10000 & percent.mt < 15)
pbmc.pool <- subset(pbmc.pool, idents=c("300BCG064", "300BCG096", "300BCG197"))

pbmc.pool$time="3m_after_vaccination"
#pbmc.pool$time[pbmc.pool1$Sample_assign=="300BCG096"]=""
#pbmc.pool$time[pbmc.pool1$Sample_assign=="300BCG197"]=""

pbmc.pool$stim="LPS"
pbmc.pool$stim[pbmc.pool$Sample_assign=="300BCG096"]="RPMI"
#pbmc.pool$stim[pbmc.pool$Sample_assign=="300BCG197"]="S.aureus"

pbmc.pool$pool = "6"
pbmc.pool6 <- pbmc.pool



pbmc.pools= merge(x=pbmc.pool1, y=c(pbmc.pool2,pbmc.pool3,pbmc.pool4,pbmc.pool5,pbmc.pool6),project = "10X", add.cell.ids = c("p1","p2","p3","p4","p5","p6"))

pbmc.pools <- NormalizeData(object = pbmc.pools)
pbmc.pools <- FindVariableFeatures(object = pbmc.pools, selection.method = "vst", nfeatures = 2000)
pbmc.pools <- ScaleData(object = pbmc.pools, features = VariableFeatures(pbmc.pools))
pbmc.pools <- RunPCA(pbmc.pools, verbose = FALSE)
pbmc.pools <- RunUMAP(pbmc.pools, dims = 1:20)
pbmc.pools <- FindNeighbors(object = pbmc.pools, reduction = "pca", dims = 1:20) #used to be 18
pbmc.pools <- FindClusters(pbmc.pools, resolution = 0.4)



mks <- c("CCR7","TCF7","IL7R","GZMK",
         "CD8A","NKG7","GZMB","CD14","LYZ","FCGR3A","CST3",
         "CD86","CD163","TNF","IL1B","CD79A","CD27","SDC1","PPBP")
get.violin.data1 <- function(seurat, genes) {
  #  output = data.frame(gene = character(0), value= numeric(0), ident = character(0))
  #####################################################################################################
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
      ############################################
      #data.melt$tech=seurat$stim # ???What was used to stimulate the cells?
      ############################################
     # if(any(data.melt$value != 0)) noise = rnorm(length(data.melt$value))/100000 else noise = 0
      data.melt$value=as.numeric(as.character(data.melt$value))#+noise
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

violin.plot.data <- get.violin.data1(pbmc.pools, mks)
#pdf("marker_cluster2.pdf",width=12,height=6)
pv <- ggplot(violin.plot.data, aes(factor(ident),value)) +
  geom_violin(scale="width",adjust=1,trim=TRUE,aes(fill=factor(ident)),show.legend = F) + 
  ylab("") + xlab("") +
  coord_flip() +facet_wrap(~ gene,scales = "free_x", ncol = length(levels(violin.plot.data$gene))) +
  theme(strip.text.x = element_text(size=12, angle=80),
        axis.text.y = element_text(size=15),
        strip.background = element_blank(),
        panel.spacing.x = unit(c(-0.2), "lines"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0.5,0.4,0.5,1.3), "cm")) +
  scale_x_discrete(limits = rev(levels(violin.plot.data$ident)), position = "top")
#dev.off()

p1 <- DimPlot(pbmc.pools,label=T)


px <- plot_grid(p0,p1)
dev.off()
plot_grid(px,pv,ncol=1)




cluster.markers.wilcox <- FindAllMarkers(object = pbmc.pools,
only.pos = TRUE,
min.pct = 0.2,
min.diff.pct = 0.2,
logfc.threshold = 0.25,
test.use = "wilcox")


top5 <- cluster.markers.wilcox %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
a <- DotPlot(pbmc.pools,features=rev(unique(top5$gene)),cols="RdBu")+coord_flip()
p2 <- a + theme(legend.position = "top",legend.title = element_blank())+xlab("")+ylab("")

px <- plot_grid(p1,pv,ncol=1)


Idents(pbmc.pools) <- "seurat_clusters"
pbmc.pools <- RenameIdents(pbmc.pools,
`0` = "CD4 T",
`1` = "NK",
`2` = "B cells",
`3` = "CD8 T",
`4` = "CD4 T",
`5` = "Monocytes",
`6` = "CD4 T",
`7` = "CD4 T",
`8` = "CD4 T",
`9` = "NK",
`10` = "B cells",
`11` = "Monocytes",
`12` = "CD4 T",
`13` = "Monocytes",
`14` = "Platelet",
`15` = "Monocytes",
`16` = "Monocytes"
)
pbmc.pools$pre.celltype <- Idents(pbmc.pools)



Idents(pbmc.pools) <- "seurat_clusters"
pbmc.pools <- RenameIdents(pbmc.pools,
`0` = "0_CD4 T",
`4` = "4_CD4 T",
`6` = "6_CD4 T",
`7` = "7_CD4 T",
`8` = "8_CD4 T",
`12` = "12_CD4 T",
`3` = "3_CD8 T",
`1` = "1_NK",
`9` = "9_NK",
`2` = "2_B",
`10` = "10_B",
`5` = "5_Monocytes",
`11` = "11_Monocytes",
`13` = "13_Monocytes",
`15` = "15_Monocytes",
`16` = "16_Monocytes",
`14` = "14_Platelet"
)
pbmc.pools$celltype <- Idents(pbmc.pools)

cols <- c(
"0_CD4 T" = "#800026",
"1_NK" = "#7bccc4",
"2_B" = "#bf812d",
"3_CD8 T" = "#084081",
"4_CD4 T" = "#bd0026",
"5_Monocytes" = "#4d9221",
"6_CD4 T" = "#e31a1c",
"7_CD4 T" = "#fc4e2a",
"8_CD4 T" = "#fd8d3c",
"9_NK" = "#2b8cbe",
 "10_B"= "#8c510a",
 "11_Monocytes" = "#de77ae",
 "12_CD4 T" = "#feb24c",
 "13_Monocytes"= "#f1b6da",
 "14_Platelet"= "#878787",
 "15_Monocytes"= "#006837",
 "16_Monocytes" = "#b8e186"
)
DimPlot(pbmc.pools, split.by="time.stim",ncol=2, cols=cols)




### monocytes subsets ####
Idents(pbmc.pools) <- "seurat_clusters"
mono.pools <- subset(pbmc.pools,idents=c("5","11","13","15","16"))

features <- c("TNF","IL1B","IL6","CXCL9","CXCL10","CXCL11")
Idents(mono.pools) <- "seurat_clusters"
DefaultAssay(mono.pools) <- "RNA"



##### TI subsets assignments ####
mpc.mk <- c("IL1B"  ,   "IL8"   ,   "IL6"    ,  "PTGS2"  ,  "IL1A" ,    "CCL20"  ,  "TNF"   ,   "CXCL3" ,   "CXCL1")
mtc.mk <- c("CXCL11", "CXCL10" ,  "CXCL9"  ,  "TNFSF10" , "HLA-DQA1", "FCN1"   ,  "IGFBP4" ,  "HLA-DPB1", "HLA-DQB1", "FAM26F" ,  "RGL1","CD4")

tmp<-as.data.frame(as.matrix(mono.pools@assays$RNA@data))
filtered_matrix<-rowSums(tmp > 0)
filtered_matrix<-filtered_matrix[filtered_matrix>10]
tmp <- tmp[row.names(tmp) %in% names(filtered_matrix),]
cells_rankings <- AUCell_buildRankings(as.matrix(tmp))

##################
cells_AUC <- AUCell_calcAUC(mpc.mk, cells_rankings, aucMaxRank = ceiling(0.03 * nrow(cells_rankings)), normAUC = T)
#Visualization
tmp <- getAUC(cells_AUC)
tmp <- as.data.frame(t(tmp))
df <- data.frame(row.names = row.names(mono.pools@meta.data), cluster = mono.pools$seurat_clusters, stringsAsFactors = F)
df <- as.data.frame(merge(df, tmp, by = 0))
row.names(df) <- df$Row.names
df$Row.names <- NULL
df.umap <- as.data.frame(mono.pools@reductions$umap@cell.embeddings)
df <- merge(df, df.umap, by=0)
df$geneset<-df$geneSet
# df <- df[df$cluster != "3",]
gg1 <- ggplot(df, aes(y = UMAP_2, x = UMAP_1, color = geneset)) + geom_point(alpha=1, size=0.3)+
  #geom_point_rast(alpha = 1, size = .1) +
  scale_colour_gradient2(low = "grey50", mid = "seashell", high = "red", midpoint = 0.1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color="black", size = .3),
        axis.line.y = element_line(color="black", size = .3),
        legend.position = "none")+ggtitle("MCI signatures")

cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)
cells_assignment$geneSet$aucThr$selected
mpc_assign <- cells_assignment


###################
cells_AUC <- AUCell_calcAUC(mtc.mk, cells_rankings, aucMaxRank = ceiling(0.03 * nrow(cells_rankings)), normAUC = T)
#Visualization
tmp <- getAUC(cells_AUC)
tmp <- as.data.frame(t(tmp))
df <- data.frame(row.names = row.names(mono.pools@meta.data), cluster = mono.pools$seurat_clusters, stringsAsFactors = F)
df <- as.data.frame(merge(df, tmp, by = 0))
row.names(df) <- df$Row.names
df$Row.names <- NULL
df.umap <- as.data.frame(mono.pools@reductions$umap@cell.embeddings)
df <- merge(df, df.umap, by=0)
df$geneset<-df$geneSet
# df <- df[df$cluster != "3",]
gg0 <- ggplot(df, aes(y = UMAP_2, x = UMAP_1, color = geneset)) + geom_point(alpha=1, size=0.3)+
  #geom_point_rast(alpha = 1, size = .1) +
  scale_colour_gradient2(low = "grey50", mid = "seashell", high = "red", midpoint = 0.1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color="black", size = .3),
        axis.line.y = element_line(color="black", size = .3),
        legend.position = "none")+ggtitle("MC signatures")

cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)
cells_assignment$geneSet$aucThr$selected
mtc_assign <- cells_assignment



#mono.pools$cellnames <- colnames(mono.pools)
mono.pools$cellnames <- colnames(mono.pools)

mono.pools$trainSig <- "NT"
mono.pools$trainSig[mono.pools$cellnames %in% mpc_assign$geneSet$assignment] = "MCI"
mono.pools$trainSig[mono.pools$cellnames %in% mtc_assign$geneSet$assignment] = "MC"


