### loading and clustering ###

library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(Matrix.utils)
library(reshape2)
library(pryr)
library(RColorBrewer)
library(cowplot)



setwd("/path/to/count/projects")
ReadBowen <- function(fid,individual,stimu,time,cells,donor){
	ff <- paste("/path/to/count/table",fid,".seurat.txt",sep="")
	datas <- read.table(ff,sep="\t")
	ERCC.WT.index <- grep(pattern = "^ERCC-", x = rownames(datas), value = FALSE)
	percent.ERCC.WT <- Matrix::colSums(datas[ERCC.WT.index, ])/Matrix::colSums(datas)
	datas <- datas[-ERCC.WT.index, ]
	bsample <- CreateSeuratObject(counts = datas, project = "IMMUNE", min.cells = 3, meta.data = data.frame(percent.ercc = percent.ERCC.WT))
	bsample$bsampid <- fid
	bsample$individual <- individual
	bsample$stimu <- stimu
	bsample$time <- time
	bsample$cells <- cells
    bsample$donor <- donor
	return(bsample)
}


d1t1 <- ReadBowen("A-S1-T1-M","A-S1-T1-M","S1","T1","MONO","A")
d1t2 <- ReadBowen("B-S1-T1-M","B-S1-T1-M","S1","T1","MONO","B")
d1t3 <- ReadBowen("C-S1-T1-M","C-S1-T1-M","S1","T1","MONO","C")

d2t1 <- ReadBowen("A-S2-T1-M","A-S2-T1-M","S2","T1","MONO","A")
d2t2 <- ReadBowen("B-S2-T1-M","B-S2-T1-M","S2","T1","MONO","B")
d2t3 <- ReadBowen("C-S2-T1-M","C-S2-T1-M","S2","T1","MONO","C")

d3t1 <- ReadBowen("A-S3-T1-M","A-S3-T1-M","S3","T1","MONO","A")
d3t2 <- ReadBowen("B-S3-T1-M","B-S3-T1-M","S3","T1","MONO","B")
d3t3 <- ReadBowen("C-S3-T1-M","C-S3-T1-M","S3","T1","MONO","C")

d4t1 <- ReadBowen("A-S4-T1-M","A-S4-T1-M","S4","T1","MONO","A")
d4t2 <- ReadBowen("B-S4-T1-M","B-S4-T1-M","S4","T1","MONO","B")
d4t3 <- ReadBowen("C-S4-T1-M","C-S4-T1-M","S4","T1","MONO","C")

d5t1 <- ReadBowen("A-S5-T1-M","A-S5-T1-M","S5","T1","MONO","A")
d5t2 <- ReadBowen("B-S5-T1-M","B-S5-T1-M","S5","T1","MONO","B")
d5t3 <- ReadBowen("C-S5-T1-M","C-S5-T1-M","S5","T1","MONO","C")

d6t1 <- ReadBowen("A-S6-T1-M","A-S6-T1-M","S6","T1","MONO","A")
d6t2 <- ReadBowen("B-S6-T1-M","B-S6-T1-M","S6","T1","MONO","B")
d6t3 <- ReadBowen("C-S6-T1-M","C-S6-T1-M","S6","T1","MONO","C")


t1t1 <- ReadBowen("A-S1-T2-M","A-S1-T2-M","S1","T2","MONO","A")
t1t2 <- ReadBowen("B-S1-T2-M","B-S1-T2-M","S1","T2","MONO","B")
t1t3 <- ReadBowen("C-S1-T2-M","C-S1-T2-M","S1","T2","MONO","C")
t1t1g <- ReadBowen("A-S1-T2-GateM","A-S1-T2-M","S1","T2","MONO","A")
t1t2g <- ReadBowen("B-S1-T2-GateM","B-S1-T2-M","S1","T2","MONO","B")

t2t1 <- ReadBowen("A-S2-T2-M","A-S2-T2-M","S2","T2","MONO","A")
t2t2 <- ReadBowen("B-S2-T2-M","B-S2-T2-M","S2","T2","MONO","B")
t2t3 <- ReadBowen("C-S2-T2-M","C-S2-T2-M","S2","T2","MONO","C")
t2t1g <- ReadBowen("A-S2-T2-GateM","A-S2-T2-M","S2","T2","MONO","A")
t2t2g <- ReadBowen("B-S2-T2-GateM","B-S2-T2-M","S2","T2","MONO","B")

t3t1 <- ReadBowen("A-S3-T2-M","A-S3-T2-M","S3","T2","MONO","A")
t3t2 <- ReadBowen("B-S3-T2-M","B-S3-T2-M","S3","T2","MONO","B")
t3t3 <- ReadBowen("C-S3-T2-M","C-S3-T2-M","S3","T2","MONO","C")
t3t1g <- ReadBowen("A-S3-T2-GateM","A-S3-T2-M","S3","T2","MONO","A")
t3t2g <- ReadBowen("B-S3-T2-GateM","B-S3-T2-M","S3","T2","MONO","B")


t4t1 <- ReadBowen("A-S4-T2-M","A-S4-T2-M","S4","T2","MONO","A")
t4t2 <- ReadBowen("B-S4-T2-M","B-S4-T2-M","S4","T2","MONO","B")
t4t3 <- ReadBowen("C-S4-T2-M","C-S4-T2-M","S4","T2","MONO","C")
t4t1g <- ReadBowen("A-S4-T2-GateM","A-S4-T2-M","S4","T2","MONO","A")
t4t2g <- ReadBowen("B-S4-T2-GateM","B-S4-T2-M","S4","T2","MONO","B")


t5t1 <- ReadBowen("A-S5-T2-M","A-S5-T2-M","S5","T2","MONO","A")
t5t2 <- ReadBowen("B-S5-T2-M","B-S5-T2-M","S5","T2","MONO","B")
t5t3 <- ReadBowen("C-S5-T2-M","C-S5-T2-M","S5","T2","MONO","C")
t5t1g <- ReadBowen("A-S5-T2-GateM","A-S5-T2-M","S5","T2","MONO","A")
t5t2g <- ReadBowen("B-S5-T2-GateM","B-S5-T2-M","S5","T2","MONO","B")

t6t1 <- ReadBowen("A-S6-T2-M","A-S6-T2-M","S6","T2","MONO","A")
t6t2 <- ReadBowen("B-S6-T2-M","B-S6-T2-M","S6","T2","MONO","B")
t6t3 <- ReadBowen("C-S6-T2-M","C-S6-T2-M","S6","T2","MONO","C")
t6t1g <- ReadBowen("A-S6-T2-GateM","A-S6-T2-M","S6","T2","MONO","A")

d1p1 <- ReadBowen("A-S1-T1-P","A-S1-T1-P","S1","T1","PBMC","A")
d1p2 <- ReadBowen("B-S1-T1-P","B-S1-T1-P","S1","T1","PBMC","B")
d1p3 <- ReadBowen("C-S1-T1-P","C-S1-T1-P","S1","T1","PBMC","C")

d2p1 <- ReadBowen("A-S2-T1-P","A-S2-T1-P","S2","T1","PBMC","A")
d2p2 <- ReadBowen("B-S2-T1-P","B-S2-T1-P","S2","T1","PBMC","B")
d2p3 <- ReadBowen("C-S2-T1-P","C-S2-T1-P","S2","T1","PBMC","C")

d3p1 <- ReadBowen("A-S3-T1-P","A-S3-T1-P","S3","T1","PBMC","A")
d3p2 <- ReadBowen("B-S3-T1-P","B-S3-T1-P","S3","T1","PBMC","B")
d3p3 <- ReadBowen("C-S3-T1-P","C-S3-T1-P","S3","T1","PBMC","C")

d4p1 <- ReadBowen("A-S4-T1-P","A-S4-T1-P","S4","T1","PBMC","A")
d4p2 <- ReadBowen("B-S4-T1-P","B-S4-T1-P","S4","T1","PBMC","B")
d4p3 <- ReadBowen("C-S4-T1-P","C-S4-T1-P","S4","T1","PBMC","C")

d5p1 <- ReadBowen("A-S5-T1-P","A-S5-T1-P","S5","T1","PBMC","A")
d5p2 <- ReadBowen("B-S5-T1-P","B-S5-T1-P","S5","T1","PBMC","B")
d5p3 <- ReadBowen("C-S5-T1-P","C-S5-T1-P","S5","T1","PBMC","C")

d6p1 <- ReadBowen("A-S6-T1-P","A-S6-T1-P","S6","T1","PBMC","A")
d6p2 <- ReadBowen("B-S6-T1-P","B-S6-T1-P","S6","T1","PBMC","B")
d6p3 <- ReadBowen("C-S6-T1-P","C-S6-T1-P","S6","T1","PBMC","C")


t1p1 <- ReadBowen("A-S1-T2-P","A-S1-T2-P","S1","T2","PBMC","A")
t1p2 <- ReadBowen("B-S1-T2-P","B-S1-T2-P","S1","T2","PBMC","B")
t1p3 <- ReadBowen("C-S1-T2-P","C-S1-T2-P","S1","T2","PBMC","C")

t2p1 <- ReadBowen("A-S2-T2-P","A-S2-T2-P","S2","T2","PBMC","A")
t2p2 <- ReadBowen("B-S2-T2-P","B-S2-T2-P","S2","T2","PBMC","B")
t2p3 <- ReadBowen("C-S2-T2-P","C-S2-T2-P","S2","T2","PBMC","C")

t3p1 <- ReadBowen("A-S3-T2-P","A-S3-T2-P","S3","T2","PBMC","A")
t3p2 <- ReadBowen("B-S3-T2-P","B-S3-T2-P","S3","T2","PBMC","B")
t3p3 <- ReadBowen("C-S3-T2-P","C-S3-T2-P","S3","T2","PBMC","C")

t4p1 <- ReadBowen("A-S4-T2-P","A-S4-T2-P","S4","T2","PBMC","A")
t4p2 <- ReadBowen("B-S4-T2-P","B-S4-T2-P","S4","T2","PBMC","B")
t4p3 <- ReadBowen("C-S4-T2-P","C-S4-T2-P","S4","T2","PBMC","C")

t5p1 <- ReadBowen("A-S5-T2-P","A-S5-T2-P","S5","T2","PBMC","A")
t5p2 <- ReadBowen("B-S5-T2-P","B-S5-T2-P","S5","T2","PBMC","B")
t5p3 <- ReadBowen("C-S5-T2-P","C-S5-T2-P","S5","T2","PBMC","C")

t6p1 <- ReadBowen("A-S6-T2-P","A-S6-T2-P","S6","T2","PBMC","A")
t6p2 <- ReadBowen("B-S6-T2-P","B-S6-T2-P","S6","T2","PBMC","B")
t6p3 <- ReadBowen("C-S6-T2-P","C-S6-T2-P","S6","T2","PBMC","C")



######merge as condition####
#without LPS##
groups1 <- merge(x=d1t1,c(d1t2,d1t3,d1p1,d1p2,d1p3,
t1t1,t1t1g,t1t2,t1t3,t1p1,t1p2,t1p3
), project="RPMI",add.cell.ids=
  c("As1t1m","Bs1t1m","Cs1t1m","As1t1p","Bs1t1p","Cs1t1p",
    "As1t2m","As1t2mg","Bs1t2m","Cs1t2m","As1t2p","Bs1t2p","Cs1t2p")
)


groups3 <- merge(x=d3t1,c(d3t2,d3t3,d3p1,d3p2,d3p3,
t3t1,t3t2g,t3t2,t3t3,t3p1,t3p2,t3p3
), project="BG",add.cell.ids=
  c("As3t1m","Bs3t1m","Cs3t1m","As3t1p","Bs3t1p","Cs3t1p",
    "As3t2m","Bs3t2mg","Bs3t2m","Cs3t2m","As3t2p","Bs3t2p","Cs3t2p")
)

groups4 <- merge(x=d4t1,c(d4t2,d4t3,d4p1,d4p2,d4p3,
t4t1,t4t2g,t4t2,t4t3,t4p1,t4p2,t4p3
), project="UA",add.cell.ids=
  c("As4t1m","Bs4t1m","Cs4t1m","As4t1p","Bs4t1p","Cs4t1p",
    "As4t2m","Bs4t2mg","Bs4t2m","Cs4t2m","As4t2p","Bs4t2p","Cs4t2p")
)

groups5 <- merge(x=d5t1,c(d5t2,d5t3,d5p1,d5p2,d5p3,
t5t1,t5t2g,t5t2,t5t3,t5p1,t5p2,t5p3
), project="oxLDL",add.cell.ids=
  c("As5t1m","Bs5t1m","Cs5t1m","As5t1p","Bs5t1p","Cs5t1p",
    "As5t2m","Bs5t2mg","Bs5t2m","Cs5t2m","As5t2p","Bs5t2p","Cs5t2p")
)

groups6 <- merge(x=d6t1,c(d6t2,d6t3,d6p1,d6p2,d6p3,
t6t1,t6t2,t6t3,t6p1,t6p2,t6p3
), project="MDP",add.cell.ids=
  c("As6t1m","Bs6t1m","Cs6t1m","As6t1p","Bs6t1p","Cs6t1p",
    "As6t2m","Bs6t2m","Cs6t2m","As6t2p","Bs6t2p","Cs6t2p")
)


groups1.SCT <- subset(x = groups1, subset = nFeature_RNA > 100 & nFeature_RNA < 7000 & percent.mt < 25 & percent.ercc < 0.9)
groups1.SCT <- SCTransform(groups1.SCT, vars.to.regress = "percent.mt", verbose = FALSE)  # vars.to.regress = "percent.mt",

groups3.SCT <- subset(x = groups3, subset = nFeature_RNA > 100 & nFeature_RNA < 7000 & percent.mt < 25 & percent.ercc < 0.9)
groups3.SCT <- SCTransform(groups3.SCT, vars.to.regress = "percent.mt", verbose = FALSE)  # vars.to.regress = "percent.mt",

groups4.SCT <- subset(x = groups4, subset = nFeature_RNA > 100 & nFeature_RNA < 7000 & percent.mt < 25 & percent.ercc < 0.9)
groups4.SCT <- SCTransform(groups4.SCT, vars.to.regress = "percent.mt", verbose = FALSE)  # vars.to.regress = "percent.mt",

groups5.SCT <- subset(x = groups5, subset = nFeature_RNA > 100 & nFeature_RNA < 7000 & percent.mt < 25 & percent.ercc < 0.9)
groups5.SCT <- SCTransform(groups5.SCT, vars.to.regress = "percent.mt", verbose = FALSE)  # vars.to.regress = "percent.mt",

groups6.SCT <- subset(x = groups6, subset = nFeature_RNA > 100 & nFeature_RNA < 7000 & percent.mt < 25 & percent.ercc < 0.9)
groups6.SCT <- SCTransform(groups6.SCT, vars.to.regress = "percent.mt", verbose = FALSE)  # vars.to.regress = "percent.mt",

panlist <- list(groups1.SCT,groups3.SCT,groups4.SCT, groups5.SCT ,groups6.SCT)

monoT2.features <- SelectIntegrationFeatures(object.list = panlist, nfeatures = 2000)
panlist <- PrepSCTIntegration(object.list = panlist, anchor.features = monoT2.features, verbose = FALSE)

monoT2.anchors <- FindIntegrationAnchors(object.list = panlist, normalization.method = "SCT", anchor.features = monoT2.features, verbose = FALSE, k.filter = 90)  #k.filter was set to  90, because some groups have no more than 96 cells.
monoT2.integrated <- IntegrateData(anchorset = monoT2.anchors, normalization.method = "SCT", verbose = FALSE)
monoT2.integrated <- RunPCA(monoT2.integrated, verbose = FALSE)
monoT2.integrated <- RunUMAP(monoT2.integrated, dims = 1:20)
monoT2.integrated <- FindNeighbors(object = monoT2.integrated, reduction = "pca", dims = 1:20) #used to be 18
monoT2. mono.pbmc.raw <- FindClusters(monoT2.integrated, resolution = 0.6)


Idents(mono.pbmc.raw)<-"stimu"
mono.pbmc.raw <- RenameIdents(mono.pbmc.raw,
`S1` = "RPMI",
`S3` = "BG",
`S4` = "UA",
`S5` = "oxLDL",
`S6` = "MDP"
)
mono.pbmc.raw$condition <- Idents(mono.pbmc.raw)


allsamples <- levels(as.factor(mono.pbmc.raw$bsampid))
allsamples <- allsamples[-1] #remove A-S1-T1-M
mono.pbmc.raw.flt <- subset(mono.pbmc.raw,idents=allsamples)
mono.pbmc.raw.flt <- RunPCA(mono.pbmc.raw.flt, verbose = FALSE)
mono.pbmc.raw.flt <- RunUMAP(mono.pbmc.raw.flt, dims = 1:20)
mono.pbmc.raw.flt <- FindNeighbors(object = mono.pbmc.raw.flt, reduction = "pca", dims = 1:20) #used to be 18
mono.pbmc.raw.flt <- FindClusters(mono.pbmc.raw.flt, resolution = 0.6)

DimPlot(mono.pbmc.raw.flt, group.by="time",split.by="donor",cols=c("blue","red"))

mono.pbmc.raw.flt$time.stim <- paste(mono.pbmc.raw.flt$time, mono.pbmc.raw.flt$stimu,sep=".")
twocolor10 <- c("grey","grey",
"grey","grey",
"grey","grey",
"grey","grey",
"red","blue"
)
d1<-DimPlot(mono.pbmc.raw.flt,group.by="time.stim",cols=twocolor10,order=c("T1.S1","T2.S1"))+NoLegend()+ggtitle("RPMI")
d3<-DimPlot(mono.pbmc.raw.flt,group.by="time.stim",cols=twocolor10,order=c("T1.S3","T2.S3"))+NoLegend()+ggtitle("BG")
d4<-DimPlot(mono.pbmc.raw.flt,group.by="time.stim",cols=twocolor10,order=c("T1.S4","T2.S4"))+NoLegend()+ggtitle("UA")
d5<-DimPlot(mono.pbmc.raw.flt,group.by="time.stim",cols=twocolor10,order=c("T1.S5","T2.S5"))+NoLegend()+ggtitle("oxLDL")
d6<-DimPlot(mono.pbmc.raw.flt,group.by="time.stim",cols=twocolor10,order=c("T1.S6","T2.S6"))+NoLegend()+ggtitle("MDP")
plot_grid(d1,d3,d4,d5,d6)


DefaultAssay(mono.pbmc.raw.flt)<-"RNA"
mono.pbmc.raw.flt <- NormalizeData(mono.pbmc.raw.flt)
markers <- c("CD14","CCL2","S100A9","FCGR3A","VMO1","CD86","ITGAX","HLA-DRA", "HLA-DQB1","CD74","TNF","CXCL1","IL1A","IL6","MSR1","IL1RN","LPL","SPP1","FBP1","CXCL10","CXCL11")

markers <- c("CD14","CCL2","S100A9","FCGR3A","VMO1","CD86","ITGAX","HLA-DRA","HLA-DQB1","TNF","IL1B","IL6","IL1A","MSR1","IL1RN","LPL","SPP1","CXCL9","CXCL10","CXCL11")

get.violin.data1 <- function(seurat, genes) {
  #  output = data.frame(gene = character(0), value= numeric(0), ident = character(0))
  #####################################################################################################
  output = data.frame(gene = character(0), value= numeric(0), ident = character(0), tech=character(0))
  #####################################################################################################
  for (gene in genes) {
    data.use = data.frame(FetchData(seurat,gene))
    data.use = t(data.use)
    data.melt=data.frame(rep(gene, length(seurat@active.ident)))
    colnames(data.melt)[1]="gene"
    data.melt$value=as.numeric(data.use[1,1:length(seurat@active.ident)])
    data.melt$id=names(data.use)[1:length(seurat@active.ident)]
    data.melt$ident=seurat@active.ident
    ############################################
    data.melt$tech=seurat$stim
    ############################################
    noise = rnorm(length(data.melt$value))/100000
    data.melt$value=as.numeric(as.character(data.melt$value))+noise
    output = rbind(output, data.melt)
  }
  return(output)
}
###############
violin.plot.data <- get.violin.data1(luis, markers)
ggplot(violin.plot.data, aes(factor(ident),value)) +
  geom_violin(scale="width",adjust=1,trim=TRUE,aes(fill=factor(ident)),show.legend = F) +
  ylab("") + xlab("") +
  coord_flip() +facet_wrap(~ gene,scales = "free_x", ncol = length(levels(violin.plot.data$gene))) +
  theme(strip.text.x = element_text(size=12, angle=50),
        axis.text.y = element_text(size=15),
        strip.background = element_blank(),
        panel.spacing.x = unit(c(-0.2), "lines"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0.5,0.4,0.5,1.3), "cm")) +
  scale_x_discrete(limits = rev(levels(violin.plot.data$ident)), position = "top")

Idents(mono.pbmc.raw.flt) <- "seurat_clusters"
mono.pbmc.raw.flt <- RenameIdents(mono.pbmc.raw.flt,
`3` = "Classical Monocytes",
`0` = "Intermediate Monocytes",
`6` = "Nonclassical Monocytes",
`9` = "Monocytic DC",
`7` = "HIF-1 signaling cells",
`10` = "Antigen-signaling cells",
`5` = "Resting cells",
`8` = "UGDH-AS1+ cell",
`1` = "Macrophages-1",
`2` = "Macrophages-2",
`4` = "Unpolarized"
)
mono.pbmc.raw.flt$celltype.sort <- Idents(mono.pbmc.raw.flt)
Idents(mono.pbmc.raw.flt) <- "seurat_clusters"
mono.pbmc.raw.flt <- RenameIdents(mono.pbmc.raw.flt,
`0` = "ITM",
`1` = "HLA",
`2` = "MM2",
`3` = "CLM",
`4` = "UPO",
`5` = "UNK",
`6` = "NCM",
`7` = "PB1",
`8` = "UAC",
`9` = "MDC",
`10` = "PB2"
)
mono.pbmc.raw.flt$celltype.short <- Idents(mono.pbmc.raw.flt)


Idents(mono.pbmc.raw.flt)<-"time"
mono.pbmc.t2 <- subset(mono.pbmc.raw.flt,idents='T2')
mono.pbmc.t1 <- subset(mono.pbmc.raw.flt,idents='T1')


