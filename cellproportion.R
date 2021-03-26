library(Seurat)
library(ggplot2)

library(dplyr)
library(Matrix)
library(Matrix.utils)
library(reshape2)
library(pryr)
library(RColorBrewer)
library(cowplot)

mono.pbmc.t2 <- readRDS("mono.pbmc.t2.rds")
mono.t2 <- subset(mono.pbmc.t2, condition != "RPMI")


#######################
##  cell proportion
#######################
Idents(mono.t2) <- "newTI"
tmp<-data.frame(cell_type= mono.t2@active.ident,sample= mono.t2@meta.data$bsampid,donor= mono.t2@meta.data$cells, stim=mono.t2@meta.data$condition)
res<-matrix(NA,ncol=5,nrow=length(unique(tmp$cell_type))*length(unique(tmp$sample)))
res<-as.data.frame(res)
names(res)<-c("cell_type","freq","sample","Condition","Stimulation")
n=1;
for (d in unique(tmp$sample)){
  t1<-subset(tmp,sample==d)
  t2<-as.data.frame(table(t1$cell_type))
  res[n:(n+length(unique(tmp$cell_type))-1),1]<-as.character(t2$Var1)
  res[n:(n+length(unique(tmp$cell_type))-1),2]<-t2$Freq
  res[n:(n+length(unique(tmp$cell_type))-1),3]<-rep(as.character(d),length(unique(tmp$cell_type)))
  res[n:(n+length(unique(tmp$cell_type))-1),4]<-rep(as.character(unique(t1$donor)),length(unique(tmp$cell_type)))
  res[n:(n+length(unique(tmp$cell_type))-1),5]<-rep(as.character(unique(t1$stim)),length(unique(tmp$cell_type)))
  n=n+length(unique(tmp$cell_type));
}
x<-split(res,res$sample)
y<-lapply(x, function(w) {w=transform(w,f=freq/sum(freq))})
y2<-do.call(rbind,y)




#######################
##  dirichlet test
#######################

library(DirichletReg)
head(y2)
y2$other <- 1-y2$f
y2$Condition <- factor(y2$Condition, levels=c("MONO","PBMC"))
y2$cells <- paste(y2$cell_type, y2$Stimulation, sep=".")

dat.test2 <- NULL
  i<-'MCI'
  print(i)
  y2.i <- y2[y2$cell_type==i,]
  y2.i$Smp <- DR_data(y2.i[,6:7])
  res <- DirichReg(Smp~Condition,y2.i,model="alternative",base=2)
  x <- summary(res)
  sc <- x$coef.mat["ConditionPBMC",3:4]
  dat.test2 <- rbind(dat.test2, c(i,sc,"MONO.vs.PBMC"))
  i<-'MC'
  print(i)
  y2.i <- y2[y2$cell_type==i,]
  y2.i$Smp <- DR_data(y2.i[,6:7])
  res <- DirichReg(Smp~Condition,y2.i,model="alternative",base=2)
  x <- summary(res)
  sc <- x$coef.mat["ConditionPBMC",3:4]
  dat.test2 <- rbind(dat.test2, c(i,sc,"MONO.vs.PBMC"))
  
 dat.test2<- as.data.frame(dat.test2)
names(dat.test2) <- c("cell_type","Z","Pr","comparison")
dat.test2$Pr <- as.numeric(as.character(dat.test2$Pr))
dat.test2$Z <- as.numeric(as.character(dat.test2$Z))

dat.test2$sig <- "ns"
dat.test2$sig[dat.test2$Pr < 0.05]="*"
dat.test2$sig[dat.test2$Pr < 0.01]="**"
dat.test2$sig[dat.test2$Pr < 0.001]="***"

dat.sig <- dat.test2[dat.test2$sig != "ns",]
dat.sig$st="MONO"
dat.sig$ed="PBMC"
dat.sig$sig <- round(dat.sig$Pr,digits = 4)
dat.sig$maxs <- 0.8

  
  p <- ggplot(data=y2, aes(x=Condition,y=f))+geom_boxplot(aes(fill=Condition),outlier.shape=NA)
p2 <- p + facet_grid(.~cell_type)+geom_jitter(aes(color= Stimulation), shape=16, width = 0.2) +theme_classic()+ theme(axis.text.x=element_text(angle = 60,vjust = 1,hjust = 1))+scale_fill_manual(values=c("MONO"="darkgrey","PBMC"="lightgrey"))+
  scale_color_manual(values=c("BG"="#67AD57","MDP"="#4A7CB3","UA"="#D1382C","oxLDL"="#E2E057"))+ylab("cell proportion")
p2 + geom_signif(data=dat.sig, aes(xmin=st,xmax=ed, annotations=sig, y_position=maxs),  manual = TRUE, vjust=0,tip_length = 0.01)  
#######################
#######################




dat.test <- NULL

for(i in levels(as.factor(y2$cells))){
#for(i in c("PLA1","PLA2","PLA3")){  
  print(i)
  y2.i <- y2[y2$cells==i,]
  y2.i$Smp <- DR_data(y2.i[,6:7])
  res <- DirichReg(Smp~Condition,y2.i,model="alternative",base=2)
  x <- summary(res)
  sc <- x$coef.mat["ConditionPBMC",3:4]
  dat.test <- rbind(dat.test, c(i,sc,"MONO.vs.PBMC"))
}
dat.test<- as.data.frame(dat.test)
names(dat.test) <- c("cell_type","Z","Pr","comparison")
dat.test$Pr <- as.numeric(as.character(dat.test$Pr))
dat.test$Z <- as.numeric(as.character(dat.test$Z))

dat.test$sig <- "ns"
dat.test$sig[dat.test$Pr < 0.05]="*"
dat.test$sig[dat.test$Pr < 0.01]="**"
dat.test$sig[dat.test$Pr < 0.001]="***"


p <- ggplot(data=y2, aes(x=Condition,y=f))+geom_boxplot(aes(fill=Stimulation),outlier.shape=NA)
p2 <- p + facet_grid(.~cells)+geom_jitter(shape=16, width = 0.2) +theme_classic()+ theme(axis.text.x=element_text(angle = 60,vjust = 1,hjust = 1))

dat.test <- NULL

for(i in levels(as.factor(y2$cell_type))){
  print(i)
  y2.i <- y2[y2$cell_type==i,]
  y2.i$Smp <- DR_data(y2.i[,5:6])
  res <- DirichReg(Smp~stim,y2.i,model="alternative",base=2)
  x <- summary(res)
  sc <- x$coef.mat["stimSC",3:]
  hc <- x$coef.mat["stimHC",3:4]
  dat.test <- rbind(dat.test, c(i,sc,"UCvsSC"))
  dat.test <- rbind(dat.test, c(i,hc,"UCvsHC"))
}
dat.test<- as.data.frame(dat.test)
names(dat.test) <- c("cells","Z","Pr","comparison")
dat.test$Pr <- as.numeric(as.character(dat.test$Pr))
dat.test$Z <- as.numeric(as.character(dat.test$Z))

dat.test$sig <- "ns"
dat.test$sig[dat.test$Pr < 0.05]="*"
dat.test$sig[dat.test$Pr < 0.01]="**"
dat.test$sig[dat.test$Pr < 0.001]="***"

#dat.sig <- dat.test[dat.test$sig != "ns",]
dat.sig$st="MONO"
dat.sig$ed="PBMC"
dat.sig$sig <- round(dat.sig$Pr,digits = 4)
maxs <- NULL
for(i in dat.sig$cells){maxs <- c(maxs,max(y2$f[y2$cells == i])+0.05)}
dat.sig$maxs <- maxs

p <- ggplot(data=y2, aes(x=Condition,y=f))+geom_boxplot(aes(fill=Stimulation),outlier.shape=NA)
p2 <- p + facet_grid(.~cells)+geom_jitter(shape=16, width = 0.2) +theme_classic()+ theme(axis.text.x=element_text(angle = 60,vjust = 1,hjust = 1))+
  scale_fill_manual(values=c("BG"="#67AD57","MDP"="#4A7CB3","UA"="#D1382C","oxLDL"="#E2E057"))+ylab("cell proportion")
  dat.sig$maxs <- 0.8
p2 + geom_signif(data=dat.sig, aes(xmin=st,xmax=ed, annotations=sig, y_position=maxs), 
+  manual = TRUE, vjust=0,tip_length = 0.01)  
  
  