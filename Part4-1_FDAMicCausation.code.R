library(presto)
library(dplyr)
library(tibble)
library(ggrepel)
library(fgsea)

setwd("/Projects/deng/Alzheimer/syn18485175")
ADMic=readRDS("ADMic.rds")
Idents(ADMic)=ADMic$MicCluster
ADMicTmp=subset(ADMic,ident=c(1:10))
new.cluster.ids <- c("BAMic", "HomMic", "HomMic", "HomMic", "HomMic", "HomMic", "NorARMic", "ARMic","FDAMic","DysMic")
names(new.cluster.ids) <- levels(ADMicTmp)
ADMicTmp <- RenameIdents(ADMicTmp, new.cluster.ids)
ADMicTmp$Gender=factor(ADMicTmp$Gender,levels=c("Male","Female"))
pdf("Manuscript/Microglia/Graph/Part4/ProlieferationGeneInClusterAddGender.pdf",width=7,height=10)
VlnPlot(subset(ADMicTmp),features=c("CD81","CD74","TYROBP","CSF1R","CORO1A","RPS3","VSIG4","HLA-DPB1","HLA-A","HLA-E","HLA-DPA1","CEBPB"),pt.size=0,cols = c("RoyalBlue","Violet"),split.by="Gender",ncol=1)&theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_blank(),plot.title = element_blank())
dev.off()

gene=c("RPS19","RPLP1","RPL13","FTL","SPP1","RPL28","ACTB","APOE","AKAP13","GAB2","ANKRD17","KANSL1","RAPGEF1","ZNF609","PIK3R5")
pdf("Manuscript/Microglia/Graph/Part4/TmpGeneInCluster.pdf",width=7,height=10)
scCustomize::Stacked_VlnPlot(seurat_object = ADMicTmp, features =gene,split.by="Gender",x_lab_rotate = TRUE, plot_spacing = 0.05,colors_use = c("RoyalBlue","Violet"))&
theme(axis.title.x=element_blank())
dev.off()

ARMic=subset(ADMic,MicCluster==9)
ARMic$Group=paste0(ARMic$MicCluster,ARMic$Gender)
ClusterMarker <- wilcoxauc(ARMic, 'Group')
ClusterMarker=ClusterMarker[ClusterMarker$group=="9Female",]
hs_data=ClusterMarker
hs_data$threshold = as.factor(ifelse(hs_data$logFC<0,ifelse(hs_data$pval<0.01&abs(hs_data$logFC)>0.25,"SigDown","Down"),ifelse(hs_data$pval<0.01&abs(hs_data$logFC)>0.25,"SigUp","Up")))
table(hs_data$threshold)
hs_data=hs_data[hs_data$threshold%in%c("SigUp","SigDown"),]
write.table(hs_data,file="Manuscript/Microglia/Graph/Part4/GenderDEGInC9.txt",sep="\t",quote=F)

Proliferation=c("CD81","VSIG4","CD74","RPS3","CSF1R","HLA-A","CEBPB","HLA-E","BST2","HLA-DPA1","CORO1A","HLA-DPB1","TYROBP")
hs_data=ClusterMarker
hs_data$threshold=ifelse(hs_data$feature%in%Proliferation,"Yes","No")
table(hs_data$threshold)

hs_data$ID=hs_data$feature
t=ggplot(data = hs_data, aes(x = logFC, y = -log10(pval), colour=threshold, label =ID, size=threshold )) +
  geom_point(alpha=0.6) +
  theme_bw() + scale_size_manual(values=c(0.3,3))+
  scale_color_manual(values=c("lightgrey", "red")) +
  geom_vline(xintercept=c(-0.25,0.25),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)",title="Differential Expressed Gene") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank()) +  
    geom_text_repel(
    data = subset(hs_data, hs_data$feature%in%Proliferation),
    aes(label = ID),
    size = 3,
    max.overlaps=20,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
tiff("Manuscript/Microglia/Graph/Part4/DEGInfor4ProgeneInARMicBetGender.tiff",height=400,width=500)
print(t)
dev.off()

t=ggplot(data = hs_data, aes(x = logFC, y = -log10(pval), colour=threshold, label =ID, size=threshold )) +
  geom_point(alpha=0.6) +
  theme_bw() + scale_size_manual(values=c(0.3,3))+
  scale_color_manual(values=c("lightgrey", "red")) +
  geom_vline(xintercept=c(-0.25,0.25),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="",y="",title="") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none", legend.title = element_blank(),axis.text.x =element_blank(),axis.text.y =element_blank()) +  
  geom_text_repel(
    data = subset(hs_data, hs_data$feature%in%Proliferation),
    aes(label = ""),
    size = 3,
    max.overlaps=20,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
tiff("Manuscript/Microglia/Graph/Part4/DEGInfor4ProgeneInARMicBetGenderNoLabel.tiff",height=400,width=450)
print(t)
dev.off()


#D:\Alzheimer\syn18485175\cluster13\monocle\SCENIC\ESR1Netowrk\permutation.pl
data=read.table("D:/Alzheimer/syn18485175/cluster13/monocle/SCENIC/ESR1Netowrk/AllTFRandNumber.txt",header=T)
ggplot(data,aes(x=Number)) + geom_density(alpha=.5,fill = "Plum1")+xlim(0,40)+theme_bw()


t=DotPlot(ADMicTmp,features=c("ESR1","ESR2"))
g=ggplot(t$data, aes(factor(features.plot,levels=rev(unique(features.plot))), id, size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_color_gradient(low="blue",high = "red")+
labs(color="Expression",size="Percent",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),angle=90),axis.text.y = element_text(size=rel(1.0))) +coord_flip()
pdf("Manuscript/Microglia/Graph/Part4/ExpressionESR1.pdf",height=2,width=4)
print(g)
dev.off()


t=DotPlot(ADMicTmp,features=c("ESR1","ESR2","AR","PGR","GPER1","INSR","CYP19A1"))
g=ggplot(t$data, aes(factor(features.plot,levels=rev(unique(features.plot))), id, size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_color_gradient(low="blue",high = "red")+
labs(color="Expression",size="Percent",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),angle=90),axis.text.y = element_text(size=rel(1.0))) +coord_flip()
pdf("Manuscript/Microglia/Graph/Part4/ExpressionHormoneReceptor.pdf",height=3,width=4)
print(g)
dev.off()

t=DotPlot(ADMicTmp,features=c("ESR1","ESR2","AR","PGR","GPER1","INSR","CYP19A1","RELA","RELB","REL"))
g=ggplot(t$data, aes(factor(features.plot,levels=rev(unique(features.plot))), id, size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_color_gradient(low="blue",high = "red")+
labs(color="Expression",size="Percent",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),angle=90),axis.text.y = element_text(size=rel(1.0))) +coord_flip()
pdf("Manuscript/Microglia/Graph/Part4/ExpressionHormoneReceptorNFkB.pdf",height=3,width=4)
print(g)
dev.off()

### AUC score for each TF across the subtypes ############
data=read.table("D:/Alzheimer/syn18485175/cluster13/monocle/SCENIC/regulonActivity_byCellType.txt",header=TRUE,row.names=1,sep="\t",check.names=FALSE)
meta=read.table("D:/Alzheimer/syn18485175/cluster13/monocle/SCENIC/ESR1Netowrk/ESR1PhysicalInteractCombine.txt",header=TRUE,row.names=1)
data1=data[matrixStats::rowMaxs(as.matrix(data[,2:dim(data)[2]]))>0.05,]
gene=intersect(rownames(data1),rownames(meta))
tmp=setdiff(rownames(data1),rownames(meta))
tmp=cbind(tmp,"ESR1"="NonInteract","ESR2"="NonInteract")
rownames(tmp)=tmp[,1]
metaAll=data.frame(rbind(meta[gene,],tmp[,-1]))
data1=data1[,-1]
ann_colors = list(ESR1=c(Interact="CadetBlue1",NonInteract="WhiteSmoke"), ESR2=c(Interact ="CadetBlue1",NonInteract = "WhiteSmoke"))
pheatmap(data1,annotation_row=metaAll,annotation_colors = ann_colors,clustering_method="ward.D2",scale="row",colorRampPalette(c("white","white","white","Wheat1","Salmon1","Firebrick3"))(50),cluster_col=FALSE,show_rownames=TRUE)

data1=data[matrixStats::rowMaxs(as.matrix(data[,2:dim(data)[2]]))>0.05,]
rownames(data1)=data1$Name
data1=data1[,-1]
pheatmap(data1,clustering_method="ward.D2",scale="row",colorRampPalette(c("white","white","white","Wheat1","Salmon1","Firebrick2"))(50),cluster_col=FALSE,show_rownames=TRUE)


### expression of the target genes for each TFs across the subtypes ############
setwd("/Projects/deng/Alzheimer/syn18485175")
ADMic=readRDS("ADMic.rds")
ADMic=subset(ADMic,SubType%in%c("BAMic", "HomMic","NorARMic", "ARMic","FDAMic","DysMic"))

setwd("/Projects/deng/Alzheimer/syn18485175/cluster13/monocle/SCENIC/")
data=read.table("/Projects/deng/Alzheimer/syn18485175/cluster13/monocle/SCENIC/regulonActivity_byCellType.txt",header=TRUE,row.names=1,sep="\t",check.names=FALSE)
data1=data[matrixStats::rowMaxs(as.matrix(data[,2:dim(data)[2]]))>0.05,]
rownames(data1)=data1$Name
data1=data1[,-1]
t=pheatmap(data1,clustering_method="ward.D2",scale="row",colorRampPalette(c("white","white","white","Wheat1","Salmon1","Firebrick2"))(50),cluster_col=FALSE,show_rownames=TRUE)

data2=data[matrixStats::rowMaxs(as.matrix(data[,2:dim(data)[2]]))>0.05,]
scenicOptions=readRDS("/Projects/deng/Alzheimer/syn18485175/cluster13/monocle/SCENIC/int/scenicOptions.Rds") 
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
TFTargetHighCon=regulonTargetsInfo[regulonTargetsInfo$TF%in%rownames(data2)&regulonTargetsInfo$highConfAnnot=="TRUE",]
extendTF=setdiff(rownames(data2),TFTargetHighCon$TF) #no high confidence target
TFTargetHighConLowTarget=unique(TFTargetHighCon$TF)[table(TFTargetHighCon$TF)<10] #smaller number of high confidence target
TFTargetLowCon=regulonTargetsInfo[regulonTargetsInfo$TF%in%c(extendTF,TFTargetHighConLowTarget),]
table(TFTargetLowCon$TF)
TFTarget=unique(rbind(TFTargetHighCon,TFTargetLowCon))
table(TFTarget$TF)
TFTargetList=split(TFTarget$gene,TFTarget$TF)


ClusterMarker <- wilcoxauc(ADMic, 'SubType')
table(ClusterMarker$group)
for(cluster in unique(ADMic$SubType)){
print (cluster)
clusterCell<- ClusterMarker %>% dplyr::filter(group == cluster) %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
ranks<- deframe(clusterCell)
fgseaRes<- fgseaMultilevel(TFTargetList, stats = ranks,eps=0)
fwrite(fgseaRes, file=paste0("/Projects/deng/Alzheimer/syn18485175/Manuscript/Microglia/Graph/Part4/GSEA/GSEA4TFTarget/",cluster,".txt",sep=""), sep="\t", sep2=c("", " ", ""))
}

TFTarget=read.table("D:/Alzheimer/syn18485175/Manuscript/Microglia/Graph/Part4/GSEA/GSEA4TFTarget.txt",header=T,sep="\t")
TFTarget=TFTarget[,c("Cluster","pathway","pval","NES","size")]
data=read.table("D:/Alzheimer/syn18485175/cluster13/monocle/SCENIC/regulonActivity_byCellType.txt",header=TRUE,row.names=1,sep="\t",check.names=FALSE)
data1=data[matrixStats::rowMaxs(as.matrix(data[,2:dim(data)[2]]))>0.05,]
data1=data1[,-1]
t=pheatmap(data1,clustering_method="ward.D2",scale="row",colorRampPalette(c("white","white","white","Wheat1","Salmon1","Firebrick2"))(50),cluster_col=FALSE,show_rownames=TRUE)
treeCol=t$tree_row
geneOrder=data.frame(TFOrder=treeCol$labels[treeCol$order])
tmp=merge(geneOrder,TFTarget,by.x="TFOrder",by.y="pathway")
tmp$sig=ifelse(tmp$pval<0.01,"Sig","NotSig")
TFOrder=factor(tmp$TFOrder,levels=geneOrder$TFOrder)
subytpeOrder=factor(tmp$Cluster,levels=c("BAMic","HomMic","NorARMic","ARMic","FDAMic","DysMic"))
t=ggplot(tmp,aes(subytpeOrder,factor(TFOrder,levels=rev(geneOrder$TFOrder)),size=size,colour=NES,shape=sig))+geom_point()+
scale_color_gradient2(low="OliveDrab",mid="white",high = "red")+
scale_size(range = c(1, 4))+
scale_shape_manual(values=c(3,19))+
theme_bw()+#theme(legend.position="bottom")+
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=rel(1.0),colour = "black",angle=90, vjust = 0.5, hjust=1),axis.text.y = element_text(size=rel(1.0),colour = "black"))
pdf("D:/Alzheimer/syn18485175/Manuscript/Microglia/Graph/Part4/GSEA4TFTarget.pdf",height=10,width=3)
print(t)
dev.off()


#### expression of the activate TFs across the subytpes ######
ADMicTmp$MicCluster=factor(ADMicTmp$MicCluster,levels=c(4,1,2,3,8,6,5,9,7,10))
activeTFs=c("BPTF","SP4","ELK3","ZEB1","MEF2C","JAZF1","MEF2A","FOXN3","FOXP2","TFEC","USF2","NR3C1","IKZF1","TBL1XR1","IRF2","IRF8","RCOR1","ZFX","HMGN3","CTCF","SF1","SMARCA4","TAF1","PHF8","RUNX1","ETV6","CEBPZ","HDAC2","KDM5A","TCF12","ELF2","FOXP1","FOXO3","THRB","JUN","MAFB","KLF6","RBPJ","CREM","SREBF1","TAF7","BHLHE41","YY1","NRF1","MAX","BACH1","KDM5B","STAT3","BCL6","BDP1","CEBPB","SPI1","YBX1","POLR2A","ESRRA")
ESR1TargetActiveTFs=c("SREBF1","STAT3","KDM5B","MAX","BCL6","SMARCA4","TBL1XR1","TCF12","IRF8","HMGN3","NR3C1","MEF2C","MEF2A","JAZF1","FOXP1","KLF6","RBPJ","THRB")

t=DotPlot(ADMicTmp,features=ESR1TargetActiveTFs,group.by="MicCluster")
g=ggplot(t$data, aes(factor(features.plot,levels=rev(unique(features.plot))), id, size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_color_gradient2(low="blue",mid="white",high = "red")+
labs(color="Expression",size="Percent",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0))) +coord_flip()
pdf("Manuscript/Microglia/Graph/Part4/ExpressionESR1TargetActiveTFs.pdf",height=4,width=4)
print(g)
dev.off()

ADMicExpr <- AverageExpression(ADMicTmp,group.by="MicCluster")[["RNA"]]
activeTFsExpr=ADMicExpr[ESR1TargetActiveTFs,]
pdf("Manuscript/Microglia/Graph/Part4/ExpressionESR1TargetActiveTFs.heatmap.pdf",width=3,height=4)
pheatmap(activeTFsExpr,clustering_method="ward.D2",border=NA,scale="row",cluster_col=FALSE,,cluster_row=FALSE,show_rownames=T,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()


colors_list <- c("forestgreen","NavajoWhite","LightBlue", "SlateBlue", "Violet", "gold")
pdf("Manuscript/Microglia/Graph/Part4/ExpressionESR1TargetActiveTFs.VlnPlot.pdf",width=6,height=10)
scCustomize::Stacked_VlnPlot(seurat_object =ADMicTmp,features = ESR1TargetActiveTFs, group.by="MicCluster",ggplot_default_colors = TRUE,x_lab_rotate = TRUE, plot_spacing = 0.05)&
theme(axis.title.x=element_blank())
dev.off()

#### co-expression for SPI1 and their target proliferation associated genes ########
setwd("/Projects/deng/Alzheimer/syn18485175/cluster13/monocle/SCENIC/")
scenicOptions=readRDS("/Projects/deng/Alzheimer/syn18485175/cluster13/monocle/SCENIC/int/scenicOptions.Rds") 
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
SPI1TargetHighCon=regulonTargetsInfo[regulonTargetsInfo$TF%in%"SPI1"&regulonTargetsInfo$highConfAnnot=="TRUE",]
ProliferationGenes=c("CD81","VSIG4","CD74","RPS3","CSF1R","HLA-A","CEBPB","HLA-E","BST2","HLA-DPA1","CORO1A","HLA-DPB1","TYROBP")
ProliferationGenesTargetBySPI1=intersect(ProliferationGenes,SPI1TargetHighCon$gene)

ProliferationGenes=c("CD74","HLA-DRB5","CD81","HLA-DRA","HLA-DRB1") #not Proliferation assocaited genes but labled
ProliferationGenesTargetBySPI1=ProliferationGenes

GeneExprFromROSMAP=read.table("/data2/deng/Alzheimer/Microglia/Gene_Syn3388564/geneExpr.txt",header=TRUE,row.names=1,sep="\t",check.names=F)
sampleInfo=read.table("/data2/deng/Alzheimer/Microglia/Gene_Syn3388564/Phenotype.txt",header=TRUE,row.names=1,sep="\t",check.names=F)
sampleList=intersect(rownames(sampleInfo),colnames(GeneExprFromROSMAP))
ProliferationGeneExprFromROSMAP=GeneExprFromROSMAP[c("SPI1",ProliferationGenesTargetBySPI1),sampleList]
ProliferationGeneExprFromROSMAP=log2(ProliferationGeneExprFromROSMAP+1)
sampleInfo=sampleInfo[sampleList,c("msex","apoe_genotype","age_at_visit_max","braaksc","ceradsc","ceradscID","cogdx","dcfdx_lv")]
all(rownames(sampleInfo)==colnames(ProliferationGeneExprFromROSMAP)) #TRUE
ProliferationGeneInfo=cbind(sampleInfo,t(ProliferationGeneExprFromROSMAP))
ProliferationGeneInfo.df=reshape2::melt(ProliferationGeneInfo,id=c(1:9))
colnames(ProliferationGeneInfo.df)=c(colnames(sampleInfo),"SPI1","Gene","Expr")
t=ggplot(ProliferationGeneInfo.df, aes(SPI1, Expr,color=msex),size=1)+
geom_point()+
scale_color_manual(values=c("Violet","RoyalBlue")) +
theme_bw()+
scale_y_continuous(position = "right")+
ggpmisc::stat_poly_line()+
ggpmisc::stat_poly_eq()+
facet_wrap(Gene~.,scale="free_y",strip.position="left",ncol=1)
pdf("/Projects/deng/Alzheimer/syn18485175/Manuscript/Microglia/Graph/Part4/Co-expressionTmp.pdf",height=10,width=4)
print(t)
dev.off()

sp <- ggpubr::ggscatter(ProliferationGeneInfo.df, x = "SPI1", y = "Expr",
  size = 1,
  color="msex",
  add = "reg.line",
  conf.int = TRUE,
  add.params = list(color = "black",fill = "lightgray"))+theme_bw()+
  scale_color_manual(values=c("Violet","RoyalBlue")) +
  facet_wrap(Gene~.,scale="free_y",strip.position="left",ncol=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggpubr::stat_cor(method = "pearson")
pdf("/Projects/deng/Alzheimer/syn18485175/Manuscript/Microglia/Graph/Part4/Co-expressionByggscatterTmp.pdf",height=8,width=4)
print(sp)
dev.off()


setwd("/Projects/deng/Alzheimer/syn18485175/")
ProliferationGeneInfo=cbind(sampleInfo,t(ProliferationGeneExprFromROSMAP))
sampleGroupInfo=read.table("Manuscript/Microglia/Graph/Part1/anno.txt",header=T)
all(rownames(sampleGroupInfo)==rownames(ProliferationGeneInfo)) #TRUE
data.df=cbind(sampleGroupInfo,ProliferationGeneInfo)
data.df$State=ifelse(data.df$Group=="G2","Ctrl","AD")
data.df$State=factor(data.df$State,levels=c("Ctrl","AD"))
data.df$Gender=factor(data.df$Gender,levels=c("Male","Female"))
t=ggplot(data.df, aes(x=State, y=CD74,color=Gender))+
geom_boxplot()+labs(title="",x="", y = "")+ #ylim(c(1.5,3.5))+
theme (strip.text= element_text(size=15)) + 
scale_color_manual(values=c("RoyalBlue","Violet")) +
theme_bw()+
theme(legend.position="right")+
geom_point(position=position_jitterdodge())+
theme(axis.title.y = element_text(size=rel(1.5)),axis.text.x = element_text(size=rel(1.5),vjust=0.5,hjust=1),axis.text.y = element_text(size=rel(1.0)))
pdf("Manuscript/Microglia/Graph/Part4/CD74ExpreInROSMAP.pdf",width=5,height=4)
print(t)
dev.off()

wilcox.test(data.df[data.df$State=="AD"&data.df$Gender=="Female","SPI1"],data.df[data.df$State=="Ctrl"&data.df$Gender=="Female","SPI1"]) #0.004595
wilcox.test(data.df[data.df$State=="AD"&data.df$Gender=="Female","SPI1"],data.df[data.df$State=="AD"&data.df$Gender=="Male","SPI1"]) #p-value = 0.0235
wilcox.test(data.df[data.df$State=="AD"&data.df$Gender=="Male","SPI1"],data.df[data.df$State=="AD"&data.df$Gender=="Male","SPI1"]) #p-value = 1
wilcox.test(data.df[data.df$State=="Ctrl"&data.df$Gender=="Male","SPI1"],data.df[data.df$State=="Ctrl"&data.df$Gender=="Female","SPI1"]) #p-value = 0.4498

wilcox.test(data.df[data.df$State=="AD"&data.df$Gender=="Female","CD81"],data.df[data.df$State=="Ctrl"&data.df$Gender=="Female","CD81"])$p.value #0.4318857
wilcox.test(data.df[data.df$State=="AD"&data.df$Gender=="Female","CD81"],data.df[data.df$State=="AD"&data.df$Gender=="Male","CD81"])$p.value #0.001806208
wilcox.test(data.df[data.df$State=="AD"&data.df$Gender=="Male","CD81"],data.df[data.df$State=="AD"&data.df$Gender=="Male","CD81"])$p.value  #1
wilcox.test(data.df[data.df$State=="Ctrl"&data.df$Gender=="Male","CD81"],data.df[data.df$State=="Ctrl"&data.df$Gender=="Female","CD81"])$p.value #0.4495116

wilcox.test(data.df[data.df$State=="AD"&data.df$Gender=="Female","CD74"],data.df[data.df$State=="Ctrl"&data.df$Gender=="Female","CD74"])$p.value #0.009800488
wilcox.test(data.df[data.df$State=="AD"&data.df$Gender=="Female","CD74"],data.df[data.df$State=="AD"&data.df$Gender=="Male","CD74"])$p.value #0.650136
wilcox.test(data.df[data.df$State=="AD"&data.df$Gender=="Male","CD74"],data.df[data.df$State=="AD"&data.df$Gender=="Male","CD74"])$p.value #1
wilcox.test(data.df[data.df$State=="Ctrl"&data.df$Gender=="Male","CD74"],data.df[data.df$State=="Ctrl"&data.df$Gender=="Female","CD74"])$p.value #0.2477014

wilcox.test(data.df[data.df$State=="AD"&data.df$Gender=="Female","TYROBP"],data.df[data.df$State=="Ctrl"&data.df$Gender=="Female","TYROBP"])$p.value #0.002056596
wilcox.test(data.df[data.df$State=="AD"&data.df$Gender=="Female","TYROBP"],data.df[data.df$State=="AD"&data.df$Gender=="Male","TYROBP"])$p.value #0.04000875
wilcox.test(data.df[data.df$State=="AD"&data.df$Gender=="Male","TYROBP"],data.df[data.df$State=="AD"&data.df$Gender=="Male","TYROBP"])$p.value #1
wilcox.test(data.df[data.df$State=="Ctrl"&data.df$Gender=="Male","TYROBP"],data.df[data.df$State=="Ctrl"&data.df$Gender=="Female","TYROBP"])$p.value #0.5654708


Cluster7DEG=read.table("D:/Alzheimer/syn18485175/cluster13/monocle/ADMicSubCluster7Marker.txt",header=T,row.names=1)
Cluster9DEG=read.table("D:/Alzheimer/syn18485175/cluster13/monocle/ADMicSubCluster9Marker.txt",header=T,row.names=1)
Cluster10DEG=read.table("D:/Alzheimer/syn18485175/cluster13/monocle/ADMicSubCluster9Marker.txt",header=T,row.names=1)
Cluster7910DEG=unique(c(rownames(Cluster7DEG),rownames(Cluster9DEG),rownames(Cluster10DEG)))
tableSubset <- regulonTargetsInfo[TF=="SMARCA4"&highConfAnnot=="TRUE"]
#gene=intersect(rownames(Cluster7DEG),tableSubset$gene)
DoHeatmap(ADMic,features=tableSubset$gene)+ theme(axis.text.y = element_text(size = 0))
viewMotifs(tableSubset, options=list(pageLength=100)) 

tableSubset <- regulonTargetsInfo[TF=="SMARCA4"]
DoHeatmap(ADMic,features=tableSubset$gene)+ theme(axis.text.y = element_text(size = 0))

t=enrichGO(gene,org.Hs.eg.db, keyType = "SYMBOL",ont = "BP",pvalueCutoff = 0.01)
t1 <- clusterProfiler::simplify(t,cutoff = 0.5,select_fun = min,measure = "Wang")
dotplot(t1,showCategory=15)

for(n in rownames(data)){
  tableSubset <- regulonTargetsInfo[TF==n ]
  gene=tableSubset$gene
  result=t(c(n,gene))
  write.table(result, file="test.txt", sep="\t", append=TRUE, row.names=FALSE, col.names=FALSE,quote=F)
}

library(clusterProfiler)
DEG97=FindMarkers(ADMic,ident.1=9,ident.2=7,pct.min=0.25)
DEG=DEG97
DEGTmp=DEG[DEG$p_val<0.01,]
dim(DEGTmp)
up=DEGTmp[DEGTmp$avg_logFC> 0,]
down=DEGTmp[DEGTmp$avg_logFC< 0,]
DegFunction <- enrichGO(rownames(DEGTmp), 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.01, keyType= "SYMBOL")
DegFunction1 <- simplify(DegFunction,cutoff = 0.5,select_fun = min,measure = "Wang")
dotplot(DegFunction1,showCategory=15)

data=read.table("D:/Alzheimer/syn18485175/cluster13/monocle/SCENIC/TFTarget/C7ProliferHiCTFRelated.txt",header=T)
pset=array()
for(i in 1:dim(data)[1]){pset[i]=dhyper((data[i,6]-1),data[i,4],(data[i,3]-data[i,4]),data[i,5])[1]} #followd by manual integrated
as.matrix(pset)
data=read.table("D:/Alzheimer/syn18485175/cluster13/monocle/SCENIC/TFTarget/ActiveTFSetInfoPvalue8dhyper.txt",header=T)
TF=factor(data$Gene,levels=rev(data$Gene))
ggplot(data=data, aes(x=TF, y=-log10(Pvalue4C7ProliferHighTF),fill=Set)) +
  geom_bar(stat="identity")+
  coord_flip()+
  theme_minimal()


#effect of menopouse (deficiency of estrogen) on FDAMic signature#
#effect of SPI1 on FDAMic signature#
setwd("/Projects/deng/Alzheimer/syn18485175")
tiff("Manuscript/Microglia/Graph/Part4/SPI1cellType.tiff",height=400,width=400)
FeaturePlot(AD,c("SPI1"),cols =c("lightgrey","red"))&NoLegend()&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

ARMic=subset(ADMic,idents=9)
Idents(ARMic)=factor(ARMic$Gender,levels=c("Male","Female"))
ARMic$Statues=factor(ARMic$Statues,levels=c("Control","Alzheimer"))
pdf("Manuscript/Microglia/Graph/Part4/SPI1InGender.pdf",height=3,width=6)
VlnPlot(ARMic,features=c("SPI1"),pt.size=0.1,cols = c("Blue","Violet"),split.by="Statues",ncol=1)&theme(axis.title.x=element_blank(),axis.title.y=element_blank(),plot.title = element_blank())
dev.off()


###effect of methylation on SPI1 ######
D:\Alzheimer\Microglia\Mathylation_syn3157275\CombineGeneExpr.code.R

###effect of methylation on FDAMic ######
D:\Alzheimer\Microglia\Mathylation_syn3157275\Mathylation_syn3157275.code.R
D:\Alzheimer\Microglia\Mathylation_syn3157275\DML\DML8Pval.code.R

### Similir methylation between female AD and estrogen deficiency ######
## Extrogen associated methylation code D:\Alzheimer\Microglia\Sklias_Estrogen_GSE132513\Sklias_Estrogen_GSE132513.code.R
Estrogen_E2D_Ctrl=read.csv("/data2/deng/Alzheimer/Microglia/Sklias_Estrogen_GSE132513/DML8limma_E2D-CTR.sigLoci.txt",header=T,row.names=1,sep="\t")
Estrogen_E2D_CtrlGene=unique(strsplit(paste0(unique(Estrogen_E2D_Ctrl$Gene),collapse = ";"),";")[[1]])
length(Estrogen_E2D_CtrlGene)
#8812
HyperInEstrogen_E2D=Estrogen_E2D_Ctrl[Estrogen_E2D_Ctrl$logFC>0,]
HypoInEstrogen_E2D=Estrogen_E2D_Ctrl[Estrogen_E2D_Ctrl$logFC<0,]
HyperGeneInEstrogen_E2D=unique(strsplit(paste0(unique(HyperInEstrogen_E2D$Gene),collapse = ";"),";")[[1]])
HypoGeneInEstrogen_E2D=unique(strsplit(paste0(unique(HypoInEstrogen_E2D$Gene),collapse = ";"),";")[[1]])
length(HyperGeneInEstrogen_E2D)
#7364
length(HypoGeneInEstrogen_E2D)
#2240
write.table(HyperGeneInEstrogen_E2D,file="DMLRelatedGene/HyperGeneInEstrogen_E2D.txt",quote=F,col.names=F,row.names=F)
write.table(HypoGeneInEstrogen_E2D,file="DMLRelatedGene/HypoGeneInEstrogen_E2D.txt",quote=F,col.names=F,row.names=F)


Estrogen_ReSt_E2D=read.csv("/data2/deng/Alzheimer/Microglia/Sklias_Estrogen_GSE132513/DML8limma_ReSt-E2D.sigLoci.txt",header=T,row.names=1,sep="\t")
Estrogen_ReSt_E2DGene=unique(strsplit(paste0(unique(Estrogen_ReSt_E2D$Gene),collapse = ";"),";")[[1]])
length(Estrogen_ReSt_E2DGene)
#4263
HyperInEstrogen_ReSt_E2D=Estrogen_ReSt_E2D[Estrogen_ReSt_E2D$logFC>0,]
HypoInEstrogen_ReSt_E2D=Estrogen_ReSt_E2D[Estrogen_ReSt_E2D$logFC<0,]
HyperGeneInEstrogen_ReSt_E2D=unique(strsplit(paste0(unique(HyperInEstrogen_ReSt_E2D$Gene),collapse = ";"),";")[[1]])
HypoGeneInEstrogen_ReSt_E2D=unique(strsplit(paste0(unique(HypoInEstrogen_ReSt_E2D$Gene),collapse = ";"),";")[[1]])
length(HyperGeneInEstrogen_ReSt_E2D)
#1074
length(HypoGeneInEstrogen_ReSt_E2D)
#3442
write.table(HyperGeneInEstrogen_ReSt_E2D,file="DMLRelatedGene/HyperGeneInEstrogen_ReSt_E2D.txt",quote=F,col.names=F,row.names=F)
write.table(HypoGeneInEstrogen_ReSt_E2D,file="DMLRelatedGene/HypoGeneInEstrogen_ReSt_E2D.txt",quote=F,col.names=F,row.names=F)

Estrogen_E2D_Ctrl_List=split(rownames(Estrogen_E2D_Ctrl),Estrogen_E2D_Ctrl$threshold)
names(Estrogen_E2D_Ctrl_List)=c("HyperInE2D","HypoInE2D")
Estrogen_ReSt_E2D_List=split(rownames(Estrogen_ReSt_E2D),Estrogen_ReSt_E2D$threshold)
names(Estrogen_ReSt_E2D_List)=c("HyperInReSt","HypoInReSt")
Estrogen_list=c(Estrogen_E2D_Ctrl_List,Estrogen_ReSt_E2D_List)
t=upset(fromList(Estrogen_list),nsets = 4,
  sets=rev(c("HyperInE2D","HypoInE2D","HyperInReSt","HypoInReSt")),
  keep.order=TRUE,
  main.bar.color=c("Orange","SlateBlue","Orange","SlateBlue","Red","lightgrey","lightgrey","lightgrey"),
  matrix.color=c("black"),
  sets.bar.color=rev(c("Orange","SlateBlue","Orange","SlateBlue"))
  )
pdf("Manuscript/Microglia/Graph/Part4/Estrogen_methylated_loci.pdf")
print(t)
dev.off()

cgMeta=read.table("/data2/deng/Alzheimer/Microglia/Sklias_Estrogen_GSE132513/cgMeta.txt",sep="\t",row.names=1,header=T)
cgMetaInfo=cgMeta[,c("Name","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_CpG_Islands_Name","Relation_to_UCSC_CpG_Island")]
dim(cgMetaInfo)
#864935      6
targetCgId=intersect(Estrogen_list$HyperInE2D,Estrogen_list$HypoInReSt) #1799
length(Estrogen_list$HyperInE2D)
#17638
length(Estrogen_list$HypoInReSt)
#6125
dhyper(1799,17638,864935-17638,6125)  #estrogen deficiency and restimulation shared significant methylated genes
pvalue=0

targetCgIdInfo=cgMetaInfo[targetCgId,]
targetCgIdInfoGene=unique(strsplit(paste0(unique(targetCgIdInfo$UCSC_RefGene_Name),collapse = ";"),";")[[1]])
length(targetCgIdInfoGene) #defined as the estrogen associated genes
#1126
write.table(targetCgIdInfoGene,file="/data2/deng/Alzheimer/Microglia/Mathylation_syn3157275/DML/EstrogenAssociatedMethylatedGene.txt",quote=F,row.names=F,col.names=F)


FemaleADMethylation=read.csv("/data2/deng/Alzheimer/Microglia/Mathylation_syn3157275/DML/SigDML_Female_G2vsG1_Symbol.txt",header=T,row.names=1,sep="\t")
FemaleADMethylationGene=unique(strsplit(paste0(unique(FemaleADMethylation$RefGene),collapse = ";"),";")[[1]])
length(FemaleADMethylationGene)
#10384 not #7448
write.table(FemaleADMethylationGene,file="/data2/deng/Alzheimer/Microglia/Mathylation_syn3157275/DML/SigDML_Female_AD_Gene.txt",quote=F,row.names=F,col.names=F)

overlap=intersect(targetCgIdInfoGene,FemaleADMethylationGene)
length(overlap) #609 estrogen associated genes shared by female AD
#731 not 609
dhyper(609,1126,(25609-1126),7448) #3.02e-72
dhyper(731,1126,(25609-1126),10384) #5.061779e-64
dhyper(162,1125,(25609-1126),2399) #P_hypo=7.07e-09
dhyper(682,1125,(25609-1125),9032) #P_hyper=2.60e-70

overlap=overlap[-1] #remove the blank
overlapMappedInfo=targetCgIdInfo[grepl(paste(overlap, collapse="|"),targetCgIdInfo$UCSC_RefGene_Name),]
dim(overlapMappedInfo)
1799

#expression of the estrogen associated cg id across different conditions
EstrogenMethy=read.csv("/data2/deng/Alzheimer/Microglia/Sklias_Estrogen_GSE132513/betaValue.txt",header=T,sep="\t",row.names=1)
EstrogenMethySample=read.csv("/data2/deng/Alzheimer/Microglia/Sklias_Estrogen_GSE132513/sampleInfo.txt",header=T,sep="\t",row.names=1)
EstrogenCgIDBeta=EstrogenMethy[targetCgId,]
dim(EstrogenCgIDBeta) #1799 cgId==1126 gene
all(colnames(EstrogenCgIDBeta)==rownames(EstrogenMethySample))
EstrogenCgIDBetaInfo=cbind(EstrogenMethySample,t(EstrogenCgIDBeta))
EstrogenCgIDBetaInfo.df=reshape2::melt(EstrogenCgIDBetaInfo,id=c(1:4))
colnames(EstrogenCgIDBetaInfo.df)=c(colnames(EstrogenMethySample),"cgId","betaValue")
EstrogenCgIDBetaInfo.df$FemaleAD=ifelse(EstrogenCgIDBetaInfo.df$cgId%in%rownames(overlapMappedInfo),"FemaleADLoci","OtherLoci")

t=ggplot(EstrogenCgIDBetaInfo.df, aes(x=Title, y=betaValue))+
geom_boxplot(aes(fill=Title),alpha=0.4)+labs(title="",x="", y = "")+ #ylim(0,0.2)+
geom_jitter(size=0.001,alpha=0.2,aes(color=FemaleAD))+
geom_line(aes(group = cgId,color=FemaleAD),linewidth=0.001,alpha=0.1,linetype="dotdash")+
scale_color_manual(values=c("Tomato","lightblue"))+
scale_fill_manual(values=c(rep("#5AB4AC",9),rep("#966223",6),rep("#3D8882",3)))+
theme (strip.text= element_text(size=10,face="bold")) + 
theme(legend.position="none")+ theme_bw()+
theme(axis.title.y = element_text(size=rel(1.0)),axis.text.x = element_text(angle = 90,vjust=0.3,hjust=1),axis.text.y = element_text(size=rel(1.0)))
pdf("Manuscript/Microglia/Graph/Part4/EstrogenAssociatedCGID.pdf",height=5,width=6)
print(t)
dev.off()

EstrogenMethySample$Time=factor(EstrogenMethySample$Time,levels=c("d0","d4","d14"))
EstrogenMethySample=EstrogenMethySample[order(EstrogenMethySample$Group,EstrogenMethySample$Time),]
EstrogenCgIDBeta=EstrogenCgIDBeta[,rownames(EstrogenMethySample)]
cgGroup=data.frame(Type=ifelse(rownames(EstrogenCgIDBeta)%in%rownames(overlapMappedInfo),"FemaleADLoci","OtherLoci"))
rownames(cgGroup)=rownames(EstrogenCgIDBeta)
ann_colors = list(
    Time = c(d0="Beige",d4="Khaki",d14="Goldenrod"),
    Group = c(CTR = "#5AB4AC", E2D = "#966223",ReSt = "#3D8882"),
    Type = c(FemaleADLoci = "firebrick", OtherLoci = "lightgray")
)
pdf("Manuscript/Microglia/Graph/Part4/EstrogenAssociatedCGID_heatmap.pdf",height=6,width=5)
pheatmap(EstrogenCgIDBeta,cluster_col=F,annotation_col=EstrogenMethySample[,c("Group","Time")],annotation_row=cgGroup,
  annotation_colors = ann_colors,color = viridis(10,option = "A"),
  clustering_method="ward.D2",scale="row",show_rownames=F,show_colnames=F)
dev.off()



MethylateGeneInMic=read.table("/data2/deng/Alzheimer/Microglia/Mathylation_syn3157275/DML/CellTypeSpecificPval/MicGeneDMC.txt",header=T,row.names=1,sep="\t")
dim(MethylateGeneInMic) #827   2
overlap=intersect(rownames(MethylateGeneInMic),targetCgIdInfoGene)
length(overlap)
#70

#25609 as the background gene list from GPL21145  #/data2/deng/Alzheimer/Microglia/Sklias_Estrogen_GSE132513
dhyper(70-1,8812,(25609-8812),827)
p=1.31e-296

overlap=intersect(rownames(MethylateGeneInMic),FemaleADMethylationGene)

#----------DEG between AD and Ctrl in each cluster----------------------------
ADMic=readRDS("ADMic.rds")
Idents(ADMic)=ADMic$MicCluster
ADMicTmp=subset(ADMic,ident=c(1:10))
ident=c(4,1,2,3,8,6,5,9,7,10)
ADMicTmp@active.ident <- factor(ADMic@active.ident,levels=ident)
new.cluster.ids <- c("BAMic", "HomMic", "HomMic", "HomMic", "HomMic", "HomMic", "NorARMic", "ARMic","FDAMic","DysMic")
names(new.cluster.ids) <- levels(ADMicTmp)
ADMicTmp <- RenameIdents(ADMicTmp, new.cluster.ids)
ADMicTmp$cellType=Idents(ADMicTmp)

for(cluster in unique(ADMicTmp$cellType)){
  tmpCluster=subset(ADMicTmp,cellType%in%cluster)
  Idents(tmpCluster)=tmpCluster$Statues
  Female=subset(tmpCluster,Gender=="Female")
  DEGInFemale=FindMarkers(Female,ident.1="Alzheimer",ident.2="Control",min.pct = 0.25)
  DEGInFemale=DEGInFemale[DEGInFemale$p_val<0.01,]
  DEGInFemale$pattern=ifelse(DEGInFemale$avg_log2FC>0,"Up","Dn")
  write.table(DEGInFemale,file=paste0("Manuscript/Microglia/Graph/Part4/DEG/Female/",cluster,".txt",sep=""),quote=F,sep="\t")
  
  Male=subset(tmpCluster,Gender=="Male")
  DEGInMale=FindMarkers(Male,ident.1="Alzheimer",ident.2="Control",min.pct = 0.25)
  DEGInMale=DEGInMale[DEGInMale$p_val<0.01,]
  DEGInMale$pattern=ifelse(DEGInMale$avg_log2FC>0,"Up","Dn")
  write.table(DEGInMale,file=paste0("Manuscript/Microglia/Graph/Part4/DEG/Male/",cluster,".txt",sep=""),quote=F,sep="\t")
}


Idents(ARMic)=factor(ARMic$Statues,levels=c("Control","Alzheimer"))
FemaleARMic=subset(ARMic,Gender=="Female")
DEG=FindMarkers(FemaleARMic,ident.1="Alzheimer")
DEG=DEG[DEG$p_val<0.01,]
DEG$pattern=ifelse(DEG$avg_log2FC>0,"Up","Dn")
table(DEG$pattern)
up=DEG[DEG$avg_log2FC>0,]
down=DEG[DEG$avg_log2FC<0,]
upId=bitr(rownames(up),"SYMBOL","ENTREZID",OrgDb="org.Hs.eg.db")[,2]
downId=bitr(rownames(down),"SYMBOL","ENTREZID",OrgDb="org.Hs.eg.db")[,2]
DegId=list("Up"=upId,"Down"=downId)
DegFunction <- compareCluster(geneCluster = DegId, fun = "enrichGO",ont="BP",OrgDb='org.Hs.eg.db',minGSSize = 5,pvalueCutoff=0.01)
write.table(DEG,file="Manuscript/Microglia/Graph/Part4/DEGbetStatuesInFemale.txt",sep="\t",quote=F)
y <- setReadable(DegFunction, 'org.Hs.eg.db', 'ENTREZID')
write.table(y,file="Manuscript/Microglia/Graph/Part4/DEGGOBPbetStatuesInFemale.txt",sep="\t",quote=F)

pdf("Manuscript/Microglia/Graph/Part4/DEGGOBPbetStatuesInFemale.pdf",width=8,height=10)
dotplot(DegFunction,showCategory=25,label_format = 100)
dev.off()

pdf("Manuscript/Microglia/Graph/Part4/DEGbetStatuesInFemale.pdf",height=7,width=6)
VlnPlot(ARMic,features=c("HLA-DRB1","HIF1A"),pt.size=0.1,cols = c("Blue","Violet"),split.by="Statues",ncol=1)&theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),plot.title = element_blank())
dev.off()
#p_val avg_log2FC pct.1 pct.2 p_val_adj
#SPI1 0.1457659  0.8070626 0.359 0.244         1

FDAMic=read.table("/Projects/deng/Alzheimer/syn18485175/cluster13/monocle/ADMicSubCluster7Marker.txt",header=T)
FDAMicMarker=FDAMic[FDAMic$p_val<0.01,]
FDAMicMarker$Pattern=ifelse(FDAMicMarker$avg_logFC>0,"Up","Dn")
FDAMicMarker$Symbol=rownames(FDAMicMarker)
fgsea_sets=split(FDAMicMarker$Symbol,FDAMicMarker$Pattern)

Idents(ARMic)=factor(ARMic$Statues,levels=c("Control","Alzheimer"))
FemaleARMic=subset(ARMic,Gender=="Female")
FemaleARMicDEG <- wilcoxauc(FemaleARMic, 'Statues')
clusterCell<- FemaleARMicDEG %>% dplyr::filter(group == "Alzheimer") %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
ranks<- deframe(clusterCell)
fgseaRes<- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0)
pdf("Manuscript/Microglia/Graph/Part4/GSEAPlotOfFDAMarkerUpInARMicFemaleAD.pdf",height=3,width=4)
plotEnrichment(fgsea_sets[["Up"]],ranks) 
dev.off()
pdf("Manuscript/Microglia/Graph/Part4/GSEAPlotOfFDAMarkerDnInARMicFemaleAD.pdf",height=3,width=4)
plotEnrichment(fgsea_sets[["Dn"]],ranks) 
dev.off()
fwrite(fgseaRes, file="Manuscript/Microglia/Graph/Part4/GSEAPlotOfFDAMarkerDnInARMicFemaleAD.txt", sep="\t", sep2=c("", " ", ""))

