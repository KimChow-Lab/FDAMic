setwd("/Projects/deng/Alzheimer/syn18485175")
library(Seurat)
library(monocle3)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

Mic.cds=readRDS("ADMicMono.rds")
t=plot_cells(Mic.cds, color_cells_by="APOEGenotype",group_label_size = 5,cell_size =1,label_cell_groups = FALSE)+scale_color_gradient(low="gray90",high = "blue")+
theme(axis.title.x = element_blank(),axis.title.y = element_blank()) 

pData(Mic.cds)$APOEGroup=as.character(pData(Mic.cds)$APOEGenotype)
tiff("Manuscript/Microglia/Graph/Part5/APOEGenotype.tiff",width=500,height=400)
plot_cells(Mic.cds, color_cells_by="APOEGroup",group_label_size = 5,cell_size =1,label_cell_groups = FALSE)+
scale_color_brewer(palette="Spectral",direction=-1)+
theme(legend.position="none",axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())&
theme(panel.border = element_rect(fill=NA,color="black", linewidth=1.5, linetype="solid"))
dev.off()

tiff("Manuscript/Microglia/Graph/Part5/APOEGenotypeLabel.tiff",width=600,height=400)
plot_cells(Mic.cds, color_cells_by="APOEGroup",group_label_size = 5,cell_size =1,label_cell_groups = FALSE)+
scale_color_brewer(palette="Spectral",direction=-1)+
theme(panel.border = element_rect(fill=NA,color="black", linewidth=1.5, linetype="solid"))
dev.off()

ADMic=readRDS("ADMic.rds")
Idents(ADMic)=ADMic$MicCluster
ADMicTmp=subset(ADMic,ident=c(1:10))
ident=c(4,1,2,3,8,6,5,9,7,10)
ADMicTmp@active.ident <- factor(ADMic@active.ident,levels=ident)
new.cluster.ids <- c("BAMic", "HomMic", "HomMic", "HomMic", "HomMic", "HomMic", "NorARMic", "ARMic","FDAMic","DysMic")
names(new.cluster.ids) <- levels(ADMicTmp)
ADMicTmp <- RenameIdents(ADMicTmp, new.cluster.ids)

cellType=c("BAMic","HomMic", "NorARMic", "ARMic","FDAMic","DysMic")
tmp=paste0(Idents(ADMicTmp),"_",ADMicTmp$APOEGenotype,"_",sep="")
tmp=as.matrix(table(tmp))
TmpInfo=as.matrix(do.call(rbind, strsplit(as.character(rownames(tmp)),'_')))
result=data.frame("Cluster"=paste0(TmpInfo[,1],sep=""),"APOEGenotype"=TmpInfo[,2],"Number"=tmp[,1])
Cluster_order=factor(result$Cluster,levels=paste0(cellType,sep=""))
t=ggplot(result, aes(Cluster_order, Number, fill=APOEGenotype)) +
  geom_bar(stat="identity",position="fill") +
  scale_fill_brewer(palette="Spectral",direction=-1)+
  theme_bw()+
  scale_y_continuous(expand=c(0,0))+
  labs(size="",x="",y="Cell Number Percent",title="")+
  theme(axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1))
pdf("Manuscript/Microglia/Graph/Part5/APOEGenotypeInEachSubytpe.pdf",width=4,height=3)
print(t)
dev.off()
pg <- ggplot_build(t)
data=result
data$ScaledValue=pg$data[[1]]$y
write.table(data,file="Manuscript/Microglia/Graph/Part5/APOEGenotypeInEachSubytpe_ScaleData.txt",sep="\t",row.names=F,quote=F)


marker=c("APOE")
cellTypeCol=c("forestgreen","NavajoWhite", "LightBlue", "SlateBlue", "Violet", "gold")
FemaleADMic=subset(ADMicTmp,Gender=="Female")
FemaleADMicNor=subset(FemaleADMic,APOEGenotype %in% c("23","33","22"))
p1=VlnPlot(FemaleADMic,features=marker,pt.size=1,cols=cellTypeCol)&theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),plot.title = element_blank())
p2=VlnPlot(FemaleADMicNor,features=marker,pt.size=1,cols=cellTypeCol)&theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),plot.title = element_blank())
pdf("Manuscript/Microglia/Graph/Part5/APOEExpression.pdf",width=5,height=3)
print(p1+p2)
dev.off()

e1=data.frame(AverageExpression(FemaleADMic,assays = "RNA", slot = "counts"))
e1["APOE",]
UpperRow(AllCells) RNA.BAMic RNA.HomMic RNA.NorARMic RNA.ARMic RNA.FDAMic RNA.DysMic
APOE 0.3176471  0.6400778    0.7567568       3.8   2.108696   1.507246
e2=data.frame(AverageExpression(FemaleADMicNor,assays = "RNA", slot = "counts"))
e2["APOE",]
BottomRow(CellsRemovedE44)      RNA.BAMic RNA.HomMic RNA.NorARMic RNA.ARMic RNA.FDAMic RNA.DysMic
APOE 0.3863636  0.6414474    0.7567568  3.986301   1.757576       1.55


VlnPlot(ADMicTmp,features=marker,pt.size=0,ncol=1,split.by="APOEGenotype",cols=c("#A3A500","#00BF7D","#00B0F6","#E76BF3"))&theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),plot.title = element_blank())

MaleExpr8APOE3Group=read.table("D:/Alzheimer/Syn3388564/Male/MaleExpr8APOE3Group.txt",header=T,row.names=1,check.names=F)
FeMaleExpr8APOE3Group=read.table("D:/Alzheimer/Syn3388564/Female/FeMaleExpr8APOE3Group.txt",header=T,row.names=1,check.names=F)
Fdata=read.table("D:/Alzheimer/Syn3388564/Female/STEM8APOE3Group/mainTable.txt",header=T,check.names=F,row.names=1)
FPro=Fdata[Fdata$Profile %in% c("12","15","13","2","3","0"),] #12, 15, 13 increased, 2,3 0 decreased
FG=intersect(rownames(FeMaleExpr8APOE3Group),rownames(FPro))
length(FG) #3269
pdf("ContinuesGeneAcrossAPOEInFemale8StemHeatmap.pdf",height=6,width=4)
F=pheatmap(log2(FeMaleExpr8APOE3Group[FG,]+1),scale="row",cluster_col=FALSE,border=FALSE,show_rownames=FALSE,clustering_method="ward.D2",treeheight_row=30,colorRampPalette(c("navy", "white","pink","firebrick3"))(50))
dev.off()

Mdata=read.table("D:/Alzheimer/Syn3388564/Male/STEM8APOE3Group/mainTable.txt",header=TRUE,check.names=FALSE,row.names=1)
MPro=Mdata[Mdata$Profile %in% c("15","12","13","2","3","0"),]#12, 15, 13 increased, 2,3 0 decreased
MG=intersect(rownames(MaleExpr8APOE3Group),rownames(MPro))
length(MG) #1909
pdf("ContinuesGeneAcrossAPOEInMale8StemHeatmap.pdf",height=4,width=3)
M=pheatmap(log2(MaleExpr8APOE3Group[MG,]+1),scale="row",cluster_col=FALSE,border=FALSE,show_rownames=FALSE,clustering_method="ward.D2",treeheight_row=20,colorRampPalette(c("navy", "white","pink","firebrick3"))(50))
dev.off()

t1=cutree(F$tree_row,k=2)
pheatmap(log2(FeMaleExpr8APOE3Group[FG,]+1),scale="row",cluster_col=FALSE,border=FALSE,show_rownames=FALSE,clustering_method="ward.D2",treeheight_row=15,annotation_row=data.frame(t1))
FGene=names(t1[t1==2])
t2=cutree(M$tree_row,k=2)
pheatmap(log2(MaleExpr8APOE3Group[MG,]+1),scale="row",clustering_method="ward.D2",cluster_col=FALSE,border=FALSE,show_rownames=FALSE,treeheight_row=15,annotation_row=data.frame(t2))
MGene=names(t2[t2==2])
Deg=list("Female"=FGene,"Male"=MGene)
DegFunction <- compareCluster(geneCluster = Deg, fun = "enrichGO",ont="BP",keyType = "SYMBOL", OrgDb='org.Hs.eg.db',minGSSize = 5,pvalueCutoff=0.01)
DegFunction1 <- clusterProfiler::simplify(DegFunction,cutoff = 0.5,select_fun = min,measure = "Wang")
pdf("GOBP4CUp.pdf",width=6,height=4)
dotplot(DegFunction1,showCategory=10)
dev.off()

write.table(data.frame(DegFunction1),file="Female1Male3ModuleFunction.txt",sep="\t",quote=FALSE)
write.table(FGene,file="FemaleAPOEContinueGene.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(MGene,file="MaleAPOEContinueGene.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

FDAMicMarker=read.table("D:/Alzheimer/syn18485175/cluster13/monocle/ADMicSubCluster7Marker.txt",header=T,row.names=1)
FDAMicMarker=FDAMicMarker[FDAMicMarker$p_val<0.01,]
FG=read.table("D:/Alzheimer/syn18485175/Manuscript/Microglia/Graph/Part5/FemaleAPOEContinueGene.txt")
intersect(FG[,1],rownames(FDAMicMarker)) #183
dim(FDAMicMarker) #368
dhyper(183,368,14251-368,length(FG)) 


FDAMicMarkerList=rownames(FDAMicMarker[FDAMicMarker$avg_logFC>0,])
length(FDAMicMarkerList)

MG.df=data.frame(Symbol=MG)
MG.df$Gender="Male"
FG.df=data.frame(Symbol=FG)
FG.df$Gender="Female"
data.df=rbind(MG.df,FG.df)

MaleExpr=log2(MaleExpr8APOE3Group[MG,]+1)
FemaleExpr=log2(FeMaleExpr8APOE3Group[FG,]+1)
rownames(MaleExpr)=paste0("Male",rownames(MaleExpr8APOE3Group[MG,]))
rownames(FemaleExpr)=paste0("Female",rownames(FeMaleExpr8APOE3Group[FG,]))

MaleGene.df=data.frame(Name=rownames(MaleExpr),symbol=MG,Gender="Male")
FemaleGene.df=data.frame(Name=rownames(FemaleExpr),symbol=FG,Gender="Female")
Gene.df=rbind(MaleGene.df,FemaleGene.df)
rownames(Gene.df)=Gene.df$Name
#Gene.df$FDAMic=ifelse(Gene.df$symbol%in%rownames(FDAMicMarker),"FDAMicMarker","OtherGenes")
ProliferationGenes=c("SPI1","CD81","CD74","TYROBP","CSF1R","CORO1A","RPS3","VSIG4","HLA-DPB1","HLA-A","HLA-E","HLA-DPA1","CEBPB")
ProliferationGenesIndex=na.omit(match(ProliferationGenes,Gene.df$symbol))
ha = rowAnnotation(foo = anno_mark(at = ProliferationGenesIndex, labels = Gene.df$symbol[ProliferationGenesIndex]))
data.df.expr=rbind(MaleExpr,FemaleExpr)
data.df.expr = t(scale(t(data.df.expr)))
geneAnno = rowAnnotation(Gender = Gene.df[,"Gender"],
    col = list(
          Gender = c("Female"="Violet","Male"="RoyalBlue"))
)
pdf("APOEGenotypeInfluence.pdf",height=4,width=4)
Heatmap(data.df.expr, cluster_rows = TRUE, left_annotation =geneAnno,right_annotation = ha,cluster_columns = FALSE, show_row_names=FALSE,row_names_side = "left")
dev.off()

Gene.df$Name=NULL
Gene.df$symbol=NULL
ann_colors = list(
    Gender = c(Female="Violet",Male="RoyalBlue")
)
pheatmap(data.df.expr,scale="row",annotation_colors = ann_colors,annotation_row=Gene.df,cluster_col=FALSE,border=FALSE,show_rownames=FALSE,clustering_method="ward.D2",treeheight_row=20,colorRampPalette(c("navy", "white","pink","firebrick3"))(50))

tmp=data.frame(data.df.expr,check.names=F)
tmp$Gender=Gene.df$Gender
tmp$Symbol=Gene.df$symbol
tmp$Pattern=ifelse(tmp[,3]-tmp[,1]>0,"Increase","Decrease")
tmp=tmp[order(tmp$Gender,tmp$Pattern),]
write.table(tmp,file="HeatmapData.txt",row.names=F,sep="\t",quote=F)


GeneExpr=read.table("D:/Alzheimer/Syn3388564/geneExpr.txt",header=TRUE,row.names=1,sep="\t",check.names=F)
SampleInfo=read.table("D:/Alzheimer/Syn3388564/Phenotype.txt",header=TRUE,row.names=1,sep="\t")
Samples=intersect(colnames(GeneExpr),rownames(SampleInfo))
GeneExpr=GeneExpr[rowSums(GeneExpr)>0,Samples]
SampleInfo=SampleInfo[Samples,]
all(colnames(GeneExpr)==rownames(SampleInfo))
BigData=cbind(SampleInfo,t(GeneExpr))

BigData=BigData[which(BigData$apoe_genotype!="NA"),]
gene=c("CD74","CEBPB","CORO1A","CSF1R","RPS3","SPI1","TYROBP","VSIG4") 
gene=intersect(gene,colnames(BigData))
tmp=BigData[,c("msex","apoe_genotype",gene)]
tmp=tmp[!tmp$apoe_genotype=="24",]
library(reshape2)
tmp1=melt(tmp,id=c(1,2))
colnames(tmp1)=c("Gender","APOEType","Gene","Expression")
APOEGeneType=factor(tmp1$APOEType,levels=rev(c("22","23","33","34","44")))
#APOEGeneType=factor(tmp1$APOEType,levels=c("22|23","33","34|44"))
t=ggplot(tmp1, aes(x= APOEGeneType, y=log2(Expression+1), fill=as.character(APOEType))) + geom_violin(trim = FALSE,scale ="width") +
geom_boxplot(width = 0.3,outlier.shape = NA)+
coord_flip() + 
scale_fill_brewer(palette="Spectral",direction=-1)+
facet_grid(factor(Gender,levels=c("Male","Female"))~Gene, scales="free_x")+
theme_bw() + 
theme (strip.background = element_blank(), strip.text.x = element_text(size=rel(1.0)),strip.text.y = element_text(size=rel(1.0))) + 
theme (axis.text.x=element_text(size=rel(1.0))) + 
guides(fill=FALSE) +
scale_x_discrete(expand = c(0, 0)) + 
scale_y_continuous(expand = c(0, 0))
pdf("ProAssoGeneAcrossAPOE.pdf",width=7,height=4)
print(t)
dev.off()


APOEGeneType=factor(tmp1$APOEType,levels=c("22","23","33","34","44"))
#APOEGeneType=factor(tmp1$APOEType,levels=c("22|23","33","34|44"))
t=ggplot(tmp1, aes(x= APOEGeneType, y=log2(Expression+1), fill=as.character(APOEType))) + geom_violin(trim = FALSE,scale ="width") +
geom_boxplot(width = 0.3,outlier.shape = NA)+
ggpubr::stat_compare_means()+
scale_fill_brewer(palette="Spectral",direction=-1)+
facet_grid(Gene~factor(Gender,levels=c("Male","Female")) , scales="free_y")+
theme_bw() + 
#theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
theme (strip.background = element_blank(), strip.text.x = element_text(size=rel(1.0)),strip.text.y = element_text(size=rel(1.0))) + 
theme (axis.text.x=element_text(size=rel(1.0))) + 
guides(fill=FALSE) +
scale_x_discrete(expand = c(0, 1)) + 
scale_y_continuous(expand = c(0, 1))
pdf("ProAssoGeneAcrossAPOE_hor.pdf",width=4,height=7)
print(t)
dev.off()


############validation the marker genes in FDAMic in 5XFAD Mouse#########################
setwd("D:/Alzheimer/syn18485175")
group=read.table("D:/Alzheimer/Microglia/MouseFAD_GSE163857/rawData/phenotype.txt",header=T,row.names=1,sep="\t")
group$APOE=ifelse(grepl("APOE3",group$Genotype),"APOE3","APOE4")
group$Sex=group$Sex
group$Statues=ifelse(grepl("FAD",group$Genotype),"Alzheimer","Control")
group$Genotype=NULL
group$Age=NULL
group=group[,c("APOE","Sex","Statues")]
GeneExpr=read.table("D:/Alzheimer/Microglia/MouseFAD_GSE163857/rawData/GeneExpr.txt",header=T,row.names=1)
GeneExpr=GeneExpr[rowSums(GeneExpr)>10,]
FADMic=read.table("D:/Alzheimer/syn18485175/cluster13/monocle/ADMicSubCluster7Marker.txt",header=T,row.names=1)
FADMarker=FADMic[FADMic$avg_logFC>0&FADMic$p_val_adj<0.01,]
dim(FADMarker)
library("homologene")
FADMarker2Mouse=unique(human2mouse(rownames(FADMarker))$mouseGene)
library(clusterProfiler)
FADMarker2MouseEnsemble=bitr(FADMarker2Mouse,"SYMBOL","ENSEMBL",OrgDb="org.Mm.eg.db")[,2]
gene=intersect(FADMarker2MouseEnsemble,rownames(GeneExpr))
ann_colors = list(
    Sex = c(Female = "Magenta", Male = "CornflowerBlue"),
    APOE=c(APOE3 = "lightblue",APOE4 = "  DarkOrange"),
    Statues=c(Alzheimer= "Violet",Control="Blue")
)
t=pheatmap(log2(GeneExpr[gene,]+1),clustering_method="ward.D2",show_colnames=F,show_rownames=F,scale="row",annotation=group,annotation_colors = ann_colors,color = colorRampPalette(c("navy", "white", "firebrick3"))(50),border="white")
pdf("Manuscript/Microglia/Graph/Part5/FDAMicMarkerIn5xFADNames.pdf",height=4,width=5)
print(t)
dev.off()

FADMarker2MouseEnsemble=bitr(FADMarker2Mouse,"SYMBOL","ENSEMBL",OrgDb="org.Mm.eg.db")
GeneExprTmp=GeneExpr
GeneExprTmp$ENSEMBL=rownames(GeneExprTmp)
GeneExprTmpExpr=merge(GeneExprTmp,FADMarker2MouseEnsemble,by="ENSEMBL")
rownames(GeneExprTmpExpr)=GeneExprTmpExpr$SYMBOL
GeneExprTmpExpr$SYMBOL=NULL
GeneExprTmpExpr$ENSEMBL=NULL
t=pheatmap(log2(GeneExprTmpExpr+1),clustering_method="ward.D2",show_colnames=F,show_rownames=T,scale="row",annotation=group,annotation_colors = ann_colors,color = colorRampPalette(c("navy", "white", "firebrick3"))(50),border="white")
pdf("Manuscript/Microglia/Graph/Part5/FDAMicMarkerIn5xFADNames.pdf",height=9,width=10)
print(t)
dev.off()



############ check the APOE dependent microglia number #########################
ratio=read.table("D:/Alzheimer/syn18485175/Scaden/TrainingData/ratio.txt",header=T,row.names=1) #data from ROSMAP, script in D:\Alzheimer\syn18485175\Cibersortx\code.R
SampleInfo=read.table("D:/Alzheimer/Syn3388564/Phenotype.txt",header=TRUE,row.names=1,sep="\t")
Samples=intersect(rownames(SampleInfo),rownames(ratio))
length(Samples)
ratio=ratio[Samples,]
SampleInfo=SampleInfo[Samples,]
all(rownames(ratio)==rownames(SampleInfo))
SampleInfo=SampleInfo[,c("apoe_genotype","age_at_visit_max")]
BigData=cbind(SampleInfo,ratio)
BigData=BigData[which(BigData$apoe_genotype!="NA"),]
Sex=factor(BigData$msex,levels=c("Male","Female"))
Statues=factor(BigData$Condition,levels=c("Control","Alzheimer"))
APOE=factor(BigData$apoe_genotype,levels=c("22","23","24","33","34","44"))
t=ggplot(BigData, aes(x=APOE, y=Mic*100, fill=Sex)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_dotplot(binaxis='y', stackdir='centerwhole',  position=position_dodge(0.7),dotsize=3,binwidth=0.05)+
  scale_fill_manual(values = c("RoyalBlue","Violet"))+
  theme_bw()+ylim(1,12)+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0)))
pdf("Manuscript/Microglia/Graph/Part5/MicCellTypeFromROSMAP8ScadenSplit8APOE.pdf",height=4,width=6)
print(t)
dev.off()

t.test(BigData[BigData$msex=="Female"&BigData$apoe_genotype=="44","Mic"],BigData[BigData$msex=="Male"&BigData$apoe_genotype=="44","Mic"]) #p-value =  0.02662
t.test(BigData[BigData$msex=="Female"&BigData$apoe_genotype=="34","Mic"],BigData[BigData$msex=="Male"&BigData$apoe_genotype=="34","Mic"]) #p-value = 0.139
t.test(BigData[BigData$msex=="Female"&BigData$apoe_genotype=="33","Mic"],BigData[BigData$msex=="Male"&BigData$apoe_genotype=="33","Mic"]) #p-value = 0.0467
t.test(BigData[BigData$msex=="Female"&BigData$apoe_genotype=="23","Mic"],BigData[BigData$msex=="Male"&BigData$apoe_genotype=="23","Mic"]) #p-value = 0.4721

t.test(BigData[BigData$msex=="Female"&BigData$apoe_genotype=="44","Mic"],BigData[BigData$msex=="Female"&BigData$apoe_genotype=="34","Mic"]) #p-value = 0.01803
t.test(BigData[BigData$msex=="Female"&BigData$apoe_genotype=="44","Mic"],BigData[BigData$msex=="Female"&BigData$apoe_genotype=="33","Mic"]) #p-value = 0.02066
t.test(BigData[BigData$msex=="Female"&BigData$apoe_genotype=="44","Mic"],BigData[BigData$msex=="Female"&BigData$apoe_genotype=="23","Mic"]) #p-value = 0.01145
t.test(BigData[BigData$msex=="Female"&BigData$apoe_genotype=="44","Mic"],BigData[BigData$msex=="Female"&BigData$apoe_genotype=="22","Mic"]) #p-value = 0.01145

t.test(BigData[BigData$msex=="Female"&BigData$apoe_genotype=="34","Mic"],BigData[BigData$msex=="Female"&BigData$apoe_genotype=="33","Mic"]) #p-value = 0.02066
t.test(BigData[BigData$msex=="Female"&BigData$apoe_genotype=="34","Mic"],BigData[BigData$msex=="Female"&BigData$apoe_genotype=="23","Mic"]) #p-value = 0.01145
t.test(BigData[BigData$msex=="Female"&BigData$apoe_genotype=="34","Mic"],BigData[BigData$msex=="Female"&BigData$apoe_genotype=="22","Mic"]) #p-value = 0.01145


t.test(BigData[BigData$msex=="Male"&BigData$apoe_genotype=="44","Mic"],BigData[BigData$msex=="Male"&BigData$apoe_genotype=="34","Mic"]) #p-value = 0.8372
t.test(BigData[BigData$msex=="Male"&BigData$apoe_genotype=="44","Mic"],BigData[BigData$msex=="Male"&BigData$apoe_genotype=="33","Mic"]) #p-value = 0.9646
t.test(BigData[BigData$msex=="Male"&BigData$apoe_genotype=="44","Mic"],BigData[BigData$msex=="Male"&BigData$apoe_genotype=="23","Mic"]) #p-value = 0.7905


wilcox.test(BigData[BigData$msex=="Female"&BigData$apoe_genotype=="34","Mic"],BigData[BigData$msex=="Female"&BigData$apoe_genotype=="33","Mic"]) #p-value = 0.01145

