library(Seurat)
library(clusterProfiler)
library(SCENIC)
library(ggrepel)
library(ggpubr)
library(presto)
library(dplyr)
library(msigdbr)
library(tibble)
library(fgsea)
library(homologene)
library(scCustomize)


setwd("/Projects/deng/Alzheimer/syn18485175")
ADMic=readRDS("ADMic.rds")
Idents(ADMic)=ADMic$SubType
#BAMic   HomMic NorARMic    ARMic   FDAMic   DysMic   OpcMic
#158      998      154      136      148      117       24
DEG=FindMarkers(ADMic,ident.1="BAMic",pct.min=0.25)
DEGTmp=DEG[DEG$p_val<0.01,]
write.table(DEGTmp,file="Graph/Part3/Marker/BAMic.txt",sep="\t",quote=F)
DEG=FindMarkers(ADMic,ident.1="HomMic",pct.min=0.25)
DEGTmp=DEG[DEG$p_val<0.01,]
write.table(DEGTmp,file="Graph/Part3/Marker/HomMic.txt",sep="\t",quote=F)
DEG=FindMarkers(ADMic,ident.1="NorARMic",pct.min=0.25)
DEGTmp=DEG[DEG$p_val<0.01,]
write.table(DEGTmp,file="Graph/Part3/Marker/NorARMic.txt",sep="\t",quote=F)
DEG=FindMarkers(ADMic,ident.1="ARMic",pct.min=0.25)
DEGTmp=DEG[DEG$p_val<0.01,]
write.table(DEGTmp,file="Graph/Part3/Marker/ARMic.txt",sep="\t",quote=F)
DEG=FindMarkers(ADMic,ident.1="FDAMic",pct.min=0.25)
DEGTmp=DEG[DEG$p_val<0.01,]
write.table(DEGTmp,file="Graph/Part3/Marker/FDAMic.txt",sep="\t",quote=F)
DEG=FindMarkers(ADMic,ident.1="DysMic",pct.min=0.25)
DEGTmp=DEG[DEG$p_val<0.01,]
write.table(DEGTmp,file="Graph/Part3/Marker/DysMic.txt",sep="\t",quote=F)

BAMic=read.table("Graph/Part3/Marker/BAMic.txt",header=T,row.names = NULL)
BAMic$Pattern=ifelse(BAMic$avg_log2FC>0,"UpInBAMic","DnInBAMic")
HomMic=read.table("Graph/Part3/Marker/HomMic.txt",header=T,row.names = NULL)
HomMic$Pattern=ifelse(HomMic$avg_log2FC>0,"UpInHomMic","DnInHomMic")
NorARMic=read.table("Graph/Part3/Marker/NorARMic.txt",header=T,row.names = NULL)
NorARMic$Pattern=ifelse(NorARMic$avg_log2FC>0,"UpInNorARMic","DnInNorARMic")
FDAMic=read.table("Graph/Part3/Marker/FDAMic.txt",header=T,row.names = NULL)
FDAMic$Pattern=ifelse(FDAMic$avg_log2FC>0,"UpInFDAMic","DnInFDAMic")
ARMic=read.table("Graph/Part3/Marker/ARMic.txt",header=T,row.names = NULL)
ARMic$Pattern=ifelse(ARMic$avg_log2FC>0,"UpInARMic","DnInARMic")
DysMic=read.table("Graph/Part3/Marker/DysMic.txt",header=T,row.names = NULL)
DysMic$Pattern=ifelse(DysMic$avg_log2FC>0,"UpInDysMic","DnInDysMic")
DEGTotal=rbind(HomMic,BAMic,NorARMic,FDAMic,ARMic,DysMic)
DEGList=split(DEGTotal$`row.names`,DEGTotal$Pattern)
setList=c("UpInBAMic","UpInHomMic","UpInNorARMic","UpInARMic","UpInFDAMic","UpInDysMic","DnInBAMic","DnInHomMic","DnInNorARMic","DnInARMic","DnInFDAMic","DnInDysMic")
t=upset(fromList(DEGList),nsets = 12,
  sets=rev(setList),
  nintersects = 41,
  keep.order=TRUE,
  main.bar.color=c(rep("Violet",11),rep("lightgrey",85)),
  matrix.color=c("lightgrey"),
  sets.bar.color=rev(c(rep("firebrick3",6),rep("SlateBlue",6)))
  )
pdf("Graph/Part3/OverlapInARMicBranchAll.pdf",height=5)
print(t)
dev.off()
UniqueUpInFDAMic=setdiff(DEGList$UpInFDAMic,c(DEGList$UpInARMic,DEGList$UpInDysMic,DEGList$DnInHomMic,DEGList$DnInNorARMic,DEGList$DnInBAMic))
length(UniqueUpInFDAMic) #31
write.table(UniqueUpInFDAMic,file="Graph/Part3/31_UniqueUpInFDAMic.txt",sep="\t",row.names=F,quote=F,col.names=F)
UniqueDnInFDAMic=setdiff(DEGList$DnInFDAMic,c(DEGList$UpInBAMic,DEGList$UpInHomMic,DEGList$UpInNorARMic,DEGList$UpInARMic,DEGList$UpInFDAMic,DEGList$UpInDysMic,DEGList$DnInBAMic,DEGList$DnInHomMic,DEGList$DnInNorARMic,DEGList$DnInARMic,DEGList$DnInDysMic))
length(UniqueDnInFDAMic) #72
write.table(UniqueDnInFDAMic,file="Graph/Part3/72_UniqueDnInFDAMic.txt",sep="\t",row.names=F,quote=F,col.names=F)

tmp=intersect(DEGList$UpInFDAMic,DEGList$DnInHomMic)
UniqueUpInFDAMicDnInHomMic=setdiff(tmp,c(DEGList$UpInBAMic,DEGList$UpInHomMic,DEGList$UpInNorARMic,DEGList$UpInARMic,DEGList$UpInDysMic,DEGList$DnInBAMic,DEGList$DnInNorARMic,DEGList$DnInARMic,DEGList$DnInDysMic))
length(UniqueUpInFDAMicDnInHomMic)
tiff("Graph/Part3/31_UniqueUpInFDAMic.tiff",height=2000,width=1000)
VlnPlot(ADMic,group.by="SubType",features=UniqueUpInFDAMic)&NoLegend()
dev.off()
tiff("Graph/Part3/6_UniqueGene.tiff",height=600,width=500)
VlnPlot(ADMic,group.by="SubType",features=UniqueUpInFDAMicDnInHomMic)&NoLegend()
dev.off()
pdf("Graph/Part3/31_UniqueUpInFDAMic.pdf",height=3,width=9)
DotPlot(ADMic,group.by="SubType",features=UniqueUpInFDAMic)&RotatedAxis()
dev.off()
pdf("Graph/Part3/6_UniqueGene.pdf",height=3,width=5)
DotPlot(ADMic,group.by="SubType",features=UniqueUpInFDAMicDnInHomMic)&RotatedAxis()
dev.off()
pdf("Graph/Part3/72_UniqueDnInFDAMic.pdf",height=3,width=18)
DotPlot(ADMic,group.by="SubType",features=UniqueDnInFDAMic)&RotatedAxis()
dev.off()

gene=c("AKAP13","GAB2","ANKRD17","TCF12","KANSL1","RAPGEF1","LDLRAD4","ZNF609","PIK3R5","ARRB2","USP15","USP39")
pdf("Graph/Part3/targetGene.pdf")
VlnPlot(ADMic,group.by="SubType",features=gene)&NoLegend()
dev.off()
colors_list <- c("forestgreen","NavajoWhite","LightBlue", "SlateBlue", "Violet", "gold")
pdf("Graph/Part3/targetGeneVln.pdf",height=8,width=4.5)
scCustomize::Stacked_VlnPlot(seurat_object = ADMic, features = gene, x_lab_rotate = TRUE, plot_spacing = 0.05,colors_use = colors_list)&
theme(axis.title.x=element_blank())
dev.off()

SharedUpGene=Reduce(intersect,list(DEGList$UpInC9,DEGList$UpInC7,DEGList$UpInC10))
write.table(SharedUpGene,file="Graph/Part3/52_SharedGene.txt",sep="\t",row.names=F,quote=F,col.names=F)

UniqueUpInC7=setdiff(DEGList$UpInC7,c(DEGList$UpInC9,DEGList$UpInC10))
targetGeneInfo=C7DEG[C7DEG$row.names%in%UniqueUpInC7,c("row.names","pct.1","pct.2","p_val")]
targetGeneInfo=targetGeneInfo[order(targetGeneInfo$p_val),]
write.table(targetGeneInfo,file="Graph/Part3/UniqueUpInC7.txt",sep="\t",row.names=F,quote=F)
UniqueDnInC7=setdiff(DEGList$DnInC7,c(DEGList$DnInC9,DEGList$DnInC10,DEGList$UpInC9,DEGList$UpInC10))
targetGeneInfo=C7DEG[C7DEG$row.names%in%UniqueDnInC7,c("row.names","pct.1","pct.2","p_val")]
targetGeneInfo=targetGeneInfo[order(targetGeneInfo$p_val),]
write.table(targetGeneInfo,file="Graph/Part3/UniqueDnInC7.txt",sep="\t",row.names=F,quote=F)

UniqueUpInC9=setdiff(DEGList$UpInC9,c(DEGList$UpInC7,DEGList$UpInC10,DEGList$DnInC7,DEGList$DnInC10))
length(UniqueUpInC9)
targetGeneInfo=C9DEG[C9DEG$row.names%in%UniqueUpInC9,c("row.names","pct.1","pct.2","p_val")]
targetGeneInfo=targetGeneInfo[order(targetGeneInfo$p_val),]
write.table(targetGeneInfo,file="Graph/Part3/UniqueUpInC9.txt",sep="\t",row.names=F,quote=F)
UniqueDnInC9=setdiff(DEGList$DnInC9,c(DEGList$DnInC7,DEGList$DnInC10,DEGList$UpInC10))
length(UniqueDnInC9)
targetGeneInfo=C9DEG[C9DEG$row.names%in%UniqueDnInC9,c("row.names","pct.1","pct.2","p_val")]
targetGeneInfo=targetGeneInfo[order(targetGeneInfo$p_val),]
write.table(targetGeneInfo,file="Graph/Part3/UniqueDnInC9.txt",sep="\t",row.names=F,quote=F)

UniqueUpInC10=setdiff(DEGList$UpInC10,c(DEGList$UpInC7,DEGList$UpInC9,DEGList$DnInC7,DEGList$DnInC9))
length(UniqueUpInC10)
targetGeneInfo=C10DEG[C10DEG$row.names%in%UniqueUpInC10,c("row.names","pct.1","pct.2","p_val")]
targetGeneInfo=targetGeneInfo[order(targetGeneInfo$p_val),]
write.table(targetGeneInfo,file="Graph/Part3/UniqueUpInC10.txt",sep="\t",row.names=F,quote=F)
UniqueDnInC10=setdiff(DEGList$DnInC10,c(DEGList$DnInC7,DEGList$DnInC9,DEGList$UpInC9))
length(UniqueDnInC10)
targetGeneInfo=C10DEG[C10DEG$row.names%in%UniqueDnInC10,c("row.names","pct.1","pct.2","p_val")]
targetGeneInfo=targetGeneInfo[order(targetGeneInfo$p_val),]
write.table(targetGeneInfo,file="Graph/Part3/UniqueDnInC10.txt",sep="\t",row.names=F,quote=F)

targetGene=UniqueUpInC7
targetGene=c("HCK","ITPR2","PAK1","PLD1","ELMO1","CYFIP1","ABR","RAPGEF1","DOCK4","RIN2","DOCK8","RASGEF1C","MAPKAP1","TRIO","GNG7")
targetGeneInfo=C7DEG[C7DEG$row.names%in%targetGene,c("row.names","pct.1","pct.2","p_val")]
targetGeneInfo=targetGeneInfo[order(targetGeneInfo$p_val),]
targetGeneInfo=targetGeneInfo[1:10,]
targetGeneInfo.df=reshape2::melt(targetGeneInfo,id=c(1,4))
colnames(targetGeneInfo.df)=c("Gene","Pval","pctGroup","pctValue")
targetGeneInfo.df$logP=-log10(targetGeneInfo.df$Pval)
targetGeneInfo.df$Gene=factor(targetGeneInfo.df$Gene,levels=unique(targetGeneInfo.df$Gene))
t=ggdotchart(targetGeneInfo.df, x = "Gene", y ="logP",
           color = "pctGroup", group="pctGroup", # here it is
           palette = c('Violet','lightgrey'), 
           size = "pctValue", 
           add = "segment", 
           title = NULL,
           sorting="none",
           add.params = list(color = "lightgray", size = 0.5),
           position = position_dodge(1),
           ggtheme = theme_bw())&scale_size(range = c(0, 3))&
           theme(axis.title.x = element_blank(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              legend.position="none")
pdf("Graph/Part3/topUpInC7.pdf",height=2,width=3)
print(t)
dev.off()


ADMicTmp=subset(ADMic,ident=c(1:11))
new.cluster.ids <- c("BAMic", "HomMic", "HomMic", "HomMic", "HomMic", "HomMic", "NorARMic", "ARMic","FDAMic","DysMic"," ")
names(new.cluster.ids) <- levels(ADMicTmp)
ADMicTmp <- RenameIdents(ADMicTmp, new.cluster.ids)

DEGFDAMicVsDysMic=FindMarkers(ADMicTmp,ident.1=c("FDAMic"),ident.2=("DysMic"),pct.min=0,logfc.threshold =0)
hs_data=DEGFDAMicVsDysMic
#Down SigDown   SigUp      Up
#    806     175      40     751

DEGFDAMicVsARM=FindMarkers(ADMicTmp,ident.1=c("FDAMic"),ident.2=("ARMic"),pct.min=0,logfc.threshold =0)
hs_data=DEGFDAMicVsARM
#Down SigDown   SigUp      Up
#    910     822      13    1167

DEGFDAMic=FindMarkers(ADMicTmp,ident.1=c("FDAMic"),pct.min=0,pct.min=0.25)
hs_data=DEGFDAMic
#Down SigDown   SigUp      Up
#418     190     183     226

hs_data$threshold = as.factor(ifelse(hs_data$avg_log2FC<0,ifelse(hs_data$p_val<0.01&abs(hs_data$avg_log2FC)>0.25,"SigDown","Down"),ifelse(hs_data$p_val<0.01&abs(hs_data$avg_log2FC)>0.25,"SigUp","Up")))
table(hs_data$threshold)
#Down NoSignificant            Up
#46         17853            27
hs_data$ID=rownames(hs_data)
t=ggplot(data = hs_data, aes(x = avg_log2FC, y = -log10(p_val), colour=threshold, label =ID )) +
  geom_point(alpha=0.6, size=1) +
  theme_bw() + 
  scale_color_manual(values=c("lightblue","blue", "red"," Salmon")) +
  geom_vline(xintercept=c(-0.25,0.25),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)",title="Differential Expressed Gene") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank()) +  
    geom_text_repel(
    data = subset(hs_data, hs_data$p_val < 0.01 & abs(hs_data$avg_log2FC) >= 0.25),
    aes(label = ID),
    size = 3,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
pdf("Graph/Part3/DEGFDAMicSymbol.pdf",height=10,width=12)
print(t)
dev.off()


DEGFDAMic_order=DEGFDAMic[order(DEGFDAMic$p_val),]
colors_list <- c("forestgreen","NavajoWhite","LightBlue", "SlateBlue", "Violet", "gold")
pdf("Graph/Part3/TopUpGeneInFDAMicExpr.pdf",height=8,width=3.5)
VlnPlot(subset(ADMic,MicCluster%in%c(1:10)),features=rownames(DEGFDAMic_order)[1:10],group.by="SubType",col=colors_list,pt.size=0,ncol=1)&theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_blank(),plot.title = element_blank())
dev.off()

pdf("Graph/Part3/TopUpGeneInFDAMicExpr.pdf",height=8,width=4.5)
scCustomize::Stacked_VlnPlot(seurat_object = subset(ADMic,MicCluster%in%c(1:10)), features = rownames(DEGFDAMic_order)[1:10], x_lab_rotate = TRUE, plot_spacing = 0.05,colors_use = colors_list)&
theme(axis.title.x=element_blank())
dev.off()


DEGARM=FindMarkers(ADMicTmp,ident.1="ARMic",pct.min=0.25)
DEGFDAMic=FindMarkers(ADMicTmp,ident.1="FDAMic",pct.min=0.25)
DEGDysMic=FindMarkers(ADMicTmp,ident.1="DysMic",pct.min=0.25)

DEGARMVsFDAMic=FindMarkers(ADMicTmp,ident.1=c("ARMic"),ident.2=("FDAMic"),pct.min=0.25)
DEGDysVsFDAMic=FindMarkers(ADMicTmp,ident.1=c("DysMic"),ident.2=("FDAMic"),pct.min=0.25)

DEGARM=DEGARM[DEGARM$p_val<0.01&DEGARM$avg_log2FC> 0,]
DEGFDAMic=DEGFDAMic[DEGFDAMic$p_val<0.01&DEGFDAMic$avg_log2FC> 0,]
DEGDysMic=DEGDysMic[DEGDysMic$p_val<0.01&DEGDysMic$avg_log2FC> 0,]

DEGDysVsFDAMicTmp=DEGDysVsFDAMic[DEGDysVsFDAMic$p_val<0.01&DEGDysVsFDAMic$avg_log2FC> 0,]
DEGFDAMicVsDys=DEGDysVsFDAMic[DEGDysVsFDAMic$p_val<0.01&DEGDysVsFDAMic$avg_log2FC< 0,]
DEGARMVsFDAMic=DEGARMVsFDAMic[DEGARMVsFDAMic$p_val<0.01&DEGARMVsFDAMic$avg_log2FC> 0,]

DegALL=list(
  "ARMic"=rownames(DEGARM),"FDAMic"=rownames(DEGFDAMic),"DysMic"=rownames(DEGDysMic),
  "ARMicvsFDAMic"=rownames(DEGARMVsFDAMic),
  "DysVsFDAMic"=rownames(DEGDysVsFDAMicTmp),
  "FDAMicVSDys"=rownames(DEGFDAMicVsDys))

DegFunction <- compareCluster(geneCluster = DegALL, fun = "enrichGO",ont="BP",OrgDb='org.Hs.eg.db',minGSSize = 5,pvalueCutoff=0.01,keyType = "SYMBOL")
DegFunction1 <- clusterProfiler::simplify(DegFunction,cutoff = 0.5,select_fun = min,measure = "Wang")
pdf("Graph/Part3/ARMicFDAMicCompairsionGOBP.pdf",width=10,height=10)
dotplot(DegFunction1,showCategory=20,label_format = 100)
dev.off()

write.table(DegFunction1,file="Graph/Part3/ARMicFDAMicCompairsionGOBP.txt",sep="\t",quote=F)

ADMicTmp=subset(ADMic,SubType%in%c("BAMic","HomMic","NorARMic","ARMic","FDAMic","DysMic"))
Proliferation=c("CSF1R","CD81","CD74","VSIG4","HLA-A","HLA-E","CEBPB","BST2","HLA-DPA1","CORO1A","HLA-DPB1","TYROBP","RPS3")
t=DotPlot(ADMicTmp,features=Proliferation)
g=ggplot(t$data, aes(id, factor(features.plot,levels=rev(unique(features.plot))), size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_colour_viridis_c()+
labs(color="avg.exp,scaled",size="pct.exp",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),angle=90,vjust=0.5,hjust=1),axis.text.y = element_text(size=rel(1.0))) 
pdf("Graph/Part3/ProliferationGene.pdf",width=4,height=4)
print(g)
dev.off()

NeuronDeath=c("APOE","GPX1","GRN","CEBPB","HSPA5","CORO1A","TREM2","ARRB2","SNCA","HSP90AB1","TYROBP","RB1")
t=DotPlot(ADMicTmp,features=NeuronDeath)
g=ggplot(t$data, aes(id, factor(features.plot,levels=rev(unique(features.plot))), size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_colour_viridis_c()+
labs(color="avg.exp.scaled",size="avg.exp",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),angle=90,vjust=0.5,hjust=1),axis.text.y = element_text(size=rel(1.0))) 
pdf("Graph/Part3/NeuronDeathGene.pdf",width=4,height=4)
print(g)
dev.off()

Antigen=c("HLA-B","CD74","HLA-DRA","CALR","HLA-A","HLA-C","HLA-E","HLA-DRB5","HLA-DPA1","B2M","HLA-DPB1","TREM2","TAPBP")
t=DotPlot(ADMicTmp,features=Antigen)
g=ggplot(t$data, aes(id, factor(features.plot,levels=rev(unique(features.plot))), size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_colour_viridis_c()+
labs(color="avg.exp.scaled",size="avg.exp",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),angle=90,vjust=0.5,hjust=1),axis.text.y = element_text(size=rel(1.0))) 
pdf("Graph/Part3/AntigenGene.pdf",width=4,height=4)
print(g)
dev.off()

MacrophageActivation=c("VSIG4","CD74","C1QA","GRN","HAMP","CEBPA","TREM2","SNCA","TYROBP")
t=DotPlot(ADMicTmp,features=MacrophageActivation)
g=ggplot(t$data, aes(id, factor(features.plot,levels=rev(unique(features.plot))), size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_colour_viridis_c()+
labs(color="avg.exp.scaled",size="avg.exp",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),angle=90),axis.text.y = element_text(size=rel(1.0))) 
pdf("Graph/Part3/MacrophageActivationGene.pdf",width=4,height=4)
print(g)
dev.off()

phagocytosis=c("CD14","FCGR3A","FCGR1A","SLC11A1","HCK","CLCN3","RAB39A","ITGB1","PRKCD","PIK3CA","PAK1","HMGB1","UNC13D","ARHGAP12","PIK3R1","P2RY6","PIP5K1A")
t=DotPlot(ADMicTmp,features=phagocytosis)
g=ggplot(t$data, aes(id, factor(features.plot,levels=rev(unique(features.plot))), size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_colour_viridis_c()+
labs(color="avg.exp.scaled",size="pct.exp",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),angle=90,vjust=0.5,hjust=1),axis.text.y = element_text(size=rel(1.0))) 
pdf("Graph/Part3/PhagocytosisGene.pdf",width=4,height=4)
print(g)
dev.off()

Insulin=c("GRB2","PIK3R1","PYGL","PIK3R5","FLOT1","MKNK1","RPTOR","SOS1","ACACA","PIK3CD","CBL","RHEB","HK2","PTPN1","RAF1","PRKACB","PPP1CB","IRS2","RHOQ","PRKAG2","PRKAG1")
t=DotPlot(ADMicTmp,features=Insulin)
g=ggplot(t$data, aes(id, factor(features.plot,levels=rev(unique(features.plot))), size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_colour_viridis_c()+
labs(color="avg.exp.scaled",size="pct.exp",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),angle=90,vjust=0.5,hjust=1),axis.text.y = element_text(size=rel(1.0))) 
pdf("Graph/Part3/InsulinGene.pdf",width=4,height=4)
print(g)
dev.off()

UpInFDA=read.table("Graph/Part3/KEGG/UpInFDANetwork_edge.txt",sep="\t")
UpInFDANumber=data.frame(sort(table(UpInFDA[,2])))
colnames(UpInFDANumber)=c("Gene","Number")
g=ggplot(UpInFDANumber, aes(Gene, y=Number,fill=Number)) +geom_bar(stat="identity")+
scale_fill_viridis_c(option="B")+coord_flip()+
labs(color="Number",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0))) 
pdf("Graph/Part3/KEGG/UpInFDAShareNumer.pdf",width=3,height=7)
print(g)
dev.off()
t=DotPlot(ADMicTmp,features=unique(UpInFDA[,2]))
g=ggplot(t$data, aes(id, factor(features.plot,levels=unique(UpInFDANumber[,1])), size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_colour_viridis_c()+
labs(color="avg.exp.scaled",size="avg.exp",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),angle=90,vjust=0.5,hjust=1,color="black"),axis.text.y = element_text(size=rel(1.0),face="italic",color="black")) 
pdf("Graph/Part3/KEGG/UpInFDA.pdf",width=4,height=8)
print(g)
dev.off()

DnInFDA=read.table("Graph/Part3/KEGG/DnInFDANetwork_edge.txt",sep="\t")
colnames(DnInFDA)=c("pathway","gene")
pathwayList=unique(DnInFDA$pathway)
for(i in 1:length(pathwayList)){
  geneList=DnInFDA[DnInFDA$pathway==pathwayList[i],"gene"]
  t=DotPlot(ADMicTmp,features=geneList)
  g=ggplot(t$data, aes(id, features.plot, size= pct.exp,color=avg.exp.scaled)) +geom_point()+
  scale_colour_viridis_c()+
  labs(color="avg.exp.scaled",size="avg.exp",x="",y="",title="")+
  theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),angle=90,vjust=0.5,hjust=1,color="black"),axis.text.y = element_text(size=rel(1.0),face="italic",color="black")) 
  h=max(length(geneList)/6,3)
  pdf(paste0("Graph/Part3/KEGG/DnInFDA_",pathwayList[i],".pdf",sep=""),width=4,height=h)
  print(g)
  dev.off()
}


##############GSEA analysis on KEGG pathways##########################
m_df<- msigdbr(species = "Homo sapiens", category = "C2",subcategory="CP:KEGG") 
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
ADdeg <- wilcoxauc(ADMicTmp, 'SubType')
cluster="FDAMic"
cluster="ARMic"
clusterCell<- ADdeg %>% dplyr::filter(group == cluster) %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
ranks<- deframe(clusterCell)
fgseaRes<- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0)
targetPathways=c("KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",
  "KEGG_FC_EPSILON_RI_SIGNALING_PATHWAY",
  "KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS",
  "KEGG_INSULIN_SIGNALING_PATHWAY")

pdf("Graph/Part3/SigPathyInARMic.pdf",height=4,width=9)
plotGseaTable(
  fgsea_sets[targetPathways],
  ranks,
  fgseaRes,
  gseaParam = 1,
  colwidths = c(7, 3, 0.8, 1.0, 1.0),
  render = TRUE
)
dev.off()


for(cluster in unique(ADMicTmp$Type)){
print (cluster)
clusterCell<- ADdeg %>% dplyr::filter(group == cluster) %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
ranks<- deframe(clusterCell)
fgseaRes<- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0)
fwrite(fgseaRes, file=paste0(cluster,".txt",sep=""), sep="\t", sep2=c("", " ", ""))
}

## manually order
data=read.table("D:/Alzheimer/syn18485175/cluster13/monocle/GSEA/cellType/combineDetail.txt",header=T,sep="\t")
data$sig=ifelse(data$pval<0.05,"Sig","NonSig")
order=read.table("Graph/Part4/FunctionHeatmapOrder.txt",sep="\t",header=TRUE)
Pathway=factor(data$pathway,levels=order[,1])
cellType=factor(data$CellType,levels=c("BAMic","HomMic","NorARMic","ARMic","FDAMic","DysMic"))
t=ggplot(data,aes(cellType,Pathway,size=-1*log(pval),color=NES,fill=NES,shape=sig))+geom_point()+
scale_color_gradient2(low="blue",mid="white",high = "red")+
scale_fill_gradient2(low="blue",mid="white",high = "red")+
scale_shape_manual(values=c(3,19))+
theme_bw()+theme(axis.text.x = element_text(size=rel(1.0),angle = 90),axis.text.y = element_text(size=rel(1.0)))

##order by kegg category
dataT=read.table("D:/Alzheimer/syn18485175/cluster13/monocle/GSEA/cellType/combineDetail.txt",header=T,sep="\t")
dataT$sig=ifelse(dataT$pval<0.05,"Sig","NonSig")
pathwayCategory=read.table("D:/Public/GSEA/KEGGPathwayCategory.txt",header=T,sep="\t")
all(dataT$pathway%in%pathwayCategory$Pathway)
dataT$Pathway=dataT$pathway
dataT=merge(dataT,pathwayCategory,by="Pathway")
dplyr::arrange(dataT,Group, SubGroup, desc(NES) )
dataT=dplyr::arrange(dataT,Group, SubGroup, desc(NES) )
dataT$PathwayName=factor(dataT$PathwayName,levels=rev(unique(dataT$PathwayName)))
cellType=factor(dataT$CellType,levels=c("BAMic","HomMic","NorARMic","ARMic","FDAMic","DysMic"))
t=ggplot(dataT,aes(cellType,PathwayName,size=-1*log(pval),color=NES,fill=NES,shape=sig))+geom_point()+
scale_color_gradient2(low="blue",mid="white",high = "red")+
scale_fill_gradient2(low="blue",mid="white",high = "red")+
scale_shape_manual(values=c(3,19))+
theme_bw()+
theme(
  axis.text.x = element_text(size=rel(1.0),colour ="black",angle = 90,vjust=0.5,hjust=1),
  axis.text.y = element_text(size=rel(1.0),colour ="black"),
  axis.title.x = element_blank(),
  axis.title.y = element_blank()) 
pdf("Graph/Part3/GSEA4MicTypeDetailOrder8Group.pdf",height=9,width=7)
print(t)
dev.off()
write.table(dataT,file="Graph/Part3/GSEA4MicTypeDetailOrder8Group.txt",sep="\t",quote=F,row.names=F)

data=data[,1:8]
data$sig=ifelse(data$pval<0.05,"Sig","NonSig")
PathwayList=factor(data$pathway,levels=rev(c("ARMMic","AgedMic","DAMic","DystrophicMic","VMAMic","RiskGene","LPS","MGnD","Neurodegeneration","Proliferation","Neutrophil","Macrophage","Interferon")))
ClusterList=factor(data$Cluster,levels=c("C4","C1","C2","C3","C8","C6","C5","C9","C7","C10"))
t=ggplot(data,aes(ClusterList,PathwayList,size=-1*log(pval),colour=NES,shape=sig))+geom_point()+
scale_color_gradient2(low="blue",mid="white",high = "red")+
scale_shape_manual(values=c(3,19))+
theme_bw()+#theme(legend.position="bottom")+
scale_size(range = c(3, 6))+
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=rel(1.0),colour = "black",angle=90, vjust = 0.5, hjust=1),axis.text.y = element_text(size=rel(1.0),colour = "black"))


#----Compaired to privious validated subtypes------------
ADMic=readRDS("ADMic.rds")
Idents(ADMic)=ADMic$MicCluster
ADMicTmp=subset(ADMic,ident=c(1:10))
ident=c(4,1,2,3,8,6,5,9,7,10)
ADMicTmp@active.ident <- factor(ADMic@active.ident,levels=ident)
new.cluster.ids <- c("BAMic", "HomMic", "HomMic", "HomMic", "HomMic", "HomMic", "NorARMic", "ARMic","FDAMic","DysMic")
names(new.cluster.ids) <- levels(ADMicTmp)
ADMicTmp <- RenameIdents(ADMicTmp, new.cluster.ids)

RISKData=read.table("/Projects/deng/Alzheimer/syn18485175/cluster13/GWAS.txt",sep="\t")
RiskGene=unique(RISKData[,1])
DAMData=read.table("/Projects/deng/Alzheimer/syn18485175/cluster13/DAM8Keren.txt",header=T,row.names=1,sep="\t")
DAMData=DAMData[DAMData[,5]>0,]
DAMGene=unique(mouse2human(rownames(DAMData))$humanGene)
AgedData=read.table("/Projects/deng/Alzheimer/syn18485175/cluster13/AgedMic8Olah.txt",header=T,row.names=1,sep="\t")
AgedGene=unique(AgedData$Genesymbol)
DystroData=read.table("/Projects/deng/Alzheimer/syn18485175/cluster13/DystrophicMic8Nguyen.txt",header=T,sep="\t")
DystroData=DystroData[DystroData$logfoldchanges>0,]
DystroGene=DystroData$names
ArmMicData=read.table("/Projects/deng/Alzheimer/syn18485175/cluster13/ARMDeg8Frigerio.txt",header=T,sep="\t")
ArmMicData=ArmMicData[ArmMicData$avg_logFC>0,]
ARMGene=unique(mouse2human(rownames(ArmMicData))$humanGene)
VMAData=read.table("/Projects/deng/Alzheimer/syn18485175/cluster13/VMA8Safaiyan.txt",header=T,row.names=1,sep="\t")
VMAGene=unique(mouse2human(rownames(VMAData))$humanGene)
MGnDData=read.table("/Projects/deng/Alzheimer/syn18485175/cluster13/MGnD8Krasemann.txt")
MGnDData=MGnDData[MGnDData[,2]=="Up",]
MGnDene=unique(mouse2human(MGnDData[,1])$humanGene)

Interferon=read.table("/Projects/deng/Alzheimer/syn18485175/cluster13/geneModule/Interferon.txt")
LPS=read.table("/Projects/deng/Alzheimer/syn18485175/cluster13/geneModule/LPS.txt")
Macrophage=read.table("/Projects/deng/Alzheimer/syn18485175/cluster13/geneModule/Macrophage.txt")
Neurodegeneration=read.table("/Projects/deng/Alzheimer/syn18485175/cluster13/geneModule/Neurodegeneration.txt")
Neutrophil=read.table("/Projects/deng/Alzheimer/syn18485175/cluster13/geneModule/Neutrophil.txt")
Proliferation=read.table("/Projects/deng/Alzheimer/syn18485175/cluster13/geneModule/Proliferation.txt")

dats=list("MGnD"=MGnDene,"RiskGene"=RiskGene,"VMAMic"=VMAGene,"DAMic"=DAMGene,"AgedMic"=AgedGene,"DystrophicMic"=DystroGene,"ARMMic"=ARMGene,"Interferon"=Interferon[,1],
  "LPS"=LPS[,1],"Macrophage"=Macrophage[,1],"Neurodegeneration"=Neurodegeneration[,1],"Neutrophil"=Neutrophil[,1],"Proliferation"=Proliferation[,1])
fgsea_sets=dats


m_df<- msigdbr(species = "Homo sapiens", category = "C2",subcategory="CP:KEGG")  #CP:KEGG
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

ClusterMarker <- wilcoxauc(ADMicTmp, 'MicCluster')
table(ClusterMarker$group)
for(cluster in unique(ADMicTmp$MicCluster)){
print (cluster)
clusterCell<- ClusterMarker %>% dplyr::filter(group == cluster) %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
ranks<- deframe(clusterCell)
fgseaRes<- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0)
fwrite(fgseaRes, file=paste0("/Projects/deng/Alzheimer/syn18485175/Graph/Part4/GSEA/GSEA4Mivaluepe/C",cluster,".txt",sep=""), sep="\t", sep2=c("", " ", ""))
}

data=read.table("Microglia/Graph/Part4/GSEA4PreMictype.txt",header=T,sep="\t")
#data=data[order(data$NES,decreasing=T),]
data=data[,1:8]
data$sig=ifelse(data$pval<0.05,"Sig","NonSig")
PathwayList=factor(data$pathway,levels=rev(c("ARMMic","AgedMic","DAMic","DystrophicMic","VMAMic","RiskGene","LPS","MGnD","Neurodegeneration","Proliferation","Neutrophil","Macrophage","Interferon")))
ClusterList=factor(data$Cluster,levels=c("C4","C1","C2","C3","C8","C6","C5","C9","C7","C10"))
t=ggplot(data,aes(ClusterList,PathwayList,size=-1*log(pval),colour=NES,shape=sig))+geom_point()+
scale_color_gradient2(low="blue",mid="white",high = "red")+
scale_shape_manual(values=c(3,19))+
theme_bw()+#theme(legend.position="bottom")+
scale_size(range = c(3, 6))+
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=rel(1.0),colour = "black",angle=90, vjust = 0.5, hjust=1),axis.text.y = element_text(size=rel(1.0),colour = "black"))
pdf("Graph/Part3/GSEA4MictypeGraph.pdf",height=3,width=5)
print(t)
dev.off()

### expression of risk genes across the subtype of microglia ###
fgsea_sets=dats
Idents(ADMicTmp)=ADMicTmp$MicCluster
ADdeg <- wilcoxauc(ADMicTmp, 'SubType')
clusterCell<- ADdeg %>% dplyr::filter(group == "FDAMic") %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
ranks<- deframe(clusterCell)
fgseaRes<- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0)
pdf("Graph/Part3/GSEAPlotOfRiskGeneInFDAMic.pdf",height=2.5,width=4)
plotEnrichment(dats[["RiskGene"]],ranks) 
dev.off()
fgseaRes[fgseaRes$pathway=="RiskGene",]
#ARMic         pval         padj   log2err        ES      NES size
#RiskGene 9.870551e-05 0.0002138619 0.5384341 0.7187046 1.331341  101
#FDAMic       pval      padj    log2err         ES       NES size
#RiskGene 0.3544669 0.3544669 0.07998588 -0.4146372 -1.061779  101

gene=fgseaRes[fgseaRes$pathway=="RiskGene",]$leadingEdge[[1]]
g=DotPlot(ADMicTmp,features=gene[1:20],group.by="SubType")
t=ggplot(g$data, aes(factor(features.plot,levels=unique(features.plot)),id, size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_colour_viridis_c()+
labs(color="avg.exp.scaled",size="pct.exp",x="",y="",title="")+
theme_bw()+
theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1)) 
pdf("Graph/Part3/GSEAPlotOfRiskGeneDetailInC9.pdf",height=3,width=5)
print(t)
dev.off()

