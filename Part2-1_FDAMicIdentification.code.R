library(monocle3)
library(dplyr)
library(Seurat)
library(presto)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)

setwd("/Projects/deng/Alzheimer/syn18485175")
MahysDataset=readRDS("MathysDataset.code.R") 
ADMic=subset(MahysDataset,ident=13)
data <- GetAssayData(object = ADMic, slot = "counts")
celldata <- as.data.frame(ADMic@meta.data)
genedata <- as.data.frame(x = row.names(ADMic), row.names = row.names(ADMic))
colnames(genedata) <- "gene_short_name"
Mic.cds <- new_cell_data_set(data, cell_metadata = celldata, gene_metadata = genedata)

#Mic.cds <- as.cell_data_set(ADMic)
#Mic.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(ADMic[["RNA"]])
colData(Mic.cds)
Mic.cds <- preprocess_cds(Mic.cds, num_dim = 20)
plot_pc_variance_explained(Mic.cds)
Mic.cds <- reduce_dimension(Mic.cds,preprocess_method = 'PCA')
Mic.cds = cluster_cells(cds = Mic.cds, resolution=1e-2,reduction_method ='UMAP')
plot_cells(Mic.cds,group_label_size = 5,cell_size = 0.5)
colData(Mic.cds)$MicCluster <- as.character(clusters(Mic.cds))
saveRDS(Mic.cds,file="ADMicMono.rds")

ADMic$MicCluster=as.character(clusters(Mic.cds))
Idents(ADMic)=ADMic$MicCluster
ADMicExpr <- AverageExpression(ADMic)[["RNA"]]
pdf("Graph/Part2/clusterCorShip.pdf",width=6,height=5)
plot(hclust(as.dist(1-cor(ADMicExpr)),method="ward.D2"),hang=-1)
dev.off()

#-----------Marker genes for each cluster--------------
ident=c(4,1,2,3,8,6,5,9,7,10,11)
ADMic@active.ident <- factor(ADMic@active.ident,levels=ident)
ADMic=subset(ADMic,ident=c(1:10))
ADMic.markers <- FindAllMarkers(ADMic, only.pos = TRUE, min.pct = 0.25)
top5 <- ADMic.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
pdf("Graph/Part2/MarkerGeneHeatmap.pdf",height=7,width=6)
DoHeatmap(ADMic, features = top5$gene) + NoLegend()
dev.off()

markerGene=c("F13A1","MRC1","SEPP1","P2RY12","CX3CR1","FCGR3A","VSIG4","CD14","RPS19","SPP1","C1QB","ARMC9","PLEKHA7","RPS24","FTH1","FTL")
pdf("Graph/Part2/MarkerGeneVlnPlot.pdf",height=9,width=6)
VlnPlot(ADMic,features=markerGene,ncol=1,pt.size=0)&
theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())
dev.off()

p=VlnPlot(ADMic,features=markerGene,ncol=1,pt.size=0)
g=ggplot_build(p)

markerGene=c("F13A1","MRC1","SEPP1","P2RY12","CX3CR1","FCGR3A","VSIG4","CD14","RPS19","SPP1","C1QB","ARMC9","PLEKHA7","RPS24","FTH1","FTL")
pdf("Graph/Part2/Mathys_MarkerGeneVlnPlot.pdf",height=9,width=6)
scCustomize::Stacked_VlnPlot(seurat_object = ADMic, features = markerGene, colors_use=unique(g$data[[1]]$fill),x_lab_rotate = TRUE, plot_spacing = 0.05)&
theme(axis.text.x=element_blank(),axis.title.x=element_blank())
dev.off()

#-----------cell type identification for each cluster--------------
new.cluster.ids <- c("BAMic", "HomMic", "HomMic", "HomMic", "HomMic", "HomMic", "NorARMic", "ARMic","FDAMic","DysMic","OliMic")
names(new.cluster.ids) <- levels(ADMic)
ADMic <- RenameIdents(ADMic, new.cluster.ids)
colData(Mic.cds)$Type <- Idents(ADMic)
colData(Mic.cds)$Type <- factor(colData(Mic.cds)$Type,levels=c("BAMic","HomMic","NorARMic","ARMic","FDAMic","DysMic","OliMic"))
colors_list <- c("forestgreen","NavajoWhite", "LightBlue", "SlateBlue", "Violet", "gold", "gray")
tiff("Graph/Part2/MonoMicType.tiff",width=500,height=400)
plot_cells(Mic.cds, color_cells_by="Type",label_cell_groups = FALSE,group_label_size = 5,cell_size =1.5)+
scale_color_manual(values = colors_list)+
theme(legend.position="none",axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())&
theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

saveRDS(ADMic,file="ADMic.rds")

#-----------Female AD number for each cluster--------------
ident=c(4,1,2,3,8,6,5,9,7,10,11)
Idents(ADMic)=ADMic$MicCluster
ADMic@active.ident <- factor(ADMic@active.ident,levels=ident)
ADMic=subset(ADMic,ident=c(1:10))
metaInfo=ADMic@meta.data
count=data.frame(table(metaInfo$MicCluster,metaInfo$Gender,metaInfo$Statues))
colnames(count)=c("Cluster","Gender","Group","Number")
count$Category=paste0(count$Gender,count$Group,sep="")
total=data.frame(table(paste0(metaInfo$Gender,metaInfo$Statues,sep="")))
colnames(total)=c("Category","Total")
graphData=merge(count,total,by="Category")
graphData$Ratio=graphData$Number/graphData$Total
graphData$Group=factor(graphData$Group,levels=c("Control","Alzheimer"))
graphData$Cluster=factor(graphData$Cluster,levels=c(4,1,2,3,8,6,5,9,7,10,11))
graphData$Gender=factor(graphData$Gender,levels=c("Male","Female"))
t=ggplot(graphData, aes(Cluster, Ratio, fill=Group)) +
  geom_bar(stat="identity",position="fill") +
  scale_fill_manual(values = c("SlateBlue", "Coral"))+
  scale_y_continuous(expand=c(0,0))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0))+
  labs(size="",x="",y="Cell ratio",title="")+
  facet_grid(Gender~.)
pdf("Graph/Part2/Mic_Monocle_GroupDis.pdf",width=6,height=4)
print(t)
dev.off()

pg <- ggplot_build(t)
data=graphData
data$ScaledValue=pg$data[[1]]$y
data=data[order(data$Cluster,data$Gender),]
write.table(data,file="Graph/Part2/Mic_Monocle_GroupDis.txt",sep="\t",row.names=F,quote=F)

#-----------sample distribution for each cluster to remove the sample bias--------------
ADMic=subset(ADMic,ident=c(1:10))
sampleInCluster=as.data.frame.matrix(table(ADMic$ProjectID,ADMic$MicCluster))
ident=as.character(c(4,1,2,3,8,6,5,9,7,10))
sampleInCluster=sampleInCluster[,ident]
colnames(sampleInCluster)=paste0("C",colnames(sampleInCluster),sep="")
SampleInfo=data.frame(ADMic$ProjectID,ADMic$Statues,ADMic$Gender)
SampleInfo=unique(SampleInfo)
colnames(SampleInfo)=c("ID","Statues","Sex")
rownames(SampleInfo)=SampleInfo$ID
SampleInfo$ID=NULL
ann_colors = list(
    Statues = c(Alzheimer="Coral", Control="SkyBlue"),
    Sex = c(Male="RoyalBlue", Female="Violet")
)
pdf("Graph/Part2/sampleDisInCluster.pdf",height=6,width=5)
pheatmap(sampleInCluster,clustering_method="ward.D2",border=NA,annotation_colors=ann_colors,annotation_row=SampleInfo,color = colorRampPalette(c("white", "white", "firebrick3"))(50),cluster_col=F,scale="row",show_rownames=F)
dev.off()

#-----------differentiation potential for each cluster with different method-------------------------
setwd("D:/Alzheimer/syn18485175")
data=read.table("D:/Alzheimer/syn18485175/cluster13/monocle/Differentiation/SCENTResult.txt",header=T)
ident=c("C4","C1","C2","C3","C8","C6","C5","C9","C7","C10","C11")
Cluster=factor(data$CellType,levels=ident)
t=ggplot(data, aes(x=Cluster, y=Entropy, color=Cluster))+
geom_violin(scale="width")+labs(title="",x="", y = "")+
theme_classic()+
theme(legend.position="none")+
geom_dotplot(binaxis='y', stackdir='centerwhole',  position=position_dodge(0.8),dotsize=0.5,binwidth=0.003)+
theme(axis.title.y = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),angle=90),axis.text.y = element_text(size=rel(1.0)))
pdf("Graph/Part2/DiffPotential8CytoSCENT.pdf",height=3,width=6)
print(t)
dev.off()

data=read.table("D:/Alzheimer/syn18485175/cluster13/monocle/Differentiation/CytoTRACE_results.txt",header=T,row.names=1)
Cluster=factor(data$Phenotype,levels=ident)
t=ggplot(data, aes(x=Cluster, y=CytoTRACE, color=Cluster))+
geom_boxplot()+labs(title="",x="", y = "")+
theme_classic()+
theme(legend.position="none")+
geom_dotplot(binaxis='y', stackdir='centerwhole', position=position_dodge(0.8),dotsize=0.5,binwidth=0.02)+
theme(axis.title.y = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),angle=90),axis.text.y = element_text(size=rel(1.0)))
pdf("Graph/Part2/DiffPotential8CytoTRACE.pdf",height=3,width=6)
print(t)
dev.off()


#-----------visualization of pesudotime and other phenotypes-------------------------
Mic.cds=readRDS("ADMicMono.rds")

metaData=colData(Mic.cds)
rootCells <- rownames(metaData[metaData$MicCluster%in%c(3,9),])
Mic.cds <- order_cells(Mic.cds,root_cells =rootCells)

tiff("Manuscript/Microglia/Graph/Part2/MonoMicPseudotime.tiff",width=500,height=400)
plot_cells(Mic.cds, color_cells_by = "pseudotime",group_label_size = 5,cell_size =1.5)+
theme(legend.position="none",axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())&
theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

ADMic=readRDS("ADMic.rds")
ADMicTmp=subset(ADMic,ident=c(1:11))
new.cluster.ids <- c("BAMic", "HomMic", "HomMic", "HomMic", "HomMic", "HomMic", "NorARMic", "ARMic","FDAMic","DysMic"," ")
names(new.cluster.ids) <- levels(ADMicTmp)
ADMicTmp <- RenameIdents(ADMicTmp, new.cluster.ids)
colData(Mic.cds)$Type <- Idents(ADMicTmp)
colData(Mic.cds)$Type <- factor(colData(Mic.cds)$Type,levels=c("BAMic","HomMic","NorARMic","ARMic","FDAMic","DysMic"," "))
colors_list <- c("forestgreen","NavajoWhite", "LightBlue", "SlateBlue", "Violet", "gold", "gray")
tiff("Manuscript/Microglia/Graph/Part2/MonoMicType.tiff",width=500,height=400)
plot_cells(Mic.cds, color_cells_by="Type",label_cell_groups = FALSE,group_label_size = 5,cell_size =1.5)+
scale_color_manual(values = colors_list)+
theme(legend.position="none",axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())&
theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

tiff("Manuscript/Microglia/Graph/Part2/MonoMicCluster.tiff",width=500,height=400)
plot_cells(Mic.cds, color_cells_by="cluster",label_cell_groups = FALSE, show_trajectory_graph = TRUE,group_label_size = 5,cell_size =1.5)&
theme(legend.position="none",axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())&
theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
tiff("Manuscript/Microglia/Graph/Part2/MonoMicStatues.tiff",width=500,height=400)
plot_cells(Mic.cds, color_cells_by="Statues",label_cell_groups = FALSE,group_label_size = 5,cell_size =1.5)+
scale_color_manual(values = c("Coral", "SlateBlue"))+
theme(legend.position="none",axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())&
theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
tiff("Manuscript/Microglia/Graph/Part2/MonoMicGender.tiff",width=500,height=400)
plot_cells(Mic.cds, color_cells_by="Gender",label_cell_groups = FALSE,group_label_size = 5,cell_size =1.5)+
scale_color_manual(values = c("Violet","RoyalBlue"))+
theme(legend.position="none",axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())&
theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

tiff("Manuscript/Microglia/Graph/Part2/MonoMicNFT.tiff",width=500,height=400)
plot_cells(Mic.cds, color_cells_by="NFT",label_cell_groups = FALSE,group_label_size = 5,cell_size =1.5)+
#scale_color_gradient(low="gray90",high = "firebrick3")+
scale_color_distiller(palette="Spectral")+
theme(legend.position="none",axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())&
theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
tiff("Manuscript/Microglia/Graph/Part2/MonoMicPlaq.tiff",width=500,height=400)
plot_cells(Mic.cds, color_cells_by="Plaq",label_cell_groups = FALSE,group_label_size = 5,cell_size =1.5)+
#scale_color_gradient(low="gray90",high = "firebrick3")+
scale_color_distiller(palette="Spectral")+
theme(legend.position="none",axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())&
theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
tiff("Manuscript/Microglia/Graph/Part2/MonoMicCogdx.tiff",width=500,height=400)
plot_cells(Mic.cds, color_cells_by="Cogdx",label_cell_groups = FALSE,group_label_size = 5,cell_size =1.5)+
#scale_color_gradient(low="gray90",high = "firebrick3")+
scale_color_distiller(palette="Spectral")+
theme(legend.position="none",axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())&
theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
tiff("Manuscript/Microglia/Graph/Part2/MonoMicPlaqLegend.tiff",width=600,height=400)
plot_cells(Mic.cds, color_cells_by="Plaq",label_cell_groups = FALSE,group_label_size = 5,cell_size =1.5)+
#scale_color_gradient(low="gray90",high = "firebrick3")+
scale_color_distiller(palette="Spectral")+
theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())&
theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

##Mean value for each cluster
ADMic=readRDS("ADMic.rds")
ident=c("4","1","2","3","8","6","5","9","7","10","11")
Idents(ADMic)=ADMic$MicCluster
ADMicTmp=subset(ADMic,ident=c(1:11))
metaInfo=ADMicTmp@meta.data
#not so good by boxplot, too many outliers
Mean=metaInfo %>% group_by(MicCluster) %>% summarize(Mean = mean(abs(cogn_global_lv), na.rm=TRUE)) %>% data.frame()
Cluster=factor(Mean$MicCluster,levels=rev(ident))
t=ggplot(Mean, aes(Cluster, y=Mean,group=1))+
geom_area(stat = "identity", fill="BlueViolet",alpha=0.6)+labs(title="",x="", y = "")+
theme_classic()+
sapply(Cluster, function(xint) 
    geom_vline(aes(xintercept = xint),lwd=0.6,lty=4,col="white")
    )+
theme(legend.position="none")+coord_flip()+
theme(axis.title.y = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0)))
pdf("Manuscript/Microglia/Graph/Part2/cogn_global_lv_Dist.pdf",height=4,width=2)
print(t)
dev.off()

tiff("Manuscript/Microglia/Graph/Part2/MonoMicAge.tiff",width=500,height=400)
plot_cells(Mic.cds, color_cells_by="Age",label_cell_groups = FALSE,group_label_size = 5,cell_size =1.5)+
theme(legend.position="none",axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())&
theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
##Mean value for Age in each cluster
Mean=metaInfo %>% group_by(MicCluster) %>% summarize(Mean = mean(Age, na.rm=TRUE)) %>% data.frame()
Cluster=factor(Mean$MicCluster,levels=rev(ident))
t=ggplot(Mean, aes(Cluster, y=abs(Mean),group=1))+
geom_area(stat = "identity", fill="BlueViolet",alpha=0.6)+labs(title="",x="", y = "")+
theme_classic()+
sapply(Cluster, function(xint) 
    geom_vline(aes(xintercept = xint),lwd=0.6,lty=4,col="white")
    )+
theme(legend.position="none")+coord_flip(ylim=c(70,95))+
theme(axis.title.y = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0)))
pdf("Manuscript/Microglia/Graph/Part2/Age_Dist.pdf",height=4,width=2)
print(t)
dev.off()

NFT=metaInfo %>% group_by(MicCluster) %>% summarize(Mean = mean(abs(NFT), na.rm=TRUE)) %>% data.frame()
Plaq=metaInfo %>% group_by(MicCluster) %>% summarize(Mean = mean(abs(Plaq), na.rm=TRUE)) %>% data.frame()
cogn_global_lv=metaInfo %>% group_by(MicCluster) %>% summarize(Mean = mean(abs(cogn_global_lv), na.rm=TRUE)) %>% data.frame()
Age=metaInfo %>% group_by(MicCluster) %>% summarize(Mean = mean(Age, na.rm=TRUE)) %>% data.frame()
AverageMetaInfo=data.frame(cluster=NFT$MicCluster,NFT=NFT$Mean,Plaq=Plaq$Mean,cogn_global_lv=cogn_global_lv$Mean,age=Age$Mean)
AverageMetaInfo$cluster=factor(AverageMetaInfo$cluster,levels=ident)
AverageMetaInfo=AverageMetaInfo[order(AverageMetaInfo$cluster),]
write.table(AverageMetaInfo,file="Manuscript/Microglia/Graph/Part2/AverageMetaInfoForEachCluster.txt",row.names=F,quote=F,sep="\t")




############ significant difference between AD and Ctrl after removing FDAMic ############
AD=readRDS("MathysDataset.code.R")
ADMic=readRDS("ADMic.rds")
ADMicRemoveC7=subset(ADMic,MicCluster%in%c(1:6,8:11))
#ADMicRemoveC7=subset(ADMic,MicCluster%in%c(1:11))
C13N=data.frame(table(ADMicRemoveC7$ProjectID))
TotalN=data.frame(table(AD$ProjectID))
all(C13N$Var1==TotalN$Var1) #should be TRUE
cellNumber=cbind(TotalN,C13N[,2])
colnames(cellNumber)=c("ProjectID","TotalNumber","C13Number")
ADInfo=AD@meta.data
sampleInfo=unique(ADInfo[,c("ProjectID","Gender","pathology.group")])
cellNumberInfo=merge(cellNumber,sampleInfo,by="ProjectID")
cellNumberInfo$Ratio=cellNumberInfo$C13Number/cellNumberInfo$TotalNumber
Sex=factor(cellNumberInfo$Gender,levels=c("Male","Female"))
Group=factor(cellNumberInfo$pathology.group,levels=c("no-pathology","early-pathology","late-pathology"))
t=ggplot(cellNumberInfo, aes(x=Group, y=Ratio*100, fill=Sex)) +
  geom_boxplot(position=position_dodge(0.8)) +  
  geom_dotplot(binaxis='y', stackdir='center',  position=position_dodge(0.8))+
  scale_fill_manual(values = c("RoyalBlue","Violet"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0)))
pdf("Manuscript/Microglia/Graph/Part2/cellNumberVariationInBeforeRemoveC7_pathology.group.pdf",width=6,height=6)
print(t)
dev.off()

wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$pathology.group=="late-pathology","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$pathology.group=="late-pathology","Ratio"]) #p-value = 0.1905  #0.5556
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$pathology.group=="early-pathology","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$pathology.group=="early-pathology","Ratio"]) #p-value = 0.7789 #0.6943
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$pathology.group=="no-pathology","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$pathology.group=="no-pathology","Ratio"]) #p-value = 0.03872 #p-value = 0.02842
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$pathology.group=="no-pathology","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$pathology.group=="late-pathology","Ratio"]) #p-value = 0.02967  #p-value = 0.2121
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$pathology.group=="no-pathology","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$pathology.group=="late-pathology","Ratio"]) #p-value = 0.3827  #p-value = 0.3284

#--------------check the DEGs affected by FDAMic--------------#
setwd("/Projects/deng/Alzheimer/syn18485175")
ADMic=readRDS("ADMic.rds")
ADdeg <- wilcoxauc(subset(ADMic,Gender=="Female"), 'Statues')
DEGBRmC7 <- ADdeg %>% dplyr::filter(group == "Alzheimer") %>% arrange(desc(logFC)) %>% dplyr::select(feature, c(logFC,pval))
DEGBRmC7$threshold = as.factor(ifelse(DEGBRmC7$pval < 0.01 & abs(DEGBRmC7$logFC) >= 0.25, ifelse(DEGBRmC7$logFC >= 0.25 ,'Up','Down'),'NoSignificant'))
table(DEGBRmC7$threshold)
#Down NoSignificant            Up
#46         17853            27
write.table(DEGBRmC7[DEGBRmC7$threshold%in%c('Up','Down'),],file="Manuscript/Microglia/Graph/Part2/DEGBRmC7.txt",sep="\t",quote=F,row.names=F)

ADdeg <- wilcoxauc(subset(ADMic,Gender=="Female"& MicCluster %in% c(1:6,8:11)), 'Statues')
DEGARmC7 <- ADdeg %>% dplyr::filter(group == "Alzheimer") %>% arrange(desc(logFC)) %>% dplyr::select(feature, c(logFC,pval))
DEGARmC7$threshold = as.factor(ifelse(DEGARmC7$pval < 0.01 & abs(DEGARmC7$logFC) >= 0.25, ifelse(DEGARmC7$logFC >= 0.25 ,'Up','Down'),'NoSignificant'))
table(DEGARmC7$threshold)
#Down NoSignificant            Up
# 33         17871            22
write.table(DEGARmC7[DEGARmC7$threshold%in%c('Up','Down'),],file="Manuscript/Microglia/Graph/Part2/DEGARmC7.txt",sep="\t",quote=F,row.names=F)

DEGBRmC7$Group="DEGBRmC7"
DEGARmC7$Group="DEGARmC7"
sigGeneB=DEGBRmC7[DEGBRmC7$threshold%in%c("Up","Down"),"feature"] #73
sigGeneA=DEGARmC7[DEGARmC7$threshold%in%c("Up","Down"),"feature"] #55
sigGeneList=unique(c(sigGeneB,sigGeneA)) #combine 84 genes
sigGeneListInfo=rbind(DEGBRmC7[DEGBRmC7$feature%in%sigGeneList,],DEGARmC7[DEGARmC7$feature%in%sigGeneList,])

FocusedGene=c("LINGO1","RASGEF1B","SRGN","RPS19","ACSL1","HIF1A","HSPA1A","APOE","SPP1","CD81","ARMC9")

t=ggplot(data = sigGeneListInfo, aes(x = logFC, y = -log10(pval), colour=Group, label =feature )) +
  geom_point(alpha=0.4, size=3.5) +
  theme_bw() + 
  scale_color_manual(values=c("blue","red")) +
  xlim(c(-1, 1)) +
  geom_vline(xintercept=c(-0.25,0.25),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)",title="Differential Expressed Gene") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank()) +  
    geom_text_repel(
    data = subset(sigGeneListInfo, sigGeneListInfo$feature%in%FocusedGene),
    aes(label = feature),
    size = 3,
    max.overlaps=20,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
pdf("Manuscript/Microglia/Graph/Part2/DEGbetBARmC7.pdf",height=5,width=7)
print(t)
dev.off()

ADMic_Female_BeforeRmC7=subset(ADMic,Gender=="Female")
ADMic_Female_AfterRmC7=subset(ADMic_Female_BeforeRmC7,MicCluster %in% c(1:6,8:11))
ADMic_Female_BeforeRmC7$Group="BeforeRmC7"
ADMic_Female_AfterRmC7$Group="AfterRmC7"
ADMic_Female_Tmp=merge(ADMic_Female_BeforeRmC7,ADMic_Female_AfterRmC7)

ADMic_Female_Tmp$Statues=factor(ADMic_Female_Tmp$Statues,levels=c("Control","Alzheimer"))
ADMic_Female_Tmp$Group=factor(ADMic_Female_Tmp$Group,levels=c("BeforeRmC7","AfterRmC7"))
pdf("Manuscript/Microglia/Graph/Part2/sigGeneAfterRmC7.pdf",width=14,height=4)
VlnPlot(ADMic_Female_Tmp,features=c("RPS19","APOE","SPP1"),pt.size=0,ncol=3,cols = c("SlateBlue","Coral"),group.by="Group",split.by="Statues",split.plot = TRUE)&
theme(legend.position = 'right')&theme(axis.title.x=element_blank(),axis.title.y=element_blank())
dev.off()
pdf("Manuscript/Microglia/Graph/Part2/sigGeneAfterRmC7_1.pdf",width=14,height=4)
VlnPlot(ADMic_Female_Tmp,features=c("RPS19","APOE","SPP1"),ncol=3,cols = c("SlateBlue","Coral"),group.by="Group",split.by="Statues")&
theme(legend.position = 'right')&theme(axis.title.x=element_blank(),axis.title.y=element_blank())
dev.off()

DEGARmC7[DEGARmC7$feature%in%c("RPS19","APOE","SPP1"),]
feature      logFC       pval     threshold
APOE 0.21172482 0.05588783 NoSignificant
RPS19 0.14530093 0.17257429 NoSignificant
SPP1 0.02436593 0.93730403 NoSignificant


#Reviewer2: to answer "Is this to be expected when any number of cells are removed"
table(DEGBRmC7$threshold)
#Down NoSignificant            Up
#46         17853            27

table(DEGARmC7$threshold)
#Down NoSignificant            Up
# 33         17871            22

table(ADMic$SubType)
#BAMic   HomMic NorARMic    ARMic   FDAMic   DysMic   OpcMic
#158      998      154      136      148      117       24
#randomly removed 148 genes from the total cells
RandomCell=ncol(ADMic)-148

PermutationNumber=1000
DEGNumber=matrix(data = NA, nrow = PermutationNumber, ncol = 2, byrow = FALSE,dimnames = list(c(1:PermutationNumber),c("UpDEGNumber","DownDEGNumber")))
for(i in c(1:PermutationNumber)){
rm(.Random.seed)
ADMicDownSampling=ADMic[, sample(colnames(ADMic), size = RandomCell, replace=F)]
#table(ADMicDownSampling$SubType)
ADRandomlydeg <- wilcoxauc(ADMicDownSampling, 'Statues')
DEGARandomlyRm <- ADRandomlydeg %>% dplyr::filter(group == "Alzheimer") %>% arrange(desc(logFC)) %>% dplyr::select(feature, c(logFC,pval))
DEGARandomlyRm$Pattern = as.factor(ifelse(DEGARandomlyRm$pval < 0.01 & abs(DEGARandomlyRm$logFC) >= 0.25, ifelse(DEGARandomlyRm$logFC >= 0.25 ,'Up','Down'),'NoSignificant'))
DEGNumber[i,1]=nrow(DEGARandomlyRm[DEGARandomlyRm$Pattern%in%c("Up"),])
DEGNumber[i,2]=nrow(DEGARandomlyRm[DEGARandomlyRm$Pattern%in%c("Down"),])
write.table(DEGARandomlyRm[DEGARandomlyRm$Pattern%in%c('Up','Down'),],file=paste0("Manuscript/Microglia/Graph/Part2/RamdomDownSample/",i,"_pval.txt",sep=""),sep="\t",quote=F,row.names=F)
}

DEGNumber.df=data.frame(DEGNumber)
mean(DEGNumber.df$UpDEGNumber)
#19.615
mean(DEGNumber.df$DownDEGNumber)
#11.494
write.table(DEGNumber.df,file="Manuscript/Microglia/Graph/Part2/DEGNumber.df.txt",sep="\t",quote=F)


DEGBRmC7$Group="DEGBRmC7"
DEGARandomlyRm$Group="DEGARandomlyRm"
sigGeneB=DEGBRmC7[DEGBRmC7$threshold%in%c("Up","Down"),"feature"] #73
sigGeneA=DEGARandomlyRm[DEGARandomlyRm$threshold%in%c("Up","Down"),"feature"] #55
sigGeneList=unique(c(sigGeneB,sigGeneA)) #combine 84 genes
sigGeneListInfo=rbind(DEGBRmC7[DEGBRmC7$feature%in%sigGeneList,],DEGARandomlyRm[DEGARandomlyRm$feature%in%sigGeneList,])
FocusedGene=c("LINGO1","RASGEF1B","SRGN","RPS19","ACSL1","HIF1A","HSPA1A","APOE","SPP1","CD81","ARMC9")
t=ggplot(data = sigGeneListInfo, aes(x = logFC, y = -log10(pval), colour=Group, label =feature )) +
  geom_point(alpha=0.4, size=3.5) +
  theme_bw() + 
  scale_color_manual(values=c("blue","red")) +
  xlim(c(-1, 1)) +
  geom_vline(xintercept=c(-0.25,0.25),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)",title="Differential Expressed Gene") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank()) +  
    geom_text_repel(
    data = subset(sigGeneListInfo, sigGeneListInfo$feature%in%FocusedGene),
    aes(label = feature),
    size = 3,
    max.overlaps=20,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
pdf("Manuscript/Microglia/Graph/Part2/DEGbetBARandomlyRmC7.pdf",height=5,width=7)
print(t)
dev.off()


#--------------check the expression of female specific adDEGs across subtypes------------------#
GlobalDEG=read.table("/Projects/deng/Alzheimer/syn18485175/Manuscript/Microglia/Graph/Part1/DEGbetADCtrlInFemale.txt",header=T,row.names=1)
GlobalDEG=GlobalDEG[GlobalDEG$threshold%in%c("Up","Down"),c("logFC","pval","threshold")]
GeneInfo=GlobalDEG
GeneInfo$logFC=NULL
GeneInfo$pval=-log10(GeneInfo$pval)
Idents(ADMic)=ADMic$MicCluster
ADMic=subset(ADMic,idents=c(1:10))
ADMicExpr <- AverageExpression(ADMic)[["RNA"]]
GlobalDEGExpr=ADMicExpr[rownames(GlobalDEG),]
pdf("Manuscript/Microglia/Graph/Part2/ExprOfGlobalDEGInFemaleSubtypeGraph.pdf",width=9,height=5)
pheatmap(t(GlobalDEGExpr),clustering_method="ward.D2",border=NA,scale="column",show_rownames=T,annotation_col=GeneInfo,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

#--------------Up and down regulated genes in FDAMic------------------#
ClusterMarker <- wilcoxauc(ADMic, 'MicCluster')
clusterCell<- ClusterMarker %>% dplyr::filter(group == 7) %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC,pval)
ranks<- deframe(clusterCell)
GlobalDEGInC7=clusterCell[clusterCell$feature%in%rownames(GlobalDEG),]
hs_data=GlobalDEGInC7
#hs_data$threshold = as.factor(ifelse(hs_data$hs_data<0,ifelse(hs_data$p_val<0.01&abs(hs_data$avg_log2FC)>0.25,"SigDown","Down"),ifelse(hs_data$p_val<0.01&abs(hs_data$avg_log2FC)>0.25,"SigUp","Up")))
hs_data$threshold = as.factor(ifelse(hs_data$feature%in%rownames(GlobalDEG[GlobalDEG$threshold=="Up",]),"UpInFemaleAD","DnInFemaleAD"))
table(hs_data$threshold)
#DnInFemaleAD UpInFemaleAD
#          46           27
hs_data$ID=hs_data$feature
t=ggplot(data = hs_data, aes(x = logFC, y = -log10(pval), colour=threshold, label =ID )) +
  geom_point(alpha=0.6, size=3.5) +
  theme_bw() + 
  scale_color_manual(values=c("blue", "red")) +
  geom_vline(xintercept=c(-0.25,0.25),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)",title="Differential Expressed Gene") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank()) +  
    ggrepel::geom_text_repel(
    data = subset(hs_data, hs_data$pval < 0.01 & abs(hs_data$logFC) >= 0.25),
    aes(label = ID),
    size = 3,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
pdf("Manuscript/Microglia/Graph/Part2/ExprOfGlobalDEGInFemaleADInARMic_New.pdf",height=5,width=7)
print(t)
dev.off()

hs_data$sigGene = ifelse(hs_data$pval < 0.01 & abs(hs_data$logFC) >0.25, ifelse(hs_data$logFC > 0.25 ,'Up','Down'),'NoSignificant')
table(hs_data$sigGene)
#Down NoSignificant            Up
#21           294            31
write.table(hs_data[hs_data$sigGene%in%c("Down","Up"),],file="Manuscript/Microglia/Graph/Part2/ExprOfGlobalDEGInFemaleADInARMic.txt",sep="\t",quote=F,row.names=F)

#---------GSEA before and after removeing FDAMic---------------------
setwd("/Projects/deng/Alzheimer/syn18485175/")
ADMic=readRDS("ADMic.rds")

m_df<- msigdbr(species = "Homo sapiens", category = "C2",subcategory="CP:KEGG")  #CP:KEGG
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
ADMicSex=subset(ADMic,Gender=="Female")
ADMicDEG <- wilcoxauc(ADMicSex, 'Statues')
clusterCell<- ADMicDEG %>% dplyr::filter(group == "Alzheimer") %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
ranks<- deframe(clusterCell)
fgseaRes<- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0)
fwrite(fgseaRes, file=paste0("/Projects/deng/Alzheimer/syn18485175/cluster13/monocle/GSEA/RemoveC7/FemaleBeforeRmC7.txt",sep=""), sep="\t", sep2=c("", " ", ""))

ADMicSubtype=subset(ADMic,idents=c(1:6,8:10))
ADMicSex=subset(ADMicSubtype,Gender=="Female")
ADMicDEG <- wilcoxauc(ADMicSex, 'Statues')
clusterCell<- ADMicDEG %>% dplyr::filter(group == "Alzheimer") %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
ranks<- deframe(clusterCell)
fgseaRes<- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0)
fwrite(fgseaRes, file=paste0("/Projects/deng/Alzheimer/syn18485175/cluster13/monocle/GSEA/RemoveC7/FemaleAfterRmC7.txt",sep=""), sep="\t", sep2=c("", " ", ""))


FAC7=read.table("D:/Alzheimer/syn18485175/cluster13/monocle/GSEA/RemoveC7/FemaleAfterRmC7.txt",header=T,sep="\t")
FAC7$Group="FAC7"
FBC7=read.table("D:/Alzheimer/syn18485175/cluster13/monocle/GSEA/RemoveC7/FemaleBeforeRmC7.txt",header=T,sep="\t")
FBC7$Group="FBC7"
MBC7=read.table("D:/Alzheimer/syn18485175/cluster13/monocle/GSEA/RemoveC7/MaleBeforeRmC7.txt",header=T,sep="\t")
MBC7$Group="MBC7"
data=rbind(FAC7,FBC7,MBC7)
dataT=data[data$pval<0.01,]
dataT=data[data$pathway %in% unique(dataT$pathway),]
dataT=dataT[order(dataT$pval),]
t=ggplot(dataT,aes(Group,factor(pathway,levels=rev(unique(pathway))),size=-1*log(pval),color=NES,fill=NES))+geom_point()+
scale_color_gradient2(low="blue",mid="white",high = "red",midpoint = 0)+
scale_fill_gradient2(low="blue",mid="white",high = "red",midpoint = 0)+
theme_bw()+
theme(axis.title.x=element_blank(),axis.title.y=element_blank())
pdf("D:/Alzheimer/syn18485175/Manuscript/Microglia/Graph/Part2/GSEABeAfRmC7.pdf",width=6.5,height=5)
print(t)
dev.off()


m_df<- msigdbr(species = "Homo sapiens", category = "C2",subcategory="CP:KEGG")  #CP:KEGG
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

ADMic=readRDS("ADMic.rds")
#ADMictmp=subset(ADMic,pathology.group%in%c("early-pathology","no-pathology"))
ADMictmp=subset(ADMic,Gender=="Female")
ADMictmp=subset(ADMictmp,Gender=="Female"&MicCluster%in%c(1:6,8:11))
ADdeg <- wilcoxauc(ADMictmp, 'Statues')
clusterCell<- ADdeg %>% dplyr::filter(group == "Alzheimer") %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
ranks<- deframe(clusterCell)
FemaleInADCtrl_fgseaRes<- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0)

targetPathways=c("KEGG_INSULIN_SIGNALING_PATHWAY","KEGG_MAPK_SIGNALING_PATHWAY","KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS","KEGG_FC_EPSILON_RI_SIGNALING_PATHWAY")
FemaleInADCtrl_fgseaRes[FemaleInADCtrl_fgseaRes$pathway%in%targetPathways,]

pdf("Manuscript/Microglia/Graph/Part1/LauAD_KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSISInADFemale.pdf",height=2.5,width=4)
plotEnrichment(fgsea_sets[["KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS"]],ranks) 
dev.off()
#pva=0.07607192 NES=-1.322086
pdf("Manuscript/Microglia/Graph/Part1/Female_LauAD_KEGG_INSULIN_SIGNALING_PATHWAY.pdf",height=2.5,width=4)
plotEnrichment(fgsea_sets[["KEGG_INSULIN_SIGNALING_PATHWAY"]],ranks) 
dev.off()
#pval=0.33, NES=-1.0719
pdf("Manuscript/Microglia/Graph/Part1/Female_LauAD_KEGG_MAPK_SIGNALING_PATHWAY.pdf",height=2.5,width=4)
plotEnrichment(fgsea_sets[["KEGG_MAPK_SIGNALING_PATHWAY"]],ranks) 
dev.off()
pdf("Manuscript/Microglia/Graph/Part1/LauAD_KEGG_FC_EPSILON_RI_SIGNALING_PATHWAYInADFemale.pdf",height=2.5,width=4)
plotEnrichment(fgsea_sets[["KEGG_FC_EPSILON_RI_SIGNALING_PATHWAY"]],ranks) 
dev.off()


pdf("Manuscript/Microglia/Graph/Part2/PHAGOCYTOSIS_early_normal_Female_RmC7.pdf",height=2.5,width=4)
plotEnrichment(fgsea_sets[["KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS"]],ranks) 
dev.off()
FemaleInADCtrl_fgseaRes[FemaleInADCtrl_fgseaRes$pathway=="KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS",]
#Before RmC7 pval=0.02474371, NES=-1.41117
#After RmC7 pval=0.06972477, NES=-1.295288


ADMictmp=subset(ADMic,pathology.group%in%c("late-pathology","no-pathology"))
#ADMictmp=subset(ADMictmp,Gender=="Male")
ADMictmp=subset(ADMictmp,Gender=="Male"&MicCluster%in%c(1:6,8:11))
#ADdeg <- wilcoxauc(ADMictmp, 'pathology.group')
clusterCell<- ADdeg %>% dplyr::filter(group == "late-pathology") %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
ranks<- deframe(clusterCell)
FemaleInADCtrl_fgseaRes<- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0)
pdf("Manuscript/Microglia/Graph/Part2/PHAGOCYTOSIS_late_normal_Male_RmC7.pdf",height=2.5,width=4)
plotEnrichment(fgsea_sets[["KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS"]],ranks) 
dev.off()
FemaleInADCtrl_fgseaRes[FemaleInADCtrl_fgseaRes$pathway=="KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS",]
#Before RmC7 pval=9.243178e-06, NES=-1.767006
#After RmC7 pval=0.0001546292, -1.725549
fwrite(FemaleInADCtrl_fgseaRes, file=paste0("Manuscript/Microglia/Graph/Part2/RmC7/late_normal_KEGG_MaleAfter_RmC7_.txt",sep=""), sep="\t", sep2=c("", " ", ""))

