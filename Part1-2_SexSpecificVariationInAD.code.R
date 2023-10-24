setwd("/Projects/deng/Alzheimer/syn18485175/Manuscript/Microglia")
library(Seurat)
library(ggplot2)
library(viridis)
library(scCustomize)
library(dplyr)
library(ggrepel)
library(presto)
library(tibble)
library(msigdbr)
library(fgsea)
library(pheatmap)
library(UpSetR)
library(clusterProfiler)

###########recluster the dataset from Mathys###################
#Output from Part1-1_MathyDataset.code.R
MahysDataset=readRDS("MathysDataset.code.R") 

##########cell type identification and cluster defination#########################
tiff("Graph/Part1/cellypeMarker.tiff",width=400,height=800)
FeaturePlot(MahysDataset,c("NRGN","GMahysDataset1","MBP","GFAP","VCAN","CSF1R","FLT1"),cols =c("lightgrey","red"),ncol=2)&NoLegend()&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
MahysDataset$cellType=factor(MahysDataset$OldCellType,levels=c("Ex","In","Oli","Opc","Ast","Mic","End","Per"))
tiff("Graph/Part1/MahysDatasetcellType.tiff",width=450,height=400)
DimPlot(MahysDataset,group.by="cellType")&NoLegend()&NoAxes()&theme(plot.title=element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
tiff("Graph/Part1/MahysDatasetcellCluster.tiff",width=450,height=400)
DimPlot(MahysDataset)&NoLegend()&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

MahysDataset$Statues=factor(MahysDataset$Statues,levels=c("Control","Alzheimer"))
tiff("Graph/Part1/CellNumberDisInAllCellType.tiff",width=500,height=250)
DimPlot(MahysDataset,group.by="Gender",split.by="Statues",cols = c("Violet","RoyalBlue"))&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
tiff("Graph/Part1/CellNumberDisInAllCellType_Pathology.Group.tiff",width=750,height=250)
DimPlot(MahysDataset,group.by="Gender",split.by="pathology.group",cols = c("Violet","RoyalBlue"))&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
tiff("Graph/Part1/CellNumberDisInAllCellType_scCustom.tiff",width=550,height=250)
DimPlot_scCustom(seurat_object = MahysDataset,group.by="Gender",split.by="Statues",colors_use = c("Violet","RoyalBlue"))&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
tiff("Graph/Part1/CellNumberDisInAllCellType_Pathology_scCustom.tiff",width=800,height=250)
DimPlot_scCustom(seurat_object = MahysDataset,group.by="Gender",split.by="pathology.group",colors_use = c("Violet","RoyalBlue"))&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

tiff("Graph/Part1/MahysDatasetcellTypeLabel.tiff",width=450,height=400)
DimPlot(MahysDataset,group.by="cellType",label=T,label.size = 6)&NoLegend()&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
tiff("Graph/Part1/MahysDatasetcellClusterLabel.tiff",width=450,height=400)
DimPlot(MahysDataset,label=T,label.size = 6)&NoLegend()&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

##########cell number variation for each cluster (microglia showed significant changed in female AD)####################
MahysDatasettmp=subset(MahysDataset,ident=c(0:21)) #cluster 22 was removed due to the small numer
tmp=paste0(MahysDatasettmp$seurat_clusters,"_",MahysDatasettmp$Gender,"_",MahysDatasettmp$Statues,sep="")
tmp=as.matrix(table(tmp))
TmpInfo=as.matrix(do.call(rbind, strsplit(as.character(rownames(tmp)),'_')))
result=data.frame("Group"=paste0(TmpInfo[,2],"_",TmpInfo[,3],sep=""),"Cluster"=paste0("C",TmpInfo[,1],sep=""),"Gender"=TmpInfo[,2],"Statues"=TmpInfo[,3],"Number"=tmp[,1])
dim(result)
i=1 
while(i<88){t=chisq.test(matrix(c(result[i,"Number"],result[(i+1),"Number"],result[(i+2),"Number"],result[(i+3),"Number"]),ncol=2))$p.value;i=i+4;print(t)}

ClusterOrder=factor(result$Cluster,levels=c("C1","C15","C14","C4","C3","C7","C17","C11","C18","C5","C8","C9","C16","C12","C19","C2","C0","C10","C6","C13","C21","C20"))
GenderList=factor(result$Gender,levels=c("Male","Female"))
t=ggplot(result, aes(ClusterOrder, Number, fill=GenderList)) +
  geom_bar(stat="identity",position="fill") +
  scale_fill_manual(values = c("RoyalBlue","Violet"))+
  scale_y_continuous(expand=c(0,0))+
  theme(axis.text.x = element_text(angle = 90))+
  labs(size="",x="",y="Cell ratio",title="")+
  facet_grid(factor(Statues,levels=c("Control","Alzheimer"))~.)
pdf("Graph/Part1/cellNumberVariation.pdf",width=10,height=4)
print(t)
dev.off()

ClusterOrder=c("C1","C15","C14","C4","C3","C7","C17","C11","C18","C5","C8","C9","C16","C12","C19","C2","C0","C10","C6","C13","C21","C20")
MahysDatasettmp$Group=paste0(MahysDatasettmp$Gender,"_",MahysDatasettmp$Statues,sep="")
tmp=as.data.frame.matrix(table(MahysDatasettmp$seurat_clusters,MahysDatasettmp$Group))
rownames(tmp)=paste0("C",rownames(tmp))
tmp=tmp[ClusterOrder,]
pdf("Graph/Part1/cellNumberCountInEachCondition.pdf",width=3)
pheatmap(tmp,cluster_row=F,clustering_method="ward.D2",legend = T,number_format="%d",cluster_col=F,display_numbers = T,color =colorRampPalette(brewer.pal(6, "Set3"))(100))
dev.off()

#---------Pvalue for each cluster-------------------
Pvalue=array()
i=1 
while(i<88){
	t=chisq.test(matrix(c(result[i,"Number"],result[(i+1),"Number"],result[(i+2),"Number"],result[(i+3),"Number"]),ncol=2))$p.value;
	index=(i+3)/4;
	Pvalue[index]=t;
	i=i+4;
	print(t)
}
PvalueData=data.frame(cellCluster=unique(result$Cluster),LogPvalue=-log10(Pvalue))
PvalueData$cellCluster=factor(PvalueData$cellCluster,levels=rev(c("C1","C15","C14","C4","C3","C7","C17","C11","C18","C5","C8","C9","C16","C12","C19","C2","C0","C10","C6","C13","C21","C20")))
PvalueDataTmp=PvalueData[order(PvalueData$cellCluster,decreasing=T),]
write.table(PvalueDataTmp,file="Graph/Part1/cellNumberVariationPvalue.txt",sep="\t",quote=F,row.names=F)

PvalueData$Type=rev(c("Ex","Ex","Ex","Ex","Ex","Ex","Ex","Ex","Ex","Ex","In","In","In","In","In","Oli","Oli","Opc","Ast","Mic","Mic","End"))
cellType=factor(PvalueData$Type,levels=rev(unique(PvalueData$Type)))
t=ggplot(PvalueData, aes(x=cellCluster, y=LogPvalue, fill=cellType)) + geom_bar(stat="identity")+theme_minimal()+ #data from cellClusterInfo.xlsx
  #theme(axis.text.x=element_blank())+
  coord_flip()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(size="",x="",y="-log10 (Pvalue)",title="")
pdf("Graph/Part1/cellNumberVariationPvalueH.pdf",width=4,height=10)
print(t)
dev.off()


C13=subset(MahysDataset,ident=13)
C13$Statues=factor(C13$Statues,levels=c("Control","Alzheimer"))
tiff("Graph/Part1/CellNumberDisInCluster13.tiff",width=500,height=250)
DimPlot(C13,group.by="Gender",split.by="Statues",cols = c("Violet","RoyalBlue"))&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

C13$pathology.group=factor(C13$pathology.group,levels=c("no-pathology","early-pathology","late-pathology"))
tiff("Graph/Part1/CellNumberDisInCluster13_Pathology.Group.tiff",width=750,height=250)
DimPlot(C13,group.by="Gender",split.by="pathology.group",cols = c("Violet","RoyalBlue"))&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

#Trends across LOAD patholgy
ClusterOrder=c("C1","C15","C14","C4","C3","C7","C17","C11","C18","C5","C8","C9","C16","C12","C19","C2","C0","C10","C6","C13","C21","C20")
Female=subset(MahysDataset,Gender=="Female")
FemaleCount=data.frame(table(paste0(Female$seurat_clusters,"_",Female$pathology.group)))
colnames(FemaleCount)=c("Group","n")
TotalCount=data.frame(table(paste0(MahysDataset$seurat_clusters,"_",MahysDataset$pathology.group)))
colnames(TotalCount)=c("Group","N")
GraphData=merge(FemaleCount,TotalCount,by="Group")
GraphData$Ratio=GraphData$n/GraphData$N
TmpInfo=as.matrix(do.call(rbind, strsplit(as.character(GraphData$Group),'_')))
GraphData$cluster=paste0("C",TmpInfo[,1])
GraphData$pathology.group=TmpInfo[,2]
GraphData=GraphData[GraphData$cluster%in%c(paste0("C",c(0:21))),]
GraphData$cluster=factor(GraphData$cluster,levels=ClusterOrder)
GraphData$pathology.group=factor(GraphData$pathology.group,levels=c("no-pathology","early-pathology","late-pathology"))
scale_this <- function(x){x - min(x)}
scaled_data <- GraphData %>% group_by(cluster) %>%mutate(ScaledRatio = scale_this(Ratio))
g=ggplot(scaled_data,aes(x=pathology.group, y=ScaledRatio,group=1))+
  geom_area(fill="Violet",outline.type="full") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.text.y =element_blank(),axis.ticks = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(panel.spacing=unit(0.1, "lines"))+
  sapply(levels(scaled_data$pathology.group), function(xint) 
    geom_vline(aes(xintercept = xint),lwd=0.6,lty=4,col="black")
    )+
  facet_grid(cluster~.)+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),axis.line =element_blank())
pdf("Graph/Part1/FemaleCellsRation_Pathology.Group.pdf",width=4,height=14)
print(g)
dev.off()

scaled_data$cluster=factor(scaled_data$cluster,levels=ClusterOrder)
scaled_data$pathology.group=factor(scaled_data$pathology.group,levels=c("no-pathology","early-pathology","late-pathology"))
scaled_data=scaled_data[order(scaled_data$cluster,scaled_data$pathology.group),]
write.table(scaled_data,file="Graph/Part1/FemaleCellsRation_Pathology.Group.txt",sep="\t",quote=F,row.names=F)



C13N=data.frame(table(C13$ProjectID))
TotalN=data.frame(table(MahysDataset$ProjectID))
all(C13N$Var1==TotalN$Var1) #should be TRUE
cellNumber=cbind(TotalN,C13N[,2])
colnames(cellNumber)=c("ProjectID","TotalNumber","C13Number")
MahysDatasetInfo=MahysDataset@meta.data
sampleInfo=unique(MahysDatasetInfo[,c("ProjectID","Gender","Statues")])
cellNumberInfo=merge(cellNumber,sampleInfo,by="ProjectID")
cellNumberInfo$Ratio=cellNumberInfo$C13Number/cellNumberInfo$TotalNumber
Sex=factor(cellNumberInfo$Gender,levels=c("Male","Female"))
Group=factor(cellNumberInfo$Statues,levels=c("Control","Alzheimer"))
t=ggplot(cellNumberInfo, aes(x=Group, y=Ratio, fill=Sex)) +
  geom_boxplot(position=position_dodge(0.8)) +  
  geom_dotplot(binaxis='y', stackdir='center',  position=position_dodge(0.8))+
  scale_fill_manual(values = c("RoyalBlue","Violet"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0)))
pdf("Graph/Part1/cellNumberVariationInMicFromMathyData.pdf",width=6,height=6)
print(t)
dev.off()
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$Statues=="Alzheimer","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$Statues=="Control","Ratio"]) #p-value = 0.08873
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$Statues=="Alzheimer","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$Statues=="Control","Ratio"]) #p-value = 0.16
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$Statues=="Alzheimer","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$Statues=="Alzheimer","Ratio"]) #p-value = 0.197
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$Statues=="Control","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$Statues=="Control","Ratio"]) #p-value = 0.03872


## Sex dependent microglia ratio across various pathology and congnitive impairment 
C13N=data.frame(table(C13$ProjectID))
TotalN=data.frame(table(MahysDataset$ProjectID))
all(C13N$Var1==TotalN$Var1) #should be TRUE
cellNumber=cbind(TotalN,C13N[,2])
colnames(cellNumber)=c("ProjectID","TotalNumber","C13Number")
MahysDatasetInfo=MahysDataset@meta.data
MahysDatasetInfo=MahysDatasetInfo[MahysDatasetInfo$Cogdx%in%c(1,2,4),]
sampleInfo=unique(MahysDatasetInfo[,c("ProjectID","Gender","Cogdx")])
cellNumberInfo=merge(cellNumber,sampleInfo,by="ProjectID")
cellNumberInfo$Ratio=cellNumberInfo$C13Number/cellNumberInfo$TotalNumber
Sex=factor(cellNumberInfo$Gender,levels=c("Male","Female"))
Group=factor(cellNumberInfo$Cogdx,levels=c(1,2,4))
t=ggplot(cellNumberInfo, aes(x=Group, y=Ratio*100, fill=Sex)) +
  geom_boxplot(position=position_dodge(0.8)) +  
  geom_dotplot(binaxis='y', stackdir='center', size=1, position=position_dodge(0.8))+
  scale_fill_manual(values = c("RoyalBlue","Violet"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0)))
pdf("Graph/Part1/cellNumberVariationInMicFromMathyDataGroup8Cogdx.pdf",width=6,height=6)
print(t)
dev.off()
summary(cellNumberInfo$Ratio)
cellNumberInfo[cellNumberInfo$Ratio==max(cellNumberInfo$Ratio),]
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$Cogdx=="1","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$Cogdx=="1","Ratio"]) #p-value = 0.4376
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$Cogdx=="2","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$Cogdx=="2","Ratio"]) #p-value = 0.4857
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$Cogdx=="4","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$Cogdx=="4","Ratio"]) #p-value = 0.08452
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$Cogdx=="4","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$Cogdx=="1","Ratio"]) #p-value = 0.4497
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$Cogdx=="4","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$Cogdx=="1","Ratio"]) #p-value = 0.4497


C13N=data.frame(table(C13$ProjectID))
TotalN=data.frame(table(MahysDataset$ProjectID))
all(C13N$Var1==TotalN$Var1) #should be TRUE
cellNumber=cbind(TotalN,C13N[,2])
colnames(cellNumber)=c("ProjectID","TotalNumber","C13Number")
MahysDatasetInfo=MahysDataset@meta.data
sampleInfo=unique(MahysDatasetInfo[,c("ProjectID","Gender","BraakStage")])
cellNumberInfo=merge(cellNumber,sampleInfo,by="ProjectID")
cellNumberInfo$Ratio=cellNumberInfo$C13Number/cellNumberInfo$TotalNumber
cellNumberInfo$Group=ifelse(cellNumberInfo$BraakStage%in%c("Stage_1","Stage_2"),"Stage_1&2",ifelse(cellNumberInfo$BraakStage%in%c("Stage_3","Stage_4"),"Stage_3&4","Stage_5&6"))
Sex=factor(cellNumberInfo$Gender,levels=c("Male","Female"))
t=ggplot(cellNumberInfo, aes(x=Group, y=Ratio*100, fill=Sex)) +
  geom_boxplot(position=position_dodge(0.8)) +  
  geom_dotplot(binaxis='y', stackdir='center', size=1, position=position_dodge(0.8))+
  scale_fill_manual(values = c("RoyalBlue","Violet"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0)))
pdf("Graph/Part1/cellNumberVariationInMicFromMathyDataGroup8BraakStage.pdf",width=6,height=6)
print(t)
dev.off()

wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$Group=="Stage_1&2","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$Group=="Stage_1&2","Ratio"]) #p-value = 0.6667
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$Group=="Stage_3&4","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$Group=="Stage_3&4","Ratio"]) #p-value = 0.223
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$Group=="Stage_5&6","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$Group=="Stage_5&6","Ratio"]) #p-value = 0.4747
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$Group=="Stage_5&6","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$Group=="Stage_1&2","Ratio"]) #p-value = 0.9371
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$Group=="Stage_5&6","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$Group=="Stage_1&2","Ratio"]) #p-value = 0.2593


C13N=data.frame(table(C13$ProjectID))
TotalN=data.frame(table(MahysDataset$ProjectID))
all(C13N$Var1==TotalN$Var1) #should be TRUE
cellNumber=cbind(TotalN,C13N[,2])
colnames(cellNumber)=c("ProjectID","TotalNumber","C13Number")
MahysDatasetInfo=MahysDataset@meta.data
sampleInfo=unique(MahysDatasetInfo[,c("ProjectID","Gender","Ceradsc")])
cellNumberInfo=merge(cellNumber,sampleInfo,by="ProjectID")
cellNumberInfo$Ratio=cellNumberInfo$C13Number/cellNumberInfo$TotalNumber
cellNumberInfo=cellNumberInfo[cellNumberInfo$Ceradsc%in%c(1,2,4),]
cellNumberInfo$Group=factor(cellNumberInfo$Ceradsc,levels=c(4,2,1))
Sex=factor(cellNumberInfo$Gender,levels=c("Male","Female"))
t=ggplot(cellNumberInfo, aes(x=Group, y=Ratio*100, fill=Sex)) +
  geom_boxplot(position=position_dodge(0.8)) +  
  geom_dotplot(binaxis='y', stackdir='center', size=1, position=position_dodge(0.8))+
  scale_fill_manual(values = c("RoyalBlue","Violet"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0)))
pdf("Graph/Part1/cellNumberVariationInMicFromMathyDataGroup8Ceradsc.pdf",width=6,height=6)
print(t)
dev.off()

wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$Group=="4","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$Group=="4","Ratio"]) #p-value = 0.01577
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$Group=="2","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$Group=="2","Ratio"]) #p-value = 0.1905
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$Group=="1","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$Group=="1","Ratio"]) #p-value = 0.07024
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$Group=="1","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$Group=="4","Ratio"]) #p-value = 0.02948
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$Group=="1","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$Group=="4","Ratio"]) #p-value = 0.1042


C13N=data.frame(table(C13$ProjectID))
TotalN=data.frame(table(MahysDataset$ProjectID))
all(C13N$Var1==TotalN$Var1) #should be TRUE
cellNumber=cbind(TotalN,C13N[,2])
colnames(cellNumber)=c("ProjectID","TotalNumber","C13Number")
MahysDatasetInfo=MahysDataset@meta.data
sampleInfo=unique(MahysDatasetInfo[,c("ProjectID","Gender","pathology.group")])
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
pdf("Graph/Part1/cellNumberVariationInMicFromMathyDataGroup8pathology.group.pdf",width=6,height=6)
print(t)
dev.off()
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$pathology.group=="late-pathology","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$pathology.group=="late-pathology","Ratio"]) #p-value = 0.1905 #p-value = 0.02857
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$pathology.group=="early-pathology","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$pathology.group=="early-pathology","Ratio"]) #p-value = 0.7789
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$pathology.group=="no-pathology","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Male"&cellNumberInfo$pathology.group=="no-pathology","Ratio"]) #p-value = 0.03872
wilcox.test(cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$pathology.group=="no-pathology","Ratio"],cellNumberInfo[cellNumberInfo$Gender=="Female"&cellNumberInfo$pathology.group=="late-pathology","Ratio"]) #p-value = 0.03872



############validate the cell number from ROSMAP bulk RNASeq##############
setwd("D:/Alzheimer/syn18485175/Manuscript/Microglia")
#downloaded from syn3157322/ROSMAP_biospecimen_metadata.csv, ROSMAP_clinical.csv, ROSMAP_clinical_codebook.pdf, ROSMAP_assay_RNAseq_metadata.csv
phenotype=read.table("D:/Alzheimer/Syn3388564/Phenotype.txt",header=T,row.names=1,sep="\t") #the data download from synapse database
trait=phenotype[,c(13,15,16,17)]
t=pheatmap(trait,clustering_method="ward.D2",show_rownames=F)
a=data.frame(cutree(t$tree_row,k=2))
anno=data.frame(Group=paste0("G",a[,1],sep=""))
rownames(anno)=rownames(a)
all(rownames(anno)==rownames(phenotype))
anno$Gender=phenotype$msex
ann_colors = list(
    Group = c(G1="Coral", G2="SlateBlue"),
    Gender=c(Female="Violet",Male="RoyalBlue")
)
t=pheatmap(trait,clustering_method="ward.D2",annotation_row=anno,annotation_colors = ann_colors,show_rownames=F,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
pdf("Graph/Part1/TraitGroup2DefineMahysDatasetInROSMAP.pdf",height=5,width=4.5)
print(t)
dev.off()
all(rownames(anno)==rownames(phenotype))
phenotype$StatusCluster=anno$Group
phenotype$Status=ifelse(phenotype$StatusCluster=="G1","LOAD","ND")
write.table(phenotype,file="Graph/Part1/PhenotypeStatus.txt",sep="\t",quote=F)


###############DEG between AD and Ctrl in female and male using bulk RNASeq ###############
library(limma)
#D:\Alzheimer\Syn3388564\AS\Manuscript\graph\graph.code.R
setwd("D:/Alzheimer/syn18485175/Manuscript/Microglia")
phenotype=read.table("Graph/Part1/PhenotypeStatus.txt",header=T,row.names=1,sep="\t")
nrow(phenotype)
#613
#downloaded from Syn3388564/ROSMAP_RNAseq_FPKM_gene.tsv
GeneExprByFPKMFromROSMAP=read.table("D:/Alzheimer/Syn3388564/geneExpr.txt",header=TRUE,row.names=1,sep="\t",check.names=F)
sampleList=intersect(rownames(phenotype),colnames(GeneExprByFPKMFromROSMAP))
length(sampleList)
#613
GeneExprByFPKMFromROSMAP=GeneExprByFPKMFromROSMAP[,sampleList]
all(rownames(phenotype)==colnames(GeneExprByFPKMFromROSMAP))
#TRUE
range(GeneExprByFPKMFromROSMAP)

data <- log2(GeneExprByFPKMFromROSMAP+1)
range(data)

phenotype.Male=phenotype[phenotype$msex=="Male",c("Status","pmi","age_at_visit_max")]
phenotype.Male[is.na(phenotype.Male)] <- 0
phenotype.Male$Status=factor(phenotype.Male$Status,levels=c("ND","LOAD"))
data.Male=data[,rownames(phenotype.Male)]
design.Male=model.matrix(~Status+pmi+age_at_visit_max,data =phenotype.Male)
fit.Male=lmFit(data.Male,design.Male)
fit.Male=eBayes(fit.Male)
top.Male=topTable(fit.Male, coef = "StatusLOAD",number=Inf,p.value=1, lfc=0)
#these soft cutoffs were accepted to idnetifed a number of DEGs in Male
top.Male$Pattern=ifelse(top.Male$P.Value < 0.01 & abs(top.Male$logFC) >= 0, ifelse(top.Male$logFC >= 0 ,'Up','Down'),'NoSig')
SigDEG.Male=top.Male[top.Male$Pattern%in%c("Up","Down"),]
table(top.Male$Pattern)
#Down NoSig    Up 
#441 49501   232
write.table(SigDEG.Male,file="Graph/Part1/BlukDEG/SigDEGInMaleByBulkROSMAP.txt",sep="\t",quote=F)


phenotype.Female=phenotype[phenotype$msex=="Female",c("Status","pmi","age_at_visit_max")]
phenotype.Female[is.na(phenotype.Female)] <- 0
phenotype.Female$Status=factor(phenotype.Female$Status,levels=c("ND","LOAD"))
data.Female=data[,rownames(phenotype.Female)]
design.Female=model.matrix(~Status+pmi+age_at_visit_max,data =phenotype.Female)
fit.Female=lmFit(data.Female,design.Female)
fit.Female=eBayes(fit.Female)
top.Female=topTable(fit.Female, coef = "StatusLOAD",number=Inf,p.value=1, lfc=0)
top.Female$Pattern=ifelse(top.Female$P.Value < 0.01 & abs(top.Female$logFC) >= 0, ifelse(top.Female$logFC >= 0,'Up','Down'),'NoSig')
SigDEG.Female=top.Female[top.Female$Pattern%in%c("Up","Down"),]
table(top.Female$Pattern)
#Down NoSig    Up 
#4910 43811  1453
write.table(SigDEG.Female,file="Graph/Part1/BlukDEG/SigDEGInFemaleByBulkROSMAP.txt",sep="\t",quote=F)


#functional enrichment analysis of the upregulated genes (focused on gliogenesis)
gene.df <- bitr(rownames(SigDEG.Female[SigDEG.Female$Pattern=="Up",]), fromType="SYMBOL",toType="ENTREZID",OrgDb = "org.Hs.eg.db") 
dim(gene.df)
#1267
go <- enrichGO(gene = gene.df$ENTREZID, OrgDb = "org.Hs.eg.db", ont="BP",readable =T)
go=data.frame(go)
tmp=go[1:5,]
tmp$Description=factor(tmp$Description,levels=rev(tmp$Description))
t=ggplot(tmp, aes(Description, -log10(p.adjust), fill=-log10(p.adjust))) +
  geom_bar(stat="identity") +
  scale_fill_viridis_c()+
  theme_bw()+
  coord_flip(ylim = c(5, 10))+
  labs(size="",x="",y="-log(p.ajdust)",title="")
pdf("Graph/Part1/BlukDEG/DEGByLimmaEnrichGO8ROSMAP.pdf",height=3,width=6)
print(t)
dev.off()
write.table(go,file="Graph/Part1/BlukDEG/DEGEnrichGO8ROSMAP.txt",sep="\t",quote=F)


SigDEG.Male=read.table("Graph/Part1/BlukDEG/SigDEGInMaleByBulkROSMAP.txt",header=T,row.names=1)
SigDEG.Female=read.table("Graph/Part1/BlukDEG/SigDEGInFemaleByBulkROSMAP.txt",header=T,row.names=1)
SigDEG.Male.list=split(rownames(SigDEG.Male),paste0(SigDEG.Male$Pattern,"InMale",sep=""))
SigDEG.Female.list=split(rownames(SigDEG.Female),paste0(SigDEG.Female$Pattern,"InFemale",sep=""))
DEGCount=c(SigDEG.Male.list,SigDEG.Female.list)
t=upset(fromList(DEGCount),sets =rev(c("DownInFemale","UpInFemale","DownInMale","UpInMale")),keep.order = T,sets.bar.color =c("#EE82EE","#FF0099","#4169E1","#CCCCFF"))
pdf("Graph/Part1/BlukDEG/DEGEnrichGO8ROSMAP_Overlap.pdf",height=4,width=6)
print(t)
dev.off()



##check the overlap between different DEG methods.
FemaleDEGByLimma=read.table("Graph/Part1/BlukDEG/SigDEGInFemaleByBulkROSMAP.txt",header=T,row.names=1)
FemaleDEGByLimma$Pattern=ifelse(FemaleDEGByLimma$logFC>0,"UpByLimma_Female","DownByLimma_Female")

#t test was used in our first version and replaced by Limma, but similiar DEGs were obtained between these two methods
FemaleDEGByTtest=read.table("D:/Alzheimer/Syn3388564/AS/FemaleDEG.txt",header=T,row.names=1)
FemaleDEGByTtest=FemaleDEGByTtest[FemaleDEGByTtest$pvalue<0.01&abs(FemaleDEGByTtest$logFC)>log2(1.2),]
FemaleDEGByTtest$Pattern=ifelse(FemaleDEGByTtest$logFC>log2(1.2),"UpByTtest_Female","DnByTtest_Female")
DEGList=c(split(rownames(FemaleDEGByLimma),FemaleDEGByLimma$Pattern),split(rownames(FemaleDEGByTtest),FemaleDEGByTtest$Pattern))
pdf("Graph/Part1/BlukDEG/OverlapBetweenTtestAndLimma_Female.pdf",width=6,height=4)
upset(fromList(DEGList),main.bar.color=c(rep("black",4),"Blue","Red"),matrix.color="black",sets.bar.color=c(rep("RoyalBlue",2),rep("Orange",2)))
dev.off()
write.table(FemaleDEGByTtest,"Graph/Part1/BlukDEG/FemaleDEGByTtest.txt",sep="\t",quote=F)
write.table(FemaleDEGByLimma,"Graph/Part1/BlukDEG/FemaleDEGByLimma.txt",sep="\t",quote=F)

##Enrichment analysis showed simliar result to single DEGs identification method
length(intersect(DEGList$UpByTtest_Female,DEGList$UpByLimma_Female)) #757
gene.df <- bitr(intersect(DEGList$UpByTtest_Female,DEGList$UpByLimma_Female), fromType="SYMBOL",toType="ENTREZID",OrgDb = "org.Hs.eg.db") 
dim(gene.df)
#626
go <- enrichGO(gene = gene.df$ENTREZID, OrgDb = "org.Hs.eg.db", ont="BP",readable =T)
go=data.frame(go)
tmp=go[1:30,]
tmp$Description=factor(tmp$Description,levels=rev(tmp$Description))
t=ggplot(tmp, aes(Description, -log10(pvalue), fill=-log10(pvalue))) +
  geom_bar(stat="identity") +
  scale_fill_viridis_c()+
  theme_bw()+
  coord_flip(ylim = c(5, 10))+
  labs(size="",x="",y="-log(pvalue)",title="")
pdf("Graph/Part1/BlukDEG/SharedUpGeneFemaleEnrichment.pdf",height=8,width=8)
print(t)
dev.off()
write.table(go,file="Graph/Part1/BlukDEG/SharedUpGeneFemaleEnrichment.txt",sep="\t",quote=F)


MaleDEGByLimma=read.table("Graph/Part1/BlukDEG/SigDEGInMaleByBulkROSMAP.txt",header=T,row.names=1)
MaleDEGByLimma$Pattern=ifelse(MaleDEGByLimma$logFC>0,"UpByLimma_Male","DownByLimma_Male")
MaleDEGByTtest=read.table("D:/Alzheimer/Syn3388564/AS/MaleDEG.txt",header=T,row.names=1)
MaleDEGByTtest=MaleDEGByTtest[MaleDEGByTtest$pvalue<0.01&abs(MaleDEGByTtest$logFC)>log2(1.2),]
MaleDEGByTtest$Pattern=ifelse(MaleDEGByTtest$logFC>log2(1.2),"UpByTtest_Male","DnByTtest_Male")
DEGList=c(split(rownames(MaleDEGByLimma),MaleDEGByLimma$Pattern),split(rownames(MaleDEGByTtest),MaleDEGByTtest$Pattern))
pdf("Graph/Part1/BlukDEG/OverlapBetweenTtestAndLimma_Male.pdf",width=6,height=4)
upset(fromList(DEGList),main.bar.color=c(rep("black",4),"Blue","Red"),matrix.color="black",sets.bar.color=c(rep("RoyalBlue",2),rep("Orange",2)))
dev.off()
write.table(MaleDEGByTtest,"Graph/Part1/BlukDEG/MaleDEGByTtest.txt",sep="\t",quote=F)
write.table(MaleDEGByLimma,"Graph/Part1/BlukDEG/MaleDEGByLimma.txt",sep="\t",quote=F)



##Downsampling the female sample to verify the affection of sample number on DEGs identification (As more female samples to male) ###################
table(phenotype$msex,phenotype$Status)
#         LOAD  ND
#  Female  173 219
#  Male     88 133
PermutationNumber=1000
DEGNumber=matrix(data = NA, nrow = PermutationNumber, ncol = 2, byrow = FALSE,dimnames = list(c(1:PermutationNumber),c("UpDEGNumber","DownDEGNumber")))
for(i in c(1:PermutationNumber)){
LOADDownSample=sample(rownames(phenotype[phenotype$msex=="Female"&phenotype$Status=="LOAD",]),88,replace=F) #randomly select 88 samples from 173 LOAD sample
NDDownSample=sample(rownames(phenotype[phenotype$msex=="Female"&phenotype$Status=="ND",]),133,replace=F) #randomly select 133 samples from 219 LOAD sample
DownSampleList=c(LOADDownSample,NDDownSample)
length(DownSampleList)
#221
phenotype.Female.DownSampling=phenotype[DownSampleList,c("Status","pmi","age_at_visit_max")]
phenotype.Female.DownSampling[is.na(phenotype.Female.DownSampling)] <- 0
phenotype.Female.DownSampling$Status=factor(phenotype.Female.DownSampling$Status,levels=c("ND","LOAD"))
data.Female.DownSampling=data[,rownames(phenotype.Female.DownSampling)]
all(rownames(phenotype.Female.DownSampling)==colnames(data.Female.DownSampling))
#TRUE
design.Female.DownSampling=model.matrix(~Status+pmi+age_at_visit_max,data =phenotype.Female.DownSampling)
fit.Female.DownSampling=lmFit(data.Female.DownSampling,design.Female.DownSampling)
fit.Female.DownSampling=eBayes(fit.Female.DownSampling)
top.Female.DownSampling=topTable(fit.Female.DownSampling, coef = "StatusLOAD",number=Inf,p.value=1, lfc=0)
top.Female.DownSampling$Pattern=ifelse(top.Female.DownSampling$P.Value < 0.01 & abs(top.Female.DownSampling$logFC) >= 0, ifelse(top.Female.DownSampling$logFC >= 0,'Up','Down'),'NoSig')
DEGNumber[i,1]=nrow(top.Female.DownSampling[top.Female.DownSampling$Pattern%in%c("Up"),])
DEGNumber[i,2]=nrow(top.Female.DownSampling[top.Female.DownSampling$Pattern%in%c("Down"),])
}
DEGNumber.df=data.frame(DEGNumber)
mean(DEGNumber.df$UpDEGNumber)
#888.957
mean(DEGNumber.df$DownDEGNumber)
#2606.542
write.table(DEGNumber.df,file="D:/Alzheimer/syn18485175/Manuscript/Microglia/Graph/Part1/BlukDEG/DEGNumberPermutation.txt",sep="\t",quote=F)

#Another 20 random permutation to check the enrichment analysis
PermutationNumber=20
for(i in c(1:PermutationNumber)){
LOADDownSample=sample(rownames(phenotype[phenotype$msex=="Female"&phenotype$Status=="LOAD",]),88,replace=F) #randomly select 88 samples from 173 LOAD sample
NDDownSample=sample(rownames(phenotype[phenotype$msex=="Female"&phenotype$Status=="ND",]),133,replace=F) #randomly select 133 samples from 219 LOAD sample
DownSampleList=c(LOADDownSample,NDDownSample)
length(DownSampleList)
#221
phenotype.Female.DownSampling=phenotype[DownSampleList,c("Status","pmi","age_at_visit_max")]
phenotype.Female.DownSampling[is.na(phenotype.Female.DownSampling)] <- 0
phenotype.Female.DownSampling$Status=factor(phenotype.Female.DownSampling$Status,levels=c("ND","LOAD"))
data.Female.DownSampling=data[,rownames(phenotype.Female.DownSampling)]
all(rownames(phenotype.Female.DownSampling)==colnames(data.Female.DownSampling))
#TRUE
design.Female.DownSampling=model.matrix(~Status+pmi+age_at_visit_max,data =phenotype.Female.DownSampling)
fit.Female.DownSampling=lmFit(data.Female.DownSampling,design.Female.DownSampling)
fit.Female.DownSampling=eBayes(fit.Female.DownSampling)
top.Female.DownSampling=topTable(fit.Female.DownSampling, coef = "StatusLOAD",number=Inf,p.value=1, lfc=0)
top.Female.DownSampling$Pattern=ifelse(top.Female.DownSampling$P.Value < 0.01 & abs(top.Female.DownSampling$logFC) >= 0, ifelse(top.Female.DownSampling$logFC >= 0,'Up','Down'),'NoSig')
SigDEG.Female=top.Female.DownSampling
#focused on the upregualted genes to decipher their functional (the downregualted genes should be derived from neuron)
gene.df <- bitr(rownames(SigDEG.Female[SigDEG.Female$Pattern=="Up",]), fromType="SYMBOL",toType="ENTREZID",OrgDb = "org.Hs.eg.db") 
dim(gene.df)
#1267
go <- enrichGO(gene = gene.df$ENTREZID, OrgDb = "org.Hs.eg.db", ont="BP",readable =T)
go=data.frame(go)
go$UpDEGNumber=nrow(top.Female.DownSampling[top.Female.DownSampling$Pattern%in%c("Up"),])
go$DownDEGNumber=nrow(top.Female.DownSampling[top.Female.DownSampling$Pattern%in%c("Down"),])
N=min(dim(go)[1],50)
tmp=go[1:N,]
write.table(tmp,file=paste0("D:/Alzheimer/syn18485175/Manuscript/Microglia/Graph/Part1/BlukDEG/GOBPEnrich/EnrichAt_",i,"_Top50GoBP.txt",sep=""),sep="\t",quote=F)
}



### DEGs identified by T test, replaced by limma using the linear model under suggestions by reviewers
#DEG=read.table("D:/Alzheimer/Syn3388564/AS/FemaleDEG.txt",header=T,row.names=1)
#DEG=DEG[DEG$pvalue<0.01,]
#DownGene=DEG[DEG$logFC < -log2(1.2),] #3523
#UpGene=DEG[DEG$logFC > log2(1.2),] #1120
#DEGCount=list()
#DEGCount$DnInF=rownames(DownGene)
#DEGCount$UpInF=rownames(UpGene)
#DEG=read.table("D:/Alzheimer/Syn3388564/AS/MaleDEG.txt",header=T,row.names=1)
#DEG=DEG[DEG$pvalue<0.01,]
#DownGene=DEG[DEG$logFC < -log2(1.2),] #3523
#UpGene=DEG[DEG$logFC > log2(1.2),] #1120
#DEGCount$DnInM=rownames(DownGene)
#DEGCount$UpInM=rownames(UpGene)
#t=upset(fromList(DEGCount))
#pdf("Graph/Part1/DEGEnrichGO8ROSMAP_Overlap.pdf",height=4,width=6)
#print(t)
#dev.off()
#DEGCount.df=rbind(
#  data.frame(Symbol=DEGCount$DnInF,Group="DownInFemale"),
#  data.frame(Symbol=DEGCount$UpInF,Group="UpInFemale"),
#  data.frame(Symbol=DEGCount$DnInM,Group="DownInMale"),
#  data.frame(Symbol=DEGCount$UpInM,Group="UpInMale")
 # )
#write.table(DEGCount.df,file="Graph/Part1/DEGByROSMAPBulk.txt",sep="\t",quote=F,row.names=F)

#DEG=read.table("D:/Alzheimer/Syn3388564/AS/FemaleDEG.txt",header=T,row.names=1)
#DEG=DEG[DEG$pvalue<0.01,]
#DownGene=DEG[DEG$logFC < -log2(1.2),] #3523
#UpGene=DEG[DEG$logFC > log2(1.2),] #1120
#gene.df <- bitr(rownames(DownGene), fromType="SYMBOL",toType="ENTREZID",OrgDb = "org.Hs.eg.db") #2982
#gene.df <- bitr(rownames(UpGene), fromType="SYMBOL",toType="ENTREZID",OrgDb = "org.Hs.eg.db") #846
#go <- enrichGO(gene = gene.df$ENTREZID, OrgDb = "org.Hs.eg.db", ont="BP",readable =T)
#go=data.frame(go)
#tmp=go[1:6,]
#tmp$Description=factor(tmp$Description,levels=rev(tmp$Description))
#t=ggplot(tmp, aes(Description, -log10(p.adjust), fill=-log10(p.adjust))) +
#  geom_bar(stat="identity") +
#  scale_fill_viridis_c()+
#  theme_bw()+
#  coord_flip(ylim = c(5, 10))+
#  labs(size="",x="",y="-log(p.ajdust)",title="")
#pdf("Graph/Part1/DEGEnrichGO8ROSMAP.pdf",height=3,width=6)
#print(t)
#dev.off()
#write.table(go,file="Graph/Part1/DEGEnrichGO8ROSMAP.txt",sep="\t",quote=F)



############DEGs between Alzheimer and ND in snRNA##############
setwd("/Projects/deng/Alzheimer/syn18485175/Manuscript/Microglia")
MahysDataset=readRDS("/Projects/deng/Alzheimer/syn18485175/MahysDataset.rds") 
#Add addition infor on seurat object
Info=read.table("/Projects/deng/Alzheimer/syn18485175/Phenotype4Cell.txt",header=TRUE,row.names=1,sep="\t")
MicInfo=Info[colnames(MahysDataset),]
all(rownames(MicInfo)==rownames(MahysDataset@meta.data))
MahysDataset$pmi=ifelse(is.na(MicInfo$Pmi),0,MicInfo$Pmi)
MahysDataset$Age=MicInfo$age_death

C13=subset(MahysDataset,ident=13)
Idents(C13)=C13$Statues


##DEGs between LOAD and ND from male
ADdegByWilcox <- wilcoxauc(subset(C13,Gender=="Male"), 'Statues')
ADdegByWilcox <- ADdegByWilcox %>% dplyr::filter(group == "Alzheimer") %>% arrange(desc(logFC)) %>% dplyr::select(feature, c(logFC,pval))
hs_data=ADdegByWilcox
hs_data$threshold = as.factor(ifelse(hs_data$pval < 0.01 & abs(hs_data$logFC) >= 0, ifelse(hs_data$logFC >= 0 ,'Up','Down'),'NoSignificant'))
table(hs_data$threshold)
#ADdegByWilcox
#Down NoSignificant            Up
#30         17895            287
write.table(hs_data,file="Graph/Part1/DEGbetAD_CtrlInMale_wilcox.txt",sep="\t",row.names=F,quote=F)
MaleListByWilcox=split(hs_data$feature,hs_data$threshold)
MaleListByWilcox$NoSignificant=NULL

#we also detect DEGs by MAST method (suggested by reviewers)
ADdegByMAST=FindMarkers(subset(C13,Gender=="Male"),ident.1="Alzheimer",ident.2="Control",test.use = "MAST",latent.vars =c("Age","pmi"),logfc.threshold = 0,min.pct = 0)
hs_data=ADdegByMAST
colnames(hs_data)=c("pval","logFC","pct.1","pct.2","padj") #rename colnames to map wilcox result
hs_data$threshold = as.factor(ifelse(hs_data$pval < 0.01 & abs(hs_data$logFC) >= 0, ifelse(hs_data$logFC >= 0 ,'Up','Down'),'NoSignificant'))
table(hs_data$threshold)
#ADdegByMAST
#Down NoSignificant          Up
#42         17810            74
write.table(hs_data,file="Graph/Part1/DEGbetAD_CtrlInMale_MAST.txt",sep="\t",quote=F)
MaleListByMAST=split(rownames(hs_data),hs_data$threshold)
MaleListByMAST$NoSignificant=NULL

names(MaleListByWilcox)=c("DnIn_MaleAD_Wilcox","UpIn_MaleAD_Wilcox")
names(MaleListByMAST)=c("DnIn_MaleAD_MAST","UpIn_MaleAD_MAST")
MaleDEGList=c(MaleListByWilcox,MaleListByMAST)


#check the overlap between different DGE algorithms (wilcox and MAST) and got the same results from male
t=upset(fromList(MaleDEGList),
    keep.order = T,
    sets = rev(c("DnIn_MaleAD_Wilcox", "UpIn_MaleAD_Wilcox", "DnIn_MaleAD_MAST","UpIn_MaleAD_MAST")),
    main.bar.color="RoyalBlue",matrix.color="RoyalBlue",sets.bar.color="LimeGreen")
pdf("Graph/Part1/MaleDEGListOverlapBetweenWilcoxAndMAST.pdf",height=4,width=5)
print(t)
dev.off()

length(intersect(MaleDEGList$DnIn_MaleAD_Wilcox,MaleDEGList$DnIn_MaleAD_MAST))
#18
length(MaleDEGList$DnIn_MaleAD_Wilcox)
#30
length(MaleDEGList$DnIn_MaleAD_MAST)
#42
1-phyper(18-1,30,17926-30,)
#0

length(intersect(MaleDEGList$UpIn_MaleAD_Wilcox,MaleDEGList$UpIn_MaleAD_MAST))
#70
length(MaleDEGList$UpIn_MaleAD_Wilcox)
#287
length(MaleDEGList$UpIn_MaleAD_MAST)
#74
1-phyper(70-1,287,17926-287,74)
#0

##viocano plots for DEGs
hs_data$ID=hs_data$feature
t=ggplot(data = hs_data, aes(x = logFC, y = -log10(pval), colour=threshold, label =ID )) +
  geom_point(alpha=0.4, size=3.5) +
  theme_bw() + 
  scale_color_manual(values=c("blue", "grey","red")) +
  xlim(c(-1, 1)) + ylim(0,12)+
  geom_vline(xintercept=c(-0.25,0.25),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)",title="Differential Expressed Gene") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank()) +  
    geom_text_repel(
    data = subset(hs_data, hs_data$pval < 0.01),
    aes(label = ID),
    max.overlaps = 10,
    size = 3,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
tiff("Graph/Part1/DEGbetMahysDatasetCtrlInMale_Gene.tiff",height=1000,width=1200)
print(t)
dev.off()


##DEGs between LOAD and ND from female
ADdegByWilcox <- wilcoxauc(subset(C13,Gender=="Female"), 'Statues')
ADdegByWilcox <- ADdegByWilcox %>% dplyr::filter(group == "Alzheimer") %>% arrange(desc(logFC)) %>% dplyr::select(feature, c(logFC,pval))
hs_data=ADdegByWilcox
hs_data$threshold = as.factor(ifelse(hs_data$pval < 0.01 & abs(hs_data$logFC) >= 0, ifelse(hs_data$logFC >= 0 ,'Up','Down'),'NoSignificant'))
table(hs_data$threshold)
#ADdegByWilcox
#Down NoSignificant            Up
#284         17580            62
write.table(hs_data,file="Graph/Part1/DEGbetAD_CtrlInFemale_wilcox.txt",sep="\t",row.names=F,quote=F)
FemaleListByWilcox=split(hs_data$feature,hs_data$threshold)
FemaleListByWilcox$NoSignificant=NULL

ADdegByMAST=FindMarkers(subset(C13,Gender=="Female"),ident.1="Alzheimer",ident.2="Control",test.use = "MAST",latent.vars =c("Age","pmi"),logfc.threshold = 0,min.pct = 0)
hs_data=ADdegByMAST
colnames(hs_data)=c("pval","logFC","pct.1","pct.2","padj") #rename colnames to map wilcox result
hs_data$threshold = as.factor(ifelse(hs_data$pval < 0.01 & abs(hs_data$logFC) >= 0, ifelse(hs_data$logFC >= 0 ,'Up','Down'),'NoSignificant'))
table(hs_data$threshold)
#ADdegByMAST
#Down NoSignificant            Up
#133         17690           103
write.table(hs_data,file="Graph/Part1/DEGbetAD_CtrlInFemale_MAST.txt",sep="\t",quote=F)
FemaleListByMAST=split(rownames(hs_data),hs_data$threshold)
FemaleListByMAST$NoSignificant=NULL

names(FemaleListByWilcox)=c("DnIn_FemaleAD_Wilcox","UpIn_FemaleAD_Wilcox")
names(FemaleListByMAST)=c("DnIn_FemaleAD_MAST","UpIn_FemaleAD_MAST")
FemaleDEGList=c(FemaleListByWilcox,FemaleListByMAST)

##check the overlap between different DGE algorithms (wilcox and MAST) and got the same results from female
t=upset(fromList(FemaleDEGList),
    keep.order = T,
    sets = rev(c("DnIn_FemaleAD_Wilcox", "UpIn_FemaleAD_Wilcox", "DnIn_FemaleAD_MAST","UpIn_FemaleAD_MAST")),
    main.bar.color="RoyalBlue",matrix.color="RoyalBlue",sets.bar.color="LimeGreen")
pdf("Graph/Part1/FemaleDEGListOverlapBetweenWilcoxAndMAST.pdf",height=4,width=5)
print(t)
dev.off()


length(intersect(FemaleDEGList$DnIn_FemaleAD_Wilcox,FemaleDEGList$DnIn_FemaleAD_MAST))
#123
length(FemaleDEGList$DnIn_FemaleAD_Wilcox)
#284
length(FemaleDEGList$DnIn_FemaleAD_MAST)
#133
1-phyper(123-1,284,17926-284,133)
#0

length(intersect(FemaleDEGList$UpIn_FemaleAD_Wilcox,FemaleDEGList$UpIn_FemaleAD_MAST))
#50
length(FemaleDEGList$UpIn_FemaleAD_Wilcox)
#62
length(FemaleDEGList$UpIn_FemaleAD_MAST)
#103
1-phyper(50-1,62,17926-62,103)
#0


hs_data$ID=hs_data$feature
t=ggplot(data = hs_data, aes(x = logFC, y = -log10(pval), colour=threshold, label =ID )) +
  geom_point(alpha=0.4, size=3.5) +
  theme_bw() + 
  scale_color_manual(values=c("blue", "grey","red")) +
  xlim(c(-1, 1)) + ylim(0,12)+
  geom_vline(xintercept=c(-0.25,0.25),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)",title="Differential Expressed Gene") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank()) +  
    geom_text_repel(
    data = subset(hs_data, hs_data$pval < 0.01),
    aes(label = ID),
    max.overlaps = 10,
    size = 3,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
tiff("Graph/Part1/DEGbetMahysDatasetCtrlInFemale_Gene.tiff",height=2000,width=2200)
print(t)
dev.off()



DEGList=c(MaleListByWilcox,FemaleListByWilcox)
library(UpSetR)
t=upset(fromList(DEGList),
    keep.order = T,
    sets = rev(c("DnIn_FemaleAD_Wilcox", "UpIn_FemaleAD_Wilcox", "UpIn_MaleAD_Wilcox","DnIn_MaleAD_Wilcox")),
    main.bar.color="RoyalBlue",matrix.color="RoyalBlue",sets.bar.color="LimeGreen")
pdf("Graph/Part1/DEGListOverlap.pdf",height=4,width=5)
print(t)
dev.off()
intersect(DEGList$UpInMahysDataset_Male,DEGList$DnInMahysDataset_Female) 
"HSPA1A"    "SLC11A1"   "MSR1"      "TG"        "DDIT4"     "HIST2H2BF" "HSPA1B"    "ZFP36"     "JAK3"      "HSDL2"
intersect(DEGList$UpInMahysDataset_Male,DEGList$UpInMahysDataset_Female) 
#"SPP1"    "APOE"    "CACNA1A" "PTPRG"   "RPS2"    "PTMA"    "APOC1"
intersect(DEGList$DnInMahysDataset_Male,DEGList$DnInMahysDataset_Female) 
#"EMID1" "CHN2"

t=rbind(paste0(DEGList$DnInMahysDataset_Female, collapse=","),paste0(DEGList$UpInMahysDataset_Female, collapse=","),paste0(DEGList$UpInMahysDataset_Male, collapse=","),paste0(DEGList$DnInMahysDataset_Male, collapse=","))
rownames(t)=c("DnInMahysDataset_Female", "UpInMahysDataset_Female", "UpInMahysDataset_Male","DnInMahysDataset_Male")
write.table(t,file="Graph/Part1/DEGlist.txt",col.names=F,quote=F,sep="\t")

## example of the specific genes derived from the upset plot: deleted suggested by reviewers
#C13$Statues=factor(C13$Statues,levels=c("Control","Alzheimer"))
#C13$Gender=factor(C13$Gender,levels=c("Male","Female"))
#Idents(C13)=C13$Gender
#pdf("Graph/Part1/DEGCluster13BetGender.pdf",width=10,height=6)
#VlnPlot(C13,features=c("HSPA1A","MSR1","APOE","SPP1"),pt.size=0,ncol=2,cols = c("SlateBlue","Coral"),group.by="Gender",split.by="Statues",split.plot = TRUE)&theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())
#dev.off()

pdf("Graph/Part1/ProliferationCluster13BetGender.pdf",width=8,height=2)
VlnPlot(C13,features=c("CSF1R","CD81"),split.plot = FALSE,pt.size=0,ncol=2,cols = c("SlateBlue","Coral"),group.by="Gender",split.by="Statues")&theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())
dev.off()



#####not show in the manuscript: CSF1R and CD81 between Alzheimer and Ctrl in ROSMAP bulk RNAseq--validation, other genes like MSR1 not good, might by affected by neurons inhibition##############
#GeneExpr=read.table("D:/Alzheimer/Syn3388564/geneExpr.txt",header=TRUE,row.names=1,sep="\t",check.names=F)
#SampleInfo=read.table("D:/Alzheimer/Syn3388564/Phenotype.txt",header=TRUE,row.names=1,sep="\t")
#Samples=intersect(colnames(GeneExpr),rownames(SampleInfo))
#GeneExpr=GeneExpr[rowSums(GeneExpr)>0,Samples]
#SampleInfo=SampleInfo[Samples,]
#all(colnames(GeneExpr)==rownames(SampleInfo))
#BigData=cbind(SampleInfo,t(GeneExpr))
#BigDataCogdx=BigData[BigData$cogdx %in% c(1,2,4,5),]
#BigDataCogdx$msex=factor(BigDataCogdx$msex,levels=c("Male","Female"))
#t=ggplot(BigDataCogdx, aes(x=as.character(cogdx), y=log2(`FCGR2A`+1),color=as.character(cogdx)))+
#geom_boxplot()+labs(title="",x="", y = "")+ #ylim(c(1.5,3.5))+
#theme (strip.text= element_text(size=15)) + 
#theme_bw()+
#stat_compare_means(comparisons =list(c("1","5")))+
#scale_color_manual(values = c("SlateBlue","DarkOrchid","SandyBrown","Coral"))+
#theme(legend.position="none")+ geom_jitter(shape=1, position=position_jitter(0.2))+
#theme(axis.title.y = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0)))+
#facet_grid(~msex)
#pdf("Graph/Part1/CSF1RInROSMAP.pdf",width=8,height=4)
#print(t)
#dev.off()


#-----------GSEA between MahysDataset and Control---------------------
setwd("/Projects/deng/Alzheimer/syn18485175/Manuscript/Microglia")
MahysDataset=readRDS("MahysDataset.rds") #MahysDataset.rds from D:\Alzheimer\syn18485175\Analysis.R
C13=subset(MahysDataset,ident=13)
m_df<- msigdbr(species = "Homo sapiens", category = "C2",subcategory="CP:KEGG")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

MaleInMahysDatasetCtrl <- wilcoxauc(subset(C13,Gender=="Male"), 'Statues')
MahysDatasetdegMahysDataset<- MaleInMahysDatasetCtrl %>% dplyr::filter(group == "Alzheimer") %>% arrange(desc(logFC)) %>% dplyr::select(feature, c(logFC,pval))
ranks<- deframe(MahysDatasetdegMahysDataset)
MaleInMahysDatasetCtrl_fgseaRes<- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0)
FemaleInMahysDatasetCtrl <- wilcoxauc(subset(C13,Gender=="Female"), 'Statues')
MahysDatasetdegMahysDataset<- FemaleInMahysDatasetCtrl %>% dplyr::filter(group == "Alzheimer") %>% arrange(desc(logFC)) %>% dplyr::select(feature, c(logFC,pval))
ranks<- deframe(MahysDatasetdegMahysDataset)
FemaleInMahysDatasetCtrl_fgseaRes<- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0)

##pathway cateogry was manully downloaded from KEGG database
pathwayCategory=read.table("/Projects/deng/Public/GSEA/KEGGPathwayCategory.txt",header=T,sep="\t")

all(MaleInMahysDatasetCtrl_fgseaRes$pathway==FemaleInMahysDatasetCtrl_fgseaRes$pathway) #TRUE
MaleInMahysDatasetCtrl_fgseaRes$Gender="Male"
FemaleInMahysDatasetCtrl_fgseaRes$Gender="Female"
dataT=rbind(MaleInMahysDatasetCtrl_fgseaRes,FemaleInMahysDatasetCtrl_fgseaRes)
sigPath=dataT[dataT$pval<0.01,]
dataT=dataT[dataT$pathway %in% unique(sigPath$pathway),]
dataT=dataT[order(dataT$pval),]
dataT$pathway=tolower(dataT$pathway)
all(dataT$pathway%in%pathwayCategory$Pathway)
dataT$Pathway=dataT$pathway
dataT=merge(dataT,pathwayCategory,by="Pathway")
dataT=dataT[order(dataT$Group,dataT$SubGroup,dataT$pval),]
dataT$PathwayName=factor(dataT$PathwayName,levels=rev(unique(dataT$PathwayName)))
t=ggplot(dataT,aes(Gender,PathwayName,size=-1*log(pval),color=NES,fill=NES))+geom_point()+
scale_color_gradient2(low="blue",mid="white",high = "red",midpoint = 0)+
scale_fill_gradient2(low="blue",mid="white",high = "red",midpoint = 0)+
theme_bw()+
theme(axis.title.x=element_blank(),axis.title.y=element_blank())
pdf("Graph/Part1/KEGGEnrichBetMahysDatasetCtrlSplit8SexTmp.pdf",width=5,height=5)
print(t)
dev.off()
dataT$leadingEdge=as.character(dataT$leadingEdge)
fwrite(dataT, file="Graph/Part1/KEGGEnrichBetMahysDatasetCtrlSplit8Sex.txt", sep="\t", sep2=c(c("","|","")))


#-----------Validate by independent datasets for some pathway---------------------
setwd("/Projects/deng/Alzheimer/syn18485175/Manuscript/Microglia")
LauMahysDatasetMic=readRDS("/Projects/deng/Alzheimer/Lau/CombineMathys/LauMahysDatasetMic.rds")
m_df<- msigdbr(species = "Homo sapiens", category = "C2",subcategory="CP:KEGG")  #CP:KEGG
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
LauMahysDatasetMictmp=subset(LauMahysDatasetMic,Gender=="Female")
MahysDatasetdeg <- wilcoxauc(LauMahysDatasetMictmp, 'Statues')
clusterCell<- MahysDatasetdeg %>% dplyr::filter(group == "Alzheimer") %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
ranks<- deframe(clusterCell)
fgseaRes<- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0)
fwrite(fgseaRes, file=paste0("GSEA/KEGG/LauMic/LauMic_MahysDatasetvsCtrlInFemale.txt",sep=""), sep="\t", sep2=c("", " ", ""))


targetPathways=c("KEGG_INSULIN_SIGNALING_PATHWAY","KEGG_MAPK_SIGNALING_PATHWAY","KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS","KEGG_FC_EPSILON_RI_SIGNALING_PATHWAY")
fgseaRes[fgseaRes$pathway%in%targetPathways,]

pdf("Graph/Part1/LauMahysDataset_KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSISInMahysDatasetFemale.pdf",height=2.5,width=4)
plotEnrichment(fgsea_sets[["KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS"]],ranks) 
dev.off()
#pva=0.07607192 NES=-1.322086
pdf("Graph/Part1/Female_LauMahysDataset_KEGG_INSULIN_SIGNALING_PATHWAY.pdf",height=2.5,width=4)
plotEnrichment(fgsea_sets[["KEGG_INSULIN_SIGNALING_PATHWAY"]],ranks) 
dev.off()
#pval=0.33, NES=-1.0719
pdf("Graph/Part1/Female_LauMahysDataset_KEGG_MAPK_SIGNALING_PATHWAY.pdf",height=2.5,width=4)
plotEnrichment(fgsea_sets[["KEGG_MAPK_SIGNALING_PATHWAY"]],ranks) 
dev.off()
pdf("Graph/Part1/LauMahysDataset_KEGG_FC_EPSILON_RI_SIGNALING_PATHWAYInMahysDatasetFemale.pdf",height=2.5,width=4)
plotEnrichment(fgsea_sets[["KEGG_FC_EPSILON_RI_SIGNALING_PATHWAY"]],ranks) 
dev.off()

pdf("Graph/Part1/LauMahysDataset_ValidatePathwaysInMahysDatasetFemale.pdf",height=6,width=8)
plotGseaTable(fgsea_sets[targetPathways], ranks, fgseaRes, gseaParam=0.5)
dev.off()
#although inhibite, not significant


#-----------GSEA between LOAD (early and late pathology) and Control ---------------------
setwd("/Projects/deng/Alzheimer/syn18485175/Manuscript/Microglia")
MahysDataset=readRDS("MahysDataset.rds")
m_df<- msigdbr(species = "Homo sapiens", category = "C2",subcategory="CP:KEGG")  #CP:KEGG
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
C13=subset(MahysDataset,idents=13)
C13tmp=subset(C13,pathology.group%in%c("early-pathology","no-pathology"))
C13tmp=subset(C13tmp,Gender=="Female")
MahysDatasetdeg <- wilcoxauc(C13tmp, 'pathology.group')
clusterCell<- MahysDatasetdeg %>% dplyr::filter(group == "early-pathology") %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
FemaleInMahysDatasetCtrl_ranks<- deframe(clusterCell)
FemaleInMahysDatasetCtrl_fgseaRes<- fgseaMultilevel(fgsea_sets, stats = FemaleInMahysDatasetCtrl_ranks,eps=0)

C13tmp=subset(C13,pathology.group%in%c("early-pathology","no-pathology"))
C13tmp=subset(C13tmp,Gender=="Male")
MahysDatasetdeg <- wilcoxauc(C13tmp, 'pathology.group')
clusterCell<- MahysDatasetdeg %>% dplyr::filter(group == "early-pathology") %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
MaleInMahysDatasetCtrl_ranks<- deframe(clusterCell)
MaleInMahysDatasetCtrl_fgseaRes<- fgseaMultilevel(fgsea_sets, stats = MaleInMahysDatasetCtrl_ranks,eps=0)


pathwayCategory=read.table("/Projects/deng/Public/GSEA/KEGGPathwayCategory.txt",header=T,sep="\t")
all(MaleInMahysDatasetCtrl_fgseaRes$pathway==FemaleInMahysDatasetCtrl_fgseaRes$pathway) #TRUE
MaleInMahysDatasetCtrl_fgseaRes$Gender="Male"
FemaleInMahysDatasetCtrl_fgseaRes$Gender="Female"
dataT=rbind(MaleInMahysDatasetCtrl_fgseaRes,FemaleInMahysDatasetCtrl_fgseaRes)
sigPath=dataT[dataT$pval<0.01,]
dataT=dataT[dataT$pathway %in% unique(sigPath$pathway),]
dataT=dataT[order(dataT$pval),]
dataT$pathway=tolower(dataT$pathway)
all(dataT$pathway%in%pathwayCategory$Pathway)
dataT$Pathway=dataT$pathway
dataT=merge(dataT,pathwayCategory,by="Pathway")
dataT=dataT[order(dataT$Group,dataT$SubGroup,dataT$pval),]
dataT$PathwayName=factor(dataT$PathwayName,levels=rev(unique(dataT$PathwayName)))
t=ggplot(dataT,aes(Gender,PathwayName,size=-1*log(pval),color=NES,fill=NES))+geom_point()+
scale_color_gradient2(low="blue",mid="white",high = "red",midpoint = 0)+
scale_fill_gradient2(low="blue",mid="white",high = "red",midpoint = 0)+
theme_bw()+
theme(axis.title.x=element_blank(),axis.title.y=element_blank())
pdf("Graph/Part1/KEGGEnrichBetEarlyNoPathologySplit8Sex.pdf",width=5.5,height=5)
print(t)
dev.off()
fwrite(dataT, file="Graph/Part1/KEGGEnrichBetEarlyNoPathologySplit8Sex.txt", sep="\t", sep2=c(c("","|","")))


#-----------visualization of some importent genes ---------------------
MahysDataset=readRDS("MahysDataset.rds") 
#https://www.genome.jp/pathway/hsa04666
data=read.table("GSEA/Cluster13KEGGpathway.txt",header=T,sep="\t")
FCGamaPhagocytosis=data[data$Gender=="Female"&data$pathway=="KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS",]
FCGamaPhagocytosisGene=unlist(strsplit(FCGamaPhagocytosis$leadingEdge," ")[[1]])
t=DotPlot(MahysDataset,features=FCGamaPhagocytosisGene,group.by="OldCellType")
graphData=t$data
Ex=graphData[graphData$id=="Ex",]
Ast=graphData[graphData$id=="Ast",]
Oli=graphData[graphData$id=="Oli",]
End=graphData[graphData$id=="End",]
Ref=cbind(Ex,Ast,Oli,End)
Ref$pct=Ref[,2]+Ref[,7]+Ref[,12]+Ref[,17]
Ref=Ref[order(Ref$pct),]
g=ggplot(graphData, aes(factor(features.plot,levels=rownames(Ref)),id,size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_colour_viridis_c()+
labs(color="avg.exp,scaled",size="pct.exp",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),angle=90,hjust=1,vjust=0.5),axis.text.y = element_text(size=rel(1.0))) 
pdf("Graph/Part1/FCGamaPhagocytosisGeneCellType.pdf",width=8,height=3.5)
print(g)
dev.off()


MicSpePhago=c("PRKCA", "INPP5D", "LIMK2", "LYN", "HCK", "GAB2", "PAK1", "DNM2", "PIK3R5", "DOCK2", "PTPRC", "ASAP1", "PRKCB",  "PRKCE", "PLD1", "FCGR2A", "MAP2K1", "ARPC3", "PLCG2", "FCGR3A", "PIP5K1A", "SCIN", "PLD2", "RAF1", "RPS6KB1", "ARPC2", "MYO10", "CDC42", "FCGR1A", "WASL",  "DNM1L")
C13=subset(MahysDataset,ident=13)
C13=subset(C13,Gender=="Male")
C13$pathology.group=factor(C13$pathology.group,levels=rev(c("no-pathology","early-pathology","late-pathology")))
Idents(C13)=C13$pathology.group
#C13Expr=AverageExpression(C13)[["RNA"]]
#MicSpePhagoExpr=C13Expr[MicSpePhago,]
#pdf("Graph/Part1/FCGGamaDotBetweenPathologyPhaseHeatmap.pdf",width=3,height=4) #
#pheatmap(MicSpePhagoExpr,clustering_method="ward.D2",border=NA,scale="row",cluster_col=F)
#dev.off()
t=DotPlot(C13,features=MicSpePhago)
g=ggplot(t$data, aes(factor(features.plot,levels=rev(unique(features.plot))),id,size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_colour_viridis_c()+
labs(color="avg.exp,scaled",size="pct.exp",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),angle=90,hjust=1,vjust=0.5),axis.text.y = element_text(size=rel(1.0))) 
pdf("Graph/Part1/FCGGamaDotBetweenPathologyPhaseMale.pdf",width=8,height=3) #
print(g)
dev.off()

MapkGene=c("CD14","HSPA1A","CACNA1A","MEF2C","STK4","HSPA1B","MAP2K1","TAOK1","STK3","RPS6KA1","TAB2","FOS","TGFBR1","RPS6KA3","CACNA1D","RAP1B","MAP3K8","MAP4K4","ARRB1","MKNK1","RASA1","RPS6KA2","MAPK1","MAP3K13","GMahysDatasetD45B","JUN","MAP2K5","PAK1","MAP4K2","RPS6KA4","CDC42","HSPA8","PLA2G4A","MECOM","MAPK14","MAP3K3","DUSP1","TGFB1","ATF2","PPM1B","TNFRSF1A","PPP3CB","ARRB2","NLK","PTPN7","AKT1","GNA12","MAP2K4","PPP3CC","TAOK3","RASGRP3","TP53","HSPA6 ","MAPK11","CACNA1E","CHUK","MAP3K5","IKBKG","ECSIT","MKNK2","PAK2","DDIT3","HSPB1","MAP4K3","PRKCA","GRB2","MAPK8IP3","NF1","STMN1","IKBKB","NFATC2","AKT3","HSPA6","RPS6KA5","RASA2","MAPK10","MAP2K2","CRK","RAPGEF2","IL1B","PRKCB","DUSP16","KRAS","MAP2K6","BRAF","RAF1","MAP4K1","CACNG8","MAP3K1","PPM1A","NR4A1")
C13=subset(MahysDataset,ident=13)
C13=subset(C13,Gender=="Male")
C13$pathology.group=factor(C13$pathology.group,levels=c("no-pathology","early-pathology","late-pathology"))
Idents(C13)=C13$pathology.group
t=DotPlot(C13,features=MapkGene)
g=ggplot(t$data, aes(id,factor(features.plot,levels=rev(unique(features.plot))),size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_colour_viridis_c()+
labs(color="avg.exp,scaled",size="pct.exp",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),angle=90,hjust=1,vjust=0.5),axis.text.y = element_text(size=rel(1.0))) 
pdf("Graph/Part1/Female_MapkGeneDotBetweenPathologyPhase.pdf",width=3,height=12) #
print(g)
dev.off()
data=t$data
data$features.plot=factor(data$features.plot,levels=(unique(data$features.plot)))
data=data[order(data$features.plot),]
write.table(data,file="Graph/Part1/Male_MAPK_Expr_Pathology.txt",sep="\t",quote=F,row.names=F)


Insulin=unique(c("INPP5D","RAPGEF1","GRB2","ACACA","IKBKB","MKNK1","PDE3B","IRS2","FLOT1","PIK3R5","PRKAG2","FLOT2","PIK3CA","CBLB","MAP2K2","MAP2K1","PRKCI","CALM3","EIF4E","SHC3","RAF1","RPTOR","MAPK10","RPS6KB1","RHOQ","PPP1CB","FOXO1","BRAF","PDPK1","ACACB","INPP5K","PRKAG1","PDE3A","PRKAR1A","GSK3B","EXOC7","TSC2","HK3","RHEB","KRAS","HK1","HK2","CRK","MAPK1","PPP1R3B","PIK3CG","PTPN1","SORBS1","AKT2","RAPGEF1","PDE3B","GRB2","INPP5D","PIK3R5","ACACA","PIK3CA","IRS2","FLOT1","MKNK1","PRKAG2","IKBKB","AKT3","CBLB","PDPK1","MAPK1","RHOQ","MAPK10","EIF4E","MAP2K2","FLOT2","RPS6KB1","CRK","PHKB","PRKCI","SHC3","CALM3","FOXO1","TSC2","GSK3B","PYGB","KRAS","PPP1CB","PRKAR1A","BRAF","RAF1","PRKAG1","HK3","ACACB","INPP5K","PDE3A"))
C13=subset(MahysDataset,ident=13)
C13=subset(C13,Gender=="Male")
C13$pathology.group=factor(C13$pathology.group,levels=c("no-pathology","early-pathology","late-pathology"))
Idents(C13)=C13$pathology.group
t=DotPlot(C13,features=Insulin)
g=ggplot(t$data, aes(id,factor(features.plot,levels=rev(unique(features.plot))),size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_colour_viridis_c()+
labs(color="avg.exp,scaled",size="pct.exp",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),angle=90,hjust=1,vjust=0.5),axis.text.y = element_text(size=rel(1.0))) 
pdf("Graph/Part1/Male_InsulinDotBetweenPathologyPhase.pdf",width=3,height=8) #
print(g)
dev.off()
data=t$data
data$features.plot=factor(data$features.plot,levels=(unique(data$features.plot)))
data=data[order(data$features.plot),]
write.table(data,file="Graph/Part1/Male_Insulin_Expr_Pathology.txt",sep="\t",quote=F,row.names=F)

C13=subset(MahysDataset,idents=13)
C13tmp=subset(C13,pathology.group%in%c("early-pathology","no-pathology"))
C13tmp=subset(C13tmp,Gender=="Female")
MahysDatasetdeg <- wilcoxauc(C13tmp, 'pathology.group')
InsulinDEG_early=MahysDatasetdeg[MahysDatasetdeg$feature%in%Insulin&MahysDatasetdeg$group == "early-pathology",]
InsulinDEG_early$feature=factor(InsulinDEG_early$feature,levels=Insulin)
InsulinDEG_early=InsulinDEG_early[order(InsulinDEG_early$feature),]

C13tmp=subset(C13,pathology.group%in%c("late-pathology","no-pathology"))
C13tmp=subset(C13tmp,Gender=="Female")
MahysDatasetdeg <- wilcoxauc(C13tmp, 'pathology.group')
InsulinDEG_late=MahysDatasetdeg[MahysDatasetdeg$feature%in%Insulin&MahysDatasetdeg$group == "late-pathology",]
InsulinDEG_late$feature=factor(InsulinDEG_late$feature,levels=Insulin)
InsulinDEG_late=InsulinDEG_late[order(InsulinDEG_late$feature),]

C13tmp=subset(C13,pathology.group%in%c("late-pathology","early-pathology"))
C13tmp=subset(C13tmp,Gender=="Female")
MahysDatasetdeg <- wilcoxauc(C13tmp, 'pathology.group')
InsulinDEG_late_early=MahysDatasetdeg[MahysDatasetdeg$feature%in%Insulin&MahysDatasetdeg$group == "late-pathology",]
InsulinDEG_late_early$feature=factor(InsulinDEG_late_early$feature,levels=Insulin)
InsulinDEG_late_early=InsulinDEG_late_early[order(InsulinDEG_late_early$feature),]

all(InsulinDEG_early$feature==InsulinDEG_late$feature)
InsulinDEG=data.frame("Symbol"=InsulinDEG_early$feature,"early vs no"=InsulinDEG_early$pval,"late vs no"=InsulinDEG_late$pval,"late vs early"=InsulinDEG_late_early$pval)
write.table(InsulinDEG,file="Graph/Part1/Female_Insulin_Pathology.txt",sep="\t",quote=F)
