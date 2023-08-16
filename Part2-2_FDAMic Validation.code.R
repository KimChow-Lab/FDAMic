###########validate the FDAMic using independent datasets (Morabito et al)###########################
library(Seurat)
library(harmony) #0.1.1
library(plyr)
library(RColorBrewer)
library(ggplot2)
library(scCustomize)
library(dplyr)
library(viridis)
library(ggrepel)
library(monocle3)

setwd("/data2/deng/Aging/Glia/scRNAseq/Morabito_GSE174367")
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174367
Morabito.data=Read10X_h5("rawData/GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5")
cellMetaInfo=read.csv("rawData/GSE174367_snRNA-seq_cell_meta.csv",header=T,row.names=1)
InvidualInfo=read.table("sampleInfo.txt",header=T)
InvidualInfo=InvidualInfo[InvidualInfo$SampleID%in%unique(cellMetaInfo$SampleID)&InvidualInfo$seqType%in%c("snRNA-seq"),]
cellMetaInfo$cellId=rownames(cellMetaInfo)
cellMetaInfo=merge(cellMetaInfo,InvidualInfo[,c("SampleID","Title","seqType")],by="SampleID")
rownames(cellMetaInfo)=cellMetaInfo$cellId

cells=intersect(rownames(cellMetaInfo),colnames(Morabito.data)) #61472
MorabitoDataSet <- CreateSeuratObject(counts = Morabito.data[,cells], project = "MorabitoDataSet") #58721 features across 61472 samples within 1 assay
MorabitoDataSet <- PercentageFeatureSet(MorabitoDataSet, "^MT-", col.name = "percent_mito")
selected_c <- WhichCells(MorabitoDataSet, expression = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent_mito < 20)
selected_f <- rownames(MorabitoDataSet)[Matrix::rowSums(MorabitoDataSet) > 50]
MorabitoDataSet <- subset(MorabitoDataSet, features = selected_f, cells = selected_c)
dim(MorabitoDataSet)
#26903 60589

UniqueGene=read.table("/Projects/deng/Public/Ensemble/gencode.v38.annotation.removeDupliSymbol.geneType.txt",header=T)
Protein_coding_Gene=UniqueGene[UniqueGene$Type%in%"protein_coding",] 
dim(Protein_coding_Gene)#19920     7
GeneList=intersect(Protein_coding_Gene$Symbol,rownames(MorabitoDataSet)) #18446
length(GeneList) #15745
#remove the mitochodrial, ribosome and noncoding gene
MorabitoDataSet <- MorabitoDataSet[!grepl("^MT-", rownames(MorabitoDataSet)), ] #remove mitochdrial genes
MorabitoDataSet <- MorabitoDataSet[!grepl("^HB[^(P)]", rownames(MorabitoDataSet)), ] #remove homoglobin genes
MorabitoDataSet <- MorabitoDataSet[GeneList, ] # Only counts associated with protein-coding genes were considered; 
dim(MorabitoDataSet) #  15725 60589

cells=intersect(colnames(MorabitoDataSet),rownames(cellMetaInfo)) #60589
cellMetaInfo=cellMetaInfo[cells,]
all(rownames(cellMetaInfo)==colnames(MorabitoDataSet)) #TRUE

MorabitoDataSet$Sample=cellMetaInfo$Title
MorabitoDataSet$Age=cellMetaInfo$Age
MorabitoDataSet$Gender=ifelse(cellMetaInfo$Sex%in%c("F"),"Female","Male")
MorabitoDataSet$Region="PFC"
MorabitoDataSet$Statues=cellMetaInfo$Diagnosis
MorabitoDataSet$cellType=cellMetaInfo$Cell.Type
MorabitoDataSet$TangStage=cellMetaInfo$Tangle.Stage
MorabitoDataSet$PlaqueStage=cellMetaInfo$Plaque.Stage
MorabitoDataSet$Batch=cellMetaInfo$Batch
MorabitoDataSet$cluster=cellMetaInfo$cluster
MorabitoDataSet$PMI=cellMetaInfo$PMI
MorabitoDataSet$RIN=cellMetaInfo$RIN

as.matrix(table(MorabitoDataSet$Sample,MorabitoDataSet$Statues))
MorabitoDataSet <- NormalizeData(MorabitoDataSet, normalization.method = "LogNormalize", scale.factor = 10000)
MorabitoDataSet <- FindVariableFeatures(MorabitoDataSet, selection.method = "vst", nfeatures = 2000)
MorabitoDataSet <- ScaleData(MorabitoDataSet)
MorabitoDataSet <- RunPCA(MorabitoDataSet, features = VariableFeatures(object = MorabitoDataSet))
MorabitoDataSet <- RunHarmony(MorabitoDataSet, group.by.vars=c("Sample"))
MorabitoDataSet<- RunUMAP(MorabitoDataSet, reduction = "harmony", dims = 1:50)
MorabitoDataSet<- FindNeighbors(MorabitoDataSet, reduction = "harmony", dims = 1:50)
MorabitoDataSet<- FindClusters(MorabitoDataSet, resolution = 0.5)

setwd("/data2/deng/Alzheimer/Microglia/Morabito_GSE174367")
pdf("Morabito_snRNA/Morabito_snRNA_Label.pdf",width=9)
DimPlot(MorabitoDataSet,reduction ="umap",label=TRUE)&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
tiff("Morabito_snRNA/Morabito_snRNA_CellTypeMarker.tiff",width=1000,height=1000)
FeaturePlot(MorabitoDataSet,reduction ="umap",features=c("CAMK2A","NRGN","RBFOX3","GAD1","GAD2","AQP4","GFAP","MBP","MOBP","PLP1","CSF1R","CD74","VCAN","FLT1","CLDN5","TAGLN"))&NoLegend()&theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
#endothelial (CLDN5), 
#mural (TAGLN), 
#fibroblGlia (CEMIP), 
#ependymal (TTR), 
#Gliaroglia (GFAP), 
#neurons (CAMK2A), 
#oligodendrocytes (PLP1), 
#microglial (CSF1R)

pal <- viridis(n = 10, option = "D")
pdf("Morabito_snRNA/Morabito_snRNA_Microglia_Marker.pdf",height=10,width=10)
Plot_Density_Custom(seurat_object = MorabitoDataSet, features = c("CSF1R","P2RY12","CD74","MRC1"), custom_palette = pal)&NoLegend()&
theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())&
theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

pdf("Morabito_snRNA/Morabito_snRNA_Sample.pdf",width=10)
DimPlot(MorabitoDataSet, reduction ="umap",group.by="Sample")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("Morabito_snRNA/Morabito_snRNA_Sample.pdf",width=10)
DimPlot(MorabitoDataSet, reduction ="umap",group.by="Sample")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("Morabito_snRNA/Morabito_snRNA_CellType.pdf",width=10)
DimPlot(MorabitoDataSet, reduction ="umap",group.by="cellType")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("Morabito_snRNA/Morabito_snRNA_cluster.pdf",width=10)
DimPlot(MorabitoDataSet, reduction ="umap",group.by="cluster")&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
saveRDS(MorabitoDataSet,file="Morabito_snRNA/Morabito_snRNA.rds")

MorabitoMic=subset(MorabitoDataSet,idents=c(3,20)) #15725 features across 4196 samples within 1 assay
MorabitoMic<- RunUMAP(MorabitoMic, reduction = "pca", dims = 1:20)
MorabitoMic<- FindNeighbors(MorabitoMic, reduction = "pca", dims = 1:20)
MorabitoMic<- FindClusters(MorabitoMic, resolution = 0.5)
pdf("Morabito_snRNA/Morabito_snRNA_Mic_Cluster.pdf")
DimPlot(MorabitoMic,label=TRUE)+NoLegend()
dev.off()
pdf("Morabito_snRNA/Morabito_snRNA_Mic_Sample.pdf")
DimPlot(MorabitoMic,group.by="Sample")+NoLegend()
dev.off()
tiff("Morabito_snRNA/Morabito_snRNA_Mic_BrainCellTypeMarker.tiff",width=1600,height=1600)
FeaturePlot(MorabitoMic,c("NRGN","CAMK2A","NRXN1","GAD1","CTNNA3","DIAPH3","MBP","GFAP","VCAN","CD74","FLT1","PDGFRA","CLDN5","FLT1","ABCB1","EBF1"))
dev.off()
tiff("Morabito_snRNA/Morabito_snRNA_Mic_MicSubTypeMarker.tiff",width=1000,height=1000)
FeaturePlot(MorabitoMic,c("MRC1","F13A1","P2RY12","VSIG4","FCGR3A","RPS19","SPP1","PLEKHA7"))
dev.off()

saveRDS(MorabitoMic,file="Morabito_snRNA/Morabito_snRNA_Mic.rds")


setwd("/Projects/deng/Alzheimer/syn18485175")
MorabitoDataSet=readRDS("/data2/deng/Alzheimer/Microglia/Morabito_GSE174367/Morabito_snRNA/Morabito_snRNA.rds")
tiff("Manuscript/Microglia/Graph/Part2/Morabito_cellType.tiff",width=450,height=400)
DimPlot(MorabitoDataSet,group.by="cellType")&NoLegend()&NoAxes()&theme(plot.title=element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
tiff("Manuscript/Microglia/Graph/Part2/Morabito_snRNA_Mic_MicSubTypeMarker.tiff",width=600,height=600)
FeaturePlot(MorabitoDataSet,c("CX3CR1","CD74","P2RY12","CSF1R"))&NoLegend()&NoAxes()&theme(plot.title=element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

cellTypeList=c("BAMic","HomMic","NorARMic","ARMic","FDAMic","DysMic")
Morabito_snRNA_Mic_Seurat=readRDS("/data2/deng/Alzheimer/Microglia/Morabito_GSE174367/Morabito_snRNA/Monocle/Morabito_snRNA_Mic.rds")
Morabito_snRNA_Mic_Seurat=subset(Morabito_snRNA_Mic_Seurat,SubType%in%cellTypeList)
pal <- viridis(n = 10, option = "D")
tiff("Manuscript/Microglia/Graph/Part2/Morabito_FDAMicMarker_Density_legend.tiff",height=400,width=600)
Plot_Density_Custom(seurat_object = Morabito_snRNA_Mic_Seurat, features = c("VSIG4","SPP1","RPS19","C1QB"), custom_palette = pal)&NoLegend()&
theme(plot.title=element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())&
theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

tiff("Manuscript/Microglia/Graph/Part2/Morabito_FDAMicMarker_Density_legend.tiff",height=400,width=500)
Plot_Density_Custom(seurat_object = Morabito_snRNA_Mic_Seurat, features = c("VSIG4","SPP1","RPS19","C1QB"), custom_palette = pal)&
#NoLegend()&
theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())&
theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

colors_list <- c("forestgreen","NavajoWhite","LightBlue", "SlateBlue", "Violet", "gold")
markerGene=c("F13A1","MRC1","SPP1","P2RY12","CX3CR1","FCGR3A","VSIG4","CD14","RPS19","SPP1","C1QB","ARMC9","PLEKHA7","RPS24","FTH1","FTL")
pdf("Manuscript/Microglia/Graph/Part2/Morabito_FDAMicMarker_VlnPlot.pdf",height=9,width=5)
scCustomize::Stacked_VlnPlot(seurat_object = Morabito_snRNA_Mic_Seurat, features = markerGene, x_lab_rotate = TRUE, plot_spacing = 0.05,colors_use = colors_list)&
theme(axis.title.x=element_blank())
dev.off()


MorabitoMic=Morabito_snRNA_Mic_Seurat
MorabitoMic$SubType="Unsign"
MorabitoMic$BraakStage=MorabitoMic$TangStage
MorabitoMic$MicCluster="Unsign"
MorabitoMic@meta.data=MorabitoMic@meta.data[,c("nCount_RNA","nFeature_RNA","percent_mito","MicCluster","Sample","Statues","SubType","Age","Gender","BraakStage")]
data <- GetAssayData(object = MorabitoMic, slot = "counts")
celldata <- as.data.frame(MorabitoMic@meta.data)
genedata <- as.data.frame(x = row.names(MorabitoMic), row.names = row.names(MorabitoMic))
colnames(genedata) <- "gene_short_name"
MorabitoMic.cds <- new_cell_data_set(data, cell_metadata = celldata, gene_metadata = genedata)
colData(MorabitoMic.cds)
MorabitoMic.cds <- preprocess_cds(MorabitoMic.cds, num_dim = 20)
pdf("Morabito_snRNA/Monocle/variance.pdf")
plot_pc_variance_explained(MorabitoMic.cds)
dev.off()
MorabitoMic.cds <- align_cds(MorabitoMic.cds, alignment_group = "Sample") #significantly batch effect existed in
#remember to cite Haghverdi L, Lun ATL, Morgan MD, Marioni JC (2018). 'Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors.' Nat. Biotechnol., 36(5), 421-427. doi: 10.1038/nbt.4091
MorabitoMic.cds <- reduce_dimension(MorabitoMic.cds,preprocess_method = 'Aligned')
#MorabitoMic.cds <- reduce_dimension(MorabitoMic.cds,preprocess_method = 'PCA')
MorabitoMic.cds = cluster_cells(cds = MorabitoMic.cds, resolution=0.005,reduction_method ='UMAP')
MorabitoMic.cds$cluster=clusters(MorabitoMic.cds)

#mapping UMAP embedding to seurat object
all(colnames(MorabitoMic)==colnames(MorabitoMic.cds)) #TRUE
MorabitoMic$MicCluster=clusters(MorabitoMic.cds)
monocle_umap <- MorabitoMic.cds@int_colData@listData$reducedDims$UMAP
colnames(monocle_umap) <- c('UMAP_1', 'UMAP_2')
MorabitoMic@reductions$umap@cell.embeddings=monocle_umap
pal <- viridis(n = 10, option = "D")
geneList=c("F13A1","P2RY12","CX3CR1","VSIG4","CD14","FCGR3A","RPS19","SPP1","MOG","PLEKHA7")
tiff("Morabito_snRNA/Monocle/Mic_Seurat__CellTypeMarkerDensity.tiff",height=1000,width=1000)
Plot_Density_Custom(seurat_object = MorabitoMic, features = geneList, custom_palette = pal)&NoLegend()&
theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())&
theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
pdf("Morabito_snRNA/Monocle/Mic_Seurat_Order_subtypeMarker.pdf",width=5,height=6)
Stacked_VlnPlot(seurat_object =MorabitoMic,group.by="MicCluster",plot_spacing = 0,features = geneList,x_lab_rotate = TRUE)
dev.off()


### defination the subytpe name for each subcluster ##############
MorabitoMic$MicCluster=factor(MorabitoMic$MicCluster,levels=c(10,2,12,7,8,11,9,4,1,5,3,6,15,14,13))
Idents(MorabitoMic)=MorabitoMic$MicCluster
new.cluster.ids <- c("BAMic","HomMic","HomMic","HomMic","HomMic","HomMic","HomMic","HomMic","NorARMic","NorARMic","ARMic","FDAMic","FDAMic","DysMic","OliMic")
names(new.cluster.ids) <- levels(MorabitoMic)
MorabitoMic <- RenameIdents(MorabitoMic, new.cluster.ids)
MorabitoMic$SubType=Idents(MorabitoMic)
all(colnames(MorabitoMic)==colnames(MorabitoMic.cds)) 
MorabitoMic.cds$subType=MorabitoMic$SubType
saveRDS(MorabitoMic.cds,file="Morabito_snRNA/Monocle/Morabito_snRNA_Mic.monocle3.cds.rds")



### visualization of the subtype of microgia from Morabito ##############
Morabito_snRNA_Mic=readRDS("/data2/deng/Alzheimer/Microglia/Morabito_GSE174367/Morabito_snRNA/Monocle/Morabito_snRNA_Mic.monocle3.cds.rds")
table(Morabito_snRNA_Mic$subType)
#BAMic   HomMic NorARMic    ARMic   FDAMic   DysMic   OliMic
#239     2002      762      402      423      166      202
metaInfo=colData(Morabito_snRNA_Mic)
valid_cells=rownames(metaInfo[metaInfo$MicCluster%in%c(1:12,14,15),]) #remove OliMic
MorabitoMic.sub.cds <- Morabito_snRNA_Mic[,valid_cells]
table(colData(MorabitoMic.sub.cds)$subType)
#BAMic   HomMic NorARMic    ARMic   FDAMic   DysMic   OliMic
#239     2002      762      402      423      166        0
colData(MorabitoMic.sub.cds)$subType <- factor(colData(MorabitoMic.sub.cds)$subType,levels=c("BAMic","HomMic","NorARMic","ARMic","FDAMic","DysMic"))
colors_list <- c("forestgreen","NavajoWhite","LightBlue", "SlateBlue", "Violet", "gold")
tiff("Manuscript/Microglia/Graph/Part2/Morabito_MonoMicSubType.tiff",width=500,height=400)
plot_cells(MorabitoMic.sub.cds, color_cells_by="subType",label_cell_groups = FALSE,group_label_size = 5,cell_size =1)+
scale_color_manual(values = colors_list)+
theme(legend.position="none",axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())&
theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

count=data.frame(table(MorabitoMic.sub.cds$subType,MorabitoMic.sub.cds$Gender,MorabitoMic.sub.cds$Statues))
colnames(count)=c("subType","Gender","Group","Number")
count$Category=paste0(count$Gender,count$Group,sep="")
total=data.frame(table(paste0(MorabitoMic.sub.cds$Gender,MorabitoMic.sub.cds$Statues,sep="")))
colnames(total)=c("Category","Total")
graphData=merge(count,total,by="Category")
graphData$Ratio=graphData$Number/graphData$Total
graphData$Gender=factor(graphData$Gender,levels=c("Male","Female"))
graphData$Group=factor(graphData$Group,levels=c("Control","AD"))
t=ggplot(graphData, aes(subType, Ratio, fill=Group)) +
  geom_bar(stat="identity",position="fill") +
  scale_fill_manual(values = c("SlateBlue", "Coral"))+
  scale_y_continuous(expand=c(0,0))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0))+
  labs(size="",x="",y="Cell ratio",title="")+
  theme(axis.text.x = element_text(size=rel(1.0),angle=90))+
  facet_grid(Gender~.)
pdf("Manuscript/Microglia/Graph/Part2/Morabito_Mic_Monocle_cellType_GroupDis.pdf",width=4,height=5)
print(t)
dev.off()

pg <- ggplot_build(t)
data=graphData
data$ScaledValue=pg$data[[1]]$y
write.table(graphData,file="Manuscript/Microglia/Graph/Part2/Morabito_Mic_Monocle_cellType_GroupDis.txt",sep="\t",row.names=F,quote=F)


tiff("Manuscript/Microglia/Graph/Part2/Morabito_MonoMicBraakStage.tiff",width=400,height=400)
plot_cells(MorabitoMic.sub.cds, color_cells_by="BraakStage",label_cell_groups = FALSE,group_label_size = 5,cell_size =0.8)+
#scale_color_gradient(low="gray90",high = "firebrick3")+
scale_color_brewer(palette="Spectral",direction = -1)+
theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())&
theme(legend.position="none",panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

order=c("BAMic","HomMic","NorARMic","ARMic","FDAMic","DysMic")
Braak=data.frame(table(MorabitoMic.sub.cds$subType,MorabitoMic.sub.cds$BraakStage))
colnames(Braak)=c("subType","Stage","Number")
total=data.frame(table(MorabitoMic.sub.cds$BraakStage))
colnames(total)=c("Stage","Total")
graphData=merge(Braak,total,by="Stage")
graphData$Ratio=graphData$Number/graphData$Total
t=ggplot(graphData, aes(factor(subType,levels=rev(order)),Ratio, fill=Stage)) +
  geom_bar(stat="identity",position="fill") +
  scale_fill_brewer(palette="Spectral",direction = -1)+
  scale_y_continuous(expand=c(0,0))+
  theme_bw()+coord_flip()+
  theme(axis.text.x = element_text(angle = 0),axis.title.x = element_blank(),axis.title.y = element_blank())
pdf("Manuscript/Microglia/Graph/Part2/Morabito_Mic_Monocle_cellType_BraakDis.pdf",width=4,height=4)
print(t)
dev.off()


