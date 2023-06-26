library(Seurat)

setwd("/Projects/deng/Alzheimer/syn18485175")
data <- Read10X(data.dir = "data/")
#data was downloaded from syn18485175 using filtered count matrix
MathysDataset <- CreateSeuratObject(counts = data, project = "MathysDataset")
Info=read.table("/Projects/deng/Alzheimer/syn18485175/Phenotype4Cell.txt",header=TRUE,row.names=1,sep="\t")
MathysDataset$Gender=Info$Gender
MathysDataset$Sex=Info$Sex
MathysDataset$OldCellType=Info$broad.cell.type
MathysDataset$OldCluster=Info$pre.cluster
MathysDataset$Statues=Info$Statues
MathysDataset$OriginalSubCluster=Info$Subcluster
MathysDataset$ProjectID=Info$ProjectID
MathysDataset$BraakStage=paste0("Stage_",Info$BraakStage,sep="")
MathysDataset$Gpath=Info$gpath
MathysDataset$Amyloid=Info$amyloid
MathysDataset$Plaq=Info$plaq_n
MathysDataset$Cogdx=Info$cogdx
MathysDataset$Ceradsc=Info$ceradsc
MathysDataset$NFT=Info$nft
MathysDataset$cogn_global_lv=Info$cogn_global_lv
MathysDataset$gpath_3neocort=Info$gpath_3neocort
MathysDataset$amyloid.group=Info$amyloid.group
MathysDataset$caa_4gp=Info$caa_4gp
MathysDataset$APOEScore=Info$apoeScore
MathysDataset$Age=Info$age_death
MathysDataset$APOEGenotype=Info$APOE
MathysDataset$pathology.group=Info$pathology.group

MathysDataset <- NormalizeData(MathysDataset, normalization.method = "LogNormalize", scale.factor = 10000)
MathysDataset <- FindVariableFeatures(MathysDataset, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MathysDataset)
MathysDataset <- ScaleData(MathysDataset, features = all.genes)
MathysDataset <- RunPCA(MathysDataset, features = VariableFeatures(object = MathysDataset))
DimHeatmap(MathysDataset , dims = 20:45, cells = 500, balanced = TRUE)
ElbowPlot(MathysDataset,ndims = 50)
MathysDataset<- RunUMAP(MathysDataset, reduction = "pca", dims = 1:30)
MathysDataset<- FindNeighbors(MathysDataset, reduction = "pca", dims = 1:30)
MathysDataset<- FindClusters(MathysDataset, resolution = 0.5)
# Visualization
DimPlot(MathysDataset, reduction = "umap", label=TRUE)
DimPlot(MathysDataset, reduction = "umap", group.by="OldCellType")
DimPlot(MathysDataset, reduction = "umap", label = TRUE,group.by="OriginalSubCluster")
DimPlot(MathysDataset, reduction = "umap", group.by="Statues",cols = c("firebrick3","green"))
DimPlot(MathysDataset, reduction = "umap", group.by="Statues",split.by="Gender",cols = c("firebrick3","green"))
DimPlot(MathysDataset, reduction = "umap", group.by="Gender",split.by="Statues")

FeaturePlot(MathysDataset,features=c("nFeature_RNA","BraakStageID","Gpath","Amyloid","Plaq","Cogdx","Ceradsc","NFT","cogn_global_lv","APOEScore","AgeDeath"))

saveRDS(MathysDataset,"MathysDataset.rds")
