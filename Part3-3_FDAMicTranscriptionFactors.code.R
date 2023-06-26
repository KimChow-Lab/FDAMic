setwd("D:/Alzheimer/syn18485175/cluster13/monocle/SCENIC")
library(SCENIC)
library(pheatmap)

org="hgnc"
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir="databases",  dbs=dbs, datasetTitle="AD", nCores=1) #check the certutil -hashfile databases/hg19-tss-centered-10kb-7species.mc9nr.feather SHA256
ADMic=readRDS("D:/Alzheimer/syn18485175/cluster13/monocle/ADMicMono.rds")
ADMic=subset(ADMic,ident=c(1,2,3,4,5,6,7,8,9,10))
exprMat  <-  as.matrix(ADMic@assays$RNA@counts)
dim(exprMat)
exprMat[1:4,1:4] 
cellInfo <-  ADMic@meta.data
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,minCountsPerGene=3*.01*ncol(exprMat),minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered, scenicOptions)

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 
library(doParalle)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat)
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
scenicOptions=readRDS("int/scenicOptions.Rds")

# Cell-type specific regulators (RSS): 
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "MicCluster"], )
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)

cellInfo <- data.frame(seuratCluster=ADMic$MicCluster)
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster), function(cells) rowMeans(getAUC(regulonAUC)[,cells]))

write.table(regulonActivity_byCellType,file="regulonActivity_byCellType.txt",sep="\t",quote=F)

data1=regulonActivity_byCellType[matrixStats::rowMaxs(as.matrix(regulonActivity_byCellType[,2:dim(regulonActivity_byCellType)[2]]))>0.05,]
gene=intersect(rownames(data1),rownames(meta))
tmp=setdiff(rownames(data1),rownames(meta))
tmp=cbind(tmp,"ESR1"="NonInteract","ESR2"="NonInteract")
rownames(tmp)=tmp[,1]
metaAll=data.frame(rbind(meta[gene,],tmp[,-1]))
data1=data1[,-1]
ann_colors = list(ESR1=c(Interact="CadetBlue1",NonInteract="WhiteSmoke"), ESR2=c(Interact ="CadetBlue1",NonInteract = "WhiteSmoke"))
pheatmap(data1,annotation_row=metaAll,annotation_colors = ann_colors,clustering_method="ward.D2",scale="row",colorRampPalette(c("white","white","white","Wheat1","Salmon1","Firebrick3"))(50),cluster_col=FALSE,show_rownames=TRUE)


