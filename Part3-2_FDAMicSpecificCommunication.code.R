library(Seurat)
library(ggplot2)
library(scCustomize)

#################Combine Yang et al and Mathys et al #########################
setwd("/data2/deng/Alzheimer/Microglia/Yang_Vascular")
#https://cells-test.gi.ucsc.edu/?ds=brain-vasc-atlas
data <- Read10X(data.dir="rawData/")
Yang <- CreateSeuratObject(counts = data, project = data)

CellMetaInfo=read.table("rawData/meta.tsv",header=T,row.names=1)
all(rownames(CellMetaInfo)==colnames(Yang)) #TRUE
Yang$Statues=CellMetaInfo$Treat
Yang$Gender=CellMetaInfo$Gender
Yang$Age=CellMetaInfo$Age
Yang$Sample=CellMetaInfo$Sample
Yang$Region=CellMetaInfo$Region
Yang$CellType=CellMetaInfo$Cell_Type
Yang$APOE4=CellMetaInfo$APOE4
Yang$APOE34=CellMetaInfo$APOE34
Yang$Batch=CellMetaInfo$Batch

Yang <- NormalizeData(Yang, normalization.method = "LogNormalize", scale.factor = 10000)
Yang <- FindVariableFeatures(Yang, selection.method = "vst", nfeatures = 2000)
Yang <- ScaleData(Yang, verbose = FALSE)
Yang <- RunPCA(Yang, features = VariableFeatures(object = Yang))
Yang <- RunHarmony(Yang, group.by.vars=c("Sample"))
Yang<- RunUMAP(Yang, reduction = "harmony", dims = 1:50)
Yang<- FindNeighbors(Yang, reduction = "pca", dims = 1:50)
Yang<- FindClusters(Yang, resolution = 0.5)

tiff("cellType.tiff",width=450,height=400)
DimPlot(Yang,group.by="CellType",label=T,label.size = 6)&NoLegend()&NoAxes()&theme(plot.title=element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
tiff("cellTypeMarker.tiff",width=400,height=800)
FeaturePlot(Yang,c("NRGN","GAD1","MBP","GFAP","VCAN","CSF1R","FLT1"),cols =c("lightgrey","red"),ncol=2)&NoLegend()&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

Yang$Group=paste0(Yang$CellType,"_",Yang$Gender,"_",Yang$Statues,sep="")
phagocytosis=c("CD14","FCGR3A","FCGR1A","SLC11A1","HCK","CLCN3","RAB39A","ITGB1","PRKCD","PIK3CA","PAK1","HMGB1","UNC13D","ARHGAP12","PIK3R1","P2RY6","PIP5K1A")
ImmuneInYang=subset(Yang,CellType%in%c("Microglia/Mφ"))
t=DotPlot(ImmuneInYang,features=phagocytosis,group.by="Group")
g=ggplot(t$data, aes(id, factor(features.plot,levels=rev(unique(features.plot))), size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_colour_viridis_c()+
labs(color="avg.exp.scaled",size="pct.exp",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),angle=90,vjust=0.5,hjust=1),axis.text.y = element_text(size=rel(1.0))) 
pdf("PhagocytosisGene.pdf",width=4,height=6)
print(g)
dev.off()
saveRDS(Yang,file="YangData.rds")

gene=c("ADAM28","ALOX5AP","APOC1","APOE","B2M","BST2","C1QA","C1QB","CALR","CD14","GRN","HLA-A","HLA-B","HLA-C","HLA-E","PSAP","PTMA","SELPLG","SPP1","ST6GAL1","SYTL3","WNT2B","TREM2","ITGA9","ALOX5","PLD1","RPSA","CD74","BMPR1A","CD81","KCNQ3","CCR1","CSF1R","CSF2RA","CSF3R","SLC1A5","HAVCR2","LINGO1","NPTN","CD68","TYROBP","PLXNC1","ADGRG1","CD63","LYVE1")
pdf("Ligand.pdf",width=15,height=5)
DotPlot(Yang,features=gene,group.by="CellType")+RotatedAxis()
dev.off()

Yang$Source="Yang"
Mathys=readRDS("/Projects/deng/Alzheimer/syn18485175/ADPhenotype.rds") 
Mathys$Source="Mathys"
Info=read.table("/Projects/deng/Alzheimer/syn18485175/Phenotype4Cell.txt",header=TRUE,row.names=1,sep="\t")
Mathys$Gender=Info$Gender
Mathys$CellType=Info$broad.cell.type
Mathys$APOE=Info$APOE
Mathys$Age=Info$age_death
Mathys$Sample=Mathys$projectID

CombineSCRNA=merge(Yang,Mathys)
CombineSCRNA@meta.data=CombineSCRNA@meta.data[,c("nCount_RNA","nFeature_RNA","Source","Sample","Statues","Gender","Age","CellType")]
CombineSCRNA <- NormalizeData(CombineSCRNA, normalization.method = "LogNormalize", scale.factor = 10000)
CombineSCRNA <- FindVariableFeatures(CombineSCRNA, selection.method = "vst", nfeatures = 2000)
CombineSCRNA <- ScaleData(CombineSCRNA, verbose = FALSE)
CombineSCRNA <- RunPCA(CombineSCRNA, features = VariableFeatures(object = CombineSCRNA))
CombineSCRNA <- RunHarmony(CombineSCRNA, group.by.vars=c("Source"))
CombineSCRNA<- RunUMAP(CombineSCRNA, reduction = "harmony", dims = 1:50)
CombineSCRNA<- FindNeighbors(CombineSCRNA, reduction = "pca", dims = 1:50)
CombineSCRNA<- FindClusters(CombineSCRNA, resolution = 0.5)

tiff("CombinedCellType.tiff",width=450,height=400)
DimPlot(CombineSCRNA,group.by="CellType",label=T,label.size = 6)&NoLegend()&NoAxes()&theme(plot.title=element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
tiff("CombinedCellTypeMarker.tiff",width=400,height=800)
FeaturePlot(CombineSCRNA,c("NRGN","GAD1","MBP","GFAP","VCAN","CSF1R","FLT1"),cols =c("lightgrey","red"),ncol=2)&NoLegend()&NoAxes()&theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

gene=c("ADAM28","ALOX5AP","APOC1","APOE","B2M","BST2","C1QA","C1QB","CALR","CD14","GRN","HLA-A","HLA-B","HLA-C","HLA-E","PSAP","PTMA","SELPLG","SPP1","ST6GAL1","SYTL3","WNT2B","TREM2","ITGA9","ALOX5","PLD1","RPSA","CD74","BMPR1A","CD81","KCNQ3","CCR1","CSF1R","CSF2RA","CSF3R","SLC1A5","HAVCR2","LINGO1","NPTN","CD68","TYROBP","PLXNC1","ADGRG1","CD63","LYVE1")
pdf("Ligand8Combine.pdf",width=15,height=5)
DotPlot(CombineSCRNA,features=gene,group.by="CellType")+RotatedAxis()
dev.off()

gene=c("ITGA4","ALOX5","VLDLR","LRP6","LRP2","LDLR","LRP5","CHRNA4","TREM2","LSR","SORL1","ABCA1","SCARB1","SDC2","LRP8","LRP4","LRP1","HLA-F","CD1A","CD3D","LILRB1","KIR2DL1","KIR2DL3","TFRC","CD247","CD1B","KLRC2","KLRD1","KIR3DL1","LILRB2","HFE","CD3G","KLRC1","LILRA4","CD93","CSPG4","CR1","C1QBP","MPL","TSHR","ITGA3","ITGAV","ITGA2B","SCARF1","MTNR1A","ITGB1","TLR6","TLR9","ITGB2","TLR4","TLR1","RIPK1","TNFRSF1A","EGFR","SORT1","TNFRSF1B","CLEC4M","LILRA1","APLP2","ERBB2","KIR3DL2","CANX","NOTCH4","KLRK1","GPR37","GPR37L1","AR","VIPR1","ITGAM","ESAM","SELL","SELP","SELE","CD44","ITGB3","ITGA9","ITGA5","ITGB5","PTGER4","CD22","NRXN1","FZD4")
pdf("Ligand8CombineReceptor.pdf",width=20,height=5)
DotPlot(CombineSCRNA,features=gene,group.by="CellType")+RotatedAxis()
dev.off()

CombineSCRNA$CombinecellType=ifelse(CombineSCRNA$CellType%in%c("Astrocyte","Ast"),"Ast",CombineSCRNA$CellType)
CombineSCRNA$CombinecellType=ifelse(CombineSCRNA$CellType%in%c("Arterial","Capillary","Veinous","End"),"End",CombineSCRNA$CombinecellType)
CombineSCRNA$CombinecellType=ifelse(CombineSCRNA$CellType%in%c("Microglia/Mφ","Mic"),"Mic",CombineSCRNA$CombinecellType)
CombineSCRNA$CombinecellType=ifelse(CombineSCRNA$CellType%in%c("M. Fibro","P. Fibro"),"Fibro",CombineSCRNA$CombinecellType)
CombineSCRNA$CombinecellType=ifelse(CombineSCRNA$CellType%in%c("OPC","Opc"),"Opc",CombineSCRNA$CombinecellType)
CombineSCRNA$CombinecellType=ifelse(CombineSCRNA$CellType%in%c("Oli","Oligo"),"Oli",CombineSCRNA$CombinecellType)
CombineSCRNA$CombinecellType=ifelse(CombineSCRNA$CellType%in%c("Pericyte","Per"),"Per",CombineSCRNA$CombinecellType)
tiff("CombinedCellType_allType.tiff",width=450,height=400)
DimPlot(CombineSCRNA,group.by="CombinecellType",label=T,label.size = 6)&NoLegend()&NoAxes()&theme(plot.title=element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
tiff("CombinedCellType_Cluster.tiff",width=450,height=400)
DimPlot(CombineSCRNA,label=T,label.size = 6)&NoLegend()&NoAxes()&theme(plot.title=element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
CombineSCRNA$Statues=ifelse(CombineSCRNA$Statues=="AD","Alzheimer",CombineSCRNA$Statues)

new.cluster.ids <- c("End","Per","Ast","Ex","Oli","Oli","Oli","End","SMC","In","In","Ex","Fibro","Opc","Mic","Ast","Ex","Ex","Opc","Ex","Ex","Mic","Ependymal","Ex","In","T","End")
names(new.cluster.ids) <- levels(CombineSCRNA)
CombineSCRNA <- RenameIdents(CombineSCRNA, new.cluster.ids)
tiff("CombinedCellType_myType.tiff",width=450,height=400)
DimPlot(CombineSCRNA,label=T,label.size = 6)&NoLegend()&NoAxes()&theme(plot.title=element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
tiff("CombinedCellType_myType_NoLabel.tiff",width=450,height=400)
DimPlot(CombineSCRNA,label=F,label.size = 6)&NoLegend()&NoAxes()&theme(plot.title=element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
CombineSCRNA$MyType=Idents(CombineSCRNA)
saveRDS(CombineSCRNA,file="CombineYangAndMathys.rds")


#########################Obtain the FDAMic mediated communication###################
setwd("/Projects/deng/Alzheimer/syn18485175")
CombineSCRNA=readRDS("/data2/deng/Alzheimer/Microglia/Yang_Vascular/CombineYangAndMathys.rds")
typeOrder=c("End","Per","SMC","Fibro","Ependymal","In","Ex","Oli","Opc","Ast","Mic","T")
CombineSCRNA$MyType=factor(CombineSCRNA$MyType,levels=typeOrder)
colorList=colorRampPalette(RColorBrewer::brewer.pal(9,name = 'Paired'))(12)
tiff("Graph/Part3/CombineSCRNACellType.tiff",width=450,height=400)
DimPlot(CombineSCRNA,group.by="MyType",cols=colorList)&NoLegend()&NoAxes()&theme(plot.title=element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

CombineSCRNAExpr=AverageExpression(CombineSCRNA,group.by ="MyType")[["RNA"]]
CombineSCRNAExpr=data.frame(CombineSCRNAExpr)
CombineSCRNAExpr$Symbol=rownames(CombineSCRNAExpr)

maxGene <- function(x){
  ifelse(max(x)==0,"Unsign",names(which.max(x)))
}

setwd("/Projects/deng/Alzheimer/syn18485175")
ADMic=readRDS("ADMic.rds")
cellType=c("BAMic","HomMic","NorARMic","ARMic","FDAMic","DysMic")
ADMic=subset(ADMic,SubType%in%cellType)
MicSubtypeExpr=AverageExpression(ADMic)[["RNA"]]

ADdeg <- wilcoxauc(ADMic, 'SubType')
LR_pair=read.table("/Projects/deng/Public/CellCommunication/CellTalkDB_human_lr_pair.txt",header=T)
LigandDEG=ADdeg[ADdeg$feature%in%LR_pair$ligand_gene_symbol&ADdeg$pval<0.01,]
LigandDEG=LigandDEG[LigandDEG$pct_in>20|LigandDEG$pct_out>20,] #at least expressed in 20%
length(unique(LigandDEG$feature)) #39 sig and highly expressed ligands
LigandDEG.df=MicSubtypeExpr[unique(LigandDEG$feature),]
pdf("Graph/Part3/LigandDEG.p20.df.pdf",width=4,height=7) #
t=pheatmap(LigandDEG.df,show_rownames=T,clustering_method="ward.D2",scale="row",border=NA,cluster_col=F,color = viridis(10))
dev.off()

geneAnno=cutree(t$tree_row,k=5)
geneAnno=data.frame(cutree(t$tree_row,k=5))
geneAnno$Group=paste0("G",geneAnno[,1])
colnames(geneAnno)=c("Cluster","Group")
#G1 FDAMic, G2:DysMic, G3:ARMic, G4:HomMic, G5:BAMic
ann_colors = list(
    Group = c(G1="Violet",G2="Gold",G3="SlateBlue",G4="NavajoWhite", G5="forestgreen")
)
pdf("Graph/Part3/LigandDEG.p20.Group.pdf",width=6,height=7) #
t=pheatmap(LigandDEG.df,annotation_row=geneAnno,annotation_colors = ann_colors,clustering_method="ward.D2",show_rownames=T,scale="row",border=NA,cluster_col=F,color = colorRampPalette(c("navy", "white","firebrick3"))(50))
dev.off()


#obtain the receptors for these 39 ligands
geneTree=t$tree_row
geneOrder=geneTree$labels[geneTree$order]
LR_Pair_FromSigLigand=LR_pair[LR_pair$ligand_gene_symbol%in%geneOrder,c("lr_pair","ligand_gene_symbol","receptor_gene_symbol")]
dim(LR_Pair_FromSigLigand) #292 LR_pair generated
#merge LR_repair and gene expression matrix by receptor to explore their expression across the brain cell types
ReceptorExpr_FromSigLigand=merge(LR_Pair_FromSigLigand,CombineSCRNAExpr,by.x="receptor_gene_symbol",by.y="Symbol")
#order receptor expression matrix according the order of Ligand in microglia heatmap
ReceptorExpr_FromSigLigand$ligand_gene_symbol=factor(ReceptorExpr_FromSigLigand$ligand_gene_symbol,levels=geneOrder)
ReceptorExpr_FromSigLigand=ReceptorExpr_FromSigLigand[order(ReceptorExpr_FromSigLigand$ligand_gene_symbol),]
#obtain the expression matrix
rownames(ReceptorExpr_FromSigLigand)=ReceptorExpr_FromSigLigand$lr_pair
ReceptorExpr_FromSigLigand_Matrix=ReceptorExpr_FromSigLigand[,levels(CombineSCRNA)] #only extract the cell type from brain cell type
index=MatrixGenerics::rowMaxs(as.matrix(ReceptorExpr_FromSigLigand_Matrix))>1 #mark the low expression receptors
ReceptorExpr_FromSigLigand_Matrix=ReceptorExpr_FromSigLigand_Matrix[index,]
dim(ReceptorExpr_FromSigLigand_Matrix)
#138 receptors were obtainned

MarkerGroup=apply(ReceptorExpr_FromSigLigand_Matrix, 1, maxGene)
ReceptorExpr_FromSigLigand=ReceptorExpr_FromSigLigand[ReceptorExpr_FromSigLigand$lr_pair%in%names(MarkerGroup),]
all(names(MarkerGroup)==ReceptorExpr_FromSigLigand$lr_pair)
ReceptorExpr_FromSigLigand$receptorType=MarkerGroup

geneAnno$symbol=rownames(geneAnno)
ReceptorExpr_FromSigLigand_Type=merge(ReceptorExpr_FromSigLigand,geneAnno,by.x="ligand_gene_symbol",by.y="symbol")
write.table(ReceptorExpr_FromSigLigand_Type,file="Graph/Part3/ReceptorTypeExpr4sigLigand.txt",quote=F,sep="\t")
recetproCellTypeBy_sigLigandGroup=data.frame(table(ReceptorExpr_FromSigLigand_Type$Group,ReceptorExpr_FromSigLigand_Type$receptorType))
colnames(recetproCellTypeBy_sigLigandGroup)=c("LigandGroup","ReceptorCellType","Number")
recetproCellTypeBy_sigLigandGroup=recetproCellTypeBy_sigLigandGroup[order(recetproCellTypeBy_sigLigandGroup$Number,decreasing=T),]
recetproCellTypeBy_sigLigandGroup$ReceptorCellType=factor(recetproCellTypeBy_sigLigandGroup$ReceptorCellType,levels=rev(unique(recetproCellTypeBy_sigLigandGroup$ReceptorCellType)))
t=ggplot(recetproCellTypeBy_sigLigandGroup, aes(ReceptorCellType,Number, fill=LigandGroup)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=rev(c("Violet", "Gold", "SlateBlue","NavajoWhite","forestgreen")))+
  #scale_y_continuous(expand=c(0,0))+
  theme_bw()+coord_flip()+
  theme(axis.text.x = element_text(angle = 0),axis.title.x = element_blank(),axis.title.y = element_blank())
pdf("Graph/Part3/ReceptorCellTypeNumber4sigLigandGroup.pdf",width=4,height=4)
print(t)
dev.off()
#G1 FDAMic, G2:DysMic, G3:ARMic, G4:HomMic, G5:BAMic

typeOrder=c("End","Per","SMC","Fibro","Ependymal","In","Ex","Oli","Opc","Ast","Mic","T")
ReceptorExpr_FromSigLigand_G1=ReceptorExpr_FromSigLigand[ReceptorExpr_FromSigLigand$ligand_gene_symbol%in%rownames(geneAnno[geneAnno$Group=="G1",]),]
ReceptorExpr_FromSigLigand_G1.df=reshape2::melt(ReceptorExpr_FromSigLigand_G1,id=c(1:3,16))
colnames(ReceptorExpr_FromSigLigand_G1.df)=c("receptor_gene_symbol","lr_pair","ligand_gene_symbol","receptorType","BrainType","ExprInBrainType")
ReceptorExpr_FromSigLigand_G1.df$BrainType=factor(ReceptorExpr_FromSigLigand_G1.df$BrainType,levels=typeOrder)
t=ggplot(ReceptorExpr_FromSigLigand_G1.df, aes(factor(lr_pair,levels=rev(unique(lr_pair))),ExprInBrainType, fill=BrainType)) +
  geom_bar(stat="identity",position="fill") +
  scale_fill_manual(values=colorRampPalette(RColorBrewer::brewer.pal(9,name = 'Paired'))(12))+
  #scale_y_continuous(expand=c(0,0))+
  theme_bw()+coord_flip()+
  theme(axis.text.x = element_text(angle = 0),axis.title.x = element_blank(),axis.title.y = element_blank())
pdf("Graph/Part3/ReceptorCellTypeExpression4sigLigandGroup_G3.pdf",width=6,height=4)
print(t)
dev.off()

ReceptorExpr_FromSigLigand_G1$receptorType=factor(ReceptorExpr_FromSigLigand_G1$receptorType,levels=typeOrder)
t=ggplot(ReceptorExpr_FromSigLigand_G1, aes(factor(lr_pair,levels=rev(unique(lr_pair))),1, fill=receptorType)) +
  geom_bar(stat="identity",position="fill") +
  scale_fill_manual(values=colorRampPalette(RColorBrewer::brewer.pal(9,name = 'Paired'))(12))+
  #scale_y_continuous(expand=c(0,0))+
  theme_bw()+coord_flip()+
  theme(axis.text.x = element_text(angle = 0),axis.title.x = element_blank(),axis.title.y = element_blank())
pdf("Graph/Part3/ReceptorCellTypeExpression4sigLigandGroup_G1_label.pdf",width=3,height=12)
print(t)
dev.off()


setwd("D:/Alzheimer/syn18485175")
ReceptorTypeExpr4sigLigand=read.table("Graph/Part3/ReceptorTypeExpr4sigLigand.txt",header=T)
RecepterLostInFDAMicFromTcell=ReceptorTypeExpr4sigLigand[ReceptorTypeExpr4sigLigand$receptorType=="T"&ReceptorTypeExpr4sigLigand$Group=="G_1",]
geneId=bitr(RecepterLostInFDAMicFromTcell$receptor_gene_symbol, "SYMBOL", "ENTREZID", OrgDb = "org.Hs.eg.db", drop = TRUE)
dim(geneId) #8
keggPaths=enrichKEGG(
  geneId$ENTREZID,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)
dim(data.frame(keggPaths))

#GeneFC=RecepterLostInFDAMicFromTcell$End/rowMeans(RecepterLostInFDAMicFromTcell[,c(5:15)])
GeneFC=log2(RecepterLostInFDAMicFromTcell$T/rowMeans(RecepterLostInFDAMicFromTcell[,c(4:14)]))
names(GeneFC)=RecepterLostInFDAMicFromTcell$receptor_gene_symbol
keggPathsSymbol <- setReadable(keggPaths, 'org.Hs.eg.db', 'ENTREZID')
pdf("Graph/Part3/KEGGPathwayBy_G1_Recepter_T_centplot.pdf",height=7,width=7)
cnetplot(keggPathsSymbol, 
  showCategory = 10,
  color.params = list(foldChange = GeneFC,category = "#EE82EE"),
  cex.params = list(category_node = 1.5, gene_node = 1, category_label = 1, gene_label = 1)
  )
dev.off()

write.table(data.frame(keggPathsSymbol),file="Graph/Part3/KEGGPathwayBy_G1_Recepter_T.txt",sep="\t",quote=F)

keggPaths=data.frame(keggPaths[1:10,])
keggPaths$logP=-log10(keggPaths$pvalue)
keggPaths$Description=factor(keggPaths$Description,levels=rev(keggPaths$Description))
pdf("Graph/Part3/KEGGPathwayBy_G1_Recepter_Mic.pdf",height=3,width=7)
ggplot(data=keggPaths, aes(x=Description, y=Count, fill=logP)) +
  geom_bar(stat="identity")+
  coord_flip()+
  scale_fill_gradient(low = "#472FED",high = "#EE82EE")+
  theme_minimal()
dev.off()

colors_list <- c("forestgreen","NavajoWhite","LightBlue", "SlateBlue", "Violet", "gold")
pdf("Graph/Part3/FDAMicLigand_VlnPlot.pdf",height=8,width=4.5)
scCustomize::Stacked_VlnPlot(seurat_object = ADMic, features = c("B2M","C1QA","C1QB","CALM1","APOC1","PKM","GNAI2","HLA-E","PSAP","SPP1","CALR","TGFB1","APOE","HLA-C","HLA-A","HLA-B"), x_lab_rotate = TRUE, plot_spacing = 0.05,colors_use = colors_list)&
theme(axis.title.x=element_blank())
dev.off()

receptor=c("CD4","TLR2","TGFBR1","TLR1","TREM2","TYROBP","CD74","CD63","PTPRC","ADCY7","ITGB2","CANX","ADGRG1","MSN","LRP4","LINGO1","CSF1R","RPSA","CD81")
pdf("Graph/Part3/FDAMicReceptor_VlnPlot.pdf",height=8,width=4.5)
scCustomize::Stacked_VlnPlot(seurat_object = ADMic, features = receptor, x_lab_rotate = TRUE, plot_spacing = 0.05,colors_use = colors_list)&
theme(axis.title.x=element_blank())
dev.off()




ADdeg <- wilcoxauc(ADMic, 'SubType')
LR_pair=read.table("/Projects/deng/Public/CellCommunication/CellTalkDB_human_lr_pair.txt",header=T)
ReceptorDEG=ADdeg[ADdeg$feature%in%LR_pair$receptor_gene_symbol&ADdeg$pval<0.01,]
ReceptorDEG=ReceptorDEG[ReceptorDEG$pct_in>20|ReceptorDEG$pct_out>20,] #at least expressed in 20%
length(unique(ReceptorDEG$feature)) #57 sig and highly expressed ligands
ReceptorDEG.df=MicSubtypeExpr[unique(ReceptorDEG$feature),]
pdf("Graph/Part3/ReceptorDEG.p20.df.pdf",width=4,height=7) #
t=pheatmap(ReceptorDEG.df,show_rownames=T,clustering_method="ward.D2",scale="row",border=NA,cluster_col=F,color = viridis(10))
dev.off()

geneAnno=cutree(t$tree_row,k=5)
geneAnno=data.frame(cutree(t$tree_row,k=5))
geneAnno$Group=paste0("G",geneAnno[,1])
colnames(geneAnno)=c("Cluster","Group")
#G1 ARMic, G2:NorARMic, G3:FDAMic, G4:HomMic, G5:BAMic
ann_colors = list(
    Group = c(G1="SlateBlue",G2="LightBlue",G3="Violet",G4="NavajoWhite", G5="forestgreen")
)
pdf("Graph/Part3/ReceptorDEG.p20.df.Group.pdf",width=6,height=9) #
t=pheatmap(ReceptorDEG.df,annotation_row=geneAnno,annotation_colors = ann_colors,clustering_method="ward.D2",show_rownames=T,scale="row",border=NA,cluster_col=F,color = colorRampPalette(c("navy", "white","firebrick3"))(50))
dev.off()

geneTree=t$tree_row
geneOrder=geneTree$labels[geneTree$order]
LR_Pair_FromSigReceptor=LR_pair[LR_pair$receptor_gene_symbol%in%geneOrder,c("lr_pair","ligand_gene_symbol","receptor_gene_symbol")]
dim(LR_Pair_FromSigReceptor) #264 LR_pair generated
#merge LR_repair and gene expression matrix by ligand to explore their expression across the brain cell types
LigandExpr_FromSigReceptor=merge(LR_Pair_FromSigReceptor,CombineSCRNAExpr,by.x="ligand_gene_symbol",by.y="Symbol")
#obtain the expression matrix
rownames(LigandExpr_FromSigReceptor)=LigandExpr_FromSigReceptor$lr_pair
LigandExpr_FromSigReceptor_Matrix=LigandExpr_FromSigReceptor[,levels(CombineSCRNA)] #only extract the cell type from brain cell type
index=MatrixGenerics::rowMaxs(as.matrix(LigandExpr_FromSigReceptor_Matrix))>1 #mark the low expression receptors
LigandExpr_FromSigReceptor_Matrix=LigandExpr_FromSigReceptor_Matrix[index,]
dim(LigandExpr_FromSigReceptor_Matrix)
#108 receptors were obtainned

MarkerGroup=apply(LigandExpr_FromSigReceptor_Matrix, 1, maxGene)
LigandExpr_FromSigReceptor=LigandExpr_FromSigReceptor[LigandExpr_FromSigReceptor$lr_pair%in%names(MarkerGroup),]
all(names(MarkerGroup)==LigandExpr_FromSigReceptor$lr_pair)
LigandExpr_FromSigReceptor$ligandType=MarkerGroup
length(unique(LigandExpr_FromSigReceptor$receptor_gene_symbol))
#37 recetpors left from 57 highly expressed in microglia (20 ligands were not mapped any receptor, might be removed by their lower expression)

geneAnno$symbol=rownames(geneAnno)
LigandExpr_FromSigReceptor_Type=merge(LigandExpr_FromSigReceptor,geneAnno,by.x="receptor_gene_symbol",by.y="symbol")
write.table(LigandExpr_FromSigReceptor_Type,file="Graph/Part3/LigandTypeExpr4sigReceptor.txt",quote=F,sep="\t")
LigandCellTypeBy_sigReceptorGroup=data.frame(table(LigandExpr_FromSigReceptor_Type$Group,LigandExpr_FromSigReceptor_Type$ligandType))
colnames(LigandCellTypeBy_sigReceptorGroup)=c("ReceptorGroup","LigandCellType","Number")
G3=LigandCellTypeBy_sigReceptorGroup[LigandCellTypeBy_sigReceptorGroup$ReceptorGroup=="G3",]
G3=G3[order(G3$Number,decreasing=T),]
LigandCellTypeBy_sigReceptorGroup$ReceptorGroup=factor(LigandCellTypeBy_sigReceptorGroup$ReceptorGroup,levels=rev(c("G3","G1","G2","G4","G5")))
LigandCellTypeBy_sigReceptorGroup$LigandCellType=factor(LigandCellTypeBy_sigReceptorGroup$LigandCellType,levels=rev(unique(G3$LigandCellType)))
t=ggplot(LigandCellTypeBy_sigReceptorGroup, aes(LigandCellType,Number, fill=ReceptorGroup)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=rev(c("Violet", "SlateBlue", "LightBlue","NavajoWhite","forestgreen")))+
  #scale_y_continuous(expand=c(0,0))+
  theme_bw()+coord_flip()+
  theme(axis.text.x = element_text(angle = 0),axis.title.x = element_blank(),axis.title.y = element_blank())
pdf("Graph/Part3/LigandCellTypeNumber4sigReceptorGroup.pdf",width=4,height=4)
print(t)
dev.off()

LigandTypeExpr4sigReceptor=read.table("Graph/Part3/LigandTypeExpr4sigReceptor.txt",header=T)
LigandLostInFDAMicFromEx=LigandTypeExpr4sigReceptor[LigandTypeExpr4sigReceptor$ligandType=="End"&LigandTypeExpr4sigReceptor$Group=="G2",]
geneId=bitr(LigandLostInFDAMicFromEx$ligand_gene_symbol, "SYMBOL", "ENTREZID", OrgDb = "org.Hs.eg.db", drop = TRUE)
keggPaths=enrichKEGG(
  geneId$ENTREZID,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.2,
  use_internal_data = FALSE
)
dim(data.frame(keggPaths))
keggPaths=data.frame(keggPaths[1:10,])
keggPaths$logP=-log10(keggPaths$pvalue)
keggPaths$Description=factor(keggPaths$Description,levels=rev(keggPaths$Description))
pdf("Graph/Part3/KEGGPathwayBy_G2_Ligand_End.pdf",height=3,width=7)
ggplot(data=keggPaths, aes(x=Description, y=Count, fill=logP)) +
  geom_bar(stat="identity")+
  coord_flip()+
  scale_fill_gradient(low = "LightSteelBlue",high = "BlueViolet")+
  theme_minimal()
dev.off()
