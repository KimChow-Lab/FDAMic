#Green et al., https://doi.org/10.1101/2023.03.07.531493
#aged prefrontal cortex from 1.64 million single-nucleus RNA-seq profiles for 424 aging individuals

setwd("D:/Alzheimer/ROSMAP/syn52293417_MIT_ROSMAP_Multiomics/AnalysisByGreen")
##Sample information was downloaded from the Table S1 (discovery cohort)
sampleInfor=read.table("D:/Alzheimer/ROSMAP/syn52293417_MIT_ROSMAP_Multiomics/AnalysisByGreen/SampleInformation.txt",header=T,row.names=1,sep="\t")
head(sampleInfor)

##using cogdx (congnitive impairment), ceradsc (amiloid plaq), braaksc (tau) to define pathology status
Pathology=sampleInfor[,c("cogdx","ceradsc","braaksc")]
t=pheatmap(Pathology,clustering_method="ward.D2",show_rownames=F)
a=data.frame(cutree(t$tree_row,k=2))
anno=data.frame(Group=paste0("G",a[,1],sep=""))
rownames(anno)=rownames(a)
all(rownames(anno)==rownames(sampleInfor))
anno$Gender=ifelse(sampleInfor$msex==0,"Female","Male")
ann_colors = list(
    Group = c(G1="Coral", G2="SlateBlue"),
    Gender=c(Female="Violet",Male="RoyalBlue")
)
t=pheatmap(Pathology,clustering_method="ward.D2",annotation_row=anno,annotation_colors = ann_colors,show_rownames=F,color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

pdf("TraitGroup2DefineMahysDatasetInROSMAP.pdf",height=5,width=4.5)
print(t)
dev.off()

all(rownames(anno)==rownames(sampleInfor))
sampleInfor$StatusCluster=anno$Group
sampleInfor$Status=ifelse(sampleInfor$StatusCluster=="G1","ND","LOAD")
sampleInfor$Gender=ifelse(sampleInfor$msex==0,"Female","Male")
write.table(sampleInfor,file="PhenotypeStatus.txt",sep="\t",quote=F)

sampleInfor$sample=rownames(sampleInfor)



##The microglia ratio for each subtype (identify by Green et al.,)
MicRatio=read.table("MicrogliaRatio.txt",header=T,row.names=1,sep="\t")
MicRatio$sample=rownames(MicRatio)

MicRatio.df=merge(sampleInfor[,c("sample","Gender","Status")],MicRatio,by="sample")
MicRatio.ldf=reshape2::melt(MicRatio.df,id=c(1:3))
colnames(MicRatio.ldf)=c("Sample","Gender","Status","MicCluster","Ratio")
MicRatio.ldf$Group=paste0(MicRatio.ldf$Gender,"_",MicRatio.ldf$Status,sep="")
MicRatio.ldf$Gender=factor(MicRatio.ldf$Gender,levels=c("Male","Female"))
MicRatio.ldf$Status=factor(MicRatio.ldf$Status,levels=c("ND","LOAD"))

t=ggplot(MicRatio.ldf, aes(x=Gender, y=as.numeric(Ratio)*100, fill=Status)) +
  geom_boxplot(position=position_dodge(0.8)) +  #,outlier.shape = NA
  #geom_point(position=position_jitterdodge(),aes(color=Status),size=0.1)+
  scale_fill_manual(values = c("RoyalBlue","Violet"))+
  scale_color_manual(values = c("RoyalBlue","Violet"))+
  theme_bw()+
  #ggpubr::stat_compare_means(method = "t.test",aes(label = sprintf("p = %4.3f", as.numeric(..p.format..))))+
  facet_wrap(.~MicCluster,scales="free_y")+
  scale_y_continuous(expand = c(0.2, 0))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0)))

pdf("cellNumberVariationFromGreen_NoPvalue.pdf",width=10,height=7)
print(t)
dev.off()


## validate the expression of marker for Mic.12 and Mic.13 across microglia subtypes
setwd("/Projects/deng/Alzheimer/syn18485175")
ADMic=readRDS("ADMic.rds")

ADMicTmp=subset(ADMic,SubType%in%c("HomMic","NorARMic", "ARMic","FDAMic","DysMic"))

TargetGenes=c("APOE","TREM2","PPARG")
t=DotPlot(ADMicTmp,features=TargetGenes,group.by="SubType")
g=ggplot(t$data, aes(id,factor(features.plot,levels=rev(unique(features.plot))),size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_colour_viridis_c()+
labs(color="avg.exp,scaled",size="pct.exp",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),angle=90,hjust=1,vjust=0.5),axis.text.y = element_text(size=rel(1.0))) 
pdf("/Projects/deng/Alzheimer/syn18485175/Manuscript/Microglia/Graph/Part1/TargetGenesByGreen.pdf",width=5,height=3) #
print(g)
dev.off()



TargetGenes=c("APOE","TREM2","PPARG","CD74","HLA-DRA","HLA-DRB1","FOXP1","SPP1","TGFB1","SMAD3","ADAM10")
t=DotPlot(ADMicTmp,features=TargetGenes,group.by="SubType")
g=ggplot(t$data, aes(factor(features.plot,levels=rev(unique(features.plot))),id,size= pct.exp,color=avg.exp.scaled)) +geom_point()+
scale_colour_viridis_c()+
labs(color="avg.exp,scaled",size="pct.exp",x="",y="",title="")+
theme_bw()+theme(title = element_text(size=rel(1.0)),axis.text.x = element_text(size=rel(1.0),angle=90,hjust=1,vjust=0.5),axis.text.y = element_text(size=rel(1.0))) 
pdf("/Projects/deng/Alzheimer/syn18485175/Manuscript/Microglia/Graph/Part1/TargetGenesByGreen_Mic12_13.pdf",width=8,height=3) #
print(g)
dev.off()















#http://compbio.mit.edu/microglia_states/
#data was downloaded from https://personal.broadinstitute.org/cboix/sun_victor_et_al_data/ not from syn52293417

ImmuneCells.meta=readRDS("rawData/ROSMAP.ImmuneCells.6regions.snRNAseq.meta.rds")
#TRUE

#only focused on the PFC (PFC batch)
SampleInfor=ImmuneCells.meta[ImmuneCells.meta$brainRegion%in%c("PFC")&ImmuneCells.meta$batch=="PFC",]
sampleClusterNumber=data.frame(table(SampleInfor$subject,SampleInfor$seurat_clusters))
colnames(sampleClusterNumber)=c("sample","cluster","SampleClusterN")
sampleNumber=data.frame(table(SampleInfor$subject))
colnames(sampleNumber)=c("sample","SampleN")
micRatio=merge(sampleClusterNumber,sampleNumber,by="sample")
micRatio$Ratio=micRatio$SampleClusterN/micRatio$SampleN
phenotype=unique(SampleInfor[,c("subject","msex","ADdiag3types")])
colnames(phenotype)=c("sample","Gender","ADdiag3types")
micRatio.df=merge(phenotype,micRatio,by="sample")
micRatio.df$Gender=factor(micRatio.df$Gender,levels=c("Male","Female"))
micRatio.df$ADdiag3types=factor(micRatio.df$ADdiag3types,levels=c("nonAD","earlyAD","lateAD"))
micRatio.df$Status=ifelse(micRatio.df$ADdiag3types=="nonAD","ND","LOAD")
micRatio.df$Status=factor(micRatio.df$Status,levels=c("ND","LOAD"))
t=ggplot(micRatio.df, aes(x=Gender, y=as.numeric(Ratio)*100, fill=Status)) +
  geom_boxplot(position=position_dodge(0.8)) +  #,outlier.shape = NA
  #geom_point(position=position_jitterdodge(),aes(fill=Status))+
  scale_fill_manual(values = c("RoyalBlue","Violet","Red"))+
  theme_bw()+
  #ggpubr::stat_compare_means(method = "t.test",aes(label = sprintf("p = %4.3f", as.numeric(..p.format..))))+
  facet_wrap(.~cluster,scales="free_y")+
  scale_y_continuous(expand = c(0.2, 0))+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0)))
pdf("cellNumberVariationFromNaSun_PFC(PFCBatchOnly)_NoPvalue.pdf",width=10,height=8)
print(t)
dev.off()

