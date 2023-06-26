setwd("/data2/deng/Alzheimer")
library(ggplot2)
library(Seurat)

load("Mathylation_syn3157275/ProcessedImputed.RData")
GeneExpr=read.table("Gene_Syn3388564/geneExpr.txt",header=TRUE,row.names=1,sep="\t",check.names=F)
sample=intersect(colnames(imputed),colnames(GeneExpr))
imputed=imputed[,sample]
GeneExpr=GeneExpr[,sample]
SampleInfo=imputed.clinical.group[sample,]
all(colnames(imputed)==rownames(SampleInfo)) #TRUE

targetGene=c("HLA-DRA","HLA-DRA","CD74","HLA-DRA","SPI1","CD14")
targetGene=c("DNMT1","DNMT3A","DNMT3B","DNMT3L") TET
all(colnames(GeneExpr)==rownames(SampleInfo)) #TRUE
targetGeneInfo=data.frame(t(GeneExpr[targetGene,]),SampleInfo)
targetGeneInfo=targetGeneInfo[,c("SexCondition",targetGene)]
targetGeneInfo.long=reshape2::melt(targetGeneInfo,id=1)
colnames(targetGeneInfo.long)=c("Condition","Gene","Expr")
targetGeneInfo.long=targetGeneInfo.long[targetGeneInfo.long$Condition%in%c("Female_G1","Female_G2","Male_G1","Male_G2"),]
t=ggplot(targetGeneInfo.long, aes(x=Condition, y=log2(Expr+1), color=Condition))+
geom_boxplot()+labs(title="",x="", y = "")+ #ylim(0,0.2)+
geom_jitter(width=0.25, alpha=0.5)+
facet_wrap(Gene~.,scale="free_y",ncol=3)+
theme (strip.text= element_text(size=10,face="bold")) + 
theme(legend.position="none")+ theme_bw()+
theme(axis.title.y = element_text(size=rel(1.0)),axis.text.x = element_text(angle = 90,vjust=0.3,hjust=1),axis.text.y = element_text(size=rel(1.0)))
pdf("Mathylation_syn3157275/CorBetMethyGene/DNAMethyltransferase_Expr.pdf",height=6,width=10)
print(t)
dev.off()


targetGene="SPI1"
all(colnames(GeneExpr)==rownames(SampleInfo)) #TRUE
targetGeneInfo=data.frame(t(GeneExpr[targetGene,]),SampleInfo)
targetGeneInfoTmp=targetGeneInfo[targetGeneInfo$HLA-DRA_genotype%in%c(22,23,33,34,44),]
targetGeneInfoTmp$HLA-DRA=as.character(ifelse(targetGeneInfoTmp$HLA-DRA_genotype%in%c(22,23),"23",ifelse(targetGeneInfoTmp$HLA-DRA_genotype%in%c(34,44),"34","33")))
t=ggplot(targetGeneInfoTmp, aes(x=SexCondition, y=log2(SPI1), color=HLA-DRA))+
geom_boxplot(outlier.size = 0,position = position_dodge(width = .95))+labs(title="",x="", y = "")+ #ylim(0,0.2)+
geom_point(position=position_jitterdodge(jitter.width = 0.2),size=0.3)+
theme (strip.text= element_text(size=10,face="bold")) + 
theme(legend.position="none")+ theme_bw()+
theme(axis.title.y = element_text(size=rel(1.0)),axis.text.x = element_text(angle = 90,vjust=0.3,hjust=1),axis.text.y = element_text(size=rel(1.0)))
pdf("Mathylation_syn3157275/CorBetMethyGene/SPI1_Expr_HLA-DRAGenotype.pdf",height=3.5,width=5)
print(t)
dev.off()



#######in scRNAseq from Mathys################
targetGene=c("HLA-DRA","HLA-DRA","CD74","HLA-DRA","SPI1","CD14")
targetGene=c("HLA-DRA","CD74","TYROBP","HLA-DRA","SPI1","VSIG4")
AD=readRDS("D:/Alzheimer/syn18485175/AD.rds")
DotPlot(AD,features=targetGene,group.by="Group")&NoLegend()&RotatedAxis()&
theme(axis.title.y = element_blank(),axis.title.x = element_blank())

ADMic=readRDS("D:/Alzheimer/syn18485175/ADMic.rds")
ADMic$Gender=factor(ADMic$Gender,levels=c("Female","Male"))


targetGene=c("C1QA","HLA-DRA","C1QC","VSIG4","CYBA")
#targetGene=c("CYBA","LAPTM5","C1QC","TYROBP","HLA-DRA","NCF4","ARHGAP30","FCER1G","HCK","STXBP2","LAIR1","C1QA","SASH3","SNX20","PTPN6","MYO1F")
VlnPlot(ADMic,features=targetGene,group.by="MicCluster",split.by="Gender",pt.size=0)&NoLegend()&
scale_fill_manual(values=c("Violet","RoyalBlue"))&
theme(axis.title.y = element_blank(),axis.title.x = element_blank())





#######loci in SPI1################
setwd("/data2/deng/Alzheimer/Microglia")
DML_Female_G2vsG1=read.table("Mathylation_syn3157275/DML/DML8limma_Female_G2vsG1.txt",header=T)
load("Mathylation_syn3157275/ProcessedImputed.RData")
GeneExpr=read.table("Gene_Syn3388564/geneExpr.txt",header=TRUE,row.names=1,sep="\t",check.names=F)
sample=intersect(colnames(imputed),colnames(GeneExpr))
imputed=imputed[,sample]
GeneExpr=GeneExpr[,sample]
SampleInfo=imputed.clinical.group[sample,]
all(colnames(imputed)==rownames(SampleInfo)) #TRUE


cgMeta=read.csv("Mathylation_syn3157275/rawData/ROSMAP_arrayMethylation_metaData.tsv",sep="\t")
SPI1Loci=cgMeta[cgMeta$RefGene%in%grep("SPI1",cgMeta$RefGene,value=T),]
SPI1Loci=[order(SPI1Loci$MAPINFO),]
targetCG=SPI1Loci$TargetID

R=array()
value=array()
for(i in 1:length(targetCG)){
	cgLoci=targetCG[i]
	R[i]=cor.test(as.numeric(imputed[cgLoci,]),as.numeric(GeneExpr["SPI1",]))$estimate[[1]]
	value[i]=cor.test(as.numeric(imputed[cgLoci,]),as.numeric(GeneExpr["SPI1",]))$p.value
}
R.df=data.frame(R=R,pval=value)
R.df$Id=targetCG
R.df$Id=factor(R.df$Id,levels=rev(R.df$Id))
#add DMC infor between Female AD vx Female Control
targetCGDMLInfo=DML_Female_G2vsG1[targetCG,]
all(R.df$Id==rownames(targetCGDMLInfo))
R.df$DMC.Pval=targetCGDMLInfo$P.Value

p<-ggplot(data=R.df, aes(x=Id, y=-log10(DMC.Pval),fill=R)) +
  geom_bar(stat="identity",width=0.6)+
  theme(axis.title.x = element_blank())+
  scale_fill_gradient2(low = "blue",mid="white",high = "red")+ theme_bw()+
  theme(axis.title.y = element_text(size=rel(1.0)),axis.title.x = element_blank(),axis.text.x = element_text(angle = 90,vjust=0.5,hjust=1),axis.text.y = element_text(size=rel(1.0)))
  #+coord_flip()

pdf("Mathylation_syn3157275/CorBetMethyGene/SPI1_corBetMethyAndExpr.pdf",height=3,width=5)
print(p)
dev.off()


#######Cor between gene and loci################
setwd("/data2/deng/Alzheimer")
library('ggpubr')
all(colnames(GeneExpr)==colnames(imputed))
targetGene=c("SPI1","HLA-DRA")
targetLoci=c("cg24234141","cg01539849","cg06784824","cg02647874","cg07698783","cg15945333","cg15982099","cg07675031","cg03106245","cg14088811","cg10435245","cg06147863","cg03565868","cg16517172","cg03301240","cg04692506","cg24307117")
df=data.frame(t(log2(GeneExpr[targetGene,]+1)),t(imputed[targetLoci,]),check.names=F)
all(rownames(df)==rownames(SampleInfo))
target.df=data.frame(SampleInfo[,c("Group","Gender")],df[,c("SPI1","HLA-DRA")],check.names=F)
target.df$States=ifelse(target.df$Group%in%c("G1","G3"),"Ctrl","AD")

sp <- ggscatter(target.df, x = "SPI1", y = "HLA-DRA",
	size = 1,
	color="Gender",
	shape="States",
	add = "reg.line",
	conf.int = TRUE,
	add.params = list(color = "black",fill = "lightgray"))+theme_bw()+
  scale_color_manual(values=c("Violet","RoyalBlue")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  stat_cor(method = "pearson")
pdf("Mathylation_syn3157275/CorBetMethyGene/SPI1_HLA-DRA.pdf",height=4,width=5)
print(sp)
dev.off()


#######Cor between SPI1 and other genes################

GeneExpr=read.table("Gene_Syn3388564/geneExpr.txt",header=TRUE,row.names=1,sep="\t",check.names=F)
SPI1Cor=matrix(NA,nrow=nrow(GeneExpr),ncol=3,dimnames = list(rownames(GeneExpr),c("Symbol","R","Pval")))
Logp1Expr=log2(GeneExpr+1)
SPI1Expr=as.numeric(GeneExpr["SPI1",])
for(i in 1:nrow(Logp1Expr)){
	singGeneExpr=as.numeric(GeneExpr[i,])
	t=cor.test(SPI1Expr,singGeneExpr)
	SPI1Cor[i,1]=rownames(Logp1Expr)[i]
	SPI1Cor[i,2]=t$estimate[[1]]
	SPI1Cor[i,3]=t$p.value
}
Tmp=data.frame(SPI1Cor)
Tmp=na.omit(Tmp)
Tmp=Tmp[order(as.numeric(Tmp$Pval)),]
write.table(Tmp,file="Mathylation_syn3157275/CorBetMethyGene/SPI1CorGene.txt",sep="\t",quote=F)

