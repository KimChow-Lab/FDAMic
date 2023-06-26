#https://www.nature.com/articles/nn.3786
#hg19 human reference genome
setwd("/data2/deng/Alzheimer/Mathylation_syn3157275/rawData")
library(data.table)
imputed=fread("ROSMAP_arrayMethylation_imputed.tsv", header = TRUE,check.names=F)
imputed=data.frame(imputed,check.names = F)
rownames(imputed)=imputed$TargetID
imputed$TargetID=NULL
#dim(imputed)
#420132    740

IDkey=read.csv("ROSMAP_IDkey.csv")
mwasID=IDkey[match(colnames(imputed),IDkey$mwas_id),c("projid","mwas_id","rnaseq_id","mirna_id","microglia_mrna_id")]
all(colnames(imputed)==mwasID$mwas_id)
#TRUE
colnames(imputed)=mwasID$projid

clinical=read.csv("ROSMAP_clinical.csv",row.names=1) #3584 samples
all(colnames(imputed)%in%rownames(clinical))
#TRUE
projidList=intersect(colnames(imputed),rownames(clinical)) #740 samples
imputed.clinical=clinical[projidList,] #740 samples
imputed.clinical$Age=ifelse(imputed.clinical$age_at_visit_max=="90+","90",imputed.clinical$age_at_visit_max)
imputed.clinical$age_at_visit_max=NULL
#msex:0 female,1 male,
#define the AD imputed.clinical
dim(imputed.clinical)
#740  17

imputed.clinical.tmp=na.omit(imputed.clinical[,c("braaksc","ceradsc","cogdx","dcfdx_lv")]) #739 samples left
dim(imputed.clinical.tmp)
pdf("clinical_739_samples.pdf",width=6,height=3)
t=pheatmap(t(imputed.clinical.tmp[,c("braaksc","ceradsc","cogdx","dcfdx_lv")]),clustering_method="ward.D2",color = colorRampPalette(c("navy", "white", "firebrick3"))(50),show_colnames=F)
dev.off()

imputed.clinical=imputed.clinical[rownames(imputed.clinical.tmp),]
Cluster=data.frame(cutree(t$tree_col,k=4))
colnames(Cluster)="Group"
Cluster$Group=paste0("G",Cluster$Group)
Cluster$Gender=ifelse(imputed.clinical$msex==1,"Male","Female")
annotation_col = list(
                    Group =c(G1 = "#0d5b26", G2 = "#2e5fa1",G3="#3fab47",G4="#52b9d8"),
                    Gender=c(Male="RoyalBlue",Female="Violet")
                )
pdf("clinical.imputed.group.pdf",width=7,height=3)
t=pheatmap(t(imputed.clinical.tmp[,c("braaksc","ceradsc","cogdx","dcfdx_lv")]),annotation_col=Cluster,annotation_colors=annotation_col,clustering_method="ward.D2",color = colorRampPalette(c("navy", "white", "firebrick3"))(50),show_colnames=F)
dev.off()

tmp=data.frame(cbind(sample=rownames(Cluster),imputed.clinical=Cluster))
tmp=tmp[order(tmp$imputed.clinical),]
imputed=imputed[,rownames(tmp)]
dim(imputed)
#420132    739

imputed.clinical=imputed.clinical[rownames(tmp),]
all(rownames(imputed.clinical)==colnames(imputed))

imputed.pca <- prcomp(t(imputed), center = TRUE,scale. = TRUE)
pcaMatrix=data.frame(imputed.pca$x)
pcaMatrixInfo=cbind(pcaMatrix,imputed.clinical,tmp)
pdf("PCA.pdf",width=7,height=6)
ggplot(pcaMatrixInfo, aes(PC1, PC2, color=imputed.clinical)) +
  geom_point(size=3) +
  scale_color_manual(values = c("#0d5b26","#2e5fa1","#3fab47","#52b9d8"))+
  theme_bw()
dev.off()

all(rownames(imputed.clinical)==rownames(tmp))
#739  17
imputed.clinical.group=cbind(imputed.clinical,tmp)
imputed.clinical.group$Age=as.numeric(imputed.clinical.group$Age)
imputed.clinical.group$Group=factor(imputed.clinical.group$imputed.clinical,levels=c("G2","G4","G1","G3"))
imputed.clinical.group$Gender=ifelse(imputed.clinical.group$msex=="1","Male","Female")
imputed.clinical.group$SexCondition=as.factor(paste0(imputed.clinical.group$Gender,"_",imputed.clinical.group$Group,sep=""))
save(imputed, imputed.clinical.group, file = "../ProcessedImputed.RData")


###DMR (methylation region) analysis using limma##############
setwd("/data2/deng/Alzheimer/Microglia/Mathylation_syn3157275/DML")
load("../ProcessedImputed.RData")
all(colnames(imputed)==rownames(imputed.clinical.group))
#TRUE
imputed.clinical.group$Age=as.numeric(imputed.clinical.group$Age)
imputed.clinical.group=imputed.clinical.group[,c("SexCondition","pmi","Study","educ","apoe_genotype","Age")]
imputed.clinical.group=na.omit(imputed.clinical.group)
imputed=imputed[,rownames(imputed.clinical.group)]
all(colnames(imputed)==rownames(imputed.clinical.group))
design <- model.matrix(~0+SexCondition+pmi+Study+educ+apoe_genotype+Age,data=imputed.clinical.group)
dim(design)
colnames(design)=c(levels(imputed.clinical.group$SexCondition),"pmi","Study","Educ","APOE","Age")
#design <- model.matrix(~0+SexCondition,data=imputed.clinical.group)
#colnames(design)=c(levels(imputed.clinical.group$SexCondition))
all(rownames(design)==rownames(imputed.clinical.group)) #TRUE
contrast.matrix <- makeContrasts(p1="Female_G2-Female_G1",p2="Male_G2-Male_G1",p3="Female_G1-Male_G1",p4="Female_G2-Male_G2",levels = design)
fit <- lmFit(imputed, design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit, trend=TRUE)
DMLInfo=topTable(fit, coef=2,n=Inf, adjust="BH")
write.table(DMLInfo,file="DML8limma_Male_G2vsG1.txt",sep="\t",quote=F)
SigInfo=DMLInfo[DMLInfo$P.Value<0.01,]
write.table(SigInfo,file="SigDML8limma_Male_G2vsG1.txt",sep="\t",quote=F)
dim(SigInfo)


DML8limma_Female_G2vsG1=read.table("DML8limma_Female_G2vsG1.txt",header=T)
DML8limma_Male_G2vsG1=read.table("DML8limma_Male_G2vsG1.txt",header=T)

SigDML_Female_G2vsG1=DML8limma_Female_G2vsG1[DML8limma_Female_G2vsG1$P.Value<0.01,]
SigDML_Male_G2vsG1=DML8limma_Male_G2vsG1[DML8limma_Male_G2vsG1$P.Value<0.01,]

cgMeta=read.csv("../rawData/ROSMAP_arrayMethylation_metaData.tsv",sep="\t")

SigDML_Female_G2vsG1$TargetID=rownames(SigDML_Female_G2vsG1)
SigDML_Female_G2vsG1_Symbol=merge(SigDML_Female_G2vsG1,cgMeta,by="TargetID")
SigDML_Female_G2vsG1_Symbol=SigDML_Female_G2vsG1_Symbol[order(SigDML_Female_G2vsG1_Symbol$P.Value),]
write.table(SigDML_Female_G2vsG1_Symbol,file="../DML/SigDML_Female_G2vsG1_Symbol.txt",sep="\t",quote=F)

SigDML_Female_G2vsG1_Symbol=read.csv("../DML/SigDML_Female_G2vsG1_Symbol.txt",header=T,row.names=1,sep="\t")
FemaleADMethylationGene=unique(strsplit(paste0(unique(SigDML_Female_G2vsG1_Symbol$RefGene),collapse = ";"),";")[[1]])
length(FemaleADMethylationGene)
write.table(FemaleADMethylationGene,file="DMLRelatedGene/FemaleADMethylationGene.txt",quote=F,col.names=F,row.names=F)

SigDML_Male_G2vsG1$TargetID=rownames(SigDML_Male_G2vsG1)
SigDML_Male_G2vsG1_Symbol=merge(SigDML_Male_G2vsG1,cgMeta,by="TargetID")
SigDML_Male_G2vsG1_Symbol=SigDML_Male_G2vsG1_Symbol[order(SigDML_Male_G2vsG1_Symbol$P.Value),]
write.table(SigDML_Male_G2vsG1_Symbol,file="../DML/SigDML_Male_G2vsG1_Symbol.txt",sep="\t",quote=F)

HypoInFAD=SigDML_Female_G2vsG1_Symbol[SigDML_Female_G2vsG1_Symbol$logFC<0,]
HyperInFAD=SigDML_Female_G2vsG1_Symbol[SigDML_Female_G2vsG1_Symbol$logFC>0,]
HypoInFADGene=unique(strsplit(paste0(unique(HypoInFAD$RefGene),collapse = ";"),";")[[1]])
HyperInFADGene=unique(strsplit(paste0(unique(HyperInFAD$RefGene),collapse = ";"),";")[[1]])
length(HypoInFADGene)
length(HyperInFADGene)
head(sort(table(HypoInFAD$Gene_Feature),decreasing=T),n=10)
write.table(HypoInFADGene,file="DMLRelatedGene/HypoInFADGene.txt",quote=F,col.names=F,row.names=F)
write.table(HyperInFADGene,file="DMLRelatedGene/HyperInFADGene.txt",quote=F,col.names=F,row.names=F)


HypoInFAD_Body=HypoInFAD[HypoInFAD$Gene_Feature%in%grep("Body",HypoInFAD$Gene_Feature,value=T),]
head(sort(table(HypoInFAD_Body$Gene_Feature),decreasing=T),n=10)
HypoInFAD_Tss=HypoInFAD[HypoInFAD$Gene_Feature%in%grep("TSS",HypoInFAD$Gene_Feature,value=T),]
head(sort(table(HypoInFAD_Tss$Gene_Feature),decreasing=T),n=10)
HyperInFAD_Body=HyperInFAD[HyperInFAD$Gene_Feature%in%grep("Body",HyperInFAD$Gene_Feature,value=T),]
head(sort(table(HyperInFAD_Body$Gene_Feature),decreasing=T),n=10)
HyperInFAD_Tss=HyperInFAD[HyperInFAD$Gene_Feature%in%grep("TSS",HyperInFAD$Gene_Feature,value=T),]
head(sort(table(HyperInFAD_Tss$Gene_Feature),decreasing=T),n=10)
BodyHypoInFAD=unique(strsplit(paste0(unique(HypoInFAD_Body$RefGene),collapse = ";"),";")[[1]])
TssHypoInFAD=unique(strsplit(paste0(unique(HypoInFAD_Tss$RefGene),collapse = ";"),";")[[1]])
BodyHyperInFAD=unique(strsplit(paste0(unique(HyperInFAD_Body$RefGene),collapse = ";"),";")[[1]])
TssHyperInFAD=unique(strsplit(paste0(unique(HyperInFAD_Tss$RefGene),collapse = ";"),";")[[1]])


HypoInMAD=SigDML_Male_G2vsG1_Symbol[SigDML_Male_G2vsG1_Symbol$logFC<0,]
HyperInMAD=SigDML_Male_G2vsG1_Symbol[SigDML_Male_G2vsG1_Symbol$logFC>0,]
HypoInMADGene=unique(strsplit(paste0(unique(HypoInMAD$RefGene),collapse = ";"),";")[[1]])
HyperInMADGene=unique(strsplit(paste0(unique(HyperInMAD$RefGene),collapse = ";"),";")[[1]])
write.table(HypoInMADGene,file="DMLRelatedGene/HypoInMADGene.txt",quote=F,col.names=F,row.names=F)
write.table(HyperInMADGene,file="DMLRelatedGene/HyperInMADGene.txt",quote=F,col.names=F,row.names=F)

HypoInMAD_Body=HypoInMAD[HypoInMAD$Gene_Feature%in%grep("Body",HypoInMAD$Gene_Feature,value=T),]
head(sort(table(HypoInMAD_Body$Gene_Feature),decreasing=T),n=10)
HypoInMAD_Tss=HypoInMAD[HypoInMAD$Gene_Feature%in%grep("TSS",HypoInMAD$Gene_Feature,value=T),]
head(sort(table(HypoInMAD_Tss$Gene_Feature),decreasing=T),n=10)
HyperInMAD_Body=HyperInMAD[HyperInMAD$Gene_Feature%in%grep("Body",HyperInMAD$Gene_Feature,value=T),]
head(sort(table(HyperInMAD_Body$Gene_Feature),decreasing=T),n=10)
HyperInMAD_Tss=HyperInMAD[HyperInMAD$Gene_Feature%in%grep("TSS",HyperInMAD$Gene_Feature,value=T),]
head(sort(table(HyperInMAD_Tss$Gene_Feature),decreasing=T),n=10)
HypoInMAD_Body_Gene=unique(strsplit(paste0(unique(HypoInMAD_Body$RefGene),collapse = ";"),";")[[1]])
HypoInMAD_Tss_Gene=unique(strsplit(paste0(unique(HypoInMAD_Tss$RefGene),collapse = ";"),";")[[1]])
HyperInMAD_Body_Gene=unique(strsplit(paste0(unique(HyperInMAD_Body$RefGene),collapse = ";"),";")[[1]])
HyperInMAD_Tss_Gene=unique(strsplit(paste0(unique(HyperInMAD_Tss$RefGene),collapse = ";"),";")[[1]])


DMLList=c(BodyHypoInFAD=list(BodyHypoInFAD),TssHypoInFAD=list(TssHypoInFAD),
	      BodyHyperInFAD=list(BodyHyperInFAD),TssHyperInFAD=list(TssHyperInFAD),
	      BodyHypoInMAD=list(HypoInMAD_Body_Gene),TssHypoInMAD=list(HypoInMAD_Tss_Gene),
	      BodyHyperInMAD=list(HyperInMAD_Body_Gene),TssHyperInMAD=list(HyperInMAD_Tss_Gene)
	      )
t=upset(fromList(DMLList),
  nsets = 8,
  keep.order = T,
  order.by = c("freq", "degree")
  )
pdf("overlap.pdf",width=12)
print(t)
dev.off()

t=rbind(paste0(DMLList$BodyHypoInFAD, collapse=","),paste0(DMLList$TssHypoInFAD, collapse=","),paste0(DMLList$BodyHyperInFAD, collapse=","),paste0(DMLList$TssHyperInFAD, collapse=","))
rownames(t)=c("BodyHypoInFAD", "TssHypoInFAD", "BodyHyperInFAD","TssHyperInFAD")
write.table(t,file="DMCInFAD.txt",col.names=F,quote=F,sep="\t")




##### cell type specific ##############
Brain_scRNAMathys=readRDS("/Projects/deng/Alzheimer/syn18485175/AD.rds")
library(presto)
library(dplyr)
library(tibble)
library(fgsea)
BrinTypeMarker <- wilcoxauc(Brain_scRNAMathys, 'OldCellType')
fgsea_sets=DMLList
for(cluster in unique(Brain_scRNAMathys$OldCellType)){
print (cluster)
singleTypeMarker <- BrinTypeMarker %>% dplyr::filter(group == cluster) %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
ranks<- deframe(singleTypeMarker)
fgseaRes<- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0)
fwrite(fgseaRes, file=paste0("CellTypeSpecific/",cluster,".txt",sep=""), sep="\t", sep2=c("", " ", ""))
}


data=read.table("D:/Alzheimer/Microglia/Mathylation_syn3157275/DML/CellTypeSpecificDML.txt",header=T,sep="\t")
data=data[,1:8]
data$sig=ifelse(data$pval<0.05,"Sig","NonSig")
PathwayList=factor(data$pathway,levels=rev(c("AgingGene8Peter","DnInLLI","Dn_G6","Dn_G7","Dn_G17","Dn_G4","Dn_G31","Dn_G21","Dn_G9","Dn_G14","Dn_G34","UpInLLI","Up_G5","Up_G10","Up_G16","Up_G22","Up_G28","Up_G13")))
ClusterList=factor(data$Cluster,levels=c("CD14 Mono ","CD16 Mono ","Doublets","Neutro Mono transitional cell","DC","pDC","Neutrophil","Developing neutrophil","Basophil","Eosinophils","Naive CD4 T","CD4 CTLs","CD4 Treg","MAIT","CD8 T","Platelet","Red blood","NK","CD265 NK","Prolif lymph","Naive B","Memory B","IgG B"))
t=ggplot(data,aes(Cluster,pathway,size=-1*log(pval),colour=NES,shape=sig))+geom_point()+
scale_colour_viridis_c()+
scale_shape_manual(values=c(3,19))+
theme_bw()+#theme(legend.position="bottom")+
theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(size=rel(1.0),colour = "black",angle=90, vjust = 0.5, hjust=1),axis.text.y = element_text(size=rel(1.0),colour = "black"))
pdf("CellTypeSpecificDML.pdf",height=5,width=7)
print(t)
dev.off()

setwd("/data2/deng/Alzheimer/Mathylation_syn3157275")
DML_Female_G2vsG1=read.table("DML/DML8limma_Female_G2vsG1.txt",header=T)
DML_Male_G2vsG1=read.table("DML/DML8limma_Male_G2vsG1.txt",header=T)
cgMeta=read.csv("rawData/ROSMAP_arrayMethylation_metaData.tsv",sep="\t")
DML_Female_G2vsG1$TargetID=rownames(DML_Female_G2vsG1)
DML_Female_G2vsG1_Symbol=merge(DML_Female_G2vsG1,cgMeta,by="TargetID")
DML_Female_G2vsG1_Symbol=DML_Female_G2vsG1_Symbol[order(DML_Female_G2vsG1_Symbol$P.Value),]
write.table(DML_Female_G2vsG1_Symbol,file="DML/DML_Female_G2vsG1_Symbol.txt",sep="\t",quote=F)
DML_Male_G2vsG1$TargetID=rownames(DML_Male_G2vsG1)
DML_Male_G2vsG1_Symbol=merge(DML_Male_G2vsG1,cgMeta,by="TargetID")
DML_Male_G2vsG1_Symbol=DML_Male_G2vsG1_Symbol[order(DML_Male_G2vsG1_Symbol$P.Value),]
write.table(DML_Male_G2vsG1_Symbol,file="DML/DML_Male_G2vsG1_Symbol.txt",sep="\t",quote=F)

###Visualization for some specific loci##############
setwd("/data2/deng/Alzheimer/Mathylation_syn3157275")
load("ProcessedImputed.RData")
all(colnames(imputed)==rownames(imputed.clinical.group))
targetPhenotype=c("apoe_genotype","braaksc","ceradsc","cogdx","Gender","SexCondition","Group")
targetPhenotypeInfo=imputed.clinical.group[,targetPhenotype]

DML_Female_G2vsG1=read.table("DML/DML8limma_Female_G2vsG1.txt",header=T,sep="\t")
DML_Male_G2vsG1=read.table("DML/DML8limma_Male_G2vsG1.txt",header=T,sep="\t")

DML_Female_G2vsG1_Symbol=read.table("DML/DML_Female_G2vsG1_Symbol.txt",header=T,sep="\t")
targetGene="SPP1"
tmp=DML_Female_G2vsG1_Symbol[DML_Female_G2vsG1_Symbol$RefGene%in%grep(targetGene,DML_Female_G2vsG1_Symbol$RefGene,value=T),]
tmp[order(tmp$P.Value),c("TargetID","logFC","P.Value","adj.P.Val","RefGene")]

targetLoci=c("cg00718409","cg02539671","cg04456284","cg10725937","cg20095587","cg25748868","cg01980222")
DML_Female_G2vsG1[targetLoci,]
DML_Male_G2vsG1[targetLoci,]

targetLociValue=imputed[targetLoci,]
GraphData=cbind(targetPhenotypeInfo,t(targetLociValue))
GraphDataLong=reshape2::melt(GraphData,id=(1:7))
colnames(GraphDataLong)=c(colnames(targetPhenotypeInfo),"cgID","betaValue")
t=ggplot(GraphDataLong, aes(x=SexCondition, y=`betaValue`, color=SexCondition))+
geom_boxplot()+labs(title="",x="", y = "")+ #ylim(0,0.2)+
facet_wrap(cgID~.,scales="free_y")+
theme (strip.text= element_text(size=10,face="bold")) + 
theme(legend.position="none")+ geom_jitter(shape=1, position=position_jitter(0.2))+theme_bw()+
theme(axis.title.y = element_text(size=rel(1.5)),axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0)))

pdf("Loci/TREM2_AllLoci.pdf",width=15,height=10)
print(t)
dev.off()

DML_Female_G2vsG1_Symbol[DML_Female_G2vsG1_Symbol$TargetID%in%"cg15460348",]