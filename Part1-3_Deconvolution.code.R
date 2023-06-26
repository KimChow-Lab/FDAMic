library(Seurat)
library(pheatmap)
library(ggplot2)

#-----------Deconvolution by Cibersortx-------------------
setwd("D:/Alzheimer/syn18485175/Cibersortx")
set.seed(12)
ADSmall <- AD[, sample(colnames(AD), size = 5000, replace=F))] #random sample 5000 cells from the single cell data
DimPlot(ADSmall,group.by="cellType") #check the cell type result based on this 5000 cells
data  <-  as.matrix(ADSmall@assays$RNA@counts)
all(colnames(data)==rownames(ADSmall@meta.data))
colnames(data)=ADSmall$cellType 
write.table(data,file="CortexCellTypeSignature.txt",sep="\t",quote=F)

phenotype=read.table("D:/Alzheimer/Syn3388564/Phenotype.txt",header=T,row.names=1,sep="\t")
trait=phenotype[,c(13,15,16,17)]
t=pheatmap(trait,clustering_method="ward.D2",show_rownames=F)
a=cutree(t$tree_row,k=2)
result=data.frame(cbind(phenotype,"Label"=a))
ratio=read.table("CIBERSORTx_Job16_Adjusted.txt",header=T,row.names=1)
sample=intersect(rownames(ratio),rownames(result))
ratio=ratio[sample,]
resultTmp=result[sample,c("msex","braaksc","ceradsc","cogdx","dcfdx_lv","Label")]
all(rownames(ratio)==rownames(resultTmp))
resultAll=data.frame("ID"=rownames(resultTmp),resultTmp,ratio)
write.table(resultAll,file="ratio.txt",quote=F,sep="\t")

ratio=read.table("D:/Alzheimer/syn18485175/Cibersortx/ratio.txt",header=T)
Sex=factor(ratio$msex,levels=c("Male","Female"))
Statues=factor(ratio$Condition,levels=c("Control","Alzheimer"))
t=ggplot(ratio, aes(x=Statues, y=Mic*100, fill=Sex)) +
  geom_boxplot(position=position_dodge(0.8),outlier.shape = NA) +  
  geom_dotplot(binaxis='y', stackdir='centerwhole',  position=position_dodge(0.8),dotsize=0.6,binwidth=0.4)+
  scale_fill_manual(values = c("RoyalBlue","Violet"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0)))
pdf("Graph/Part1/MicCellTypeFromROSMAP8Cibersoftx.pdf",height=4,width=4.5)
print(t)
dev.off()
wilcox.test(ratio[ratio$Group=="Female-Alzheimer","Mic"],ratio[ratio$Group=="Female-Control","Mic"]) #p-value = 0.01871
wilcox.test(ratio[ratio$Group=="Male-Alzheimer","Mic"],ratio[ratio$Group=="Male-Control","Mic"]) #p-value = 0.9837
wilcox.test(ratio[ratio$Group=="Female-Alzheimer","Mic"],ratio[ratio$Group=="Male-Alzheimer","Mic"]) #p-value = 0.06688
wilcox.test(ratio[ratio$Group=="Female-Control","Mic"],ratio[ratio$Group=="Male-Control","Mic"]) #p-value = 0.846



#-----------Deconvolution by Scaden-------------------
ADSmallCount=as.matrix(GetAssayData(ADSmall$RNA, slot = "counts")) 
ADSmallCPM=RelativeCounts(ADSmallCount, scale.factor = 1000000)
ADSmallCPM=as.matrix(ADSmallCPM)
ADSmallTmp=ADSmallCPM[apply(ADSmallCPM,1,sd)>0.05,]
write.table(t(ADSmallCPM),file="Cortex_counts.txt",sep="\t",quote=F) #output the normolized matrix for scaden


conda create -n scaden python=3.6
conda activate scaden
pip install scaden
cd TrainingData


scaden simulate --cells 100 --n_samples 1000
scaden process data.h5ad ../geneExprFromROSMAPR0.txt
scaden train processed.h5ad --model_dir model/
scaden predict ../geneExprFromROSMAPR0.txt --model_dir model/



setwd("D:/Alzheimer/syn18485175/Scaden/TrainingData/")
ratio=read.table("scaden_predictions.txt",header=T,row.names=1)
sample=intersect(rownames(ratio),rownames(result))
ratio=ratio[sample,]
resultTmp=result[sample,c("msex","braaksc","ceradsc","cogdx","dcfdx_lv","Label")]
all(rownames(ratio)==rownames(resultTmp))
resultAll=data.frame("ID"=rownames(resultTmp),resultTmp,ratio)
write.table(resultAll,file="ratio.txt",sep="\t",quote=F)

ratio=read.table("D:/Alzheimer/syn18485175/Scaden/TrainingData/ratio.txt",header=T) #data from ROSMAP, script in D:\Alzheimer\syn18485175\Cibersortx\code.R
Sex=factor(ratio$msex,levels=c("Male","Female"))
Statues=factor(ratio$Condition,levels=c("Control","Alzheimer"))
t=ggplot(ratio, aes(x=Statues, y=Mic*100, fill=Sex)) +
  geom_boxplot(position=position_dodge(0.8),outlier.shape = NA) +  
  geom_dotplot(binaxis='y', stackdir='centerwhole',  position=position_dodge(0.8),dotsize=0.6,binwidth=0.2)+
  scale_fill_manual(values = c("RoyalBlue","Violet"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0)))
pdf("Graph/Part1/MicCellTypeFromROSMAP8Scaden.pdf",height=4,width=4)
print(t)
dev.off()
wilcox.test(ratio[ratio$Group=="Female-Alzheimer","Mic"],ratio[ratio$Group=="Female-Control","Mic"]) #p-value = 0.003703
wilcox.test(ratio[ratio$Group=="Male-Alzheimer","Mic"],ratio[ratio$Group=="Male-Control","Mic"]) #p-value = 0.4656
wilcox.test(ratio[ratio$Group=="Female-Alzheimer","Mic"],ratio[ratio$Group=="Male-Alzheimer","Mic"]) #p-value = 0.001345
wilcox.test(ratio[ratio$Group=="Female-Control","Mic"],ratio[ratio$Group=="Male-Control","Mic"]) #p-value = 0.706


ratio=read.table("D:/Alzheimer/syn18485175/Scaden/TrainingData/ratio.txt",header=T) #data from ROSMAP, script in D:\Alzheimer\syn18485175\Cibersortx\code.R
ratio=ratio[ratio$cogdx%in%c(1,2,4,5),]
dim(ratio)
Sex=factor(ratio$msex,levels=c("Male","Female"))
Phenotype=factor(ratio$cogdx,levels=c(1,2,4,5))
t=ggplot(ratio, aes(x=Phenotype, y=Mic*100, fill=Sex)) +
  geom_boxplot(position=position_dodge(0.8),outlier.shape = NA) +  
  geom_dotplot(binaxis='y', stackdir='centerwhole',  position=position_dodge(0.8),dotsize=0.6,binwidth=0.2)+
  scale_fill_manual(values = c("RoyalBlue","Violet"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0)))
pdf("Graph/Part1/MicCellTypeFromROSMAP8Scaden4Cogdx.pdf",height=4,width=5)
print(t)
dev.off()
wilcox.test(ratio[ratio$cogdx=="5"&ratio$msex=="Female","Mic"],ratio[ratio$cogdx=="5"&ratio$msex=="Male","Mic"]) #p-value = 0.01472
wilcox.test(ratio[ratio$cogdx=="4"&ratio$msex=="Female","Mic"],ratio[ratio$cogdx=="4"&ratio$msex=="Male","Mic"]) #p-value = 0.02134
wilcox.test(ratio[ratio$cogdx=="2"&ratio$msex=="Female","Mic"],ratio[ratio$cogdx=="2"&ratio$msex=="Male","Mic"]) #p-value = 0.6274
wilcox.test(ratio[ratio$cogdx=="1"&ratio$msex=="Female","Mic"],ratio[ratio$cogdx=="1"&ratio$msex=="Male","Mic"]) #p-value = 0.3743
wilcox.test(ratio[ratio$cogdx=="5"&ratio$msex=="Female","Mic"],ratio[ratio$cogdx=="1"&ratio$msex=="Female","Mic"]) #p-value = 0.009597
wilcox.test(ratio[ratio$cogdx=="5"&ratio$msex=="Male","Mic"],ratio[ratio$cogdx=="1"&ratio$msex=="Male","Mic"]) #p-value = 0.5744


ratio=read.table("D:/Alzheimer/syn18485175/Scaden/TrainingData/ratio.txt",header=T) #data from ROSMAP, script in D:\Alzheimer\syn18485175\Cibersortx\code.R
count=as.data.frame.matrix(table(ratio$braaksc,ratio$msex))
candidate=rownames(count[count$Female>10&count$Male>10,])
ratio=ratio[ratio$braaksc%in%candidate,]
dim(ratio)
table(ratio$braaksc)
Sex=factor(ratio$msex,levels=c("Male","Female"))
t=ggplot(ratio, aes(x=as.character(braaksc), y=Mic*100, fill=Sex)) +
  geom_boxplot(position=position_dodge(0.8),outlier.shape = NA) +  
  geom_dotplot(binaxis='y', stackdir='centerwhole',  position=position_dodge(0.8),dotsize=0.6,binwidth=0.2)+
  scale_fill_manual(values = c("RoyalBlue","Violet"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0)))
pdf("Graph/Part1/MicCellTypeFromROSMAP8Scaden4braaksc.pdf",height=4,width=5)
print(t)
dev.off()
wilcox.test(ratio[ratio$braaksc=="5"&ratio$msex=="Female","Mic"],ratio[ratio$braaksc=="5"&ratio$msex=="Male","Mic"]) #p-value = 0.07788
wilcox.test(ratio[ratio$braaksc=="4"&ratio$msex=="Female","Mic"],ratio[ratio$braaksc=="4"&ratio$msex=="Male","Mic"]) #p-value = 0.501
wilcox.test(ratio[ratio$braaksc=="3"&ratio$msex=="Female","Mic"],ratio[ratio$braaksc=="3"&ratio$msex=="Male","Mic"]) #p-value = 0.04163
wilcox.test(ratio[ratio$braaksc=="2"&ratio$msex=="Female","Mic"],ratio[ratio$braaksc=="2"&ratio$msex=="Male","Mic"]) #p-value = 0.6432
wilcox.test(ratio[ratio$braaksc=="1"&ratio$msex=="Female","Mic"],ratio[ratio$braaksc=="1"&ratio$msex=="Male","Mic"]) #p-value = 0.7871
wilcox.test(ratio[ratio$braaksc=="5"&ratio$msex=="Female","Mic"],ratio[ratio$braaksc=="1"&ratio$msex=="Female","Mic"]) #p-value = 0.05829
wilcox.test(ratio[ratio$braaksc=="5"&ratio$msex=="Male","Mic"],ratio[ratio$braaksc=="1"&ratio$msex=="Male","Mic"]) #p-value = 0.7104


ratio=read.table("D:/Alzheimer/syn18485175/Scaden/TrainingData/ratio.txt",header=T) #data from ROSMAP, script in D:\Alzheimer\syn18485175\Cibersortx\code.R
table(ratio$ceradsc)
count=as.data.frame.matrix(table(ratio$ceradsc,ratio$msex))
candidate=rownames(count[count$Female>10&count$Male>10,])
ratio=ratio[ratio$ceradsc%in%candidate,]
dim(ratio)
table(ratio$ceradsc)
Sex=factor(ratio$msex,levels=c("Male","Female"))
t=ggplot(ratio, aes(x=factor(ceradsc,levels=c(4:1)), y=Mic*100, fill=Sex)) +
  geom_boxplot(position=position_dodge(0.8),outlier.shape = NA) +  
  geom_dotplot(binaxis='y', stackdir='centerwhole',  position=position_dodge(0.8),dotsize=0.6,binwidth=0.2)+
  scale_fill_manual(values = c("RoyalBlue","Violet"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0)))
pdf("Graph/Part1/MicCellTypeFromROSMAP8Scaden4ceradsc.pdf",height=4,width=5)
print(t)
dev.off()
wilcox.test(ratio[ratio$ceradsc=="1"&ratio$msex=="Female","Mic"],ratio[ratio$ceradsc=="1"&ratio$msex=="Male","Mic"]) #p-value = 0.1347
wilcox.test(ratio[ratio$ceradsc=="2"&ratio$msex=="Female","Mic"],ratio[ratio$ceradsc=="2"&ratio$msex=="Male","Mic"]) #p-value = 0.1551
wilcox.test(ratio[ratio$ceradsc=="3"&ratio$msex=="Female","Mic"],ratio[ratio$ceradsc=="3"&ratio$msex=="Male","Mic"]) #p-value = 0.4826
wilcox.test(ratio[ratio$ceradsc=="4"&ratio$msex=="Female","Mic"],ratio[ratio$ceradsc=="4"&ratio$msex=="Male","Mic"]) #p-value = 0.3948
wilcox.test(ratio[ratio$ceradsc=="4"&ratio$msex=="Female","Mic"],ratio[ratio$ceradsc=="1"&ratio$msex=="Female","Mic"]) #p-value = 0.04766
wilcox.test(ratio[ratio$ceradsc=="4"&ratio$msex=="Male","Mic"],ratio[ratio$ceradsc=="1"&ratio$msex=="Male","Mic"]) #p-value = 0.7471

ratio=read.table("D:/Alzheimer/syn18485175/Scaden/TrainingData/ratio.txt",header=T) #data from ROSMAP, script in D:\Alzheimer\syn18485175\Cibersortx\code.R
phenotype=read.table("D:/Alzheimer/Syn3388564/Phenotype.txt",header=T,row.names=1,sep="\t") #the data download from synapse database
trait=phenotype[,c(13,15,16,17)]
t=pheatmap(trait,clustering_method="ward.D2",show_rownames=F)
a=data.frame(cutree(t$tree_row,k=4))
anno=data.frame(Group=paste0("G",a[,1],sep=""))
rownames(anno)=rownames(a)
anno$ID=rownames(anno)
ratioPathologyDegree=merge(anno,ratio,by="ID")
Sex=factor(ratioPathologyDegree$msex,levels=c("Male","Female"))
Statues=factor(ratioPathologyDegree$Group.x,levels=c("G4","G2","G3","G1"))
t=ggplot(ratioPathologyDegree, aes(x=Statues, y=Mic*100, fill=Sex)) +
  geom_boxplot(position=position_dodge(0.8),outlier.shape = NA) +  
  geom_dotplot(binaxis='y', stackdir='centerwhole',  position=position_dodge(0.8),dotsize=0.6,binwidth=0.2)+
  scale_fill_manual(values = c("RoyalBlue","Violet"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  theme(axis.text.x = element_text(size=rel(1.0)),axis.text.y = element_text(size=rel(1.0)))
pdf("Graph/Part1/MicCellTypeFromROSMAP8ScadenBetADDegree.pdf",height=4,width=5)
print(t)
dev.off()

ratioPathologyDegree$Group=paste0(ratioPathologyDegree$Group.x,ratioPathologyDegree$msex,sep="")
wilcox.test(ratioPathologyDegree[ratioPathologyDegree$Group=="G1Female","Mic"],ratioPathologyDegree[ratioPathologyDegree$Group=="G1Male","Mic"]) #p-value = 0.004165
wilcox.test(ratioPathologyDegree[ratioPathologyDegree$Group=="G3Female","Mic"],ratioPathologyDegree[ratioPathologyDegree$Group=="G3Male","Mic"]) #p-value = 0.2813
wilcox.test(ratioPathologyDegree[ratioPathologyDegree$Group=="G2Female","Mic"],ratioPathologyDegree[ratioPathologyDegree$Group=="G2Male","Mic"]) #p-value = 0.8847
