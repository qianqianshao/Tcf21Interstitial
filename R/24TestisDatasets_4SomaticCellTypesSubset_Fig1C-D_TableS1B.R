### R script for zoomed-in view of 4 intertitial cell types in 24 testis datasets by Qianyi
### Related to Figure 1C-D and Table S1B

# Used somatic cell type assignments from analysis of merged 25 testis datasets in Green et al, Cell 2018 paper. Please refer to:  
# https://github.com/qianqianshao/Drop-seq_ST/blob/master/R/FocusedClustering_SomaticCells/Somatic_7SomaticCellTypes.R
# Removed INT6 dataset from merged 25 testis datasets to avoid batch effect of Tcf21+ dataset -> 24 testis datasets 
# removed InnateLymphoid which are very few in the 24 testis datasets


home="/scratch/junzli_root/junzli/qzm/Dropseq_analysis/"
setwd(home)

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
redblue100<-rgb(read.table('data_DGE/redblue100.txt',sep='\t',row.names=1,header=T))
source(paste0(home,"data_DGE/Rcode_multiplot.R"))


dgefile=dgename=""

### 4 somatic cell types for 24 ST datasets - kept cell type annotation from 25 ST datasets
sets=c("Endothelial","Myoid","Unknown","Leydig")
setsname=setslabel=c("NoINT6_Endo-Myo-Unk-Ley")


### 4 colors from ggplot
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
myBrewerPalette=gg_color_hue(4)


### load all soma of 25 ST datasets 
load(file ="data_DGE/MouseAdultST25Somatic.Robj")
dgeall=dge


### Extract Endo-Myo-Unk-Ley cells from 25 ST datasets, then remove INT6 dataset
sets=c("Endothelial","Myoid","Unknown","Leydig")
setsname=setslabel=c("NoINT6_Endo-Myo-Unk-Ley")
dge=dgeall
table(dge@ident)
cells.use1=names(dge@ident[which(dge@ident %in% sets)])
table(dge@ident[cells.use1])
   Endothelial InnateLymphoid         Leydig     Macrophage          Myoid 
           179              0            314              0             49 
       Sertoli        Unknown 
             0           2205 

cells.use=cells.use1[!grepl("INT6",cells.use1)]
table(gsub("_.*","",cells.use))
#INT1 INT2 INT3 INT4 INT5 SER2 SER3 SER4 SER5 SER6 SER7 SER8 SPG1 SPG3  ST1  ST2 
# 100  131   66  483  145    1    2   25   30    4   15    3    1   88   43   52 
# ST3  ST4  ST5  ST6  ST7  ST8 
#  58   77   56   40   19   20 
table(dge@ident[cells.use])
   Endothelial InnateLymphoid         Leydig     Macrophage          Myoid 
           151              0            292              0             49 
       Sertoli        Unknown 
             0            967


dge=dgeall
dgedata=dge@raw.data[,cells.use]
length(cells.use) # 1459 
dim(dgedata)         
### Filter for genes (Non-0 genes)
nCellperGene <- rowSums(dgedata>0)
dgedata2=dgedata[which(nCellperGene >= 1),]
dim(dgedata2)   # [1] 20568  1459
### Setup Seurat object
dge <- CreateSeuratObject(raw.data = dgedata2, project=setsname, min.cells=1, min.genes=1)
mito.genes <- grep(pattern = "^MT-|^mt-", x = rownames(x = dge@data), value = TRUE)
percent.mito <- Matrix::colSums(dge@raw.data[mito.genes, ])/Matrix::colSums(dge@raw.data)
dge <- AddMetaData(object = dge, metadata = percent.mito, col.name = "percent.mito")
### Normalize data
dge <- NormalizeData(object = dge)
### Highly variable genes
pdf(paste0(dgename,"HVG_x0.1_y0.1.pdf"),height=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
dge <- FindVariableGenes(object = dge, x.low.cutoff = 0.1, x.high.cutoff = 10, y.cutoff = 0.1, do.plot = TRUE, col.use=rgb(0,0,0,0.8))
dev.off()
print(length(dge@var.genes)) #2547 # 2443
### Scale data 
dge <- ScaleData(object = dge)
### PCA
Sys.time()  
dge <- RunPCA(dge, pc.genes = dge@var.genes,pcs.compute = 40,  do.print = TRUE, pcs.print = 5, genes.print = 5)
dge <- ProjectPCA(object = dge, do.print = FALSE)

### Scree Plot for single PCA
numPCs=10    # ST25-NoINT6_Endo-Myo-Unk-Ley
i=1

pdf(paste(dgefile,"dge_PCA_Variablel_variation.pdf",sep=""),height=3.5,width=7)
par(mfrow=c(1,2),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(dge@dr$pca@sdev[1:40],type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
abline(v=numPCs[i]+0.5,lty=2)
eigenvalue=dge@dr$pca@sdev[numPCs[i]]
text(numPCs[i],eigenvalue+5,col="red",paste(numPCs[i],"PCs"))
eigenvalue=dge@dr$pca@sdev[numPCs[i]]
print(eigenvalue)
plot(density(dge@dr$pca@sdev),col="red",lwd=2,xlab="Eigenvalue",main="",cex.lab=1.5)
polygon(density(dge@dr$pca@sdev),col="black")
lines(density(dge@dr$pca@sdev),col="red",lwd=2,xlab="Eigenvalue")
abline(v=eigenvalue,col="red",lwd=2,lty=2)
text(eigenvalue+0.2,1,col="red",paste(numPCs[i],"PCs"))
dev.off()


  dge <- RunTSNE(dge, dims.use = 1:numPCs[i], do.fast = T)
  dge <- RunUMAP(dge, dims.use = 1:numPCs[i])

ident=factor(dgeall@ident[cells.use],levels=sets)
table(ident)
#Macrophage Endothelial       Myoid     Unknown      Leydig     Sertoli 
#         43         151          49         967         292        2120
#Endothelial       Myoid     Unknown      Leydig 
#        151          49         967         292 
dge=AddMetaData(dge,ident,"celltype")
dge <- SetAllIdent(dge, id = "celltype")
dge@ident=factor(dge@ident,levels=sets)
table(dge@ident)
Endothelial       Myoid     Unknown      Leydig 
        151          49         967         292 

save(dge,file=paste0("data_DGE/MouseAdultST25",setsname,".Robj"))



### visualization for NoINT6_Endo-Ley-Unknown-Myoid
pdf(paste(dgefile,"PCA_tSNE_UMAP.pdf",sep=""),height=5.8,width=8.5)
plot2=PCAPlot(dge,do.return=TRUE,pt.size=0.8,do.label=TRUE)
plot3=PCAPlot(dge,1,3,do.return=TRUE,pt.size=0.8,do.label=TRUE)
plotu=DimPlot(dge,reduction.use="umap",pt.size=0.8,do.return = TRUE,do.label=TRUE,)
plott=TSNEPlot(dge,do.return=TRUE,pt.size=0.8,do.label=TRUE,label.size=5)
multiplot(list(plot2,plot3,plott,plotu),cols=2)
plot2=PCAPlot(dge,do.return=TRUE,pt.size=0.8)
plot3=PCAPlot(dge,1,3,do.return=TRUE,pt.size=0.8)
plotu=DimPlot(dge,reduction.use="umap",pt.size=0.8,do.return = TRUE)
plott=TSNEPlot(dge,do.return=TRUE,pt.size=0.8)
multiplot(list(plot2,plot3,plott,plotu),cols=2)
dev.off()


markersall=FindAllMarkers(dge,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
table(markersall$cluster)
#Endothelial       Myoid     Unknown      Leydig 
#        218         117         105         206
write.table(markersall,paste0("data_DGE/MouseAdultST25_MarkersAll_",setsname,"_pct0.2_diffpct0.2_thresh2fold_1.27.2020.txt"),col.names=T,row.names=T,quote=F,sep="\t")
# saved as Table S1B

######## Heatmap for markers across NoINT6_Endo-Myo-Unk-Ley
###### 1. Heatmap for markers in all cells in NoINT6_Endo-Myo-Unk-Ley
redblue100<-rgb(read.table(paste0(home,"data_DGE/redblue100.txt"),sep='\t',row.names=1,header=T))

dge # 20568 genes across 1459 samples

dge=SetAllIdent(dge,id="celltype")
dge@ident=factor(dge@ident,levels=sets)
celltype=ident=dge@ident

levels=levels(celltype)
cells=sort(ident)
cells.use=NULL
for(i in 1:length(levels)){
   set.seed(i)
   tmp=cells[which(cells == levels[i])]
   if(length(tmp)>0){
      tmpname=sample(names(tmp),length(tmp),replace=FALSE)
      cells.use=c(cells.use,tmpname)
   }
}
cells.ident=as.factor(cells)
names(cells.ident)=cells.use
levels(cells.ident)=levels

minmax=function(data,min,max) {
  data2=data
  data2[data2>max]=max
  data2[data2<min]=min
  return(data2)
}
disp.min=-2.5;disp.max=2.5;draw.line=TRUE;
            data.use2=dge@scale.data[markersall$gene,cells.use]
            data.use2=minmax(data.use2,min=disp.min,max=disp.max)

            lab2=rep("",length(cells.use))
            lab2[round(cumsum(table(cells.ident)[levels(cells.ident)])-table(cells.ident)[levels(cells.ident)]/2)+15]=levels(cells.ident)

            row.lab2=gsub(".*_","",lab2)

            orig.ident=factor(gsub("_.*","",cells.ident),levels=unique(gsub("_.*","",cells.ident)))
            col.lab2=rep("",length(cells.use))
            col.lab2[round(cumsum(table(orig.ident)[levels(orig.ident)])-table(orig.ident)[levels(orig.ident)]/2)+5]=levels(orig.ident)

            colsep.use2=cumsum(table(cells.ident)[levels(cells.ident)])
            colsep.use2=cumsum(table(orig.ident)[levels(orig.ident)]) # draw a line between datasets

sidecol2=do.call(rbind,strsplit(as.character(cells.ident),"_"))
sidecol2=cbind(sidecol2,sidecol2)
for(rep in 1:length(unique(sidecol2[,1]))){
a=unique(sidecol2[,1])[rep]
sidecol2[which(sidecol2[,1]==a),2]<-rep
}

rlab2=rbind(rep("white",length(levels))[as.numeric(sidecol2[,2])],myBrewerPalette[as.numeric(sidecol2[,2])])
clab2=cbind(rlab2[2,],rlab2[1,])
rownames(rlab2)=c("","Cell Type")
colnames(clab2)=c("Cell Type","")

col.use2=col.use

library(gplots)
library(devtools)

jpeg(file=paste(dgename,"allcells_markersall_clab.jpeg",sep=""),height=3000,width=2800,res=300)
par(mar=c(10,4,1,2),mgp=c(2.5, 1, 0))
heatmap.3(data.use2,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use2,colsep = colsep.use2,sepcolor="black",sepwidth=c(0.01,0.01),ColSideColors=clab2,labCol=col.lab2,cexCol=1.5,cexRow=0.5,ColSideColorsSize = 3,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F,scale="none",margins=c(7,5))                    # symm=F,symkey=F,symbreaks=F,
dev.off()

# saved as Figure 1D

