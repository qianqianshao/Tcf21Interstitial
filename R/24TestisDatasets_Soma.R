### R script for visualization of somatic cell types in 24 testis datasets using cell type classification of merged 25 testis datasets by Qianyi
### Related to Figure 1A and 1B

# Used somatic cell type assignments from analysis of merged 25 testis datasets in Green et al, Cell 2018 paper. Please refer to:  
# https://github.com/qianqianshao/Drop-seq_ST/blob/master/R/FocusedClustering_SomaticCells/Somatic_7SomaticCellTypes.R
# Remove INT6 dataset from merged 25 testis datasets -> 24 testis datasets 
# remove InnateLymphoid which are very few in the 24 testis datasets

home="/scratch/junzli_root/junzli/qzm/Dropseq_analysis/"
setwd(home)

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
redblue100<-rgb(read.table('data_DGE/redblue100.txt',sep='\t',row.names=1,header=T))
source(paste0(home,"data_DGE/Rcode_multiplot.R"))


dgefile=dgename="figJan2020_MouseAdultST25NoINT6_Endo-Myo-Unk-Ley/"

### 6 somatic cell types for 24 ST datasets - assigned from 25 ST datasets
sets=c("Macrophage","Endothelial","Myoid","Unknown","Leydig","Sertoli")
setsname=setslabel=c("ST24NoINT6_NoInnateLymph")

load(paste0("data_DGE/MouseAdultST25",setsname,".Robj"))
dge #  22565 genes across 3622 samples


### 4 colors from ggplot
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
myBrewerPalette=gg_color_hue(4)


######### PCA for soma of 24 ST datasets (removing INT6 from ST25)
dge=dgeall
cells.use=names(dge@ident[!grepl("INT6",names(dge@ident))])
table(gsub("_.*","",cells.use))
table(dge@ident[cells.use])
   Endothelial InnateLymphoid         Leydig     Macrophage          Myoid 
           151              6            292             43             49 
       Sertoli        Unknown 
          2120            967
# also remove innate lymphoid as there are very few
cells.use=cells.use[which(dge@ident[cells.use]!="InnateLymphoid")]
sets=c("Macrophage","Endothelial","Myoid","Unknown","Leydig","Sertoli")
setsname=setslabel=c("ST24NoINT6_NoInnateLymph")

dge=dgeall
dgedata=dge@raw.data[,cells.use]
length(cells.use) # 3622 
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
numPCs=10     # ST25-NoINT6_NoInnateLymphoid
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
dge=AddMetaData(dge,ident,"celltype")
dge <- SetAllIdent(dge, id = "celltype")
dge@ident=factor(dge@ident,levels=sets)
table(dge@ident)

save(dge,file=paste0("data_DGE/MouseAdultST25",setsname,".Robj"))


### Visualization for 24ST-NoINT6_NoInnateLymph
library(RColorBrewer)
display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE, colorblindFriendly=TRUE)
myBrewerPalette=c(brewer.pal(12,"Paired")[3],gg_color_hue(4),brewer.pal(12,"Paired")[1])
dge@dr$pca@cell.embeddings=as.matrix(dgeall@pca.rot[cells.use,])
pdf(paste(dgefile,setslabel,"PCA_tSNE_UMAP.pdf",sep=""),height=6.4,width=14)
plot2=PCAPlot(dge,do.return=TRUE,cols.use=myBrewerPalette,pt.size=0.8,do.label=TRUE)
plot3=PCAPlot(dge,1,3,do.return=TRUE,cols.use=myBrewerPalette,pt.size=0.8,do.label=TRUE)
plot4=PCAPlot(dge,1,4,do.return=TRUE,cols.use=myBrewerPalette,pt.size=0.8,do.label=TRUE)
plot5=PCAPlot(dge,1,5,do.return=TRUE,cols.use=myBrewerPalette,pt.size=0.8,do.label=TRUE)
plot6=PCAPlot(dge,1,6,do.return=TRUE,cols.use=myBrewerPalette,pt.size=0.8,do.label=TRUE)
plot7=PCAPlot(dge,1,7,do.return=TRUE,cols.use=myBrewerPalette,pt.size=0.8,do.label=TRUE)
plotu=DimPlot(dge,reduction.use="umap",cols.use=myBrewerPalette,pt.size=0.8,do.return = TRUE,do.label=TRUE,)
plott=TSNEPlot(dge,do.return=TRUE,colors.use=myBrewerPalette,pt.size=0.8,do.label=TRUE,label.size=5)
multiplot(list(plot2,plot3,plot4,plot5,plott,plotu),cols=3)
plot2=PCAPlot(dge,do.return=TRUE,cols.use=myBrewerPalette,pt.size=0.8)
plot3=PCAPlot(dge,1,3,do.return=TRUE,cols.use=myBrewerPalette,pt.size=0.8)
plot4=PCAPlot(dge,1,4,do.return=TRUE,cols.use=myBrewerPalette,pt.size=0.8)
plot5=PCAPlot(dge,1,5,do.return=TRUE,cols.use=myBrewerPalette,pt.size=0.8)
plot6=PCAPlot(dge,1,6,do.return=TRUE,cols.use=myBrewerPalette,pt.size=0.8)
plot7=PCAPlot(dge,1,7,do.return=TRUE,cols.use=myBrewerPalette,pt.size=0.8)
plotu=DimPlot(dge,reduction.use="umap",cols.use=myBrewerPalette,pt.size=0.8,do.return = TRUE)
plott=TSNEPlot(dge,do.return=TRUE,colors.use=myBrewerPalette,pt.size=0.8)
multiplot(list(plot2,plot3,plot4,plot5,plott,plotu),cols=3)
dev.off()
# saved as Figure 1A

### Seriation for 24ST-NoINT6_NoInnateLymph
genecountsall=log(AverageExpression(dge)+1)
levels=sets
library(seriation)
n=ncluster=length(levels)
nbatch=1 # nbatch=length(dataset)
bb=1
tmp=genecountsall[,levels]
tmp=tmp
colnames(tmp)=gsub(".*_","",colnames(tmp))
da <- dist(t(as.matrix(tmp)), method = "euclidean")
# note: dist calculate distance between each row
length(da) # 91
da
do=seriate(da,method="OLO")
print(get_order(do))

###### plot with seriation
pdf(file=paste0(dgefile,setslabel,"Centroid_norm_Seriation.pdf"))
 hmap(da,margin=c(6,6)) # default method="OLO"
dev.off()
# saved as Figure 1B



markersall=FindAllMarkers(dge,only.pos=TRUE,test.use="bimod",logfc.threshold = log(2),min.diff.pct=0.2,do.print = TRUE)
table(markersall$cluster)
#Macrophage Endothelial       Myoid     Unknown      Leydig     Sertoli 
#        224         284         169         274         242         370 
write.table(markersall,paste0("data_DGE/MouseAdultST25_MarkersAll_",setsname,"_pct0.2_diffpct0.2_thresh2fold_1.27.2020.txt"),col.names=T,row.names=T,quote=F,sep="\t")

