# scRNA-seq analysis of in-vitro Leydig differentiation at 3 different time points by Qianyi on 10.18.2019 
# clustering for each of the 3 time points of in-vitro leydig differentiation separately 


library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
redblue100<-rgb(read.table(paste0(home,'data_DGE/redblue100.txt'),sep='\t',row.names=1,header=T))
library(RColorBrewer)
myBrewerPalette=c(brewer.pal(7,"Set2"),brewer.pal(12,"Paired")[c(10,12)])

infile=read.table("topcellS",stringsAsFactors=F)
nbatch=nrow(infile)
infile
dataset=infile[,3]
subject=unique(gsub("-B","",gsub("-A","",dataset)))
#140460  3300  LM4day
#140461  2400  LM7day
#140462  1300  LM14day



########## Analysis for individual dataset
setlist=list()
dgelist=list()
countcell=countgene=NULL
for(i in 1:length(dataset)){
### read in data
dgedata=read.table(paste0(infile[i,1],"/out_",infile[i,2],"cells_gene_exon_tagged_clean.dge.txt.gz"),header=T,row.names=1)
colnames(dgedata)=paste(infile[i,3],colnames(dgedata),sep="_")
### Filter for cells (>500 genes and <10% MT)
nGeneperCell <- colSums(dgedata>0)
dgedata.tmp=dgedata[,nGeneperCell>500]
mito.genes <- grep("^MT-|^mt-", rownames(dgedata.tmp), value = T) 
percent.mito <- colSums(dgedata.tmp[mito.genes, ])/colSums(dgedata.tmp)
dgedata.tmp=dgedata[,nGeneperCell>500][,percent.mito<0.1] 
print(c(ncol(dgedata),ncol(dgedata[,nGeneperCell>500]),ncol(dgedata.tmp)))
countcell=rbind(countcell,c(ncol(dgedata),ncol(dgedata[,nGeneperCell>500]),ncol(dgedata.tmp)))
### Filter for genes (Non-0 genes)
nCellperGene <- rowSums(dgedata.tmp>0)
dgedata2=dgedata.tmp[which(nCellperGene >= 1),]
print(c(nrow(dgedata),nrow(dgedata2)))
countgene=rbind(countgene,c(nrow(dgedata),nrow(dgedata2)))
print(summary(rowSums(dgedata2)))
setlist[[i]]=dgedata2
### Setup Seurat object
dge <- CreateSeuratObject(raw.data = dgedata2, project=infile[i,3], min.cells=1, min.genes=1)
mito.genes <- grep(pattern = "^MT-|^mt-", x = rownames(x = dge@data), value = TRUE)
percent.mito <- Matrix::colSums(dge@raw.data[mito.genes, ])/Matrix::colSums(dge@raw.data)
dge <- AddMetaData(object = dge, metadata = percent.mito, col.name = "percent.mito")
### Normalize data
dge <- NormalizeData(object = dge)
### Highly variable genes
pdf(paste0(dgefile,infile[i,3],"_HVG_x0.2_y0.2.pdf"),height=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
dge <- FindVariableGenes(object = dge, x.low.cutoff = 0.2, x.high.cutoff = 10, y.cutoff = 0.2, do.plot = TRUE, col.use=rgb(0,0,0,0.8))
dev.off()
print(length(dge@var.genes))
### Scale data 
dge <- ScaleData(object = dge)
### PCA
Sys.time()  
dge <- RunPCA(dge, pc.genes = dge@var.genes,pcs.compute = 40,  do.print = TRUE, pcs.print = 5, genes.print = 5)
dge <- ProjectPCA(object = dge, do.print = FALSE)
dgelist[[i]]=dge
save(dge,file=paste0(infile[i,3],".Robj"))
}

### Initial number of cells and genes
countcell
countgene
### Number of detected genes and filtered cells
for(i in 1:length(dataset)){
	print(dim(setlist[[i]]))
}
for(i in 1:length(dataset)){
  print(dim(dgelist[[i]]@data))
} ## Make sure the same 
### Number of HVG
for(i in 1:length(dataset)){
  print(length(dgelist[[i]]@var.genes))
}
### Per-cell attributes: Depth, Average number of genes, UMIs, %MT per cell
#pdf(file=paste0(dgefile,"NGeneUMImt.pdf"))
for(i in 1:length(dataset)){
dge=dgelist[[i]]
#VlnPlot(dge, features.plot=c("nGene", "nUMI", "percent.mito"), nCol = 3)
print(c(mean(dge@meta.data$nGene),mean(dge@meta.data$nUMI),mean(dge@meta.data$percent.mito)))
}
#dev.off()

###### Determine top PCs
numPCs= c(8,9,8)

pdf(paste(dgefile,"dge_PCA_Variablel_variation.pdf",sep=""),height=18,width=18)
par(mfrow=c(6,6),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
for(i in 1:length(dataset)){
dge=dgelist[[i]]
plot(dge@dr$pca@sdev[1:40],type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
legend("topright",legend=dataset[i],cex=1.5)
abline(v=numPCs[i]+0.5,lty=2)
}


###### Find clusters using top PCs
for(i in 1:length(dataset)){
dge=dgelist[[i]]
### tSNE
	dge <- RunTSNE(dge, dims.use = 1:numPCs[i], do.fast = T)
  dge <- RunUMAP(dge, dims.use = 1:numPCs[i])
### Louvain-Jaccard Clustering
dge <- FindClusters(dge, reduction.type = "pca", dims.use = 1:numPCs[i], resolution = seq(0.1,3,by=0.1), save.SNN = T)
dgelist[[i]]=dge
save(dge,file=paste0(infile[i,3],".Robj"))
}


for(i in 1:length(dataset)){
dge=dgelist[[i]]
print(c( length(unique(dge@meta.data$res.0.1)),length(unique(dge@meta.data$res.0.2)),length(unique(dge@meta.data$res.0.3)),length(unique(dge@meta.data$res.0.4)),length(unique(dge@meta.data$res.0.5)),length(unique(dge@meta.data$res.0.6)),length(unique(dge@meta.data$res.0.7)),length(unique(dge@meta.data$res.0.8)),length(unique(dge@meta.data$res.0.9)),length(unique(dge@meta.data$res.1)),length(unique(dge@meta.data$res.1.1)),length(unique(dge@meta.data$res.1.2)) ))
}



###### order clusters for each dataset
res=paste0("res.0.",c(5,5,7))

## order cell by cluster ID and randomly shuffle cells within each batch
levelss=list()
for(resi in 1:length(dataset)){
dge=dgelist[[resi]]
dge=SetAllIdent(dge,id=res[resi])
print(c(dataset[resi],length(unique(dge@ident))))
levels=levels(dge@ident)
ident=factor(dge@ident,levels=levels)

### randomly shuffling cells within each cluster
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

### for each cluster, calculate average normalized expression of each gene
tmpdge=data.frame(t(as.matrix(dge@data[,cells.use])))
# make sure same order for cells.ident and dge before combining
which(names(cells.ident)!=colnames(dge@data)[cells.use])  # nothing
mouseclustersall=data.frame(ident=cells.ident,t(as.matrix(dge@data[,cells.use])))
genecountsall=matrix(,dim(dge@data)[1],length(unique(mouseclustersall$ident)))
rownames(genecountsall)=rownames(dge@data)
colnames(genecountsall)=unique(mouseclustersall$ident)
for(i in unique(mouseclustersall$ident)){
    genecountsall[,i]=apply(mouseclustersall[mouseclustersall$ident==i,-1],2,function(x) ExpMean(as.numeric(x))) # log(mean(exp(x) - 1) + 1)
    print(i)
}

### Reordering cluster centroid using dissimilarity matrix
library(seriation)
n=ncluster=length(levels)
nbatch=1 # nbatch=length(dataset)
bb=1
tmp=genecountsall[,levels]
tmp=tmp
colnames(tmp)=gsub(".*_","",colnames(tmp))
da <- dist(t(as.matrix(tmp)), method = "euclidean")
pdf(file=paste0(dgename,"Centroid_norm_Seriation_",infile[resi,3],"_",res[resi],".pdf"))
 dissplot(da, method="OLO",options = list(main = paste("Dissimilarity with seriation OLO")))
 hmap(da) # default method="OLO"
dev.off()

### get order of seriation
 do=seriate(da,method="OLO")
levelss[[resi]]=get_order(seriate(da,method="OLO"))
# levelss=levels[get_order(seriate(da,method="OLO"))]
levelss[[resi]]
levels=levelss[[resi]]
if(resi>2){
  levels=rev(levelss[[resi]])
}

### Reordered clusters for all cells
cells.use=colnames(dge@data)
# random shuffling cells within ordered clusters
ident=factor(dge@ident,levels=levels-1)

cells=sort(ident)
cells.use=NULL
for(i in 1:length(levels)){
   set.seed(i)
   tmp=cells[which(cells == levels[i]-1)]
   if(length(tmp)>0){
      tmpname=sample(names(tmp),length(tmp),replace=FALSE)
      cells.use=c(cells.use,tmpname)
   }
}
cells.ident.ordered=factor(as.numeric(cells),ordered=TRUE)
names(cells.ident.ordered)=cells.use

### save ordered cluster ID in dge object
which(unique(cells.ident)!=get_order(do)) # integer(0)

ordered="ordered"

dge=AddMetaData(dge,cells.ident.ordered,ordered)
dge@meta.data[,ordered]=factor(dge@meta.data[,ordered])
dge=SetAllIdent(dge,ordered)
# save the dge file
dgelist[[resi]]=dge
#save(dge,file=paste0("data_DGE/",infile[resi,3],".Robj"))
}
# go above to re-plot PCA and tSNE

### save dge with ordered clusters
for(i in 1:length(dataset)){
dge=dgelist[[i]]
save(dge,file=paste0(infile[i,3],".Robj"))
}


## double-check if I need to reverse the cluster ID orders
### put SPG (Zbtb16) / Leydig (Star) as beginning, and Myoid (Myh11) as end
knownmarkers=c("Mki67","Nr2f2","Pdgfra","Tcf21","Cyp17a1","Star","Acta2","Myh11","Clu","Sox9","Zbtb16","Kit","Stra8","Acrv1","Prm1")
plotlist=list()
for(i in 1:length(dataset)){
gene=knownmarkers[which(knownmarkers %in% rownames(dgelist[[i]]@data))]
plotlist[[i]]=VlnPlot(dgelist[[i]],gene,nCol=4)
}
pdf("clusters_ordered0_knownmarkers_Violin.pdf",height=8,width=12)
print(plotlist)
dev.off()
pdf("clusters_ordered0_LM_knownmarkers_Violin.pdf",height=10,width=12)
print(plotlist)
dev.off()
pdf("clusters_ordered0_knownmarkers_Feature.pdf",height=6,width=8)
for(i in 1:length(dataset)){
gene=knownmarkers[which(knownmarkers %in% rownames(dgelist[[i]]@data))]
  FeaturePlot(dgelist[[i]],gene,cols.use=c("gray80","red"))
}
dev.off()
#FeaturePlot(dgelist[[i]],"PRM1")
#TSNEPlot(dgelist[[i]],do.return=T)
# download the tSNE plot below to check the cluster IDs

## double-check if I need to flip PCs, or flip tSNE
### PCA and tSNE for ordered clusters of each individual replicate
plotlistt=plotlist2=plotlist3=plotlist4=plotlist5=list()
for(i in 1:length(dataset)){
dge=dgelist[[i]]
plotlist2[[i]]=PCAPlot(dge,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlist3[[i]]=PCAPlot(dge,1,3,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlist4[[i]]=PCAPlot(dge,1,4,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlist5[[i]]=PCAPlot(dge,1,5,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlistt[[i]]=TSNEPlot(dge,pt.size=1,do.return=TRUE,do.label=TRUE,label.size=4)
}
pdf("clusters_ordered0.pdf",height=8,width=18)
MultiPlotList(plotlist2,cols = 5)
MultiPlotList(plotlist3,cols = 5)
MultiPlotList(plotlist4,cols = 5)
MultiPlotList(plotlist5,cols = 5)
MultiPlotList(plotlistt,cols = 5)
dev.off()

### flip PC1
for(i in c(2)){
dge=dgelist[[i]]
dge@dr$pca@cell.embeddings[,1]=-dge@dr$pca@cell.embeddings[,1]
dgelist[[i]]=dge
}
### flip PC2
for(i in c(1)){
dge=dgelist[[i]]
dge@dr$pca@cell.embeddings[,2]=-dge@dr$pca@cell.embeddings[,2]
dgelist[[i]]=dge
}
### flip PC3
for(i in c(1)){
dge=dgelist[[i]]
dge@dr$pca@cell.embeddings[,3]=-dge@dr$pca@cell.embeddings[,3]
dgelist[[i]]=dge
}
### flip tSNE2
for(i in c(1)){
dge=dgelist[[i]]
dge@dr$tsne@cell.embeddings[,2]=-dge@dr$tsne@cell.embeddings[,2]
dgelist[[i]]=dge
}
### switch tSNE1 and tSNE2
for(i in c(2)){
dge=dgelist[[i]]
tmp=dge@dr$tsne@cell.embeddings[,2]
dge@dr$tsne@cell.embeddings[,2]=dge@dr$tsne@cell.embeddings[,1]
dge@dr$tsne@cell.embeddings[,1]=tmp
dgelist[[i]]=dge
}


## after reverse cluster order and flip cooridnates, repeat the above plots
# change ordered0 in the file name to order1
pdf("clusters_order1_PRM1_Violin.pdf",height=4,width=8)
dev.off()
#FeaturePlot(dgelist[[i]],"PRM1")
#TSNEPlot(dgelist[[i]],do.return=T)
# download the tSNE plot below to check the cluster IDs

## double-check if I need to flip PCs, or flip tSNE
### PCA and tSNE for ordered clusters of each individual replicate
plotlistt=plotlist2=plotlist3=plotlist4=plotlist5=list()
for(i in 1:length(dataset)){
dge=dgelist[[i]]
plotlist2[[i]]=PCAPlot(dge,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlist3[[i]]=PCAPlot(dge,1,3,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlist4[[i]]=PCAPlot(dge,1,4,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlist5[[i]]=PCAPlot(dge,1,5,pt.size=1,do.return=TRUE,do.label=TRUE)
plotlistt[[i]]=TSNEPlot(dge,pt.size=1,do.return=TRUE,do.label=TRUE,label.size=4)
}
pdf("clusters_order1.pdf",height=8,width=18)
multiplot(plotlist2,cols = 5)
multiplot(plotlist3,cols = 5)
multiplot(plotlist4,cols = 5)
multiplot(plotlist5,cols = 5)
multiplot(plotlistt,cols = 5)
dev.off()

### Visualize PC1-3 and tSNE for each time point in one plot
plotlist2=plotlist3=plotlistt=plotlistu=list()
for(i in 1:length(dataset)){
  rep=dataset[i]
dge=dgelist[[rep]]
plotlist2[[i]]=PCAPlot(dge,pt.size=1,cols.use=myBrewerPalette,do.return=TRUE,do.label=TRUE)
plotlist3[[i]]=PCAPlot(dge,1,3,pt.size=1,cols.use=myBrewerPalette,do.return=TRUE,do.label=TRUE)
plotlistt[[i]]=TSNEPlot(dge,pt.size=1,colors.use=myBrewerPalette,do.return=TRUE,do.label=TRUE,label.size=4)
plotlistu[[i]]=DimPlot(dge,pt.size=1,reduction.use="umap",cols.use=myBrewerPalette,do.return=TRUE,do.label=TRUE)
}
pdf(paste0("subject_rep_clusters_order1.pdf"),height=2.7*2,width=3.5*4)
multiplot(c(plotlist2,plotlist3,plotlistt,plotlistu),cols = 4)
dev.off()


### saved Seurat R object for each time point of in-vitro leydig differentiation
for(i in 1:length(dataset)){
dge=dgelist[[i]]
save(dge,file=paste0(infile[i,3],".Robj"))
}

