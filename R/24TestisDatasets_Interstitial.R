### R script for focused analysis of all Interstitial Cells in merged 24 ST batches in Jul 2017 by Qianyi


home="/scratch/junzli_flux/qzm/Dropseq_analysis/"
redblue100<-rgb(read.table('/scratch/junzli_flux/qzm/Dropseq_analysis/data_DGE/redblue100.txt',sep='\t',row.names=1,header=T))

library(Seurat)
library(dplyr)
library(Matrix)

load(file =paste0(home,"data_DGE/MouseAdultST24genesUMI20cell15.Robj"))
dge=SetAllIdent(dge,id="clusters31_GeneFilter1")


sets=list(c(6:10))
setsname=c("6-10")
setslabel=c("Interstitial")


######### PCA for interstitial cells only
j=1
print(j)
### Subset of cells in neighboring clusters
dgetmp=SubsetData(dge,ident.use=sets[[j]])
### genes expressed in the subset of cells of neighboring clusters
nCellperGene <- rowSums(as.matrix(dgetmp@data)>0)
genes.use=names(nCellperGene[which(nCellperGene!=0)])
dgetmp@data=dgetmp@data[genes.use,]
print(table(dgetmp@ident))
print(dim(dgetmp@data))
print(Sys.time())
dgetmp <- PCA(dgetmp, pc.genes = rownames(dgetmp@data), do.print = TRUE, pcs.print = 5, genes.print = 5)
print(Sys.time())
save(dgetmp,file =paste0(home,"data_DGE/MouseAdultST24genesUMI20cell15_31clusters_cluster",setsname[j],".Robj"))
plotall=list()
plotall[[1]]=PCAPlot(dgetmp,1,2,do.return = TRUE,pt.size = 1,do.label=T)
plotall[[2]]=PCAPlot(dgetmp,1,3,do.return = TRUE,pt.size = 1,do.label=T)
pdf(paste(dgefile,"ZoomedInReDo_PCA_NeighborClusters_",setsname[j],".pdf",sep=""),height=10,width=6)
MultiPlotList(plotall,cols = 1)
dev.off()

dge=dgetmp
#check components representing greatest variation
plot(dge@pca.obj[[1]]$sdev[1:120],type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) 
numPCs=14

pdf(paste(dgefile,"dge_PCA_all2.pdf",sep=""),height=7.5,width=16)
plot1=PCAPlot(dge,1,4,do.return = TRUE,pt.size = dge@data.info$nUMIperCell2)
plot2=PCAPlot(dge,1,5,do.return = TRUE,pt.size = dge@data.info$nUMIperCell2)
MultiPlotList(list(plot1,plot2),cols = 2)
dev.off()


dge=RunTSNE(dge,dims.use = 1:numPCs[j],do.fast=T)    # max_iter=2000
dge <- FindClusters(dge, pc.use = 1:numPCs[j], resolution = seq(0.1,3,0.1), print.output = 0, save.SNN = T)

res="res.0.2";i=1
dge <- SetAllIdent(dge, id = res[i])
dge <- BuildClusterTree(dge, do.reorder = T, reorder.numeric = T,pcs.use=1:numPCs[i])

library(seriation)
levels=levels(dge@ident)
n=ncluster=length(levels)
nbatch=1 # nbatch=length(dataset)
bb=1
### Caculate Dissimilarity matrix distance and Do Seriation
tmp=genecountsall[,levels][,(sum(ncluster[-(bb:nbatch)])+1):(sum(ncluster[1:bb]))]
tmp=tmp
colnames(tmp)=gsub(".*_","",colnames(tmp))
da <- dist(t(as.matrix(tmp)), method = "euclidean")
# note: dist calculate distance between each row
length(da) # 91
da
# order by seriation
print(get_order(seriate(da,method="OLO")))
 # plot with seriation
pdf(file=paste(dgename,"Centroid_norm_Seriation_Dissimilarity.pdf"))
 hmap(da) # default method="OLO"
dev.off()

save(dge, file =paste0(home,"data_DGE/MouseAdultST24genesUMI20cell15_31clusters_cluster",setsname[j],".Robj"))



