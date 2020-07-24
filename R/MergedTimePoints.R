# scRNA-seq analysis of in-vitro Leydig differentiation at 3 different time points by Qianyi on 10.18.2019 
# merge the INT4 soma with 3 time points of in-vitro leydig differentiation data and do clustering


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
#140460  3300  LM4day
#140461  2400  LM7day
#140462  1300  LM14day
dataset=infile[,3]
subject="LM"
exp="OldSca1LM"


### merge somatic cells of Old Adult Sca1 (INT4) together with in-vitro LM datasets
Sca1record=read.table("Sca1.record",stringsAsFactors=F)
indiv=1
Sca1=gsub("LM","",exp[indiv])
prefix=Sca1record[which(Sca1record[,1]==Sca1),3]
id=Sca1record[which(Sca1record[,1]==Sca1),2]

set=setlist=list()
###### Load raw gene expression matrix for old Sca1 (INT4) 
### Sca1 dataset with all cells and all genes
Sca1data=read.table(paste0("../data_DGE/",prefix,"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"), header=T,row.names=1)
print(dim(Sca1data))
# [1] 22758   653 
colnames(Sca1data)=paste(id,colnames(Sca1data),sep="_")
###### Load Adult Mouse Testis Merged Somatic Data - To extract somatic cell types from old Sca1 (INT4)
### merged somatic cell types with filtered genes
load(file =paste0("../data_DGE/MouseAdultST25Somatic.Robj"))
dgeSomatic=dge
dim(dgeSomatic@data) # 22734 genes across 5081 samples
### Extract somatic cell type from old Adult Sca1 datasets (INT4)
table(gsub("_.*","",names(dgeSomatic@ident)))
#INT1 INT2 INT3 INT4 INT5 INT6 SER1 SER2 SER3 SER4 SER5 SER6 SER7 SER8 SPG1 SPG2 
# 106  134  102  495  164 1453   24   33  118  218  518  531  198  411    3    1 
#SPG3  ST1  ST2  ST3  ST4  ST5  ST6  ST7  ST8 
#  91   48   69   71   91   73   55   32   42 
table(dgeSomatic@ident)
#   Endothelial InnateLymphoid         Leydig     Macrophage          Myoid 
#           179             64            314            139             49 
#       Sertoli        Unknown 
#          2131           2205
cells.Sca1=grep(id,names(dgeSomatic@ident),value=T)
print(length(cells.Sca1)) # 495 # 1453
adultSca1 <- Sca1data[,cells.Sca1]
nCellperGene<-rowSums(adultSca1>0) 
adultSca1 <- adultSca1[which(nCellperGene>0),]
print(dim(adultSca1)) 
# [1] 17997   495
# [1] 25187  1453
print(c(mean(dgeSomatic@data.info[cells.Sca1,]$nGene),mean(dgeSomatic@data.info[cells.Sca1,]$nUMI),mean(dgeSomatic@data.info[cells.Sca1,]$percent.mito)))
setlist[[1]]=adultSca1
table(dgeSomatic@ident[cells.Sca1])
#   Endothelial InnateLymphoid         Leydig     Macrophage          Myoid 
#            29              5             13              1              0 
#       Sertoli        Unknown 
#             6            441 
#   Endothelial InnateLymphoid         Leydig     Macrophage          Myoid 
#            28             58             22             96              0 
#       Sertoli        Unknown 
#            11           1238 

###### load Seurat object of merged all time points for in-vitro leydig differentiation
load(file="LM.Robj")
table(dgeall@ident)
#LM14day  LM4day  LM7day 
#    973    2980    2171 
setlist[[2]]=dgeall@raw.data
for(i in 1:2){
  set[[i]]=data.frame(GENE=rownames(setlist[[i]]),setlist[[i]])
}
dgedataall=Reduce(function(x,y) merge(x,y,all=TRUE), set)
dgedataall[is.na(dgedataall)] <- 0
row.names(dgedataall)=dgedataall[,1]
dgedataall=dgedataall[,-1]
dim(dgedataall)             
names(dgedataall)=gsub("^X","",names(dgedataall))
names(dgedataall)=gsub("\\.","-",names(dgedataall))
nCellperGene <- rowSums(dgedataall>0)
length(which(nCellperGene==0))
nCellperGene[which(nCellperGene==0)]

dge <- CreateSeuratObject(raw.data = dgedataall, project=exp[indiv], min.cells=1, min.genes=1)
mito.genes <- grep(pattern = "^MT-|^mt-", x = rownames(x = dge@data), value = TRUE)
percent.mito <- Matrix::colSums(dge@raw.data[mito.genes, ])/Matrix::colSums(dge@raw.data)
dge <- AddMetaData(object = dge, metadata = percent.mito, col.name = "percent.mito")
### Normalize data
dge <- NormalizeData(object = dge)
### Highly variable genes
pdf(paste0(dgefile,exp[indiv],"_HVG_x0.15_y0.2.pdf"),height=6)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
dge <- FindVariableGenes(object = dge, x.low.cutoff = 0.15, x.high.cutoff = 10, y.cutoff = 0.2, do.plot = TRUE, col.use=rgb(0,0,0,0.8))
dev.off()
print(length(dge@var.genes))
### Scale data 
dge <- ScaleData(object = dge)
### PCA
Sys.time()  
dge <- RunPCA(dge, pc.genes = dge@var.genes,pcs.compute = 40,  do.print = TRUE, pcs.print = 5, genes.print = 5)
dge <- ProjectPCA(object = dge, do.print = FALSE)
dgeall=dge
dgealllist[[indiv]]=dgeall
save(dgeall,file=paste0(exp[indiv],".Robj"))
table(gsub("_.*","",names(dge@ident)))
#   INT4 LM14day  LM4day  LM7day 
#    495     973    2980    2171
  save(dgeall,file=paste0(exp[i],".Robj"))

numPCs=6;i=1
pdf(paste(dgefile,"dgeall_PCA_Variablel_variation.pdf",sep=""),height=6,width=12)
par(mfrow=c(2,4),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(dgeall@dr$pca@sdev[1:40],type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
legend("topright",legend=exp[i],cex=1.5)
abline(v=numPCs[i]+0.5,lty=2)
dev.off()

### tSNE
  dgeall <- RunTSNE(dgeall, dims.use = 1:numPCs[i], do.fast = T)
  dgeall <- RunUMAP(dgeall, dims.use = 1:numPCs[i])
### Louvain-Jaccard Clustering
dgeall <- FindClusters(dgeall, reduction.type = "pca", dims.use = 1:numPCs[i], resolution = seq(0.1,3,by=0.1), save.SNN = T)
print(c( length(unique(dgeall@meta.data$res.0.1)),length(unique(dgeall@meta.data$res.0.2)),length(unique(dgeall@meta.data$res.0.3)),length(unique(dgeall@meta.data$res.0.4)),length(unique(dgeall@meta.data$res.0.5)),length(unique(dgeall@meta.data$res.0.6)),length(unique(dgeall@meta.data$res.0.7)),length(unique(dgeall@meta.data$res.0.8)),length(unique(dgeall@meta.data$res.0.9)),length(unique(dgeall@meta.data$res.1)),length(unique(dgeall@meta.data$res.1.1)),length(unique(dgeall@meta.data$res.1.2)) ))

###### order clusters for each exp
res=paste0("res.0.",7)
length(res)

## order cell by cluster ID and randomly shuffle cells within each batch
resi=1
dgeall=SetAllIdent(dgeall,id=res[resi])
print(c(exp[resi],length(unique(dgeall@ident))))
levels=levels(dgeall@ident)
ident=factor(dgeall@ident,levels=levels)

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
tmpdgeall=data.frame(t(as.matrix(dgeall@data[,cells.use])))
# make sure same order for cells.ident and dgeall before combining
which(names(cells.ident)!=colnames(dgeall@data)[cells.use])  # nothing
mouseclustersall=data.frame(ident=cells.ident,t(as.matrix(dgeall@data[,cells.use])))
genecountsall=matrix(,dim(dgeall@data)[1],length(unique(mouseclustersall$ident)))
rownames(genecountsall)=rownames(dgeall@data)
colnames(genecountsall)=unique(mouseclustersall$ident)
for(i in unique(mouseclustersall$ident)){
    genecountsall[,i]=apply(mouseclustersall[mouseclustersall$ident==i,-1],2,function(x) ExpMean(as.numeric(x))) # log(mean(exp(x) - 1) + 1)
    print(i)
}

### Reordering cluster centroid using dissimilarity matrix
library(seriation)
n=ncluster=length(levels)
nbatch=1 # nbatch=length(exp)
bb=1
tmp=genecountsall[,levels]
tmp=tmp
colnames(tmp)=gsub(".*_","",colnames(tmp))
da <- dist(t(as.matrix(tmp)), method = "euclidean")
pdf(file=paste0(dgefile,"Centroid_norm_Seriation_",exp[resi],"_",res[resi],".pdf"))
 dissplot(da, method="OLO",options = list(main = paste("Dissimilarity with seriation OLO")))
 hmap(da) # default method="OLO"
dev.off()

### get order of seriation
 do=seriate(da,method="OLO")
levels=get_order(seriate(da,method="OLO"))

### Reordered clusters for all cells
cells.use=colnames(dgeall@data)
# random shuffling cells within ordered clusters
ident=factor(dgeall@ident,levels=levels-1)

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

### save ordered cluster ID in dgeall object
which(unique(cells.ident)!=get_order(do)) # integer(0)

ordered=paste0(res[resi],"order")

dgeall=AddMetaData(dgeall,cells.ident.ordered,ordered)
dgeall@meta.data[,ordered]=factor(dgeall@meta.data[,ordered])
dgeall=SetAllIdent(dgeall,ordered)

# OldSca1+LM res.0.7order merge clusters4-5 together, merge clusters6-8 together
  dge=SetAllIdent(dge,"res.0.7order")
  ident=dge@ident
  table(dge@ident)
  ident[which(ident==5)]<-4
  ident[which(ident %in% c(6:8))]<-5
  ident[which(ident==9)]<-6
  ident[which(ident==10)]<-7
  ident=factor(ident,levels=1:7)
  table(ident)
  dge=AddMetaData(dge,ident,"res.0.7merge")
  dge=SetAllIdent(dge,"res.0.7merge")
  table(dge@ident)
  dgeall=dge
  dgealllist[[i]]=dgeall
# 1.31.2020 further merge clusters3 and 5 together
  dge=SetAllIdent(dge,"res.0.7merge")
  ident=dge@ident
  table(dge@ident)
  ident[which(ident==5)]<-3
  ident[which(ident==6)]<-5
  ident[which(ident==7)]<-6
  ident=factor(ident,levels=1:6)
  table(ident)
  dge=AddMetaData(dge,ident,"res.0.7merge2")
  dge=SetAllIdent(dge,"res.0.7merge2")
  table(dge@ident)
  dgeall=dge
  dgealllist[[i]]=dgeall

### save dgeall with ordered and merged clusters
save(dgeall,file=paste0(exp[i],".Robj"))

