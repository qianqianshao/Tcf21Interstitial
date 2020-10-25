###### Comparison between in-vitro leydig differentiation data and Adult Mouse Testis Soma by Qianyi on 3/13/2020
### Related to Figure S3A

### Adults Mouse Testis Soma (25 batches) -> 7 Somatic cell types by Qianyi in 2017, published in Green, Ma, et al., DevCell 2018
# Drop-seq
# data format: log(mean(counts-per-10k)+1), 26694 detected genes, 7 cell type centroids

### 7 clusters for OldSca1 (time point 0) + in-vitro leydig regeneration (3 time points)
# data format: log(mean(counts-per-10k)+1), 24698 detected genes, 7 cluster centroids

###### load our adult mouse testis somatic data by Qianyi in 2017
load(file ="data_DGE/MouseAdultST25Somatic.Robj")
dge1=dge
length(dge1@var.genes) # 2464

### load 7 cell type centroids for adult mouse testis somatic data
#data format: log(mean(counts-per-10k)+1), 26694 detected genes, 7 cell type centroids
soma=read.table("data_DGE/MouseAdultST25Somatic_7celltypes_centroid.txt",header=T,stringsAsFactors=F)
#soma=read.table("somaGonad4ClusterCentroids_RPKM.txt",header=T)
dim(soma) # 22734     7

### load markers for 7 cell types of adult mouse testis somatic data
somamarkers=read.table("data_DGE/MouseAdultST25Somatic_MarkersAll_7celltypes_pct0.2_diffpct0.2_thresh2fold_1.8.2018.txt",header=T,stringsAsFactors=F)
#soma=read.table("somaGonad4ClusterCentroids_RPKM.txt",header=T)
dim(somamarkers) # 1797     7
### number of markers for each cluster
table(somamarkers$cluster)
#   Endothelial InnateLymphoid         Leydig     Macrophage          Myoid 
#           226            388            218            248            116 
#       Sertoli        Unknown 
#           284            317 

###### load our OldSca1 + in-vitro leydig regeneration data by Qianyi in Dec 2019
load(file="DropseqSue_Oct2019/OldSca1LM.Robj")
dge2=dgeall
length(dge2@var.genes) # 2344

### load 7 cluster centroids for OldSca1+LM in-vitro leydig regeneration
# data format: log(mean(counts-per-10k)+1), 24698 genes, 7 cluster centroids (columns)
adult=read.table("DropseqSue_Oct2019/OldSca1LM_res.0.7merge_Centroid.txt",header=T)
dim(adult)   # [1] 24698     7
names(adult)=c("1.IntProgLy6a","2.IntProgTcf21","3.ProlifIntProg","4.KidneyCell","5.DiffIntProg","6.ImmLeydig","7.Leydig")

### load markers for 7 cell types of OldSca1+LM in-vitro leydig regeneration
adultmarkers=read.table("DropseqSue_Oct2019/exp_OldSca1LM_res.0.7merge_mindiff0.2_logfc2fold_11.2019.txt",header=T,stringsAsFactors=F)
#adult=read.table("adultGonad4ClusterCentroids_RPKM.txt",header=T)
dim(adultmarkers) # 617     7
### number of markers for each cluster
table(adultmarkers$cluster)
#  1   2   3   4   5   6   7 
#225 125  73  11   7  21 155 

###### Choose gene sets for cross-tabulation
markerlist=list()
### 1. use union of Highly-variable genes
markerlist[[1]]=unique(c(dge1@var.genes,dge2@var.genes))
### 2. use intersect of Highly-variable genes
markerlist[[2]]=intersect(dge1@var.genes,dge2@var.genes)
### 3. use union of markers 
markerlist[[3]]=unique(c(somamarkers$gene,adultmarkers$gene))

### Match overlapping genes
for(i in 1:3){
  markers=markerlist[[i]]
  markers=markers[which(markers %in% rownames(soma))]
  markers=markers[which(markers %in% rownames(adult))]
  markerlist[[i]]=markers
}

for(i in 1:3){
  print(length(markerlist[[i]]))
}
#[1] 4126
#[1] 662
#[1] 1769

### Cross-tabulate using rank correlation with each gene set
# Rank Correlation between the 7 cell type centroids for OldSca1 + in-vitro leydig differentiation and Adult Mouse Testis Soma (7 cell type centroids in 25 ST batches) 
rholist=list()
for(i in 1:3){
  markers=markerlist[[i]]
  all=cbind(soma[markers,],adult[markers,])
  rho=cor(all,method="spearman")
  rholist[[i]]=rho
  print(rho[8:14,1:7])
}
# saved as FigS3A


