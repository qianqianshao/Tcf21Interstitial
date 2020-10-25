# Pseudotemporal ordering of in-vitro Leydig differentiation data using Monocle3 by Qianyi in Dec 2019 
# Related to Figure 3D-E

setwd("/scratch/junzli_root/junzli/qzm/Dropseq_analysis/DropseqSue_Oct2019")
home="/scratch/junzli_root/junzli/qzm/Dropseq_analysis/"
dgefile=dgename="exp_"
exp="OldSCa1LM"
  load(file=paste0(exp[indiv],".Robj"))
  resSeurat="res.0.7merge"

library(monocle3)
library(ggplot2)
library(dplyr)
dge=dgeall
cells.use=names(dge@ident)
print(table(dge@ident[cells.use]))

umi_matrix=as.matrix(dge@raw.data[,cells.use])
cell_metadata=dge@meta.data[cells.use,]
gene_metadata=data.frame(gene_short_name=rownames(dge@data))
rownames(gene_metadata)=rownames(dge@data)

#Store data in a cell_data_set object
cds <- new_cell_data_set(as(umi_matrix, "sparseMatrix"),
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)

levels=as.numeric(sort(unique(pData(cds)[resSeurat]))[,1])
ntime=nrow(unique(pData(cds)["time2"]))

## Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 50) # default
pdf(paste0(exp,"_screeplot.pdf"),height=4,width=4.5)
plot_pc_variance_explained(cds)
dev.off()

## Reduce the dimensions using UMAP
cds <- reduce_dimension(cds,preprocess_method = 'PCA')
pdf(paste0(exp[indiv],"_umap.pdf"),height=2.5,width=3.5)
# visualize clusters from Seurat
plot_cells(cds, color_cells_by=resSeurat[indiv],
           label_cell_groups = FALSE,
           show_trajectory_graph=FALSE,
           graph_label_size=3) + scale_color_manual(values=myBrewerPalette[levels])
dev.off()

## Cluster the cells
res=1e-6
	cds = cluster_cells(cds, resolution=res)
pdf(paste0(exp,"_cluster.pdf"),height=2.5,width=2.5)
# visualize clusters from Monocle3
plot_cells(cds, color_cells_by="cluster", group_cells_by="cluster",
           show_trajectory_graph=FALSE,
           group_label_size=5)
dev.off()

## learn trajectory
cds <- learn_graph(cds)
pdf(paste0(exp,"_graph.pdf"),height=2.5,width=2.5)
plot_cells(cds,
           color_cells_by = resSeurat,
           label_groups_by_cluster=FALSE,
           group_label_size=5,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5) + scale_color_manual(values=myBrewerPalette[levels])
dev.off()
# saved as Fig3D

## Order cells
cds <- order_cells(cds)
pdf(paste0(exp,"_pseudotime.pdf"),height=2.5,width=3.5)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=3)
dev.off()

### plot known markers in pseudotime
knownmarkers=c("Pdgfra","Tcf21","Ly6a","Hsd3b1","Cyp17a1","Cyp11a1","Star","Insl3","Thra","Thrb","Acta2","Myh11","Clu","Kit","Mki67","Nr2f2")
knownmarkers_cds <- cds[rowData(cds)$gene_short_name %in% knownmarkers,
                       ]
jpeg(paste0(exp,"_pseudotime_knownmarkers.jpeg"),res=300,height=1200,width=2500)
plot_genes_in_pseudotime(knownmarkers_cds,
                         color_cells_by=resSeurat[indiv],
                         min_expr=0.5,ncol=4) + scale_color_manual(values=myBrewerPalette[levels])
dev.off()
## visualize more markers - 3/13/2020
knownmarkers=read.table("3.12.20plotmarkers.txt",stringsAsFactors=F)[,1]
length(knownmarkers) # 37
knownmarkers[which(!(knownmarkers %in% rowData(cds)$gene_short_name))]
knownmarkers_cds <- cds[rowData(cds)$gene_short_name %in% knownmarkers,
                       ]
jpeg(paste0(exp,"_pseudotime_more_markers.jpeg"),res=300,height=1800,width=3750)
plot_genes_in_pseudotime(knownmarkers_cds,
                         color_cells_by=resSeurat[indiv],
                         min_expr=0.5,ncol=6) + scale_color_manual(values=myBrewerPalette[levels])
dev.off()
# saved in Fig3E
