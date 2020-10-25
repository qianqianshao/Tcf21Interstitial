###### Comparison between in-vitro leydig differentiation data and Fetal Gonad by Qianyi on 1/15/2020
### Related to Figure S3B

### Stevant Cell Reports 2018 for Nr5a1+ Male Fetal Gonad 
# Nr5a1-GFP+ cells, Fluidigm C1 single-cell RNA-seq
# https://www.sciencedirect.com/science/article/pii/S2211124718300755?via%3Dihub

### 7 clusters for OldSca1 (time point 0) + in-vitro leydig regeneration (3 time points)
# data format: log(mean(counts-per-10k)+1), 24698 detected genes, 7 cluster centroids

home="C:/Users/qzm.UMHS/Desktop/Leydig-Regeneration/fig_LeyReg/IndivExp_all3exp_SomaMergedTimeTreat/LM_CrossTabulation/"
setwd(home)
setwd("ComparisonWithLM_Stevant_Gonad_2018/")

### load Stevant 2018 Fetal Gonad data
# data format: mean(RPKM), 4435 markers (rows), 6 cluster centroids (columns)
fetal=read.csv("FetalGonad6ClusterCentroids_mmc2.csv",header=T,stringsAsFactors=F)
#fetal=read.table("FetalGonad4ClusterCentroids_RPKM.txt",header=T)
dim(fetal) # 2802 11
anyDuplicated(fetal[,1]) # duplicted marker for each cluster
fetal[anyDuplicated(fetal[,1]),1] 
rownames(fetal)=fetal[,1]
fetal=fetal[,-1]
names(fetal)[1:7]=c("C1.Endothelial","C2.EarlyProg","C3.IntProg","C4.Pre.Sertoli","C5.FetalLeydig","C6.Sertoli","ClusterID")
### number of markers for each cluster
table(fetal$ClusterID)
# C1  C2  C3  C4  C5  C6 
#377 602 625 389 329 480 
dim(fetal) # [1] 2802   10
fetal[1:2,]
#        C1.Endothelial C2.EarlyProg C3.IntProg C4.Pre.Sertoli C5.FetalLeydig C6.Sertoli
#Afap1l1       56.46611    0.5203821 0.07745985      0.1990415              0 0.09825918
#Akr1c14       53.25356    0.6311708 0.02809066      0.0000000              0 0.27971679
#        ClusterID num_cells_expressed pval qval
#Afap1l1        C1                  14    0    0
#Akr1c14        C1                  11    0    0

summary(c(as.matrix(fetal[,1:6])))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#   0.000    0.382    5.093   32.527   20.009 7764.770 

###### load the 7 cluster centroids for OldSca1 + in-vitro leydig regeneration data
### load cluster centroids
# data format: log(mean(counts-per-10k)+1), 24698 genes, 7 cluster centroids (columns)
adult=read.table("OldSca1LM_res.0.7merge_Centroid.txt",header=T)
dim(adult)   # [1] 24698     7
summary(c(as.matrix(adult)))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.00000 0.01171 0.18291 0.18601 5.30190
### change data format to mean(counts-per-1M) for our 7 clusters
adult<-expm1(adult)*100
summary(c(as.matrix(adult)))
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#    0.000     0.000     1.178    40.489    20.444 19971.728 
names(adult)=c("1.IntProgLy6a","2.IntProgTcf21","3.ProlifIntProg","4.KidneyCell","5.DiffIntProg","6.ImmLeydig","7.Leydig")
adult[1:2,]
#              1.IntProgLy6a 2.IntProgTcf21 3.ProlifIntProg 4.KidneyCell 5.DiffIntProg
#0610007P14Rik       0.00000       78.40848        89.45929     96.79621      83.65371
#0610009B22Rik      33.01358       40.15885        52.14396     58.87424      42.74908
#              6.ImmLeydig  7.Leydig
#0610007P14Rik   105.48432 211.00790
#0610009B22Rik    57.20171  87.74432

### Match overlapped markers
overlapped=rownames(fetal[rownames(fetal) %in% rownames(adult),])
length(overlapped) # [1] 2692
fetal<-fetal[overlapped,]
table(fetal$ClusterID)
# C1  C2  C3  C4  C5  C6 
#373 567 604 367 323 458
adult<-adult[overlapped,]
all(rownames(adult) == rownames(fetal)) #TRUE
summary(adult)
 1.IntProgLy6a     2.IntProgTcf21      3.ProlifIntProg     4.KidneyCell      
 Min.   :   0.00   Min.   :    0.000   Min.   :   0.000   Min.   :    0.000  
 1st Qu.:   0.00   1st Qu.:    3.612   1st Qu.:   3.433   1st Qu.:    2.694  
 Median :  20.37   Median :   20.701   Median :  18.515   Median :   14.910  
 Mean   :  79.07   Mean   :   91.399   Mean   :  87.678   Mean   :   86.986  
 3rd Qu.:  70.81   3rd Qu.:   63.817   3rd Qu.:  60.968   3rd Qu.:   52.235  
 Max.   :8435.01   Max.   :11434.148   Max.   :9125.539   Max.   :12832.019  
 5.DiffIntProg        6.ImmLeydig           7.Leydig        
 Min.   :    0.000   Min.   :    0.000   Min.   :    0.000  
 1st Qu.:    2.606   1st Qu.:    3.226   1st Qu.:    2.177  
 Median :   14.827   Median :   16.030   Median :   12.658  
 Mean   :   88.204   Mean   :   88.541   Mean   :   91.050  
 3rd Qu.:   53.345   3rd Qu.:   53.157   3rd Qu.:   44.404  
 Max.   :10382.483   Max.   :14423.992   Max.   :13416.779 
summary(fetal[,1:4])
 C1.Endothelial      C2.EarlyProg        C3.IntProg       C4.Pre.Sertoli     
 Min.   :   0.000   Min.   :   0.000   Min.   :   0.000   Min.   :   0.0000  
 1st Qu.:   0.000   1st Qu.:   2.064   1st Qu.:   2.469   1st Qu.:   0.7385  
 Median :   0.057   Median :   7.647   Median :   8.997   Median :   5.8290  
 Mean   :  29.683   Mean   :  26.548   Mean   :  28.617   Mean   :  31.5718  
 3rd Qu.:  11.979   3rd Qu.:  20.191   3rd Qu.:  23.451   3rd Qu.:  21.4563  
 Max.   :6527.126   Max.   :4514.734   Max.   :2904.171   Max.   :2892.7538  

###### calculate rank correlation between Stevant fetal and in-vitro leydig differentiation cluster centroids using overlapped markers
all=cbind(fetal[,1:6],adult)
all[1:2,]
#        C1.Endothelial C2.EarlyProg C3.IntProg C4.Pre.Sertoli C5.FetalLeydig C6.Sertoli
#Afap1l1       56.46611    0.5203821 0.07745985      0.1990415              0 0.09825918
#Akr1c14       53.25356    0.6311708 0.02809066      0.0000000              0 0.27971679
#        1.IntProgLy6a 2.IntProgTcf21 3.ProlifIntProg 4.KidneyCell 5.DiffIntProg 6.ImmLeydig
#Afap1l1     119.74328       3.469532        4.423441    0.1462779      1.179371     0.00000
#Akr1c14      59.32737      59.258125       10.515709    2.4109137     19.846068    37.26798
#        7.Leydig
#Afap1l1 0.000000
#Akr1c14 4.312793

rho=cor(all,method="spearman")
summary(c(rho[7:13,1:6]))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3452  0.4087  0.4557  0.4615  0.4936  0.6412 
rho[7:13,1:6]
                C1.Endothelial C2.EarlyProg C3.IntProg C4.Pre.Sertoli C5.FetalLeydig C6.Sertoli
1.IntProgLy6a        0.4921033    0.3489701  0.3867597      0.3476159      0.3452332  0.3599595
2.IntProgTcf21       0.4151423    0.4075939  0.5612188      0.3950201      0.4585931  0.3935706
3.ProlifIntProg      0.4579267    0.5392240  0.6412484      0.3917628      0.4751843  0.4522093
4.KidneyCell         0.4577483    0.4940439  0.5754167      0.4264225      0.4524019  0.4118348
5.DiffIntProg        0.4612894    0.4885170  0.5976840      0.4132407      0.4532031  0.4033228
6.ImmLeydig          0.4333398    0.4642210  0.6148058      0.4543325      0.5282606  0.4569868
7.Leydig             0.3944326    0.4182602  0.5459248      0.4883753      0.5834696  0.4954776
# saved as FigS3B
