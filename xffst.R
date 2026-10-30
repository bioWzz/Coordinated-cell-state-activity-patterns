install.packages("clustree")

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(clustree)
library(patchwork)

#测试一套SCAR_ST_000007
SCAR_ST_000007 <- readRDS("D:/细胞状态2/论文/空间/paad/Pancreatic ductal adenocarcinoma/SCAR_ST_000007/other/SCAR_ST_000007.rds")
SCAR_ST_000007 <- readRDS("/wangzhenzenlab/slb/st_data/Pancreatic ductal adenocarcinoma/SCAR_ST_000007/other/SCAR_ST_000007.rds")


SpatialDimPlot(SCAR_ST_000007,alpha=1)

plot1 <- VlnPlot(SCAR_ST_000007, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(SCAR_ST_000007, features = "nCount_Spatial") + theme(legend.position = "right")
p = wrap_plots(plot1, plot2)

#降维
SCAR_ST_000007 <- RunPCA(SCAR_ST_000007, assay = "SCT", verbose = FALSE)
SCAR_ST_000007 <- FindNeighbors(SCAR_ST_000007, reduction = "pca", dims = 1:30)
SCAR_ST_000007 <- FindClusters(SCAR_ST_000007,resolution = 1, verbose = FALSE)
SCAR_ST_000007 <- RunUMAP(SCAR_ST_000007, reduction = "pca", dims = 1:20)

ElbowPlot(SCAR_ST_000007)

DimPlot(SCAR_ST_000007, reduction = "umap", label = TRUE)
SpatialDimPlot(SCAR_ST_000007, label = TRUE, label.size = 3, alpha =1, pt.size = 1)
  #scale_colour_manual(values=c("stromal"="gray60","Cancer"="#d9352a"))

#将位置信息合并入单细胞seurat对象
position<- SCAR_ST_000007@images[["slice1"]]@coordinates
position=select(position,imagecol,imagerow)
head(position)
colnames(position)=paste0("Spatial_",1:ncol(position))

SCAR_ST_000007[["spatial"]]<-CreateDimReducObject(embeddings=as.matrix(position),
                                             key="Spatial_",assay="Spatial")
DimPlot(SCAR_ST_000007,reduction="spatial",pt.size=2)


#提取表达量变化最高的10个基因
top10<- head(VariableFeatures(SCAR_ST_000007),10)
top10
p5 <- VariableFeaturePlot(SCAR_ST_000007,cols = c( "gray60", "red"))+NoLegend()
p6 <- LabelPoints(plot = p5,points = top10, repel = TRUE,xnudge=0,ynudge=0)
p5 | p6


ST_pbmc.markers <-  FindAllMarkers(SCAR_ST_000007, only.pos = TRUE,min.pct = 0.25)
top10 <- ST_pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

top10$gene

BiocManager::install("SingleR")
BiocManager::install("celldex")

library(SingleR)
library(celldex)

ref_data <- celldex::HumanPrimaryCellAtlasData()
ref_us <- celldex::MouseRNAseqData()

expression_data <- SCAR_ST_000007[["Spatial"]]@data
expression_data[1:3,1:3]

singleR_results <- SingleR(
  test = expression_data,  # 空间转录组数据的表达矩阵
  ref = ref_data,          # 参考数据集
  labels = ref_data$label.main , # 选择细胞类型标签
  clusters = SCAR_ST_000007$seurat_clusters,            #细胞聚类信息（只有method = "cluster"才使用）
  method = "cluster"
)

#SingleR注释
#pred <- SingleR(test=as.matrix(pbmc@assays$RNA@data),       #输入表达矩阵
##                ref=ref_use,                                #参考数据
#                labels=ref_use$label.fine,                  #标签列
#                clusters = pbmc$seurat_clusters,            #细胞聚类信息（只有method = "cluster"才使用）
#                method = "cluster")                         #注释类型，按照cluster还是按照cell进行注释

SCAR_ST_000007$SingleR_labels <- singleR_results$labels

SpatialDimPlot(SCAR_ST_000007, group.by = "SingleR_labels", label = TRUE, label.size = 3)









#测试二套
#setwd("D:/细胞状态2/论文/空间/paad/nfs")
setwd("/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/nfs")
seurat_data=Load10X_Spatial(data.dir="./SS1905133_processed",
                          filename="filtered_feature_bc_matrix.h5",
                          assay="Spatial",
                          slice="test")
class(seurat_data)

seurat_data <- NormalizeData(seurat_data, normalization.method = "LogNormalize", scale.factor = 10000)
#seurat_data <- NormalizeData(seurat_data)
#识别高变基因
seurat_data <- FindVariableFeatures(seurat_data, selection.method = "vst", nfeatures = 2000)
#找出排在前10个的高变基因
top10 <- head(VariableFeatures(seurat_data), 10)
top10
#变量特征可视化
#VariableFeaturePlot(seurat_data)
#LabelPoints(plot = plot1, points = top10, repel = TRUE)

#数据归一化
all.genes <- rownames(seurat_data)
seurat_data <- ScaleData(seurat_data, features = all.genes)
#线性降维PCA主成分分析
seurat_data <- RunPCA(seurat_data, features = VariableFeatures(object = seurat_data))
ElbowPlot(seurat_data)

#--用KNN方法对细胞聚类----------#
seurat_data <- FindNeighbors(seurat_data, dims = 1:20)
seurat_data <- FindClusters(seurat_data, resolution =0.6)
head(Idents(seurat_data), 5)
#非线性降维(UMAP/tSNE)
pdf("/wangzhenzenlab/fsh/xff/st2/umap/umap.pdf")
seurat_data <- RunUMAP(seurat_data, dims = 1:20)
DimPlot(seurat_data, reduction = "umap",label = TRUE, pt.size = 1)
dev.off()
seurat_data <- RunTSNE(seurat_data, dims = 1:20)
DimPlot(seurat_data, reduction = "tsne",label = TRUE, pt.size = 1)

position=read.csv("./SS1905133_processed/spatial/tissue_positions_list.csv",header=F,row.names=1)
head(position)
position=select(position,V5,V6)
colnames(position)=paste0("Spatial_",1:ncol(position))
#seurat_data[["spatial"]]<-CreateDimReducObject(embeddings=as.matrix(position),
#                                             key="Spatial",assay="RNA")


cells_select <- intersect(rownames(position), colnames(seurat_data))
dim(seurat_data)
dim(position)
position2 <- position[cells_select,]

seurat_data[["spatial"]]<-CreateDimReducObject(embeddings=as.matrix(position2),
                                               key="Spatial_",assay="RNA")
DimPlot(seurat_data,reduction="spatial",pt.size=1)
SpatialFeaturePlot(seurat_data, features = "nCount_Spatial") + theme(legend.position = "right")

head(seurat_data@meta.data)
SpatialDimPlot(seurat_data,group.by = "celltype",alpha=1)

#---数据的预处理-----#
VlnPlot(seurat_data, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
SpatialFeaturePlot(seurat_data, features = "nCount_Spatial") + theme(legend.position = "right")
#-------------#
VlnPlot(seurat_data,  features = "nFeature_Spatial",pt.size = 0.1,cols = "tomato",group.by = "orig.ident") + NoLegend()
FeatureScatter(seurat_data, feature1 = "nCount_Spatial",feature2= "nFeature_Spatial",group.by = "orig.ident")+ NoLegend()
#--------------------#
BiocManager::install("glmGamPoi")
install.packages("future.apply")
library(future)
library(glmGamPoi)

seurat_data <- SCTransform(seurat_data, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
#提取表达量变化最高的10个基因
top10<- head(VariableFeatures(seurat_data),10)
top10
p5 <- VariableFeaturePlot(seurat_data,cols = c( "gray60", "red"))+NoLegend()
p6 <- LabelPoints(plot = p5,points = top10, repel = TRUE,xnudge=0,ynudge=0)
p5 | p6




seurat_data.markers <-  FindAllMarkers(seurat_data, only.pos = TRUE,min.pct = 0.25)
top10 <- seurat_data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10$gene
table(top10$cluster)
write.csv(top10,file = "/wangzhenzenlab/fsh/xff/st2/top10.csv")

getwd()
save(seurat_data,file = "/wangzhenzenlab/fsh/xff/st2/seurat_data.RData")
load()


genes <- c("EPCAM", "KRT19", "CDH1", "KRT18")            #Epithelial
genes <- c("CD3D", "CD3E", "CD3G", "TRAC")               #T
genes <- c("CD79A", "IGHM", "IGHG3", "IGHA2")            #B
genes <- c("CD68", "MARCO", "FCGR3A", "CD16", "LYZ")     #Myeloid
genes <- c("NCAM1","CD56", "NKG7", "GNLY", "KLRD1")      #NK
genes <- c("KIT", "MS4A2","GATA2")                       #MAST
genes <- c("DCN", "COL1A1", "COL1A2", "THY1")            #Fibroblasts
genes <- c("PECAM1", "CLDN5", "FLT1", "RAMP2")           #Endothelial
genes <- c("EGFR", "KRAS", "TP53", "CDKN2A","VEGFA")
genes <- c("SMAD4", "MUC1", "CA19-9", "VEGF","CXCR4")

pdf("/wangzhenzenlab/fsh/xff/st2/Spatial/Spatial.pdf")
SpatialDimPlot(seurat_data,alpha=1)
dev.off()

seurat_data@meta.data
seurat_data$seurat_clusters
getwd()
#绘制umap
pdf("/wangzhenzenlab/fsh/xff/st2/umap/Endothelial.pdf")
FeaturePlot(seurat_data, features = genes, reduction = "umap")
dev.off()
# 绘制气泡图
pdf("/wangzhenzenlab/fsh/xff/st2/Endothelial.pdf")
DotPlot(seurat_data, features = genes, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
pdf("/wangzhenzenlab/fsh/xff/st2/Spatial/Endothelial.pdf")
SpatialFeaturePlot(seurat_data, features = genes)
dev.off()
0:Endothelial/Fibroblasts
1:Epithelial
2:Epithelial
3:Myeloid
4:B
5:Epithelial
6:T/NK
7:Fibroblasts
8:Epithelial
9:Fibroblasts
10:Epithelial

#GSVA
package.version('GSVA')
library(GSVA)
gene_sets <- read.csv("D:/细胞状态2/cell_state.csv",header = T)
head(gene_sets)
list=split(gene_sets$GeneName,gene_sets$Signature )

expression_matrix <- GetAssayData(seurat_data, assay = "Spatial", slot = "data")
dim(expression_matrix)
deal <- as.data.frame(expression_matrix)
deal[1:3,1:3]

dealT <- t(deal)
dealT2 <- dealT[,c("ACTA2" , "AEBP1" , "ANXA2","ATF3","BGN","CALD1","COL1A1","CTSK",
                   "EGR1","EGR3","FN1","HTRA1","INHBA","MMP11","MMP14","MMP2","NR4A1","PDGFA",
                   "S100A4","SPP1","VCAN")]
head(dealT2)
write.csv(dealT2,file = "D:/细胞状态2/论文/空间/paad/dealT2.csv")
dealT2 <- read.csv("D:/细胞状态2/论文/空间/paad/paaddealT2 (2).csv")


ES <- gsva(expr=deal, 
           gset.idx.list=list,
           kcdf=c("Gaussian"),
           parallel.sz=32)


load("D:/细胞状态2/论文/空间/st2/seurat_data.RData")
load("D:/细胞状态2/论文/空间/st/EST.RData")
load("/wangzhenzenlab/fsh/xff/st2/seurat_data.RData")
load("/wangzhenzenlab/fsh/xff/st2/seurat_data.RData")
head(EST)
EST <- as.data.frame(EST)
EST$cell <- rownames(EST)

stmeta.data<- seurat_data@meta.data
head(stmeta.data)
stmeta.data$cell <- rownames(stmeta.data)
stmeta.data <- merge(stmeta.data,EST,by="cell")
stmeta.data$celltype <- ifelse(stmeta.data$seurat_clusters==0,"Endothelial",
                               ifelse(stmeta.data$seurat_clusters%in%c(1,2,5,8,10),"Epithelial",
                                      ifelse(stmeta.data$seurat_clusters==3,"Myeloid",
                                             ifelse(stmeta.data$seurat_clusters==4,"B",
                                                    ifelse(stmeta.data$seurat_clusters==6,"T/NK",
                                                           ifelse(stmeta.data$seurat_clusters%in%c(7,9),"Fibroblasts","NA"))))))
table(stmeta.data$seurat_clusters)
table(stmeta.data$celltype)




singlemt <- stmeta.data[,7:20]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
my_colors <- colorRampPalette(c("steelblue", "white", "#c11e21"))(100)
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)


stmeta.data$model1 <- rowMeans(stmeta.data[,c("Angiogenesis","EMT","Invasion","Metastasis")])
stmeta.data$group1 <- ifelse(stmeta.data$model1>median(stmeta.data$model1),"high","low")
table(stmeta.data$group1,stmeta.data$celltype)
stmeta.data$region <- ifelse(stmeta.data$celltype=="Epithelial","tumor region","other region")
table(stmeta.data$group1,stmeta.data$region)
head(stmeta.data)
rownames(stmeta.data)<- stmeta.data$cell

dim(stmeta.data)
dim(seurat_data@meta.data)

seurat_cells <- rownames(seurat_data@meta.data)
model1_cells <- rownames(stmeta.data)

setdiff(seurat_cells, model1_cells)
setdiff(model1_cells, seurat_cells)

table(seurat_data$group1)

seurat_data@meta.data <- stmeta.data
table(seurat_data$group1)
pdf("/wangzhenzenlab/fsh/xff/st2/group1.pdf")
SpatialDimPlot(seurat_data,alpha=1,group.by = "group1",cols = c("high"="#c11e21","low"="steelblue"))
dev.off()
pdf("/wangzhenzenlab/fsh/xff/st2/region.pdf")
SpatialDimPlot(seurat_data,alpha=1,group.by = "region")
dev.off()
write.csv(stmeta.data,file = "D:/细胞状态2/论文/空间/st2/stmeta.data.csv")
stmeta.data <- read.csv("D:/细胞状态2/论文/空间/st2/stmeta.data.csv")

head(seurat_data)
SpatialDimPlot(seurat_data,alpha=1,group.by = "group1")+ stat_contour(aes(z = seurat_data$model1), color = "black", size = 1)


pdf("/wangzhenzenlab/fsh/xff/st2/paadcelltype.pdf")
SpatialDimPlot(seurat_data,group.by = "celltype",alpha=1)
dev.off()

head(dealT2)
head(stmeta.data)
dealT3 <- merge(dealT2,stmeta.data,by.x="X",by.y="cell")
head(dealT3)

singlemt <- dealT3[,c("EGR3","EGR1","ANXA2","ATF3","PDGFA","model1")]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
library(corrplot)
my_colors <- colorRampPalette(c("steelblue", "white", "#c11e21"))(100)
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)










#STAD第一套
#setwd("D:/细胞状态2/论文/空间/stad")
setwd("/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/stad")
stadseurat_data=Load10X_Spatial(data.dir="./21_00732_LI_SING",
                                filename="21_00732_LI_SING_filtered_feature_bc_matrix.h5",
                                assay="Spatial",
                                slice="test")
class(stadseurat_data)

head(stadseurat_data@meta.data)

stadseurat_data <- NormalizeData(stadseurat_data, normalization.method = "LogNormalize", scale.factor = 10000)

#识别高变基因
stadseurat_data <- FindVariableFeatures(stadseurat_data, selection.method = "vst", nfeatures = 2000)
#找出排在前10个的高变基因
top10 <- head(VariableFeatures(stadseurat_data), 10)
top10




#数据归一化
all.genes <- rownames(stadseurat_data)
stadseurat_data <- ScaleData(stadseurat_data, features = all.genes)
#线性降维PCA主成分分析
stadseurat_data <- RunPCA(stadseurat_data, features = VariableFeatures(object = stadseurat_data))
ElbowPlot(stadseurat_data)
#seurat_data
#--用KNN方法对细胞聚类----------#
stadseurat_data <- FindNeighbors(stadseurat_data, dims = 1:20)
stadseurat_data <- FindClusters(stadseurat_data, resolution =0.6)
head(Idents(stadseurat_data), 5)
#非线性降维(UMAP/tSNE)
pdf("/wangzhenzenlab/fsh/xff/stad/umap/umap.pdf")
stadseurat_data <- RunUMAP(stadseurat_data, dims = 1:20)
DimPlot(stadseurat_data, reduction = "umap",label = TRUE, pt.size = 1)
dev.off()
#stadseurat_data <- RunTSNE(stadseurat_data, dims = 1:20)
#DimPlot(stadseurat_data, reduction = "tsne",label = TRUE, pt.size = 1)

position=read.csv("./21_00732_LI_SING/spatial/tissue_positions_list.csv",header=F,row.names=1)
head(position)
class(position)
colnames(position)

position <- position[ ,c("V5", "V6")]
colnames(position)=paste0("Spatial_",1:ncol(position))

cells_select <- intersect(rownames(position), colnames(stadseurat_data))
dim(stadseurat_data)
dim(position)
position2 <- position[cells_select,]

stadseurat_data[["spatial"]]<-CreateDimReducObject(embeddings=as.matrix(position2),
                                               key="Spatial_",assay="RNA")
DimPlot(stadseurat_data,reduction="spatial",pt.size=1)
SpatialFeaturePlot(stadseurat_data, features = "nCount_Spatial") + theme(legend.position = "right")

head(stadseurat_data@meta.data)
SpatialDimPlot(stadseurat_data,alpha=1)



BiocManager::install("glmGamPoi")
install.packages("future.apply")
library(future)
library(glmGamPoi)


stadseurat_data <- SCTransform(stadseurat_data, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)



stadseurat_data.markers <-  FindAllMarkers(stadseurat_data, only.pos = TRUE,min.pct = 0.25)
top10 <- stadseurat_data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10$gene
table(top10$cluster)
write.csv(top10,file = "/wangzhenzenlab/fsh/xff/stad/top10.csv")

stadseurat_data0732 <- stadseurat_data
save(stadseurat_data,file = "/wangzhenzenlab/fsh/xff/stad/stadseurat_data0732.RData")
load("/wangzhenzenlab/fsh/xff/stad/stadseurat_data.RData")


genes <- c("EPCAM", "KRT19", "CDH1", "KRT18")            #Epithelial
genes <- c("CD3D", "CD3E", "CD3G", "TRAC")               #T
genes <- c("CD79A", "IGHM", "IGHG3", "IGHA2")            #B
genes <- c("CD68", "MARCO", "FCGR3A", "CD16", "LYZ")     #Myeloid
genes <- c("NCAM1","CD56", "NKG7", "GNLY", "KLRD1")      #NK
genes <- c("KIT", "MS4A2","GATA2")                       #MAST
genes <- c("DCN", "COL1A1", "COL1A2", "THY1")            #Fibroblasts
genes <- c("PECAM1", "CLDN5", "FLT1", "RAMP2")           #Endothelial
genes <- c("EGFR", "KRAS", "TP53", "CDKN2A","CD44","CDH1")    #cancer


pdf("/wangzhenzenlab/fsh/xff/stad/Spatial/Spatial.pdf")
SpatialDimPlot(stadseurat_data,alpha=1)
dev.off()

stadseurat_data@meta.data
stadseurat_data$seurat_clusters
getwd()
#绘制umap
pdf("/wangzhenzenlab/fsh/xff/stad/umap/Epithelial.pdf")
FeaturePlot(stadseurat_data, features = genes, reduction = "umap")

dev.off()
# 绘制气泡图
pdf("/wangzhenzenlab/fsh/xff/stad/Epithelial.pdf")
DotPlot(stadseurat_data, features = genes, group.by = "seurat_clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
pdf("/wangzhenzenlab/fsh/xff/stad/Spatial/Epithelial.pdf")
SpatialFeaturePlot(stadseurat_data, features = genes)
dev.off()
0:B,Myeloid
1:mast
2:Epithelial
3:Epithelial
4:Epithelial
5:Fibroblasts
6:Endothelial
7:Epithelial
8:Epithelial
9:T
10:Fibroblasts
11:B,T


#GSVA
package.version('GSVA')
library(GSVA)
gene_sets <- read.csv("D:/细胞状态2/cell_state.csv",header = T)
head(gene_sets)
list=split(gene_sets$GeneName,gene_sets$Signature )

expression_matrix <- GetAssayData(stadseurat_data, assay = "Spatial", slot = "data")
dim(expression_matrix)
deal <- as.data.frame(expression_matrix)
deal[1:3,1:3]



ES <- gsva(expr=deal, 
           gset.idx.list=list,
           kcdf=c("Gaussian"),
           parallel.sz=32)


load("D:/细胞状态2/论文/空间/stad/stadseurat_data0732.RData")
load("D:/细胞状态2/论文/空间/stad/stadEST0732.RData")
load("/wangzhenzenlab/fsh/xff/stad/stadseurat_data.RData")
load("/wangzhenzenlab/fsh/xff/stad/stadEST.RData")
head(stadEST)
stadEST <- as.data.frame(stadEST)
stadEST$cell <- rownames(stadEST)

stmeta.data<- stadseurat_data@meta.data
head(stmeta.data)
stmeta.data$cell <- rownames(stmeta.data)
stmeta.data <- merge(stmeta.data,stadEST,by="cell")
stmeta.data$celltype <- ifelse(stmeta.data$seurat_clusters==0,"B/Myeloid",
                               ifelse(stmeta.data$seurat_clusters%in%c(2,3,4,7,8),"Epithelial",
                                      ifelse(stmeta.data$seurat_clusters==1,"Mast",
                                             ifelse(stmeta.data$seurat_clusters==6,"Endothelial",
                                                    ifelse(stmeta.data$seurat_clusters%in%c(9,11),"T",
                                                           ifelse(stmeta.data$seurat_clusters%in%c(5,10),"Fibroblasts","NA"))))))
table(stmeta.data$seurat_clusters)
table(stmeta.data$celltype)

SpatialDimPlot(stadseurat_data,alpha=1,group.by = "celltype")

0:B,Myeloid
1:mast
2:Epithelial
3:Epithelial
4:Epithelial
5:Fibroblasts
6:Endothelial
7:Epithelial
8:Epithelial
9:T
10:Fibroblasts
11:B,T


singlemt <- stmeta.data[,9:22]
head(singlemt)
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
my_colors <- colorRampPalette(c("steelblue", "white", "#c11e21"))(100)
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)


stmeta.data$model1 <- rowMeans(stmeta.data[,c("Angiogenesis","EMT","Invasion","Metastasis")])
stmeta.data$group1 <- ifelse(stmeta.data$model1>median(stmeta.data$model1),"high","low")
table(stmeta.data$group1,stmeta.data$celltype)
stmeta.data$region <- ifelse(stmeta.data$celltype=="Epithelial","tumor region","other region")
table(stmeta.data$group1,stmeta.data$region)
head(stmeta.data)
rownames(stmeta.data)<- stmeta.data$cell

dim(stmeta.data)
dim(stadseurat_data@meta.data)

seurat_cells <- rownames(stadseurat_data@meta.data)
model1_cells <- rownames(stmeta.data)

setdiff(seurat_cells, model1_cells)
setdiff(model1_cells, seurat_cells)

stadseurat_data@meta.data <- stmeta.data
table(stadseurat_data$group1)

dim(stadseurat_data@meta.data)
dim(stmeta.data)

pdf("/wangzhenzenlab/fsh/xff/stad/group1.pdf")
SpatialDimPlot(stadseurat_data,alpha=1,group.by = "group1",cols = c("high"="#c11e21","low"="steelblue"))
dev.off()
pdf("/wangzhenzenlab/fsh/xff/stad/region.pdf")
SpatialDimPlot(stadseurat_data,alpha=1,group.by = "region")
dev.off()
write.csv(stmeta.data,file = "D:/细胞状态2/论文/空间/stad/stmeta.data.csv")

head(stadseurat_data)
SpatialDimPlot(stadseurat_data,alpha=1,group.by = "group1")+ stat_contour(aes(z = stadseurat_data$model1), color = "black", size = 1)


dealT <- t(deal)
dealT2 <- dealT[,c("SELE","CXCL12","EGR1","PDPN","CD44","COL6A2","LGALS1","FBLN5")]
head(dealT2)
write.csv(dealT2,file = "D:/细胞状态2/论文/空间/stad/dealT2.csv")

save(stadseurat_data,file = "D:/细胞状态2/论文/空间/stad/stadseurat_data.RData")
load("/wangzhenzenlab/fsh/xff/stad/stadseurat_data.RData")

pdf("/wangzhenzenlab/fsh/xff/stad/stadcelltype.pdf")
SpatialDimPlot(stadseurat_data,group.by = "celltype",alpha=1)
dev.off()










#STAD
#setwd("D:/细胞状态2/论文/空间/stad")
setwd("/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/stad")
stadseurat_data0331=Load10X_Spatial(data.dir="./20_00331_LI_SING",
                                filename="20_00331_LI_SING_filtered_feature_bc_matrix.h5",
                                assay="Spatial",
                                slice="test")

stadseurat_data0731=Load10X_Spatial(data.dir="./21_00731_LI_SING",
                                filename="21_00731_LI_SING_filtered_feature_bc_matrix.h5",
                                assay="Spatial",
                                slice="test")

stadseurat_data0733=Load10X_Spatial(data.dir="./21_00733_LI_SING",
                                filename="21_00733_LI_SING_filtered_feature_bc_matrix.h5",
                                assay="Spatial",
                                slice="test")

stadseurat_data0734=Load10X_Spatial(data.dir="./21_00734_LI_SING",
                                filename="21_00734_LI_SING_filtered_feature_bc_matrix.h5",
                                assay="Spatial",
                                slice="test")

stadseurat_data1251=Load10X_Spatial(data.dir="./21_01251_LI_SING",
                                filename="21_01251_LI_SING_filtered_feature_bc_matrix.h5",
                                assay="Spatial",
                                slice="test")

stadseurat_data1252=Load10X_Spatial(data.dir="./21_01252_LI_SING",
                                filename="21_01252_LI_SING_filtered_feature_bc_matrix.h5",
                                assay="Spatial",
                                slice="test")

stadseurat_data1253=Load10X_Spatial(data.dir="./21_01253_LI_SING",
                                filename="21_01253_LI_SING_filtered_feature_bc_matrix.h5",
                                assay="Spatial",
                                slice="test")

stadseurat_data1254=Load10X_Spatial(data.dir="./21_01254_LI_SING",
                                filename="21_01254_LI_SING_filtered_feature_bc_matrix.h5",
                                assay="Spatial",
                                slice="test")

stadseurat_data1675=Load10X_Spatial(data.dir="./21_01675_LI_SING",
                                filename="21_01675_LI_SING_filtered_feature_bc_matrix.h5",
                                assay="Spatial",
                                slice="test")


stmetadata0331
#GSVA
library(GSVA)
gene_sets <- read.csv("D:/细胞状态2/cell_state.csv",header = T)
head(gene_sets)
list=split(gene_sets$GeneName,gene_sets$Signature )

expression_matrix0331 <- GetAssayData(stadseurat_data0331, assay = "Spatial", slot = "count")
expression_matrix0731 <- GetAssayData(stadseurat_data0731, assay = "Spatial", slot = "count")
expression_matrix0733 <- GetAssayData(stadseurat_data0733, assay = "Spatial", slot = "count")
expression_matrix0734 <- GetAssayData(stadseurat_data0734, assay = "Spatial", slot = "count")
expression_matrix1251 <- GetAssayData(stadseurat_data1251, assay = "Spatial", slot = "count")
expression_matrix1252 <- GetAssayData(stadseurat_data1252, assay = "Spatial", slot = "count")
expression_matrix1253 <- GetAssayData(stadseurat_data1253, assay = "Spatial", slot = "count")
expression_matrix1254 <- GetAssayData(stadseurat_data1254, assay = "Spatial", slot = "count")
expression_matrix1675 <- GetAssayData(stadseurat_data1675, assay = "Spatial", slot = "count")
expression_matrix0331[1:3,1:3]

getwd()
save(expression_matrix0331,expression_matrix0731,expression_matrix0733,expression_matrix0734,
     expression_matrix1251,expression_matrix1252,expression_matrix1253,expression_matrix1254,
     expression_matrix1675,file = "stadexp.RData")

dim(expression_matrix)
deal <- as.data.frame(expression_matrix)
deal[1:3,1:3]



ES <- gsva(expr=deal, 
           gset.idx.list=list,
           kcdf=c("Gaussian"),
           parallel.sz=32)

load("D:/细胞状态2/论文/空间/stad/ES.RData")

EST0331 <- as.data.frame(t(ES0331))
EST0731 <- as.data.frame(t(ES0731))
EST0733 <- as.data.frame(t(ES0733))
EST0734 <- as.data.frame(t(ES0734))
EST1251 <- as.data.frame(t(ES1251))
EST1252 <- as.data.frame(t(ES1252))
EST1253 <- as.data.frame(t(ES1253))
EST1254 <- as.data.frame(t(ES1254))
EST1675 <- as.data.frame(t(ES1675))


stadseurat_data <- stadseurat_data0331
stadseurat_data <- stadseurat_data0731
stadseurat_data <- stadseurat_data0733
stadseurat_data <- stadseurat_data0734
stadseurat_data <- stadseurat_data1251
stadseurat_data <- stadseurat_data1252
stadseurat_data <- stadseurat_data1253
stadseurat_data <- stadseurat_data1254
stadseurat_data <- stadseurat_data1675

class(stadseurat_data)
head(stadseurat_data@meta.data)
stadseurat_data <- NormalizeData(stadseurat_data, normalization.method = "LogNormalize", scale.factor = 10000)
#识别高变基因
stadseurat_data <- FindVariableFeatures(stadseurat_data, selection.method = "vst", nfeatures = 2000)
#找出排在前10个的高变基因
top10 <- head(VariableFeatures(stadseurat_data), 10)
top10
#数据归一化
all.genes <- rownames(stadseurat_data)
stadseurat_data <- ScaleData(stadseurat_data, features = all.genes)
#线性降维PCA主成分分析
stadseurat_data <- RunPCA(stadseurat_data, features = VariableFeatures(object = stadseurat_data))
ElbowPlot(stadseurat_data)
#seurat_data
#--用KNN方法对细胞聚类----------#
stadseurat_data <- FindNeighbors(stadseurat_data, dims = 1:20)
stadseurat_data <- FindClusters(stadseurat_data, resolution =0.6)
head(Idents(stadseurat_data), 5)
#非线性降维(UMAP/tSNE)
stadseurat_data <- RunUMAP(stadseurat_data, dims = 1:20)
DimPlot(stadseurat_data, reduction = "umap",label = TRUE, pt.size = 1)

position=read.csv("./20_00331_LI_SING/spatial/tissue_positions_list.csv",header=F,row.names=1)
position=read.csv("./21_00731_LI_SING/spatial/tissue_positions_list.csv",header=F,row.names=1)
position=read.csv("./21_00733_LI_SING/spatial/tissue_positions_list.csv",header=F,row.names=1)
position=read.csv("./21_00734_LI_SING/spatial/tissue_positions_list.csv",header=F,row.names=1)
position=read.csv("./21_01251_LI_SING/spatial/tissue_positions_list.csv",header=F,row.names=1)
position=read.csv("./21_01252_LI_SING/spatial/tissue_positions_list.csv",header=F,row.names=1)
position=read.csv("./21_01253_LI_SING/spatial/tissue_positions_list.csv",header=F,row.names=1)
position=read.csv("./21_01254_LI_SING/spatial/tissue_positions_list.csv",header=F,row.names=1)
position=read.csv("./21_01675_LI_SING/spatial/tissue_positions_list.csv",header=F,row.names=1)
head(position)
class(position)
colnames(position)

position <- position[ ,c("V5", "V6")]
colnames(position)=paste0("Spatial_",1:ncol(position))
cells_select <- intersect(rownames(position), colnames(stadseurat_data))
dim(stadseurat_data)
dim(position)
position2 <- position[cells_select,]
stadseurat_data[["spatial"]]<-CreateDimReducObject(embeddings=as.matrix(position2),
                                                   key="Spatial_",assay="RNA")
DimPlot(stadseurat_data,reduction="spatial",pt.size=1)

SpatialFeaturePlot(stadseurat_data, features = "nCount_Spatial") + theme(legend.position = "right")

head(stadseurat_data@meta.data)
SpatialDimPlot(stadseurat_data,alpha=1)







BiocManager::install("glmGamPoi")
install.packages("future.apply")
library(future)
library(glmGamPoi)
library(dplyr)


genes <- c("EGFR", "KRAS", "TP53", "CDKN2A","CD44","CDH1")    #cancer

pdf("/wangzhenzenlab/fsh/xff/stad/Spatial/Epithelial.pdf")
custom_colors <- c("gray90", "#c11e21")
SpatialFeaturePlot(stadseurat_data, features = genes,image.alpha = 0,alpha = c(0.5,1),pt.size.factor = 3)&
  scale_fill_gradientn(
    colours = custom_colors, # 使用自定义颜色
    na.value = "gray90"      # 修改图例标题
  ) 
#dev.off()

SpatialFeaturePlot(stadseurat_data, features = c("CXCL6","CXCL8","MMP3"),image.alpha = 0,alpha = c(0.5,1),pt.size.factor = 3)&
  scale_fill_gradientn(
    colours = custom_colors, # 使用自定义颜色
    na.value = "gray90"      # 修改图例标题
  ) 
SpatialFeaturePlot(stadseurat_data, features = c("FCRLA","CR2","FDCSP","CXCL13","MS4A1","TCL1A","VPREB3"),image.alpha = 0,alpha = c(0.5,1),pt.size.factor = 3)&
  scale_fill_gradientn(
    colours = custom_colors, # 使用自定义颜色
    na.value = "gray90"      # 修改图例标题
  ) 
SpatialFeaturePlot(stadseurat_data, features = c("PRSS21","DPEP3","COL9A3","C20orf204"),image.alpha = 0,alpha = c(0.5,1),pt.size.factor = 3)&
  scale_fill_gradientn(
    colours = custom_colors, # 使用自定义颜色
    na.value = "gray90"      # 修改图例标题
  ) 



0:B,Myeloid
1:mast
2:Epithelial
3:Epithelial
4:Epithelial
5:Fibroblasts
6:Endothelial
7:Epithelial
8:Epithelial
9:T
10:Fibroblasts
11:B,T




# load("D:/细胞状态2/论文/空间/stad/stadseurat_data.RData")
# load("D:/细胞状态2/论文/空间/stad/stadEST.RData")
# load("/wangzhenzenlab/fsh/xff/stad/stadseurat_data.RData")
# load("/wangzhenzenlab/fsh/xff/stad/stadEST.RData")
# head(stadEST)
# stadEST <- as.data.frame(stadEST)

EST0331$cell <- rownames(EST0331)
EST0731$cell <- rownames(EST0731)
EST0733$cell <- rownames(EST0733)
EST0734$cell <- rownames(EST0734)
EST1251$cell <- rownames(EST1251)
EST1252$cell <- rownames(EST1252)
EST1253$cell <- rownames(EST1253)
EST1254$cell <- rownames(EST1254)
EST1675$cell <- rownames(EST1675)

stmetadata0331 <- stadseurat_data0331@meta.data
stmetadata0731 <- stadseurat_data0731@meta.data
stmetadata0733 <- stadseurat_data0733@meta.data
stmetadata0734 <- stadseurat_data0734@meta.data
stmetadata1251 <- stadseurat_data1251@meta.data
stmetadata1252 <- stadseurat_data1252@meta.data
stmetadata1253 <- stadseurat_data1253@meta.data
stmetadata1254 <- stadseurat_data1254@meta.data
stmetadata1675 <- stadseurat_data1675@meta.data


stmetadata0331$cell <- rownames(stmetadata0331)
stmetadata0731$cell <- rownames(stmetadata0731)
stmetadata0733$cell <- rownames(stmetadata0733)
stmetadata0734$cell <- rownames(stmetadata0734)
stmetadata1251$cell <- rownames(stmetadata1251)
stmetadata1252$cell <- rownames(stmetadata1252)
stmetadata1253$cell <- rownames(stmetadata1253)
stmetadata1254$cell <- rownames(stmetadata1254)
stmetadata1675$cell <- rownames(stmetadata1675)

stmetadata0331 <- merge(stmetadata0331, EST0331, by = "cell")
stmetadata0731 <- merge(stmetadata0731, EST0731, by = "cell")
stmetadata0733 <- merge(stmetadata0733, EST0733, by = "cell")
stmetadata0734 <- merge(stmetadata0734, EST0734, by = "cell")
stmetadata1251 <- merge(stmetadata1251, EST1251, by = "cell")
stmetadata1252 <- merge(stmetadata1252, EST1252, by = "cell")
stmetadata1253 <- merge(stmetadata1253, EST1253, by = "cell")
stmetadata1254 <- merge(stmetadata1254, EST1254, by = "cell")
stmetadata1675 <- merge(stmetadata1675, EST1675, by = "cell")



stmetadata0331$model1 <- rowMeans(stmetadata0331[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
stmetadata0331$group1 <- ifelse(stmetadata0331$model1 > median(stmetadata0331$model1), "high", "low")
rownames(stmetadata0331) <- stmetadata0331$cell

stmetadata0731$model1 <- rowMeans(stmetadata0731[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
stmetadata0731$group1 <- ifelse(stmetadata0731$model1 > median(stmetadata0731$model1), "high", "low")
rownames(stmetadata0731) <- stmetadata0731$cell

stmetadata0733$model1 <- rowMeans(stmetadata0733[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
stmetadata0733$group1 <- ifelse(stmetadata0733$model1 > median(stmetadata0733$model1), "high", "low")
rownames(stmetadata0733) <- stmetadata0733$cell

stmetadata0734$model1 <- rowMeans(stmetadata0734[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
stmetadata0734$group1 <- ifelse(stmetadata0734$model1 > median(stmetadata0734$model1), "high", "low")
rownames(stmetadata0734) <- stmetadata0734$cell

stmetadata1251$model1 <- rowMeans(stmetadata1251[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
stmetadata1251$group1 <- ifelse(stmetadata1251$model1 > median(stmetadata1251$model1), "high", "low")
rownames(stmetadata1251) <- stmetadata1251$cell

stmetadata1252$model1 <- rowMeans(stmetadata1252[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
stmetadata1252$group1 <- ifelse(stmetadata1252$model1 > median(stmetadata1252$model1), "high", "low")
rownames(stmetadata1252) <- stmetadata1252$cell

stmetadata1253$model1 <- rowMeans(stmetadata1253[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
stmetadata1253$group1 <- ifelse(stmetadata1253$model1 > median(stmetadata1253$model1), "high", "low")
rownames(stmetadata1253) <- stmetadata1253$cell

stmetadata1254$model1 <- rowMeans(stmetadata1254[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
stmetadata1254$group1 <- ifelse(stmetadata1254$model1 > median(stmetadata1254$model1), "high", "low")
rownames(stmetadata1254) <- stmetadata1254$cell

stmetadata1675$model1 <- rowMeans(stmetadata1675[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
stmetadata1675$group1 <- ifelse(stmetadata1675$model1 > median(stmetadata1675$model1), "high", "low")
rownames(stmetadata1675) <- stmetadata1675$cell

save(stmetadata0331,stmetadata0731,stmetadata0733,stmetadata0734,
     stmetadata1251,stmetadata1252,stmetadata1253,stmetadata1254,
     stmetadata1675,file = "stmetadata.RData")

stadseurat_data0331@meta.data <- stmetadata0331
stadseurat_data0731@meta.data <- stmetadata0731
stadseurat_data0733@meta.data <- stmetadata0733
stadseurat_data0734@meta.data <- stmetadata0734
stadseurat_data1251@meta.data <- stmetadata1251
stadseurat_data1252@meta.data <- stmetadata1252
stadseurat_data1253@meta.data <- stmetadata1253
stadseurat_data1254@meta.data <- stmetadata1254
stadseurat_data1675@meta.data <- stmetadata1675


save(stadseurat_data0331,stadseurat_data0731,stadseurat_data0733,stadseurat_data0734,
     stadseurat_data1251,stadseurat_data1252,stadseurat_data1253,stadseurat_data1254,
     stadseurat_data1675,file = "stadseurat_data.RData")

SpatialDimPlot(stadseurat_data0331, alpha = 1, group.by = "group1", cols = c("high" = "#c11e21", "low" = "steelblue"), image.alpha = 0)
SpatialDimPlot(stadseurat_data0731, alpha = 1, group.by = "group1", cols = c("high" = "#c11e21", "low" = "steelblue"), image.alpha = 0)
SpatialDimPlot(stadseurat_data0733, alpha = 1, group.by = "group1", cols = c("high" = "#c11e21", "low" = "steelblue"), image.alpha = 0)
SpatialDimPlot(stadseurat_data0734, alpha = 1, group.by = "group1", cols = c("high" = "#c11e21", "low" = "steelblue"), image.alpha = 0)
SpatialDimPlot(stadseurat_data1251, alpha = 1, group.by = "group1", cols = c("high" = "#c11e21", "low" = "steelblue"), image.alpha = 0)
SpatialDimPlot(stadseurat_data1252, alpha = 1, group.by = "group1", cols = c("high" = "#c11e21", "low" = "steelblue"), image.alpha = 0)
SpatialDimPlot(stadseurat_data1253, alpha = 1, group.by = "group1", cols = c("high" = "#c11e21", "low" = "steelblue"), image.alpha = 0)
SpatialDimPlot(stadseurat_data1254, alpha = 1, group.by = "group1", cols = c("high" = "#c11e21", "low" = "steelblue"), image.alpha = 0)
SpatialDimPlot(stadseurat_data1675, alpha = 1, group.by = "group1", cols = c("high" = "#c11e21", "low" = "steelblue"), image.alpha = 0)

SpatialDimPlot(stadseurat_data, alpha = 1, group.by = "group1", cols = c("high" = "#c11e21", "low" = "steelblue"), image.alpha = 0)

head(stadseurat_data0331)
SpatialDimPlot(stadseurat_data0331,alpha=1,group.by = "group1")+ stat_contour(aes(z = stadseurat_data0331$model1), color = "black", size = 1)

SpatialFeaturePlot(stadseurat_data0331, features = c("model1"), image.alpha = 0)
SpatialFeaturePlot(stadseurat_data0731, features = c("model1"), image.alpha = 0)
SpatialFeaturePlot(stadseurat_data0733, features = c("model1"), image.alpha = 0)
SpatialFeaturePlot(stadseurat_data0734, features = c("model1"), image.alpha = 0)
SpatialFeaturePlot(stadseurat_data1251, features = c("model1"), image.alpha = 0)
SpatialFeaturePlot(stadseurat_data1252, features = c("model1"), image.alpha = 0)
SpatialFeaturePlot(stadseurat_data1253, features = c("model1"), image.alpha = 0)
SpatialFeaturePlot(stadseurat_data1254, features = c("model1"), image.alpha = 0)
SpatialFeaturePlot(stadseurat_data1675, features = c("model1"), image.alpha = 0)


SpatialFeaturePlot(stadseurat_data, features = c("model1"), image.alpha = 0)










dealT <- t(deal)
dealT2 <- dealT[,c("SELE","CXCL12","EGR1","PDPN","CD44","COL6A2","LGALS1","FBLN5")]
head(dealT2)
write.csv(dealT2,file = "D:/细胞状态2/论文/空间/stad/dealT2.csv")

save(stadseurat_data,file = "D:/细胞状态2/论文/空间/stad/stadseurat_data.RData")
load("/wangzhenzenlab/fsh/xff/stad/stadseurat_data.RData")

pdf("/wangzhenzenlab/fsh/xff/stad/stadcelltype.pdf")
SpatialDimPlot(stadseurat_data,group.by = "celltype",alpha=1)
dev.off()














#PAAD
#setwd("D:/细胞状态2/论文/空间/paad/nfs")
setwd("/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/paad/")
paadseurat_data1905=Load10X_Spatial(data.dir="./SS1905133_processed",
                                    filename="filtered_feature_bc_matrix.h5",
                                    assay="Spatial",
                                    slice="test")

paadseurat_data1923=Load10X_Spatial(data.dir="./SS1923404_processed",
                                    filename="filtered_feature_bc_matrix.h5",
                                    assay="Spatial",
                                    slice="test")

paadseurat_data1945=Load10X_Spatial(data.dir="./SS1945070_processed",
                                    filename="filtered_feature_bc_matrix.h5",
                                    assay="Spatial",
                                    slice="test")

paadseurat_data1960=Load10X_Spatial(data.dir="./SS1960050_processed",
                                    filename="filtered_feature_bc_matrix.h5",
                                    assay="Spatial",
                                    slice="test")

paadseurat_data2002=Load10X_Spatial(data.dir="./SS2002154_processed",
                                    filename="filtered_feature_bc_matrix.h5",
                                    assay="Spatial",
                                    slice="test")

paadseurat_data2005=Load10X_Spatial(data.dir="./SS2005094_processed",
                                    filename="filtered_feature_bc_matrix.h5",
                                    assay="Spatial",
                                    slice="test")

paadseurat_data2021=Load10X_Spatial(data.dir="./SS2021309_processed",
                                    filename="filtered_feature_bc_matrix.h5",
                                    assay="Spatial",
                                    slice="test")



#计算细胞状态的gsva
library(GSVA)
gene_sets <- read.csv("D:/细胞状态2/cell_state.csv",header = T)
head(gene_sets)
list=split(gene_sets$GeneName,gene_sets$Signature )

expression_matrix1905 <- GetAssayData(paadseurat_data1905, assay = "Spatial", slot = "count")
expression_matrix1923 <- GetAssayData(paadseurat_data1923, assay = "Spatial", slot = "count")
expression_matrix1945 <- GetAssayData(paadseurat_data1945, assay = "Spatial", slot = "count")
expression_matrix1960 <- GetAssayData(paadseurat_data1960, assay = "Spatial", slot = "count")
expression_matrix2002 <- GetAssayData(paadseurat_data2002, assay = "Spatial", slot = "count")
expression_matrix2005 <- GetAssayData(paadseurat_data2005, assay = "Spatial", slot = "count")
expression_matrix2021 <- GetAssayData(paadseurat_data2021, assay = "Spatial", slot = "count")

expression_matrix1905[1:5,1:5]

getwd()
save(expression_matrix1905,expression_matrix1923,expression_matrix1945,expression_matrix1960,
     expression_matrix2002,expression_matrix2005,expression_matrix2021,file = "paadexp.RData")

dim(expression_matrix)
deal <- as.data.frame(expression_matrix)
deal[1:3,1:3]



ES <- gsva(expr=deal, 
           gset.idx.list=list,
           kcdf=c("Gaussian"),
           parallel.sz=32)

#load("D:/细胞状态2/论文/挖掘/独立验证集/PAAD/paadES.RData")
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/挖掘/特征/paad/paadES.RData")
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/挖掘/特征/paad/ES1905.RData")


EST1905 <- as.data.frame(t(ES1905))
EST1923 <- as.data.frame(t(ES1923)) 
EST1945 <- as.data.frame(t(ES1945))
EST1960 <- as.data.frame(t(ES1960))
EST2002 <- as.data.frame(t(ES2002))
EST2005 <- as.data.frame(t(ES2005))
EST2021 <- as.data.frame(t(ES2021))

stadseurat_data <- paadseurat_data1905
stadseurat_data <- paadseurat_data1923
stadseurat_data <- paadseurat_data1945
stadseurat_data <- paadseurat_data1960
stadseurat_data <- paadseurat_data2002
stadseurat_data <- paadseurat_data2005
stadseurat_data <- paadseurat_data2021



position=read.csv("./SS1905133_processed/spatial/tissue_positions_list.csv",header=F,row.names=1)
position=read.csv("./SS1923404_processed/spatial/tissue_positions_list.csv",header=F,row.names=1)
position=read.csv("./SS1945070_processed/spatial/tissue_positions_list.csv",header=F,row.names=1)
position=read.csv("./SS1960050_processed/spatial/tissue_positions_list.csv",header=F,row.names=1)
position=read.csv("./SS2002154_processed/spatial/tissue_positions_list.csv",header=F,row.names=1)
position=read.csv("./SS2005094_processed/spatial/tissue_positions_list.csv",header=F,row.names=1)
position=read.csv("./SS2021309_processed/spatial/tissue_positions_list.csv",header=F,row.names=1)


position <- position[ ,c("V5", "V6")]
colnames(position)=paste0("Spatial_",1:ncol(position))
cells_select <- intersect(rownames(position), colnames(stadseurat_data))
dim(stadseurat_data)
dim(position)
position2 <- position[cells_select,]
stadseurat_data[["spatial"]]<-CreateDimReducObject(embeddings=as.matrix(position2),
                                                   key="Spatial_",assay="RNA")
DimPlot(stadseurat_data,reduction="spatial",pt.size=1)
SpatialFeaturePlot(stadseurat_data, features = "nCount_Spatial") + theme(legend.position = "right")
head(stadseurat_data@meta.data)
SpatialDimPlot(stadseurat_data,alpha=1)







EST1905$cell <- rownames(EST1905)
EST1923$cell <- rownames(EST1923)
EST1945$cell <- rownames(EST1945)
EST1960$cell <- rownames(EST1960)
EST2002$cell <- rownames(EST2002)
EST2005$cell <- rownames(EST2005)
EST2021$cell <- rownames(EST2021)

paadstmetadata1905 <- paadseurat_data1905@meta.data
paadstmetadata1923 <- paadseurat_data1923@meta.data
paadstmetadata1945 <- paadseurat_data1945@meta.data
paadstmetadata1960 <- paadseurat_data1960@meta.data
paadstmetadata2002 <- paadseurat_data2002@meta.data
paadstmetadata2005 <- paadseurat_data2005@meta.data
paadstmetadata2021 <- paadseurat_data2021@meta.data


paadstmetadata1905$cell <- rownames(paadstmetadata1905)
paadstmetadata1923$cell <- rownames(paadstmetadata1923)
paadstmetadata1945$cell <- rownames(paadstmetadata1945)
paadstmetadata1960$cell <- rownames(paadstmetadata1960)
paadstmetadata2002$cell <- rownames(paadstmetadata2002)
paadstmetadata2005$cell <- rownames(paadstmetadata2005)
paadstmetadata2021$cell <- rownames(paadstmetadata2021)

paadstmetadata1905 <- merge(paadstmetadata1905, EST1905, by = "cell")
paadstmetadata1923 <- merge(paadstmetadata1923, EST1923, by = "cell")
paadstmetadata1945 <- merge(paadstmetadata1945, EST1945, by = "cell")
paadstmetadata1960 <- merge(paadstmetadata1960, EST1960, by = "cell")
paadstmetadata2002 <- merge(paadstmetadata2002, EST2002, by = "cell")
paadstmetadata2005 <- merge(paadstmetadata2005, EST2005, by = "cell")
paadstmetadata2021 <- merge(paadstmetadata2021, EST2021, by = "cell")



paadstmetadata1905$model1 <- rowMeans(paadstmetadata1905[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
paadstmetadata1905$group1 <- ifelse(paadstmetadata1905$model1 > median(paadstmetadata1905$model1), "high", "low")
rownames(paadstmetadata1905) <- paadstmetadata1905$cell

paadstmetadata1923$model1 <- rowMeans(paadstmetadata1923[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
paadstmetadata1923$group1 <- ifelse(paadstmetadata1923$model1 > median(paadstmetadata1923$model1), "high", "low")
rownames(paadstmetadata1923) <- paadstmetadata1923$cell

paadstmetadata1945$model1 <- rowMeans(paadstmetadata1945[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
paadstmetadata1945$group1 <- ifelse(paadstmetadata1945$model1 > median(paadstmetadata1945$model1), "high", "low")
rownames(paadstmetadata1945) <- paadstmetadata1945$cell

paadstmetadata1960$model1 <- rowMeans(paadstmetadata1960[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
paadstmetadata1960$group1 <- ifelse(paadstmetadata1960$model1 > median(paadstmetadata1960$model1), "high", "low")
rownames(paadstmetadata1960) <- paadstmetadata1960$cell

paadstmetadata2002$model1 <- rowMeans(paadstmetadata2002[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
paadstmetadata2002$group1 <- ifelse(paadstmetadata2002$model1 > median(paadstmetadata2002$model1), "high", "low")
rownames(paadstmetadata2002) <- paadstmetadata2002$cell

paadstmetadata2005$model1 <- rowMeans(paadstmetadata2005[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
paadstmetadata2005$group1 <- ifelse(paadstmetadata2005$model1 > median(paadstmetadata2005$model1), "high", "low")
rownames(paadstmetadata2005) <- paadstmetadata2005$cell

paadstmetadata2021$model1 <- rowMeans(paadstmetadata2021[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
paadstmetadata2021$group1 <- ifelse(paadstmetadata2021$model1 > median(paadstmetadata2021$model1), "high", "low")
rownames(paadstmetadata2021) <- paadstmetadata2021$cell

save(paadstmetadata1905, paadstmetadata1923, paadstmetadata1945, paadstmetadata1960, 
     paadstmetadata2002, paadstmetadata2005, paadstmetadata2021, 
     file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/paad/paadstmetadata.RData")






#BRCA
setwd("/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/BRCA")

brcaseurat_data7599=Load10X_Spatial(data.dir="./GSM6177599",
                                    filename="GSM6177599.h5",
                                    assay="Spatial",
                                    slice="test")

brcaseurat_data7601=Load10X_Spatial(data.dir="./GSM6177601",
                                    filename="GSM6177601.h5",
                                    assay="Spatial",
                                    slice="test")

brcaseurat_data7603=Load10X_Spatial(data.dir="./GSM6177603",
                                    filename="GSM6177603.h5",
                                    assay="Spatial",
                                    slice="test")


expression_matrix7599 <- GetAssayData(brcaseurat_data7599, assay = "Spatial", slot = "count")
expression_matrix7601 <- GetAssayData(brcaseurat_data7601, assay = "Spatial", slot = "count")
expression_matrix7603 <- GetAssayData(brcaseurat_data7603, assay = "Spatial", slot = "count")

save(expression_matrix7599,expression_matrix7601,expression_matrix7603,file = "BRCAexp.RData")



ES7599 <- gsva(expr=expression_matrix7599, 
           gset.idx.list=list,
           kcdf=c("Gaussian"),
           parallel.sz=32)

ES7601 <- gsva(expr=expression_matrix7601, 
               gset.idx.list=list,
               kcdf=c("Gaussian"),
               parallel.sz=32)

ES7603 <- gsva(expr=expression_matrix7603, 
               gset.idx.list=list,
               kcdf=c("Gaussian"),
               parallel.sz=32)

save(ES7599,ES7601,ES7603,file = "/wangzhenzenlab/shenjie/xff/st/BRCAES.RData")
load("/wangzhenzenlab/shenjie/xff/st/BRCAES.RData")

EST7599 <- as.data.frame(t(ES7599))
EST7601 <- as.data.frame(t(ES7601)) 
EST7603 <- as.data.frame(t(ES7603))


position7599=read.csv("./GSM6177599/spatial/tissue_positions_list.csv",header=F,row.names=1)
position7601=read.csv("./GSM6177601/spatial/tissue_positions_list.csv",header=F,row.names=1)
position7603=read.csv("./GSM6177603/spatial/tissue_positions_list.csv",header=F,row.names=1)

position7599 <- position7599[ ,c("V5", "V6")]
colnames(position7599)=paste0("Spatial_",1:ncol(position7599))
cells_select <- intersect(rownames(position7599), colnames(brcaseurat_data7599))
dim(brcaseurat_data7599)
dim(position7599)
position7599_2 <- position7599[cells_select,]
brcaseurat_data7599[["spatial"]]<-CreateDimReducObject(embeddings=as.matrix(position7599_2),
                                                       key="Spatial_",assay="RNA")

position7601 <- position7601[ ,c("V5", "V6")]
colnames(position7601)=paste0("Spatial_",1:ncol(position7601))
cells_select <- intersect(rownames(position7601), colnames(brcaseurat_data7601))
dim(brcaseurat_data7601)
dim(position7601)
position7601_2 <- position7601[cells_select,]
brcaseurat_data7601[["spatial"]]<-CreateDimReducObject(embeddings=as.matrix(position7601_2),
                                                       key="Spatial_",assay="RNA")

position7603 <- position7603[ ,c("V5", "V6")]
colnames(position7603)=paste0("Spatial_",1:ncol(position7603))
cells_select <- intersect(rownames(position7603), colnames(brcaseurat_data7603))
dim(brcaseurat_data7603)
dim(position7603)
position7603_2 <- position7603[cells_select,]
brcaseurat_data7603[["spatial"]]<-CreateDimReducObject(embeddings=as.matrix(position7603_2),
                                                       key="Spatial_",assay="RNA")

EST7599$cell <- rownames(EST7599)
EST7601$cell <- rownames(EST7601)
EST7603$cell <- rownames(EST7603)

brcastmetadata7599 <- brcaseurat_data7599@meta.data
brcastmetadata7601 <- brcaseurat_data7601@meta.data
brcastmetadata7603 <- brcaseurat_data7603@meta.data

brcastmetadata7599$cell <- rownames(brcastmetadata7599)
brcastmetadata7601$cell <- rownames(brcastmetadata7601)
brcastmetadata7603$cell <- rownames(brcastmetadata7603)


brcastmetadata7599 <- merge(brcastmetadata7599, EST7599, by = "cell")
brcastmetadata7601 <- merge(brcastmetadata7601, EST7601, by = "cell")
brcastmetadata7603 <- merge(brcastmetadata7603, EST7603, by = "cell")


brcastmetadata7599$model1 <- rowMeans(brcastmetadata7599[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
brcastmetadata7599$group1 <- ifelse(brcastmetadata7599$model1 > median(brcastmetadata7599$model1), "high", "low")
rownames(brcastmetadata7599) <- brcastmetadata7599$cell

brcastmetadata7601$model1 <- rowMeans(brcastmetadata7601[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
brcastmetadata7601$group1 <- ifelse(brcastmetadata7601$model1 > median(brcastmetadata7601$model1), "high", "low")
rownames(brcastmetadata7601) <- brcastmetadata7601$cell

brcastmetadata7603$model1 <- rowMeans(brcastmetadata7603[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
brcastmetadata7603$group1 <- ifelse(brcastmetadata7603$model1 > median(brcastmetadata7603$model1), "high", "low")
rownames(brcastmetadata7603) <- brcastmetadata7603$cell

save(brcastmetadata7599, brcastmetadata7601, brcastmetadata7603,
     file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/BRCA/brcastmetadata.RData")


#KIRC
setwd("/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/KIRC")

kircseurat_data4030=Load10X_Spatial(data.dir="./GSM5924030",
                                    filename="GSM5924030.h5",
                                    assay="Spatial",
                                    slice="test")

kircseurat_data4031=Load10X_Spatial(data.dir="./GSM5924031",
                                    filename="GSM5924031.h5",
                                    assay="Spatial",
                                    slice="test")

kircseurat_data4032=Load10X_Spatial(data.dir="./GSM5924032",
                                    filename="GSM5924032.h5",
                                    assay="Spatial",
                                    slice="test")


expression_matrix4030 <- GetAssayData(kircseurat_data4030, assay = "Spatial", slot = "count")
expression_matrix4031 <- GetAssayData(kircseurat_data4031, assay = "Spatial", slot = "count")
expression_matrix4032 <- GetAssayData(kircseurat_data4032, assay = "Spatial", slot = "count")

save(expression_matrix4030,expression_matrix4031,expression_matrix4032,file = "KIRCexp.RData")



ES4030 <- gsva(expr=expression_matrix4030, 
               gset.idx.list=list,
               kcdf=c("Gaussian"),
               parallel.sz=32)

ES4031 <- gsva(expr=expression_matrix4031, 
               gset.idx.list=list,
               kcdf=c("Gaussian"),
               parallel.sz=32)

ES4032 <- gsva(expr=expression_matrix4032, 
               gset.idx.list=list,
               kcdf=c("Gaussian"),
               parallel.sz=32)

save(ES4030,ES4031,ES4032,file = "/wangzhenzenlab/shenjie/xff/st/KIRCES.RData")

load("/wangzhenzenlab/shenjie/xff/st/KIRCES.RData")

EST4030 <- as.data.frame(t(ES4030))
EST4031 <- as.data.frame(t(ES4031)) 
EST4032 <- as.data.frame(t(ES4032))


position4030=read.csv("./GSM5924030/spatial/tissue_positions_list.csv",header=F,row.names=1)
position4031=read.csv("./GSM5924031/spatial/tissue_positions_list.csv",header=F,row.names=1)
position4032=read.csv("./GSM5924032/spatial/tissue_positions_list.csv",header=F,row.names=1)

position4030 <- position4030[ ,c("V5", "V6")]
colnames(position4030)=paste0("Spatial_",1:ncol(position4030))
cells_select <- intersect(rownames(position4030), colnames(kircseurat_data4030))
dim(kircseurat_data4030)
dim(position4030)
position4030_2 <- position4030[cells_select,]
kircseurat_data4030[["spatial"]]<-CreateDimReducObject(embeddings=as.matrix(position4030_2),
                                                       key="Spatial_",assay="RNA")

position4031 <- position4031[ ,c("V5", "V6")]
colnames(position4031)=paste0("Spatial_",1:ncol(position4031))
cells_select <- intersect(rownames(position4031), colnames(kircseurat_data4031))
dim(kircseurat_data4031)
dim(position4031)
position4031_2 <- position4031[cells_select,]
kircseurat_data4031[["spatial"]]<-CreateDimReducObject(embeddings=as.matrix(position4031_2),
                                                       key="Spatial_",assay="RNA")

position4032 <- position4032[ ,c("V5", "V6")]
colnames(position4032)=paste0("Spatial_",1:ncol(position4032))
cells_select <- intersect(rownames(position4032), colnames(kircseurat_data4032))
dim(kircseurat_data4032)
dim(position4032)
position4032_2 <- position4032[cells_select,]
kircseurat_data4032[["spatial"]]<-CreateDimReducObject(embeddings=as.matrix(position4032_2),
                                                       key="Spatial_",assay="RNA")

EST4030$cell <- rownames(EST4030)
EST4031$cell <- rownames(EST4031)
EST4032$cell <- rownames(EST4032)

kircstmetadata4030 <- kircseurat_data4030@meta.data
kircstmetadata4031 <- kircseurat_data4031@meta.data
kircstmetadata4032 <- kircseurat_data4032@meta.data

kircstmetadata4030$cell <- rownames(kircstmetadata4030)
kircstmetadata4031$cell <- rownames(kircstmetadata4031)
kircstmetadata4032$cell <- rownames(kircstmetadata4032)


kircstmetadata4030 <- merge(kircstmetadata4030, EST4030, by = "cell")
kircstmetadata4031 <- merge(kircstmetadata4031, EST4031, by = "cell")
kircstmetadata4032 <- merge(kircstmetadata4032, EST4032, by = "cell")


kircstmetadata4030$model1 <- rowMeans(kircstmetadata4030[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
kircstmetadata4030$group1 <- ifelse(kircstmetadata4030$model1 > median(kircstmetadata4030$model1), "high", "low")
rownames(kircstmetadata4030) <- kircstmetadata4030$cell

kircstmetadata4031$model1 <- rowMeans(kircstmetadata4031[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
kircstmetadata4031$group1 <- ifelse(kircstmetadata4031$model1 > median(kircstmetadata4031$model1), "high", "low")
rownames(kircstmetadata4031) <- kircstmetadata4031$cell

kircstmetadata4032$model1 <- rowMeans(kircstmetadata4032[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
kircstmetadata4032$group1 <- ifelse(kircstmetadata4032$model1 > median(kircstmetadata4032$model1), "high", "low")
rownames(kircstmetadata4032) <- kircstmetadata4032$cell

save(kircstmetadata4030, kircstmetadata4031, kircstmetadata4032,
     file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/KIRC/kircstmetadata.RData")




#SKCM
setwd("/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/SKCM")

skcmseurat_data7624=Load10X_Spatial(data.dir="./GSM6177624",
                                    filename="GSM6177624.h5",
                                    assay="Spatial",
                                    slice="test")

skcmseurat_data7625=Load10X_Spatial(data.dir="./GSM6177625",
                                    filename="GSM6177625.h5",
                                    assay="Spatial",
                                    slice="test")

expression_matrix7624 <- GetAssayData(skcmseurat_data7624, assay = "Spatial", slot = "count")
expression_matrix7625 <- GetAssayData(skcmseurat_data7625, assay = "Spatial", slot = "count")

rownames(expression_matrix7624) <- sub(".*?-", "", rownames(expression_matrix7624))
rownames(expression_matrix7625) <- sub(".*?-", "", rownames(expression_matrix7625))

save(expression_matrix7624,expression_matrix7625,file = "SKCMexp.RData")



ES7624 <- gsva(expr=expression_matrix7624, 
               gset.idx.list=list,
               kcdf=c("Gaussian"),
               parallel.sz=32)

ES7625 <- gsva(expr=expression_matrix7625, 
               gset.idx.list=list,
               kcdf=c("Gaussian"),
               parallel.sz=32)


save(ES7624,ES7625,file = "/wangzhenzenlab/shenjie/xff/st/SKCMES.RData")

load("/wangzhenzenlab/shenjie/xff/st/SKCMES.RData")

EST7624 <- as.data.frame(t(ES7624))
EST7625 <- as.data.frame(t(ES7625)) 


position7624=read.csv("./GSM6177624/spatial/tissue_positions_list.csv",header=F,row.names=1)
position7625=read.csv("./GSM6177625/spatial/tissue_positions_list.csv",header=F,row.names=1)

position7624 <- position7624[ ,c("V5", "V6")]
colnames(position7624)=paste0("Spatial_",1:ncol(position7624))
cells_select <- intersect(rownames(position7624), colnames(skcmseurat_data7624))
dim(skcmseurat_data7624)
dim(position7624)
position7624_2 <- position7624[cells_select,]
skcmseurat_data7624[["spatial"]]<-CreateDimReducObject(embeddings=as.matrix(position7624_2),
                                                       key="Spatial_",assay="RNA")

position7625 <- position7625[ ,c("V5", "V6")]
colnames(position7625)=paste0("Spatial_",1:ncol(position7625))
cells_select <- intersect(rownames(position7625), colnames(skcmseurat_data7625))
dim(skcmseurat_data7625)
dim(position7625)
position7625_2 <- position7625[cells_select,]
skcmseurat_data7625[["spatial"]]<-CreateDimReducObject(embeddings=as.matrix(position7625_2),
                                                       key="Spatial_",assay="RNA")



EST7624$cell <- rownames(EST7624)
EST7625$cell <- rownames(EST7625)

skcmstmetadata7624 <- skcmseurat_data7624@meta.data
skcmstmetadata7625 <- skcmseurat_data7625@meta.data

skcmstmetadata7624$cell <- rownames(skcmstmetadata7624)
skcmstmetadata7625$cell <- rownames(skcmstmetadata7625)


skcmstmetadata7624 <- merge(skcmstmetadata7624, EST7624, by = "cell")
skcmstmetadata7625 <- merge(skcmstmetadata7625, EST7625, by = "cell")


skcmstmetadata7624$model1 <- rowMeans(skcmstmetadata7624[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
skcmstmetadata7624$group1 <- ifelse(skcmstmetadata7624$model1 > median(skcmstmetadata7624$model1), "high", "low")
rownames(skcmstmetadata7624) <- skcmstmetadata7624$cell

skcmstmetadata7625$model1 <- rowMeans(skcmstmetadata7625[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
skcmstmetadata7625$group1 <- ifelse(skcmstmetadata7625$model1 > median(skcmstmetadata7625$model1), "high", "low")
rownames(skcmstmetadata7625) <- skcmstmetadata7625$cell



save(skcmstmetadata7624, skcmstmetadata7625,
     file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/SKCM/skcmstmetadata.RData")

#分组函数

fenzufun <- function(stmetadata){
  stmetadata2 <- stmetadata%>%
    mutate(level = cut(model1,
                       breaks = quantile(model1, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                       labels = c("low", "medium", "high"),
                       include.lowest = TRUE))
  return(stmetadata2)
}


brcastmetadata7599_2 <- fenzufun(brcastmetadata7599)
brcastmetadata7601_2 <- fenzufun(brcastmetadata7601)
brcastmetadata7603_2 <- fenzufun(brcastmetadata7603)
kircstmetadata4030_2 <- fenzufun(kircstmetadata4030)
kircstmetadata4031_2 <- fenzufun(kircstmetadata4031)
kircstmetadata4032_2 <- fenzufun(kircstmetadata4032)
skcmstmetadata7624_2 <- fenzufun(skcmstmetadata7624)
skcmstmetadata7625_2 <- fenzufun(skcmstmetadata7625)


save(brcastmetadata7599_2, brcastmetadata7601_2, brcastmetadata7603_2,
     file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/BRCA/newbrcametadata.RData")

save(kircstmetadata4030_2, kircstmetadata4031_2, kircstmetadata4032_2,
     file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/KIRC/newkircmetadata.RData")

save(skcmstmetadata7624_2, skcmstmetadata7625_2,
     file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/SKCM/newskcmmetadata.RData")




load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/BRCA/newbrcametadata.RData")
brcaseurat_data7599@meta.data <- brcastmetadata7599_2
brcaseurat_data7601@meta.data <- brcastmetadata7601_2
brcaseurat_data7603@meta.data <- brcastmetadata7603_2
save(brcaseurat_data7599,brcaseurat_data7601,brcaseurat_data7603,
     file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/BRCA/newbrcaseurat_data.RData")

#7599
pdf("/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/BRCA/7599.pdf")
SpatialDimPlot(brcaseurat_data7599, alpha = 1, group.by = "level", 
               cols = c("high" = "#c11e21", "medium"= "#99bca0","low" = "steelblue"),
               image.alpha = 0,pt.size.factor = 2.5)
dev.off()
# 7601
pdf("/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/BRCA/7601.pdf")
SpatialDimPlot(brcaseurat_data7601, alpha = 1, group.by = "level",
               cols = c("high" = "#c11e21", "medium"= "#99bca0","low" = "steelblue"), 
               image.alpha = 0,pt.size.factor = 2.5)
dev.off()
# 7603
pdf("/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/BRCA/7603.pdf")
SpatialDimPlot(brcaseurat_data7603, alpha = 1, group.by = "level", 
               cols = c("high" = "#c11e21", "medium"= "#99bca0","low" = "steelblue"), 
               image.alpha = 0,pt.size.factor = 2.5)
dev.off()



load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/KIRC/newkircmetadata.RData")
kircseurat_data4030@meta.data <- kircstmetadata4030_2
kircseurat_data4031@meta.data <- kircstmetadata4031_2
kircseurat_data4032@meta.data <- kircstmetadata4032_2
save(kircseurat_data4030,kircseurat_data4031,kircseurat_data4032,
     file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/KIRC/newkircseurat_data.RData")
# 4030
pdf("/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/KIRC/4030.pdf")
SpatialDimPlot(kircseurat_data4030, alpha = 1, group.by = "level", 
               cols = c("high" = "#c11e21", "medium"= "#99bca0","low" = "steelblue"), 
               image.alpha = 0,pt.size.factor = 2)
dev.off()

# 4031
pdf("/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/KIRC/4031.pdf")
SpatialDimPlot(kircseurat_data4031, alpha = 1, group.by = "level", 
               cols = c("high" = "#c11e21", "medium"= "#99bca0","low" = "steelblue"), 
               image.alpha = 0,pt.size.factor = 2)
dev.off()

# 4032
pdf("/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/KIRC/4032.pdf")
SpatialDimPlot(kircseurat_data4032, alpha = 1, group.by = "level", 
               cols = c("high" = "#c11e21", "medium"= "#99bca0","low" = "steelblue"), 
               image.alpha = 0,pt.size.factor = 2)
dev.off()



load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/SKCM/newskcmmetadata.RData")
skcmseurat_data7624@meta.data <- skcmstmetadata7624_2
skcmseurat_data7625@meta.data <- skcmstmetadata7625_2
save(skcmseurat_data7624,skcmseurat_data7625,
     file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/SKCM/newskcmseurat_data.RData")
# 7624
pdf("/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/SKCM/7624.pdf")
SpatialDimPlot(skcmseurat_data7624, alpha = 1, group.by = "level", 
               cols = c("high" = "#c11e21", "medium"= "#99bca0","low" = "steelblue"), 
               image.alpha = 0,pt.size.factor = 2.2)
dev.off()

#7625
pdf("/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/SKCM/7625.pdf")
SpatialDimPlot(skcmseurat_data7625, alpha = 1, group.by = "level", 
               cols = c("high" = "#c11e21", "medium"= "#99bca0","low" = "steelblue"), 
               image.alpha = 0,pt.size.factor = 2.7)
dev.off()













