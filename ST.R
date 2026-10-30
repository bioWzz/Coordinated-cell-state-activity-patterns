library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(clustree)
library(patchwork)
#InstallData("stxBrain")
#remotes::install_github('satijalab/seurat-data')
# 读取10x数据，data.dir参数指定存放文件的路径
#setwd("F:/fsh/ST_matrix2/ST/ProstateAcinarCancer_10x")
#setwd("F:/fsh/ST_matrix2/ST/SkinMelanoma_10x/ST_data")
# 创建Seurat对象
#seurat_data <- Read10X(data.dir = "./ST_data")
#读取.h5
fs=list.files(pattern = '.h5')
seurat_data <-Read10X_h5(fs) #读取文件
seurat_data<- CreateSeuratObject(counts = seurat_data, project = "pbmc3k", min.cells = 3, min.features = 200)
#-------#
#QC质量控制来筛选高质量的细胞
seurat_data[["percent.mt"]] <- PercentageFeatureSet(seurat_data, pattern = "^MT-")
#将QC可视化为小提琴图
VlnPlot(seurat_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol =3)
#将特征间关系可视化
FeatureScatter(seurat_data, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(seurat_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#QC指标过滤细胞
seurat_data <- subset(seurat_data, subset = nFeature_RNA > 0 & nFeature_RNA < 11000 & percent.mt < 20)
#数据标准化
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
#确定数据维度
seurat_data <- JackStraw(seurat_data, num.replicate = 100)
seurat_data <- ScoreJackStraw(seurat_data, dims = 1:20)
#对PC打分可视化
JackStrawPlot(seurat_data, dims = 1:15)
#查看resolution
seurat_data <- FindClusters(
  object = seurat_data,
  resolution = c(seq(0,1,.2))
)
clustree(seurat_data@meta.data,prefix='RNA_snn_res.')+ coord_flip()
#--用KNN方法对细胞聚类----------#
seurat_data <- FindNeighbors(seurat_data, dims = 1:10)
seurat_data <- FindClusters(seurat_data, resolution =0.2)
head(Idents(seurat_data), 5)
#非线性降维(UMAP/tSNE)
seurat_data <- RunUMAP(seurat_data, dims = 1:30)
DimPlot(seurat_data, reduction = "umap",label = TRUE, pt.size = 1)
seurat_data <- RunTSNE(seurat_data, dims = 1:30)
DimPlot(seurat_data, reduction = "tsne",label = TRUE, pt.size = 1)
#-----------#
#write.csv(pbmc.markers,"F:/fsh/ST_matrix2/ST/Lung2/ST_data/pbmc.markers.csv")#导出结果
pbmc.markers <- FindAllMarkers(seurat_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
# VlnPlot(pbmc, features = c("IDH1", "IDH2"), slot = "counts", log = TRUE)
# FeaturePlot(pbmc, features = c("IDH1", "IDH2"))
top20 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#DoHeatmap(seurat_data, features = top20$gene) + NoLegend()
FeaturePlot(seurat_data, features = c("IL2RA", "IDO1"), reduction = "tsne")
VlnPlot(seurat_data, features = c("CD8A", "LEF1","TCF7","KLF2","TXK","BACH2","CD8B"))## NAIVE 
VlnPlot(seurat_data, features = c("CD8A", "SLC4A10","ZBTB16","NCR3")) ## MAIT
VlnPlot(seurat_data, features = c("CD8A","CD3D","GNLY","NKG7","GZMB","GZMA","GZMM"))## Effector
VlnPlot(seurat_data, features = c("CD8A","CD3D","GPR183","S100A4","GZMA","GZMK"))     ## Memory 
##基因注释
cellmarker <- read.csv("F:/fsh/ST_matrix2/ST/BRCA/ST_data/cellmarker_h.csv")
cellmarker <- cellmarker[cellmarker$tissue_class=="Lung",]
cellmarker <- merge(cellmarker,top20,by.x="Symbol",by.y ="gene")
cellmarker <- cellmarker[,c(1,3,6,8,10,22,25,26)]
cellmarker <- cellmarker[!duplicated(cellmarker),]

#marker<-cellmarker[cellmarker$cluster==6,]
type0 <-"FIB"#"CAF cell" 
type1 <-"macrophage" #"EPI cell"
type2 <-"CD8T" #"EPI cell" 
type3 <-"EPI" #"EPI cell" 
type4 <-"EPI" #"T cell" 
type5 <-"FlB + macrophage" #"T cell" 
type6 <-"Alveolar epithelial cells" #"END cell" 
type7 <-"EPI" #"EPI cell" 
type8 <-"Alveolar epithelial cells" #"EPI cell" 
type9 <-"SLC16A7+ cell" #"MYE cell"
type10 <-"macrophage"
type11 <-"Ciliated epithelial cells"
type12 <-" Smooth muscle cells"
type13 <-"EPI"
type14 <-"CD8T"
#------------#
type0 <-"T cell"#"CAF cell" 
type1 <-"Other cell" #"EPI cell"
type2 <-"Other cell" #"EPI cell" 
type3 <-"Other cell" #"EPI cell" 
type4 <-"Other cell" #"T cell" 
type5 <-"Other cell" #"T cell" 
type6 <-"Other cell" #"END cell" 
type7 <-"Other cell" #"EPI cell" 
#------------#
# VlnPlot(pbmc, features = c("'SLIT2'"), slot = "counts", log = TRUE)
# FeaturePlot(pbmc, features = c("'SLIT2'"))
# top20$gene[top20$cluster=="7"]
# FeaturePlot(pbmc, features = top20$gene[top20$cluster=="0"])
new.cluster.ids <- c(type0,type1,type2,type3,type4,type5,type6,type7,type8,type9,
                     type10,type11,type12,type13,type14)
new.cluster.ids <- c(type0,type1,type2,type3,type4,type5,type6,type7)
names(new.cluster.ids) <- levels(seurat_data)
seurat_data <- RenameIdents(seurat_data, new.cluster.ids)
#pdf("plot/Treated_umap注释信息.pdf",width = 8,height = 6)
DimPlot(seurat_data, reduction = "tsne",label = TRUE, pt.size = 1)
dev.off()
#----------#
seurat_data@meta.data$cell_type<-seurat_data@active.ident
#--------------#
seurat_data
#------------#
gene_data<-seurat_data[["RNA"]]$counts
meta.data<-seurat_data@meta.data
save(seurat_data,file = "ST_gene_data.RData")

#------#
gene_exp<-seurat_data[["RNA"]]$counts
#gene_exp<-data.frame(gene_exp)
#----------#
mate<-seurat_data@meta.data
T_cell<-mate[mate$cell_type=="T cell",]
gene_exp1<-gene_exp[,T_cell$cell_type]
gene_exp1<-data.frame(gene_exp1)

#------#
save(gene_exp1,file = "gene_exp.RData")
#----GSVA结果-----#
GSVA<- read.csv("F:/fsh/ST_matrix2/ST/Lung2/ST_data/GSVA.csv")
rownames(GSVA)<-GSVA$X
GSVA<-GSVA[,-1]
GSVA<-data.frame(t(GSVA))
#分组
GSVA_fenzu <- GSVA %>%
  mutate(exclude_category = case_when(
    exclude > 0 ~ "高",
    exclude < 0 ~ "低"
  ),
  exhaustion_category = case_when(
    exhaustion > 0 ~ "高",
    exhaustion <0 ~ "低"
  ),
  senility_category = case_when(
    senility > 0 ~ "高",
    senility <0 ~ "低"
  ),
  Tcell.dysfunction_category = case_when(
    Tcell.dysfunction > 0 ~ "高",
    Tcell.dysfunction < 0 ~ "低"
  ))
#分两组
df_grouped <- GSVA_fenzu %>%
  mutate(group = ifelse(exclude_category == "高" & exhaustion_category == "低" & senility_category == "低" & Tcell.dysfunction_category == "低", "Moderate", "Disorderly"))
table(df_grouped$group)
df_grouped$X<-rownames(df_grouped)
df_grouped$X <- gsub("\\.", "-", df_grouped$X)
#-------#
M_cell<-df_grouped[df_grouped$group=="Moderate",]
M_cell<-M_cell[,c(10,9)]
names(M_cell)<-c("cell","group")

seurat_data@meta.data<-seurat_data@meta.data
seurat_data@meta.data$cell<-rownames(seurat_data@meta.data)
seurat_data@meta.data <- seurat_data@meta.data %>%
  left_join(M_cell, by = "cell") %>%  # 将两个数据框按"cell"列进行左连接
  mutate(cell_type = if_else(!is.na(group), paste0(cell_type, "_M"), cell_type)) 

rownames(seurat_data@meta.data)<-seurat_data@meta.data$cell
meta.data<-seurat_data@meta.data
#--------#
# 假设 seurat_data 是一个 Seurat 对象
# 假设我们要检查的列是 cell_type
# 获取 cell_type 不在 active.ident 中的细胞
not_in_active <- seurat_data@meta.data$cell_type[!seurat_data@meta.data$cell_type %in% seurat_data@active.ident]
# 找到这些细胞在 seurat_data 中的索引
rows_to_update <- which(!seurat_data@meta.data$cell_type %in% seurat_data@active.ident)
# 确保 not_in_active 中的所有值都是 active.ident 因子的水平的一部分
# 如果 not_in_active 包含新的水平，我们需要先添加这些水平
new_levels <- setdiff(not_in_active, levels(seurat_data@active.ident))
if (length(new_levels) > 0) {
  levels(seurat_data@active.ident) <- c(levels(seurat_data@active.ident), new_levels)
}
# 更新 active.ident
seurat_data@active.ident[rows_to_update] <- not_in_active
# 检查结果
table(seurat_data@active.ident)
#------------------------#
seurat_data
#------对单细胞数据进行 SCTransform转换-------#
allen_reference <- SCTransform(seurat_data, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
#---可视化注释结果----#
DimPlot(allen_reference, group.by = "cell_type", label = TRUE)
DimPlot(seurat_data, reduction = "umap",label = TRUE, pt.size = 0.8)+ 
  scale_colour_manual(values=c("#c1c1c1","#c1c1c1"  , "#80B1D3" ,"#c1c1c1" ,"#c1c1c1" , "#c1c1c1",
                               "#c1c1c1" ,"#c1c1c1","#c1c1c1"  ,"#d9352a" ))

#----------#
#读入空间信息
test_data=Load10X_Spatial(data.dir="./ST_data",
                            filename="filtered_feature_bc_matrix.h5",
                            assay="Spatial",
                            slice="test")
head(test_data@meta.data)
SpatialDimPlot(test_data,alpha=0)
#-----------#
##获取切片图像坐标：
img<-GetTissueCoordinates(test_data)
img$imagerow=540-img$imagerow
write.csv(img,file="./ST_data/position_information.csv")
head(img)
#读入单细胞的位置信息
position=read.csv("./ST_data/position_information.csv",header=T,row.names=1)
head(position)
position=select(position,imagecol,imagerow)
#将位置信息合并入单细胞seurat对象
colnames(position)=paste0("Spatial_",1:ncol(position))
test_data[["spatial"]]<-CreateDimReducObject(embeddings=as.matrix(position),
                                                key="Spatial",assay="RNA")
DimPlot(test_data,reduction="spatial",pt.size=1)
#---数据的预处理-----#
VlnPlot(test_data, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
SpatialFeaturePlot(test_data, features = "nCount_Spatial") + theme(legend.position = "right")
#-------------#
VlnPlot(test_data,  features = "nFeature_Spatial",pt.size = 0.1,cols = "tomato",group.by = "orig.ident") + NoLegend()
FeatureScatter(test_data, feature1 = "nCount_Spatial",feature2= "nFeature_Spatial",group.by = "orig.ident")+ NoLegend()
#--------------------#
test_data <- SCTransform(test_data, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
#提取表达量变化最高的10个基因
top10<- head(VariableFeatures(test_data),10)
top10
##[1] "Ttr"     "Hbb-bs"  "Plp1"    "Mbp"     "Hba-a1"  "Ptgds"   "Penk"    "Hba-a2"  "S100a5"  "Ppp1r1b"

#类似常规转录组，我们也可以展示高变基因；
VariableFeaturePlot(test_data,cols = c( "gray60", "red"))+NoLegend()
LabelPoints(plot = p5,points = top10, repel = TRUE,xnudge=0,ynudge=0)
#-----基因表达可视化-----------#
SpatialFeaturePlot(test_data, features = c("IGKC", "ZG16"))
SpatialFeaturePlot(test_data, features = "IGKC", pt.size.factor = 1)
SpatialFeaturePlot(test_data, features = "IGKC", alpha = c(0.1, 1))
#--------降维聚类及可视化------------------#
test_data <- RunPCA(test_data, assay = "SCT", verbose = FALSE)
test_data <- FindNeighbors(test_data, reduction = "pca", dims = 1:30)
test_data <- FindClusters(test_data, verbose = FALSE)
test_data <- RunUMAP(test_data, reduction = "pca", dims = 1:30)
#------------#
DimPlot(test_data, reduction = "umap", label = TRUE)
SpatialDimPlot(test_data, label = TRUE, label.size = 3, alpha =0.5)+
  scale_colour_manual(values=c("stromal"="gray60","Cancer"="#d9352a"))
#----#
test_data
test_data@meta.data$cluster<-test_data@active.ident
#meta.data<-test_data@meta.data
#-----计算差异基因做细胞注释--------#
#de_markers <- FindMarkers(test_data, ident.1 = 5, ident.2 = 6)
#head(de_markers)

ST_pbmc.markers <-  FindAllMarkers(test_data, only.pos = TRUE,min.pct = 0.25)
#pbmc.markers1<-pbmc.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top20 <- ST_pbmc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
#-------#
cellmarker <- read.csv("F:/fsh/ST_matrix2/ST/BRCA/ST_data/cellmarker_h.csv")
cellmarker <- cellmarker[cellmarker$tissue_class=="Breast",]
cellmarker <- merge(cellmarker,top20,by.x="Symbol",by.y ="gene")
cellmarker <- cellmarker[,c(1,3,6,8,10,22,25,26)]
cellmarker <- cellmarker[!duplicated(cellmarker),]
type0 <- "stromal"
type1 <- "stromal"
type2 <- "stromal" 
type3 <- "Cancer" 
type4 <- "stromal" 
type5 <- "stromal" 
type6 <- "stromal" 
type7 <- "Cancer" 
type8 <- "stromal"
type9 <- "Cancer"
type10 <- "stromal"
type11 <- "stromal"
type12 <- "Cancer"
type13 <- "stromal"
type14 <- "stromal"

#--------#
#差异表达基因可视化
#VlnPlot(test_data, features = c("CD4", "CD8A", "CD8B"))
new.cluster.ids <- c(type0,type1,type2,type3,type4,type5,type6,
                     type7,type8,type9,type10,type11,type12,type13,type14)
names(new.cluster.ids) <- levels(test_data)
test_data <- RenameIdents(test_data, new.cluster.ids)
#pdf("plot/Treated_umap注释信息.pdf",width = 8,height = 6)
DimPlot(test_data, reduction = "umap",label = TRUE, pt.size = 1)
dev.off()
test_data@meta.data$region<-test_data@active.ident
#------------#
#-------#
#library(celldex)
#library(SingleR)
# 安装依赖包
#install.packages(c("outliers", "doFuture"))
#BiocManager::install("singscore")
#BiocManager::install("SingleR")
#BiocManager::install("celldex")
#remove.packages("celldex")
#install.packages("C:/Users/Administrator/Downloads/celldex_1.0.0.tar.gz", repos = NULL, type = "source")
#save.image(file = "F:/fsh/ST_matrix2/ST/BRCA/ST_data/ST_pbmc.RData")
save(seurat_data,test_data,cortex,summarized_progeny_scores_df,MIA.results1,file = "F:/fsh/ST_matrix2/ST/Lung2/ST_data/ST_lung.RData")
#-------------#
#install.packages("scCATCH")
#devtools::install_github('ZJUFanLab/scCATCH')
#------对于取子集后的空间数据-------#
cortex  <- SCTransform(test_data, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
#------整合空间数据和单细胞数据-----#

anchors <- FindTransferAnchors(reference = allen_reference, query = cortex, normalization.method = "SCT")

predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$cell_type, prediction.assay = TRUE, weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay
#---将细胞类型映射到空间区域-----#
rownames(cortex)
DefaultAssay(cortex) <- "predictions"
#-----#
#-----#
p<-SpatialFeaturePlot(cortex, features = c("T cell","T cell-M"), pt.size.factor = 1.6, ncol = 2, crop = TRUE, alpha = c(0.1, 1))#

DimPlot(cortex, group.by = "cell_type", label = TRUE)
#-----------#
## 利用预测的细胞cluster，寻找高变特征（类似常规的高变基因 ）
cortex <- FindSpatiallyVariableFeatures(cortex, assay = "predictions", selection.method = "moransi", 
                                        features = rownames(cortex), r.metric = 5, slot = "data")

DECs = SpatiallyVariableFeatures(cortex,method = "moransi",
                                 assay = NULL,
                                 decreasing = TRUE)

top.clusters <- head(DECs, 4)
SpatialPlot(object = cortex, features = top.clusters, ncol = 2)
cortex@active.ident<-seurat_data@active.ident
## 2.1 获取切片图像坐标和Key：
img<- GetTissueCoordinates(cortex)
head(img)
## 2.2 分群信息；
ident<- Idents(cortex)
head(ident)
SpatialDimPlot(cortex, label = TRUE, label.size = 3)
## 2.3 合并切片图像坐标、分群信息；
imgId<- data.frame(img,ident)
head(imgId)
#------细胞注释-------#
library(scCATCH)
library(magrittr)
scCATCH_obj <- scCATCH::createscCATCH(seurat_data@assays$RNA$data, cluster = as.character(seurat_data@meta.data$seurat_clusters))
cellmatch<-cellmatch
cellmatch_new <- cellmatch[cellmatch$cancer %in% c("Prostate Cancer"),]
clu_markers1 <- scCATCH::findmarkergene(scCATCH_obj, species = "Human",
                                        cluster = "All", 
                                        cancer ="Prostate" ,
                                        if_use_custom_marker =  TRUE,
                                        marker = cellmatch,
                                        tissue = c("Prostate","Blood"),
                                        use_method = "1", cell_min_pct = 0.25,
                                        logfc = 0.25, pvalue = 0.25,
                                        verbose = TRUE)

clu_ann <- scCATCH::findcelltype(clu_markers1)
scCATCH_ann_data = merge(clu_ann@meta, clu_ann@celltype,by = "cluster") %>% tibble::column_to_rownames("cell")
cell<-data.frame(scCATCH_ann_data$cluster,scCATCH_ann_data$cell_type)
cell<-unique(cell)
#--singleR-------#
library(SingleR)
library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
hpca.se<-human_reference
meta=seurat_data@meta.data #pbmc的meta文件，包含了seurat的聚类结果
head(meta)
#进行singleR注释
pbmc_for_SingleR<-seurat_data[["RNA"]]$counts ##获取标准化矩阵
pbmc.hesc <- SingleR(test = pbmc_for_SingleR, ref = hpca.se, labels = hpca.se$label.main)#
pbmc.hesc

#seurat 和 singleR的table表
table(pbmc.hesc$labels,meta$seurat_clusters)
#------------------------#
### 找差异基因
# 控制三个阈值：logfc.threshold p_val_adj d
maintype_marker=FindAllMarkers(seurat_data,logfc.threshold = 0.5,only.pos = T)
maintype_marker=maintype_marker%>%filter(p_val_adj < 1e-04)
maintype_marker$d=maintype_marker$pct.1 - maintype_marker$pct.2
maintype_marker=maintype_marker%>%filter(d > 0.2)
maintype_marker=maintype_marker%>%arrange(cluster,desc(avg_log2FC))
maintype_marker=as.data.frame(maintype_marker)
table(maintype_marker$cluster)
### 整合坐标、region
#（此处把坐标类比单细胞转录组分析中的注释信息）
# 
#region.info=read.table("PDAC_A.region.txt",header = T,sep = "\t",stringsAsFactors = F)
test_data@meta.data$SB=test_data@active.ident
library(Seurat)
library(tidyverse)
### 找region特异基因
#Idents(test_data)="region"
region_marker=FindAllMarkers(test_data,logfc.threshold = 0,only.pos = T)
region_marker=region_marker%>%filter(p_val_adj < 1e-04)
region_marker$d=region_marker$pct.1 - region_marker$pct.2
region_marker=region_marker%>%filter(d > 0.05)
region_marker=region_marker%>%arrange(cluster,desc(avg_log2FC))
region_marker=as.data.frame(region_marker)
table(region_marker$cluster)
### 画图看看
top20 <- ST_pbmc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
#-------#
cellmarker <- read.csv("F:/fsh/ST_matrix2/ST/BRCA/ST_data/cellmarker_h.csv")
cellmarker <- cellmarker[cellmarker$tissue_class=="Lung",]
cellmarker <- merge(cellmarker,top20,by.x="Symbol",by.y ="gene")
cellmarker <- cellmarker[,c(1,3,6,8,10,22,25,26)]
cellmarker <- cellmarker[!duplicated(cellmarker),]
type0 <-"stromal" #"stromal"
type1 <-"stromal" #"Cancer"
type2 <-"stromal" #"Cancer" 
type3 <-"stromal" #"Cancer" 
type4 <-"stromal" #"stromal" 
type5 <-"stromal" #"Cancer" 
type6 <-"Cancer" ##"stromal" 
type7 <-"stromal" #"stromal" 
type8 <-"stromal" ##"stromal"
type9 <-"Cancer" #"stromal"
type10 <-"stromal" #"stromal"
type11 <-"Cancer" #"stromal"
#type12 <-"stromal" #"stromal"
#type13 <-"stromal" ##"stromal"
#type14 <-"stromal" #"stromal"
#type15 <-"stromal" #"Cancer"
#--------#
#差异表达基因可视化
#VlnPlot(test_data, features = c("CD4", "CD8A", "CD8B"))
new.cluster.ids <- c(type0,type1,type2,type3,type4,type5,type6,
                     type7,type8,type9,type10,type11)
names(new.cluster.ids) <- levels(test_data)
test_data <- RenameIdents(test_data, new.cluster.ids)
#pdf("plot/Treated_umap注释信息.pdf",width = 8,height = 6)
test_data@meta.data$region<-test_data@active.ident
SpatialDimPlot(test_data, label = TRUE, label.size = 3)
#----区域marker----#
region_marker=FindAllMarkers(test_data,logfc.threshold = 0.5,only.pos = T)
top20 <- region_marker %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
SpatialFeaturePlot(test_data, features = c("MUC2","S100A7","KLK5","KLK11","CHRNA9"), alpha = c(0.1, 1))
#-----TCGA临床验证----------#
BRCA_exp<-read.delim("F:/fsh/test/sample/TCGA-BRCA.htseq_counts.tsv", header=T)
clinical<-read.delim("F:/fsh/test/sample/clinical-all.txt", header=T)
clinical<-clinical[,c(4,5,16)]
clinical$bcr_sample_barcode<- gsub("-", ".", clinical$bcr_sample_barcode)

mart<-fread("F:/fsh/test/sample/mart_export.txt", header=T)
mart<-mart[,c(1,5)]
#----------------#
BRCA_exp$Ensembl_ID<- sub("\\..*$", "", BRCA_exp$Ensembl_ID)
# 将数据框转换为矩阵，然后取前出现的每个id
BRCA_exp <- BRCA_exp[!duplicated(BRCA_exp$Ensembl_ID), ]
#------------#
BRCA_exp_1<-merge(mart,BRCA_exp,by.x="Gene stable ID",by.y="Ensembl_ID")
BRCA_exp_1 <- BRCA_exp_1[!duplicated(BRCA_exp_1$`Gene name`), ]
BRCA_exp_1<-BRCA_exp_1[!BRCA_exp_1$`Gene name`=="",]
BRCA_exp_1<-data.frame(BRCA_exp_1)
rownames(BRCA_exp_1)<-BRCA_exp_1$Gene.name
BRCA_exp_1<-BRCA_exp_1[,-c(1,2)]
#-----------#
MUC2<-BRCA_exp_1["CHRNA9",]
# 假设expr_df是基因表达的数据框，行是基因，列是样本
# 首先计算每个基因表达的中位数
median_expr <- apply(MUC2, 1, median, na.rm = TRUE)

# 初始化一个新的数据框来存储分组信息
group_df <- data.frame(sample = colnames(MUC2), expression_group = character(length = ncol(MUC2)))

# 将每个样本的表达水平与中位数进行比较，确定分组
for (i in 1:ncol(MUC2)) {
  # 计算每个样本的表达水平是高于还是低于中位数
  group_df$expression_group[i] <- ifelse(rowMeans(MUC2[, i, drop = FALSE]) > median_expr, "High", "Low")
}

# 查看分组信息
table(group_df$expression_group)
#-----------------#

gene_clinical<-merge(clinical,group_df,by.x="bcr_sample_barcode",by.y="sample")
#-----生存分析-------#
library(survival)
library(survminer)
names(gene_clinical)<-c("sample","OS.time","OS","Group")
gene_clinical <- gene_clinical[!is.na(gene_clinical$OS.time), ]
sfit <- survfit(Surv(OS.time, OS)~Group, data=gene_clinical)
print(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)
