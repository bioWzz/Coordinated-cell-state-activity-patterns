library(monocle3)
library(dplyr)
library(Seurat)
library(patchwork)
library(dplyr)
library(CytoTRACE)
library(Matrix)
library(corrplot)
my_colors <- colorRampPalette(c("steelblue", "white", "#c11e21"))(100)
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/threecancersingle.RData")
#注释
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/skcm/skcmCD8Tpbmc.RData")
load("D:/细胞状态2/论文/分化/fig4d/skcm/skcmCD8Tpbmc.RData")


VlnPlot(skcmCD8Tpbmc, features = c("CD8A","HLA-E","CXCL10","IDO1","HLA-DRA","CXCL9","STAT1","GZMB","GZMH"))
VlnPlot(skcmCD8Tpbmc, features = c("IL2RA","CD3E","PRF1","CD2","VIM","NKG7","CD3D","LAPTM5","PTPRC","CD3G","GZMB","CD28","CD8A","GZMA"))
VlnPlot(skcmCD8Tpbmc, features = c("NKG7","CD8A","CCR7","CD3D","CD3G"))
VlnPlot(skcmCD8Tpbmc, features = c("HLA-DRB1","ENTPD1","KRT86","LY6E","IFI6","IFI27","HLA-DQA2","IFI44","OAS1","PRF1","HAVCR2","CXCL13","HLA-DRB5","HLA-DRA","ALOX5AP","MX1","GZMB","ISG15","CD82","IFIT3","IFITM3","BST2","GNLY","IFI44L","CCL3"))
VlnPlot(skcmCD8Tpbmc, features = c("CXCR6","CCR5"))
VlnPlot(skcmCD8Tpbmc, features = c("CCR7", "LEF1"))
VlnPlot(skcmCD8Tpbmc, features = c("NCR1","CD3D","GNLY","KLRD1"))
VlnPlot(skcmCD8Tpbmc, features = c("CXCL13","CD200","BCL6","IL21","PDCD1"))
VlnPlot(skcmCD8Tpbmc, features = c("ZEB2","CD4"))


head(skcmCD8Tpbmc)
skcmCD8Tpbmc
FindMarkers(skcmCD8Tpbmc, ident.1 = 2, min.pct = 0.001)
FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)

genes_to_check <- c("CD8A","HLA-E","CXCL10","IDO1","HLA-DRA","CXCL9","STAT1","GZMB","GZMH",
                    "IL2RA","CD3E","PRF1","CD2","VIM","NKG7","CD3D","LAPTM5","PTPRC","CD3G","GZMB","CD28","CD8A","GZMA",
                    "NKG7","CD8A","CCR7","CD3D","CD3G",
                    "HLA-DRB1","ENTPD1","KRT86","LY6E","IFI6","IFI27","HLA-DQA2","IFI44","OAS1","PRF1","HAVCR2","CXCL13","HLA-DRB5","HLA-DRA","ALOX5AP","MX1","GZMB","ISG15","CD82","IFIT3","IFITM3","BST2","GNLY","IFI44L","CCL3",
                    "CXCR6","CCR5",
                    "CCR7","LEF1",
                    "NCR1","CD3D","GNLY","KLRD1",
                    "CXCL13","CD200","BCL6","IL21","PDCD1",
                    "ZEB2","CD4")
unigene <- unique(genes_to_check)
genes_to_check <- as.factor(genes_to_check)
genes_to_check <- factor(genes_to_check,levels = unigene)
skcmCD8Tpbmc$seurat_clusters <- as.factor(skcmCD8Tpbmc$seurat_clusters)
DotPlot(skcmCD8Tpbmc,group.by = 'seurat_clusters', features = unique(genes_to_check),cluster.idents = T) + coord_flip()

DotPlot(skcmCD8Tpbmc, features = unique(genes_to_check),group.by = "seurat_clusters")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")

markers <- FindMarkers(skcmCD8Tpbmc, ident.1 = "1", min.pct = 0.25)

clusters <- unique(Idents(skcmCD8Tpbmc))
marker_genes <- list()
for (cluster in clusters) {
  markers <- FindMarkers(skcmCD8Tpbmc, ident.1 = cluster, min.pct = 0.25)
  marker_genes[[cluster]] <- markers
}







#自动注释
library(SingleR)
library(celldex)

#加载人类数据
hpca.se <- celldex::HumanPrimaryCellAtlasData()
hpca.se
hpca.se@metadata

bpe.se<-celldex::BlueprintEncodeData() 
bpe.se@metadata

scobj4SingleR <- GetAssayData(skcmCD8Tpbmc, slot = "data")

pred.scobj <- SingleR(test = scobj4SingleR, ref = hpca.se, assay.type.test=1,
                      labels = hpca.se$label.main)







#SKCMTprolif
#无注释
exprMat<-skcmTprolifpbmc[["RNA"]]$counts
exprMat <- exprMat[apply(exprMat > 0,1,sum) >= 5,]

meta.data<-skcmTprolifpbmc@meta.data
meta.data$cell<-rownames(meta.data)

meta.data1 <- meta.data

exprMat2 <- data.frame(exprMat)

meta.data1$CELL <- gsub("[-@]", ".",meta.data1$cell)
dim(exprMat2)
#-------CytoTRACE---------#
results <- CytoTRACE(exprMat2)
emb <- skcmTprolifpbmc@reductions[["umap"]]@cell.embeddings

rownames(emb)<-gsub("[-@]", ".",rownames(emb))

pheno <- meta.data1$seurat_clusters
pheno <- as.character(pheno)
names(pheno) <- rownames(meta.data1)
#names(pheno) <- gsub("[^a-zA-Z0-9]", ".",names(pheno))
names(pheno) <- gsub("[-@]", ".",names(pheno))

setdiff(names(pheno) ,colnames(exprMat2))
setdiff(colnames(exprMat2),names(pheno))
#------结果--------#
setwd("/wangzhenzenlab/fsh/xff/fenhua/zhushi/skcm/Tprolif")
plotCytoTRACE(results, phenotype = pheno, emb = emb)
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

#----monocle3-----#
f_data_cd4 <- data.frame(gene_short_name = row.names(exprMat2) ,row.names = row.names(exprMat2))
#exprMat<-skcmCD8Tpbmc[["RNA"]]$counts
#exprMat[1:3,1:3]
#exprMat2 <- exprMat[,!colnames(exprMat)%in%diffcell]
sorted_column_names <- colnames(exprMat2)
rownames(meta.data1) <- meta.data1$CELL
meta.data2 <- meta.data1[sorted_column_names, ]
head(meta.data2)
cds <- new_cell_data_set(as.matrix(exprMat2),
                         cell_metadata = meta.data2,
                         gene_metadata = f_data_cd4)
class(exprMat)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = 'PCA')
#这里是umap的藏身之地!!!!!!!!!!!!
#emb <- skcmCD8Tpbmc@reductions[["umap"]]@cell.embeddings
#emb <- emb[!rownames(emb)%in%diffcell,]
SingleCellExperiment::reducedDims(cds)[["UMAP"]]=emb
cds <- cluster_cells(cds, reduction_method ="UMAP")
#---------#
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition")

#table(colData(cds)$CELLTYPE)
#----自定义root-----#
## 定义函数 get_earliest_principal_node ：根据指定的簇（seurat_clusters）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin="6"){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "seurat_clusters"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}
nodes_vec <- c(get_earliest_principal_node(cds,"6"),
               get_earliest_principal_node(cds,"5"),
               get_earliest_principal_node(cds,"3")
)
cds = order_cells(cds, root_pr_nodes=nodes_vec,reduction_method = "UMAP")
#----------------------#
#--轨迹图-#
pdf("plot_cells.pdf")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1)
dev.off()

pseudotime1 <- as.data.frame(cds@principal_graph_aux$UMAP$pseudotime)
head(pseudotime1)  
pseudotime1$cell <- rownames(pseudotime1)
colnames(pseudotime1) <- c("pseudotime","cell")

pseudotime2 <- merge(pseudotime1,SKCMsingle,by.x="cell",by.y="Cell")
head(pseudotime2)
pseudotime2 <-pseudotime2[,c("cell","pseudotime","model4","model6","celltype")]
#rownames(pseudotime2) <- pseudotime2$cell
pseudotime2$cell <- NULL
pseudotime2 <- pseudotime2[!is.infinite(pseudotime2$pseudotime),]
table(pseudotime2$celltype)
save(pseudotime2,file = "pseudotime2.RData")
#计算相关性
singlemt <- pseudotime2[pseudotime2$celltype=="Tprolif",c("pseudotime","model4")]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
#0.077
pdf("cor4.pdf")
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
dev.off()














#SKCMMali
#无注释
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/skcm/skcmMalipbmc.RData")
exprMat<-skcmMalipbmc[["RNA"]]$counts
exprMat <- exprMat[apply(exprMat > 0,1,sum) >= 5,]

meta.data<-skcmMalipbmc@meta.data
meta.data$cell<-rownames(meta.data)

meta.data1 <- meta.data

exprMat2 <- data.frame(exprMat)

meta.data1$CELL <- gsub("[-@]", ".",meta.data1$cell)
dim(exprMat2)
#-------CytoTRACE---------#
results <- CytoTRACE(exprMat2)
emb <- skcmMalipbmc@reductions[["umap"]]@cell.embeddings

rownames(emb)<-gsub("[-@]", ".",rownames(emb))

pheno <- meta.data1$seurat_clusters
pheno <- as.character(pheno)
names(pheno) <- rownames(meta.data1)
#names(pheno) <- gsub("[^a-zA-Z0-9]", ".",names(pheno))
names(pheno) <- gsub("[-@]", ".",names(pheno))

setdiff(names(pheno) ,colnames(exprMat2))
setdiff(colnames(exprMat2),names(pheno))
#------结果--------#
setwd("/wangzhenzenlab/fsh/xff/fenhua/zhushi/skcm/Mali")
plotCytoTRACE(results, phenotype = pheno, emb = emb)
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

#----monocle3-----#
f_data_cd4 <- data.frame(gene_short_name = row.names(exprMat2) ,row.names = row.names(exprMat2))
#exprMat<-skcmCD8Tpbmc[["RNA"]]$counts
#exprMat[1:3,1:3]
#exprMat2 <- exprMat[,!colnames(exprMat)%in%diffcell]
sorted_column_names <- colnames(exprMat2)
rownames(meta.data1) <- meta.data1$CELL
meta.data2 <- meta.data1[sorted_column_names, ]
head(meta.data2)
cds <- new_cell_data_set(as.matrix(exprMat2),
                         cell_metadata = meta.data2,
                         gene_metadata = f_data_cd4)
class(exprMat)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = 'PCA')
#这里是umap的藏身之地!!!!!!!!!!!!
#emb <- skcmCD8Tpbmc@reductions[["umap"]]@cell.embeddings
#emb <- emb[!rownames(emb)%in%diffcell,]
SingleCellExperiment::reducedDims(cds)[["UMAP"]]=emb
cds <- cluster_cells(cds, reduction_method ="UMAP")
#---------#
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition")

#table(colData(cds)$CELLTYPE)
#----自定义root-----#
## 定义函数 get_earliest_principal_node ：根据指定的簇（seurat_clusters）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin="4"){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "seurat_clusters"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}
nodes_vec <- c(get_earliest_principal_node(cds,"4"),
               get_earliest_principal_node(cds,"2")
)
cds = order_cells(cds, root_pr_nodes=nodes_vec,reduction_method = "UMAP")
#----------------------#
#--轨迹图-#
pdf("plot_cells.pdf")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1)
dev.off()

pseudotime1 <- as.data.frame(cds@principal_graph_aux$UMAP$pseudotime)
head(pseudotime1)  
pseudotime1$cell <- rownames(pseudotime1)
colnames(pseudotime1) <- c("pseudotime","cell")

pseudotime2 <- merge(pseudotime1,SKCMsingle,by.x="cell",by.y="Cell")
head(pseudotime2)
pseudotime2 <-pseudotime2[,c("cell","pseudotime","model4","model6","celltype")]
#rownames(pseudotime2) <- pseudotime2$cell
pseudotime2$cell <- NULL
pseudotime2 <- pseudotime2[!is.infinite(pseudotime2$pseudotime),]
table(pseudotime2$celltype)
save(pseudotime2,file = "pseudotime2.RData")
#计算相关性
singlemt <- pseudotime2[pseudotime2$celltype=="Malignant",c("pseudotime","model6")]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
#0.127
pdf("cor6.pdf")
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
dev.off()












#SKCMFib
#无注释
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/skcm/skcmFibpbmc.RData")
exprMat<-skcmFibpbmc[["RNA"]]$counts
exprMat <- exprMat[apply(exprMat > 0,1,sum) >= 5,]

meta.data<-skcmFibpbmc@meta.data
meta.data$cell<-rownames(meta.data)

meta.data1 <- meta.data

exprMat2 <- data.frame(exprMat)

meta.data1$CELL <- gsub("[-@]", ".",meta.data1$cell)
dim(exprMat2)
#-------CytoTRACE---------#
results <- CytoTRACE(exprMat2)
emb <- skcmFibpbmc@reductions[["umap"]]@cell.embeddings

rownames(emb)<-gsub("[-@]", ".",rownames(emb))

pheno <- meta.data1$seurat_clusters
pheno <- as.character(pheno)
names(pheno) <- rownames(meta.data1)
#names(pheno) <- gsub("[^a-zA-Z0-9]", ".",names(pheno))
names(pheno) <- gsub("[-@]", ".",names(pheno))

setdiff(names(pheno) ,colnames(exprMat2))
setdiff(colnames(exprMat2),names(pheno))
#------结果--------#
setwd("/wangzhenzenlab/fsh/xff/fenhua/zhushi/skcm/Fib")
plotCytoTRACE(results, phenotype = pheno, emb = emb)
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

#----monocle3-----#
f_data_cd4 <- data.frame(gene_short_name = row.names(exprMat2) ,row.names = row.names(exprMat2))
#exprMat<-skcmCD8Tpbmc[["RNA"]]$counts
#exprMat[1:3,1:3]
#exprMat2 <- exprMat[,!colnames(exprMat)%in%diffcell]
sorted_column_names <- colnames(exprMat2)
rownames(meta.data1) <- meta.data1$CELL
meta.data2 <- meta.data1[sorted_column_names, ]
head(meta.data2)
cds <- new_cell_data_set(as.matrix(exprMat2),
                         cell_metadata = meta.data2,
                         gene_metadata = f_data_cd4)
class(exprMat)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = 'PCA')
#这里是umap的藏身之地!!!!!!!!!!!!
#emb <- skcmCD8Tpbmc@reductions[["umap"]]@cell.embeddings
#emb <- emb[!rownames(emb)%in%diffcell,]
SingleCellExperiment::reducedDims(cds)[["UMAP"]]=emb
cds <- cluster_cells(cds, reduction_method ="UMAP")
#---------#
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition")

#table(colData(cds)$CELLTYPE)
#----自定义root-----#
## 定义函数 get_earliest_principal_node ：根据指定的簇（seurat_clusters）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin="2"){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "seurat_clusters"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}
nodes_vec <- c(get_earliest_principal_node(cds,"2"),
               get_earliest_principal_node(cds,"3"),
               get_earliest_principal_node(cds,"7"),
               get_earliest_principal_node(cds,"8"),
               get_earliest_principal_node(cds,"1"),
               get_earliest_principal_node(cds,"0")
)
cds = order_cells(cds, root_pr_nodes=nodes_vec,reduction_method = "UMAP")
#----------------------#
#--轨迹图-#
pdf("plot_cells.pdf")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1)
dev.off()
save(cds,file = "cds.RData")
pseudotime1 <- as.data.frame(cds@principal_graph_aux$UMAP$pseudotime)
head(pseudotime1)  
pseudotime1$cell <- rownames(pseudotime1)
colnames(pseudotime1) <- c("pseudotime","cell")

pseudotime2 <- merge(pseudotime1,SKCMsingle,by.x="cell",by.y="Cell")
head(pseudotime2)
pseudotime2 <-pseudotime2[,c("cell","pseudotime","model4","model6","celltype")]
#rownames(pseudotime2) <- pseudotime2$cell
pseudotime2$cell <- NULL
pseudotime2 <- pseudotime2[!is.infinite(pseudotime2$pseudotime),]
table(pseudotime2$celltype)
save(pseudotime2,file = "pseudotime2.RData")
#计算相关性
singlemt <- pseudotime2[pseudotime2$celltype=="Fibroblasts",c("pseudotime","model6")]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
#-0.228
pdf("cor6.pdf")
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
dev.off()














#PAADFib
#无注释
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/paad/paadFibpbmc.RData")
exprMat<-paadFibpbmc[["RNA"]]$counts
exprMat <- exprMat[apply(exprMat > 0,1,sum) >= 5,]

meta.data<-paadFibpbmc@meta.data
meta.data$cell<-rownames(meta.data)

meta.data1 <- meta.data

exprMat2 <- data.frame(exprMat)
dim(exprMat2)
meta.data1$CELL <- gsub("[-@]", ".",meta.data1$cell)

#-------CytoTRACE---------#
results <- CytoTRACE(exprMat2)
emb <- paadFibpbmc@reductions[["umap"]]@cell.embeddings

rownames(emb)<-gsub("[-@]", ".",rownames(emb))

pheno <- meta.data1$seurat_clusters
pheno <- as.character(pheno)
names(pheno) <- rownames(meta.data1)
#names(pheno) <- gsub("[^a-zA-Z0-9]", ".",names(pheno))
names(pheno) <- gsub("[-@]", ".",names(pheno))

setdiff(names(pheno) ,colnames(exprMat2))
setdiff(colnames(exprMat2),names(pheno))
#------结果--------#
setwd("/wangzhenzenlab/fsh/xff/fenhua/zhushi/paad/Fib")
plotCytoTRACE(results, phenotype = pheno, emb = emb)
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

#----monocle3-----#
f_data_cd4 <- data.frame(gene_short_name = row.names(exprMat2) ,row.names = row.names(exprMat2))
#exprMat<-skcmCD8Tpbmc[["RNA"]]$counts
#exprMat[1:3,1:3]
#exprMat2 <- exprMat[,!colnames(exprMat)%in%diffcell]
sorted_column_names <- colnames(exprMat2)
rownames(meta.data1) <- meta.data1$CELL
meta.data2 <- meta.data1[sorted_column_names, ]
head(meta.data2)
cds <- new_cell_data_set(as.matrix(exprMat2),
                         cell_metadata = meta.data2,
                         gene_metadata = f_data_cd4)
class(exprMat)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = 'PCA')
#这里是umap的藏身之地!!!!!!!!!!!!
#emb <- skcmCD8Tpbmc@reductions[["umap"]]@cell.embeddings
#emb <- emb[!rownames(emb)%in%diffcell,]
SingleCellExperiment::reducedDims(cds)[["UMAP"]]=emb
cds <- cluster_cells(cds, reduction_method ="UMAP")
#---------#
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition")

#table(colData(cds)$CELLTYPE)
#----自定义root-----#
## 定义函数 get_earliest_principal_node ：根据指定的簇（seurat_clusters）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin="7"){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "seurat_clusters"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}
nodes_vec <- c(get_earliest_principal_node(cds,"7")
)
cds = order_cells(cds, root_pr_nodes=nodes_vec,reduction_method = "UMAP")
#----------------------#
#--轨迹图-#
pdf("plot_cells.pdf")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1)
dev.off()
save(cds,file = "cds.RData")

pseudotime1 <- as.data.frame(cds@principal_graph_aux$UMAP$pseudotime)
head(pseudotime1)  
pseudotime1$cell <- rownames(pseudotime1)
pseudotime1$cell <- gsub("[-_@]", ".",meta.data1$cell)
colnames(pseudotime1) <- c("pseudotime","cell")

pseudotime2 <- merge(pseudotime1,PAADsingle,by.x="cell",by.y="cell")
head(pseudotime2)
pseudotime2 <-pseudotime2[,c("cell","pseudotime","model1","celltype")]
#rownames(pseudotime2) <- pseudotime2$cell
pseudotime2$cell <- NULL
pseudotime2 <- pseudotime2[!is.infinite(pseudotime2$pseudotime),]
table(pseudotime2$celltype)
save(pseudotime2,file = "pseudotime2.RData")
#计算相关性
singlemt <- pseudotime2[pseudotime2$celltype=="Fibroblasts",c("pseudotime","model1")]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
#-0.481
#0.246
pdf("Fibcor1.pdf")
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
dev.off()
















#PAADmyoFib
#无注释
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/paad/paadmyoFibpbmc.RData")
#load("/wangzhenzenlab/fsh/xff/fenhua/paad/paadmyoFibpbmc.RData")
exprMat<-paadmyoFibpbmc[["RNA"]]$counts
exprMat <- exprMat[apply(exprMat > 0,1,sum) >= 5,]

meta.data<-paadmyoFibpbmc@meta.data
meta.data$cell<-rownames(meta.data)

meta.data1 <- meta.data

exprMat2 <- data.frame(exprMat)
dim(exprMat2)
meta.data1$CELL <- gsub("[-@]", ".",meta.data1$cell)

#-------CytoTRACE---------#
results <- CytoTRACE(exprMat2)
emb <- paadmyoFibpbmc@reductions[["umap"]]@cell.embeddings

rownames(emb)<-gsub("[-@]", ".",rownames(emb))

pheno <- meta.data1$seurat_clusters
pheno <- as.character(pheno)
names(pheno) <- rownames(meta.data1)
#names(pheno) <- gsub("[^a-zA-Z0-9]", ".",names(pheno))
names(pheno) <- gsub("[-@]", ".",names(pheno))

setdiff(names(pheno) ,colnames(exprMat2))
setdiff(colnames(exprMat2),names(pheno))
#------结果--------#
setwd("/wangzhenzenlab/fsh/xff/fenhua/zhushi/paad/Myofib")
plotCytoTRACE(results, phenotype = pheno, emb = emb)
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

#----monocle3-----#
f_data_cd4 <- data.frame(gene_short_name = row.names(exprMat2) ,row.names = row.names(exprMat2))
#exprMat<-skcmCD8Tpbmc[["RNA"]]$counts
#exprMat[1:3,1:3]
#exprMat2 <- exprMat[,!colnames(exprMat)%in%diffcell]
sorted_column_names <- colnames(exprMat2)
rownames(meta.data1) <- meta.data1$CELL
meta.data2 <- meta.data1[sorted_column_names, ]
head(meta.data2)
cds <- new_cell_data_set(as.matrix(exprMat2),
                         cell_metadata = meta.data2,
                         gene_metadata = f_data_cd4)
class(exprMat)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = 'PCA')
#这里是umap的藏身之地!!!!!!!!!!!!
#emb <- skcmCD8Tpbmc@reductions[["umap"]]@cell.embeddings
#emb <- emb[!rownames(emb)%in%diffcell,]
SingleCellExperiment::reducedDims(cds)[["UMAP"]]=emb
cds <- cluster_cells(cds, reduction_method ="UMAP")
#---------#
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition")

#table(colData(cds)$CELLTYPE)
#----自定义root-----#
## 定义函数 get_earliest_principal_node ：根据指定的簇（seurat_clusters）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin="8"){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "seurat_clusters"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}
nodes_vec <- c(get_earliest_principal_node(cds,"8")
)
cds = order_cells(cds, root_pr_nodes=nodes_vec,reduction_method = "UMAP")
#----------------------#
#--轨迹图-#
pdf("plot_cells.pdf")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1)
dev.off()
save(cds,file = "cds.RData")

pseudotime1 <- as.data.frame(cds@principal_graph_aux$UMAP$pseudotime)
head(pseudotime1)  
pseudotime1$cell <- rownames(pseudotime1)
pseudotime1$cell <- gsub("[-_@]", ".",meta.data1$cell)
colnames(pseudotime1) <- c("pseudotime","cell")

pseudotime2 <- merge(pseudotime1,PAADsingle,by.x="cell",by.y="cell")
head(pseudotime2)
pseudotime2 <-pseudotime2[,c("cell","pseudotime","model1","celltype")]
#rownames(pseudotime2) <- pseudotime2$cell
pseudotime2$cell <- NULL
pseudotime2 <- pseudotime2[!is.infinite(pseudotime2$pseudotime),]
table(pseudotime2$celltype)
save(pseudotime2,file = "pseudotime2.RData")
#计算相关性
singlemt <- pseudotime2[pseudotime2$celltype=="Myofibroblasts",c("pseudotime","model1")]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
#-0.481
#0.246
pdf("Fibcor1.pdf")
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
dev.off()

















#STADEnd
#无注释
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/stad/stadEndpbmc.RData")
exprMat<-stadEndpbmc[["RNA"]]$counts
exprMat <- exprMat[apply(exprMat > 0,1,sum) >= 5,]

meta.data<-stadEndpbmc@meta.data
meta.data$cell<-rownames(meta.data)

meta.data1 <- meta.data

exprMat2 <- data.frame(exprMat)
dim(exprMat2)
meta.data1$CELL <- gsub("[-@]", ".",meta.data1$cell)

#-------CytoTRACE---------#
results <- CytoTRACE(exprMat2)
emb <- stadEndpbmc@reductions[["umap"]]@cell.embeddings

rownames(emb)<-gsub("[-@]", ".",rownames(emb))

pheno <- meta.data1$seurat_clusters
pheno <- as.character(pheno)
names(pheno) <- rownames(meta.data1)
#names(pheno) <- gsub("[^a-zA-Z0-9]", ".",names(pheno))
names(pheno) <- gsub("[-@]", ".",names(pheno))

setdiff(names(pheno) ,colnames(exprMat2))
setdiff(colnames(exprMat2),names(pheno))
#------结果--------#
setwd("/wangzhenzenlab/fsh/xff/fenhua/zhushi/stad/End")
plotCytoTRACE(results, phenotype = pheno, emb = emb)
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

#----monocle3-----#
f_data_cd4 <- data.frame(gene_short_name = row.names(exprMat2) ,row.names = row.names(exprMat2))
sorted_column_names <- colnames(exprMat2)
rownames(meta.data1) <- meta.data1$CELL
meta.data2 <- meta.data1[sorted_column_names, ]
head(meta.data2)
cds <- new_cell_data_set(as.matrix(exprMat2),
                         cell_metadata = meta.data2,
                         gene_metadata = f_data_cd4)
class(exprMat)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = 'PCA')
#这里是umap的藏身之地!!!!!!!!!!!!
#emb <- skcmCD8Tpbmc@reductions[["umap"]]@cell.embeddings
#emb <- emb[!rownames(emb)%in%diffcell,]
SingleCellExperiment::reducedDims(cds)[["UMAP"]]=emb
cds <- cluster_cells(cds, reduction_method ="UMAP")
#---------#
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition")

#table(colData(cds)$CELLTYPE)
#----自定义root-----#
## 定义函数 get_earliest_principal_node ：根据指定的簇（seurat_clusters）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin="0"){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "seurat_clusters"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}
nodes_vec <- c(get_earliest_principal_node(cds,"0")
)
cds = order_cells(cds, root_pr_nodes=nodes_vec,reduction_method = "UMAP")
#----------------------#
#--轨迹图-#
pdf("plot_cells.pdf")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1)
dev.off()
save(cds,file = "cds.RData")

pseudotime1 <- as.data.frame(cds@principal_graph_aux$UMAP$pseudotime)
head(pseudotime1)  
pseudotime1$cell <- rownames(pseudotime1)
pseudotime1$cell <- gsub("[-_@]", ".",meta.data1$cell)
colnames(pseudotime1) <- c("pseudotime","cell")

pseudotime2 <- merge(pseudotime1,STADsingle,by.x="cell",by.y="cell")
head(pseudotime2)
pseudotime2 <-pseudotime2[,c("cell","pseudotime","model1","model3","model5","celltype")]
#rownames(pseudotime2) <- pseudotime2$cell
pseudotime2$cell <- NULL
pseudotime2 <- pseudotime2[!is.infinite(pseudotime2$pseudotime),]
table(pseudotime2$celltype)
save(pseudotime2,file = "pseudotime2.RData")
#计算相关性
singlemt <- pseudotime2[pseudotime2$celltype=="Endothelial",c("pseudotime","model5")]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
#0.216
#0.223
pdf("cor5.pdf")
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
dev.off()











#STADFib
#无注释
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/stad/stadFibpbmc.RData")
exprMat<-stadFibpbmc[["RNA"]]$counts
exprMat <- exprMat[apply(exprMat > 0,1,sum) >= 5,]

meta.data<-stadFibpbmc@meta.data
meta.data$cell<-rownames(meta.data)

meta.data1 <- meta.data

exprMat2 <- data.frame(exprMat)
dim(exprMat2)
meta.data1$CELL <- gsub("[-@]", ".",meta.data1$cell)

#-------CytoTRACE---------#
results <- CytoTRACE(exprMat2)
emb <- stadFibpbmc@reductions[["umap"]]@cell.embeddings

rownames(emb)<-gsub("[-@]", ".",rownames(emb))

pheno <- meta.data1$seurat_clusters
pheno <- as.character(pheno)
names(pheno) <- rownames(meta.data1)
#names(pheno) <- gsub("[^a-zA-Z0-9]", ".",names(pheno))
names(pheno) <- gsub("[-@]", ".",names(pheno))

setdiff(names(pheno) ,colnames(exprMat2))
setdiff(colnames(exprMat2),names(pheno))
#------结果--------#
setwd("/wangzhenzenlab/fsh/xff/fenhua/zhushi/stad/Fib")
plotCytoTRACE(results, phenotype = pheno, emb = emb)
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

#----monocle3-----#
f_data_cd4 <- data.frame(gene_short_name = row.names(exprMat2) ,row.names = row.names(exprMat2))
sorted_column_names <- colnames(exprMat2)
rownames(meta.data1) <- meta.data1$CELL
meta.data2 <- meta.data1[sorted_column_names, ]
head(meta.data2)
cds <- new_cell_data_set(as.matrix(exprMat2),
                         cell_metadata = meta.data2,
                         gene_metadata = f_data_cd4)
class(exprMat)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = 'PCA')
#这里是umap的藏身之地!!!!!!!!!!!!
#emb <- skcmCD8Tpbmc@reductions[["umap"]]@cell.embeddings
#emb <- emb[!rownames(emb)%in%diffcell,]
SingleCellExperiment::reducedDims(cds)[["UMAP"]]=emb
cds <- cluster_cells(cds, reduction_method ="UMAP")
#---------#
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition")

#table(colData(cds)$CELLTYPE)
#----自定义root-----#
## 定义函数 get_earliest_principal_node ：根据指定的簇（seurat_clusters）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin="0"){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "seurat_clusters"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}
nodes_vec <- c(get_earliest_principal_node(cds,"2"),
               get_earliest_principal_node(cds,"1")
)
cds = order_cells(cds, root_pr_nodes=nodes_vec,reduction_method = "UMAP")
#----------------------#
#--轨迹图-#
pdf("plot_cells.pdf")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1)
dev.off()
save(cds,file = "cds.RData")


pseudotime1 <- as.data.frame(cds@principal_graph_aux$UMAP$pseudotime)
head(pseudotime1)  
pseudotime1$cell <- rownames(pseudotime1)
pseudotime1$cell <- gsub("[-_@]", ".",meta.data1$cell)
colnames(pseudotime1) <- c("pseudotime","cell")

pseudotime2 <- merge(pseudotime1,STADsingle,by.x="cell",by.y="cell")
head(pseudotime2)
pseudotime2 <-pseudotime2[,c("cell","pseudotime","model1","model3","model5","celltype")]
#rownames(pseudotime2) <- pseudotime2$cell
pseudotime2$cell <- NULL
pseudotime2 <- pseudotime2[!is.infinite(pseudotime2$pseudotime),]
table(pseudotime2$celltype)
save(pseudotime2,file = "pseudotime2.RData")
#计算相关性
singlemt <- pseudotime2[pseudotime2$celltype=="Fibroblasts",c("pseudotime","model1")]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
#0.080
#-0.299
pdf("cor1.pdf")
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
dev.off()












#STADMono
#无注释
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/stad/stadMonopbmc.RData")
exprMat<-stadMonopbmc[["RNA"]]$counts
exprMat <- exprMat[apply(exprMat > 0,1,sum) >= 5,]

meta.data<-stadMonopbmc@meta.data
meta.data$cell<-rownames(meta.data)

meta.data1 <- meta.data

exprMat2 <- data.frame(exprMat)
dim(exprMat2)
meta.data1$CELL <- gsub("[-@]", ".",meta.data1$cell)

#-------CytoTRACE---------#
results <- CytoTRACE(exprMat2)
emb <- stadMonopbmc@reductions[["umap"]]@cell.embeddings

rownames(emb)<-gsub("[-@]", ".",rownames(emb))

pheno <- meta.data1$seurat_clusters
pheno <- as.character(pheno)
names(pheno) <- rownames(meta.data1)
#names(pheno) <- gsub("[^a-zA-Z0-9]", ".",names(pheno))
names(pheno) <- gsub("[-@]", ".",names(pheno))

setdiff(names(pheno) ,colnames(exprMat2))
setdiff(colnames(exprMat2),names(pheno))
#------结果--------#
setwd("/wangzhenzenlab/fsh/xff/fenhua/zhushi/stad/Mono")
plotCytoTRACE(results, phenotype = pheno, emb = emb)
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

#----monocle3-----#
f_data_cd4 <- data.frame(gene_short_name = row.names(exprMat2) ,row.names = row.names(exprMat2))
sorted_column_names <- colnames(exprMat2)
rownames(meta.data1) <- meta.data1$CELL
meta.data2 <- meta.data1[sorted_column_names, ]
head(meta.data2)
cds <- new_cell_data_set(as.matrix(exprMat2),
                         cell_metadata = meta.data2,
                         gene_metadata = f_data_cd4)
class(exprMat)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = 'PCA')
#这里是umap的藏身之地!!!!!!!!!!!!
#emb <- skcmCD8Tpbmc@reductions[["umap"]]@cell.embeddings
#emb <- emb[!rownames(emb)%in%diffcell,]
SingleCellExperiment::reducedDims(cds)[["UMAP"]]=emb
cds <- cluster_cells(cds, reduction_method ="UMAP")
#---------#
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition")

#table(colData(cds)$CELLTYPE)
#----自定义root-----#
## 定义函数 get_earliest_principal_node ：根据指定的簇（seurat_clusters）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin="2"){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "seurat_clusters"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}
nodes_vec <- c(get_earliest_principal_node(cds,"2")
)
cds = order_cells(cds, root_pr_nodes=nodes_vec,reduction_method = "UMAP")
#----------------------#
#--轨迹图-#
pdf("plot_cells.pdf")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1)
dev.off()

pseudotime1 <- as.data.frame(cds@principal_graph_aux$UMAP$pseudotime)
head(pseudotime1)  
pseudotime1$cell <- rownames(pseudotime1)
pseudotime1$cell <- gsub("[-_@]", ".",meta.data1$cell)
colnames(pseudotime1) <- c("pseudotime","cell")

pseudotime2 <- merge(pseudotime1,STADsingle,by.x="cell",by.y="cell")
head(pseudotime2)
pseudotime2 <-pseudotime2[,c("cell","pseudotime","model1","model3","model5","celltype")]
#rownames(pseudotime2) <- pseudotime2$cell
pseudotime2$cell <- NULL
pseudotime2 <- pseudotime2[!is.infinite(pseudotime2$pseudotime),]
table(pseudotime2$celltype)
save(pseudotime2,file = "pseudotime2.RData")
#计算相关性
singlemt <- pseudotime2[pseudotime2$celltype=="Mono/Macro",c("pseudotime","model5")]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
#0.047
#-0.079
#0.176
pdf("cor5.pdf")
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
dev.off()














#STADDC
#无注释
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/stad/stadDCpbmc.RData")
exprMat<-stadDCpbmc[["RNA"]]$counts
exprMat <- exprMat[apply(exprMat > 0,1,sum) >= 5,]

meta.data<-stadDCpbmc@meta.data
meta.data$cell<-rownames(meta.data)

meta.data1 <- meta.data

exprMat2 <- data.frame(exprMat)
dim(exprMat2)
meta.data1$CELL <- gsub("[-@]", ".",meta.data1$cell)

#-------CytoTRACE---------#
results <- CytoTRACE(exprMat2)
emb <- stadDCpbmc@reductions[["umap"]]@cell.embeddings

rownames(emb)<-gsub("[-@]", ".",rownames(emb))

pheno <- meta.data1$seurat_clusters
pheno <- as.character(pheno)
names(pheno) <- rownames(meta.data1)
#names(pheno) <- gsub("[^a-zA-Z0-9]", ".",names(pheno))
names(pheno) <- gsub("[-@]", ".",names(pheno))

setdiff(names(pheno) ,colnames(exprMat2))
setdiff(colnames(exprMat2),names(pheno))
#------结果--------#
setwd("/wangzhenzenlab/fsh/xff/fenhua/zhushi/stad/DC")
plotCytoTRACE(results, phenotype = pheno, emb = emb)
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

#----monocle3-----#
f_data_cd4 <- data.frame(gene_short_name = row.names(exprMat2) ,row.names = row.names(exprMat2))
sorted_column_names <- colnames(exprMat2)
rownames(meta.data1) <- meta.data1$CELL
meta.data2 <- meta.data1[sorted_column_names, ]
head(meta.data2)
cds <- new_cell_data_set(as.matrix(exprMat2),
                         cell_metadata = meta.data2,
                         gene_metadata = f_data_cd4)
class(exprMat)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = 'PCA')
#这里是umap的藏身之地!!!!!!!!!!!!
#emb <- skcmCD8Tpbmc@reductions[["umap"]]@cell.embeddings
#emb <- emb[!rownames(emb)%in%diffcell,]
SingleCellExperiment::reducedDims(cds)[["UMAP"]]=emb
cds <- cluster_cells(cds, reduction_method ="UMAP")
#---------#
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition")

#table(colData(cds)$CELLTYPE)
#----自定义root-----#
## 定义函数 get_earliest_principal_node ：根据指定的簇（seurat_clusters）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin="2"){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "seurat_clusters"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}
nodes_vec <- c(get_earliest_principal_node(cds,"3"),
               get_earliest_principal_node(cds,"6")
)
cds = order_cells(cds, root_pr_nodes=nodes_vec,reduction_method = "UMAP")
#----------------------#
#--轨迹图-#
pdf("plot_cells.pdf")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1)
dev.off()

pseudotime1 <- as.data.frame(cds@principal_graph_aux$UMAP$pseudotime)
head(pseudotime1)  
pseudotime1$cell <- rownames(pseudotime1)
pseudotime1$cell <- gsub("[-_@]", ".",meta.data1$cell)
colnames(pseudotime1) <- c("pseudotime","cell")

pseudotime2 <- merge(pseudotime1,STADsingle,by.x="cell",by.y="cell")
head(pseudotime2)
pseudotime2 <-pseudotime2[,c("cell","pseudotime","model1","model3","model5","celltype")]
#rownames(pseudotime2) <- pseudotime2$cell
pseudotime2$cell <- NULL
pseudotime2 <- pseudotime2[!is.infinite(pseudotime2$pseudotime),]
table(pseudotime2$celltype)
save(pseudotime2,file = "pseudotime2.RData")
#计算相关性
singlemt <- pseudotime2[pseudotime2$celltype=="DC",c("pseudotime","model5")]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
#-0.127
pdf("cor5.pdf")
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
dev.off()









#STADMast
#无注释
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/stad/stadMastpbmc.RData")
exprMat<-stadMastpbmc[["RNA"]]$counts
exprMat <- exprMat[apply(exprMat > 0,1,sum) >= 5,]

meta.data<-stadMastpbmc@meta.data
meta.data$cell<-rownames(meta.data)

meta.data1 <- meta.data

exprMat2 <- data.frame(exprMat)
dim(exprMat2)
meta.data1$CELL <- gsub("[-@]", ".",meta.data1$cell)

#-------CytoTRACE---------#
results <- CytoTRACE(exprMat2)
emb <- stadMastpbmc@reductions[["umap"]]@cell.embeddings

rownames(emb)<-gsub("[-@]", ".",rownames(emb))

pheno <- meta.data1$seurat_clusters
pheno <- as.character(pheno)
names(pheno) <- rownames(meta.data1)
#names(pheno) <- gsub("[^a-zA-Z0-9]", ".",names(pheno))
names(pheno) <- gsub("[-@]", ".",names(pheno))

setdiff(names(pheno) ,colnames(exprMat2))
setdiff(colnames(exprMat2),names(pheno))
#------结果--------#
setwd("/wangzhenzenlab/fsh/xff/fenhua/zhushi/stad/Mast")
plotCytoTRACE(results, phenotype = pheno, emb = emb)
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

#----monocle3-----#
f_data_cd4 <- data.frame(gene_short_name = row.names(exprMat2) ,row.names = row.names(exprMat2))
sorted_column_names <- colnames(exprMat2)
rownames(meta.data1) <- meta.data1$CELL
meta.data2 <- meta.data1[sorted_column_names, ]
head(meta.data2)
cds <- new_cell_data_set(as.matrix(exprMat2),
                         cell_metadata = meta.data2,
                         gene_metadata = f_data_cd4)
class(exprMat)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = 'PCA')
#这里是umap的藏身之地!!!!!!!!!!!!
#emb <- skcmCD8Tpbmc@reductions[["umap"]]@cell.embeddings
#emb <- emb[!rownames(emb)%in%diffcell,]
SingleCellExperiment::reducedDims(cds)[["UMAP"]]=emb
cds <- cluster_cells(cds, reduction_method ="UMAP")
#---------#
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition")

#table(colData(cds)$CELLTYPE)
#----自定义root-----#
## 定义函数 get_earliest_principal_node ：根据指定的簇（seurat_clusters）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin="1"){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "seurat_clusters"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}
nodes_vec <- c(get_earliest_principal_node(cds,"1")
)
cds = order_cells(cds, root_pr_nodes=nodes_vec,reduction_method = "UMAP")
#----------------------#
#--轨迹图-#
pdf("plot_cells.pdf")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1)
dev.off()
save(cds,file = "cds.RData")

pseudotime1 <- as.data.frame(cds@principal_graph_aux$UMAP$pseudotime)
head(pseudotime1)  
pseudotime1$cell <- rownames(pseudotime1)
pseudotime1$cell <- gsub("[-_@]", ".",meta.data1$cell)
colnames(pseudotime1) <- c("pseudotime","cell")

pseudotime2 <- merge(pseudotime1,STADsingle,by.x="cell",by.y="cell")
head(pseudotime2)
pseudotime2 <-pseudotime2[,c("cell","pseudotime","model1","model3","model5","celltype")]
#rownames(pseudotime2) <- pseudotime2$cell
pseudotime2$cell <- NULL
pseudotime2 <- pseudotime2[!is.infinite(pseudotime2$pseudotime),]
table(pseudotime2$celltype)
save(pseudotime2,file = "pseudotime2.RData")
#计算相关性
singlemt <- pseudotime2[pseudotime2$celltype=="Mast",c("pseudotime","model1")]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
#-0.136
#0.045
pdf("cor1.pdf")
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
dev.off()













#UCECEndothelial
#无注释
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/ucec/ucecEndothelialpbmc.RData")
exprMat<-ucecEndothelialpbmc[["RNA"]]$counts
exprMat <- exprMat[apply(exprMat > 0,1,sum) >= 5,]

meta.data<-ucecEndothelialpbmc@meta.data
meta.data$cell<-rownames(meta.data)

meta.data1 <- meta.data

exprMat2 <- data.frame(exprMat)
dim(exprMat2)
meta.data1$CELL <- gsub("[-@]", ".",meta.data1$cell)

#-------CytoTRACE---------#
results <- CytoTRACE(exprMat2)
emb <- ucecEndothelialpbmc@reductions[["umap"]]@cell.embeddings

rownames(emb)<-gsub("[-@]", ".",rownames(emb))

pheno <- meta.data1$seurat_clusters
pheno <- as.character(pheno)
names(pheno) <- rownames(meta.data1)
#names(pheno) <- gsub("[^a-zA-Z0-9]", ".",names(pheno))
names(pheno) <- gsub("[-@]", ".",names(pheno))

setdiff(names(pheno) ,colnames(exprMat2))
setdiff(colnames(exprMat2),names(pheno))
#------结果--------#
setwd("/wangzhenzenlab/fsh/xff/fenhua/zhushi/ucec/Endothelial")
plotCytoTRACE(results, phenotype = pheno, emb = emb)
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

#----monocle3-----#
f_data_cd4 <- data.frame(gene_short_name = row.names(exprMat2) ,row.names = row.names(exprMat2))
sorted_column_names <- colnames(exprMat2)
rownames(meta.data1) <- meta.data1$CELL
meta.data2 <- meta.data1[sorted_column_names, ]
head(meta.data2)
cds <- new_cell_data_set(as.matrix(exprMat2),
                         cell_metadata = meta.data2,
                         gene_metadata = f_data_cd4)
class(exprMat)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = 'PCA')
#这里是umap的藏身之地!!!!!!!!!!!!
#emb <- skcmCD8Tpbmc@reductions[["umap"]]@cell.embeddings
#emb <- emb[!rownames(emb)%in%diffcell,]
SingleCellExperiment::reducedDims(cds)[["UMAP"]]=emb
cds <- cluster_cells(cds, reduction_method ="UMAP")
#---------#
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition")

#table(colData(cds)$CELLTYPE)
#----自定义root-----#
## 定义函数 get_earliest_principal_node ：根据指定的簇（seurat_clusters）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin="1"){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "seurat_clusters"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}
nodes_vec <- c(get_earliest_principal_node(cds,"1")
)
cds = order_cells(cds, root_pr_nodes=nodes_vec,reduction_method = "UMAP")
#----------------------#
#--轨迹图-#
pdf("plot_cells.pdf")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1)
dev.off()
save(cds,file = "cds.RData")

pseudotime1 <- as.data.frame(cds@principal_graph_aux$UMAP$pseudotime)
head(pseudotime1)  
pseudotime1$cell <- rownames(pseudotime1)
pseudotime1$cell <- gsub("[-_@]", ".",meta.data1$cell)
colnames(pseudotime1) <- c("pseudotime","cell")

pseudotime2 <- merge(pseudotime1,UCECsingle,by.x="cell",by.y="cell")
head(pseudotime2)
pseudotime2 <-pseudotime2[,c("cell","pseudotime","model1","model3","model5","celltype")]
#rownames(pseudotime2) <- pseudotime2$cell
pseudotime2$cell <- NULL
pseudotime2 <- pseudotime2[!is.infinite(pseudotime2$pseudotime),]
table(pseudotime2$celltype)
save(pseudotime2,file = "pseudotime2.RData")
#计算相关性
singlemt <- pseudotime2[pseudotime2$celltype=="Endothelial",c("pseudotime","model1")]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
#-0.136
#0.045
pdf("cor1.pdf")
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
dev.off()

























#umap~model
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/threecancersingle.RData")

#skcmFib
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/skcm/skcmFibpbmc.RData")

skcmFibpbmc
skcmFibumap<- skcmFibpbmc@reductions$umap@cell.embeddings
skcmFibumap <- as.data.frame(skcmFibumap)
skcmFibumap$cell <- rownames(skcmFibumap)
head(skcmFibumap)
dim(skcmFibumap)
head(SKCMsingle)

skcmFibumap2 <- merge(skcmFibumap,SKCMsingle,by.x="cell",by.y="Cell")
dim(skcmFibumap2)
max(skcmFibumap2$model6)
min(skcmFibumap2$model6)
median(skcmFibumap2$model6)

ggplot(skcmFibumap2, aes(x = umap_1, y = umap_2, color = model6)) +
  geom_point(size = 1,alpha = 1) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_gradient2(low = "steelblue", mid = "white", high = "#c11e21",
                       midpoint = 0, limit = c(-0.2,0.2),oob = scales::squish)


skcmFibumap2$group6 <- ifelse(skcmFibumap2$model6>median(skcmFibumap2$model6),"high","low")
pdf("/wangzhenzenlab/fsh/xff/fenhua/zhushi/skcm/Fib/model6umap2.pdf",height = 4,width = 4)
ggplot(skcmFibumap2, aes(x = umap_1, y = umap_2, color = group6)) +
  geom_point(size = 1,alpha = 1) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_manual(values=c("low" = "steelblue", "high" = "#c11e21"))
dev.off()




#paadFib
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/paad/paadFibpbmc.RData")

paadFibpbmc
paadFibumap<- paadFibpbmc@reductions$umap@cell.embeddings
paadFibumap <- as.data.frame(paadFibumap)
paadFibumap$cell <- rownames(paadFibumap)
head(paadFibumap)
dim(paadFibumap)
head(PAADsingle)

paadFibumap2 <- merge(paadFibumap,PAADsingle,by.x="cell",by.y="Cell")
dim(paadFibumap2)
max(paadFibumap2$model1)
min(paadFibumap2$model1)
median(paadFibumap2$model1)

pdf("/wangzhenzenlab/fsh/xff/fenhua/zhushi/paad/Fib/model1umap.pdf",height = 4,width = 4)
ggplot(paadFibumap2, aes(x = umap_1, y = umap_2, color = model1)) +
  geom_point(size = 1,alpha = 1) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_gradient2(low = "#f0f623",mid = "#cd4a75", high = "#040081",
                        midpoint = 0.15, limit = c(-0.1,0.3),oob = scales::squish)
dev.off()


paadFibumap2$group1 <- ifelse(paadFibumap2$model1>median(paadFibumap2$model1),"high","low")
pdf("/wangzhenzenlab/fsh/xff/fenhua/zhushi/paad/Fib/model1umap2.pdf",height = 4,width = 4)
ggplot(paadFibumap2, aes(x = umap_1, y = umap_2, color = group1)) +
  geom_point(size = 1,alpha = 1) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_manual(values=c("low" = "steelblue", "high" = "#c11e21"))
dev.off()





#stadFib
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/stad/stadFibpbmc.RData")

stadFibpbmc
stadFibumap<- stadFibpbmc@reductions$umap@cell.embeddings
stadFibumap <- as.data.frame(stadFibumap)
stadFibumap$cell <- rownames(stadFibumap)
head(stadFibumap)
dim(stadFibumap)
head(STADsingle)

stadFibumap2 <- merge(stadFibumap,STADsingle,by.x="cell",by.y="Cell")
dim(stadFibumap2)
max(stadFibumap2$model5)
min(stadFibumap2$model5)
median(stadFibumap2$model5)

pdf("/wangzhenzenlab/fsh/xff/fenhua/zhushi/stad/Fib/model5umap.pdf",height = 4,width = 4)
ggplot(stadFibumap2, aes(x = umap_1, y = umap_2, color = model5)) +
  geom_point(size = 2,alpha = 1) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_gradient2(low = "#f0f623",mid = "#cd4a75", high = "#040081",
                        midpoint = -0.15, limit = c(-0.3,0),oob = scales::squish)
dev.off()


stadFibumap2$group5 <- ifelse(stadFibumap2$model5>median(stadFibumap2$model5),"high","low")
pdf("/wangzhenzenlab/fsh/xff/fenhua/zhushi/stad/Fib/model5umap2.pdf",height = 4,width = 4)
ggplot(stadFibumap2, aes(x = umap_1, y = umap_2, color = group5)) +
  geom_point(size = 1,alpha = 1) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_manual(values=c("low" = "steelblue", "high" = "#c11e21"))
dev.off()





#stadEnd
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/stad/stadEndpbmc.RData")

stadEndpbmc
stadEndumap<- stadEndpbmc@reductions$umap@cell.embeddings
stadEndumap <- as.data.frame(stadEndumap)
stadEndumap$cell <- rownames(stadEndumap)
head(stadEndumap)
dim(stadEndumap)
head(STADsingle)

stadEndumap2 <- merge(stadEndumap,STADsingle,by.x="cell",by.y="Cell")
dim(stadEndumap2)
max(stadEndumap2$model1)
min(stadEndumap2$model1)
median(stadEndumap2$model1)

pdf("/wangzhenzenlab/fsh/xff/fenhua/zhushi/stad/End/model1umap.pdf",height = 4,width = 4)
ggplot(stadEndumap2, aes(x = umap_1, y = umap_2, color = model1)) +
  geom_point(size = 2,alpha = 1) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_gradient2(low = "#f0f623",mid = "white", high = "#040081",
                        midpoint = 0.15, limit = c(0.05,0.25),oob = scales::squish)
dev.off()


stadEndumap2$group1 <- ifelse(stadEndumap2$model1>median(stadEndumap2$model1),"high","low")
pdf("/wangzhenzenlab/fsh/xff/fenhua/zhushi/stad/End/model1umap2.pdf",height = 4,width = 4)
ggplot(stadEndumap2, aes(x = umap_1, y = umap_2, color = group1)) +
  geom_point(size = 1,alpha = 1) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_manual(values=c("low" = "steelblue", "high" = "#c11e21"))
dev.off()
