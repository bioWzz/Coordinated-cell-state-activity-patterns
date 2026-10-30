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

#CESC
load("D:/细胞状态2/论文/fig3/fig3d/cesc/cescpbmc.RData")

exprMat<-cesccd8tpbmc[["RNA"]]$counts
exprMat<-cescmalipbmc[["RNA"]]$counts
exprMat<-cescmonopbmc[["RNA"]]$counts

exprMat <- exprMat[apply(exprMat > 0,1,sum) >= 5,]
# 读取注释文件
annotation_data <- CESCsingle
table(annotation_data$celltype)
rownames(annotation_data)<-annotation_data$Cell
meta.data<-cescmonopbmc@meta.data
meta.data$cell<-rownames(meta.data)
#查看是否有相同的细胞类型
cell<-intersect(rownames(meta.data),annotation_data$Cell)
#提取相同的细胞类型
cell_inf1<-annotation_data[cell,]
cell_inf2<-data.frame(cell_inf1$Cell,cell_inf1$celltype)
names(cell_inf2)<-c("cell","celltype")
#提取相同细胞类型的meta.data文件
meta.data1<-merge(meta.data,cell_inf2,by.x = "cell",by.y = "cell")
rownames(meta.data1)<-meta.data1$cell
meta.data1[1:3,1:3]
table(meta.data1$celltype)
exprMat2 <- data.frame(exprMat)

#-------CytoTRACE---------#
results <- CytoTRACE(exprMat2)
emb <- cescmonopbmc@reductions[["umap"]]@cell.embeddings
rownames(emb)<-gsub("[^a-zA-Z0-9]", ".",rownames(emb))

pheno <- meta.data1$celltype
pheno <- as.character(pheno)
names(pheno) <- rownames(meta.data1)
names(pheno) <- gsub("[^a-zA-Z0-9]", ".",names(pheno))
#------结果--------#
setwd("D:/细胞状态2/论文/fig3/fig3d/cesc/Mono/Macro")
plotCytoTRACE(results, phenotype = pheno, emb = emb)
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

#----monocle3-----#
f_data_cd4 <- data.frame(gene_short_name = row.names(cescmonopbmc) ,row.names = row.names(cescmonopbmc))
exprMat<-cescmonopbmc[["RNA"]]$counts
exprMat[1:3,1:3]
sorted_column_names <- colnames(exprMat)
meta.data2 <- meta.data1[sorted_column_names, ]
head(meta.data2)
cds <- new_cell_data_set(exprMat,
                         cell_metadata = meta.data2,
                         gene_metadata = f_data_cd4)
class(exprMat)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = 'PCA')
#这里是umap的藏身之地!!!!!!!!!!!!
SingleCellExperiment::reducedDims(cds)[["UMAP"]]=emb
cds <- cluster_cells(cds, reduction_method ="UMAP")
#---------#
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition")

#----自定义root-----#
## 定义函数 get_earliest_principal_node ：根据指定的细胞类型（celltype）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin="Mono/Macro"){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "celltype"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}
nodes_vec <- c(get_earliest_principal_node(cds,"Mono/Macro"))
cds = order_cells(cds, root_pr_nodes=nodes_vec,reduction_method = "UMAP")
#----------------------#
#--轨迹图-#
pdf("D:/细胞状态2/论文/fig3/fig3d/cesc/Mono/Macro/plot_cells.pdf")
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

pseudotime2 <- merge(pseudotime1,CESCsingle,by.x="cell",by.y="Cell")
head(pseudotime2)
pseudotime2 <- pseudotime2[,c("cell","pseudotime","model1","model3","celltype")]
rownames(pseudotime2) <- pseudotime2$cell
pseudotime2$cell <- NULL
pseudotime2 <- pseudotime2[!is.infinite(pseudotime2$pseudotime),]
table(pseudotime2$celltype)
save(pseudotime2,file = "D:/细胞状态2/论文/fig3/fig3d/cesc/Mono/Macro/pseudotime2.RData")
#计算相关性
singlemt <- pseudotime2[,c("pseudotime","model3")]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson

corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)













#ESCA
load("D:/细胞状态2/论文/fig3/fig3d/esca/escapbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/esca/escaEndpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/esca/escaEpipbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/esca/escaFibpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/esca/escaMastpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/esca/escaMonopbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/esca/escaCD8Tpbmc.RData")

exprMat<-escaBpbmc[["RNA"]]$counts
exprMat<-escaCD4Tpbmc[["RNA"]]$counts
exprMat<-escaCD8Tpbmc[["RNA"]]$counts
exprMat<-escaDCpbmc[["RNA"]]$counts
exprMat<-escaEndpbmc[["RNA"]]$counts
exprMat<-escaEpipbmc[["RNA"]]$counts
exprMat<-escaFibpbmc[["RNA"]]$counts
exprMat<-escaMastpbmc[["RNA"]]$counts
exprMat<-escaMonopbmc[["RNA"]]$counts


#escaBpbmc       Mono/Macro
exprMat <- exprMat[apply(exprMat > 0,1,sum) >= 5,]
# 读取注释文件
annotation_data <- ESCAsingle

duplicated_row_names <- duplicated(annotation_data$Cell) | annotation_data$Cell %in% duplicated(annotation_data$Cell)
# 保留第一个出现的行，删除其余重复的行
annotation_data <- annotation_data[!duplicated_row_names & !is.na(annotation_data$Cell), ]
# 更新行名
rownames(annotation_data) <- annotation_data$Cell

table(annotation_data$celltype)
#rownames(annotation_data)<-annotation_data$Cell
meta.data<-escaMonopbmc@meta.data
meta.data$cell<-rownames(meta.data)
#查看是否有相同的细胞类型
cell<-intersect(rownames(meta.data),annotation_data$Cell)
#提取相同的细胞类型
cell_inf1<-annotation_data[cell,]
cell_inf2<-data.frame(cell_inf1$Cell,cell_inf1$celltype)
names(cell_inf2)<-c("cell","celltype")
#提取相同细胞类型的meta.data文件
meta.data1<-merge(meta.data,cell_inf2,by.x = "cell",by.y = "cell")
rownames(meta.data1)<-meta.data1$cell
meta.data1[1:3,]
table(meta.data1$celltype)
exprMat2 <- data.frame(exprMat)
#-------CytoTRACE---------#
results <- CytoTRACE(exprMat2)
emb <- escaMonopbmc@reductions[["umap"]]@cell.embeddings
rownames(emb)<-gsub("[^a-zA-Z0-9]", ".",rownames(emb))

#meta.data$celltype <- "Mono/Macro"
pheno <- meta.data$celltype
pheno <- as.character(pheno)
names(pheno) <- rownames(meta.data)
names(pheno) <- gsub("[^a-zA-Z0-9]", ".",names(pheno))

colnames(exprMat2)<- gsub("[^a-zA-Z0-9]", ".",colnames(exprMat2))
rownames(meta.data)<- gsub("[^a-zA-Z0-9]", ".",rownames(meta.data))
length(intersect(rownames(emb) ,sorted_column_names))
setdiff(rownames(emb) ,colnames(exprMat2))
sorted_column_names <- colnames(exprMat2)
emb <- emb[sorted_column_names,]
class(exprMat2)
meta.data
emb <- as.data.frame(emb)
emb <- as.matrix(emb)
#------结果--------#
setwd("D:/细胞状态2/论文/fig3/fig3d/esca/Mono/Macro")
plotCytoTRACE(results, phenotype = pheno, emb = emb)
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

#----monocle3-----#
f_data_cd4 <- data.frame(gene_short_name = row.names(escaMonopbmc) ,row.names = row.names(escaMonopbmc))
exprMat<-escaMonopbmc[["RNA"]]$counts
exprMat[1:3,1:3]
colnames(exprMat)<- gsub("[^a-zA-Z0-9]", ".",colnames(exprMat))
sorted_column_names <- colnames(exprMat)
meta.data2 <- meta.data[sorted_column_names, ]
rownames(meta.data1) <- gsub("[^a-zA-Z0-9]", ".",rownames(meta.data1))
head(meta.data2)
rownames(meta.data2)<- gsub("[^a-zA-Z0-9]", ".",rownames(meta.data2))
intersect(rownames(meta.data2),colnames(exprMat))
setdiff(rownames(meta.data2),colnames(exprMat))
cds <- new_cell_data_set(exprMat,
                         cell_metadata = meta.data2,
                         gene_metadata = f_data_cd4)
class(exprMat)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = 'PCA')
#这里是umap的藏身之地!!!!!!!!!!!!
SingleCellExperiment::reducedDims(cds)[["UMAP"]]=emb
cds <- cluster_cells(cds, reduction_method ="UMAP")
#---------#
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition")

#----自定义root-----#
## 定义函数 get_earliest_principal_node ：根据指定的细胞类型（celltype）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin="Mono/Macro"){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "celltype"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}
nodes_vec <- c(get_earliest_principal_node(cds,"Mono/Macro"))
cds = order_cells(cds, root_pr_nodes=nodes_vec,reduction_method = "UMAP")
#----------------------#
#--轨迹图-#
pdf("D:/细胞状态2/论文/fig3/fig3d/esca/Mono/Macro/plot_cells.pdf")
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

pseudotime2 <- merge(pseudotime1,ESCAsingle,by.x="cell",by.y="cell")
head(pseudotime2)
pseudotime2 <- pseudotime2[,c("cell","pseudotime","model4","celltype")]
rownames(pseudotime2) <- pseudotime2$cell
pseudotime2$cell <- NULL
pseudotime2 <- pseudotime2[!is.infinite(pseudotime2$pseudotime),]
table(pseudotime2$celltype)
save(pseudotime2,file = "D:/细胞状态2/论文/fig3/fig3d/esca/Mono/Macro/pseudotime2.RData")
#计算相关性
singlemt <- pseudotime2[,c("pseudotime","model4")]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson

corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)










#kirc
load("D:/细胞状态2/论文/fig3/fig3d/kirc/kircBpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/kirc/kircCD4Tpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/kirc/kircCD8Tpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/kirc/kircDCpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/kirc/kircEndpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/kirc/kircEpipbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/kirc/kircFibpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/kirc/kircMastpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/kirc/kircMonopbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/kirc/kircNKpbmc.RData")

exprMat<-kircBpbmc[["RNA"]]$counts
exprMat<-kircCD4Tpbmc[["RNA"]]$counts
exprMat<-kircCD8Tpbmc[["RNA"]]$counts
exprMat<-kircDCpbmc[["RNA"]]$counts
exprMat<-kircEndpbmc[["RNA"]]$counts
exprMat<-kircEpipbmc[["RNA"]]$counts
exprMat<-kircFibpbmc[["RNA"]]$counts
exprMat<-kircMastpbmc[["RNA"]]$counts
exprMat<-kircMonopbmc[["RNA"]]$counts
exprMat<-kircNKpbmc[["RNA"]]$counts


exprMat <- exprMat[apply(exprMat > 0,1,sum) >= 5,]
# 读取注释文件
annotation_data <- KIRCsingle

duplicated_row_names <- duplicated(annotation_data$Cell) | annotation_data$Cell %in% duplicated(annotation_data$Cell)
# 保留第一个出现的行，删除其余重复的行
annotation_data <- annotation_data[!duplicated_row_names & !is.na(annotation_data$Cell), ]
# 更新行名
rownames(annotation_data) <- annotation_data$Cell

table(annotation_data$celltype)
#rownames(annotation_data)<-annotation_data$Cell
meta.data<-kircBpbmc@meta.data
meta.data$cell<-rownames(meta.data)
#查看是否有相同的细胞类型
cell<-intersect(rownames(meta.data),annotation_data$Cell)
#提取相同的细胞类型
cell_inf1<-annotation_data[cell,]
cell_inf2<-data.frame(cell_inf1$Cell,cell_inf1$celltype)
names(cell_inf2)<-c("cell","celltype")
#提取相同细胞类型的meta.data文件
meta.data1<-merge(meta.data,cell_inf2,by.x = "cell",by.y = "cell")
rownames(meta.data1)<-meta.data1$cell
meta.data1[1:3,]
table(meta.data1$celltype)
exprMat2 <- data.frame(exprMat)
#-------CytoTRACE---------#
results <- CytoTRACE(exprMat2)
emb <- kircBpbmc@reductions[["umap"]]@cell.embeddings
rownames(emb)<-gsub("[-@]", ".",rownames(emb))


pheno <- meta.data1$celltype
pheno <- as.character(pheno)
names(pheno) <- rownames(meta.data1)
names(pheno) <- gsub("[-@]", ".",names(pheno))

setdiff(names(pheno) ,colnames(exprMat2))
#------结果--------#
setwd("D:/细胞状态2/论文/fig3/fig3d/kirc/B")
plotCytoTRACE(results, phenotype = pheno, emb = emb)
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

#----monocle3-----#
f_data_cd4 <- data.frame(gene_short_name = row.names(kircBpbmc) ,row.names = row.names(kircBpbmc))
exprMat<-kircBpbmc[["RNA"]]$counts
exprMat[1:3,1:3]
sorted_column_names <- colnames(exprMat)
meta.data2 <- meta.data1[sorted_column_names, ]
head(meta.data2)
cds <- new_cell_data_set(exprMat,
                         cell_metadata = meta.data2,
                         gene_metadata = f_data_cd4)
class(exprMat)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = 'PCA')
#这里是umap的藏身之地!!!!!!!!!!!!
emb <- kircBpbmc@reductions[["umap"]]@cell.embeddings
SingleCellExperiment::reducedDims(cds)[["UMAP"]]=emb
cds <- cluster_cells(cds, reduction_method ="UMAP")
#---------#
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition")

#----自定义root-----#
## 定义函数 get_earliest_principal_node ：根据指定的细胞类型（celltype）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin="B"){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "celltype"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}
nodes_vec <- c(get_earliest_principal_node(cds,"B"))
cds = order_cells(cds, root_pr_nodes=nodes_vec,reduction_method = "UMAP")
#----------------------#
#--轨迹图-#
pdf("D:/细胞状态2/论文/fig3/fig3d/kirc/B/plot_cells.pdf")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1)
dev.off()

save(cds,file = "D:/细胞状态2/论文/fig3/fig3d/kirc/B/cds.RData")

pseudotime1 <- as.data.frame(cds@principal_graph_aux$UMAP$pseudotime)
head(pseudotime1)  
pseudotime1$cell <- rownames(pseudotime1)
colnames(pseudotime1) <- c("pseudotime","cell")

pseudotime2 <- merge(pseudotime1,KIRCsingle,by.x="cell",by.y="Cell")
head(pseudotime2)
pseudotime2 <-pseudotime2[,c("cell","pseudotime","model1","model2","model4","model6","celltype")]
rownames(pseudotime2) <- pseudotime2$cell
pseudotime2$cell <- NULL
pseudotime2 <- pseudotime2[!is.infinite(pseudotime2$pseudotime),]
table(pseudotime2$celltype)
save(pseudotime2,file = "D:/细胞状态2/论文/fig3/fig3d/kirc/B/pseudotime2.RData")
#计算相关性
singlemt <- pseudotime2[,c("pseudotime","model6")]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson

corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)


a


















#paad
load("D:/细胞状态2/论文/fig3/fig3d/paad/paadBpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/paad/paadCD8Tpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/paad/paadDCpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/paad/paadEndpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/paad/paadEpipbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/paad/paadFibpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/paad/paadMastpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/paad/paadMonopbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/paad/paadNKpbmc.RData")
#load("D:/细胞状态2/论文/分化/fig4d/paad/paadFibpbmc.RData")
#load("D:/细胞状态2/论文/分化/fig4d/paad/paadmyoFibpbmc.RData")

exprMat<-paadBpbmc[["RNA"]]$counts
exprMat<-paadCD8Tpbmc[["RNA"]]$counts
exprMat<-paadDCpbmc[["RNA"]]$counts
exprMat<-paadEndpbmc[["RNA"]]$counts
exprMat<-paadEpipbmc[["RNA"]]$counts
exprMat<-paadFibpbmc[["RNA"]]$counts
exprMat<-paadmyoFibpbmc[["RNA"]]$counts
exprMat<-paadMastpbmc[["RNA"]]$counts
exprMat<-paadMonopbmc[["RNA"]]$counts
exprMat<-paadNKpbmc[["RNA"]]$counts


load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/paad/paadFibpbmc.RData")
#load("/wangzhenzenlab/fsh/xff/fenhua/paad/paadmyoFibpbmc.RData")
exprMat<-paadmyoFibpbmc[["RNA"]]$counts
exprMat <- exprMat[apply(exprMat > 0,1,sum) >= 5,]
# 读取注释文件
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/threecancersingle.RData")
annotation_data <- PAADsingle

duplicated_row_names <- duplicated(annotation_data$Cell) | annotation_data$Cell %in% duplicated(annotation_data$Cell)
# 保留第一个出现的行，删除其余重复的行
annotation_data <- annotation_data[!duplicated_row_names & !is.na(annotation_data$Cell), ]
# 更新行名
rownames(annotation_data) <- annotation_data$Cell

table(annotation_data$celltype)
#rownames(annotation_data)<-annotation_data$Cell
meta.data<-paadmyoFibpbmc@meta.data
meta.data$cell<-rownames(meta.data)
#查看是否有相同的细胞类型
cell<-intersect(rownames(meta.data),annotation_data$Cell)
#提取相同的细胞类型
cell_inf1<-annotation_data[cell,]
cell_inf2<-data.frame(cell_inf1$Cell,cell_inf1$celltype,cell_inf1$model1)
names(cell_inf2)<-c("cell","celltype","model1")
#提取相同细胞类型的meta.data文件
meta.data1<-merge(meta.data,cell_inf2,by.x = "cell",by.y = "cell")
rownames(meta.data1)<-meta.data1$cell
meta.data1[1:3,]
meta.data1$group1 <- ifelse(meta.data1$model1>median(meta.data1$model1),"high","low")
meta.data1$CELLTYPE <- paste(meta.data1$celltype,"_",meta.data1$group1,sep = "")
table(meta.data1$CELLTYPE)
exprMat2 <- data.frame(exprMat)
#-------CytoTRACE---------#
results <- CytoTRACE(exprMat2)
emb <- paadmyoFibpbmc@reductions[["umap"]]@cell.embeddings
#rownames(emb)<-gsub("[^a-zA-Z0-9]", ".",rownames(emb))
rownames(emb)<-gsub("[@-]", ".",rownames(emb))

pheno <- meta.data1$CELLTYPE
pheno <- as.character(pheno)
names(pheno) <- rownames(meta.data1)
#names(pheno) <- gsub("[^a-zA-Z0-9]", ".",names(pheno))
names(pheno) <- gsub("[@-]", ".",names(pheno))

setdiff(names(pheno) ,colnames(exprMat2))
#------结果--------#
setwd("/wangzhenzenlab/fsh/xff/fenhua/paad/Myofibroblasts")
plotCytoTRACE(results, phenotype = pheno, emb = emb)
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

#----monocle3-----#
f_data_cd4 <- data.frame(gene_short_name = row.names(paadmyoFibpbmc) ,row.names = row.names(paadmyoFibpbmc))
exprMat<-paadmyoFibpbmc[["RNA"]]$counts
exprMat[1:3,1:3]
sorted_column_names <- colnames(exprMat)
meta.data2 <- meta.data1[sorted_column_names, ]
head(meta.data2)
cds <- new_cell_data_set(exprMat,
                         cell_metadata = meta.data2,
                         gene_metadata = f_data_cd4)
class(exprMat)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = 'PCA')
#这里是umap的藏身之地!!!!!!!!!!!!
SingleCellExperiment::reducedDims(cds)[["UMAP"]]=emb
cds <- cluster_cells(cds, reduction_method ="UMAP")
#---------#
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition")

#----自定义root-----#
## 定义函数 get_earliest_principal_node ：根据指定的细胞类型（celltype）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin="Fibroblasts_high"){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "CELLTYPE"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}
nodes_vec <- c(get_earliest_principal_node(cds,"Fibroblasts_high"),
               get_earliest_principal_node(cds,"Fibroblasts_low"),
               get_earliest_principal_node(cds,"Myofibroblasts_high")
               )
cds = order_cells(cds, root_pr_nodes=nodes_vec,reduction_method = "UMAP")
#----------------------#
#--轨迹图-#
pdf("/wangzhenzenlab/fsh/xff/fenhua/paad/Myofibroblasts/plot_cells.pdf")
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

pseudotime2 <- merge(pseudotime1,PAADsingle,by.x="cell",by.y="Cell")
head(pseudotime2)
pseudotime2 <-pseudotime2[,c("cell","pseudotime","model1","celltype")]
rownames(pseudotime2) <- pseudotime2$cell
pseudotime2$cell <- NULL
pseudotime2 <- pseudotime2[!is.infinite(pseudotime2$pseudotime),]
table(pseudotime2$celltype)
save(pseudotime2,file = "/wangzhenzenlab/fsh/xff/fenhua/paad/Myofibroblasts/pseudotime2.RData")
#计算相关性
singlemt <- pseudotime2[pseudotime2$celltype=="Myofibroblasts",c("pseudotime","model1")]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
pdf("/wangzhenzenlab/fsh/xff/fenhua/paad/Myofibroblasts/Fibroblastscor.pdf")
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
dev.off()

a















#skcm
load("D:/细胞状态2/论文/fig3/fig3d/skcm/skcmBpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/skcm/skcmCD4Tpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/skcm/skcmCD8Tpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/skcm/skcmDCpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/skcm/skcmEndpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/skcm/skcmFibpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/skcm/skcmMalipbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/skcm/skcmMastpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/skcm/skcmMonopbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/skcm/skcmNKpbmc.RData")

exprMat<-skcmBpbmc[["RNA"]]$counts
exprMat<-skcmCD4Tpbmc[["RNA"]]$counts
exprMat<-skcmCD8Tpbmc[["RNA"]]$counts
exprMat<-skcmDCpbmc[["RNA"]]$counts
exprMat<-skcmEndpbmc[["RNA"]]$counts
exprMat<-skcmFibpbmc[["RNA"]]$counts
exprMat<-skcmMalipbmc[["RNA"]]$counts
exprMat<-skcmMastpbmc[["RNA"]]$counts
exprMat<-skcmMonopbmc[["RNA"]]$counts
exprMat<-skcmNKpbmc[["RNA"]]$counts



#CD8T
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/skcm/skcmCD8Tpbmc.RData")
exprMat<-skcmCD8Tpbmc[["RNA"]]$counts
exprMat <- exprMat[apply(exprMat > 0,1,sum) >= 5,]
# 读取注释文件
annotation_data <- SKCMsingle

duplicated_row_names <- duplicated(annotation_data$Cell) | annotation_data$Cell %in% duplicated(annotation_data$Cell)
# 保留第一个出现的行，删除其余重复的行
annotation_data <- annotation_data[!duplicated_row_names & !is.na(annotation_data$Cell), ]
# 更新行名
rownames(annotation_data) <- annotation_data$Cell

table(annotation_data$celltype)
#rownames(annotation_data)<-annotation_data$Cell
meta.data<-skcmCD8Tpbmc@meta.data
meta.data$cell<-rownames(meta.data)
#查看是否有相同的细胞类型
annotation_data <- annotation_data[annotation_data$celltype%in%c("Tprolif","CD8T","CD8Tex"),]
cell<-intersect(rownames(meta.data),annotation_data$Cell)
#提取相同的细胞类型
cell_inf1<-annotation_data[cell,]
cell_inf2<-data.frame(cell_inf1$Cell,cell_inf1$celltype,cell_inf1$model4)
names(cell_inf2)<-c("cell","celltype","model4")
#提取相同细胞类型的meta.data文件
meta.data1<-merge(meta.data,cell_inf2,by.x = "cell",by.y = "cell")
diffcell<- setdiff(meta.data$cell,meta.data1$cell)
diffcell1<- setdiff(meta.data1$cell,meta.data$cell)

rownames(meta.data1)<-meta.data1$cell
meta.data1[1:3,]
meta.data1$group4 <- ifelse(meta.data1$model4>median(meta.data1$model4),"high","low")
meta.data1$CELLTYPE <- paste(meta.data1$celltype,"_",meta.data1$group4,sep = "")
table(meta.data1$CELLTYPE)
table(meta.data1$celltype)
exprMat2 <- data.frame(exprMat)
exprMat2 <- exprMat2[,!colnames(exprMat2)%in%diffcell]
diffcell2<- setdiff(meta.data1$cell,colnames(exprMat2))
diffcell3 <- setdiff(colnames(exprMat2),meta.data1$cell)
exprMat2[1:3,1:3]
meta.data1$CELL <- gsub("[-@]", ".",meta.data1$cell)
diffcell5<- setdiff(meta.data1$CELL,colnames(exprMat2))
diffcell6 <- setdiff(colnames(exprMat2),meta.data1$CELL)
exprMat2 <- exprMat2[,!colnames(exprMat2)%in%diffcell6]
#-------CytoTRACE---------#
results <- CytoTRACE(exprMat2)
emb <- skcmCD8Tpbmc@reductions[["umap"]]@cell.embeddings
emb <- emb[!rownames(emb)%in%diffcell,]
#rownames(emb)<-gsub("[^a-zA-Z0-9]", ".",rownames(emb))
rownames(emb)<-gsub("[-@]", ".",rownames(emb))

pheno <- meta.data1$CELLTYPE
pheno <- as.character(pheno)
names(pheno) <- rownames(meta.data1)
#names(pheno) <- gsub("[^a-zA-Z0-9]", ".",names(pheno))
names(pheno) <- gsub("[-@]", ".",names(pheno))

setdiff(names(pheno) ,colnames(exprMat2))
setdiff(colnames(exprMat2),names(pheno))
#------结果--------#
setwd("/wangzhenzenlab/fsh/xff/fenhua/skcm/Tprolif")
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
## 定义函数 get_earliest_principal_node ：根据指定的细胞类型（celltype）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin="CD8Tex_high"){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "CELLTYPE"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}
nodes_vec <- c(get_earliest_principal_node(cds,"CD8Tex_high"),
               get_earliest_principal_node(cds,"CD8T_high"),
               get_earliest_principal_node(cds,"CD8T_low"),
               get_earliest_principal_node(cds,"CD8Tex_low"),
               get_earliest_principal_node(cds,"Tprolif_high")
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

pdf("cor4.pdf")
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
dev.off()





#Fib
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/skcm/skcmFibpbmc.RData")
exprMat<-skcmFibpbmc[["RNA"]]$counts
exprMat <- exprMat[apply(exprMat > 0,1,sum) >= 5,]
# 读取注释文件
annotation_data <- SKCMsingle

duplicated_row_names <- duplicated(annotation_data$Cell) | annotation_data$Cell %in% duplicated(annotation_data$Cell)
# 保留第一个出现的行，删除其余重复的行
annotation_data <- annotation_data[!duplicated_row_names & !is.na(annotation_data$Cell), ]
# 更新行名
rownames(annotation_data) <- annotation_data$Cell

table(annotation_data$celltype)
#rownames(annotation_data)<-annotation_data$Cell
meta.data<-skcmFibpbmc@meta.data
meta.data$cell<-rownames(meta.data)
#查看是否有相同的细胞类型
annotation_data <- annotation_data[annotation_data$celltype%in%c("Fibroblasts","Myofibroblasts"),]
cell<-intersect(rownames(meta.data),annotation_data$Cell)
#提取相同的细胞类型
cell_inf1<-annotation_data[cell,]
cell_inf2<-data.frame(cell_inf1$Cell,cell_inf1$celltype,cell_inf1$model6)
names(cell_inf2)<-c("cell","celltype","model6")
#提取相同细胞类型的meta.data文件
meta.data1<-merge(meta.data,cell_inf2,by.x = "cell",by.y = "cell")
diffcell<- setdiff(meta.data$cell,meta.data1$cell)
diffcell1<- setdiff(meta.data1$cell,meta.data$cell)

rownames(meta.data1)<-meta.data1$cell
meta.data1[1:3,]
meta.data1$group6 <- ifelse(meta.data1$model6>median(meta.data1$model6),"high","low")
meta.data1$CELLTYPE <- paste(meta.data1$celltype,"_",meta.data1$group6,sep = "")
table(meta.data1$CELLTYPE)
table(meta.data1$celltype)
exprMat2 <- data.frame(exprMat)
exprMat2 <- exprMat2[,!colnames(exprMat2)%in%diffcell]
#-------CytoTRACE---------#
results <- CytoTRACE(exprMat2)
emb <- skcmFibpbmc@reductions[["umap"]]@cell.embeddings
emb <- emb[!rownames(emb)%in%diffcell,]
#rownames(emb)<-gsub("[^a-zA-Z0-9]", ".",rownames(emb))
rownames(emb)<-gsub("[-@]", ".",rownames(emb))

pheno <- meta.data1$CELLTYPE
pheno <- as.character(pheno)
names(pheno) <- rownames(meta.data1)
#names(pheno) <- gsub("[^a-zA-Z0-9]", ".",names(pheno))
names(pheno) <- gsub("[-@]", ".",names(pheno))

setdiff(names(pheno) ,colnames(exprMat2))
#------结果--------#
setwd("/wangzhenzenlab/fsh/xff/fenhua/skcm/Fibroblasts")
plotCytoTRACE(results, phenotype = pheno, emb = emb)
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

#----monocle3-----#
f_data_cd4 <- data.frame(gene_short_name = row.names(exprMat2) ,row.names = row.names(exprMat2))
exprMat<-skcmFibpbmc[["RNA"]]$counts
exprMat[1:3,1:3]
exprMat2 <- exprMat[,!colnames(exprMat)%in%diffcell]
sorted_column_names <- colnames(exprMat2)
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
emb <- skcmFibpbmc@reductions[["umap"]]@cell.embeddings
emb <- emb[!rownames(emb)%in%diffcell,]
SingleCellExperiment::reducedDims(cds)[["UMAP"]]=emb
cds <- cluster_cells(cds, reduction_method ="UMAP")
#---------#
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition")

#----自定义root-----#
## 定义函数 get_earliest_principal_node ：根据指定的细胞类型（celltype）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin="Fibroblasts_low"){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "CELLTYPE"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}
nodes_vec <- c(get_earliest_principal_node(cds,"Fibroblasts_low"),
               get_earliest_principal_node(cds,"Malignant_low"))
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
singlemt <- pseudotime2[pseudotime2$celltype=="Myofibroblasts",c("pseudotime","model6")]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson

pdf("cor6.pdf")
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
dev.off()










#Mali
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/skcm/skcmMalipbmc.RData")
exprMat<-skcmMalipbmc[["RNA"]]$counts
exprMat <- exprMat[apply(exprMat > 0,1,sum) >= 5,]
# 读取注释文件
annotation_data <- SKCMsingle

duplicated_row_names <- duplicated(annotation_data$Cell) | annotation_data$Cell %in% duplicated(annotation_data$Cell)
# 保留第一个出现的行，删除其余重复的行
annotation_data <- annotation_data[!duplicated_row_names & !is.na(annotation_data$Cell), ]
# 更新行名
rownames(annotation_data) <- annotation_data$Cell

table(annotation_data$celltype)
#rownames(annotation_data)<-annotation_data$Cell
meta.data<-skcmMalipbmc@meta.data
meta.data$cell<-rownames(meta.data)
#查看是否有相同的细胞类型
annotation_data <- annotation_data[annotation_data$celltype%in%c("Malignant"),]
cell<-intersect(meta.data$cell,annotation_data$Cell)
#提取相同的细胞类型
cell_inf1<-annotation_data[cell,]
cell_inf2<-data.frame(cell_inf1$Cell,cell_inf1$celltype,cell_inf1$model6)
names(cell_inf2)<-c("cell","celltype","model6")
#提取相同细胞类型的meta.data文件
meta.data1<-merge(meta.data,cell_inf2,by.x = "cell",by.y = "cell")
diffcell<- setdiff(meta.data$cell,meta.data1$cell)
diffcell1<- setdiff(meta.data1$cell,meta.data$cell)

rownames(meta.data1)<-meta.data1$cell
meta.data1[1:3,]
meta.data1$group6 <- ifelse(meta.data1$model6>median(meta.data1$model6),"high","low")
meta.data1$CELLTYPE <- paste(meta.data1$celltype,"_",meta.data1$group6,sep = "")
table(meta.data1$CELLTYPE)
table(meta.data1$celltype)
exprMat2 <- data.frame(exprMat)
exprMat2 <- exprMat2[,!colnames(exprMat2)%in%diffcell]
#-------CytoTRACE---------#
results <- CytoTRACE(exprMat2)
emb <- skcmMalipbmc@reductions[["umap"]]@cell.embeddings
emb <- emb[!rownames(emb)%in%diffcell,]
#rownames(emb)<-gsub("[^a-zA-Z0-9]", ".",rownames(emb))
rownames(emb)<-gsub("[-@]", ".",rownames(emb))

pheno <- meta.data1$CELLTYPE
pheno <- as.character(pheno)
names(pheno) <- rownames(meta.data1)
#names(pheno) <- gsub("[^a-zA-Z0-9]", ".",names(pheno))
names(pheno) <- gsub("[-@]", ".",names(pheno))

setdiff(names(pheno) ,colnames(exprMat2))
#------结果--------#
setwd("/wangzhenzenlab/fsh/xff/fenhua/skcm/Malignant")
plotCytoTRACE(results, phenotype = pheno, emb = emb)
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

#----monocle3-----#
f_data_cd4 <- data.frame(gene_short_name = row.names(exprMat2) ,row.names = row.names(exprMat2))
#exprMat<-skcmMalipbmc[["RNA"]]$counts
#exprMat[1:3,1:3]
#exprMat2 <- exprMat[,!colnames(exprMat)%in%diffcell]
sorted_column_names <- colnames(exprMat2)
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
emb <- skcmMalipbmc@reductions[["umap"]]@cell.embeddings
emb <- emb[!rownames(emb)%in%diffcell,]
SingleCellExperiment::reducedDims(cds)[["UMAP"]]=emb
cds <- cluster_cells(cds, reduction_method ="UMAP")
#---------#
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition")

#----自定义root-----#
## 定义函数 get_earliest_principal_node ：根据指定的细胞类型（celltype）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin="Malignant_low"){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "CELLTYPE"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}
nodes_vec <- c(get_earliest_principal_node(cds,"Malignant_low"))
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

pdf("cor6.pdf")
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
dev.off()

a


















#stad
load("D:/细胞状态2/论文/fig3/fig3d/stad/stadBpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/stad/stadCD8Tpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/stad/stadDCpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/stad/stadEndpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/stad/stadEpipbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/stad/stadFibpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/stad/stadMastpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/stad/stadMonopbmc.RData")


exprMat<-stadBpbmc[["RNA"]]$counts
exprMat<-stadCD8Tpbmc[["RNA"]]$counts
exprMat<-stadDCpbmc[["RNA"]]$counts
exprMat<-stadEndpbmc[["RNA"]]$counts
exprMat<-stadEpipbmc[["RNA"]]$counts
exprMat<-stadFibpbmc[["RNA"]]$counts
exprMat<-stadMastpbmc[["RNA"]]$counts
exprMat<-stadMonopbmc[["RNA"]]$counts


#End
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/stad/stadEndpbmc.RData")
exprMat<-stadEndpbmc[["RNA"]]$counts
exprMat <- exprMat[apply(exprMat > 0,1,sum) >= 5,]
# 读取注释文件
annotation_data <- STADsingle

duplicated_row_names <- duplicated(annotation_data$Cell) | annotation_data$Cell %in% duplicated(annotation_data$Cell)
# 保留第一个出现的行，删除其余重复的行
annotation_data <- annotation_data[!duplicated_row_names & !is.na(annotation_data$Cell), ]
# 更新行名
rownames(annotation_data) <- annotation_data$Cell
table(annotation_data$celltype)

#rownames(annotation_data)<-annotation_data$Cell
meta.data<-stadEndpbmc@meta.data
meta.data$cell<-rownames(meta.data)
#查看是否有相同的细胞类型
cell<-intersect(rownames(meta.data),annotation_data$Cell)
#提取相同的细胞类型
cell_inf1<-annotation_data[cell,]
cell_inf2<-data.frame(cell_inf1$Cell,cell_inf1$celltype,cell_inf1$model1,cell_inf1$EMT,cell_inf1$Angiogenesis,cell_inf1$Invasion,cell_inf1$Metastasis)
names(cell_inf2)<-c("cell","celltype","model1","EMT","Angiogenesis","Invasion","Metastasis")
head(cell_inf2)
#提取相同细胞类型的meta.data文件
meta.data1<-merge(meta.data,cell_inf2,by.x = "cell",by.y = "cell")
rownames(meta.data1)<-meta.data1$cell
meta.data1[1:3,]
table(meta.data1$celltype)
meta.data1$group1 <- ifelse(meta.data1$model1>median(meta.data1$model1),"high","low")
meta.data1$CELLTYPE <- paste(meta.data1$celltype,"_",meta.data1$group1,sep = "")
table(meta.data1$CELLTYPE)
table(meta.data1$group1)

exprMat2 <- data.frame(exprMat)
#-------CytoTRACE---------#
results <- CytoTRACE(exprMat2)
emb <- stadEndpbmc@reductions[["umap"]]@cell.embeddings
#rownames(emb)<-gsub("[^a-zA-Z0-9]", ".",rownames(emb))
rownames(emb)<-gsub("[@-]", ".",rownames(emb))

pheno <- meta.data1$CELLTYPE
pheno <- as.character(pheno)
names(pheno) <- rownames(meta.data1)
#names(pheno) <- gsub("[^a-zA-Z0-9]", ".",names(pheno))
names(pheno) <- gsub("[@-]", ".",names(pheno))

setdiff(names(pheno) ,colnames(exprMat2))
#------结果--------#
setwd("/wangzhenzenlab/fsh/xff/fenhua/stad/Endothelial")
plotCytoTRACE(results, phenotype = pheno, emb = emb)
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

#----monocle3-----#
f_data_cd4 <- data.frame(gene_short_name = row.names(stadEndpbmc) ,row.names = row.names(stadEndpbmc))
exprMat<-stadEndpbmc[["RNA"]]$counts
exprMat[1:3,1:3]
sorted_column_names <- colnames(exprMat)
meta.data2 <- meta.data1[sorted_column_names, ]
head(meta.data2)
table(meta.data2$group1)
cds <- new_cell_data_set(exprMat,
                         cell_metadata = meta.data2,
                         gene_metadata = f_data_cd4)
class(exprMat)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = 'PCA')
#这里是umap的藏身之地!!!!!!!!!!!!
emb <- stadEndpbmc@reductions[["umap"]]@cell.embeddings
SingleCellExperiment::reducedDims(cds)[["UMAP"]]=emb
cds <- cluster_cells(cds, reduction_method ="UMAP")
#---------#
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition")

#----自定义root-----#
## 定义函数 get_earliest_principal_node ：根据指定的细胞类型（CELLTYPE）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin="Endothelial_high"){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "CELLTYPE"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}
nodes_vec <- c(get_earliest_principal_node(cds,"Endothelial_high"))
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

pseudotime2 <- merge(pseudotime1,STADsingle,by.x="cell",by.y="Cell")
head(pseudotime2)
pseudotime2 <-pseudotime2[,c("cell","pseudotime","model1","model3","model5","celltype")]
rownames(pseudotime2) <- pseudotime2$cell
pseudotime2$cell <- NULL
pseudotime2 <- pseudotime2[!is.infinite(pseudotime2$pseudotime),]
table(pseudotime2$celltype)
save(pseudotime2,file = "pseudotime2.RData")
#计算相关性
singlemt <- pseudotime2[,c("pseudotime","model1")]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
pdf("model1cor.pdf")
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
dev.off()






#Fib
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/stad/stadFibpbmc.RData")
exprMat<-stadFibpbmc[["RNA"]]$counts
exprMat <- exprMat[apply(exprMat > 0,1,sum) >= 5,]
# 读取注释文件
annotation_data <- STADsingle

duplicated_row_names <- duplicated(annotation_data$Cell) | annotation_data$Cell %in% duplicated(annotation_data$Cell)
# 保留第一个出现的行，删除其余重复的行
annotation_data <- annotation_data[!duplicated_row_names & !is.na(annotation_data$Cell), ]
# 更新行名
rownames(annotation_data) <- annotation_data$Cell
table(annotation_data$celltype)

#rownames(annotation_data)<-annotation_data$Cell
meta.data<-stadFibpbmc@meta.data
meta.data$cell<-rownames(meta.data)
#查看是否有相同的细胞类型
cell<-intersect(rownames(meta.data),annotation_data$Cell)
#提取相同的细胞类型
cell_inf1<-annotation_data[cell,]
cell_inf2<-data.frame(cell_inf1$Cell,cell_inf1$celltype,cell_inf1$model5,cell_inf1$`Cell Cycle`,cell_inf1$`DNA damage`,cell_inf1$`DNA repair`)
names(cell_inf2)<-c("cell","celltype","model5","Cell Cycle","DNA damage","DNA repair")
head(cell_inf2)
#提取相同细胞类型的meta.data文件
meta.data1<-merge(meta.data,cell_inf2,by.x = "cell",by.y = "cell")
rownames(meta.data1)<-meta.data1$cell
meta.data1[1:3,]
table(meta.data1$celltype)
meta.data1$group5 <- ifelse(meta.data1$model5>median(meta.data1$model5),"high","low")
meta.data1$CELLTYPE <- paste(meta.data1$celltype,"_",meta.data1$group5,sep = "")
table(meta.data1$group5)
table(meta.data1$CELLTYPE)

exprMat2 <- data.frame(exprMat)
#-------CytoTRACE---------#
results <- CytoTRACE(exprMat2)
emb <- stadFibpbmc@reductions[["umap"]]@cell.embeddings
#rownames(emb)<-gsub("[^a-zA-Z0-9]", ".",rownames(emb))
rownames(emb)<-gsub("[@-]", ".",rownames(emb))

pheno <- meta.data1$CELLTYPE
pheno <- as.character(pheno)
names(pheno) <- rownames(meta.data1)
#names(pheno) <- gsub("[^a-zA-Z0-9]", ".",names(pheno))
names(pheno) <- gsub("[@-]", ".",names(pheno))

setdiff(names(pheno) ,colnames(exprMat2))
#------结果--------#
setwd("/wangzhenzenlab/fsh/xff/fenhua/stad/Fibroblasts")
plotCytoTRACE(results, phenotype = pheno, emb = emb)
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

#----monocle3-----#
f_data_cd4 <- data.frame(gene_short_name = row.names(stadFibpbmc) ,row.names = row.names(stadFibpbmc))
exprMat<-stadFibpbmc[["RNA"]]$counts
exprMat[1:3,1:3]
sorted_column_names <- colnames(exprMat)
meta.data2 <- meta.data1[sorted_column_names, ]
head(meta.data2)
table(meta.data2$group5)
cds <- new_cell_data_set(exprMat,
                         cell_metadata = meta.data2,
                         gene_metadata = f_data_cd4)
class(exprMat)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = 'PCA')
#这里是umap的藏身之地!!!!!!!!!!!!
emb <- stadFibpbmc@reductions[["umap"]]@cell.embeddings
SingleCellExperiment::reducedDims(cds)[["UMAP"]]=emb
cds <- cluster_cells(cds, reduction_method ="UMAP")
#---------#
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition")

#----自定义root-----#
## 定义函数 get_earliest_principal_node ：根据指定的细胞类型（celltype）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin="Fibroblasts_low"){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "CELLTYPE"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}
nodes_vec <- c(get_earliest_principal_node(cds,"Fibroblasts_low"))
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

pseudotime2 <- merge(pseudotime1,STADsingle,by.x="cell",by.y="Cell")
head(pseudotime2)
pseudotime2 <-pseudotime2[,c("cell","pseudotime","model1","model3","model5","celltype")]
rownames(pseudotime2) <- pseudotime2$cell
pseudotime2$cell <- NULL
pseudotime2 <- pseudotime2[!is.infinite(pseudotime2$pseudotime),]
table(pseudotime2$celltype)
save(pseudotime2,file = "pseudotime2.RData")
#计算相关性
singlemt <- pseudotime2[,c("pseudotime","model1")]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson

pdf("model1cor.pdf")
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
dev.off()



#Mono
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/stad/stadMonopbmc.RData")
exprMat<-stadMonopbmc[["RNA"]]$counts
exprMat <- exprMat[apply(exprMat > 0,1,sum) >= 5,]
# 读取注释文件
annotation_data <- STADsingle

duplicated_row_names <- duplicated(annotation_data$Cell) | annotation_data$Cell %in% duplicated(annotation_data$Cell)
# 保留第一个出现的行，删除其余重复的行
annotation_data <- annotation_data[!duplicated_row_names & !is.na(annotation_data$Cell), ]
# 更新行名
rownames(annotation_data) <- annotation_data$Cell
table(annotation_data$celltype)

#rownames(annotation_data)<-annotation_data$Cell
meta.data<-stadMonopbmc@meta.data
meta.data$cell<-rownames(meta.data)
#查看是否有相同的细胞类型
cell<-intersect(rownames(meta.data),annotation_data$Cell)
#提取相同的细胞类型
cell_inf1<-annotation_data[cell,]
cell_inf2<-data.frame(cell_inf1$Cell,cell_inf1$celltype,cell_inf1$model5,cell_inf1$model1,cell_inf1$model3)
names(cell_inf2)<-c("cell","celltype","model5","model1","model3")
head(cell_inf2)
#提取相同细胞类型的meta.data文件
meta.data1<-merge(meta.data,cell_inf2,by.x = "cell",by.y = "cell")
rownames(meta.data1)<-meta.data1$cell
meta.data1[1:3,]
table(meta.data1$celltype)
meta.data1$group1 <- ifelse(meta.data1$model1>median(meta.data1$model1),"high","low")
meta.data1$group3 <- ifelse(meta.data1$model3>median(meta.data1$model3),"high","low")
meta.data1$group5 <- ifelse(meta.data1$model5>median(meta.data1$model5),"high","low")
meta.data1$CELLTYPE1 <- paste(meta.data1$celltype,"_",meta.data1$group1,sep = "")
meta.data1$CELLTYPE3 <- paste(meta.data1$celltype,"_",meta.data1$group3,sep = "")
meta.data1$CELLTYPE5 <- paste(meta.data1$celltype,"_",meta.data1$group5,sep = "")
table(meta.data1$group5)
table(meta.data1$CELLTYPE1)

exprMat2 <- data.frame(exprMat)
#-------CytoTRACE---------#
results <- CytoTRACE(exprMat2)
emb <- stadMonopbmc@reductions[["umap"]]@cell.embeddings
#rownames(emb)<-gsub("[^a-zA-Z0-9]", ".",rownames(emb))
rownames(emb)<-gsub("[@-]", ".",rownames(emb))

pheno <- meta.data1$CELLTYPE1
#pheno <- meta.data1$CELLTYPE3
#pheno <- meta.data1$CELLTYPE5
pheno <- as.character(pheno)
names(pheno) <- rownames(meta.data1)
#names(pheno) <- gsub("[^a-zA-Z0-9]", ".",names(pheno))
names(pheno) <- gsub("[@-]", ".",names(pheno))

setdiff(names(pheno) ,colnames(exprMat2))
#------结果--------#
setwd("/wangzhenzenlab/fsh/xff/fenhua/stad/Mono/Macro/model5")
plotCytoTRACE(results, phenotype = pheno, emb = emb)
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

#----monocle3-----#
f_data_cd4 <- data.frame(gene_short_name = row.names(stadMonopbmc) ,row.names = row.names(stadMonopbmc))
exprMat<-stadMonopbmc[["RNA"]]$counts
exprMat[1:3,1:3]
sorted_column_names <- colnames(exprMat)
meta.data2 <- meta.data1[sorted_column_names, ]
head(meta.data2)
table(meta.data2$group5)
cds <- new_cell_data_set(exprMat,
                         cell_metadata = meta.data2,
                         gene_metadata = f_data_cd4)
class(exprMat)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = 'PCA')
#这里是umap的藏身之地!!!!!!!!!!!!
emb <- stadMonopbmc@reductions[["umap"]]@cell.embeddings
SingleCellExperiment::reducedDims(cds)[["UMAP"]]=emb
cds <- cluster_cells(cds, reduction_method ="UMAP")
#---------#
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition")

#----自定义root-----#
## 定义函数 get_earliest_principal_node ：根据指定的细胞类型（celltype）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin="Mono/Macro_low"){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "CELLTYPE5"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}
nodes_vec <- c(get_earliest_principal_node(cds,"Mono/Macro_low"))
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

pseudotime2 <- merge(pseudotime1,STADsingle,by.x="cell",by.y="Cell")
head(pseudotime2)
pseudotime2 <-pseudotime2[,c("cell","pseudotime","model1","model3","model5","celltype")]
rownames(pseudotime2) <- pseudotime2$cell
pseudotime2$cell <- NULL
pseudotime2 <- pseudotime2[!is.infinite(pseudotime2$pseudotime),]
table(pseudotime2$celltype)
save(pseudotime2,file = "pseudotime2.RData")
#计算相关性
singlemt <- pseudotime2[,c("pseudotime","model5")]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson

pdf("model5cor.pdf")
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
dev.off()





#DC
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/stad/stadDCpbmc.RData")
exprMat<-stadDCpbmc[["RNA"]]$counts
exprMat <- exprMat[apply(exprMat > 0,1,sum) >= 5,]
# 读取注释文件
annotation_data <- STADsingle

duplicated_row_names <- duplicated(annotation_data$Cell) | annotation_data$Cell %in% duplicated(annotation_data$Cell)
# 保留第一个出现的行，删除其余重复的行
annotation_data <- annotation_data[!duplicated_row_names & !is.na(annotation_data$Cell), ]
# 更新行名
rownames(annotation_data) <- annotation_data$Cell
table(annotation_data$celltype)

#rownames(annotation_data)<-annotation_data$Cell
meta.data<-stadDCpbmc@meta.data
meta.data$cell<-rownames(meta.data)
#查看是否有相同的细胞类型
cell<-intersect(rownames(meta.data),annotation_data$Cell)
#提取相同的细胞类型
cell_inf1<-annotation_data[cell,]
cell_inf2<-data.frame(cell_inf1$Cell,cell_inf1$celltype,cell_inf1$model5,cell_inf1$`Cell Cycle`,cell_inf1$`DNA damage`,cell_inf1$`DNA repair`)
names(cell_inf2)<-c("cell","celltype","model5","Cell Cycle","DNA damage","DNA repair")
head(cell_inf2)
#提取相同细胞类型的meta.data文件
meta.data1<-merge(meta.data,cell_inf2,by.x = "cell",by.y = "cell")
rownames(meta.data1)<-meta.data1$cell
meta.data1[1:3,]
table(meta.data1$celltype)
meta.data1$group5 <- ifelse(meta.data1$model5>median(meta.data1$model5),"high","low")
meta.data1$CELLTYPE <- paste(meta.data1$celltype,"_",meta.data1$group5,sep = "")
table(meta.data1$group5)
table(meta.data1$CELLTYPE)

exprMat2 <- data.frame(exprMat)
#-------CytoTRACE---------#
results <- CytoTRACE(exprMat2)
emb <- stadDCpbmc@reductions[["umap"]]@cell.embeddings
#rownames(emb)<-gsub("[^a-zA-Z0-9]", ".",rownames(emb))
rownames(emb)<-gsub("[@-]", ".",rownames(emb))

pheno <- meta.data1$CELLTYPE
pheno <- as.character(pheno)
names(pheno) <- rownames(meta.data1)
#names(pheno) <- gsub("[^a-zA-Z0-9]", ".",names(pheno))
names(pheno) <- gsub("[@-]", ".",names(pheno))

setdiff(names(pheno) ,colnames(exprMat2))
#------结果--------#
setwd("/wangzhenzenlab/fsh/xff/fenhua/stad/DC")
plotCytoTRACE(results, phenotype = pheno, emb = emb)
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

#----monocle3-----#
f_data_cd4 <- data.frame(gene_short_name = row.names(stadDCpbmc) ,row.names = row.names(stadDCpbmc))
exprMat<-stadDCpbmc[["RNA"]]$counts
exprMat[1:3,1:3]
sorted_column_names <- colnames(exprMat)
meta.data2 <- meta.data1[sorted_column_names, ]
head(meta.data2)
table(meta.data2$group5)
cds <- new_cell_data_set(exprMat,
                         cell_metadata = meta.data2,
                         gene_metadata = f_data_cd4)
class(exprMat)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = 'PCA')
#这里是umap的藏身之地!!!!!!!!!!!!
emb <- stadDCpbmc@reductions[["umap"]]@cell.embeddings
SingleCellExperiment::reducedDims(cds)[["UMAP"]]=emb
cds <- cluster_cells(cds, reduction_method ="UMAP")
#---------#
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition")

#----自定义root-----#
## 定义函数 get_earliest_principal_node ：根据指定的细胞类型（celltype）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin="DC_high"){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "CELLTYPE"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}
nodes_vec <- c(get_earliest_principal_node(cds,"DC_high"))
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

pseudotime2 <- merge(pseudotime1,STADsingle,by.x="cell",by.y="Cell")
head(pseudotime2)
pseudotime2 <-pseudotime2[,c("cell","pseudotime","model1","model3","model5","celltype")]
rownames(pseudotime2) <- pseudotime2$cell
pseudotime2$cell <- NULL
pseudotime2 <- pseudotime2[!is.infinite(pseudotime2$pseudotime),]
table(pseudotime2$celltype)
save(pseudotime2,file = "pseudotime2.RData")
#计算相关性
singlemt <- pseudotime2[,c("pseudotime","model5")]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson

pdf("model5cor.pdf")
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
dev.off()






#Mast
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/stad/stadMastpbmc.RData")
exprMat<-stadMastpbmc[["RNA"]]$counts
exprMat <- exprMat[apply(exprMat > 0,1,sum) >= 5,]
# 读取注释文件
annotation_data <- STADsingle

duplicated_row_names <- duplicated(annotation_data$Cell) | annotation_data$Cell %in% duplicated(annotation_data$Cell)
# 保留第一个出现的行，删除其余重复的行
annotation_data <- annotation_data[!duplicated_row_names & !is.na(annotation_data$Cell), ]
# 更新行名
rownames(annotation_data) <- annotation_data$Cell
table(annotation_data$celltype)

#rownames(annotation_data)<-annotation_data$Cell
meta.data<-stadMastpbmc@meta.data
meta.data$cell<-rownames(meta.data)
#查看是否有相同的细胞类型
cell<-intersect(rownames(meta.data),annotation_data$Cell)
#提取相同的细胞类型
cell_inf1<-annotation_data[cell,]
cell_inf2<-data.frame(cell_inf1$Cell,cell_inf1$celltype,cell_inf1$model5,cell_inf1$model1)
names(cell_inf2)<-c("cell","celltype","model5","model1")
head(cell_inf2)
#提取相同细胞类型的meta.data文件
meta.data1<-merge(meta.data,cell_inf2,by.x = "cell",by.y = "cell")
rownames(meta.data1)<-meta.data1$cell
meta.data1[1:3,]
table(meta.data1$celltype)
meta.data1$group5 <- ifelse(meta.data1$model5>median(meta.data1$model5),"high","low")
meta.data1$CELLTYPE <- paste(meta.data1$celltype,"_",meta.data1$group5,sep = "")
table(meta.data1$group5)
table(meta.data1$CELLTYPE)

exprMat2 <- data.frame(exprMat)
#-------CytoTRACE---------#
results <- CytoTRACE(exprMat2)
emb <- stadMastpbmc@reductions[["umap"]]@cell.embeddings
#rownames(emb)<-gsub("[^a-zA-Z0-9]", ".",rownames(emb))
rownames(emb)<-gsub("[@-]", ".",rownames(emb))

pheno <- meta.data1$CELLTYPE
pheno <- as.character(pheno)
names(pheno) <- rownames(meta.data1)
#names(pheno) <- gsub("[^a-zA-Z0-9]", ".",names(pheno))
names(pheno) <- gsub("[@-]", ".",names(pheno))

setdiff(names(pheno) ,colnames(exprMat2))
#------结果--------#
setwd("/wangzhenzenlab/fsh/xff/fenhua/stad/Mast")
plotCytoTRACE(results, phenotype = pheno, emb = emb)
plotCytoGenes(results, numOfGenes = 30, outputDir = './')

#----monocle3-----#
f_data_cd4 <- data.frame(gene_short_name = row.names(stadMastpbmc) ,row.names = row.names(stadMastpbmc))
exprMat<-stadMastpbmc[["RNA"]]$counts
exprMat[1:3,1:3]
sorted_column_names <- colnames(exprMat)
meta.data2 <- meta.data1[sorted_column_names, ]
head(meta.data2)
table(meta.data2$group5)
cds <- new_cell_data_set(exprMat,
                         cell_metadata = meta.data2,
                         gene_metadata = f_data_cd4)
class(exprMat)
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds,preprocess_method = 'PCA')
#这里是umap的藏身之地!!!!!!!!!!!!
emb <- stadMastpbmc@reductions[["umap"]]@cell.embeddings
SingleCellExperiment::reducedDims(cds)[["UMAP"]]=emb
cds <- cluster_cells(cds, reduction_method ="UMAP")
#---------#
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "partition")

#----自定义root-----#
## 定义函数 get_earliest_principal_node ：根据指定的细胞类型（celltype）名称，找到Seurat对象中对应细胞的主成分图（UMAP）上的最近顶点
get_earliest_principal_node <- function(cds, time_bin="Mast_high"){
  # 首先找到指定celltype的 ID
  cell_ids <- which(colData(cds)[, "CELLTYPE"] == time_bin)
  # 获取主成分图（UMAP）中细胞投影到最近分支结点的信息
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 找到在指定的celltype出现次数最多的分支节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))] # igraph::V()函数获取主成分图中的所有顶点信息
  root_pr_nodes
}
nodes_vec <- c(get_earliest_principal_node(cds,"Mast_high"))
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

pseudotime2 <- merge(pseudotime1,STADsingle,by.x="cell",by.y="Cell")
head(pseudotime2)
pseudotime2 <-pseudotime2[,c("cell","pseudotime","model1","model3","model5","celltype")]
rownames(pseudotime2) <- pseudotime2$cell
pseudotime2$cell <- NULL
pseudotime2 <- pseudotime2[!is.infinite(pseudotime2$pseudotime),]
table(pseudotime2$celltype)
save(pseudotime2,file = "pseudotime2.RData")
#计算相关性
singlemt <- pseudotime2[,c("pseudotime","model5")]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson

pdf("model5cor.pdf")
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
dev.off()

























#ucec
load("D:/细胞状态2/论文/fig3/fig3d/ucec/ucecCD4Tpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/ucec/ucecCD8Tpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/ucec/ucecDCpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/ucec/ucecMastpbmc.RData")
load("D:/细胞状态2/论文/fig3/fig3d/ucec/ucecMonopbmc.RData")



exprMat<-ucecCD4Tpbmc[["RNA"]]$counts
exprMat<-ucecCD8Tpbmc[["RNA"]]$counts
exprMat<-ucecDCpbmc[["RNA"]]$counts
exprMat<-ucecMastpbmc[["RNA"]]$counts
exprMat<-ucecMonopbmc[["RNA"]]$counts


getwd()
setwd("D:/细胞状态2/论文/fig3/fig3d/STAD")
head(pseudotime2)
table(pseudotime2$celltype)
singlemt <- pseudotime2[,c("pseudotime","model6")]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson


#model值umap图
stadEndumap <- stadEndpbmc@reductions[["umap"]]@cell.embeddings
stadFibumap <- stadFibpbmc@reductions[["umap"]]@cell.embeddings

class(stadEndumap)
stadEndumap <- as.data.frame(stadEndumap)
stadFibumap <- as.data.frame(stadFibumap)

stadEndumap$cell <- rownames(stadEndumap)
stadFibumap$cell <- rownames(stadFibumap)

stadEndumap2 <- merge(stadEndumap,STADsingle,by.x="cell",by.y="Cell")#1
stadFibumap2 <- merge(stadFibumap,STADsingle,by.x="cell",by.y="Cell")#5


library(ggplot2)
max(stadFibumap2$model5)
min(stadFibumap2$model5)
ggplot(stadFibumap2, aes(x = umap_1, y = umap_2, color = model5)) +
  geom_point(size = 3,alpha = 1) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_gradient2(midpoint = -0.15, low = "#f0f623", mid = "#cd4a75", high = "#040081", 
                        breaks = c(-0.3, -0.15, 0), labels = c("-0.3", "-0.15", "0"))













#绘制umap
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/skcm/skcmCD8Tpbmc.RData")
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/stad/stadEndpbmc.RData")
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/stad/stadFibpbmc.RData")

skcmCD8Tumap <- skcmCD8Tpbmc@reductions$umap@cell.embeddings
stadEndumap <- stadEndpbmc@reductions$umap@cell.embeddings
stadFib <- stadFibpbmc@reductions$umap@cell.embeddings

head(skcmCD8Tumap)
head(SKCMsingle)
skcmCD8Tumap <- as.data.frame(skcmCD8Tumap)
skcmCD8Tumap$Cell <- rownames(skcmCD8Tumap)
skcmCD8Tumap2 <- merge(skcmCD8Tumap,SKCMsingle,by="Cell")
head(skcmCD8Tumap2)
max(skcmCD8Tumap2$model4)
min(skcmCD8Tumap2$model4)

pdf("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fenhua/skcm/skcmCD8Tumap2.pdf")
ggplot(skcmCD8Tumap2, aes(x = umap_1, y = umap_2, color = model5)) +
  geom_point(size = 0.1,alpha = 1) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_gradient2(midpoint = 0, low = "#f0f623", mid = "#cd4a75", high = "#040081", 
                        breaks = c(-0.45, 0, 0.45), labels = c("-0.45", "0", "0.45"))
dev.off()



head(stadEndumap)
head(STADsingle)
stadEndumap <- as.data.frame(stadEndumap)
stadEndumap$Cell <- rownames(stadEndumap)
stadEndumap2 <- merge(stadEndumap,STADsingle,by="Cell")
head(stadEndumap2)
max(stadEndumap2$model1)
min(stadEndumap2$model1)
pdf("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fenhua/stad/stadEndumap2.pdf")
ggplot(stadEndumap2, aes(x = umap_1, y = umap_2, color = model1)) +
  geom_point(size = 3,alpha = 1) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_gradient2(midpoint = 0.15, low = "#f0f623", mid = "#cd4a75", high = "#040081", 
                        breaks = c(0, 0.15, 0.35), labels = c("0", "0.15", "0.35"))
dev.off()




head(stadFib)
head(STADsingle)
stadFib <- as.data.frame(stadFib)
stadFib$Cell <- rownames(stadFib)
stadFib2 <- merge(stadFib,STADsingle,by="Cell")
head(stadFib2)
max(stadFib2$model5)
min(stadFib2$model5)
pdf("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fenhua/stad/stadFib2.pdf")
ggplot(stadFib2, aes(x = umap_1, y = umap_2, color = model5)) +
  geom_point(size = 3,alpha = 1) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_gradient2(midpoint = -0.15, low = "#f0f623", mid = "#cd4a75", high = "#040081", 
                        breaks = c(-0.3, -0.15, 0), labels = c("-0.3", "-0.15", "0"))
dev.off()











