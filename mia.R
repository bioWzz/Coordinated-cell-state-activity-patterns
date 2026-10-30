load("D:/细胞状态2/论文/空间/st2/seurat_data.RData")


cluster_gene_stat=FindAllMarkers(seurat_data,only.pos = TRUE,logfc.threshold = 0.9)
head(cluster_gene_stat)
table(cluster_gene_stat$cluster)
st.clusts=paste0(seq(0,10,1))
head(cluster_gene_stat)

st.marker.list=split(cluster_gene_stat,f = cluster_gene_stat$cluster)
st.marker.list = lapply(st.marker.list,rownames)
head(st.marker.list)

library(readxl)
markgene <- read_excel("D:/细胞状态2/论文/空间/markgene.xlsx")
#markgene <- read.csv("D:/细胞状态2/论文/空间/markgene.csv",header = T)
markgene <- as.data.frame(markgene)
rownames(markgene) <- markgene$markgene

sc.clusts <- markgene$celltype
sc.clusts <- unique(sc.clusts)
cellType_marker <- split(markgene,f=markgene$celltype)
cellType_marker = lapply(cellType_marker,rownames)
sc.marker.list <- cellType_marker


N <- length(st.clusts) 
M <- length(sc.clusts)
MIA.results <- matrix(0,nrow = M, ncol = N)
row.names(MIA.results) <- sc.clusts
colnames(MIA.results) <- st.clusts
head(MIA.results)



gene.universe <- length(rownames(cluster_gene_stat))



# Loop over ST clusters
for (i in 1:N) {
  # Then loop over SC clusters
  # i=1
  for (j in 1:M) {
    #j <- 1
    genes1 <- st.marker.list[[st.clusts[i]]]
    genes2 <- sc.marker.list[[sc.clusts[j]]]
    
    # Hypergeometric    
    A <- length(intersect(genes1,genes2))
    B <- length(genes1)
    C <- length(genes2)
    enr <- -log10(phyper(A, B, gene.universe-B, C, lower.tail = FALSE))
    dep <- -log10(1-phyper(A, B, gene.universe-B, C, lower.tail = FALSE))
    if (enr < dep) {
      MIA.results[j,i] = -dep
    } else {
      MIA.results[j,i] = enr
    }
  }
}


# Some results were -Inf...check why this is the case...
MIA.results[is.infinite(MIA.results)] <- 0
MIA.results


pheatmap::pheatmap(MIA.results,scale ='row' )


MIA.results2=MIA.results %>%as.matrix()

MIA.results2[ MIA.results >5 ] =5


pheatmap::pheatmap(MIA.results2,scale ='row' )


# Visualize as heatmap
library(reshape2)
library(ggplot2)
library(dplyr)
library(Seurat)
heatmap_df <- data.frame('Cell types' = melt(MIA.results)[,1],
                         'Tissue regions' = melt(MIA.results)[,2],
                         enrichment = melt(MIA.results)[,3])
range(heatmap_df$enrichment)

head(heatmap_df)

boxplot(heatmap_df$enrichment)

#根据具体数值，设置最大值和最小值
heatmap_df[abs(heatmap_df$enrichment) >7,]$enrichment=7


ggplot(data = heatmap_df, aes(x = Tissue.regions, y = Cell.types, fill = enrichment)) +
  geom_tile() + 
  scale_fill_gradient2(low = "navyblue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-7,7), space = "Lab", 
                       name="Enrichment \n -log10(p)") +
  ylim(heatmap_df$Cell.types %>% levels() %>% sort() %>% rev())+
  theme_minimal()+RotatedAxis() 


library(reshape2)
library(ggplot2)
library(dplyr)
library(Seurat)
library(scales)
heatmap_df <- data.frame('Cell types' = melt(MIA.results)[,1],
                         'Tissue regions' = melt(MIA.results)[,2],
                         enrichment = melt(MIA.results)[,3])
head(heatmap_df)

ggplot(data = heatmap_df, aes(x = Tissue.regions, y = Cell.types, fill = enrichment)) +
  geom_tile() + 
  scale_fill_gradient2(low = "blue", high = "red",mid = 'white', 
                       midpoint = 0, limit = c(-1,1), space = "Lab",oob=squish,
                       name="Enrichment \n -log10(p)") +
  ylim(heatmap_df$Cell.types %>% levels() %>% sort() %>% rev())+
  theme_minimal()+RotatedAxis() 
















#fsh
st.clusts <- Idents(test_data) %>% levels()
N <- length(st.clusts)

st.marker.list <- vector(mode = 'list', length = N)
names(st.marker.list) <- st.clusts
for(i in st.clusts) {
  st.marker.list[[i]] <- region_marker[region_marker$cluster == i,'gene']
}

# Create a list object containing the marker genes for each cell type:
sc.clusts <- Idents(seurat_data) %>% levels()
M <- length(sc.clusts)

sc.marker.list <- vector(mode = 'list', length = M)
names(sc.marker.list) <- sc.clusts

for (i in sc.clusts) {
  sc.marker.list[[i]] <- maintype_marker[maintype_marker$cluster == i,'gene']
}
#### Let's first demonstrate MIA with a single example: fibroblasts (scRNA-seq) and the cancer region (ST):
# Initialize a dataframe for us to store values in:
N <- length(st.clusts)
M <- length(sc.clusts)
MIA.results <- matrix(0,nrow = M, ncol = N)
row.names(MIA.results) <- sc.clusts
colnames(MIA.results) <- st.clusts

# Gene universe
gene.universe <- length(rownames(test_data))

# Loop over ST clusters
for (i in 1:N) {
  # Then loop over SC clusters
  for (j in 1:M) {
    genes1 <- st.marker.list[[st.clusts[i]]]
    genes2 <- sc.marker.list[[sc.clusts[j]]]
    
    # Hypergeometric    
    A <- length(intersect(genes1,genes2))
    B <- length(genes1)
    C <- length(genes2)
    enr <- -log10(phyper(A, B, gene.universe-B, C, lower.tail = FALSE))
    dep <- -log10(1-phyper(A, B, gene.universe-B, C, lower.tail = FALSE))
    if (enr < dep) {
      MIA.results[j,i] = -dep
    } else {
      MIA.results[j,i] = enr
    }
  }
}

print(MIA.results)
#----------------#
# Visualize as heatmap
library(dplyr)
heatmap_df <- data.frame('Cell types' = melt(MIA.results)[,1],
                         'Tissue regions' = melt(MIA.results)[,2],
                         enrichment = melt(MIA.results)[,3])

ggplot(data = heatmap_df, aes(x = Tissue.regions, y = Cell.types, fill = enrichment)) +
  geom_tile() + 
  scale_fill_gradient2(low = "navyblue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-10,10), space = "Lab", name="Enrichment \n -log10(p)") +
  ylim(heatmap_df$Cell.types %>% levels() %>% sort() %>% rev()) + 
  theme_minimal()
#------------#
meta.data<-test_data@meta.data
ann_col<-unique(meta.data[,c(8,9)])
names(ann_col)<-c("SB","region")
rownames(ann_col)<-ann_col$SB

# 定义属于 "Moderate" 组的列名
moderate_groups <- ann_col[ann_col$region=="Cancer",]
moderate_groups<-moderate_groups$SB

# 如果列名属于 "Moderate" 组，则排在前，不属于则排在后面
is_M_group <- MIA.results %>% colnames() %in% moderate_groups
df_sorted <- MIA.results1[, c(which(is_M_group), which(!is_M_group))]


df_sorted<-scale(df_sorted)

ann_color = list(region =c("stromal"="#fb776e","Cancer"="#00c5ca")) #定义分组颜色


MIA.results1<-MIA.results
MIA.results1<-MIA.results1[-2,]
MIA.results1[MIA.results1 == Inf] <- 0
MIA.results1<-MIA.results1+1
MIA.results1<-data.frame(MIA.results1)
names(MIA.results)<- gsub("X", "", names(MIA.results))

pheatmap(df_sorted, scale = "column", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",#表示行、列聚类使用的度量方法，默认为“euclidean”，也可为  "correlation"即按照 Pearson correlation方法进行聚类
         show_rownames = T,
         #cluster_rows = T,
         cellwidth = 15, cellheight = 15,
         clustering_method = "complete", #表示聚类方法，包括：‘ward’, ‘ward.D’, ‘ward.D2’, ‘single’, ‘complete’, ‘average’, ‘mcquitty’, ‘median’, ‘centroid’
         cluster_rows = T,cluster_cols =F , #cluster_rows表示仅对行聚类，cluster_cols表示仅对列聚类，值为TRUE或FALSE
         # cutree_rows = 2, #cutree_cols = , #若进行了行、列聚类，根据行、列聚类数量分隔热图行
         treeheight_row = 30, treeheight_col = 30,#、行、列聚类树高度调整
         annotation_row =NA , annotation_col =ann_col , #表示是否对行、列进行注释，默认NA
         annotation = NA, annotation_colors = ann_color,  #表示行注释及列注释的颜色，默认NA
         color = colorRampPalette(c('#266EAF',"white",'#B52230'))(100), #表示热图颜色,(100)表示100个等级
) 
dev.off()
# 对行和列都不聚类，自定义划分行和列的gap
pheatmap(MIA.results1, annotation_col = ann_col,
         cluster_rows = FALSE, cluster_cols = FALSE, 
         gaps_row = c(1, 2), gaps_col = c(5),
         cellwidth = 20, cellheight = 15)




