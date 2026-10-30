#RNA速率
load("/wangzhenzenlab/fsh/xff/fenhua/paad/paadmyoFibpbmc.RData")


install.packages("remotes")
remotes::install_github("velocyto-team/velocyto.R")
Sys.which("make")
Sys.getenv("PATH")


exprMat<-paadmyoFibpbmc[["RNA"]]$counts
exprMat[1:3,1:3]

umap_coords <- Embeddings(paadmyoFibpbmc, reduction = "umap")
head(umap_coords)
umap_coords <- as.data.frame(umap_coords)
umap_coords$cell <- rownames(umap_coords)

#seurat_clusters
head(paadmyoFibpbmc@meta.data)
annotation1<- paadmyoFibpbmc@meta.data
#annotation1 <- as.data.frame(annotation1)
annotation1$cell <- rownames(annotation1)
head(annotation1)
annotation1 <- annotation1[,c("seurat_clusters","cell")]

data2 <- merge(annotation1,umap_coords,by="cell")
head(data2)

write.csv(exprMat,file = "/wangzhenzenlab/fsh/xff/fenhua/RNA/PAAD/exprMat.csv")
write.csv(data2,file = "/wangzhenzenlab/fsh/xff/fenhua/RNA/PAAD/position.csv")

data2 <- read.csv("/wangzhenzenlab/fsh/xff/fenhua/RNA/PAAD/position.csv",header = T)







load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/stad/stadEndpbmc.RData")
UMAPPlot(stadEndpbmc)

exprMat<-stadEndpbmc[["RNA"]]$counts
exprMat[1:3,1:3]

umap_coords <- Embeddings(stadEndpbmc, reduction = "umap")
head(umap_coords)
umap_coords <- as.data.frame(umap_coords)
umap_coords$cell <- rownames(umap_coords)

#seurat_clusters
head(stadEndpbmc@meta.data)
annotation1<- stadEndpbmc@meta.data
#annotation1 <- as.data.frame(annotation1)
annotation1$cell <- rownames(annotation1)
head(annotation1)
annotation1 <- annotation1[,c("seurat_clusters","cell")]

data2 <- merge(annotation1,umap_coords,by="cell")
head(data2)

write.csv(exprMat,file = "/wangzhenzenlab/fsh/xff/fenhua/RNA/STAD/exprMat.csv")
write.csv(data2,file = "/wangzhenzenlab/fsh/xff/fenhua/RNA/STAD/position.csv")










#TP高低
library(ggplot2)
head(data2)
head(PAADsingle)
dim(data3)
data3 <- merge(data2,PAADsingle,by.x="cell",by.y="Cell")
head(data3)
data3$group1 <- ifelse(data3$model1>median(data3$model1),"high","low")
table(data3$group1)
pdf("/wangzhenzenlab/fsh/xff/fenhua/RNA/PAAD/highlow.pdf")
ggplot(data3, aes(x = UMAP_1.x , y = UMAP_2.x , color = group1)) +
  geom_point(size = 2,alpha = 1) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_manual(values=c("low" = "steelblue", "high" = "#c11e21"))
dev.off()
