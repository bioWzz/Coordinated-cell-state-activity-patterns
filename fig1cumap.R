#服务器
load("/wangzhenzenlab/xff/cellstate/xff/cellstate/cellstatemergedata.RData")
load("/wangzhenzenlab/xff/cellstate/xff/randomcellumap.RData")

head(randomcellumap)
class(randomcellumap)
randomcellumap <- as.data.frame(randomcellumap)
randomcellumap$Cell <- rownames(randomcellumap)
randomcellumap$Cell <- gsub("_", ".", randomcellumap$Cell)

head(cellstatemergedata)
cellstatemergedata[,17:18] <- NULL
randomcellumap2 <- merge(randomcellumap,cellstatemergedata,by.x = "Cell",by.y = "cell")
head(randomcellumap2)
randomcellumap3 <- randomcellumap2[,c(2,3,20)]
head(randomcellumap3)
table(randomcellumap3$`Celltype (major-lineage)`)

library(ggplot2)
pdf("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig1/randomcellumap.pdf")
ggplot(randomcellumap3, aes(x = UMAP_1, y = UMAP_2, color = `Celltype (major-lineage)`)) +
  geom_point(size = 0.1,alpha = 0.5) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("B"="#F8766C","CD4Tconv"="#c11e21","CD8T"="#E48600",
                                "CD8Tex"="#E7B781","DC"="#44bd28","Endothelial"="steelblue",
                                "Epithelial"="#B883FF","Fibroblasts"="#A3B5D6","Malignant"="#3C6882",
                                "Mast"="#BFD3E6","Mono/Macro"="#70C1B3","Myofibroblasts"="#73B3BE",
                                "NK"="#6D9EC3","Tprolif"="#c6c322","Treg"="#A7C66B"))


dev.off()

adjusted_colors <- c("#FB8072", "#c11e21", "#E48600", "#E7B781", "#44bd28", 
                     "#A7C66B", "#96C18C", "#BFD3E6", "#70C1B3", "#6D9EC3",
                     "#3C6882", "steelblue", "#B883FF", "#A3B5D6", "#73B3BE")





table(randomcellumap2$cancertype)

randomcellumap4 <- randomcellumap2[randomcellumap2$cancertype%in%c("BCC","BRCA","CESC",
                                                                   "CRC","ESCA","KIRC",
                                                                   "LIHC","NSCLC","OV",
                                                                   "PAAD","PRAD","SKCM",
                                                                   "STAD","THCA","UCEC"),]

library(ggplot2)
pdf("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig1/randomcellumapcancer.pdf")
ggplot(randomcellumap4, aes(x = UMAP_1, y = UMAP_2, color = cancertype)) +
  geom_point(size = 0.1,alpha = 0.5) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("BCC"="#F8766C","BRCA"="#c11e21","CESC"="#E48600",
                                "CRC"="#E7B781","ESCA"="#44bd28","KIRC"="steelblue",
                                "LIHC"="#B883FF","NSCLC"="#A3B5D6","OV"="#3C6882",
                                "PAAD"="#BFD3E6","PRAD"="#70C1B3","SKCM"="#73B3BE",
                                "STAD"="#6D9EC3","THCA"="#c6c322","UCEC"="#A7C66B"))


dev.off()


