library(ggplot2)
library(Seurat)
#umap图
load("D:/细胞状态2/论文/fig1/TUMOR3.RData")
head(TUMOR3)

#基因在umap上的表达
# 
# load("D:/细胞状态2/论文/fig4/fig4d/pbmcpaad.RData")
# head(pbmcpaad@meta.data)
# pbmcpaad[1:3,1:3]
# expr_matrix <- GetAssayData(object = pbmcpaad, assay = "RNA", slot = "data")
# expr_matrix[1:3,1:3]
# class(expr_matrix)
# expr_matrix <- as.data.frame(expr_matrix)
# head(EMTgene)
# emtmatrix <- expr_matrix[EMTgene$GeneName,]
# emtmatrix <- na.omit(emtmatrix)
# paadexpmatrix1<- t(emtmatrix)
# paadexpmatrix1 <- as.data.frame(paadexpmatrix1)
# paadexpmatrix1$cell <- rownames(paadexpmatrix1)
# 
# paadexpmatrix2 <- merge(umap3dpaadmergedata,paadexpmatrix1,by="cell")
# head(paadexpmatrix2)
# 
# ggplot(paadexpmatrix2, aes(x = umap_1, y = umap_2, color = COL1A1)) +
#   geom_point(size = 0.1,alpha = 0.5) +  # 绘制散点
#   labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
#   theme_minimal()+ 
#   theme(panel.grid = element_blank())+
#   scale_color_gradient2(low = "steelblue", mid = "white", high = "#c11e21", midpoint = 0, limits = c(-0.3, 0.3)) +
#   scale_fill_gradient2(low = "steelblue", mid = "white", high = "#c11e21", midpoint = 0, limits = c(-0.3, 0.3)) 



library(ggplot2)
library(ggprism)
library(ggpubr)
library(tidyr)

#多个基因在不同组的表达情况箱式图
# group <- levels(factor(paadexpmatrix2$group1))
# comp <- combn(group,2)
# my_comparisons <- list()
# for (i in 1:ncol(comp)) {
#   my_comparisons[[i]] <- comp[,i]
# }
# 
# paadexpmatrix3 <- paadexpmatrix2[,c(1,41,48:97)]
# head(paadexpmatrix3)
# paadexpmatrix3_long <- pivot_longer(paadexpmatrix3, cols = -c("cell","group1"), names_to = "gene", values_to = "exp")
# head(paadexpmatrix3_long)
# 
# table(paadexpmatrix2$celltype)
# emtpaadexpmatrix <- paadexpmatrix2[paadexpmatrix2$celltype=="Fibroblasts",]
# emtpaadexpmatrix2 <- emtpaadexpmatrix[,c(1,41,48:97)]
# head(emtpaadexpmatrix2)
# emtpaadexpmatrix2_long <- pivot_longer(emtpaadexpmatrix2, cols = -c("cell","group1"), names_to = "gene", values_to = "exp")
# head(emtpaadexpmatrix2_long)
# 
# ggplot(data = emtpaadexpmatrix2_long, aes(x = gene, y = exp, fill = group1)) +
#   geom_violin(position = position_dodge(0.8))+
#   geom_boxplot(width=0.5,
#                #outlier.colour = NA,
#                position = position_dodge(0.8))+  #每组的小提琴图之间的距离
#   stat_compare_means(aes(group = group1))+  #,label = 'p.signif'
#   scale_fill_manual(values = c("#c11e21", "steelblue"))+
#   theme_prism()+
#   theme(axis.text.x = element_text(size = 20,angle = 30, hjust = 0.5, vjust = 0.5))





load("D:/细胞状态2/论文/fig3/fig3d/umap3dpaad.RData")
head(umap3dpaad)
umap3dpaad <- as.data.frame(umap3dpaad)
umap3dpaad$cell <- rownames(umap3dpaad)
umap3dpaadmergedata <- merge(umap3dpaad,TUMOR3,by.x = "cell",by.y = "Cell")
head(umap3dpaadmergedata)
table(umap3dpaadmergedata$celltype)

umap3dpaadmergedata$group1 <- ifelse(umap3dpaadmergedata$model1>median(umap3dpaadmergedata$model1),"high","low")

ggplot(umap3dpaad, aes(x = umap_1, y = umap_2))+
  geom_point(size = 0.1,alpha = 0.5)
ggplot(umap3dpaadmergedata, aes(x = umap_1, y = umap_2, color = celltype)) +
  geom_point(size = 0.1,alpha = 0.5) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("B"="#F8766C","CD4Tconv"="#c11e21","CD8T"="#E48600",
                                "CD8Tex"="#E7B781","DC"="#44bd28","Endothelial"="steelblue",
                                "Epithelial"="#B883FF","Fibroblasts"="#A3B5D6","Malignant"="#3C6882",
                                "Mast"="#BFD3E6","Mono/Macro"="#70C1B3","Myofibroblasts"="#73B3BE",
                                "NK"="#6D9EC3","Tprolif"="#c6c322","Treg"="#A7C66B"))


ggplot(umap3dpaadmergedata, aes(x = umap_1, y = umap_2, color = group1)) +
  geom_point(size = 0.1,alpha = 0.5) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("high"="#c11e21","low"="steelblue"))


head(umap3dpaadmergedata)
max(umap3dpaadmergedata$EMT)
min(umap3dpaadmergedata$EMT)
ggplot(umap3dpaadmergedata, aes(x = umap_1, y = umap_2, color = Metastasis)) +
  geom_point(size = 0.1,alpha = 0.5) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_gradient2(low = "steelblue",  high = "#c11e21", midpoint = 0, limits = c(-0.4, 0.4),oob = scales::squish) +
  scale_fill_gradient2(low = "steelblue",  high = "#c11e21", midpoint = 0, limits = c(-0.4, 0.4),oob = scales::squish) 



#save(pbmcpaad,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/pbmcpaad.RData")
#save(umap3dpaad,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/umap3dpaad.RData")



#SKCM
load("D:/细胞状态2/论文/fig3/fig3d/umap3dskcm.RData")

head(umap3dskcm)
umap3dskcm <- as.data.frame(umap3dskcm)
umap3dskcm$cell <- rownames(umap3dskcm)
umap3dskcmmergedata <- merge(umap3dskcm,TUMOR3,by.x = "cell",by.y = "Cell")
head(umap3dskcmmergedata)
table(umap3dskcmmergedata$celltype)

umap3dskcmmergedata$group4 <- ifelse(umap3dskcmmergedata$model4>median(umap3dskcmmergedata$model4),"high","low")
umap3dskcmmergedata$group6 <- ifelse(umap3dskcmmergedata$model6>median(umap3dskcmmergedata$model6),"high","low")

ggplot(umap3dskcmmergedata, aes(x = umap_1, y = umap_2, color = celltype)) +
  geom_point(size = 0.1,alpha = 0.5) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("B"="#F8766C","CD4Tconv"="#c11e21","CD8T"="#E48600",
                                "CD8Tex"="#E7B781","DC"="#44bd28","Endothelial"="steelblue",
                                "Epithelial"="#B883FF","Fibroblasts"="#A3B5D6","Malignant"="#3C6882",
                                "Mast"="#BFD3E6","Mono/Macro"="#70C1B3","Myofibroblasts"="#73B3BE",
                                "NK"="#6D9EC3","Tprolif"="#c6c322","Treg"="#A7C66B"))

ggplot(umap3dskcmmergedata, aes(x = umap_1, y = umap_2, color = group6)) +
  geom_point(size = 0.1,alpha = 0.5) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("high"="#c11e21","low"="steelblue"))

head(umap3dskcmmergedata)
max(umap3dskcmmergedata$Proliferation)
min(umap3dskcmmergedata$Proliferation)
ggplot(umap3dskcmmergedata, aes(x = umap_1, y = umap_2, color = Quiescence)) +
  geom_point(size = 0.1,alpha = 0.5) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_gradient2(low = "steelblue", mid = "white", high = "#c11e21", midpoint = 0, limits = c(-0.4, 0.4),oob = scales::squish) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "#c11e21", midpoint = 0, limits = c(-0.4, 0.4),oob = scales::squish) 


#STAD
load("D:/细胞状态2/论文/fig3/fig3d/umap3dstad.RData")

head(umap3dstad)
umap3dstad <- as.data.frame(umap3dstad)
umap3dstad$cell <- rownames(umap3dstad)
umap3dstadmergedata <- merge(umap3dstad,TUMOR3,by.x = "cell",by.y = "Cell")
head(umap3dstadmergedata)
table(umap3dstadmergedata$celltype)

umap3dstadmergedata$group1 <- ifelse(umap3dstadmergedata$model1>median(umap3dstadmergedata$model1),"high","low")
umap3dstadmergedata$group3 <- ifelse(umap3dstadmergedata$model3>median(umap3dstadmergedata$model3),"high","low")
umap3dstadmergedata$group5 <- ifelse(umap3dstadmergedata$model5>median(umap3dstadmergedata$model5),"high","low")


ggplot(umap3dstadmergedata, aes(x = umap_1, y = umap_2, color = celltype)) +
  geom_point(size = 0.1,alpha = 0.5) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("B"="#F8766C","CD4Tconv"="#c11e21","CD8T"="#E48600",
                                "CD8Tex"="#E7B781","DC"="#44bd28","Endothelial"="steelblue",
                                "Epithelial"="#B883FF","Fibroblasts"="#A3B5D6","Malignant"="#3C6882",
                                "Mast"="#BFD3E6","Mono/Macro"="#70C1B3","Myofibroblasts"="#73B3BE",
                                "NK"="#6D9EC3","Tprolif"="#c6c322","Treg"="#A7C66B"))


ggplot(umap3dstadmergedata, aes(x = umap_1, y = umap_2, color = group5)) +
  geom_point(size = 0.1,alpha = 0.5) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("high"="#c11e21","low"="steelblue"))



head(umap3dstadmergedata)
max(umap3dstadmergedata$`Cell Cycle`)
min(umap3dstadmergedata$`DNA damage`)
max(umap3dstadmergedata$`DNA repair`)
min(umap3dstadmergedata$`DNA repair`)

ggplot(umap3dstadmergedata, aes(x = umap_1, y = umap_2, color = Metastasis)) +
  geom_point(size = 0.5,alpha = 0.7) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_gradient2(low = "steelblue", mid = "white", high = "#c11e21", midpoint = 0, limits = c(-0.4, 0.4),oob = scales::squish) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "#c11e21", midpoint = 0, limits = c(-0.4, 0.4),oob = scales::squish) 




#UCEC
load("D:/细胞状态2/论文/fig3/fig3d/umap3ducec.RData")

head(umap3ducec)
umap3ducec <- as.data.frame(umap3ducec)
umap3ducec$cell <- rownames(umap3ducec)
umap3ducecmergedata <- merge(umap3ducec,TUMOR3,by.x = "cell",by.y = "Cell")
head(umap3ducecmergedata)
table(umap3ducecmergedata$celltype)

umap3ducecmergedata$TPgroup <- ifelse(umap3ducecmergedata$TP>median(umap3ducecmergedata$TP),"high","low")


ggplot(umap3ducecmergedata, aes(x = umap_1, y = umap_2, color = celltype)) +
  geom_point(size = 0.1,alpha = 0.5) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("B"="#F8766C","CD4Tconv"="#c11e21","CD8T"="#E48600",
                                "CD8Tex"="#E7B781","DC"="#44bd28","Endothelial"="steelblue",
                                "Epithelial"="#B883FF","Fibroblasts"="#A3B5D6","Malignant"="#3C6882",
                                "Mast"="#BFD3E6","Mono/Macro"="#70C1B3","Myofibroblasts"="#73B3BE",
                                "NK"="#6D9EC3","Tprolif"="#c6c322","Treg"="#A7C66B"))


ggplot(umap3ducecmergedata, aes(x = umap_1, y = umap_2, color = TPgroup)) +
  geom_point(size = 0.1,alpha = 0.5) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_manual(values = c("high"="#c11e21","low"="steelblue"))



head(umap3ducecmergedata)
max(umap3ducecmergedata$Invasion)
median(umap3ducecmergedata$Invasion)


#Angiogenesis   EMT     Invasion   Metastasis
ggplot(umap3ducecmergedata, aes(x = umap_1, y = umap_2, color = Angiogenesis)) +
  geom_point(size = 0.5,alpha = 0.7) +  # 绘制散点
  labs(title = "UMAP Plot", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()+ 
  theme(panel.grid = element_blank())+
  scale_color_gradient2(low = "steelblue", mid = "white", high = "#c11e21", midpoint = 0, limits = c(-0.4, 0.4),oob = scales::squish) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "#c11e21", midpoint = 0, limits = c(-0.4, 0.4),oob = scales::squish) 





library(Seurat)
library(hdf5r)

CESCsingle <- TUMOR3[TUMOR3$cancertype=="CESC",]
ESCAsingle <- TUMOR3[TUMOR3$cancertype=="ESCA",]
KIRCsingle <- TUMOR3[TUMOR3$cancertype=="KIRC",]
PAADsingle <- TUMOR3[TUMOR3$cancertype=="PAAD",]
SKCMsingle <- TUMOR3[TUMOR3$cancertype=="SKCM",]
STADsingle <- TUMOR3[TUMOR3$cancertype=="STAD",]
UCECsingle <- TUMOR3[TUMOR3$cancertype=="UCEC",]

head(CESCsingle)
ESCAsingle
KIRCsingle
PAADsingle
SKCMsingle
STADsingle
UCECsingle


#读取数据
#PAAD
paadexp1 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/PAAD_CRA001160/PAAD_CRA001160_expression.h5")
paadexp2 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/PAAD_GSE111672/PAAD_GSE111672_expression.h5")
paadexp3 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/PAAD_GSE141017/PAAD_GSE141017_expression.h5")
paadexp4 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/PAAD_GSE148673/PAAD_GSE148673_expression.h5")
paadexp5 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/PAAD_GSE154763/PAAD_GSE154763_expression.h5")
paadexp6 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/PAAD_GSE154778/PAAD_GSE154778_expression.h5")
paadexp7 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/PAAD_GSE158356/PAAD_GSE158356_expression.h5")
paadexp8 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/PAAD_GSE162708/PAAD_GSE162708_expression.h5")
paadexp9 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/PAAD_GSE165399/PAAD_GSE165399_expression.h5")


paadinter <- Reduce(intersect,list(rownames(paadexp1),rownames(paadexp2),rownames(paadexp3),rownames(paadexp4),rownames(paadexp5),rownames(paadexp6),rownames(paadexp7),rownames(paadexp8),rownames(paadexp9)))
paadinter <- unique(paadinter)
PAADexp <- cbind(paadexp1[paadinter,],paadexp2[paadinter,],paadexp3[paadinter,],paadexp4[paadinter,],paadexp5[paadinter,],paadexp6[paadinter,],paadexp7[paadinter,],paadexp8[paadinter,],paadexp9[paadinter,])

#SKCM
skcmexp1 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/SKCM_GSE72056/SKCM_GSE72056_expression.h5")
skcmexp2 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/SKCM_GSE115978_aPD1/SKCM_GSE115978_aPD1_expression.h5")
skcmexp3 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/SKCM_GSE120575_aPD1aCTLA4/SKCM_GSE120575_aPD1aCTLA4_expression.h5")
skcmexp4 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/SKCM_GSE123139/SKCM_GSE123139_expression.h5")
skcmexp5 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/SKCM_GSE134388_aPD1/SKCM_GSE134388_aPD1_expression.h5")
skcmexp6 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/SKCM_GSE139249/SKCM_GSE139249_expression.h5")
skcmexp7 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/SKCM_GSE148190/SKCM_GSE148190_expression.h5")
skcmexp8 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/SKCM_GSE159251/SKCM_GSE159251_expression.h5")
skcmexp9 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/SKCM_GSE166181_aPD1/SKCM_GSE166181_aPD1_expression.h5")
skcmexp10 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/SKCM_GSE179373/SKCM_GSE179373_expression.h5")


skcminter <- Reduce(intersect,list(rownames(skcmexp1),rownames(skcmexp2),rownames(skcmexp3),rownames(skcmexp4),rownames(skcmexp5),rownames(skcmexp6),rownames(skcmexp7),rownames(skcmexp8),rownames(skcmexp9),rownames(skcmexp10)))
skcminter <- unique(skcminter)
SKCMexp <- cbind(skcmexp1[skcminter,],skcmexp2[skcminter,],skcmexp3[skcminter,],skcmexp4[skcminter,],skcmexp5[skcminter,],skcmexp6[skcminter,],skcmexp7[skcminter,],skcmexp8[skcminter,],skcmexp9[skcminter,],skcmexp10[skcminter,])

#STAD
stadexp1 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/STAD_GSE134520/STAD_GSE134520_expression.h5")
stadexp2 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/STAD_GSE167297/STAD_GSE167297_expression.h5")


stadinter <- Reduce(intersect,list(rownames(stadexp1),rownames(stadexp2)))
stadinter <- unique(stadinter)
STADexp <- cbind(stadexp1[stadinter,],stadexp2[stadinter,])




#差异基因分析

library(limma)


#PAAD
PAADsingle$Fibroblasts <- ifelse(PAADsingle$celltype=="Fibroblasts","yes","no")
PAADsingle$Myofibroblasts <- ifelse(PAADsingle$celltype=="Myofibroblasts","yes","no")
PAADsingle$group1 <- ifelse(PAADsingle$model1>median(PAADsingle$model1),"high","low")

table(PAADsingle$Fibroblasts)
table(PAADsingle$celltype)


head(PAADsingle)
table(PAADsingle$group1)
PAADcoldata <- PAADsingle[,c("Cell","Fibroblasts","Myofibroblasts","group1")]
save(PAADcoldata,file = "D:/细胞状态2/论文/fig3/fig3d/PAADcoldata.RData")
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/PAADcoldata.RData")

PAADexp <- PAADexp[,PAADcoldata$Cell]
rownames(PAADcoldata)<- PAADcoldata$Cell

#方法一：limma
design <- model.matrix(~ Fibroblasts, data = PAADcoldata)
fit <- lmFit(PAADexp, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, adjust="BH", number=Inf)
print(results)
paadFibresult <- results
save(paadFibresult,file="/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/paadFibresult.RData")


design <- model.matrix(~ Myofibroblasts, data = PAADcoldata)
fit <- lmFit(PAADexp, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, adjust="BH", number=Inf)
print(results)
paadMyoresult <- results
save(paadMyoresult,file="/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/paadMyoresult.RData")


design <- model.matrix(~ group1, data = PAADcoldata)
fit <- lmFit(PAADexp, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, adjust="BH", number=Inf)
print(results)
paadgroup1result <- results
save(paadgroup1result,file="/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/paadgroup1result.RData")


#方法二：DESeq(需要count数据)
#dds <- DESeqDataSetFromMatrix(countData = PAADexp, colData = PAADcoldata, design = ~ Fibroblasts)
#dds <- DESeq(dds)

#SKCM
table(SKCMsingle$celltype)
SKCMsingle$Tprolif <- ifelse(SKCMsingle$celltype=="Tprolif","yes","no")
SKCMsingle$Malignant <- ifelse(SKCMsingle$celltype=="Malignant","yes","no")
SKCMsingle$Fibroblasts <- ifelse(SKCMsingle$celltype=="Fibroblasts","yes","no")
SKCMsingle$group4 <- ifelse(SKCMsingle$model4>median(SKCMsingle$model4),"high","low")
SKCMsingle$group6 <- ifelse(SKCMsingle$model6>median(SKCMsingle$model6),"high","low")

head(SKCMsingle)
SKCMcoldata <- SKCMsingle[,c("Cell","Tprolif","Malignant","Fibroblasts","group4","group6")]
save(SKCMcoldata,file = "D:/细胞状态2/论文/fig3/fig3d/SKCMcoldata.RData")
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/SKCMcoldata.RData")

SKCMexp <- SKCMexp[,SKCMcoldata$Cell]
rownames(SKCMcoldata)<- SKCMcoldata$Cell


design <- model.matrix(~ Tprolif, data = SKCMcoldata)
fit <- lmFit(SKCMexp, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, adjust="BH", number=Inf)
print(results)
SKCMTprolifresult <- results
save(SKCMTprolifresult,file="/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/SKCMTprolifresult.RData")


design <- model.matrix(~ Malignant, data = SKCMcoldata)
fit <- lmFit(SKCMexp, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, adjust="BH", number=Inf)
print(results)
SKCMMalignantresult <- results
save(SKCMMalignantresult,file="/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/SKCMMalignantresult.RData")


design <- model.matrix(~ Fibroblasts, data = SKCMcoldata)
fit <- lmFit(SKCMexp, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, adjust="BH", number=Inf)
print(results)
SKCMFibroblastsresult <- results
save(SKCMFibroblastsresult,file="/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/SKCMFibroblastsresult.RData")


design <- model.matrix(~ group4, data = SKCMcoldata)
fit <- lmFit(SKCMexp, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, adjust="BH", number=Inf)
print(results)
SKCMgroup4result <- results
save(SKCMgroup4result,file="/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/SKCMgroup4result.RData")


design <- model.matrix(~ group6, data = SKCMcoldata)
fit <- lmFit(SKCMexp, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, adjust="BH", number=Inf)
print(results)
SKCMgroup6result <- results
save(SKCMgroup6result,file="/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/SKCMgroup6result.RData")


#STAD
table(STADsingle$celltype)
STADsingle$DC <- ifelse(STADsingle$celltype=="DC","yes","no")
STADsingle$Mast <- ifelse(STADsingle$celltype=="Mast","yes","no")
STADsingle$Mono_Macro <- ifelse(STADsingle$celltype=="Mono/Macro","yes","no")
STADsingle$Endothelial <- ifelse(STADsingle$celltype=="Endothelial","yes","no")
STADsingle$Fibroblasts <- ifelse(STADsingle$celltype=="Fibroblasts","yes","no")
STADsingle$group1 <- ifelse(STADsingle$model1>median(STADsingle$model1),"high","low")
STADsingle$group3 <- ifelse(STADsingle$model3>median(STADsingle$model3),"high","low")
STADsingle$group5 <- ifelse(STADsingle$model5>median(STADsingle$model5),"high","low")


head(STADsingle)
STADcoldata <- STADsingle[,c("Cell","DC","Mast","Mono_Macro","Endothelial","Fibroblasts","group1","group3","group5")]
save(STADcoldata,file = "D:/细胞状态2/论文/fig3/fig3d/STADcoldata.RData")
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/STADcoldata.RData")

STADexp <- STADexp[,STADcoldata$Cell]
rownames(STADcoldata)<- STADcoldata$Cell


design <- model.matrix(~ DC, data = STADcoldata)
fit <- lmFit(STADexp, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, adjust="BH", number=Inf)
print(results)
STADDCresult <- results
save(STADDCresult,file="/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/STADDCresult.RData")


design <- model.matrix(~ Mast, data = STADcoldata)
fit <- lmFit(STADexp, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, adjust="BH", number=Inf)
print(results)
STADMastresult <- results
save(STADMastresult,file="/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/STADMastresult.RData")


design <- model.matrix(~ Mono_Macro, data = STADcoldata)
fit <- lmFit(STADexp, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, adjust="BH", number=Inf)
print(results)
STADMono_Macroresult <- results
save(STADMono_Macroresult,file="/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/STADMono_Macroresult.RData")


design <- model.matrix(~ Endothelial, data = STADcoldata)
fit <- lmFit(STADexp, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, adjust="BH", number=Inf)
print(results)
STADEndothelialresult <- results
save(STADEndothelialresult,file="/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/STADEndothelialresult.RData")


design <- model.matrix(~ Fibroblasts, data = STADcoldata)
fit <- lmFit(STADexp, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, adjust="BH", number=Inf)
print(results)
STADFibroblastsresult <- results
save(STADFibroblastsresult,file="/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/STADFibroblastsresult.RData")


design <- model.matrix(~ group1, data = STADcoldata)
fit <- lmFit(STADexp, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, adjust="BH", number=Inf)
print(results)
STADgroup1result <- results
save(STADgroup1result,file="/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/STADgroup1result.RData")


design <- model.matrix(~ group3, data = STADcoldata)
fit <- lmFit(STADexp, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, adjust="BH", number=Inf)
print(results)
STADgroup3result <- results
save(STADgroup3result,file="/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/STADgroup3result.RData")


design <- model.matrix(~ group5, data = STADcoldata)
fit <- lmFit(STADexp, design)
fit <- eBayes(fit)
results <- topTable(fit, coef=2, adjust="BH", number=Inf)
print(results)
STADgroup5result <- results
save(STADgroup5result,file="/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/STADgroup5result.RData")

#dim(STADexp)
#dim(results)
#sig_genes <- rownames(results[results$adj.P.Val < 0.05, ])
#length(sig_genes)

#处理差异基因
library(ggrepel)
load("D:/细胞状态2/论文/fig3/fig3e/paadFibresult.RData")
load("D:/细胞状态2/论文/fig3/fig3e/paadMyoresult.RData")
load("D:/细胞状态2/论文/fig3/fig3e/paadgroup1result.RData")
paadFibsig_genes <- rownames(paadFibresult[paadFibresult$adj.P.Val < 0.05& abs(paadFibresult$logFC)>0.5, ])
paadMyosig_genes <- rownames(paadMyoresult[paadMyoresult$adj.P.Val < 0.05& abs(paadMyoresult$logFC)>0.5, ])
paadgroup1sig_genes <- rownames(paadgroup1result[paadgroup1result$adj.P.Val < 0.05& abs(paadgroup1result$logFC)>0.5, ])

paadFibsig_genes2 <- intersect(paadFibsig_genes,paadgroup1sig_genes)

max(paadFibresult$logFC)
min(paadFibresult$logFC)

paadFibsigresult <- paadFibresult[paadFibsig_genes,]
paadMyosigresult <- paadMyoresult[paadMyosig_genes,]

paadinter <- intersect(paadFibsig_genes,paadMyosig_genes)
paadinter2 <- intersect(paadinter,paadgroup1sig_genes)


load("D:/细胞状态2/论文/fig3/fig3e/SKCMTprolifresult.RData")
load("D:/细胞状态2/论文/fig3/fig3e/SKCMFibroblastsresult.RData")
load("D:/细胞状态2/论文/fig3/fig3e/SKCMMalignantresult.RData")
#load("D:/细胞状态2/论文/fig3/fig3e/SKCMgroup4result.RData")
#load("D:/细胞状态2/论文/fig3/fig3e/SKCMgroup6result.RData")
skcmTprsig_genes <- rownames(SKCMTprolifresult[SKCMTprolifresult$adj.P.Val < 0.05& abs(SKCMTprolifresult$logFC)>0.5, ])
skcmFibsig_genes <- rownames(SKCMFibroblastsresult[SKCMFibroblastsresult$adj.P.Val < 0.05& abs(SKCMFibroblastsresult$logFC)>0.5, ])
skcmMalsig_genes <- rownames(SKCMMalignantresult[SKCMMalignantresult$adj.P.Val < 0.05& abs(SKCMMalignantresult$logFC)>0.5, ])



load("D:/细胞状态2/论文/fig3/fig3e/STADDCresult.RData")
load("D:/细胞状态2/论文/fig3/fig3e/STADEndothelialresult.RData")
load("D:/细胞状态2/论文/fig3/fig3e/STADFibroblastsresult.RData")
load("D:/细胞状态2/论文/fig3/fig3e/STADMastresult.RData")
load("D:/细胞状态2/论文/fig3/fig3e/STADMono_Macroresult.RData")

staddcsig_genes <- rownames(STADDCresult[STADDCresult$adj.P.Val < 0.05& abs(STADDCresult$logFC)>0.5, ])
stadendsig_genes <- rownames(STADEndothelialresult[STADEndothelialresult$adj.P.Val < 0.05& abs(STADEndothelialresult$logFC)>0.5, ])
stadfibsig_genes <- rownames(STADFibroblastsresult[STADFibroblastsresult$adj.P.Val < 0.05& abs(STADFibroblastsresult$logFC)>0.5, ])
stadmastsig_genes <- rownames(STADMastresult[STADMastresult$adj.P.Val < 0.05& abs(STADMastresult$logFC)>0.5, ])
stadmonosig_genes <- rownames(STADMono_Macroresult[STADMono_Macroresult$adj.P.Val < 0.05& abs(STADMono_Macroresult$logFC)>0.5, ])



#功能富集
library(clusterProfiler)
library(org.Hs.eg.db)


#valid_keys <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
gene_list <- bitr(paadinter, fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Hs.eg.db)

go_results <- enrichGO(gene = gene_list$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "ALL",  # 可以选择 "BP", "CC", "MF"
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

library(enrichplot)
library(ggplot2)

go_results@result$pvalue2 <- -log10(go_results@result$pvalue)
ggplot(go_results, aes(x = Description, y = Count, fill = pvalue2)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "steelblue", high = "#c11e21") +
  theme_minimal() +
  labs(x = "Description", y = "Count") +
  ggtitle("GO Barplot with Gradient Fill") +
  coord_flip()+
  theme_minimal()+ 
  theme(panel.grid = element_blank())


head(go_results)
dim(go_results@result$Description)
max(go_results@result$pvalue)
min(go_results@result$pvalue)

go_resultsdataframe <- as.data.frame(go_results@result)
go_resultsdataframe <- go_resultsdataframe[go_resultsdataframe$pvalue<0.05,]
go_resultsdataframe <- go_resultsdataframe[go_resultsdataframe$p.adjust<0.05,]
go_resultsdataframe$geneID <- NULL
head(go_resultsdataframe,n=10)

paadFibmodel1go_resultsdataframe <- go_resultsdataframe



#cellstate标签基因

cellstategene <- read.csv("D:/细胞状态/cell_state.csv")
head(cellstategene)
table(cellstategene$Signature)
Angiogenesisgene <- cellstategene[cellstategene$Signature=="Angiogenesis",]
EMTgene <- cellstategene[cellstategene$Signature=="EMT",]
Invasiongene <- cellstategene[cellstategene$Signature=="Invasion",]
Metastasisgene <- cellstategene[cellstategene$Signature=="Metastasis",]

Apoptosisgene <- cellstategene[cellstategene$Signature=="Apoptosis",]

Hypoxiagene <- cellstategene[cellstategene$Signature=="Hypoxia",]

Inflammationgene <- cellstategene[cellstategene$Signature=="Inflammation",]
Quiescencegene <- cellstategene[cellstategene$Signature=="Quiescence",]

CellCyclegene <- cellstategene[cellstategene$Signature=="Cell Cycle",]
DNAdamagegene <- cellstategene[cellstategene$Signature=="DNA damage",]
DNArepairgene <- cellstategene[cellstategene$Signature=="DNA repair",]

Proliferationgene <- cellstategene[cellstategene$Signature=="Proliferation",]

Differentiationgene <- cellstategene[cellstategene$Signature=="Differentiation",]
Stemnessgene <- cellstategene[cellstategene$Signature=="Stemness",]

#paadFibsig_model1
paadFibsig_Angiogenesis_genes <- intersect(paadFibsig_genes,Angiogenesisgene$GeneName)
paadFibsig_EMT_genes <- intersect(paadFibsig_genes,EMTgene$GeneName)
paadFibsig_Invasion_genes <- intersect(paadFibsig_genes,Invasiongene$GeneName)
paadFibsig_Metastasis_genes <- intersect(paadFibsig_genes,Metastasisgene$GeneName)

head(paadFibresult)
paadFibresult[paadFibsig_Angiogenesis_genes,]
paadFibresult[paadFibsig_EMT_genes,]
paadFibresult[paadFibsig_Invasion_genes,]
paadFibresult[paadFibsig_Metastasis_genes,]

#paadMyosig_model1
paadMyosig_Angiogenesis_genes <- intersect(paadMyosig_genes,Angiogenesisgene$GeneName)
paadMyosig_EMT_genes <- intersect(paadMyosig_genes,EMTgene$GeneName)
paadMyosig_Invasion_genes <- intersect(paadMyosig_genes,Invasiongene$GeneName)
paadMyosig_Metastasis_genes <- intersect(paadMyosig_genes,Metastasisgene$GeneName)

paadMyoresult[paadMyosig_Angiogenesis_genes,]
paadMyoresult[paadMyosig_EMT_genes,]
paadMyoresult[paadMyosig_Invasion_genes,]
paadMyoresult[paadMyosig_Metastasis_genes,]

intersect(paadMyosig_Angiogenesis_genes,paadFibsig_Angiogenesis_genes)
intersect(paadMyosig_EMT_genes,paadFibsig_EMT_genes)
intersect(paadMyosig_Invasion_genes,paadFibsig_Invasion_genes)
intersect(paadMyosig_Metastasis_genes,paadFibsig_Metastasis_genes)

#paad交集
paadsig_Angiogenesis_genes <- intersect(paadinter,Angiogenesisgene$GeneName)
paadsig_EMT_genes <- intersect(paadinter,EMTgene$GeneName)
paadsig_Invasion_genes <- intersect(paadinter,Invasiongene$GeneName)
paadsig_Metastasis_genes <- intersect(paadinter,Metastasisgene$GeneName)

paadFibresult[paadsig_Angiogenesis_genes,]
paadFibresult[paadsig_EMT_genes,]
paadFibresult[paadsig_Invasion_genes,]
paadFibresult[paadsig_Metastasis_genes,]

paadMyoresult[paadsig_Angiogenesis_genes,]
paadMyoresult[paadsig_EMT_genes,]
paadMyoresult[paadsig_Invasion_genes,]
paadMyoresult[paadsig_Metastasis_genes,]

intersect(paadMyosig_Angiogenesis_genes,paadFibsig_Angiogenesis_genes)
intersect(paadMyosig_EMT_genes,paadFibsig_EMT_genes)
intersect(paadMyosig_Invasion_genes,paadFibsig_Invasion_genes)
intersect(paadMyosig_Metastasis_genes,paadFibsig_Metastasis_genes)


paadset <- list(
  paadMyosig_genes=paadMyosig_genes,
  paadFibsig_genes=paadFibsig_genes,
  Invasion=Invasiongene$GeneName
)
ggvenn(paadset,
       fill_color=c('#D2868C','#B1AD8E','#72C4D1'),
       fill_alpha = 0.4,
       stroke_color=NA
)
library(VennDiagram)
vstad <-venn.diagram(x = paadset, 
                     filename = NULL, 
                     scaled = T, # 根据比例显示大小
                     alpha=c(0.4, 0.4,0.4), # 设置每个区块的透明度
                     ## 下面是除了标题外，图形其他元素的设置参数
                     # 图形元素设置：圈
                     lwd = 1, # 圆圈线条的粗细：1 2 3 4 5 6
                     lty = 1, # 圆圈线条的类型：1为实线，2为虚线，blank为无线条
                     col = NA, # 圆圈线条颜色
                     fill = c('#D2868C','#B1AD8E','#72C4D1'),  # 圆圈颜色
                     # 图形元素设置：数字
                     cex = 1, # 数字大小
                     fontface = "bold", # 加粗
                     fonrfamily = "sans", # 数字字体
                     # 图形元素设置：标签即(category)
                     cat.cex = 1,  # 标签字体大小
                     cat.col = "black",  # 标签字体色
                     cat.fontface = "bold",  # 加粗
                     cat.default.pos = "outer",  # 标签内外位置, 在圆圈内还是圆圈外，outer 内 text 外
                     cat.pos = c(0, 0,0),  # 标签旋转位置，用圆的度数
                     cat.dist = c(0.05,0.03,0.03),  # 标签离圆圈位置,离圆的距离，如果标签与圆圈重叠，可以调整这个参数
                     cat.fontfamily = "sans",  # 标签字体
)

cowplot::plot_grid(vstad)

#skcm

#qipao_paadFib$logFC <- qipao_paadFib$log2FC
#qipao_paadFib$log2FC <- NULL
qipao_paadFib<- paadFibresult[c("COL1A1","BGN","TAGLN"),]
qipao_paadMyo <- paadMyoresult[c("COL1A1","BGN","TAGLN"),]
qipao_paadFib$group <- "PAAD_Fib"
qipao_paadMyo$group <- "PAAD_Myo"
qipao_paadFib$gene <- rownames(qipao_paadFib)
qipao_paadMyo$gene <- rownames(qipao_paadMyo)
rownames(qipao_paadFib) <- NULL
rownames(qipao_paadMyo) <- NULL

qipaotest <- rbind(qipao_paadFib,qipao_paadMyo)


ggplot(qipaotest, aes(x = gene, y = group)) +
  geom_point(aes(size = logFC, color = logFC,fill=logFC), shape = 21) +
  scale_size_continuous(range = c(1, 10)) +  # 调整气泡大小范围
  scale_color_gradient2(low = "steelblue", mid = "white", high = "#c11e21", midpoint = 0, limits = c(-3.5, 3.5)) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "#c11e21", midpoint = 0, limits = c(-3.5, 3.5)) +
  theme_minimal()+
  theme_bw() + 
  theme(panel.grid=element_blank())


