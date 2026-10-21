#单机
load("D:/细胞状态2/论文/fig1/cellstatemergedata.RData")
#fig1b 测试皮尔森
head(cellstatemergedata)

testperson <- cellstatemergedata
head(testperson)
testperson <- testperson[,2:15]
testperson <- scale(testperson)
cor_pearson <- cor(testperson, method = 'pearson')
cor_pearson
library(corrplot)
pdf("D:/细胞状态2/论文/fig1/gsvacor.pdf")
my_colors <- colorRampPalette(c("steelblue", "white", "#c11e21"))(100)
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
hc <- hclust(as.dist(1 - cor_pearson))
rect.hclust(hc, k = 7) 
dev.off()

hc <- hclust(dist(cor_pearson))


#计算foldchange值
#筛选正常与癌症
load("/boot2/xff/cellstate/cellstatemergedata.RData")
load("D:/细胞状态2/论文/fig1/cellstatemergedata.RData")
head(cellstatemergedata)
table(cellstatemergedata$Source)
table(cellstatemergedata$`Celltype (malignancy)`)
cellstatemergedata <- cellstatemergedata[!cellstatemergedata$`Celltype (malignancy)`=="Others",]

table(cellstatemergedata$`Celltype (major-lineage)`)
names(table(cellstatemergedata$`Celltype (major-lineage)`))
cellstatemergedata$celltype <- cellstatemergedata$`Celltype (major-lineage)`

cellstatemergedata$celltype[cellstatemergedata$celltype %in% c("TProlif","Tproilf") ] <- "Tprolif"
cellstatemergedata$celltype[cellstatemergedata$celltype %in% c("Epithelial ","Ciliated","Hepatocytes","Keratinocytes") ] <- "Epithelial"
cellstatemergedata$celltype[cellstatemergedata$celltype %in% c("Fibroblast") ] <- "Fibroblasts"
cellstatemergedata$celltype[cellstatemergedata$celltype %in% c("Mono.Macro") ] <- "Mono/Macro"

cellstatemergedata$Source[cellstatemergedata$celltype %in% c("Malignant")] <- "Tumor"


TUMOR <- cellstatemergedata[cellstatemergedata$Source %in% c("AML", "Ascites", "B-ALL","BMMC",
                                                             "center","cortex","Human tumor",
                                                             "Leptomeningeal metastasis","mBrain",
                                                             "Metastasis","Metastatic","mLiver",
                                                             "mFallopian_tube","nodule","Primary",
                                                             "Relapsed","T-ALL","TIL","tLN","tumor",
                                                             "Tumor"), ]
NORMAL <- cellstatemergedata[cellstatemergedata$Source %in% c("Normal", "NAT", "Adjacent normal",
                                                              "PBMC","healthy","BM","Cell culture",
                                                              "Skin","LN","Tonsil"), ]
head(TUMOR)
head(NORMAL)


#去除部分癌症与细胞
singlecan3 <- c("BCC","BRCA","CESC","CRC","ESCA",
                "KIRC","LIHC","NSCLC","OV","PAAD",
                "PRAD","SKCM","STAD","THCA","UCEC")

TUMOR2 <- TUMOR[TUMOR$cancertype %in% c("BCC","BRCA","CESC","CRC","ESCA",
                                                "KIRC","LIHC","NSCLC","OV","PAAD",
                                                "PRAD","SKCM","STAD","THCA","UCEC"),]

NORMAL2 <- NORMAL[NORMAL$cancertype %in% c("BCC","BRCA","CESC","CRC","ESCA",
                                                  "KIRC","LIHC","NSCLC","OV","PAAD",
                                                  "PRAD","SKCM","STAD","THCA","UCEC"),]
#NORMAL2[NORMAL2$celltype=="Malignant",]
table(TUMOR2$celltype)
table(NORMAL2$celltype)
write.csv(table(TUMOR2$celltype),file = "D:/细胞状态2/论文/fig1/TUMORcelltype.csv")
write.csv(table(NORMAL$celltype),file = "D:/细胞状态2/论文/fig1/NORMALcelltype.csv")
intersect(names(table(TUMOR2$celltype)),names(table(NORMAL2$celltype)))
write.csv(intersect(names(table(TUMOR2$celltype)),names(table(NORMAL2$celltype))),
          file = "D:/细胞状态2/论文/fig1/intersectcelltype.csv")



TUMOR3 <- TUMOR2[TUMOR2$`Celltype (major-lineage)` %in% c("B","CD4Tconv","CD8T","CD8Tex","DC",
                                                        "Endothelial","Epithelial","Fibroblasts","Malignant","Mast",
                                                        "Mono/Macro","Myofibroblasts","NK","Tprolif","Treg"),]
NORMAL3 <- NORMAL2[NORMAL2$`Celltype (major-lineage)` %in% c("B","CD4T","CD4Tconv","CD8T","CD8Tex","DC",
                                                           "Endothelial","Epithelial","Fibroblasts","Malignant","Mast",
                                                           "Mono/Macro","Myofibroblasts","NK","Tprolif","Treg"),]



head(TUMOR3)
head(NORMAL3)

for (i in singlecan3) {
  print(i)
  TUMORtest <- TUMOR3[TUMOR3$cancertype==i,]
  NORMALtest <- NORMAL3[NORMAL3$cancertype==i,]
  print(table(TUMORtest$celltype))
  print(table(NORMALtest$celltype))
  print(intersect(names(table(TUMORtest$celltype)),names(table(NORMALtest$celltype))))
}



singlecell3 <- c("B","CD4Tconv","CD8T","CD8Tex","DC",
                 "Endothelial","Epithelial","Fibroblasts","Malignant","Mast",
                 "Mono/Macro","Myofibroblasts","NK","Tprolif","Treg")



head(TUMOR3)
head(NORMAL3)
table(TUMOR3$celltype)
table(NORMAL3$celltype)
table(TUMOR3$celltype,TUMOR3$cancertype)
table(NORMAL3$celltype,NORMAL3$cancertype)


load("D:/细胞状态2/论文/fig1/TUMOR3.RData")
head(TUMOR3)
testperson <- TUMOR3
head(testperson)
testperson <- testperson[,2:15]
testperson <- scale(testperson)
cor_pearson <- cor(testperson, method = 'pearson')
cor_pearson
library(corrplot)
pdf("D:/细胞状态2/论文/fig1/fig1bf/normalgsvacor.pdf")
my_colors <- colorRampPalette(c("steelblue", "white", "#c11e21"))(100)
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
hc <- hclust(as.dist(1 - cor_pearson))
rect.hclust(hc, k = 4) 
dev.off()

cor_pearsondata <- as.data.frame(cor_pearson)
write.csv(cor_pearsondata,file="D:/细胞状态2/论文/fig1/cor_pearsondata.csv")

testperson <- NORMAL3
head(testperson)
testperson <- testperson[,2:15]
testperson <- scale(testperson)
cor_pearson <- cor(testperson, method = 'pearson')
cor_pearson
library(corrplot)
pdf("D:/细胞状态2/论文/fig1/fig1bf/normalgsvacor.pdf")
my_colors <- colorRampPalette(c("steelblue", "white", "#c11e21"))(100)
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
hc <- hclust(as.dist(1 - cor_pearson))
rect.hclust(hc, k = 7) 
dev.off()


save(TUMOR3,NORMAL3,file = "D:/细胞状态2/论文/fig1/TUMOR3.RData")
write.csv(TUMOR3,file = "D:/细胞状态2/论文/fig1/TUMOR3.csv")
load("D:/细胞状态2/论文/fig1/TUMOR3.RData")
head(TUMOR3)
table(TUMOR3$cancertype)
table(TUMOR3$celltype)

#每种细胞类型计算皮尔森相关系数
for (m in singlecell3) {
  TUMORtest <- TUMOR3[TUMOR3$celltype==m,]
  #i <- 1
  everycellperson <- TUMORtest[,2:15]
  everycellperson <- scale(everycellperson)
  cor_pearson <- cor(everycellperson, method = 'pearson')
  cor_pearson
  myfilename <- paste("D:/细胞状态2/论文/fig1/fig1bf/cor/everycell/",m,"cor.pdf",sep = "")
  pdf(myfilename)
  corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = FALSE, diag = FALSE,col = my_colors)
  hc <- hclust(as.dist(1 - cor_pearson))
  rect.hclust(hc, k = 7) 
  dev.off()
}
dev.new()


#每种癌症类型计算皮尔森相关系数
for (m in singlecan3) {
  TUMORtest <- TUMOR3[TUMOR3$cancertype==m,]
  #i <- 1
  everycancerperson <- TUMORtest[,2:15]
  everycancerperson <- scale(everycancerperson)
  cor_pearson <- cor(everycancerperson, method = 'pearson')
  cor_pearson
  myfilename <- paste("D:/细胞状态2/论文/fig1/fig1bf/cor/everycancer/",m,"cor.pdf",sep = "")
  pdf(myfilename)
  corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = FALSE, diag = FALSE,col = my_colors)
  hc <- hclust(as.dist(1 - cor_pearson))
  rect.hclust(hc, k = 7) 
  dev.off()
}


#查看不同model的分布情况

fig1gcell<- read.csv("D:/细胞状态2/论文/fig1/fig1bf/fig1g1.csv",header = T)
fig1gcancer <- read.csv("D:/细胞状态2/论文/fig1/fig1bf/fig1g2.csv",header = T)


library(pheatmap)
library(RColorBrewer)
custom_colors2 <- colorRampPalette(c("steelblue","#9FBED8","#e39b9c","#c11e21"))
num_breaks <- 100
breaks2 <- seq(0.25, 1, length.out = num_breaks + 1)
library(tidyr)

long_cell <- gather(fig1gcell, key = "Model", value = "Value", -X)
long_cancer <- gather(fig1gcancer, key = "Model", value = "Value", -X)


pdf("D:/细胞状态2/论文/fig1/fig1bf/cor/fig1gcell.pdf")
ggplot(long_cell, aes(x = Model, y = X, fill = Value,size =Value)) +
  geom_point(shape = 21, color = "white") +  # 使用形状21绘制气泡
  scale_size_continuous(range = c(2,6)) +  # 设置气泡大小范围
  scale_fill_gradient2(low = "steelblue",mid="white", high = "#c11e21",midpoint = 0.625, name = "Value") +  # 设置颜色梯度
  theme_minimal() +  # 使用简洁主题
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),  # 去除主要网格线
        panel.grid.minor = element_blank())  # 旋转X轴标签以提高可读性
dev.off()

dev.new()

pdf("D:/细胞状态2/论文/fig1/fig1bf/cor/fig1gcancer.pdf")
ggplot(long_cancer, aes(x = Model, y = X, fill = Value,size =Value)) +
  geom_point(shape = 21, color = "white") +  # 使用形状21绘制气泡
  scale_size_continuous(range = c(2,6)) +  # 设置气泡大小范围
  scale_fill_gradient2(low = "steelblue",mid="white", high = "#c11e21",midpoint = 0.625, name = "Value") +  # 设置颜色梯度
  theme_minimal() +  # 使用简洁主题
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),  # 去除主要网格线
        panel.grid.minor = element_blank())  # 旋转X轴标签以提高可读性
dev.off()

