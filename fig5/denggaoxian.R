library(reshape2) #数据处理相关
library(ggplot2) # 绘图相关
library(grDevices) #绘图颜色相关
library(RColorBrewer)#绘图颜色相关
library(Seurat)
library(ggdensity)
library(ggblanket)
library(ggsci)
library(sp)
library(corrplot)
library(tidyr)

#测试数据
# 创建一个示例网格数据
x <- seq(-10, 10, by = 0.1)
y <- seq(-10, 10, by = 0.1)
grid <- expand.grid(x = x, y = y)

# 定义 Z 值（例如，Z = sin(sqrt(x^2 + y^2))）
grid$z <- with(grid, sin(sqrt(x^2 + y^2)))

# 查看数据
head(grid)
class(grid)
# 绘制等高线图
ggplot(grid, aes(x = x, y = y, z = z)) +
  geom_contour() +
  labs(title = "Contour Plot") +
  theme_minimal()

str(grid$z)


#测试2
# 创建示例数据
set.seed(123)  # 为了重现性
n <- 100
data <- data.frame(
  x = rnorm(n, mean = 0, sd = 10),  # x坐标
  y = rnorm(n, mean = 0, sd = 10),  # y坐标
  gene_expression = rnorm(n, mean = 10, sd = 2)  # 基因表达量
)

# 你可以替换成你的实际数据，注意要有空间坐标和对应的基因表达量
ggplot(data, aes(x = x, y = y, z = gene_expression)) +
  geom_contour_filled() +  # 等高线图，使用填充色
  labs(title = "Spatial Transcriptomics Contour Plot", 
       x = "X Coordinate", 
       y = "Y Coordinate", 
       fill = "Gene Expression") +
  theme_minimal()


#paad
stmetadata <- read.csv("D:/细胞状态2/论文/空间/st2/stmeta.data.csv")
head(stmetadata)
position1 <- read.csv("D:/细胞状态2/论文/空间/st/SS1905133_processed/spatial/tissue_positions_list.csv",header = F)
head(position1)

stmetadata2 <- merge(stmetadata,position1,by.x="cell",by.y="V1")
head(stmetadata2)


min(stmetadata2$V6) 
range(stmetadata2$model1)


# 绘制等高线图
stmetadata3 <- stmetadata2
head(stmetadata3)

#write.csv(stmetadata3,file = "D:/细胞状态2/论文/空间/st2/stmeta.data3.csv")

install.packages("ggdensity")
install.packages("ggblanket")
library(ggdensity)
library(ggblanket)
library(ggsci)

stmetadata3$group1 <- as.factor(stmetadata3$group1)
stmetadata3$group1 <- factor(stmetadata3$group1,levels = c("low","high"))
head(stmetadata3)
denggaoxian <- ggplot(stmetadata3,aes(x = V5, y = V6, fill = group1)) +
  geom_hdr() + 
  geom_point(shape = 21) +
  scale_fill_aaas() +
  facet_wrap(vars(group1)) +
  theme_bw() + 
  theme(panel.grid=element_blank())

ggplot(stmetadata3,aes(x = V5, y = V6, fill = group1)) +
  geom_hdr() + 
  scale_fill_aaas() +
  facet_wrap(vars(group1)) 

ggplot(stmetadata3,aes(x = V5, y = V6, fill = group1))+
  geom_point(shape = 21) +
  scale_fill_manual(values = c("high"="#c11e21","low"="steelblue"))

head(stmetadata3)
dim(stmetadata3)


# 提取 geom_hdr 生成的数据
plot_data <- ggplot_build(denggaoxian)$data[[1]]  
dim(plot_data)
head(plot_data)
table(plot_data$fill)
table(plot_data$group)
table(plot_data$probs)
table(plot_data$linetype)
table(plot_data$probs,plot_data$group)
max(plot_data$x)

dim(ggplot_build(denggaoxian)$data[[2]] )
head(ggplot_build(denggaoxian)$data[[2]])



# 查看结果
head(stmetadata3)

prob50_spots <- plot_data[plot_data$probs == "50%", ]
table(prob50_spots$group)
table(prob50_spots$probs)
head(prob50_spots)
prob50_spotshigh <- plot_data[plot_data$group == "2-004", ]
prob80_spotshigh <- plot_data[plot_data$group == "2-003", ]
prob95_spotshigh <- plot_data[plot_data$group == "2-002", ]
prob99_spotshigh <- plot_data[plot_data$group == "2-001", ]

max(prob50_spots$x)
min(prob50_spots$x)
max(prob50_spots$y)
min(prob50_spots$y)



# 安装并加载 sp 包
install.packages("sp")
library(sp)

# 封闭曲线（多边形）坐标
polygon_x <- prob50_spotshigh$x  
polygon_y <- prob50_spotshigh$y  

polygon_x <- prob80_spotshigh$x
polygon_y <- prob80_spotshigh$y

dim(prob80_spotshigh[prob80_spotshigh$y<17000,])
testmetadata<- prob80_spotshigh[prob80_spotshigh$y<17000,]
polygon_x <- testmetadata$x
polygon_y <- testmetadata$y



# 待检测点的坐标
points_x <- stmetadata3$V5
points_y <- stmetadata3$V6


# 将封闭曲线和点坐标转为矩阵
polygon_coords <- cbind(polygon_x, polygon_y)
points_coords <- cbind(points_x, points_y)

# 创建 SpatialPolygons对象（封闭曲线多边形）
polygon_sp <- Polygon(polygon_coords)
polygons_sp <- Polygons(list(polygon_sp), "polygon1")
sp_poly <- SpatialPolygons(list(polygons_sp))

# 使用 polygon() 函数绘制多边形
plot(1, 1, type="n", xlim=c(4500, 55000), ylim=c(4500, 50000), xlab="X", ylab="Y")
polygon(polygon_x, polygon_y, col="skyblue", border="black")


# 使用 point.in.polygon 函数检查每个点是否在多边形内
result <- point.in.polygon(points_x, points_y, polygon_x, polygon_y)

# 输出结果：结果为1表示在多边形内，0表示不在多边形内
print(result)
length(result)
stmetadata3$result50 <- result
stmetadata3$result80 <- result
#stmetadata3$result95 <- result
#stmetadata3$result99 <- result
table(stmetadata3$result50)


stmetadata3$result <- ifelse(stmetadata3$result50==1,"50%",
                             ifelse(stmetadata3$result80==1,"80%","95%"))
table(stmetadata3$result)

head(stmetadata3)

table(stmetadata3$result,stmetadata3$group1)
table(stmetadata3$result,stmetadata3$celltype)
head(stmeta.data)
head(stmetadata3)
write.csv(stmetadata3,file = "D:/细胞状态2/论文/空间/paad/denggaoxian.csv")
paaddenggaoxian <- read.csv("D:/细胞状态2/论文/空间/manypaad/denggaoxian.csv")
head(paaddenggaoxian)
table(paaddenggaoxian$celltype)

stmetadata3 <- pivot_longer(paaddenggaoxian, 
                            cols = c("EMT","Angiogenesis","Invasion","Metastasis","model1"), 
                            names_to = "cellstate", 
                            values_to = "Value")
head(stmetadata3)
stmetadata3$result <- as.factor(stmetadata3$result)

stmetadataB <- stmetadata3[stmetadata3$celltype=="B",]
stmetadataend <- stmetadata3[stmetadata3$celltype=="Endothelial",]
stmetadataepi <- stmetadata3[stmetadata3$celltype=="Epithelial",]
stmetadatafib <- stmetadata3[stmetadata3$celltype=="Fibroblasts",]
stmetadatamyeloid <- stmetadata3[stmetadata3$celltype=="Myeloid",]
stmetadataT <- stmetadata3[stmetadata3$celltype=="T/NK",]


library(ggpubr)
ggplot(data = stmetadataend, aes(x = cellstate, y = Value, fill = result)) +
  geom_violin(position = position_dodge(0.8))+
  geom_boxplot(width=0.2,
               position = position_dodge(0.8))+  #每组的小提琴图之间的距离
  stat_compare_means(aes(group = result))+ 
  scale_fill_manual(values = c("#e99593", "#99bca0","steelblue"))+
  theme(axis.text.x = element_text(size = 20,angle = 30, hjust = 0.5, vjust = 0.5))+
  theme_minimal()+ 
  theme(panel.grid = element_blank())

pd <- position_dodge2(width = 0.8, preserve = "single")

ggplot(stmetadataB, aes(x = cellstate, y = Value, fill = result)) +
  geom_violin(position = pd, trim = FALSE) +
  geom_boxplot(width = 0.2, position = pd, outlier.shape = NA) +
  stat_compare_means(aes(group = result)) +
  scale_fill_manual(values = c("#e99593", "#99bca0", "steelblue")) +
  theme_minimal() +
  theme(panel.grid = element_blank())


a








#STAD
position$V1 <- rownames(position)
stmetadata2 <- merge(stmeta.data,position,by.x="cell",by.y="V1")

head(stmetadata2)
stmetadata2$group1 <- as.factor(stmetadata2$group1)
stmetadata2$group1 <- factor(stmetadata2$group1,levels = c("low","high"))
head(stmetadata2)
denggaoxian <- ggplot(stmetadata2,aes(x = Spatial_1, y = Spatial_2, fill = group1)) +
  geom_hdr() + 
  geom_point(shape = 21) +
  scale_fill_aaas() +
  facet_wrap(vars(group1)) +
  theme_bw() + 
  theme(panel.grid=element_blank())

ggplot(stmetadata2,aes(x = Spatial_1, y = Spatial_2, fill = group1)) +
  geom_hdr() + 
  geom_point(shape = 21) +
  scale_fill_aaas() +
  facet_wrap(vars(group1)) 

ggplot(stmetadata2,aes(x = Spatial_1, y = Spatial_2, fill = group1))+
  geom_point(shape = 21) +
  scale_fill_manual(values = c("high"="#c11e21","low"="steelblue"))






# 提取 geom_hdr 生成的数据
plot_data <- ggplot_build(denggaoxian)$data[[1]]  
dim(plot_data)
head(plot_data)
table(plot_data$fill)
table(plot_data$group)
table(plot_data$probs)
table(plot_data$linetype)
table(plot_data$probs,plot_data$group)
max(plot_data$x)

dim(ggplot_build(denggaoxian)$data[[2]] )
head(ggplot_build(denggaoxian)$data[[2]])



# 查看结果
head(stmetadata3)
head(stmetadata2)

prob50_spots <- plot_data[plot_data$probs == "50%", ]
table(prob50_spots$group)
table(prob50_spots$probs)
head(prob50_spots)
prob50_spotshigh <- plot_data[plot_data$group == "2-004", ]
prob80_spotshigh <- plot_data[plot_data$group == "2-003", ]
prob95_spotshigh <- plot_data[plot_data$group == "2-002", ]
prob99_spotshigh <- plot_data[plot_data$group == "2-001", ]

max(prob50_spots$x)
min(prob50_spots$x)
max(prob50_spots$y)
min(prob50_spots$y)



# 安装并加载 sp 包
install.packages("sp")
library(sp)

# 封闭曲线（多边形）坐标
polygon_x <- prob50_spotshigh$x  
polygon_y <- prob50_spotshigh$y  

polygon_x <- prob80_spotshigh$x
polygon_y <- prob80_spotshigh$y

dim(prob80_spotshigh[prob80_spotshigh$y<17000,])
testmetadata<- prob80_spotshigh[prob80_spotshigh$y<17000,]
polygon_x <- testmetadata$x
polygon_y <- testmetadata$y


polygon_x <- prob95_spotshigh$x
polygon_y <- prob95_spotshigh$y

polygon_x <- prob99_spotshigh$x
polygon_y <- prob99_spotshigh$y

# 待检测点的坐标
points_x <- stmetadata2$Spatial_1
points_y <- stmetadata2$Spatial_2


# 将封闭曲线和点坐标转为矩阵
polygon_coords <- cbind(polygon_x, polygon_y)
points_coords <- cbind(points_x, points_y)

# 创建 SpatialPolygons对象（封闭曲线多边形）
polygon_sp <- Polygon(polygon_coords)
polygons_sp <- Polygons(list(polygon_sp), "polygon1")
sp_poly <- SpatialPolygons(list(polygons_sp))

# 使用 polygon() 函数绘制多边形
plot(1, 1, type="n", xlim=c(4000, 25000), ylim=c(3000, 25000), xlab="X", ylab="Y")
polygon(polygon_x, polygon_y, col="skyblue", border="black")


# 使用 point.in.polygon 函数检查每个点是否在多边形内
result <- point.in.polygon(points_x, points_y, polygon_x, polygon_y)

# 输出结果：结果为1表示在多边形内，0表示不在多边形内
print(result)
length(result)
stmetadata2$result50 <- result
stmetadata2$result80_1 <- result
stmetadata2$result80_2 <- result
stmetadata2$result95 <- result
stmetadata2$result99 <- result
table(stmetadata2$result50)
table(stmetadata2$result80_2)


stmetadata2$result <- ifelse(stmetadata2$result50==1,"50%",
                             ifelse(stmetadata2$result80_1==1,"80%",
                                    ifelse(stmetadata2$result80_2%in%c(1,2),"80%","95%")))
table(stmetadata2$result)

head(stmetadata2)

table(stmetadata2$result,stmetadata2$group1)
stmetadata2$celltype <- stmeta.data$celltype
table(stmetadata2$result,stmetadata2$celltype)
head(stmeta.data)
head(stmetadata2)
write.csv(stmetadata2,file = "D:/细胞状态2/论文/空间/stad/stad/denggaoxian.csv")

stmetadata2 <- read.csv("D:/细胞状态2/论文/空间/stad/0732stadzhushi/denggaoxian.csv")

ggplot(data = stmetadata2, aes(x = result, y = EMT, fill = result)) +
  geom_violin(position = position_dodge(0.8))+
  geom_boxplot(width=0.5,
               #outlier.colour = NA,
               position = position_dodge(0.8))+  #每组的小提琴图之间的距离
  scale_fill_manual(values = c("#e99593", "#99bca0","steelblue"))+
  theme(axis.text.x = element_text(size = 20,angle = 30, hjust = 0.5, vjust = 0.5))

EMT
Angiogenesis
Invasion
Metastasis

stmetadata3 <- pivot_longer(stmetadata2, 
             cols = c("EMT","Angiogenesis","Invasion","Metastasis","model1"), 
             names_to = "cellstate", 
             values_to = "Value")

stmetadataB <- stmetadata3[stmetadata3$celltype=="B/Myeloid",]
stmetadataend <- stmetadata3[stmetadata3$celltype=="Endothelial",]
stmetadataepi <- stmetadata3[stmetadata3$celltype=="Epithelial",]
stmetadatafib <- stmetadata3[stmetadata3$celltype=="Fibroblasts",]
stmetadatamast <- stmetadata3[stmetadata3$celltype=="Mast",]
stmetadataT <- stmetadata3[stmetadata3$celltype=="T",]



ggplot(data = stmetadataepi, aes(x = cellstate, y = Value, fill = result)) +
  geom_violin(position = position_dodge(0.8))+
  geom_boxplot(width=0.2,
               position = position_dodge(0.8))+  #每组的小提琴图之间的距离
  stat_compare_means(aes(group = result))+ 
  scale_fill_manual(values = c("#e99593", "#99bca0","steelblue"))+
  theme(axis.text.x = element_text(size = 20,angle = 30, hjust = 0.5, vjust = 0.5))+
  theme_minimal()+ 
  theme(panel.grid = element_blank())




#PAAD
data <- data.frame(
  group = c("high", "low"),
  value = c(354, 1164)
)
pie(data$value,labels=data$group,col = c("#c11e21","steelblue"))




data <- data.frame(
  category = c("B", "End", "Epi","Fib","Myeloid","T/NK","B", "End", "Epi","Fib","Myeloid","T/NK","B", "End", "Epi","Fib","Myeloid","T/NK"),
  value = c(0, 9, 670, 198,329,119,1,297,569,88,107,42,309,643,490,15,55,6),
  prob = c("50","50","50","50","50","50","80","80","80","80","80","80","95","95","95","95","95","95")
)

ggplot(data, aes(x = prob, y = value, fill = category)) +
  geom_bar(stat = "identity") +
  labs(title = "Stacked Bar Chart", x = "Category", y = "Value") +
  theme_bw() + 
  theme(panel.grid=element_blank())



#STAD
data <- data.frame(
  group = c("high", "low"),
  value = c(201,1035)
)
pie(data$value,labels=data$group,col = c("#c11e21","steelblue"))



data <- data.frame(
  category = c("B/Myeloid", "End", "Epi","Fib","Mast","T","B/Myeloid", "End", "Epi","Fib","Mast","T","B/Myeloid", "End", "Epi","Fib","Mast","T"),
  value = c(0, 93, 270, 256,273,28,12,122,422,80,104,147,457,23,698,6,2,50),
  prob = c("50","50","50","50","50","50","80","80","80","80","80","80","95","95","95","95","95","95")
)

ggplot(data, aes(x = prob, y = value, fill = category)) +
  geom_bar(stat = "identity") +
  labs(title = "Stacked Bar Chart", x = "Category", y = "Value") +
  theme_bw() + 
  theme(panel.grid=element_blank())









#奇特的基因表达
staddenggaoxian <- read.csv("D:/细胞状态2/论文/空间/stad/stad/denggaoxian.csv")
head(staddenggaoxian)
table(staddenggaoxian$group1)
table(staddenggaoxian$result)
table(staddenggaoxian$group1,staddenggaoxian$result)


load("D:/细胞状态2/论文/空间/stad/stadseurat_data.RData")
load("D:/细胞状态2/论文/空间/stad/stadEST.RData")
head(stadseurat_data)
stadEST[1:3,]
stadexpression_matrix <- GetAssayData(stadseurat_data, assay = "Spatial", slot = "data")
stadexpression_matrix[1:3,1:3]
stadexpression_matrix2 <- as.data.frame(t(stadexpression_matrix))
stadexpression_matrix2[1:3,1:3]
stadexpression_matrix2$cell <- rownames(stadexpression_matrix2)

stadexpression_matrix3 <- merge(staddenggaoxian,stadexpression_matrix2,by="cell")
dim(stadexpression_matrix2)
dim(stadexpression_matrix3)
dim(staddenggaoxian)

i <- "CXCL12"
midmore <- data.frame(row.names = "test")
midmore$highmean <- 0
midmore$midmean <- 0
midmore$lowmean <- 0
midmore$FC1 <- 0
midmore$FC2 <- 0

for (i in rownames(stadexpression_matrix)) {
  highmean<- mean(stadexpression_matrix3[stadexpression_matrix3$result=="50%",i])
  midmean <- mean(stadexpression_matrix3[stadexpression_matrix3$result=="80%",i])
  lowmean <- mean(stadexpression_matrix3[stadexpression_matrix3$result=="95%",i])
  if(midmean>lowmean&midmean>highmean){
    print(i)
    genemore <- data.frame(row.names = i)
    genemore$highmean <- highmean
    genemore$midmean <- midmean
    genemore$lowmean <- lowmean
    genemore$FC1 <- midmean/highmean
    genemore$FC2 <- midmean/lowmean
    midmore <- rbind(midmore,genemore)
  }
}

dim(midmore)
head(midmore)

midmore["CXCL12",]
dim(midmore[midmore$FC1>2& midmore$midmean>0.1,])
dim(midmore[midmore$FC2>2& midmore$midmean>0.1,])
dim(midmore[midmore$FC2>2 & midmore$FC1>2& midmore$midmean>0.1,])
head(midmore[midmore$FC2>2 & midmore$FC1>2& midmore$midmean>0.1,])
head(midmore[midmore$FC1>2& midmore$midmean>0.1,])
midhigh <- midmore[midmore$FC2>2 & midmore$FC1>2& midmore$midmean>0.1,]


stadexpression_matrix3[1:3,1:40]
rownames(stadexpression_matrix3) <- stadexpression_matrix3$cell
midhighexp<- stadexpression_matrix3[,rownames(midhigh)]
head(midhighexp)
write.csv(midhighexp,file = "D:/细胞状态2/论文/挖掘/中间高/midhighexp.csv")


library(corrplot)
singlemt <- midhighexp
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
my_colors <- colorRampPalette(c("steelblue", "white", "#c11e21"))(100)
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)

cor_pearson2 <- as.data.frame(cor_pearson)
cor_pearson2$gene1 <- rownames(cor_pearson2)

library(tidyr)
cor_pearsonlong <- cor_pearson2 %>%
  pivot_longer(
    cols = -gene1, # 指定要转换的列（以 "score" 开头的列）
    names_to = "gene2",     # 新生成的列，表示原列名
    values_to = "corvalue"    # 新生成的列，表示原列的值
  )
head(cor_pearsonlong)

cor_pearsonlong2 <- cor_pearsonlong[cor_pearsonlong$corvalue>0.6&cor_pearsonlong$corvalue!=1,]
unique(cor_pearsonlong2$gene1)








#趋化因子
quhuayinzi <- stadexpression_matrix3[,c("CCL1","CCL2","CCL3","CCL3L1","CCL4","CCL4L1",
                                        "CCL5","CCL7","CCL8","CCL11","CCL13","CCL14",
                                        "CCL15","CCL16","CCL17","CCL18","CCL19","CCL20",
                                        "CCL21","CCL22","CCL23","CCL24","CCL25","CCL26",
                                        "CCL27","CCL28","CXCL1","CXCL2","CXCL3","CXCL4",
                                        "CXCL4L1","CXCL5","CXCL6","CXCL7","CXCL8","CXCL9",
                                        "CXCL10","CXCL11","CXCL12","CXCL13","CXCL14","CXCL16",
                                        "CXCL17","CX3CL1","XCL1","XCL2")]

all_columns <- colnames(stadexpression_matrix3)

# 筛选出所有趋化因子相关列（包含"CCL"或"CXCL"）
chemokine_cols <- grep("CCL|CXCL|CX3CL|XCL", all_columns, value = TRUE, ignore.case = TRUE)

# 打印可用列
print(chemokine_cols)
quhuayinzi <- stadexpression_matrix3[,chemokine_cols]
quhuayinzi$CCL24 <- NULL
library(corrplot)
singlemt <- quhuayinzi
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
my_colors <- colorRampPalette(c("steelblue", "white", "#c11e21"))(100)
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)


#趋化因子下游共表达
quhuayinzi2 <- stadexpression_matrix3[,c("CXCL6","CXCL8","ACKR1","CXCR1","CXCR2")]
head(quhuayinzi2)
singlemt <- quhuayinzi2
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
my_colors <- colorRampPalette(c("steelblue", "white", "#c11e21"))(100)
corrplot(cor_pearson, method = 'number', number.cex = 0.8, diag = FALSE, tl.cex = 0.8,col = my_colors)
corrplot(cor_pearson, add = TRUE, type = 'upper', method = 'pie', diag = FALSE, tl.pos = 'n', cl.pos = 'n',col = my_colors)


quhuayinzi3 <- stadexpression_matrix3[,c("CXCL13","ACKR1","ACKR4","CXCR3","CXCR5")]
head(quhuayinzi3)
singlemt <- quhuayinzi3
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
my_colors <- colorRampPalette(c("steelblue", "white", "#c11e21"))(100)
corrplot(cor_pearson, method = 'number', number.cex = 0.8, diag = FALSE, tl.cex = 0.8,col = my_colors)
corrplot(cor_pearson, add = TRUE, type = 'upper', method = 'pie', diag = FALSE, tl.pos = 'n', cl.pos = 'n',col = my_colors)





#middown基因二维图
load("/wangzhenzenlab/fsh/xff/stad/stadseurat_data.RData")
genes <- c("CXCL8","CXCL6","MMP3",
           "FCRLA","CR2","FDCSP","CXCL13","MS4A1","TCL1A","VPREB3",
           "PRSS21","DPEP3","COL9A3","C20orf204")
pdf("/wangzhenzenlab/xff/cellstate/xff/cellstate2/wajue/gene.pdf")
SpatialFeaturePlot(stadseurat_data, features = genes)
dev.off()

head(stadseurat_data)

pdf("/wangzhenzenlab/xff/cellstate/xff/cellstate2/wajue/model1.pdf")
SpatialFeaturePlot(stadseurat_data, features = "model1")
SpatialFeaturePlot(stadseurat_data, features = c("EMT","Invasion","Metastasis","Angiogenesis"))

dev.off()




#非0区域共表达情况
fei0mat <- stadexpression_matrix3[,c("model1","EMT","Angiogenesis","Invasion","Metastasis",
                                     "CXCL8","CXCL6","MMP3",
                                     "FCRLA","CR2","FDCSP","CXCL13","MS4A1","TCL1A","VPREB3",
                                     "PRSS21","DPEP3","COL9A3","C20orf204")]
head(fei0mat)
colnames(fei0mat)[-1:-5]

i <- "XCL1"
my_colors <- colorRampPalette(c("steelblue", "white", "#c11e21"))(100)
for (i in colnames(fei0mat)[-1:-5]) {
  fei0data <- fei0mat[,c(i,"model1","EMT","Angiogenesis","Invasion","Metastasis")]
  fei0area <- fei0data[fei0data[,i]!=0,]
  singlemt <- fei0area
  singlemt <- scale(singlemt)
  cor_pearson <- cor(singlemt, method = 'pearson')
  myfilename <- paste("D:/细胞状态2/论文/挖掘/非0区域相关性/",i,".pdf",sep = "")
  pdf(myfilename)
  corrplot(cor_pearson, method = 'number', number.cex = 0.8, diag = FALSE, tl.cex = 0.8,col = my_colors)
  corrplot(cor_pearson, add = TRUE, type = 'upper', method = 'pie', diag = FALSE, tl.pos = 'n', cl.pos = 'n',col = my_colors)
  dev.off()
}

head(fei0data)
head(fei0area)

stadexpression_matrix3[,"EMT"]

for (i in colnames(fei0mat)[-1:-5]) {
  fei0data <- fei0mat[,c(i,"model1","EMT","Angiogenesis","Invasion","Metastasis")]
  fei0area <- fei0data[fei0data[,i]>=0.5,]
  singlemt <- fei0area
  singlemt <- scale(singlemt)
  cor_pearson <- cor(singlemt, method = 'pearson')
  myfilename <- paste("D:/细胞状态2/论文/挖掘/非0区域相关性/大于0.1/",i,".pdf",sep = "")
  pdf(myfilename)
  corrplot(cor_pearson, method = 'number', number.cex = 0.8, diag = FALSE, tl.cex = 0.8,col = my_colors)
  corrplot(cor_pearson, add = TRUE, type = 'upper', method = 'pie', diag = FALSE, tl.pos = 'n', cl.pos = 'n',col = my_colors)
  dev.off()
}




#cox
library(stringr)
library(survival)
library(survminer)
setwd("D:/细胞状态2/论文/分化/fenhua/cox")
#Ensembl_ID转基因数据
genedata <- read.table("mart_export.txt", sep=",", header=TRUE)
head(genedata)
#临床数据
clin_all <- read.table("D:/细胞状态2/bulk/clinical/cli_all.txt",header = T)
head(clin_all)


stadbulk <- read.table("D:/细胞状态2/bulk/expression-counts/TCGA-STAD.htseq_counts.tsv", sep="\t",header=TRUE)
stadbulk[1:3,1:3]

stadbulk$Ensembl_ID <- str_split_fixed(stadbulk$Ensembl_ID,'[.]',2)[,1]
stadnew_data <- merge(stadbulk,genedata,by.x = "Ensembl_ID",by.y = "Gene.stable.ID")
dim(stadnew_data)
col<- c("Ensembl_ID", "Gene.stable.ID.version", "Transcript.stable.ID.version", "Transcript.stable.ID") 
stadnew_df <- stadnew_data[, !(colnames(stadnew_data)%in%col)]
stadunique_df <- stadnew_df[!duplicated(stadnew_df$Gene.name),]
rownames(stadunique_df) <- stadunique_df[,"Gene.name"]
stadunique_df$Gene.name <- NULL

stadunique_df[1:3,1:3]

# 进行log转化，避免零值影响
log_stadcounts <- log2(stadunique_df + 1)
log_stadcounts[1:3,1:3]
max(log_stadcounts)


log_stadcountsT <- t(log_stadcounts)
log_stadcountsT[1:3,1:3]
log_stadcountsT <- as.data.frame(log_stadcountsT)
log_stadcountsT$sample <- rownames(log_stadcountsT)
log_stadcountsT$patient <- sub("\\.[^\\.]*$", "", log_stadcountsT$sample)
head(log_stadcountsT$patient)


clin_stad <- merge(log_stadcountsT,clin_all,by.x="patient",by.y="PATIENT")
clin_stad[1:3,1:3]


#"CXCL8","CXCL6","MMP3","FCRLA","CR2","FDCSP","CXCL13","MS4A1","TCL1A","VPREB3","PRSS21","DPEP3","COL9A3","C20orf204"
cox_model <- coxph(Surv(OS, EVENT) ~ CXCL8, data = clin_stad)
cox_model <- coxph(Surv(OS, EVENT) ~ CXCL6, data = clin_stad)
cox_model <- coxph(Surv(OS, EVENT) ~ MMP3, data = clin_stad)
cox_model <- coxph(Surv(OS, EVENT) ~ FCRLA, data = clin_stad)
cox_model <- coxph(Surv(OS, EVENT) ~ CR2, data = clin_stad)
cox_model <- coxph(Surv(OS, EVENT) ~ FDCSP, data = clin_stad)
cox_model <- coxph(Surv(OS, EVENT) ~ CXCL13, data = clin_stad)
cox_model <- coxph(Surv(OS, EVENT) ~ MS4A1, data = clin_stad)
cox_model <- coxph(Surv(OS, EVENT) ~ TCL1A, data = clin_stad)
cox_model <- coxph(Surv(OS, EVENT) ~ VPREB3, data = clin_stad)
cox_model <- coxph(Surv(OS, EVENT) ~ PRSS21, data = clin_stad)
cox_model <- coxph(Surv(OS, EVENT) ~ DPEP3, data = clin_stad)
cox_model <- coxph(Surv(OS, EVENT) ~ COL9A3, data = clin_stad)
cox_model <- coxph(Surv(OS, EVENT) ~ C20orf204, data = clin_stad)

summary(cox_model)
# 绘制森林图
ggforest(cox_model, data = clin_stad)


1






































#独立验证
load("D:/细胞状态2/论文/空间/stad/stadexp.RData")
load("D:/细胞状态2/论文/空间/stad/stadseurat_data.RData")
load("D:/细胞状态2/论文/空间/stad/stmetadata.RData")
load("D:/细胞状态2/论文/空间/stad/ES.RData")

position0331=read.csv("./20_00331_LI_SING/spatial/tissue_positions_list.csv",header=F,row.names=1)
position0731=read.csv("./21_00731_LI_SING/spatial/tissue_positions_list.csv",header=F,row.names=1)
position0733=read.csv("./21_00733_LI_SING/spatial/tissue_positions_list.csv",header=F,row.names=1)
position0734=read.csv("./21_00734_LI_SING/spatial/tissue_positions_list.csv",header=F,row.names=1)
position1251=read.csv("./21_01251_LI_SING/spatial/tissue_positions_list.csv",header=F,row.names=1)
position1252=read.csv("./21_01252_LI_SING/spatial/tissue_positions_list.csv",header=F,row.names=1)
position1253=read.csv("./21_01253_LI_SING/spatial/tissue_positions_list.csv",header=F,row.names=1)
position1254=read.csv("./21_01254_LI_SING/spatial/tissue_positions_list.csv",header=F,row.names=1)
position1675=read.csv("./21_01675_LI_SING/spatial/tissue_positions_list.csv",header=F,row.names=1)

position0331$V1 <- rownames(position0331)
position0731$V1 <- rownames(position0731)
position0733$V1 <- rownames(position0733)
position0734$V1 <- rownames(position0734)
position1251$V1 <- rownames(position1251)
position1252$V1 <- rownames(position1252)
position1253$V1 <- rownames(position1253)
position1254$V1 <- rownames(position1254)
position1675$V1 <- rownames(position1675)

head(stmetadata0331)
head(position0331)

stmetadata2_0331 <- merge(stmetadata0331, position0331, by.x = "cell", by.y = "V1")
stmetadata2_0731 <- merge(stmetadata0731, position0731, by.x = "cell", by.y = "V1")
stmetadata2_0733 <- merge(stmetadata0733, position0733, by.x = "cell", by.y = "V1")
stmetadata2_0734 <- merge(stmetadata0734, position0734, by.x = "cell", by.y = "V1")
stmetadata2_1251 <- merge(stmetadata1251, position1251, by.x = "cell", by.y = "V1")
stmetadata2_1252 <- merge(stmetadata1252, position1252, by.x = "cell", by.y = "V1")
stmetadata2_1253 <- merge(stmetadata1253, position1253, by.x = "cell", by.y = "V1")
stmetadata2_1254 <- merge(stmetadata1254, position1254, by.x = "cell", by.y = "V1")
stmetadata2_1675 <- merge(stmetadata1675, position1675, by.x = "cell", by.y = "V1")


head(stmetadata2_0331)


stmetadata2_0331$group1 <- as.factor(stmetadata2_0331$group1)
stmetadata2_0331$group1 <- factor(stmetadata2_0331$group1, levels = c("low", "high"))
stmetadata2_0731$group1 <- as.factor(stmetadata2_0731$group1)
stmetadata2_0731$group1 <- factor(stmetadata2_0731$group1, levels = c("low", "high"))
stmetadata2_0733$group1 <- as.factor(stmetadata2_0733$group1)
stmetadata2_0733$group1 <- factor(stmetadata2_0733$group1, levels = c("low", "high"))
stmetadata2_0734$group1 <- as.factor(stmetadata2_0734$group1)
stmetadata2_0734$group1 <- factor(stmetadata2_0734$group1, levels = c("low", "high"))
stmetadata2_1251$group1 <- as.factor(stmetadata2_1251$group1)
stmetadata2_1251$group1 <- factor(stmetadata2_1251$group1, levels = c("low", "high"))
stmetadata2_1252$group1 <- as.factor(stmetadata2_1252$group1)
stmetadata2_1252$group1 <- factor(stmetadata2_1252$group1, levels = c("low", "high"))
stmetadata2_1253$group1 <- as.factor(stmetadata2_1253$group1)
stmetadata2_1253$group1 <- factor(stmetadata2_1253$group1, levels = c("low", "high"))
stmetadata2_1254$group1 <- as.factor(stmetadata2_1254$group1)
stmetadata2_1254$group1 <- factor(stmetadata2_1254$group1, levels = c("low", "high"))
stmetadata2_1675$group1 <- as.factor(stmetadata2_1675$group1)
stmetadata2_1675$group1 <- factor(stmetadata2_1675$group1, levels = c("low", "high"))


# 0331
denggaoxian0331 <- ggplot(stmetadata2_0331,aes(x = V3, y = V4, fill = group1)) +
  geom_hdr() +
  scale_fill_aaas() +
  facet_wrap(vars(group1)) +
  theme_bw() +
  theme(panel.grid = element_blank())

# 0731
denggaoxian0731 <- ggplot(stmetadata2_0731,aes(x = V3, y = V4, fill = group1)) +
  geom_hdr() +
  scale_fill_aaas() +
  facet_wrap(vars(group1)) +
  theme_bw() +
  theme(panel.grid = element_blank())

# 0733
denggaoxian0733 <- ggplot(stmetadata2_0733,aes(x = V3, y = V4, fill = group1)) +
  geom_hdr() +
  scale_fill_aaas() +
  facet_wrap(vars(group1)) +
  theme_bw() +
  theme(panel.grid = element_blank())

# 0734
denggaoxian0734 <- ggplot(stmetadata2_0734,aes(x = V3, y = V4, fill = group1)) +
  geom_hdr() +
  scale_fill_aaas() +
  facet_wrap(vars(group1)) +
  theme_bw() +
  theme(panel.grid = element_blank())

# 1251
denggaoxian1251 <- ggplot(stmetadata2_1251,aes(x = V3, y = V4, fill = group1)) +
  geom_hdr() +
  scale_fill_aaas() +
  facet_wrap(vars(group1)) +
  theme_bw() +
  theme(panel.grid = element_blank())

# 1252
denggaoxian1252 <- ggplot(stmetadata2_1252,aes(x = V3, y = V4, fill = group1)) +
  geom_hdr() +
  scale_fill_aaas() +
  facet_wrap(vars(group1)) +
  theme_bw() +
  theme(panel.grid = element_blank())

# 1253
denggaoxian1253 <- ggplot(stmetadata2_1253,aes(x = V3, y = V4, fill = group1)) +
  geom_hdr() +
  scale_fill_aaas() +
  facet_wrap(vars(group1)) +
  theme_bw() +
  theme(panel.grid = element_blank())

# 1254
denggaoxian1254 <- ggplot(stmetadata2_1254,aes(x = V3, y = V4, fill = group1)) +
  geom_hdr() +
  scale_fill_aaas() +
  facet_wrap(vars(group1)) +
  theme_bw() +
  theme(panel.grid = element_blank())

# 1675
denggaoxian1675 <- ggplot(stmetadata2_1675,aes(x = V3, y = V4, fill = group1)) +
  geom_hdr() +
  scale_fill_aaas() +
  facet_wrap(vars(group1)) +
  theme_bw() +
  theme(panel.grid = element_blank())


ggplot(stmetadata2_0331,aes(x = V3, y = V4, fill = group1))+
  geom_point(shape = 21) +
  scale_fill_manual(values = c("high"="#c11e21","low"="steelblue"))



#0331给点划分区域
plot_data0331 <- ggplot_build(denggaoxian0331)$data[[1]]  
prob50_spotshigh0331 <- plot_data0331[plot_data0331$group == "2-004", ]
prob80_spotshigh0331 <- plot_data0331[plot_data0331$group == "2-003", ]
prob95_spotshigh0331 <- plot_data0331[plot_data0331$group == "2-002", ]
prob99_spotshigh0331 <- plot_data0331[plot_data0331$group == "2-001", ]


max(prob99_spotshigh0331$x)
min(prob99_spotshigh0331$x)
max(prob99_spotshigh0331$y)
min(prob99_spotshigh0331$y)


# 封闭曲线（多边形）坐标
polygon_x0331_50 <- prob50_spotshigh0331$x  
polygon_y0331_50 <- prob50_spotshigh0331$y  
polygon_x0331_80 <- prob80_spotshigh0331$x
polygon_y0331_80 <- prob80_spotshigh0331$y
polygon_x0331_95 <- prob95_spotshigh0331$x
polygon_y0331_95 <- prob95_spotshigh0331$y
polygon_x0331_99 <- prob99_spotshigh0331$x
polygon_y0331_99 <- prob99_spotshigh0331$y

# 待检测点的坐标
points_x0331 <- stmetadata2_0331$V3
points_y0331 <- stmetadata2_0331$V4

# 将封闭曲线和点坐标转为矩阵
polygon_coords0331_50 <- cbind(polygon_x0331_50, polygon_y0331_50)
polygon_coords0331_80 <- cbind(polygon_x0331_80, polygon_y0331_80)
polygon_coords0331_95 <- cbind(polygon_x0331_95, polygon_y0331_95)
polygon_coords0331_99 <- cbind(polygon_x0331_99, polygon_y0331_99)
points_coords0331 <- cbind(points_x0331, points_y0331)

# 创建 SpatialPolygons对象（封闭曲线多边形）
polygon_sp0331_50 <- Polygon(polygon_coords0331_50)
polygon_sp0331_80 <- Polygon(polygon_coords0331_80)
polygon_sp0331_95 <- Polygon(polygon_coords0331_95)
polygon_sp0331_99 <- Polygon(polygon_coords0331_99)
polygons_sp0331_50 <- Polygons(list(polygon_sp0331_50), "polygon1")
polygons_sp0331_80 <- Polygons(list(polygon_sp0331_80), "polygon1")
polygons_sp0331_95 <- Polygons(list(polygon_sp0331_95), "polygon1")
polygons_sp0331_99 <- Polygons(list(polygon_sp0331_99), "polygon1")
sp_poly0331_50 <- SpatialPolygons(list(polygons_sp0331_50))
sp_poly0331_80 <- SpatialPolygons(list(polygons_sp0331_80))
sp_poly0331_95 <- SpatialPolygons(list(polygons_sp0331_95))
sp_poly0331_99 <- SpatialPolygons(list(polygons_sp0331_99))
# 使用 polygon() 函数绘制多边形
plot(1, 1, type="n", xlim=c(0, 130), ylim=c(0, 130), xlab="X", ylab="Y")
polygon(polygon_x0331_50, polygon_y0331_50, col="skyblue", border="black")
polygon(polygon_x0331_80, polygon_y0331_80, col="skyblue", border="black")
polygon(polygon_x0331_95, polygon_y0331_95, col="skyblue", border="black")
polygon(polygon_x0331_99, polygon_y0331_99, col="skyblue", border="black")

# 使用 point.in.polygon 函数检查每个点是否在多边形内
result0331_50 <- point.in.polygon(points_x0331, points_y0331, polygon_x0331_50, polygon_y0331_50)
result0331_80 <- point.in.polygon(points_x0331, points_y0331, polygon_x0331_80, polygon_y0331_80)
result0331_95 <- point.in.polygon(points_x0331, points_y0331, polygon_x0331_95, polygon_y0331_95)
result0331_99 <- point.in.polygon(points_x0331, points_y0331, polygon_x0331_99, polygon_y0331_99)

# 输出结果：结果为1表示在多边形内，0表示不在多边形内
stmetadata2_0331$result50 <- result0331_50
stmetadata2_0331$result80 <- result0331_80
stmetadata2_0331$result95 <- result0331_95
stmetadata2_0331$result99 <- result0331_99


stmetadata2_0331$result <- ifelse(stmetadata2_0331$result50==1,"50%",
                             ifelse(stmetadata2_0331$result80==1,"80%","95%"))
table(stmetadata2_0331$result)
write.csv(stmetadata2_0331,file = "D:/细胞状态2/论文/挖掘/独立验证集/等高线/denggaoxian0331.csv")


ggplot(data = stmetadata2_0331, aes(x = result, y = model1, fill = result)) +
  geom_violin(position = position_dodge(0.8))+
  geom_boxplot(width=0.5,
               #outlier.colour = NA,
               position = position_dodge(0.8))+  #每组的小提琴图之间的距离
  scale_fill_manual(values = c("#e99593", "#99bca0","steelblue"))+
  theme(axis.text.x = element_text(size = 20,angle = 30, hjust = 0.5, vjust = 0.5))

EMT
Angiogenesis
Invasion
Metastasis









#0731给点划分区域
plot_data0731 <- ggplot_build(denggaoxian0731)$data[[1]]  
prob50_spotshigh0731 <- plot_data0731[plot_data0731$group == "2-004", ]
prob80_spotshigh0731 <- plot_data0731[plot_data0731$group == "2-003", ]
prob95_spotshigh0731 <- plot_data0731[plot_data0731$group == "2-002", ]
prob99_spotshigh0731 <- plot_data0731[plot_data0731$group == "2-001", ]


max(prob99_spotshigh0731$x)
min(prob99_spotshigh0731$x)
max(prob99_spotshigh0731$y)
min(prob99_spotshigh0731$y)


# 封闭曲线（多边形）坐标
polygon_x0731_50 <- prob50_spotshigh0731$x  
polygon_y0731_50 <- prob50_spotshigh0731$y  
polygon_x0731_80 <- prob80_spotshigh0731$x
polygon_y0731_80 <- prob80_spotshigh0731$y
polygon_x0731_95 <- prob95_spotshigh0731$x
polygon_y0731_95 <- prob95_spotshigh0731$y
polygon_x0731_99 <- prob99_spotshigh0731$x
polygon_y0731_99 <- prob99_spotshigh0731$y

# 待检测点的坐标
points_x0731 <- stmetadata2_0731$V3
points_y0731 <- stmetadata2_0731$V4

# 将封闭曲线和点坐标转为矩阵
polygon_coords0731_50 <- cbind(polygon_x0731_50, polygon_y0731_50)
polygon_coords0731_80 <- cbind(polygon_x0731_80, polygon_y0731_80)
polygon_coords0731_95 <- cbind(polygon_x0731_95, polygon_y0731_95)
polygon_coords0731_99 <- cbind(polygon_x0731_99, polygon_y0731_99)
points_coords0731 <- cbind(points_x0731, points_y0731)

# 创建 SpatialPolygons对象（封闭曲线多边形）
polygon_sp0731_50 <- Polygon(polygon_coords0731_50)
polygon_sp0731_80 <- Polygon(polygon_coords0731_80)
polygon_sp0731_95 <- Polygon(polygon_coords0731_95)
polygon_sp0731_99 <- Polygon(polygon_coords0731_99)
polygons_sp0731_50 <- Polygons(list(polygon_sp0731_50), "polygon1")
polygons_sp0731_80 <- Polygons(list(polygon_sp0731_80), "polygon1")
polygons_sp0731_95 <- Polygons(list(polygon_sp0731_95), "polygon1")
polygons_sp0731_99 <- Polygons(list(polygon_sp0731_99), "polygon1")
sp_poly0731_50 <- SpatialPolygons(list(polygons_sp0731_50))
sp_poly0731_80 <- SpatialPolygons(list(polygons_sp0731_80))
sp_poly0731_95 <- SpatialPolygons(list(polygons_sp0731_95))
sp_poly0731_99 <- SpatialPolygons(list(polygons_sp0731_99))
# 使用 polygon() 函数绘制多边形
plot(1, 1, type="n", xlim=c(0, 130), ylim=c(0, 130), xlab="X", ylab="Y")
polygon(polygon_x0731_50, polygon_y0731_50, col="skyblue", border="black")
polygon(polygon_x0731_80, polygon_y0731_80, col="skyblue", border="black")
polygon(polygon_x0731_95, polygon_y0731_95, col="skyblue", border="black")
polygon(polygon_x0731_99, polygon_y0731_99, col="skyblue", border="black")

# 使用 point.in.polygon 函数检查每个点是否在多边形内
result0731_50 <- point.in.polygon(points_x0731, points_y0731, polygon_x0731_50, polygon_y0731_50)
result0731_80 <- point.in.polygon(points_x0731, points_y0731, polygon_x0731_80, polygon_y0731_80)
result0731_95 <- point.in.polygon(points_x0731, points_y0731, polygon_x0731_95, polygon_y0731_95)
result0731_99 <- point.in.polygon(points_x0731, points_y0731, polygon_x0731_99, polygon_y0731_99)

# 输出结果：结果为1表示在多边形内，0表示不在多边形内
stmetadata2_0731$result50 <- result0731_50
stmetadata2_0731$result80 <- result0731_80
stmetadata2_0731$result95 <- result0731_95
stmetadata2_0731$result99 <- result0731_99


stmetadata2_0731$result <- ifelse(stmetadata2_0731$result50==1,"50%",
                                  ifelse(stmetadata2_0731$result80==1,"80%","95%"))
table(stmetadata2_0731$result)
write.csv(stmetadata2_0731,file = "D:/细胞状态2/论文/挖掘/独立验证集/等高线/denggaoxian0731.csv")


ggplot(data = stmetadata2_0731, aes(x = result, y = model1, fill = result)) +
  geom_violin(position = position_dodge(0.8))+
  geom_boxplot(width=0.5,
               #outlier.colour = NA,
               position = position_dodge(0.8))+  #每组的小提琴图之间的距离
  scale_fill_manual(values = c("#e99593", "#99bca0","steelblue"))+
  theme(axis.text.x = element_text(size = 20,angle = 30, hjust = 0.5, vjust = 0.5))

EMT
Angiogenesis
Invasion
Metastasis










#0733给点划分区域
plot_data0733 <- ggplot_build(denggaoxian0733)$data[[1]]  
prob50_spotshigh0733 <- plot_data0733[plot_data0733$group == "2-004", ]
prob80_spotshigh0733 <- plot_data0733[plot_data0733$group == "2-003", ]
prob95_spotshigh0733 <- plot_data0733[plot_data0733$group == "2-002", ]
prob99_spotshigh0733 <- plot_data0733[plot_data0733$group == "2-001", ]


max(prob99_spotshigh0733$x)
min(prob99_spotshigh0733$x)
max(prob99_spotshigh0733$y)
min(prob99_spotshigh0733$y)


# 封闭曲线（多边形）坐标
polygon_x0733_50 <- prob50_spotshigh0733$x  
polygon_y0733_50 <- prob50_spotshigh0733$y  
polygon_x0733_80 <- prob80_spotshigh0733$x
polygon_y0733_80 <- prob80_spotshigh0733$y
polygon_x0733_95 <- prob95_spotshigh0733$x
polygon_y0733_95 <- prob95_spotshigh0733$y
polygon_x0733_99 <- prob99_spotshigh0733$x
polygon_y0733_99 <- prob99_spotshigh0733$y

# 待检测点的坐标
points_x0733 <- stmetadata2_0733$V3
points_y0733 <- stmetadata2_0733$V4

# 将封闭曲线和点坐标转为矩阵
polygon_coords0733_50 <- cbind(polygon_x0733_50, polygon_y0733_50)
polygon_coords0733_80 <- cbind(polygon_x0733_80, polygon_y0733_80)
polygon_coords0733_95 <- cbind(polygon_x0733_95, polygon_y0733_95)
polygon_coords0733_99 <- cbind(polygon_x0733_99, polygon_y0733_99)
points_coords0733 <- cbind(points_x0733, points_y0733)

# 创建 SpatialPolygons对象（封闭曲线多边形）
polygon_sp0733_50 <- Polygon(polygon_coords0733_50)
polygon_sp0733_80 <- Polygon(polygon_coords0733_80)
polygon_sp0733_95 <- Polygon(polygon_coords0733_95)
polygon_sp0733_99 <- Polygon(polygon_coords0733_99)
polygons_sp0733_50 <- Polygons(list(polygon_sp0733_50), "polygon1")
polygons_sp0733_80 <- Polygons(list(polygon_sp0733_80), "polygon1")
polygons_sp0733_95 <- Polygons(list(polygon_sp0733_95), "polygon1")
polygons_sp0733_99 <- Polygons(list(polygon_sp0733_99), "polygon1")
sp_poly0733_50 <- SpatialPolygons(list(polygons_sp0733_50))
sp_poly0733_80 <- SpatialPolygons(list(polygons_sp0733_80))
sp_poly0733_95 <- SpatialPolygons(list(polygons_sp0733_95))
sp_poly0733_99 <- SpatialPolygons(list(polygons_sp0733_99))
# 使用 polygon() 函数绘制多边形
plot(1, 1, type="n", xlim=c(0, 130), ylim=c(0, 130), xlab="X", ylab="Y")
polygon(polygon_x0733_50, polygon_y0733_50, col="skyblue", border="black")
polygon(polygon_x0733_80, polygon_y0733_80, col="skyblue", border="black")
polygon(polygon_x0733_95, polygon_y0733_95, col="skyblue", border="black")
polygon(polygon_x0733_99, polygon_y0733_99, col="skyblue", border="black")

# 使用 point.in.polygon 函数检查每个点是否在多边形内
result0733_50 <- point.in.polygon(points_x0733, points_y0733, polygon_x0733_50, polygon_y0733_50)
result0733_80 <- point.in.polygon(points_x0733, points_y0733, polygon_x0733_80, polygon_y0733_80)
result0733_95 <- point.in.polygon(points_x0733, points_y0733, polygon_x0733_95, polygon_y0733_95)
result0733_99 <- point.in.polygon(points_x0733, points_y0733, polygon_x0733_99, polygon_y0733_99)

# 输出结果：结果为1表示在多边形内，0表示不在多边形内
stmetadata2_0733$result50 <- result0733_50
stmetadata2_0733$result80 <- result0733_80
stmetadata2_0733$result95 <- result0733_95
stmetadata2_0733$result99 <- result0733_99


stmetadata2_0733$result <- ifelse(stmetadata2_0733$result50==1,"50%",
                                  ifelse(stmetadata2_0733$result80==1,"80%","95%"))
table(stmetadata2_0733$result)
write.csv(stmetadata2_0733,file = "D:/细胞状态2/论文/挖掘/独立验证集/等高线/denggaoxian0733.csv")


ggplot(data = stmetadata2_0733, aes(x = result, y = model1, fill = result)) +
  geom_violin(position = position_dodge(0.8))+
  geom_boxplot(width=0.5,
               #outlier.colour = NA,
               position = position_dodge(0.8))+  #每组的小提琴图之间的距离
  scale_fill_manual(values = c("#e99593", "#99bca0","steelblue"))+
  theme(axis.text.x = element_text(size = 20,angle = 30, hjust = 0.5, vjust = 0.5))

EMT
Angiogenesis
Invasion
Metastasis











#0734给点划分区域
plot_data0734 <- ggplot_build(denggaoxian0734)$data[[1]]  
prob50_spotshigh0734 <- plot_data0734[plot_data0734$group == "2-004", ]
prob80_spotshigh0734 <- plot_data0734[plot_data0734$group == "2-003", ]
prob95_spotshigh0734 <- plot_data0734[plot_data0734$group == "2-002", ]
prob99_spotshigh0734 <- plot_data0734[plot_data0734$group == "2-001", ]


max(prob99_spotshigh0734$x)
min(prob99_spotshigh0734$x)
max(prob99_spotshigh0734$y)
min(prob99_spotshigh0734$y)


# 封闭曲线（多边形）坐标
polygon_x0734_50 <- prob50_spotshigh0734$x  
polygon_y0734_50 <- prob50_spotshigh0734$y  
polygon_x0734_80 <- prob80_spotshigh0734$x
polygon_y0734_80 <- prob80_spotshigh0734$y
polygon_x0734_95 <- prob95_spotshigh0734$x
polygon_y0734_95 <- prob95_spotshigh0734$y
polygon_x0734_99 <- prob99_spotshigh0734$x
polygon_y0734_99 <- prob99_spotshigh0734$y

# 待检测点的坐标
points_x0734 <- stmetadata2_0734$V3
points_y0734 <- stmetadata2_0734$V4

# 将封闭曲线和点坐标转为矩阵
polygon_coords0734_50 <- cbind(polygon_x0734_50, polygon_y0734_50)
polygon_coords0734_80 <- cbind(polygon_x0734_80, polygon_y0734_80)
polygon_coords0734_95 <- cbind(polygon_x0734_95, polygon_y0734_95)
polygon_coords0734_99 <- cbind(polygon_x0734_99, polygon_y0734_99)
points_coords0734 <- cbind(points_x0734, points_y0734)

# 创建 SpatialPolygons对象（封闭曲线多边形）
polygon_sp0734_50 <- Polygon(polygon_coords0734_50)
polygon_sp0734_80 <- Polygon(polygon_coords0734_80)
polygon_sp0734_95 <- Polygon(polygon_coords0734_95)
polygon_sp0734_99 <- Polygon(polygon_coords0734_99)
polygons_sp0734_50 <- Polygons(list(polygon_sp0734_50), "polygon1")
polygons_sp0734_80 <- Polygons(list(polygon_sp0734_80), "polygon1")
polygons_sp0734_95 <- Polygons(list(polygon_sp0734_95), "polygon1")
polygons_sp0734_99 <- Polygons(list(polygon_sp0734_99), "polygon1")
sp_poly0734_50 <- SpatialPolygons(list(polygons_sp0734_50))
sp_poly0734_80 <- SpatialPolygons(list(polygons_sp0734_80))
sp_poly0734_95 <- SpatialPolygons(list(polygons_sp0734_95))
sp_poly0734_99 <- SpatialPolygons(list(polygons_sp0734_99))
# 使用 polygon() 函数绘制多边形
plot(1, 1, type="n", xlim=c(0, 130), ylim=c(0, 130), xlab="X", ylab="Y")
polygon(polygon_x0734_50, polygon_y0734_50, col="skyblue", border="black")
polygon(polygon_x0734_80, polygon_y0734_80, col="skyblue", border="black")
polygon(polygon_x0734_95, polygon_y0734_95, col="skyblue", border="black")
polygon(polygon_x0734_99, polygon_y0734_99, col="skyblue", border="black")

# 使用 point.in.polygon 函数检查每个点是否在多边形内
result0734_50 <- point.in.polygon(points_x0734, points_y0734, polygon_x0734_50, polygon_y0734_50)
result0734_80 <- point.in.polygon(points_x0734, points_y0734, polygon_x0734_80, polygon_y0734_80)
result0734_95 <- point.in.polygon(points_x0734, points_y0734, polygon_x0734_95, polygon_y0734_95)
result0734_99 <- point.in.polygon(points_x0734, points_y0734, polygon_x0734_99, polygon_y0734_99)

# 输出结果：结果为1表示在多边形内，0表示不在多边形内
stmetadata2_0734$result50 <- result0734_50
stmetadata2_0734$result80 <- result0734_80
stmetadata2_0734$result95 <- result0734_95
stmetadata2_0734$result99 <- result0734_99


stmetadata2_0734$result <- ifelse(stmetadata2_0734$result50==1,"50%",
                                  ifelse(stmetadata2_0734$result80==1,"80%","95%"))
table(stmetadata2_0734$result)
write.csv(stmetadata2_0734,file = "D:/细胞状态2/论文/挖掘/独立验证集/等高线/denggaoxian0734.csv")


ggplot(data = stmetadata2_0734, aes(x = result, y = model1, fill = result)) +
  geom_violin(position = position_dodge(0.8))+
  geom_boxplot(width=0.5,
               #outlier.colour = NA,
               position = position_dodge(0.8))+  #每组的小提琴图之间的距离
  scale_fill_manual(values = c("#e99593", "#99bca0","steelblue"))+
  theme(axis.text.x = element_text(size = 20,angle = 30, hjust = 0.5, vjust = 0.5))

EMT
Angiogenesis
Invasion
Metastasis













#1251给点划分区域
plot_data1251 <- ggplot_build(denggaoxian1251)$data[[1]]  
prob50_spotshigh1251 <- plot_data1251[plot_data1251$group == "2-004", ]
prob80_spotshigh1251 <- plot_data1251[plot_data1251$group == "2-003", ]
prob95_spotshigh1251 <- plot_data1251[plot_data1251$group == "2-002", ]
prob99_spotshigh1251 <- plot_data1251[plot_data1251$group == "2-001", ]


max(prob99_spotshigh1251$x)
min(prob99_spotshigh1251$x)
max(prob99_spotshigh1251$y)
min(prob99_spotshigh1251$y)


# 封闭曲线（多边形）坐标
polygon_x1251_50 <- prob50_spotshigh1251$x  
polygon_y1251_50 <- prob50_spotshigh1251$y  
polygon_x1251_80 <- prob80_spotshigh1251$x
polygon_y1251_80 <- prob80_spotshigh1251$y
polygon_x1251_95 <- prob95_spotshigh1251$x
polygon_y1251_95 <- prob95_spotshigh1251$y
polygon_x1251_99 <- prob99_spotshigh1251$x
polygon_y1251_99 <- prob99_spotshigh1251$y

# 待检测点的坐标
points_x1251 <- stmetadata2_1251$V3
points_y1251 <- stmetadata2_1251$V4

# 将封闭曲线和点坐标转为矩阵
polygon_coords1251_50 <- cbind(polygon_x1251_50, polygon_y1251_50)
polygon_coords1251_80 <- cbind(polygon_x1251_80, polygon_y1251_80)
polygon_coords1251_95 <- cbind(polygon_x1251_95, polygon_y1251_95)
polygon_coords1251_99 <- cbind(polygon_x1251_99, polygon_y1251_99)
points_coords1251 <- cbind(points_x1251, points_y1251)

# 创建 SpatialPolygons对象（封闭曲线多边形）
polygon_sp1251_50 <- Polygon(polygon_coords1251_50)
polygon_sp1251_80 <- Polygon(polygon_coords1251_80)
polygon_sp1251_95 <- Polygon(polygon_coords1251_95)
polygon_sp1251_99 <- Polygon(polygon_coords1251_99)
polygons_sp1251_50 <- Polygons(list(polygon_sp1251_50), "polygon1")
polygons_sp1251_80 <- Polygons(list(polygon_sp1251_80), "polygon1")
polygons_sp1251_95 <- Polygons(list(polygon_sp1251_95), "polygon1")
polygons_sp1251_99 <- Polygons(list(polygon_sp1251_99), "polygon1")
sp_poly1251_50 <- SpatialPolygons(list(polygons_sp1251_50))
sp_poly1251_80 <- SpatialPolygons(list(polygons_sp1251_80))
sp_poly1251_95 <- SpatialPolygons(list(polygons_sp1251_95))
sp_poly1251_99 <- SpatialPolygons(list(polygons_sp1251_99))
# 使用 polygon() 函数绘制多边形
plot(1, 1, type="n", xlim=c(0, 130), ylim=c(0, 130), xlab="X", ylab="Y")
polygon(polygon_x1251_50, polygon_y1251_50, col="skyblue", border="black")
polygon(polygon_x1251_80, polygon_y1251_80, col="skyblue", border="black")
polygon(polygon_x1251_95, polygon_y1251_95, col="skyblue", border="black")
polygon(polygon_x1251_99, polygon_y1251_99, col="skyblue", border="black")

# 使用 point.in.polygon 函数检查每个点是否在多边形内
result1251_50 <- point.in.polygon(points_x1251, points_y1251, polygon_x1251_50, polygon_y1251_50)
result1251_80 <- point.in.polygon(points_x1251, points_y1251, polygon_x1251_80, polygon_y1251_80)
result1251_95 <- point.in.polygon(points_x1251, points_y1251, polygon_x1251_95, polygon_y1251_95)
result1251_99 <- point.in.polygon(points_x1251, points_y1251, polygon_x1251_99, polygon_y1251_99)

# 输出结果：结果为1表示在多边形内，0表示不在多边形内
stmetadata2_1251$result50 <- result1251_50
stmetadata2_1251$result80 <- result1251_80
stmetadata2_1251$result95 <- result1251_95
stmetadata2_1251$result99 <- result1251_99


stmetadata2_1251$result <- ifelse(stmetadata2_1251$result50==1,"50%",
                                  ifelse(stmetadata2_1251$result80==1,"80%","95%"))
table(stmetadata2_1251$result)
write.csv(stmetadata2_1251,file = "D:/细胞状态2/论文/挖掘/独立验证集/等高线/denggaoxian1251.csv")


ggplot(data = stmetadata2_1251, aes(x = result, y = model1, fill = result)) +
  geom_violin(position = position_dodge(0.8))+
  geom_boxplot(width=0.5,
               #outlier.colour = NA,
               position = position_dodge(0.8))+  #每组的小提琴图之间的距离
  scale_fill_manual(values = c("#e99593", "#99bca0","steelblue"))+
  theme(axis.text.x = element_text(size = 20,angle = 30, hjust = 0.5, vjust = 0.5))

EMT
Angiogenesis
Invasion
Metastasis












#1252给点划分区域
plot_data1252 <- ggplot_build(denggaoxian1252)$data[[1]]  
prob50_spotshigh1252 <- plot_data1252[plot_data1252$group == "2-004", ]
prob80_spotshigh1252 <- plot_data1252[plot_data1252$group == "2-003", ]
prob95_spotshigh1252 <- plot_data1252[plot_data1252$group == "2-002", ]
prob99_spotshigh1252 <- plot_data1252[plot_data1252$group == "2-001", ]


max(prob99_spotshigh1252$x)
min(prob99_spotshigh1252$x)
max(prob99_spotshigh1252$y)
min(prob99_spotshigh1252$y)


# 封闭曲线（多边形）坐标
polygon_x1252_50 <- prob50_spotshigh1252$x  
polygon_y1252_50 <- prob50_spotshigh1252$y  
polygon_x1252_80 <- prob80_spotshigh1252$x
polygon_y1252_80 <- prob80_spotshigh1252$y
polygon_x1252_95 <- prob95_spotshigh1252$x
polygon_y1252_95 <- prob95_spotshigh1252$y
polygon_x1252_99 <- prob99_spotshigh1252$x
polygon_y1252_99 <- prob99_spotshigh1252$y

# 待检测点的坐标
points_x1252 <- stmetadata2_1252$V3
points_y1252 <- stmetadata2_1252$V4

# 将封闭曲线和点坐标转为矩阵
polygon_coords1252_50 <- cbind(polygon_x1252_50, polygon_y1252_50)
polygon_coords1252_80 <- cbind(polygon_x1252_80, polygon_y1252_80)
polygon_coords1252_95 <- cbind(polygon_x1252_95, polygon_y1252_95)
polygon_coords1252_99 <- cbind(polygon_x1252_99, polygon_y1252_99)
points_coords1252 <- cbind(points_x1252, points_y1252)

# 创建 SpatialPolygons对象（封闭曲线多边形）
polygon_sp1252_50 <- Polygon(polygon_coords1252_50)
polygon_sp1252_80 <- Polygon(polygon_coords1252_80)
polygon_sp1252_95 <- Polygon(polygon_coords1252_95)
polygon_sp1252_99 <- Polygon(polygon_coords1252_99)
polygons_sp1252_50 <- Polygons(list(polygon_sp1252_50), "polygon1")
polygons_sp1252_80 <- Polygons(list(polygon_sp1252_80), "polygon1")
polygons_sp1252_95 <- Polygons(list(polygon_sp1252_95), "polygon1")
polygons_sp1252_99 <- Polygons(list(polygon_sp1252_99), "polygon1")
sp_poly1252_50 <- SpatialPolygons(list(polygons_sp1252_50))
sp_poly1252_80 <- SpatialPolygons(list(polygons_sp1252_80))
sp_poly1252_95 <- SpatialPolygons(list(polygons_sp1252_95))
sp_poly1252_99 <- SpatialPolygons(list(polygons_sp1252_99))
# 使用 polygon() 函数绘制多边形
plot(1, 1, type="n", xlim=c(0, 130), ylim=c(0, 130), xlab="X", ylab="Y")
polygon(polygon_x1252_50, polygon_y1252_50, col="skyblue", border="black")
polygon(polygon_x1252_80, polygon_y1252_80, col="skyblue", border="black")
polygon(polygon_x1252_95, polygon_y1252_95, col="skyblue", border="black")
polygon(polygon_x1252_99, polygon_y1252_99, col="skyblue", border="black")

# 使用 point.in.polygon 函数检查每个点是否在多边形内
result1252_50 <- point.in.polygon(points_x1252, points_y1252, polygon_x1252_50, polygon_y1252_50)
result1252_80 <- point.in.polygon(points_x1252, points_y1252, polygon_x1252_80, polygon_y1252_80)
result1252_95 <- point.in.polygon(points_x1252, points_y1252, polygon_x1252_95, polygon_y1252_95)
result1252_99 <- point.in.polygon(points_x1252, points_y1252, polygon_x1252_99, polygon_y1252_99)

# 输出结果：结果为1表示在多边形内，0表示不在多边形内
stmetadata2_1252$result50 <- result1252_50
stmetadata2_1252$result80 <- result1252_80
stmetadata2_1252$result95 <- result1252_95
stmetadata2_1252$result99 <- result1252_99


stmetadata2_1252$result <- ifelse(stmetadata2_1252$result50==1,"50%",
                                  ifelse(stmetadata2_1252$result80==1,"80%","95%"))
table(stmetadata2_1252$result)
write.csv(stmetadata2_1252,file = "D:/细胞状态2/论文/挖掘/独立验证集/等高线/denggaoxian1252.csv")


ggplot(data = stmetadata2_1252, aes(x = result, y = model1, fill = result)) +
  geom_violin(position = position_dodge(0.8))+
  geom_boxplot(width=0.5,
               #outlier.colour = NA,
               position = position_dodge(0.8))+  #每组的小提琴图之间的距离
  scale_fill_manual(values = c("#e99593", "#99bca0","steelblue"))+
  theme(axis.text.x = element_text(size = 20,angle = 30, hjust = 0.5, vjust = 0.5))

EMT
Angiogenesis
Invasion
Metastasis












#1253给点划分区域
plot_data1253 <- ggplot_build(denggaoxian1253)$data[[1]]  
prob50_spotshigh1253 <- plot_data1253[plot_data1253$group == "2-004", ]
prob80_spotshigh1253 <- plot_data1253[plot_data1253$group == "2-003", ]
prob95_spotshigh1253 <- plot_data1253[plot_data1253$group == "2-002", ]
prob99_spotshigh1253 <- plot_data1253[plot_data1253$group == "2-001", ]


max(prob99_spotshigh1253$x)
min(prob99_spotshigh1253$x)
max(prob99_spotshigh1253$y)
min(prob99_spotshigh1253$y)


# 封闭曲线（多边形）坐标
polygon_x1253_50 <- prob50_spotshigh1253$x  
polygon_y1253_50 <- prob50_spotshigh1253$y  
polygon_x1253_80 <- prob80_spotshigh1253$x
polygon_y1253_80 <- prob80_spotshigh1253$y
polygon_x1253_95 <- prob95_spotshigh1253$x
polygon_y1253_95 <- prob95_spotshigh1253$y
polygon_x1253_99 <- prob99_spotshigh1253$x
polygon_y1253_99 <- prob99_spotshigh1253$y

# 待检测点的坐标
points_x1253 <- stmetadata2_1253$V3
points_y1253 <- stmetadata2_1253$V4

# 将封闭曲线和点坐标转为矩阵
polygon_coords1253_50 <- cbind(polygon_x1253_50, polygon_y1253_50)
polygon_coords1253_80 <- cbind(polygon_x1253_80, polygon_y1253_80)
polygon_coords1253_95 <- cbind(polygon_x1253_95, polygon_y1253_95)
polygon_coords1253_99 <- cbind(polygon_x1253_99, polygon_y1253_99)
points_coords1253 <- cbind(points_x1253, points_y1253)

# 创建 SpatialPolygons对象（封闭曲线多边形）
polygon_sp1253_50 <- Polygon(polygon_coords1253_50)
polygon_sp1253_80 <- Polygon(polygon_coords1253_80)
polygon_sp1253_95 <- Polygon(polygon_coords1253_95)
polygon_sp1253_99 <- Polygon(polygon_coords1253_99)
polygons_sp1253_50 <- Polygons(list(polygon_sp1253_50), "polygon1")
polygons_sp1253_80 <- Polygons(list(polygon_sp1253_80), "polygon1")
polygons_sp1253_95 <- Polygons(list(polygon_sp1253_95), "polygon1")
polygons_sp1253_99 <- Polygons(list(polygon_sp1253_99), "polygon1")
sp_poly1253_50 <- SpatialPolygons(list(polygons_sp1253_50))
sp_poly1253_80 <- SpatialPolygons(list(polygons_sp1253_80))
sp_poly1253_95 <- SpatialPolygons(list(polygons_sp1253_95))
sp_poly1253_99 <- SpatialPolygons(list(polygons_sp1253_99))
# 使用 polygon() 函数绘制多边形
plot(1, 1, type="n", xlim=c(0, 130), ylim=c(0, 130), xlab="X", ylab="Y")
polygon(polygon_x1253_50, polygon_y1253_50, col="skyblue", border="black")
polygon(polygon_x1253_80, polygon_y1253_80, col="skyblue", border="black")
polygon(polygon_x1253_95, polygon_y1253_95, col="skyblue", border="black")
polygon(polygon_x1253_99, polygon_y1253_99, col="skyblue", border="black")

# 使用 point.in.polygon 函数检查每个点是否在多边形内
result1253_50 <- point.in.polygon(points_x1253, points_y1253, polygon_x1253_50, polygon_y1253_50)
result1253_80 <- point.in.polygon(points_x1253, points_y1253, polygon_x1253_80, polygon_y1253_80)
result1253_95 <- point.in.polygon(points_x1253, points_y1253, polygon_x1253_95, polygon_y1253_95)
result1253_99 <- point.in.polygon(points_x1253, points_y1253, polygon_x1253_99, polygon_y1253_99)

# 输出结果：结果为1表示在多边形内，0表示不在多边形内
stmetadata2_1253$result50 <- result1253_50
stmetadata2_1253$result80 <- result1253_80
stmetadata2_1253$result95 <- result1253_95
stmetadata2_1253$result99 <- result1253_99


stmetadata2_1253$result <- ifelse(stmetadata2_1253$result50==1,"50%",
                                  ifelse(stmetadata2_1253$result80==1,"80%","95%"))
table(stmetadata2_1253$result)
write.csv(stmetadata2_1253,file = "D:/细胞状态2/论文/挖掘/独立验证集/等高线/denggaoxian1253.csv")


ggplot(data = stmetadata2_1253, aes(x = result, y = model1, fill = result)) +
  geom_violin(position = position_dodge(0.8))+
  geom_boxplot(width=0.5,
               #outlier.colour = NA,
               position = position_dodge(0.8))+  #每组的小提琴图之间的距离
  scale_fill_manual(values = c("#e99593", "#99bca0","steelblue"))+
  theme(axis.text.x = element_text(size = 20,angle = 30, hjust = 0.5, vjust = 0.5))

EMT
Angiogenesis
Invasion
Metastasis










#1254给点划分区域
plot_data1254 <- ggplot_build(denggaoxian1254)$data[[1]]  
prob50_spotshigh1254 <- plot_data1254[plot_data1254$group == "2-004", ]
prob80_spotshigh1254 <- plot_data1254[plot_data1254$group == "2-003", ]
prob95_spotshigh1254 <- plot_data1254[plot_data1254$group == "2-002", ]
prob99_spotshigh1254 <- plot_data1254[plot_data1254$group == "2-001", ]


max(prob99_spotshigh1254$x)
min(prob99_spotshigh1254$x)
max(prob99_spotshigh1254$y)
min(prob99_spotshigh1254$y)


# 封闭曲线（多边形）坐标
polygon_x1254_50 <- prob50_spotshigh1254$x  
polygon_y1254_50 <- prob50_spotshigh1254$y  
polygon_x1254_80 <- prob80_spotshigh1254$x
polygon_y1254_80 <- prob80_spotshigh1254$y
polygon_x1254_95 <- prob95_spotshigh1254$x
polygon_y1254_95 <- prob95_spotshigh1254$y
polygon_x1254_99 <- prob99_spotshigh1254$x
polygon_y1254_99 <- prob99_spotshigh1254$y

# 待检测点的坐标
points_x1254 <- stmetadata2_1254$V3
points_y1254 <- stmetadata2_1254$V4

# 将封闭曲线和点坐标转为矩阵
polygon_coords1254_50 <- cbind(polygon_x1254_50, polygon_y1254_50)
polygon_coords1254_80 <- cbind(polygon_x1254_80, polygon_y1254_80)
polygon_coords1254_95 <- cbind(polygon_x1254_95, polygon_y1254_95)
polygon_coords1254_99 <- cbind(polygon_x1254_99, polygon_y1254_99)
points_coords1254 <- cbind(points_x1254, points_y1254)

# 创建 SpatialPolygons对象（封闭曲线多边形）
polygon_sp1254_50 <- Polygon(polygon_coords1254_50)
polygon_sp1254_80 <- Polygon(polygon_coords1254_80)
polygon_sp1254_95 <- Polygon(polygon_coords1254_95)
polygon_sp1254_99 <- Polygon(polygon_coords1254_99)
polygons_sp1254_50 <- Polygons(list(polygon_sp1254_50), "polygon1")
polygons_sp1254_80 <- Polygons(list(polygon_sp1254_80), "polygon1")
polygons_sp1254_95 <- Polygons(list(polygon_sp1254_95), "polygon1")
polygons_sp1254_99 <- Polygons(list(polygon_sp1254_99), "polygon1")
sp_poly1254_50 <- SpatialPolygons(list(polygons_sp1254_50))
sp_poly1254_80 <- SpatialPolygons(list(polygons_sp1254_80))
sp_poly1254_95 <- SpatialPolygons(list(polygons_sp1254_95))
sp_poly1254_99 <- SpatialPolygons(list(polygons_sp1254_99))
# 使用 polygon() 函数绘制多边形
plot(1, 1, type="n", xlim=c(0, 130), ylim=c(0, 130), xlab="X", ylab="Y")
polygon(polygon_x1254_50, polygon_y1254_50, col="skyblue", border="black")
polygon(polygon_x1254_80, polygon_y1254_80, col="skyblue", border="black")
polygon(polygon_x1254_95, polygon_y1254_95, col="skyblue", border="black")
polygon(polygon_x1254_99, polygon_y1254_99, col="skyblue", border="black")

# 使用 point.in.polygon 函数检查每个点是否在多边形内
result1254_50 <- point.in.polygon(points_x1254, points_y1254, polygon_x1254_50, polygon_y1254_50)
result1254_80 <- point.in.polygon(points_x1254, points_y1254, polygon_x1254_80, polygon_y1254_80)
result1254_95 <- point.in.polygon(points_x1254, points_y1254, polygon_x1254_95, polygon_y1254_95)
result1254_99 <- point.in.polygon(points_x1254, points_y1254, polygon_x1254_99, polygon_y1254_99)

# 输出结果：结果为1表示在多边形内，0表示不在多边形内
stmetadata2_1254$result50 <- result1254_50
stmetadata2_1254$result80 <- result1254_80
stmetadata2_1254$result95 <- result1254_95
stmetadata2_1254$result99 <- result1254_99


stmetadata2_1254$result <- ifelse(stmetadata2_1254$result50==1,"50%",
                                  ifelse(stmetadata2_1254$result80==1,"80%","95%"))
table(stmetadata2_1254$result)
write.csv(stmetadata2_1254,file = "D:/细胞状态2/论文/挖掘/独立验证集/等高线/denggaoxian1254.csv")


ggplot(data = stmetadata2_1254, aes(x = result, y = model1, fill = result)) +
  geom_violin(position = position_dodge(0.8))+
  geom_boxplot(width=0.5,
               #outlier.colour = NA,
               position = position_dodge(0.8))+  #每组的小提琴图之间的距离
  scale_fill_manual(values = c("#e99593", "#99bca0","steelblue"))+
  theme(axis.text.x = element_text(size = 20,angle = 30, hjust = 0.5, vjust = 0.5))

EMT
Angiogenesis
Invasion
Metastasis













#1675给点划分区域
plot_data1675 <- ggplot_build(denggaoxian1675)$data[[1]]  
prob50_spotshigh1675 <- plot_data1675[plot_data1675$group == "2-004", ]
prob80_spotshigh1675 <- plot_data1675[plot_data1675$group == "2-003", ]
prob95_spotshigh1675 <- plot_data1675[plot_data1675$group == "2-002", ]
prob99_spotshigh1675 <- plot_data1675[plot_data1675$group == "2-001", ]


max(prob99_spotshigh1675$x)
min(prob99_spotshigh1675$x)
max(prob99_spotshigh1675$y)
min(prob99_spotshigh1675$y)


# 封闭曲线（多边形）坐标
polygon_x1675_50 <- prob50_spotshigh1675$x  
polygon_y1675_50 <- prob50_spotshigh1675$y  
polygon_x1675_80 <- prob80_spotshigh1675$x
polygon_y1675_80 <- prob80_spotshigh1675$y
polygon_x1675_95 <- prob95_spotshigh1675$x
polygon_y1675_95 <- prob95_spotshigh1675$y
polygon_x1675_99 <- prob99_spotshigh1675$x
polygon_y1675_99 <- prob99_spotshigh1675$y

# 待检测点的坐标
points_x1675 <- stmetadata2_1675$V3
points_y1675 <- stmetadata2_1675$V4

# 将封闭曲线和点坐标转为矩阵
polygon_coords1675_50 <- cbind(polygon_x1675_50, polygon_y1675_50)
polygon_coords1675_80 <- cbind(polygon_x1675_80, polygon_y1675_80)
polygon_coords1675_95 <- cbind(polygon_x1675_95, polygon_y1675_95)
polygon_coords1675_99 <- cbind(polygon_x1675_99, polygon_y1675_99)
points_coords1675 <- cbind(points_x1675, points_y1675)

# 创建 SpatialPolygons对象（封闭曲线多边形）
polygon_sp1675_50 <- Polygon(polygon_coords1675_50)
polygon_sp1675_80 <- Polygon(polygon_coords1675_80)
polygon_sp1675_95 <- Polygon(polygon_coords1675_95)
polygon_sp1675_99 <- Polygon(polygon_coords1675_99)
polygons_sp1675_50 <- Polygons(list(polygon_sp1675_50), "polygon1")
polygons_sp1675_80 <- Polygons(list(polygon_sp1675_80), "polygon1")
polygons_sp1675_95 <- Polygons(list(polygon_sp1675_95), "polygon1")
polygons_sp1675_99 <- Polygons(list(polygon_sp1675_99), "polygon1")
sp_poly1675_50 <- SpatialPolygons(list(polygons_sp1675_50))
sp_poly1675_80 <- SpatialPolygons(list(polygons_sp1675_80))
sp_poly1675_95 <- SpatialPolygons(list(polygons_sp1675_95))
sp_poly1675_99 <- SpatialPolygons(list(polygons_sp1675_99))
# 使用 polygon() 函数绘制多边形
plot(1, 1, type="n", xlim=c(0, 130), ylim=c(0, 130), xlab="X", ylab="Y")
polygon(polygon_x1675_50, polygon_y1675_50, col="skyblue", border="black")
polygon(polygon_x1675_80, polygon_y1675_80, col="skyblue", border="black")
polygon(polygon_x1675_95, polygon_y1675_95, col="skyblue", border="black")
polygon(polygon_x1675_99, polygon_y1675_99, col="skyblue", border="black")

# 使用 point.in.polygon 函数检查每个点是否在多边形内
result1675_50 <- point.in.polygon(points_x1675, points_y1675, polygon_x1675_50, polygon_y1675_50)
result1675_80 <- point.in.polygon(points_x1675, points_y1675, polygon_x1675_80, polygon_y1675_80)
result1675_95 <- point.in.polygon(points_x1675, points_y1675, polygon_x1675_95, polygon_y1675_95)
result1675_99 <- point.in.polygon(points_x1675, points_y1675, polygon_x1675_99, polygon_y1675_99)

# 输出结果：结果为1表示在多边形内，0表示不在多边形内
stmetadata2_1675$result50 <- result1675_50
stmetadata2_1675$result80 <- result1675_80
stmetadata2_1675$result95 <- result1675_95
stmetadata2_1675$result99 <- result1675_99


stmetadata2_1675$result <- ifelse(stmetadata2_1675$result50==1,"50%",
                                  ifelse(stmetadata2_1675$result80==1,"80%","95%"))
table(stmetadata2_1675$result)
write.csv(stmetadata2_1675,file = "D:/细胞状态2/论文/挖掘/独立验证集/等高线/denggaoxian1675.csv")


ggplot(data = stmetadata2_1675, aes(x = result, y = model1, fill = result)) +
  geom_violin(position = position_dodge(0.8))+
  geom_boxplot(width=0.5,
               #outlier.colour = NA,
               position = position_dodge(0.8))+  #每组的小提琴图之间的距离
  scale_fill_manual(values = c("#e99593", "#99bca0","steelblue"))+
  theme(axis.text.x = element_text(size = 20,angle = 30, hjust = 0.5, vjust = 0.5))

EMT
Angiogenesis
Invasion
Metastasis












#奇特的基因表达
# load("D:/细胞状态2/论文/空间/stad/stadexp.RData")
# load("D:/细胞状态2/论文/空间/stad/stadseurat_data.RData")
# load("D:/细胞状态2/论文/空间/stad/stmetadata.RData")
# load("D:/细胞状态2/论文/空间/stad/ES.RData")

load("/wangzhenzenlab/fsh/xff/stdata/duli/stadexp.RData")
load("/wangzhenzenlab/fsh/xff/stdata/duli/stadseurat_data.RData")
load("/wangzhenzenlab/fsh/xff/stdata/duli/stmetadata.RData")
load("/wangzhenzenlab/fsh/xff/stdata/duli/ES.RData")

# stmetadata2_0331 <- read.csv("D:/细胞状态2/论文/挖掘/独立验证集/等高线/denggaoxian0331.csv")
# stmetadata2_0731 <- read.csv("D:/细胞状态2/论文/挖掘/独立验证集/等高线/denggaoxian0731.csv")
# stmetadata2_0733 <- read.csv("D:/细胞状态2/论文/挖掘/独立验证集/等高线/denggaoxian0733.csv")
# stmetadata2_0734 <- read.csv("D:/细胞状态2/论文/挖掘/独立验证集/等高线/denggaoxian0734.csv")
# stmetadata2_1251 <- read.csv("D:/细胞状态2/论文/挖掘/独立验证集/等高线/denggaoxian1251.csv")
# stmetadata2_1252 <- read.csv("D:/细胞状态2/论文/挖掘/独立验证集/等高线/denggaoxian1252.csv")
# stmetadata2_1253 <- read.csv("D:/细胞状态2/论文/挖掘/独立验证集/等高线/denggaoxian1253.csv")
# stmetadata2_1254 <- read.csv("D:/细胞状态2/论文/挖掘/独立验证集/等高线/denggaoxian1254.csv")
# stmetadata2_1675 <- read.csv("D:/细胞状态2/论文/挖掘/独立验证集/等高线/denggaoxian1675.csv")

stmetadata2_0331 <- read.csv("/wangzhenzenlab/fsh/xff/stdata/duli/denggaoxian0331.csv")
stmetadata2_0731 <- read.csv("/wangzhenzenlab/fsh/xff/stdata/duli/denggaoxian0731.csv")
stmetadata2_0733 <- read.csv("/wangzhenzenlab/fsh/xff/stdata/duli/denggaoxian0733.csv")
stmetadata2_0734 <- read.csv("/wangzhenzenlab/fsh/xff/stdata/duli/denggaoxian0734.csv")
stmetadata2_1251 <- read.csv("/wangzhenzenlab/fsh/xff/stdata/duli/denggaoxian1251.csv")
stmetadata2_1252 <- read.csv("/wangzhenzenlab/fsh/xff/stdata/duli/denggaoxian1252.csv")
stmetadata2_1253 <- read.csv("/wangzhenzenlab/fsh/xff/stdata/duli/denggaoxian1253.csv")
stmetadata2_1254 <- read.csv("/wangzhenzenlab/fsh/xff/stdata/duli/denggaoxian1254.csv")
stmetadata2_1675 <- read.csv("/wangzhenzenlab/fsh/xff/stdata/duli/denggaoxian1675.csv")

expression_matrix0331[1:3,1:3]

expression_matrix0331_2 <- as.data.frame(t(expression_matrix0331))
expression_matrix0331_2$cell <- rownames(expression_matrix0331_2)
expression_matrix0331_3 <- merge(stmetadata2_0331,expression_matrix0331_2,by="cell")
stadexpression_matrix  <- expression_matrix0331
stadexpression_matrix3 <- expression_matrix0331_3

expression_matrix0731_2 <- as.data.frame(t(expression_matrix0731))
expression_matrix0731_2$cell <- rownames(expression_matrix0731_2)
expression_matrix0731_3 <- merge(stmetadata2_0731, expression_matrix0731_2, by = "cell")
stadexpression_matrix  <- expression_matrix0731
stadexpression_matrix3 <- expression_matrix0731_3

expression_matrix0733_2 <- as.data.frame(t(expression_matrix0733))
expression_matrix0733_2$cell <- rownames(expression_matrix0733_2)
expression_matrix0733_3 <- merge(stmetadata2_0733, expression_matrix0733_2, by = "cell")
stadexpression_matrix  <- expression_matrix0733
stadexpression_matrix3 <- expression_matrix0733_3

expression_matrix0734_2 <- as.data.frame(t(expression_matrix0734))
expression_matrix0734_2$cell <- rownames(expression_matrix0734_2)
expression_matrix0734_3 <- merge(stmetadata2_0734, expression_matrix0734_2, by = "cell")
stadexpression_matrix  <- expression_matrix0734
stadexpression_matrix3 <- expression_matrix0734_3

expression_matrix1251_2 <- as.data.frame(t(expression_matrix1251))
expression_matrix1251_2$cell <- rownames(expression_matrix1251_2)
expression_matrix1251_3 <- merge(stmetadata2_1251, expression_matrix1251_2, by = "cell")
stadexpression_matrix  <- expression_matrix1251
stadexpression_matrix3 <- expression_matrix1251_3

expression_matrix1252_2 <- as.data.frame(t(expression_matrix1252))
expression_matrix1252_2$cell <- rownames(expression_matrix1252_2)
expression_matrix1252_3 <- merge(stmetadata2_1252, expression_matrix1252_2, by = "cell")
stadexpression_matrix  <- expression_matrix1252
stadexpression_matrix3 <- expression_matrix1252_3

expression_matrix1253_2 <- as.data.frame(t(expression_matrix1253))
expression_matrix1253_2$cell <- rownames(expression_matrix1253_2)
expression_matrix1253_3 <- merge(stmetadata2_1253, expression_matrix1253_2, by = "cell")
stadexpression_matrix  <- expression_matrix1253
stadexpression_matrix3 <- expression_matrix1253_3

expression_matrix1254_2 <- as.data.frame(t(expression_matrix1254))
expression_matrix1254_2$cell <- rownames(expression_matrix1254_2)
expression_matrix1254_3 <- merge(stmetadata2_1254, expression_matrix1254_2, by = "cell")
stadexpression_matrix  <- expression_matrix1254
stadexpression_matrix3 <- expression_matrix1254_3

expression_matrix1675_2 <- as.data.frame(t(expression_matrix1675))
expression_matrix1675_2$cell <- rownames(expression_matrix1675_2)
expression_matrix1675_3 <- merge(stmetadata2_1675, expression_matrix1675_2, by = "cell")
stadexpression_matrix  <- expression_matrix1675
stadexpression_matrix3 <- expression_matrix1675_3

midmore <- data.frame(row.names = "test")
midmore$highmean <- 0
midmore$midmean <- 0
midmore$lowmean <- 0
midmore$FC1 <- 0
midmore$FC2 <- 0

for (i in rownames(stadexpression_matrix)) {
  highmean<- mean(stadexpression_matrix3[stadexpression_matrix3$result=="50%",i])
  midmean <- mean(stadexpression_matrix3[stadexpression_matrix3$result=="80%",i])
  lowmean <- mean(stadexpression_matrix3[stadexpression_matrix3$result=="95%",i])
  if(midmean>lowmean&midmean>highmean){
    genemore <- data.frame(row.names = i)
    genemore$highmean <- highmean
    genemore$midmean <- midmean
    genemore$lowmean <- lowmean
    genemore$FC1 <- midmean/highmean
    genemore$FC2 <- midmean/lowmean
    midmore <- rbind(midmore,genemore)
  }
}

dim(midmore)
head(midmore)



midhigh <- midmore[midmore$FC2>2 & midmore$FC1>2& midmore$midmean>0.01,]
dim(midhigh)
head(midhigh)

stadexpression_matrix3[1:3,1:40]
rownames(stadexpression_matrix3) <- stadexpression_matrix3$cell
midhighexp<- stadexpression_matrix3[,rownames(midhigh)]
head(midhighexp)

write.csv(midhighexp, file = "/wangzhenzenlab/fsh/xff/stdata/duli/midhighexp0331.csv")
write.csv(midhighexp, file = "/wangzhenzenlab/fsh/xff/stdata/duli/midhighexp0731.csv")
write.csv(midhighexp, file = "/wangzhenzenlab/fsh/xff/stdata/duli/midhighexp0733.csv")
write.csv(midhighexp, file = "/wangzhenzenlab/fsh/xff/stdata/duli/midhighexp0734.csv")
write.csv(midhighexp, file = "/wangzhenzenlab/fsh/xff/stdata/duli/midhighexp1251.csv")
write.csv(midhighexp, file = "/wangzhenzenlab/fsh/xff/stdata/duli/midhighexp1252.csv")
write.csv(midhighexp, file = "/wangzhenzenlab/fsh/xff/stdata/duli/midhighexp1253.csv")
write.csv(midhighexp, file = "/wangzhenzenlab/fsh/xff/stdata/duli/midhighexp1254.csv")
write.csv(midhighexp, file = "/wangzhenzenlab/fsh/xff/stdata/duli/midhighexp1675.csv")


library(corrplot)
singlemt <- midhighexp
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
my_colors <- colorRampPalette(c("steelblue", "white", "#c11e21"))(100)
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)

cor_pearson2 <- as.data.frame(cor_pearson)
cor_pearson2$gene1 <- rownames(cor_pearson2)

library(tidyr)
cor_pearsonlong <- cor_pearson2 %>%
  pivot_longer(
    cols = -gene1, # 指定要转换的列（以 "score" 开头的列）
    names_to = "gene2",     # 新生成的列，表示原列名
    values_to = "corvalue"    # 新生成的列，表示原列的值
  )
head(cor_pearsonlong)

cor_pearsonlong2 <- cor_pearsonlong[cor_pearsonlong$corvalue>0.6&cor_pearsonlong$corvalue!=1,]
unique(cor_pearsonlong2$gene1)





#先前基因
pregene <- c("CXCL8","CXCL6","MMP3",
             "FCRLA","CR2","FDCSP","CXCL13","MS4A1","TCL1A","VPREB3",
             "PRSS21","DPEP3","COL9A3","C20orf204")


expression_matrix0331_2[1:3,1:40]
write.csv(expression_matrix0331_2[, pregene], file = "D:/细胞状态2/论文/挖掘/独立验证集/等高线/pregene0331.csv")
write.csv(expression_matrix0731_2[, pregene], file = "D:/细胞状态2/论文/挖掘/独立验证集/等高线/pregene0731.csv")
write.csv(expression_matrix0733_2[, pregene], file = "D:/细胞状态2/论文/挖掘/独立验证集/等高线/pregene0733.csv")
write.csv(expression_matrix0734_2[, pregene], file = "D:/细胞状态2/论文/挖掘/独立验证集/等高线/pregene0734.csv")
write.csv(expression_matrix1251_2[, pregene], file = "D:/细胞状态2/论文/挖掘/独立验证集/等高线/pregene1251.csv")
write.csv(expression_matrix1252_2[, pregene], file = "D:/细胞状态2/论文/挖掘/独立验证集/等高线/pregene1252.csv")
write.csv(expression_matrix1253_2[, pregene], file = "D:/细胞状态2/论文/挖掘/独立验证集/等高线/pregene1253.csv")
write.csv(expression_matrix1254_2[, pregene], file = "D:/细胞状态2/论文/挖掘/独立验证集/等高线/pregene1254.csv")
write.csv(expression_matrix1675_2[, pregene], file = "D:/细胞状态2/论文/挖掘/独立验证集/等高线/pregene1675.csv")





