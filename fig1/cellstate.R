load("/boot2/xff/cellstate/cellstatemergedata.RData")
load("/boot2/xff/gdh_data2/annodataALL.RData")
load("/wangzhenzenlab/xff/cellstate/xff/cellstate/cellstatemergedata.RData")

head(cellstatemergedata)
#去掉NHL
cellstatemergedata <- cellstatemergedata[!cellstatemergedata$cancertype=="NHL",]
table(cellstatemergedata$"Celltype (major-lineage)")

TGSVA <- cellstatemergedata[cellstatemergedata$"Celltype (minor-lineage)"%in%
                              c("CD4T", "CD4Teff", "CD4Tn", "Th1", "Th2", "Th17",
                                "CD8Teff", "Treg", "Tfh", "CD4Tconv", "Tprolif", 
                                "CD8T", "CD8Tex", "CD8Tcm",  "NK","CD8Tn", "CD8Tem"),]
head(TGSVA)
dim(TGSVA)
table(TGSVA$`Celltype (minor-lineage)`)
save(TGSVA,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate/TGSVA.RData")

#去掉特异性细胞
table(cellstatemergedata$`Celltype (major-lineage)`)
celltype49 <- read.csv("/boot2/xff/55celltype.csv")
cellstatemergedata <- merge(celltype49,cellstatemergedata,by.x = "Celltype",by.y = "Celltype (major-lineage)")
cellstatemergedata$Celltype <- as.character(cellstatemergedata$Celltype)
#查看细胞类型
table(cellstatemergedata$Celltype)
#查看癌症类型
table(cellstatemergedata$cancertype)
#查看来源组织
table(cellstatemergedata$Source)
names(table(cellstatemergedata$Source))
class(cellstatemergedata$Source)

#分析来源组织是正常还是癌症
test1 <- cellstatemergedata[cellstatemergedata$Source=="nodule",]
table(test1$Stage)
#Glioma Adult Pediatric  GSE131928
#LIHC  PBMC  GSE98638 Ascites Blood Lymphnode GSE140228
#SKCM blood PBMC GSE123139 nodule center cortex GSE139249
#TISCH2_185pancancer_scRNA-seq\SKCM_GSE159251\SKCM_GSE159251_CellMetainfo_table.tsv
#SKCM_GSE179373 SKCM_GSE166181_aPD1
cellstatemergedata$Source[cellstatemergedata$Source=="BMMC"&cellstatemergedata$Gender=="None",]

test2<- cellstatemergedata[cellstatemergedata$cancertype=="UVM",]
test3 <- cellstatemergedata[cellstatemergedata$Celltype=="Oligodendrocyte",]
#修改来源为normal和tumor
cellstatemergedata$Source <- ifelse(cellstatemergedata$Gender == "None" & cellstatemergedata$Source == "BMMC", "Normal", cellstatemergedata$Source)
cellstatemergedata$Source <- ifelse(cellstatemergedata$Stage %in% c("Metastatic", "Primary") & cellstatemergedata$Source == "PBMC", "Tumor", cellstatemergedata$Source)

#分类
NORMAL <- cellstatemergedata[cellstatemergedata$Source %in% c("Normal", "NAT", "Adjacent normal",
                                                              "PBMC","healthy","BM","Cell culture",
                                                              "Skin","LN","Tonsil"), ]
table(NORMAL$Source)
table(NORMAL$`Celltype (major-lineage)`)
table(NORMAL$cancertype)


TUMOR <- cellstatemergedata[cellstatemergedata$Source %in% c("AML", "Ascites", "B-ALL","BMMC",
                                                             "center","cortex","Human tumor",
                                                             "Leptomeningeal metastasis","mBrain",
                                                             "Metastasis","Metastatic","mLiver",
                                                             "mFallopian_tube","nodule","Primary",
                                                             "Relapsed","T-ALL","TIL","tLN","tumor",
                                                             "Tumor"), ]
table(TUMOR$Celltype)
table(TUMOR$cancertype)
test3 <- TUMOR[TUMOR$Celltype=="Oligodendrocyte",]
head(TUMOR)


#不同癌症与细胞状态
cancertype_cellstate <- TUMOR[,c(3:16,28)]
head(cancertype_cellstate)
cancertype_cellstatemedian <- aggregate(. ~ cancertype, data = cancertype_cellstate, mean)
rownames(cancertype_cellstatemedian) <- cancertype_cellstatemedian$cancertype
cancertype_cellstatemedian$cancertype <- NULL
pdf("/boot2/xff/cellstate2/cancertype_cellstatemean.pdf")
pheatmap(cancertype_cellstatemedian, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = T, 
         #cluster_rows表示仅对行聚类，cluster_cols表示仅对列聚类，值为TRUE或FALSE
         cutree_rows = 5, cutree_cols = NA, #若进行了行/列聚类，根据行/列聚类数量分隔热图行,cutree_rows=num分割行，cutree_cols=num分割列
         #treeheight_row = 30, treeheight_col = 30, #若行、列聚类树高度调整
         border = NA, #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         cellwidth = NA, cellheight = NA,  #表示单个单元格的宽度\高度，默认为 “NA”
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
         
)

dev.off()

NORMAL_cellstate <- NORMAL[,c(3:16,28)]
head(NORMAL_cellstate)
NORMAL_cellstatemedian <- aggregate(. ~ cancertype, data = NORMAL_cellstate, mean)
rownames(NORMAL_cellstatemedian) <- NORMAL_cellstatemedian$cancertype
NORMAL_cellstatemedian$cancertype <- NULL
pdf("/boot2/xff/cellstate2/NORMAL_cellstatemean.pdf")
pheatmap(NORMAL_cellstatemedian, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = T, 
         #cluster_rows表示仅对行聚类，cluster_cols表示仅对列聚类，值为TRUE或FALSE
         cutree_rows = 5, cutree_cols = NA, #若进行了行/列聚类，根据行/列聚类数量分隔热图行,cutree_rows=num分割行，cutree_cols=num分割列
         #treeheight_row = 30, treeheight_col = 30, #若行、列聚类树高度调整
         border = NA, #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         cellwidth = NA, cellheight = NA,  #表示单个单元格的宽度\高度，默认为 “NA”
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
         
)

dev.off()

CANCER <- "ALL"
foldchange <- function(CANCER){
  cancerTUMOR <- cancertype_cellstatemedian[rownames(cancertype_cellstatemedian)==CANCER,]
  cancerNORMAL <- NORMAL_cellstatemedian[rownames(NORMAL_cellstatemedian)==CANCER,]
  FC <- log2((cancerTUMOR+1)/(cancerNORMAL+1))
  return(FC)
}
ALLFC <- foldchange("ALL")
allFC <- data.frame()
for (i in rownames(NORMAL_cellstatemedian)) {
  print(i)
  if(i!="PBMC"){
    singleFC <- foldchange(i)
    allFC <- rbind(allFC,singleFC)
  }
}
allFC
allFC <- allFC[!rownames(allFC) %in% c("ALL", "AML", "BLCA", "Glioma","HNSC",
                                       "KIPAN","MCC","MF","MM","NPC","UVM"),]
custom_colors <- colorRampPalette(c("steelblue","white","#e99593"))
num_breaks <- 100
breaks <- seq(-0.1, 0.1, length.out = num_breaks + 1)
pdf("/boot2/xff/cellstate2/cancertypeheatmap/cancertypeFC.pdf")
pheatmap(allFC, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = T, 
         #cluster_rows表示仅对行聚类，cluster_cols表示仅对列聚类，值为TRUE或FALSE
         cutree_rows = 4, cutree_cols = NA, #若进行了行/列聚类，根据行/列聚类数量分隔热图行,cutree_rows=num分割行，cutree_cols=num分割列
         #treeheight_row = 30, treeheight_col = 30, #若行、列聚类树高度调整
         border = NA, #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         cellwidth = NA, cellheight = NA,  #表示单个单元格的宽度\高度，默认为 “NA”
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
         
)

dev.off()


#不同细胞与细胞状态
head(TUMOR)
#加上细胞大类
CELLTYPE <- read.csv("/boot2/xff/celltype.csv")
#cellstatemergedata$CELLTYPE.y <- NULL
TUMOR <- merge(TUMOR,CELLTYPE,by.x = "Celltype",by.y = "celltype")
TUMOR$celltype <- paste0(TUMOR$Celltype,"_", substr(TUMOR$CELLTYPE, 1, 1))
head(TUMOR)
names(table(TUMOR$celltype))

NORMAL <- merge(NORMAL,CELLTYPE,by.x = "Celltype",by.y = "celltype")
NORMAL$celltype <- paste0(NORMAL$Celltype,"_", substr(NORMAL$CELLTYPE, 1, 1))
head(NORMAL)
names(table(NORMAL$celltype))

celltype_cellstate <- TUMOR[,c(3:16,32)]
head(celltype_cellstate)
celltype_cellstatemedian <- aggregate(. ~ celltype, data = celltype_cellstate, mean)
rownames(celltype_cellstatemedian) <- celltype_cellstatemedian$celltype
celltype_cellstatemedian$celltype <- NULL
pdf("/boot2/xff/cellstate2/celltype_cellstatemedian.pdf")
pheatmap(celltype_cellstatemedian, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = T, 
         #cluster_rows表示仅对行聚类，cluster_cols表示仅对列聚类，值为TRUE或FALSE
         cutree_rows = 5, cutree_cols = NA, #若进行了行/列聚类，根据行/列聚类数量分隔热图行,cutree_rows=num分割行，cutree_cols=num分割列
         #treeheight_row = 30, treeheight_col = 30, #若行、列聚类树高度调整
         border = NA, #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         cellwidth = NA, cellheight = NA,  #表示单个单元格的宽度\高度，默认为 “NA”
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
         
)

dev.off()

NORMALcelltype_cellstate <- NORMAL[,c(3:16,32)]
head(NORMALcelltype_cellstate)
NORMALcelltype_cellstatemedian <- aggregate(. ~ celltype, data = NORMALcelltype_cellstate, mean)
rownames(NORMALcelltype_cellstatemedian) <- NORMALcelltype_cellstatemedian$celltype
NORMALcelltype_cellstatemedian$celltype <- NULL
pdf("/boot2/xff/cellstate2/NORMALcelltype_cellstatemedian.pdf")
pheatmap(NORMALcelltype_cellstatemedian, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = T, 
         #cluster_rows表示仅对行聚类，cluster_cols表示仅对列聚类，值为TRUE或FALSE
         cutree_rows = 5, cutree_cols = NA, #若进行了行/列聚类，根据行/列聚类数量分隔热图行,cutree_rows=num分割行，cutree_cols=num分割列
         #treeheight_row = 30, treeheight_col = 30, #若行、列聚类树高度调整
         border = NA, #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         cellwidth = NA, cellheight = NA,  #表示单个单元格的宽度\高度，默认为 “NA”
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
         
)

dev.off()

CELL <- "Acinar_O"
foldchangecell <- function(CELL){
  cellTUMOR <- celltype_cellstatemedian[rownames(celltype_cellstatemedian)==CELL,]
  cellNORMAL <- NORMALcelltype_cellstatemedian[rownames(NORMALcelltype_cellstatemedian)==CELL,]
  FC <- log2((cellTUMOR+1)/(cellNORMAL+1))
  return(FC)
}


Acinar_OFC <- foldchangecell("Acinar_O")
allcellFC <- data.frame()
for (i in rownames(NORMALcelltype_cellstatemedian)) {
  print(i)
  if(i!="PBMC"){
    singleFC <- foldchangecell(i)
    allcellFC <- rbind(allcellFC,singleFC)
  }
}
allcellFC

pdf("/boot2/xff/cellstate2/celltypeFC.pdf")
pheatmap(allcellFC, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = T, 
         #cluster_rows表示仅对行聚类，cluster_cols表示仅对列聚类，值为TRUE或FALSE
         cutree_rows = 4, cutree_cols = NA, #若进行了行/列聚类，根据行/列聚类数量分隔热图行,cutree_rows=num分割行，cutree_cols=num分割列
         #treeheight_row = 30, treeheight_col = 30, #若行、列聚类树高度调整
         border = NA, #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         cellwidth = NA, cellheight = NA,  #表示单个单元格的宽度\高度，默认为 “NA”
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
         
)

dev.off()

#去掉other
TUMOR2 <- TUMOR[TUMOR$CELLTYPE!="O",]
NORMAL2 <- NORMAL[NORMAL$CELLTYPE!="O",]

celltype_cellstate2 <- TUMOR2[,c(3:16,32)]
head(celltype_cellstate2)
celltype_cellstatemedian2 <- aggregate(. ~ celltype, data = celltype_cellstate2, mean)
rownames(celltype_cellstatemedian2) <- celltype_cellstatemedian2$celltype
celltype_cellstatemedian2$celltype <- NULL
pdf("/boot2/xff/cellstate2/celltype_cellstatemedian2.pdf")
pheatmap(celltype_cellstatemedian2, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = T, 
         #cluster_rows表示仅对行聚类，cluster_cols表示仅对列聚类，值为TRUE或FALSE
         cutree_rows = 5, cutree_cols = NA, #若进行了行/列聚类，根据行/列聚类数量分隔热图行,cutree_rows=num分割行，cutree_cols=num分割列
         #treeheight_row = 30, treeheight_col = 30, #若行、列聚类树高度调整
         border = NA, #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         cellwidth = NA, cellheight = NA,  #表示单个单元格的宽度\高度，默认为 “NA”
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
         
)

dev.off()

NORMALcelltype_cellstate2 <- NORMAL2[,c(3:16,32)]
head(NORMALcelltype_cellstate2)
NORMALcelltype_cellstatemedian2 <- aggregate(. ~ celltype, data = NORMALcelltype_cellstate2, mean)
rownames(NORMALcelltype_cellstatemedian2) <- NORMALcelltype_cellstatemedian2$celltype
NORMALcelltype_cellstatemedian2$celltype <- NULL
pdf("/boot2/xff/cellstate2/NORMALcelltype_cellstatemedian2.pdf")
pheatmap(NORMALcelltype_cellstatemedian2, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = T, 
         #cluster_rows表示仅对行聚类，cluster_cols表示仅对列聚类，值为TRUE或FALSE
         cutree_rows = 5, cutree_cols = NA, #若进行了行/列聚类，根据行/列聚类数量分隔热图行,cutree_rows=num分割行，cutree_cols=num分割列
         #treeheight_row = 30, treeheight_col = 30, #若行、列聚类树高度调整
         border = NA, #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         cellwidth = NA, cellheight = NA,  #表示单个单元格的宽度\高度，默认为 “NA”
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
         
)

dev.off()

allcellFC2 <- data.frame()
for (i in rownames(NORMALcelltype_cellstatemedian2)) {
  print(i)
  if(i!="PBMC"){
    singleFC <- foldchangecell(i)
    allcellFC2 <- rbind(allcellFC2,singleFC)
  }
}
allcellFC2

pdf("/boot2/xff/cellstate2/celltypeFC2.pdf")
pheatmap(allcellFC2, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = T, 
         #cluster_rows表示仅对行聚类，cluster_cols表示仅对列聚类，值为TRUE或FALSE
         cutree_rows = 4, cutree_cols = NA, #若进行了行/列聚类，根据行/列聚类数量分隔热图行,cutree_rows=num分割行，cutree_cols=num分割列
         #treeheight_row = 30, treeheight_col = 30, #若行、列聚类树高度调整
         border = NA, #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         cellwidth = NA, cellheight = NA,  #表示单个单元格的宽度\高度，默认为 “NA”
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
         
)

dev.off()


save(TUMOR,TUMOR2,NORMAL,NORMAL2,file = "/boot2/xff/cellstate2/singlecancer/tumandnor.RData")
load("/boot2/xff/cellstate2/singlecancer/tumandnor.RData")

head(TUMOR)


#单癌单细胞
a <- "THCA"
singlecancerfun <- function(a){
  singlecancer <- TUMOR[TUMOR$cancertype==a,]
  singlenormal <- NORMAL[NORMAL$cancertype==a,]
  singlecancer <- singlecancer[,c(3:16,28,32)]
  singlenormal <- singlenormal[,c(3:16,28,32)]
  head(singlecancer)
  singlecell <- names(table(singlecancer$celltype))
  singlecell2 <- names(table(singlenormal$celltype))
  singlecell3 <- intersect(singlecell,singlecell2)
  normedian <- apply(singlenormal[,1:14],2, median)
  normedian <- t(as.data.frame(normedian))
  #m <- "B_I"
  #singlecancer2 <- singlecancer[,2:11]
  singlecancer3 <- data.frame()
  for (m in singlecell) {
    if(m %in% singlecell3){
      singlenormal2 <- singlenormal[singlenormal$celltype==m,]
      normalmedian <- apply(singlenormal2[,1:14],2, median)
      normalmedian <- t(as.data.frame(normalmedian))
    }else{
      normalmedian <- normedian
    }
    singlecancer2 <- singlecancer[singlecancer$celltype==m,]
    #i <- 1
    for (i in 1:14) {
      singlecancer2[,i] <- (log2(singlecancer2[,i]+1)-log2(normalmedian[,i]+1))
    }
    singlecancer3 <- bind_rows(singlecancer3,singlecancer2)
  }
  
  head(singlecancer3)
  return(singlecancer3)
}
library(dplyr)
intersectcancer <- intersect(names(table(TUMOR$cancertype)),names(table(NORMAL$cancertype)))
#n <- "ALL"
for (n in intersectcancer[1:26]) {
  print(n)
  singlecancerheatmap <- singlecancerfun(n)
  rownames(singlecancerheatmap) <- paste0("row_", seq(nrow(singlecancerheatmap)))
  print(head(singlecancerheatmap))
  ann_row = singlecancerheatmap[,c(1,16)]
  ann_row$Angiogenesis <- NULL
  print(head(ann_row))
  row.names(ann_row) = rownames(ann_row)
  filename <- paste("/boot2/xff/cellstate2/singlecancer/",n,".pdf",sep = "")
  pdf(filename,width=18,height=12)
  pheatmap(singlecancerheatmap[,1:14], 
           scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
           cluster_rows = F,
           cluster_cols = F, 
           #cluster_rows表示仅对行聚类，cluster_cols表示仅对列聚类，值为TRUE或FALSE
           cutree_rows = NA, cutree_cols = NA, #若进行了行/列聚类，根据行/列聚类数量分隔热图行,cutree_rows=num分割行，cutree_cols=num分割列
           #treeheight_row = 30, treeheight_col = 30, #若行、列聚类树高度调整
           border = NA, #表示热图每个小的单元格边框的颜色，默认为 "grey60"
           cellwidth = NA, cellheight = NA,  #表示单个单元格的宽度\高度，默认为 “NA”
           show_rownames = F, show_colnames = T, #表示是否显示行名、列名
           main = "cellstateFC", #表示热图的标题名字
           color = custom_colors(num_breaks),
           breaks = breaks,
           angle_col = "45", #表示列标签的角度
           gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
           #gaps_col = c(1,2,3,4,5,6),  #仅在未进行列聚类时使用，表示在列方向上热图的隔断位置
           annotation_row = ann_row, 
           annotation_col = NA, #表示是否对行、列进行注释，默认NA
           annotation = NA, 
           #annotation_colors = list(Celltype = ann_color),  #表示行注释及列注释的颜色，默认NA
           annotation_legend = TRUE, #表示是否显示注释的图例信息
           annotation_names_row = TRUE, annotation_names_col = TRUE) #表示是否显示行、列注释的名称
  dev.off()
}

#按细胞类型取上述数据框的中位数
a <- n
singlecancerfun2 <- function(a){
  singlecancer <- TUMOR[TUMOR$cancertype==a,]
  singlenormal <- NORMAL[NORMAL$cancertype==a,]
  singlecancer <- singlecancer[,c(3:16,28,32)]
  singlenormal <- singlenormal[,c(3:16,28,32)]
  head(singlecancer)
  singlecell <- names(table(singlecancer$celltype))
  singlecell2 <- names(table(singlenormal$celltype))
  singlecell3 <- intersect(singlecell,singlecell2)
  normedian <- apply(singlenormal[,1:14],2, median)
  normedian <- t(as.data.frame(normedian))
  #m <- "B_I"
  #singlecancer2 <- singlecancer[,2:11]
  singlecancer3 <- data.frame()
  for (m in singlecell) {
    if(m %in% singlecell3){
      singlenormal2 <- singlenormal[singlenormal$celltype==m,]
      normalmedian <- apply(singlenormal2[,1:14],2, median)
      normalmedian <- t(as.data.frame(normalmedian))
    }else{
      normalmedian <- normedian
    }
    singlecancer2 <- singlecancer[singlecancer$celltype==m,]
    #singlecancer2==singlenormal2
    tumormedian <- apply(singlecancer2[,1:14],2, median)
    tumormedian <- t(as.data.frame(tumormedian))
    
    singlecancermedian <- (log2(tumormedian+1)-log2(normalmedian+1))
    singlecancermedian <- as.data.frame(singlecancermedian)
    rownames(singlecancermedian) <- m
    
    singlecancer3 <- bind_rows(singlecancer3,singlecancermedian)
  }
  head(singlecancer3)
  return(singlecancer3)
}

#绘制扇形图
library(dplyr)
intersectcancer <- intersect(names(table(TUMOR$cancertype)),names(table(NORMAL$cancertype)))
intersectcelltype<- names(table(NORMAL$celltype))
#n <- "ALL"
#Btest <- result_df
Btest[] <- 0
Btest <- Btest[1:5,]
rownames(Btest) <- c("-2","-1","0","1","2")
class(Btest)
cellstatename <- colnames(Btest)

for (n in intersectcancer[1:26]) {
  #print(n)
  data_df <- singlecancerfun2(n)
  result_df <- ifelse(data_df > 0.2, 2, 
                      ifelse(data_df > 0 & data_df <= 0.2, 1, 
                             ifelse(data_df == 0, 0, 
                                    ifelse(data_df > -0.2 & data_df < 0, -1, -2))))
  result_df <- as.data.frame(result_df)
  rownames(result_df) <- rownames(data_df)
  colnames(result_df) <- colnames(data_df)
  if(!all(is.na(result_df["B_I",]))){
    print(n)
    print("this have")
    #l <- 1
    for (l in 1:14) {
      if(result_df["B_I",l]==-2){
        Btest["-2",l]=Btest["-2",l]+1
      }else if(result_df["B_I",l]==-1){
        Btest["-1",l]=Btest["-1",l]+1
      }else if(result_df["B_I",l]==0){
        Btest["0",l]=Btest["0",l]+1
      }else if(result_df["B_I",l]==1){
        Btest["1",l]=Btest["1",l]+1
      }else if(result_df["B_I",l]==2){
        Btest["2",l]=Btest["2",l]+1
      }
    }
  }
}
library(ggplot2)
Btest$group <- rownames(Btest)

Btest$color <- c("steelblue", "#87d6e8", "white", "#d0adb8", "#e99593")
Btest$radius <- c(0.4, 0.7, 1, 1.3, 1.6)


Btest$`Cell Cycle`
#"Angiogenesis"    "Apoptosis"       "Cell Cycle"      "Differentiation" "DNA damage"     
#"DNA repair"      "EMT"             "Hypoxia"         "Inflammation"    "Invasion"       
#"Metastasis"      "Proliferation"   "Quiescence"      "Stemness"

pdf("/boot2/xff/cellstate2/pan/BStemness.pdf")
ggplot(Btest, aes(x = 1, y = Stemness, fill = color),
                               stat = "identity", width = 1, color = "black") +
  coord_polar(theta = "y") +
  theme_void() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = Btest$color) +
  guides(fill = guide_legend(title = "Group"))+ 
  geom_rect(aes(xmin = -radius, xmax = radius,
                       ymin = lag(cumsum(Stemness), default = 0),
                       ymax = cumsum(Stemness)),
                   color = "black", fill = Btest$color)



dev.off()

#循环画图
m <- "Acinar_O"
plots <- list()
for (m in intersectcelltype) {
  print(m)
  Btest[] <- 0
  
  for (n in intersectcancer[1:26]) {
    #print(n)
    data_df <- singlecancerfun2(n)
    result_df <- ifelse(data_df > 0.2, 2, 
                        ifelse(data_df > 0 & data_df <= 0.2, 1, 
                               ifelse(data_df == 0, 0, 
                                      ifelse(data_df > -0.2 & data_df < 0, -1, -2))))
    result_df <- as.data.frame(result_df)
    rownames(result_df) <- rownames(data_df)
    colnames(result_df) <- colnames(data_df)
    if(!all(is.na(result_df[m,]))){
      print(n)
      #l <- 1
      for (l in 1:14) {
        if(result_df[m,l]==-2){
          Btest["-2",l]=Btest["-2",l]+1
        }else if(result_df[m,l]==-1){
          Btest["-1",l]=Btest["-1",l]+1
        }else if(result_df[m,l]==0){
          Btest["0",l]=Btest["0",l]+1
        }else if(result_df[m,l]==1){
          Btest["1",l]=Btest["1",l]+1
        }else if(result_df[m,l]==2){
          Btest["2",l]=Btest["2",l]+1
        }
      }
    }
  }
  Btest$group <- rownames(Btest)
  Btest$color <- c("steelblue", "#87d6e8", "white", "#d0adb8", "#e99593")
  Btest$radius <- c(0.4, 0.7, 1, 1.3, 1.6)
  
  for (k in cellstatename) {
    filename <- paste("/boot2/xff/cellstate2/pan/",m,"/",k,".pdf",sep = "")
    pdf(filename)
    p1 <- ggplot(Btest, aes(x = 1, y = !!sym(k), fill = color)) +
      geom_rect(aes(xmin = -radius, xmax = radius,
                    ymin = lag(cumsum(!!sym(k)), default = 0),
                    ymax = cumsum(!!sym(k))),
                color = "black", fill = Btest$color) +
      coord_polar(theta = "y") +
      theme_void() +
      theme(legend.position = "bottom") +
      scale_fill_manual(values = Btest$color) +
      guides(fill = guide_legend(title = "Group"))
    print(p1)
    dev.off()
    plots[[length(plots) + 1]] <- p1
  }
}



library(gridExtra)
library(grid)

pdf("/boot2/xff/cellstate2/test/pan.pdf")
grid_arrange <- do.call(grid.arrange, c(list=plots, ncol = 14))
dev.off()

#饼图按癌型分块

cancertest <- data.frame(row.names = intersectcancer, matrix(NA, nrow = length(intersectcancer), ncol = length(cellstatename)))
colnames(cancertest) <- cellstatename
cancertest <- cancertest[!(cancertest$group %in% c("ALL", "AML", "BLCA", "Glioma", "HNSC", "KIPAN", "MCC", "MF", "MM", "NPC","UVM")), ]
cancertest2 <- cancertest
intersectcancer2 <- rownames(cancertest2)

plots <- list()

for (m in intersectcelltype) {
  print(m)
  cancertest[] <- 0
  cancertest2[] <- 0
  for (n in intersectcancer2) {
    #print(n)
    data_df <- singlecancerfun2(n)
    
    if(!all(is.na(data_df[m,]))){
      print(n)
      #l <- 1
      cancertest[n,] <- 1
      cancertest2[n,] <- data_df[m,]
    }
  }
  cancertest$group <- rownames(cancertest)
  cancertest$color <- c("rosybrown", "goldenrod", "olivedrab", "darkcyan", "darkmagenta",
                        "brown","orange","seagreen","steelblue","mediumvioletred",
                        "indianred","gold","palegreen","dodgerblue","darkblue"
                        )
  
  #min(cancertest2$Stemness)
  for (k in cellstatename) {
    min_values <- unique(sort(cancertest2[,k]))
    radius_values <- seq(0.4, length.out = length(min_values), by = 0.1)
    to_radius <- setNames(radius_values, min_values)
    cancertest2$radius <- to_radius[as.character(cancertest2[,k])]
    cancertest$radius <- cancertest2$radius
    filename <- paste("/boot2/xff/cellstate2/pan2/",m,"/",k,".pdf",sep = "")
    pdf(filename)
    p1 <- ggplot(cancertest, aes(x = 1, y = !!sym(k), fill = color)) +
      geom_rect(aes(xmin = -radius, xmax = radius,
                    ymin = lag(cumsum(!!sym(k)), default = 0),
                    ymax = cumsum(!!sym(k))),
                color = "black", fill = cancertest$color) +
      coord_polar(theta = "y") +
      theme_void() +
      theme(legend.position = "bottom") +
      scale_fill_manual(values = cancertest$color) +
      guides(fill = guide_legend(title = "Group"))
    print(p1)
    dev.off()
    plots[[length(plots) + 1]] <- p1
  }
}

pdf("/boot2/xff/cellstate2/test/pan2.pdf")
grid_arrange <- do.call(grid.arrange, c(list=plots, ncol = 14))
dev.off()



#行是癌型，列是状态
intersectcelltype2 <- c("B_I","CD4Tconv_I","CD8T_I","CD8Tex_I","Tprolif_I","Treg_I",
                        "DC_I","Mast_I","Mono_Macro_I","NK_I","Malignant_M","Endothelial_S",
                        "Epithelial_S","Fibroblasts_S","Myofibroblasts_S")
celltest <- data.frame(row.names = intersectcelltype2, matrix(NA, nrow = length(intersectcelltype2), ncol = length(cellstatename)))
colnames(celltest) <- cellstatename
celltest2 <- celltest

plots <- list()
#m <- "BCC"
for (m in intersectcancer2) {
  print(m)
  celltest[] <- 0
  celltest2[] <- 0
  #n <- "B_I"
  for (n in intersectcelltype2) {
    #print(n)
    data_df <- singlecancerfun2(m)
    
    if(!all(is.na(data_df[n,]))){
      print(n)
      #l <- 1
      celltest[n,] <- 1
      celltest2[n,] <- data_df[n,]
    }
  }
  celltest$group <- rownames(celltest)
  celltest$color <- c("rosybrown", "goldenrod", "olivedrab", "darkcyan", "darkmagenta",
                        "brown","orange","seagreen","steelblue","mediumvioletred",
                        "indianred","gold","palegreen","dodgerblue","darkblue"
  )
  
  #min(celltest2$Stemness)
  for (k in cellstatename) {
    min_values <- unique(sort(celltest2[,k]))
    radius_values <- seq(0.4, length.out = length(min_values), by = 0.1)
    to_radius <- setNames(radius_values, min_values)
    celltest2$radius <- to_radius[as.character(celltest2[,k])]
    celltest$radius <- celltest2$radius
    filename <- paste("/boot2/xff/cellstate2/pancancer_cellstate/",m,"/",k,".pdf",sep = "")
    pdf(filename)
    p1 <- ggplot(celltest, aes(x = 1, y = !!sym(k), fill = color)) +
      geom_rect(aes(xmin = -radius, xmax = radius,
                    ymin = lag(cumsum(!!sym(k)), default = 0),
                    ymax = cumsum(!!sym(k))),
                color = "black", fill = celltest$color) +
      coord_polar(theta = "y") +
      theme_void() +
      theme(legend.position = "bottom") +
      scale_fill_manual(values = celltest$color) +
      guides(fill = guide_legend(title = "Group"))
    print(p1)
    dev.off()
    plots[[length(plots) + 1]] <- p1
  }
}

pdf("/boot2/xff/cellstate2/test/pancancer_cellstate.pdf")
grid_arrange <- do.call(grid.arrange, c(list=plots, ncol = 14))
dev.off()



#计算不同状态的皮尔森相关系数
#先去除罕见癌症和特异性细胞类型
singlecancer <- TUMOR2[,c(3:16,28,32)]
singlenormal <- NORMAL2[,c(3:16,28,32)]
singlecell <- names(table(singlecancer$celltype))
singlecell2 <- names(table(singlenormal$celltype))
singlecell3 <- intersect(singlecell,singlecell2)
singlecell3 <- c("B_I","CD4Tconv_I","CD8T_I","CD8Tex_I","DC_I",
                 "Endothelial_S","Epithelial_S","Fibroblasts_S","Malignant_M","Mast_I",
                 "Mono/Macro_I","Myofibroblasts_S","NK_I","Tprolif_I","Treg_I")
singlecancer <- TUMOR2[TUMOR2$celltype %in% c("B_I","CD4Tconv_I","CD8T_I","CD8Tex_I","DC_I",
                                              "Endothelial_S","Epithelial_S","Fibroblasts_S","Malignant_M","Mast_I",
                                              "Mono/Macro_I","Myofibroblasts_S","NK_I","Tprolif_I","Treg_I"),c(3:16,28,32)]
singlenormal <- NORMAL2[NORMAL2$celltype %in% c("B_I","CD4Tconv_I","CD8T_I","CD8Tex_I","DC_I",
                                               "Endothelial_S","Epithelial_S","Fibroblasts_S","Malignant_M","Mast_I",
                                               "Mono/Macro_I","Myofibroblasts_S","NK_I","Tprolif_I","Treg_I"),c(3:16,28,32)]
singlecan <- names(table(singlecancer$cancertype))
singlecan2 <- names(table(singlenormal$cancertype))
singlecan3 <- intersect(singlecan,singlecan2)
singlecan3 <- c("BCC","BRCA","CESC","CRC","ESCA",
                "KIRC","LIHC","NSCLC","OV","PAAD",
                "PRAD","SKCM","STAD","THCA","UCEC")

singlecancer <- singlecancer[singlecancer$cancertype %in% c("BCC","BRCA","CESC","CRC","ESCA",
                                                            "KIRC","LIHC","NSCLC","OV","PAAD",
                                                            "PRAD","SKCM","STAD","THCA","UCEC"),]
singlenormal <- singlenormal[singlenormal$cancertype %in% c("BCC","BRCA","CESC","CRC","ESCA",
                                                            "KIRC","LIHC","NSCLC","OV","PAAD",
                                                            "PRAD","SKCM","STAD","THCA","UCEC"),]


normedian <- apply(singlenormal[,1:14],2, median)
normedian <- t(as.data.frame(normedian))
#m <- "B_I"
#singlecancer2 <- singlecancer[,2:11]
singlecancer3 <- data.frame()
for (m in singlecell) {
  if(m %in% singlecell3){
    singlenormal2 <- singlenormal[singlenormal$celltype==m,]
    normalmedian <- apply(singlenormal2[,1:14],2, median)
    normalmedian <- t(as.data.frame(normalmedian))
  }else{
    normalmedian <- normedian
  }
  singlecancer2 <- singlecancer[singlecancer$celltype==m,]
  #i <- 1
  for (i in 1:14) {
    singlecancer2[,i] <- (log2(singlecancer2[,i]+1)-log2(normalmedian[,i]+1))
  }
  singlecancer3 <- bind_rows(singlecancer3,singlecancer2)
}
singlemt <- singlecancer3[,1:14]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
library(corrplot)
pdf("/boot2/xff/cellstate2/cor/cor.pdf")
corrplot(cor_pearson, method = 'number', number.cex = 0.8, diag = FALSE, tl.cex = 0.8,col = my_colors)
corrplot(cor_pearson, add = TRUE, type = 'upper', method = 'pie', diag = FALSE, tl.pos = 'n', cl.pos = 'n',col = my_colors)

dev.off()

hc <- hclust(dist(cor_pearson))
pdf("/boot2/xff/cellstate2/cor/cor2.pdf")
my_colors <- colorRampPalette(c("steelblue", "white", "#e99593"))(100)
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
rect.hclust(hc, k = 6) 
rect.hclust(hc, k = 5) 
dev.off()
head(singlecancer3)

dev.new()
#每种细胞类型计算皮尔森相关系数
everycellperson <- data.frame()
for (m in singlecell3) {
  if(m %in% singlecell3){
    singlenormal2 <- singlenormal[singlenormal$celltype==m,]
    normalmedian <- apply(singlenormal2[,1:14],2, median)
    normalmedian <- t(as.data.frame(normalmedian))
  }else{
    normalmedian <- normedian
  }
  singlecancer2 <- singlecancer[singlecancer$celltype==m,]
  #i <- 1
  for (i in 1:14) {
    singlecancer2[,i] <- (log2(singlecancer2[,i]+1)-log2(normalmedian[,i]+1))
  }
  everycellperson <- singlecancer2[,1:14]
  everycellperson <- scale(everycellperson)
  cor_pearson <- cor(everycellperson, method = 'pearson')
  cor_pearson
  
  myfilename <- paste("/boot2/xff/cellstate2/cor/everycell/",m,"cor.pdf",sep = "")
  pdf(myfilename)
  corrplot(cor_pearson, method = 'number', number.cex = 0.8, diag = FALSE, tl.cex = 0.8,col = my_colors)
  corrplot(cor_pearson, add = TRUE, type = 'upper', method = 'pie', diag = FALSE, tl.pos = 'n', cl.pos = 'n',col = my_colors)
  
  dev.off()
  
}
#聚类
for (m in singlecell3) {
  if(m %in% singlecell3){
    singlenormal2 <- singlenormal[singlenormal$celltype==m,]
    normalmedian <- apply(singlenormal2[,1:14],2, median)
    normalmedian <- t(as.data.frame(normalmedian))
  }else{
    normalmedian <- normedian
  }
  singlecancer2 <- singlecancer[singlecancer$celltype==m,]
  #i <- 1
  for (i in 1:14) {
    singlecancer2[,i] <- (log2(singlecancer2[,i]+1)-log2(normalmedian[,i]+1))
  }
  everycellperson <- singlecancer2[,1:14]
  everycellperson <- scale(everycellperson)
  cor_pearson <- cor(everycellperson, method = 'pearson')
  cor_pearson
  hc <- hclust(dist(cor_pearson))
  myfilename <- paste("/boot2/xff/cellstate2/cor/everycell2/",m,"cor.pdf",sep = "")
  pdf(myfilename)
  corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = FALSE, diag = FALSE,col = my_colors)
  rect.hclust(hc, k = 6) 
  rect.hclust(hc, k = 5) 
  dev.off()
  
}

#每种癌症类型计算皮尔森相关系数
everycancerperson <- data.frame()
#n <- "BCC"
for (n in singlecan3) {
  singlecancer2 <- singlecancer[singlecancer$cancertype==n,]
  singlenormal2 <- singlenormal[singlenormal$cancertype==n,]
  normalmedian <- apply(singlenormal2[,1:14],2, median)
  normalmedian <- t(as.data.frame(normalmedian))
      #i <- 1
  for (i in 1:14) {
    singlecancer2[,i] <- (log2(singlecancer2[,i]+1)-log2(normalmedian[,i]+1))
  }
  everycancerperson <- singlecancer2[,1:14]
  everycancerperson <- scale(everycancerperson)
  cor_pearson <- cor(everycancerperson, method = 'pearson')
  cor_pearson
  myfilename <- paste("/boot2/xff/cellstate2/cor/everycan/",n,"cor.pdf",sep = "")
  pdf(myfilename)
  corrplot(cor_pearson, method = 'number', number.cex = 0.8, diag = FALSE, tl.cex = 0.8)
  corrplot(cor_pearson, add = TRUE, type = 'upper', method = 'pie', diag = FALSE, tl.pos = 'n', cl.pos = 'n')
  dev.off()
}
#聚类
for (n in singlecan3) {
  singlecancer2 <- singlecancer[singlecancer$cancertype==n,]
  singlenormal2 <- singlenormal[singlenormal$cancertype==n,]
  normalmedian <- apply(singlenormal2[,1:14],2, median)
  normalmedian <- t(as.data.frame(normalmedian))
  #i <- 1
  for (i in 1:14) {
    singlecancer2[,i] <- (log2(singlecancer2[,i]+1)-log2(normalmedian[,i]+1))
  }
  everycancerperson <- singlecancer2[,1:14]
  everycancerperson <- scale(everycancerperson)
  cor_pearson <- cor(everycancerperson, method = 'pearson')
  cor_pearson
  hc <- hclust(dist(cor_pearson))
  myfilename <- paste("/boot2/xff/cellstate2/cor/everycan2/",n,"cor.pdf",sep = "")
  pdf(myfilename)
  corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = FALSE, diag = FALSE,col = my_colors)
  rect.hclust(hc, k = 6) 
  rect.hclust(hc, k = 5) 
  dev.off()
}


#查看不同model的分布情况
fig1gcell<- read.csv("/boot2/xff/cellstate2/cor/fig1g1.csv",header = T)
fig1gcancer <- read.csv("/boot2/xff/cellstate2/cor/fig1g2.csv",header = T)
rownames(fig1gcell) <- fig1gcell$X
fig1gcell$X <- NULL
rownames(fig1gcancer) <- fig1gcancer$X
fig1gcancer$X <- NULL

library(pheatmap)
library(RColorBrewer)
custom_colors2 <- colorRampPalette(c("steelblue","#9FBED8","#FE9F79","#A5607C"))
num_breaks <- 100
breaks2 <- seq(0.25, 1, length.out = num_breaks + 1)

pdf("/boot2/xff/cellstate2/cor/fig1g1.pdf")
pheatmap(fig1gcell, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = F, 
         border = "grey60", #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         cellwidth = NA, cellheight = NA,  #表示单个单元格的宽度\高度，默认为 “NA”
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors2(num_breaks),
         breaks = breaks2,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
)

dev.off()

#绘制气泡热图
fig1gcell<- read.csv("/boot2/xff/cellstate2/cor/fig1g1.csv",header = T)
fig1gcancer <- read.csv("/boot2/xff/cellstate2/cor/fig1g2.csv",header = T)
library(tidyr)
long_cell <- gather(fig1gcell, key = "Model", value = "Value", -X)
long_cancer <- gather(fig1gcancer, key = "Model", value = "Value", -X)

pdf("/boot2/xff/cellstate2/cor/fig1g1.pdf")
ggplot(long_cancer, aes(x = Model, y = X, fill = Value,size =Value)) +
  geom_point(shape = 21, color = "white") +  # 使用形状21绘制气泡
  scale_size_continuous(range = c(2,6)) +  # 设置气泡大小范围
  scale_fill_gradient2(low = "steelblue",mid="white", high = "#A5607C",midpoint = 0.625, name = "Value") +  # 设置颜色梯度
  theme_minimal() +  # 使用简洁主题
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),  # 去除主要网格线
        panel.grid.minor = element_blank())  # 旋转X轴标签以提高可读性
dev.off()

#看每一种模式如何协调变化
cellmodel <- singlecancer3[,-15]
cancermodel <- singlecancer3[,-16]
cellmodelmean <- aggregate(. ~ celltype, data = cellmodel, mean)
cancermodelmean <- aggregate(. ~ cancertype, data = cancermodel, mean)

custom_colors <- colorRampPalette(c("steelblue","white","#e99593"))
num_breaks <- 100
breaks <- seq(-0.2, 0.2, length.out = num_breaks + 1)

#model1
cellmodel1 <- cellmodelmean[,c("celltype","EMT","Invasion","Angiogenesis","Metastasis")]
cancermodel1 <- cancermodelmean[,c("cancertype","EMT","Invasion","Angiogenesis","Metastasis")]
rownames(cellmodel1) <- cellmodel1$celltype
rownames(cancermodel1) <- cancermodel1$cancertype
cellmodel1$celltype <- NULL
cancermodel1$cancertype <- NULL

pheatmap(cellmodel1, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = F, 
         border = "grey60", #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
)



pheatmap(cancermodel1, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = F, 
         border = "grey60", #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
)

#model2
cellmodel2 <- cellmodelmean[,c("celltype","Apoptosis")]
cancermodel2 <- cancermodelmean[,c("cancertype","Apoptosis")]
rownames(cellmodel2) <- cellmodel2$celltype
rownames(cancermodel2) <- cancermodel2$cancertype
cellmodel2$celltype <- NULL
cancermodel2$cancertype <- NULL

pheatmap(cellmodel2, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = F, 
         border = "grey60", #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
)
pheatmap(cancermodel2, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = F, 
         border = "grey60", #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
)

#model3
cellmodel3 <- cellmodelmean[,c("celltype","Hypoxia")]
cancermodel3 <- cancermodelmean[,c("cancertype","Hypoxia")]
rownames(cellmodel3) <- cellmodel3$celltype
rownames(cancermodel3) <- cancermodel3$cancertype
cellmodel3$celltype <- NULL
cancermodel3$cancertype <- NULL

pheatmap(cellmodel3, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = F, 
         border = "grey60", #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
)
pheatmap(cancermodel3, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = F, 
         border = "grey60", #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
)

#model4
cellmodel4 <- cellmodelmean[,c("celltype","Inflammation","Quiescence")]
cancermodel4 <- cancermodelmean[,c("cancertype","Inflammation","Quiescence")]
rownames(cellmodel4) <- cellmodel4$celltype
rownames(cancermodel4) <- cancermodel4$cancertype
cellmodel4$celltype <- NULL
cancermodel4$cancertype <- NULL

pheatmap(cellmodel4, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = F, 
         border = "grey60", #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
)
pheatmap(cancermodel4, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = F, 
         border = "grey60", #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
)

#model5
cellmodel5 <- cellmodelmean[,c("celltype","Cell Cycle","DNA damage","DNA repair")]
cancermodel5 <- cancermodelmean[,c("cancertype","Cell Cycle","DNA damage","DNA repair")]
rownames(cellmodel5) <- cellmodel5$celltype
rownames(cancermodel5) <- cancermodel5$cancertype
cellmodel5$celltype <- NULL
cancermodel5$cancertype <- NULL

pheatmap(cellmodel5, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = F, 
         border = "grey60", #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
)
pheatmap(cancermodel5, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = F, 
         border = "grey60", #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
)

#model6
cellmodel6 <- cellmodelmean[,c("celltype","Proliferation")]
cancermodel6 <- cancermodelmean[,c("cancertype","Proliferation")]
rownames(cellmodel6) <- cellmodel6$celltype
rownames(cancermodel6) <- cancermodel6$cancertype
cellmodel6$celltype <- NULL
cancermodel6$cancertype <- NULL

pheatmap(cellmodel6, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = F, 
         border = "grey60", #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
)
pheatmap(cancermodel6, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = F, 
         border = "grey60", #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
)

#model7
cellmodel7 <- cellmodelmean[,c("celltype","Differentiation","Stemness")]
cancermodel7 <- cancermodelmean[,c("cancertype","Differentiation","Stemness")]
rownames(cellmodel7) <- cellmodel7$celltype
rownames(cancermodel7) <- cancermodel7$cancertype
cellmodel7$celltype <- NULL
cancermodel7$cancertype <- NULL

pheatmap(cellmodel7, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = F, 
         border = "grey60", #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
)
pheatmap(cancermodel7, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = F, 
         border = "grey60", #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
)

#model2,3,6
cellmodel236 <- cellmodelmean[,c("celltype","Apoptosis","Hypoxia","Proliferation")]
cancermodel236 <- cancermodelmean[,c("cancertype","Apoptosis","Hypoxia","Proliferation")]
rownames(cellmodel236) <- cellmodel236$celltype
rownames(cancermodel236) <- cancermodel236$cancertype
cellmodel236$celltype <- NULL
cancermodel236$cancertype <- NULL

pheatmap(cellmodel236, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = F, 
         border = "grey60", #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
)
pheatmap(cancermodel236, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = F, 
         border = "grey60", #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
)
cellstatename

#把状态值换为model值

cellmodel1$model1 <- rowMeans(cellmodel1)
cellmean <- cellmodel1
cellmean$EMT <- NULL
cellmean$Invasion <- NULL
cellmean$Angiogenesis <- NULL
cellmean$Metastasis <- NULL
cellmean$model2 <- cellmodel236$Apoptosis
cellmean$model3 <- cellmodel236$Hypoxia
cellmean$model4 <- rowMeans(cellmodel4)
cellmean$model5 <- rowMeans(cellmodel5)
cellmean$model6 <- cellmodel236$Proliferation
cellmean$model7 <- rowMeans(cellmodel7)

cancermodel1$model1 <- rowMeans(cancermodel1)
cancermean <- cancermodel1
cancermean$EMT <- NULL
cancermean$Invasion <- NULL
cancermean$Angiogenesis <- NULL
cancermean$Metastasis <- NULL
cancermean$model2 <- cancermodel236$Apoptosis
cancermean$model3 <- cancermodel236$Hypoxia
cancermean$model4 <- rowMeans(cancermodel4)
cancermean$model5 <- rowMeans(cancermodel5)
cancermean$model6 <- cancermodel236$Proliferation
cancermean$model7 <- rowMeans(cancermodel7)
cellmean
cancermean
pheatmap(cellmean, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = T, 
         border = "grey60", #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
)
pheatmap(cancermean, 
         scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         cluster_rows = T,
         cluster_cols = T, 
         border = "grey60", #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         show_rownames = T, show_colnames = T, #表示是否显示行名、列名
         main = "cellstate", #表示热图的标题名字
         color = custom_colors(num_breaks),
         breaks = breaks,
         angle_col = "45", #表示列标签的角度
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
)

head(singlecancer3)
table(singlecancer3$celltype)
save(singlecancer3,file = "/boot2/xff/cellstate2/model/justup/singlecancer3.RData")
