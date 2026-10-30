setwd("E:/T细胞障碍/肺癌/泛癌/药物敏感性")
#基因和细胞系关系
drug_id<-read.csv("E:/T细胞障碍/肺癌/泛癌/药物敏感性/细胞系编号.csv", header = T)
data <- read.table("E:/T细胞障碍/肺癌/泛癌/药物敏感性/Cell_line_RMA_proc_basalExp.txt", header=T, sep="\t")#读入表达数据

data1<-data[!data$GENE_SYMBOLS=="",]
rownames(data1)<-data1$GENE_SYMBOLS
data1<-data1[,-c(1,2)]
names(data1)<- sub("^.*?\\.", "",names(data1))
names(data1)<- drug_id$Line[match(names(data1),drug_id$COSMIC_ID)]

#-------------#
a<-c( "ANXA1" ,"CD3E","CD55" ,  "HLA-DRA" , "HLA-DRB1", "IL2RG"  ,  "NFKBIZ" ,  "B2M" ,     "HLA-A"  ,  "HLA-B" ,  
      "PTPRC", "EZR","GBP1" , "LCK" ,"LGALS3",   "THY1" ,    "HSPD1",    "IL6" ,     "CD3D" ,    "CD74" , 
      "CTSL", "AIF1" ,"APP"  , "CCL5" , "CORO1A",   "CXCL12",   "S100A7",   "CD24",     "HES1" ,    "HLA-DPA1",
      "HLA-DPB1","IGF1",  "IGFBP2" ,  "VCAM1"  ,  "IL7R",     "MAFB" ,    "RPS6",     "ZFP36L1" , "ZFP36L2" )
#----------------#
#------筛选免疫治疗药和靶向药---细胞系与药物得关系-------#
drug_cell<-read.csv("E:/T细胞障碍/肺癌/泛癌/药物敏感性/GDSC2_fitted_dose_response_27Oct23.csv", header = T)
drug_cell<-drug_cell[,c(5,7,9,16)]
drug_cell_1<-drug_cell[!drug_cell$TCGA_DESC=="UNCLASSIFIED",]

drug_name<-c("Lenalidomide","Venetoclax","AZD5363","MK-2206","VE821")
subset_df <- drug_cell_1[drug_cell_1$DRUG_NAME %in% drug_name,]
wide_df <-dcast(subset_df,CELL_LINE_NAME~DRUG_NAME , value.var = "LN_IC50",fun.aggregate = mean)
wide_df[is.na(wide_df)] <- 0
rownames(wide_df)<-wide_df$CELL_LINE_NAME
wide_df<-wide_df[,-1]
#------------#
#---------药物必须要一致-----------------#
b<-intersect(rownames(wide_df),colnames(data1))
wide_df<-wide_df[b,]
data1<-data1[,b]
identical(rownames(wide_df),colnames(data1))
#-----------------#
cell<-read.csv("E:/T细胞障碍/肺癌/泛癌/药物敏感性/cell_exp.csv", header = T)
rownames(cell)<-cell$X
cell<-cell[,-1]

cell<-cell[a,]
cell <- na.omit(cell)

cell<-as.matrix(cell)

#-----#
calcPhen<-calcPhenotype(trainingExprData = as.matrix(data1),
                        trainingPtype = as.matrix(wide_df),
                        testExprData = cell, # 需要matrix
                        batchCorrect = 'eb',  # 校正参数，"eb" for array, "standardize" for rnaseq
                        #IC50是对数转换的，所以表达矩阵也用对数转换过的
                        powerTransformPhenotype = F,
                        minNumSamples = 20,
                        printOutput = T,
                        removeLowVaryingGenes = 0.2,
                        removeLowVaringGenesFrom = "homogenizeData"
)
#-----------------#
pre_cell<-read.csv("E:/T细胞障碍/肺癌/泛癌/药物敏感性/calcPhenotype_Output/DrugPredictions.csv", header = T)
rownames(pre_cell)<-pre_cell$X
pre_cell<-pre_cell[,-1]
#----------------#
#
library(pheatmap)
bk <- c(seq(-15,-1,by=1),seq(1,15,by=1))
pheatmap(pre_cell, scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",#表示行、列聚类使用的度量方法，默认为“euclidean”，也可为  "correlation"即按照 Pearson correlation方法进行聚类
         
         show_rownames = T,
         
         #cluster_rows = T,
         
         #main =cancer_name[i] ,
         
         cellwidth = 15, cellheight = 15,
         
         clustering_method = "complete", #表示聚类方法，包括：‘ward’, ‘ward.D’, ‘ward.D2’, ‘single’, ‘complete’, ‘average’, ‘mcquitty’, ‘median’, ‘centroid’
         
         cluster_rows = T,cluster_cols =F , #cluster_rows表示仅对行聚类，cluster_cols表示仅对列聚类，值为TRUE或FALSE
         
         #cutree_rows = 2, #cutree_cols = , #若进行了行、列聚类，根据行、列聚类数量分隔热图行
         
         treeheight_row = 30, treeheight_col = 30,#、行、列聚类树高度调整
         
         #annotation_row = row_ann, annotation_col = NA, #表示是否对行、列进行注释，默认NA
         
         #annotation = NA, annotation_colors = ann_color,  #表示行注释及列注释的颜色，默认NA
         
         color = colorRampPalette(c('#266EAF',"white",'#B52230'))(100), #表示热图颜色,(100)表示100个等级
         #color = c(colorRampPalette(colors = c("#266EAF","white"))(length(bk)/2),colorRampPalette(colors = c("white","#B52230"))(length(bk)/2)),
         #breaks=bk
) 
