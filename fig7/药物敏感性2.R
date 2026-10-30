library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
library(TCGAbiolinks)
library(dplyr)
library(tidyverse)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)


setwd("D:/细胞状态2/论文/药物2")

drug_id<-read.csv("细胞系编号.csv", header = T)
data <- read.table("Cell_line_RMA_proc_basalExp.txt", header=T, sep="\t")#读入表达数据

data[1:3,1:3]

data1<-data[!data$GENE_SYMBOLS=="",]
rownames(data1)<-data1$GENE_SYMBOLS
data1<-data1[,-c(1,2)]
names(data1)<- sub("^.*?\\.", "",names(data1))
names(data1)<- drug_id$Line[match(names(data1),drug_id$COSMIC_ID)]
data1[1:3,1:3]

#------筛选免疫治疗药和靶向药---细胞系与药物得关系-------#
drug_cell<-read.csv("GDSC2_fitted_dose_response_27Oct23.csv", header = T)
drug_cell<-drug_cell[,c(5,7,9,16)]
drug_cell_1<-drug_cell[!drug_cell$TCGA_DESC=="UNCLASSIFIED",]
cnacer_type<-unique(drug_cell_1$TCGA_DESC)



#PAAD
#取与model1皮尔森相关系数高的基因"D:/细胞状态2/论文/空间/paad/paadcor.csv"
a<-c("COL1A1","FN1","MMP2","MMP14","BGN","INHBA")

#---选择癌症-----#
paadcancer<-drug_cell_1[drug_cell_1$TCGA_DESC=="PAAD",]
paadcancer<-paadcancer[,-2]
paadcancerwide_df <-dcast(paadcancer,CELL_LINE_NAME~DRUG_NAME , value.var = "LN_IC50",fun.aggregate = mean)
#paadcancerwide_df <-reshape2::dcast(paadcancer,CELL_LINE_NAME~DRUG_NAME , value.var = "LN_IC50",fun.aggregate = mean)
paadcancerwide_df[is.na(paadcancerwide_df)] <- 0
rownames(paadcancerwide_df)<-paadcancerwide_df$CELL_LINE_NAME
paadcancerwide_df<-paadcancerwide_df[,-1]


#---------药物必须要一致-----------------#
b<-intersect(rownames(paadcancerwide_df),colnames(data1))
paadcancerwide_df<-paadcancerwide_df[b,]
data2<-data1[,b]
identical(rownames(paadcancerwide_df),colnames(data2))

genedata <- read.table("D:/细胞状态2/论文/分化/fenhua/cox/mart_export.txt", sep=",", header=TRUE)

drug_gene<- read.table("D:/细胞状态2/bulk/expression-counts/TCGA-PAAD.htseq_counts.tsv", header=T, sep="\t")#读入表达数据
drug_gene[1:3,1:3]

drug_gene$Ensembl_ID <- str_split_fixed(drug_gene$Ensembl_ID,'[.]',2)[,1]
drugnew_data <- merge(drug_gene,genedata,by.x = "Ensembl_ID",by.y = "Gene.stable.ID")
dim(drugnew_data)
col<- c("Ensembl_ID", "Gene.stable.ID.version", "Transcript.stable.ID.version", "Transcript.stable.ID") 
drugnew_df <- drugnew_data[, !(colnames(drugnew_data)%in%col)]
drugunique_df <- drugnew_df[!duplicated(drugnew_df$Gene.name),]
rownames(drugunique_df) <- drugunique_df[,"Gene.name"]
drugunique_df$Gene.name <- NULL
drugunique_df[1:3,1:3]


drug_gene1<-drugunique_df[a,]
drug_gene_mx<-as.matrix(drug_gene1)
rownames(drug_gene_mx)
class(drug_gene_mx)
class(data2)

dim(data2)  # 检查data2的维度
dim(paadcancerwide_df)  # 检查paadcancerwide_df的维度
dim(drug_gene_mx)  # 检查drug_gene_mx的维度
sum(is.na(data2))  # 检查缺失值数量
sum(is.na(paadcancerwide_df))  # 检查缺失值数量
sum(is.na(drug_gene_mx))  # 检查缺失值数量
a
data2["MMP14",]
rownames(paadcancerwide_df)

calcPhen<-calcPhenotype(trainingExprData = as.matrix(data2),
                        trainingPtype = as.matrix(paadcancerwide_df),
                        testExprData = drug_gene_mx, # 需要matrix
                        batchCorrect = 'eb',  # 校正参数，"eb" for array, "standardize" for rnaseq
                        #IC50是对数转换的，所以表达矩阵也用对数转换过的
                        powerTransformPhenotype = F,
                        minNumSamples = 5,
                        printOutput = T,
                        removeLowVaryingGenes = 0.2,
                        removeLowVaringGenesFrom = "homogenizeData"
)

grug_pre<-read.csv("calcPhenotype_Output/paadDrugPredictions.csv", header = T)

load("D:/细胞状态2/论文/fig2/totalESTclin2.RData")
group<-totalESTclin2
head(group)
paadgroup <- totalESTclin2[totalESTclin2$cancertype=="PAAD",]
table(paadgroup$group1)
paadgroup$group1 <- ifelse(paadgroup$model1> median(paadgroup$model1),"high" ,"low")


#分组
paadgroup<-paadgroup[,c("bcr_sample_barcode","group1")]
table(paadgroup$group1)

head(paadgroup)
paadgroup$bcr_sample_barcode <- gsub("[-]", ".",paadgroup$bcr_sample_barcode)



paadDrugPred <- grug_pre
rownames(paadDrugPred)<-paadDrugPred$X
paadDrugPred<-paadDrugPred[,-1]


# 计算每列中NA值的数量
na_count <- colSums(is.na(paadDrugPred))

# 找到不含NA值的列的索引
clean_columns <- which(na_count == 0)

# 创建一个不含NA值的数据框
clean_df <- paadDrugPred[, clean_columns]
# 删除数据框中值大于8的列
clean_df <- clean_df[, !apply(clean_df, 2, function(x) any(x > 10))]

expression_matrix <- NULL

data1_df<-data.frame(t(clean_df))
#--------#
a<-intersect(colnames(data1_df),paadgroup$bcr_sample_barcode)
#------#
data1_df<-data1_df[,a]
drug_group<-paadgroup[paadgroup$bcr_sample_barcode%in% a,]


#差异分析
paaddftest <- t(data1_df)
paaddftest <- as.data.frame(paaddftest)
paaddftest$sample <- rownames(paaddftest)
paaddftest2 <- merge(paaddftest,drug_group,by.x="sample",by.y="bcr_sample_barcode")
head(paaddftest2)
dim(paaddftest2)
paaddftest4 <- paaddftest2[,c(-1,-287)]
head(paaddftest4)


paadnewFC <- data.frame(row.names = colnames(paaddftest4))
paadnewFC$FC <- NA
paadnewFC$t <- NA
paadnewFC$wilcox <- NA
paadnewFC$aov <- NA
paadnewFC$highmean <- NA


for (i in colnames(paaddftest4)) {
  print(i)
  yao <- paaddftest2[,c(i,"group1")]
  head(yao)
  yao$value <- exp(yao[,1])
  mean_high <- mean(yao[yao$group1 == "high","value"])
  mean_low <- mean(yao[yao$group1 == "low","value"])
  paadnewFC[i,"highmean"] <- mean_high
  FC <- mean_high / mean_low
  paadnewFC[i,"FC"] <- FC
  t_test_result <- t.test(value ~ group1, data = yao)
  wilcox_test_result <- wilcox.test(value ~ group1, data = yao)
  aov_result <- aov(value ~ group1, data = yao)
  paadnewFC[i,"t"] <- t_test_result$p.value
  paadnewFC[i,"wilcox"] <- wilcox_test_result$p.value
  paadnewFC[i,"aov"] <- summary(aov_result)[[1]]$"Pr(>F)"[1]
  getwd()
  myfilename <- paste("./药物预测/boxplot/paad/",i,".pdf",sep = "")
  pdf(myfilename)
  p <- ggplot(data = yao, aes(x =group1 , y =value , fill = group1)) +
    geom_violin(position = position_dodge(0.8))+
    geom_boxplot(width=0.5,
                 outlier.colour = NA,
                 position = position_dodge(0.8))+  #每组的小提琴图之间的距离
    stat_compare_means(aes(group = group1))+  #,label = 'p.signif'
    scale_fill_manual(values = c("#c11e21","steelblue") )+
    theme_bw()+theme(panel.grid=element_blank())+
    labs( x="drug", y="IC50 score")+
    theme(axis.text.x = element_text(size = 6,angle = 0, hjust = 0.5, vjust = 0.5))
  print(p)
  dev.off()
}

paadnewFC$LN_IC50 <- log(paadnewFC$highmean)
head(paadnewFC)
#paadnewFC[paadnewFC$FC>2,]
paadnewFC[paadnewFC$FC<0.5,]

dim(paadnewFC[paadnewFC$t<0.05,])
dim(paadnewFC[paadnewFC$wilcox<0.05,])
dim(paadnewFC[paadnewFC$aov<0.05,])


paadnewFC <- paadnewFC[paadnewFC$t<0.05&paadnewFC$wilcox<0.05&paadnewFC$aov<0.05,]
paadnewFC <- paadnewFC[paadnewFC$FC<1,]
paadnewFC[order(paadnewFC$highmean), ]
write.csv(paadnewFC[order(paadnewFC$highmean), ],file = "./药物预测/paad.csv")

rownames(paadnewFC[order(paadnewFC$highmean), ])


maxmindftest <- dftest4
min(colMeans(maxmindftest))
maxmindftest <- rbind(maxmindftest,colMeans(maxmindftest))
maxmindftestT <- t(maxmindftest)
sort(colMeans(maxmindftest))

#Romidepsin   Tanespimycin      Bortezomib  Dactinomycin   Paclitaxel
















#STAD
setwd("D:/细胞状态2/论文/药物2")

drug_id<-read.csv("细胞系编号.csv", header = T)
data <- read.table("Cell_line_RMA_proc_basalExp.txt", header=T, sep="\t")#读入表达数据

data[1:3,1:3]

data1<-data[!data$GENE_SYMBOLS=="",]
rownames(data1)<-data1$GENE_SYMBOLS
data1<-data1[,-c(1,2)]
names(data1)<- sub("^.*?\\.", "",names(data1))
names(data1)<- drug_id$Line[match(names(data1),drug_id$COSMIC_ID)]
data1[1:3,1:3]

#------筛选免疫治疗药和靶向药---细胞系与药物得关系-------#
drug_cell<-read.csv("GDSC2_fitted_dose_response_27Oct23.csv", header = T)
drug_cell<-drug_cell[,c(5,7,9,16)]
drug_cell_1<-drug_cell[!drug_cell$TCGA_DESC=="UNCLASSIFIED",]
cnacer_type<-unique(drug_cell_1$TCGA_DESC)



#STAD
#取与model1皮尔森相关系数高的基因"D:/细胞状态2/论文/空间/stad/stadcor.csv"
a<-c("CXCL12","LGALS1","COL6A2","PDPN","EGR1","CD44")
a <- c("FN1","COL6A3","SPARC","COL1A1","COL1A2",
       "COL3A1","LGALS1","BGN","CXCL6" ,
       "MMP3" ,"FCRLA" ,"FDCSP", "CXCL13", "MS4A1", 
       "TCL1A", "VPREB3", "PRSS21", "DPEP3", "COL9A3" 
       )
#---选择癌症-----#
stadcancer<-drug_cell_1[drug_cell_1$TCGA_DESC=="STAD",]
stadcancer<-stadcancer[,-2]
stadcancerwide_df <-dcast(stadcancer,CELL_LINE_NAME~DRUG_NAME , value.var = "LN_IC50",fun.aggregate = mean)
#stadcancerwide_df <-reshape2::dcast(stadcancer,CELL_LINE_NAME~DRUG_NAME , value.var = "LN_IC50",fun.aggregate = mean)
stadcancerwide_df[is.na(stadcancerwide_df)] <- 0
rownames(stadcancerwide_df)<-stadcancerwide_df$CELL_LINE_NAME
stadcancerwide_df<-stadcancerwide_df[,-1]


#---------药物必须要一致-----------------#
b<-intersect(rownames(stadcancerwide_df),colnames(data1))
stadcancerwide_df<-stadcancerwide_df[b,]
data2<-data1[,b]
identical(rownames(stadcancerwide_df),colnames(data2))

genedata <- read.table("D:/细胞状态2/论文/分化/fenhua/cox/mart_export.txt", sep=",", header=TRUE)

drug_gene<- read.table("D:/细胞状态2/bulk/expression-counts/TCGA-STAD.htseq_counts.tsv", header=T, sep="\t")#读入表达数据
drug_gene[1:3,1:3]

drug_gene$Ensembl_ID <- str_split_fixed(drug_gene$Ensembl_ID,'[.]',2)[,1]
drugnew_data <- merge(drug_gene,genedata,by.x = "Ensembl_ID",by.y = "Gene.stable.ID")
dim(drugnew_data)
col<- c("Ensembl_ID", "Gene.stable.ID.version", "Transcript.stable.ID.version", "Transcript.stable.ID") 
drugnew_df <- drugnew_data[, !(colnames(drugnew_data)%in%col)]
drugunique_df <- drugnew_df[!duplicated(drugnew_df$Gene.name),]
rownames(drugunique_df) <- drugunique_df[,"Gene.name"]
drugunique_df$Gene.name <- NULL
drugunique_df[1:3,1:3]


drug_gene1<-drugunique_df[a,]
drug_gene_mx<-as.matrix(drug_gene1)
rownames(drug_gene_mx)
class(drug_gene_mx)
class(data2)

dim(data2)  # 检查data2的维度
dim(stadcancerwide_df)  # 检查stadcancerwide_df的维度
dim(drug_gene_mx)  # 检查drug_gene_mx的维度
sum(is.na(data2))  # 检查缺失值数量
sum(is.na(stadcancerwide_df))  # 检查缺失值数量
sum(is.na(drug_gene_mx))  # 检查缺失值数量

a

rownames(stadcancerwide_df)
getwd()
calcPhen<-calcPhenotype(trainingExprData = as.matrix(data2),
                        trainingPtype = as.matrix(stadcancerwide_df),
                        testExprData = drug_gene_mx, # 需要matrix
                        batchCorrect = 'eb',  # 校正参数，"eb" for array, "standardize" for rnaseq
                        #IC50是对数转换的，所以表达矩阵也用对数转换过的
                        powerTransformPhenotype = F,
                        minNumSamples = 5,
                        printOutput = T,
                        removeLowVaryingGenes = 0.2,
                        removeLowVaringGenesFrom = "homogenizeData"
)

grug_pre<-read.csv("calcPhenotype_Output/stadDrugPredictions2.csv", header = T)

load("D:/细胞状态2/论文/fig2/totalESTclin2.RData")
group<-totalESTclin2
head(group)
stadgroup <- totalESTclin2[totalESTclin2$cancertype=="STAD",]
table(stadgroup$group1)
stadgroup$group1 <- ifelse(stadgroup$model1> median(stadgroup$model1),"high" ,"low")


#分组
stadgroup<-stadgroup[,c("bcr_sample_barcode","group1")]
table(stadgroup$group1)

head(stadgroup)
stadgroup$bcr_sample_barcode <- gsub("[-]", ".",stadgroup$bcr_sample_barcode)



stadDrugPred <- grug_pre
rownames(stadDrugPred)<-stadDrugPred$X
stadDrugPred<-stadDrugPred[,-1]


# 计算每列中NA值的数量
na_count <- colSums(is.na(stadDrugPred))

# 找到不含NA值的列的索引
clean_columns <- which(na_count == 0)

# 创建一个不含NA值的数据框
clean_df <- stadDrugPred[, clean_columns]
# 删除数据框中值大于8的列
clean_df <- clean_df[, !apply(clean_df, 2, function(x) any(x > 10))]


data1_df<-data.frame(t(clean_df))
#--------#
a<-intersect(colnames(data1_df),stadgroup$bcr_sample_barcode)
#------#
data1_df<-data1_df[,a]
drug_group<-stadgroup[stadgroup$bcr_sample_barcode%in% a,]


#差异分析
staddftest <- t(data1_df)
staddftest <- as.data.frame(staddftest)
staddftest$sample <- rownames(staddftest)
staddftest2 <- merge(staddftest,drug_group,by.x="sample",by.y="bcr_sample_barcode")
head(staddftest2)
dim(staddftest2)
staddftest4 <- staddftest2
staddftest4$sample <- NULL
staddftest4$group1 <- NULL

head(staddftest4)


stadnewFC <- data.frame(row.names = colnames(staddftest4))
stadnewFC$FC <- NA
stadnewFC$t <- NA
stadnewFC$wilcox <- NA
stadnewFC$aov <- NA
stadnewFC$highmean <- NA


for (i in colnames(staddftest4)) {
  print(i)
  yao <- staddftest2[,c(i,"group1")]
  head(yao)
  yao$value <- exp(yao[,1])
  mean_high <- mean(yao[yao$group1 == "high","value"])
  mean_low <- mean(yao[yao$group1 == "low","value"])
  stadnewFC[i,"highmean"] <- mean_high
  FC <- mean_high / mean_low
  stadnewFC[i,"FC"] <- FC
  t_test_result <- t.test(value ~ group1, data = yao)
  wilcox_test_result <- wilcox.test(value ~ group1, data = yao)
  aov_result <- aov(value ~ group1, data = yao)
  stadnewFC[i,"t"] <- t_test_result$p.value
  stadnewFC[i,"wilcox"] <- wilcox_test_result$p.value
  stadnewFC[i,"aov"] <- summary(aov_result)[[1]]$"Pr(>F)"[1]
  getwd()
  myfilename <- paste("./药物预测/boxplot/stad2/",i,".pdf",sep = "")
  pdf(myfilename)
  p <- ggplot(data = yao, aes(x =group1 , y =value , fill = group1)) +
    geom_violin(position = position_dodge(0.8))+
    geom_boxplot(width=0.5,
                 outlier.colour = NA,
                 position = position_dodge(0.8))+  #每组的小提琴图之间的距离
    stat_compare_means(aes(group = group1))+  #,label = 'p.signif'
    scale_fill_manual(values = c("#c11e21","steelblue") )+
    theme_bw()+theme(panel.grid=element_blank())+
    labs( x="drug", y="IC50 score")+
    theme(axis.text.x = element_text(size = 6,angle = 0, hjust = 0.5, vjust = 0.5))
  print(p)
  dev.off()
}

stadnewFC$LN_IC50 <- log(stadnewFC$highmean)
head(stadnewFC)
#stadnewFC[stadnewFC$FC>2,]
stadnewFC[stadnewFC$FC<0.5,]

dim(stadnewFC[stadnewFC$t<0.05,])
dim(stadnewFC[stadnewFC$wilcox<0.05,])
dim(stadnewFC[stadnewFC$aov<0.05,])


stadnewFC <- stadnewFC[stadnewFC$t<0.05&stadnewFC$wilcox<0.05&stadnewFC$aov<0.05,]
stadnewFC <- stadnewFC[stadnewFC$FC<1,]
stadnewFC[order(stadnewFC$highmean), ]
write.csv(stadnewFC[order(stadnewFC$highmean), ],file = "./药物预测/stad2.csv")

rownames(stadnewFC[order(stadnewFC$highmean), ])





