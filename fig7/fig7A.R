install.packages("oncoPredict")
BiocManager::install("GenomicFeatures")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

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

setwd("D:/细胞状态2/论文/药物")

drug_id<-read.csv("细胞系编号.csv", header = T)
data <- read.table("Cell_line_RMA_proc_basalExp.txt", header=T, sep="\t")#读入表达数据

data[1:3,1:3]

#PAAD
a<-c("ACTA2" , "AEBP1" , "ANXA2","ATF3","BGN","CALD1","COL1A1","CTSK",
     "EGR1","EGR3","FN1","HTRA1","INHBA","MMP11","MMP14","MMP2","NR4A1","PDGFA",
     "S100A4","SPP1","VCAN")

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

#---选择癌症-----#
cancer<-drug_cell_1[drug_cell_1$TCGA_DESC=="PAAD",]
cancer<-cancer[,-2]
wide_df <-dcast(cancer,CELL_LINE_NAME~DRUG_NAME , value.var = "LN_IC50",fun.aggregate = mean)
#wide_df <-reshape2::dcast(cancer,CELL_LINE_NAME~DRUG_NAME , value.var = "LN_IC50",fun.aggregate = mean)
wide_df[is.na(wide_df)] <- 0
rownames(wide_df)<-wide_df$CELL_LINE_NAME
wide_df<-wide_df[,-1]

#---------药物必须要一致-----------------#
b<-intersect(rownames(wide_df),colnames(data1))
wide_df<-wide_df[b,]
data2<-data1[,b]
identical(rownames(wide_df),colnames(data2))

genedata <- read.table("D:/细胞状态2/论文/分化/fenhua/cox/mart_export.txt", sep=",", header=TRUE)

#setwd("F:/fsh/bulk/BULK_表达")
#file_name <- list.files("D:/细胞状态2/bulk/expression-counts")
#cancer_strings <- file_name[grepl("PAAD", file_name)]
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

#-----#
calcPhen<-calcPhenotype(trainingExprData = as.matrix(data2),
                        trainingPtype = as.matrix(wide_df),
                        testExprData = drug_gene_mx, # 需要matrix
                        batchCorrect = 'eb',  # 校正参数，"eb" for array, "standardize" for rnaseq
                        #IC50是对数转换的，所以表达矩阵也用对数转换过的
                        powerTransformPhenotype = F,
                        minNumSamples = 10,
                        printOutput = T,
                        removeLowVaryingGenes = 0.2,
                        removeLowVaringGenesFrom = "homogenizeData"
)
#-----------------#

grug_pre<-read.csv("calcPhenotype_Output/DrugPredictions.csv", header = T)
max(grug_pre)

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


#------------#
#brca_mcp_sub2<-unique(brca_mcp_sub2)
brca_mcp_sub2 <- reshape2::melt(grug_pre,id=c("X"))
grug_sample<-merge(paadgroup,brca_mcp_sub2,by.x="bcr_sample_barcode",by.y="X")
head(grug_sample)
table(grug_sample$variable)










#STAD
a<-c("SELE","CXCL12","EGR1","PDPN","CD44","COL6A2","LGALS1","FBLN5")

# data1<-data[!data$GENE_SYMBOLS=="",]
# rownames(data1)<-data1$GENE_SYMBOLS
# data1<-data1[,-c(1,2)]
# names(data1)<- sub("^.*?\\.", "",names(data1))
# names(data1)<- drug_id$Line[match(names(data1),drug_id$COSMIC_ID)]
# data1[1:3,1:3]
# 
# 
# #------筛选免疫治疗药和靶向药---细胞系与药物得关系-------#
# drug_cell<-read.csv("GDSC2_fitted_dose_response_27Oct23.csv", header = T)
# drug_cell<-drug_cell[,c(5,7,9,16)]
# drug_cell_1<-drug_cell[!drug_cell$TCGA_DESC=="UNCLASSIFIED",]
# cnacer_type<-unique(drug_cell_1$TCGA_DESC)
#---选择癌症-----#
cancer<-drug_cell_1[drug_cell_1$TCGA_DESC=="STAD",]
cancer<-cancer[,-2]
wide_df <-dcast(cancer,CELL_LINE_NAME~DRUG_NAME , value.var = "LN_IC50",fun.aggregate = mean)
#wide_df <-reshape2::dcast(cancer,CELL_LINE_NAME~DRUG_NAME , value.var = "LN_IC50",fun.aggregate = mean)
wide_df[is.na(wide_df)] <- 0
rownames(wide_df)<-wide_df$CELL_LINE_NAME
wide_df<-wide_df[,-1]

#---------药物必须要一致-----------------#
b<-intersect(rownames(wide_df),colnames(data1))
wide_df<-wide_df[b,]
data2<-data1[,b]
identical(rownames(wide_df),colnames(data2))

genedata <- read.table("D:/细胞状态2/论文/分化/fenhua/cox/mart_export.txt", sep=",", header=TRUE)

#setwd("F:/fsh/bulk/BULK_表达")
#file_name <- list.files("D:/细胞状态2/bulk/expression-counts")
#cancer_strings <- file_name[grepl("PAAD", file_name)]
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

#-----#
calcPhen<-calcPhenotype(trainingExprData = as.matrix(data2),
                        trainingPtype = as.matrix(wide_df),
                        testExprData = drug_gene_mx, # 需要matrix
                        batchCorrect = 'eb',  # 校正参数，"eb" for array, "standardize" for rnaseq
                        #IC50是对数转换的，所以表达矩阵也用对数转换过的
                        powerTransformPhenotype = F,
                        minNumSamples = 5,
                        printOutput = T,
                        removeLowVaryingGenes = 0.2,
                        removeLowVaringGenesFrom = "homogenizeData"
)
#-----------------#

grug_pre<-read.csv("calcPhenotype_Output/DrugPredictions.csv", header = T)

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


#------------#
#brca_mcp_sub2<-unique(brca_mcp_sub2)
brca_mcp_sub2 <- reshape2::melt(grug_pre,id=c("X"))
grug_sample<-merge(stadgroup,brca_mcp_sub2,by.x="bcr_sample_barcode",by.y="X")
head(grug_sample)
table(grug_sample$variable)














#候选药
grug_name<-grug_sample[grug_sample$variable=="Palbociclib",]

grug_name$group <- as.factor(grug_name$group )
group <- levels(factor(grug_name$group ))
grug_name$group  <- factor(grug_name$group,levels = group)
comp <- combn(group,2)
my_comparisons <- list()
for (j in 1:ncol(comp)) {
  my_comparisons[[j]] <- comp[,j]
}
ggplot(data = grug_name, aes(x = group, y = value, fill = group)) +
  geom_violin(position = position_dodge(0.8))+
  geom_boxplot(width=0.5,
               outlier.colour = NA,
               position = position_dodge(0.8))+  #每组的小提琴图之间的距离
  stat_compare_means(aes(group = group),)+  #,label = 'p.signif'
  scale_fill_manual(values = c("#7493C2","#F7D1AC") )+
  theme_bw()+theme(panel.grid=element_blank())+
  labs( x="Obstacle types of T cell", y="IC50 score")+
  #ylim(0,1.15)+
  #theme_prism()+
  #scale_y_log10()+
  theme(axis.text.x = element_text(size = 20,angle = 0, hjust = 0.5, vjust = 0.5),legend.position = "none")


