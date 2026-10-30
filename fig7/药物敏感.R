#install.packages("oncoPredict")
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
setwd("E:/T细胞障碍/肺癌/泛癌/药物敏感性")
#基因和细胞系关系
drug_id<-read.csv("E:/T细胞障碍/肺癌/泛癌/药物敏感性/细胞系编号.csv", header = T)
data <- read.table("E:/T细胞障碍/肺癌/泛癌/药物敏感性/Cell_line_RMA_proc_basalExp.txt", header=T, sep="\t")#读入表达数据

data[1:3,1:3]

a<-c("CD8A","CD8B","CD4","FOSB","JUNB","IKZF1","RUNX3","DDX5","MXD1","HSP90AB1",
     "DCAF12","ARID5A","SEMA4D","VAV1")
#细胞系与基因
data1<-data[!data$GENE_SYMBOLS=="",]
rownames(data1)<-data1$GENE_SYMBOLS
data1<-data1[,-c(1,2)]
names(data1)<- sub("^.*?\\.", "",names(data1))
names(data1)<- drug_id$Line[match(names(data1),drug_id$COSMIC_ID)]

#data2<-data1[a,]
#data2<-na.omit(data2)
#data3<-exp(data1)
#------筛选免疫治疗药和靶向药---细胞系与药物得关系-------#
drug_cell<-read.csv("E:/T细胞障碍/肺癌/泛癌/药物敏感性/GDSC2_fitted_dose_response_27Oct23.csv", header = T)
drug_cell<-drug_cell[,c(5,7,9,16)]
drug_cell_1<-drug_cell[!drug_cell$TCGA_DESC=="UNCLASSIFIED",]
cnacer_type<-unique(drug_cell_1$TCGA_DESC)

#drug_name<-c("Lenalidomide","Palbociclib","Dasatinib","Erlotinib","Tamoxifen")
#subset_df <- drug_cell[drug_cell$DRUG_NAME %in% drug_name,]
#---选择癌症-----#
cancer<-drug_cell_1[drug_cell_1$TCGA_DESC=="UCEC",]
cancer<-cancer[,-2]
wide_df <-dcast(cancer,CELL_LINE_NAME~DRUG_NAME , value.var = "LN_IC50",fun.aggregate = mean)
wide_df[is.na(wide_df)] <- 0
rownames(wide_df)<-wide_df$CELL_LINE_NAME
wide_df<-wide_df[,-1]
#wide_df1<- exp(wide_df)
#---------药物必须要一致-----------------#
b<-intersect(rownames(wide_df),colnames(data1))
wide_df<-wide_df[b,]
data2<-data1[,b]
identical(rownames(wide_df),colnames(data2))
#-----#
setwd("F:/fsh/bulk/BULK_表达")
file_name <- list.files()
cancer_strings <- file_name[grepl("UCEC", file_name)]
#------#
drug_gene<- read.table(cancer_strings, header=T, sep="\t")#读入表达数据
#drug_gene<-log2(drug_gene)
drug_gene1<-drug_gene[a,]
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
grug_pre<-read.csv("E:/T细胞障碍/肺癌/泛癌/calcPhenotype_Output/DrugPredictions.csv", header = T)

group<-read.csv("E:/T细胞障碍/肺癌/泛癌/bulk/bulk_GSVA/bulk_GSVA_2_model_all.csv", header = T)
#分组
group<-group[,c(2,12)]
table(group$group)

#df_grouped1<-df_grouped[!df_grouped$sample %in% c("CHOL","OV","LUSC","LGG","UCS"),]

#table(df_grouped1$group)


#------------#
#brca_mcp_sub2<-unique(brca_mcp_sub2)
brca_mcp_sub2 <- reshape2::melt(grug_pre,id=c("X"))
grug_sample<-merge(group,brca_mcp_sub2,by.x="simple",by.y="X")

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

