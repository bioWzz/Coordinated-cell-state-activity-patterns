setwd("F:/fsh/bulk/药物预测")
first_file_name<- list.files() 
dir<-"F:/fsh/bulk/重做不同癌型药物的差异"
simple_group<-read.csv("E:/T细胞障碍/肺癌/泛癌/bulk/bulk_GSVA/bulk_gsva_2model.csv", header = T)
simple_group<-simple_group[,c(2,12)]
for (i in 1:19) {
  data<-read.csv(first_file_name[i], header =T)
  rownames(data)<-data$X
  data<-data[,-1]

  # 计算每列中NA值的数量
  na_count <- colSums(is.na(data))

  # 找到不含NA值的列的索引
  clean_columns <- which(na_count == 0)

  # 创建一个不含NA值的数据框
  clean_df <- data[, clean_columns]
  # 删除数据框中值大于8的列
  clean_df <- clean_df[, !apply(clean_df, 2, function(x) any(x > 8))]


  data1_df<-data.frame(t(clean_df))
  #--------#
  a<-intersect(colnames(data1_df),simple_group$X)
  #------#
  data1_df<-data1_df[,a]
  drug_group<-simple_group[simple_group$X%in% a,]
  #-----差异分析-----#
  design<-model.matrix(~drug_group$group)#只要分组信息
  fit=lmFit(data1_df,design)
  fit=eBayes(fit)
  deg=topTable(fit,coef=2,number = Inf)  #提取所有基因的差异分析结果
  print(colnames(deg))
  ##设置1ogFC值和P值来标记上下调基因
  # 将logFC转换为原始倍数变化
  #deg$FC <- exp(deg$logFC)
  logFC=1#设置1ogFC为1
  P.Value = 0.05#设置P.Va1ue为0.05
  type1 = (deg$P.Value < P.Value)&(deg$logFC < -logFC)#提取P.Va1ue小于0.05，1ogFC小于1的
  type2 = (deg$P.Value < P.Value)&(deg$logFC > logFC)
  deg$Group = ifelse(type1,"Down",ifelse(type2,"Up","Not-Sig"))#判断基因表达的上下调
  table(deg$Group)##统计上下调的基因数量
  deg$cancer<-sub("_.*", "",first_file_name[i])
  
  dir_name<-paste(dir,first_file_name[i],sep = "/")
  write.csv(deg,dir_name)#导出结果
}
#----------------#
setwd("F:/fsh/bulk/重做不同癌型药物的差异")
first_file_name<- list.files() 
dir<-"F:/fsh/bulk/药物预测"
simple_group<-read.csv("E:/T细胞障碍/肺癌/泛癌/bulk/bulk_GSVA/bulk_gsva_2model.csv", header = T)
simple_group<-simple_group[,c(2,12)]

for (i in 5:length(first_file_name)) {
  #差异药物
  data<-read.csv(first_file_name[19], header =T)
  rownames(data)<-data$X
  data<-data[,-1]
  UP<-data[data$Group=="Down",]
  drug_name<-rownames(UP)
  #预测药物
  dir_name<-paste(dir,first_file_name[19],sep = "/")
  pre_data<-read.csv(dir_name, header =T)
  rownames(pre_data)<-pre_data$X
  pre_data<-pre_data[,-1]
  
  subset_df <- pre_data[,drug_name]
  
  subset_df$simple<-rownames(subset_df)
  drug_group<-merge(subset_df,simple_group,by.x="simple",by.y="X")
  rownames(drug_group)<-drug_group$simple
  drug_group<-drug_group[,-1]
  
  brca_mcp_sub2 <- reshape2::melt(drug_group,id=c("group"))
  brca_mcp_sub2$group <- as.factor(brca_mcp_sub2$group)
  group <- levels(factor(brca_mcp_sub2$group ))
  brca_mcp_sub2$group  <- factor(brca_mcp_sub2$group,levels = group)
  comp <- combn(group,2)
  my_comparisons <- list()
  for (j in 1:ncol(comp)) {
    my_comparisons[[j]] <- comp[,j]
  }
  
  cancer_name <- sub("_.*", "", first_file_name[i])
  name<-paste(cancer_name,".pdf",sep = "") 
  #name_dir<-paste(dir,name,sep = "/") 
  
  pdf(file=name, height=15,width=20)
  
  #--111--#
  p<-ggplot(data = brca_mcp_sub2, aes(x =variable , y =value , fill = group)) +
    geom_violin(position = position_dodge(0.8))+
    geom_boxplot(width=0.5,
                 outlier.colour = NA,
                 position = position_dodge(0.8))+  #每组的小提琴图之间的距离
    stat_compare_means(aes(group = group),label = 'p.signif')+  #,label = 'p.signif'
    scale_fill_manual(values = c("#7493C2","#F7D1AC") )+
    theme_bw()+theme(panel.grid=element_blank())+
    labs( x="drug", y="IC50 score")+
    #ylim(0,1.15)+
    #theme_prism()+
    #scale_y_log10()+
    theme(axis.text.x = element_text(size = 6,angle = 0, hjust = 0.5, vjust = 0.5))
  print(p)
  dev.off()
}
#------------------#
setwd("F:/fsh/bulk/敏感癌型")
first_file_name<- list.files() 
for (i in 1:length(first_file_name)) {
  data<-read.csv(first_file_name[i], header =T)
  rownames(data)<-data$X
  data<-data[,-1]
  UP<-data[data$Group=="Down",]
  drug_name<-rownames(UP)
  print(first_file_name[i])
  print(drug_name)
  print("-----------")
}


#-------Melanoma---------#
library(dplyr)
setwd("E:/T细胞障碍/肺癌/免疫治疗/GSVA")
first_file_name<- list.files() 

all<-NA
for (i in 1:length(first_file_name)) {
  data<-read.csv(first_file_name[i], header = T)
  rownames(data)<-data$X
  data<-data[,-1]
  
  all<-cbind(all,data)
}
all<-all[,-1]
GSVA<-data.frame(t(all))
GSVA$simple<-rownames(GSVA)
GSVA$simple <- gsub("X", "", GSVA$simple)
#------------#
setwd("E:/T细胞障碍/肺癌/免疫治疗/临床数据")
first_file_name<- list.files() 
all<-NA
for (i in 1:length(first_file_name)) {
  data<-data <- read.table(first_file_name[i], header=T, sep="\t",fill=TRUE)
  all<-rbind(all,data)
}
all<-all[-1,]
cin<-all
cin<-cin[,c(1,9,10,13)]
names(cin)<-c("simple_id","OS_time","OS","Therapy")
table(cin$Therapy)

anti_PD_cin<-cin[cin$Therapy=="anti-CTLA-4",]
anti_PD_cin1<-cin[cin$Therapy=="anti-CTLA-4+anti-PD-1",]
#---------------#
GSVA_fenzu <- GSVA %>%
  mutate(exclude_category = case_when(
    exclude > 0 ~ "高",
    exclude < 0 ~ "低"
  ),
  exhaustion_category = case_when(
    exhaustion > 0 ~ "高",
    exhaustion <0 ~ "低"
  ),
  senility_category = case_when(
    senility > 0 ~ "高",
    senility <0 ~ "低"
  ),
  Tcell.dysfunction_category = case_when(
    Tcell.dysfunction > 0 ~ "高",
    Tcell.dysfunction < 0 ~ "低"
  ))
#分两组
library("survival")
library("survminer")
df_grouped <- GSVA_fenzu %>%
  mutate(group = ifelse(exclude_category == "高" & exhaustion_category == "低" & senility_category == "低" & Tcell.dysfunction_category == "低", "Moderate", "Disorderly")) 
#
table(df_grouped$group)
#---#
M_cin<-merge(df_grouped,anti_PD_cin,by.x = "simple",by.y = "simple_id")
#-------------#
table(M_cin$Therapy)
#-----------#

M_cin<-M_cin[,c(1,10,11,12)]
M_cin$OS <- ifelse(M_cin$OS  == "Alive", 0, ifelse(M_cin$OS == "Dead", 1,M_cin$OS ))
M_cin <- na.omit(M_cin)
M_cin$OS<-as.numeric(M_cin$OS)
sfit <- survfit(Surv(OS_time, OS)~group, data=M_cin)
print(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)
#---------------------------#
anti_PD_cin2<-rbind(anti_PD_cin,anti_PD_cin1)


PD1<-merge(df_grouped,anti_PD_cin2,by.x = "simple",by.y = "simple_id")
PD1<-PD1[,c(1,2,10)]



PD1$group <- as.factor(PD1$group)
group <- levels(factor(PD1$group ))
PD1$group  <- factor(PD1$group,levels = group)
comp <- combn(group,2)
my_comparisons <- list()
for (j in 1:ncol(comp)) {
  my_comparisons[[j]] <- comp[,j]
}


#--111--#
p<-ggplot(data = PD1, aes(x =group , y =exclude , fill = group)) +
  #geom_violin(position = position_dodge(0.8))+
  geom_boxplot(width=0.5,
               outlier.colour = NA,
               position = position_dodge(0.8))+  #每组的小提琴图之间的距离
  stat_compare_means(aes(group = group),)+  #,label = 'p.signif'
  scale_fill_manual(values = c("#F9A4A4","#B5E2EA") )+
  theme_bw()+theme(panel.grid=element_blank())+
  labs( x="Obstacle types of T cell", y="exclude score")+
  #ylim(0,1.15)+
  #theme_prism()+
  #scale_y_log10()+
  theme(axis.text.x = element_text(size = 10,angle = 0, hjust = 0.5, vjust = 0.5))
#---------筛选药物的ADME-------------#
drug_data<-read.csv("E:/T细胞障碍/肺癌/泛癌/药物敏感性/SwissADME.csv", header = T)
hi_drug_data<-drug_data[drug_data$GI.absorption=="High",]


filtered_df <- hi_drug_data[(hi_drug_data$Lipinski == 0 & hi_drug_data$Ghose == 0 & hi_drug_data$Veber == 0 & hi_drug_data$Egan == 0) |
                    (hi_drug_data$Lipinski == 0 & hi_drug_data$Ghose == 0 & hi_drug_data$Veber == 0 & hi_drug_data$Muegge == 0) |
                    (hi_drug_data$Lipinski == 0 & hi_drug_data$Ghose == 0 & hi_drug_data$Egan == 0 & hi_drug_data$Muegge == 0) |
                    (hi_drug_data$Lipinski == 0 & hi_drug_data$Veber == 0 & hi_drug_data$Egan == 0 & hi_drug_data$Muegge == 0) |
                    (hi_drug_data$Ghose == 0 & hi_drug_data$Veber == 0 & hi_drug_data$Egan == 0 & hi_drug_data$Muegge == 0), ]

filtered_df <- hi_drug_data[(hi_drug_data$Lipinski == 0 & 
                               hi_drug_data$Ghose == 0 & 
                               hi_drug_data$Veber == 0 & 
                               hi_drug_data$Egan == 0&
                               hi_drug_data$Muegge == 0),]
write.csv(filtered_df,"E:/T细胞障碍/肺癌/泛癌/药物敏感性/筛选过后的药物.csv")#导出结果
#------------#
setwd("E:/T细胞障碍/肺癌/泛癌/药物敏感性/ADME2")
first_file_name<- list.files() 
all<-NA


for (i in 1:length(first_file_name)) {
  data<-read.csv(first_file_name[i], header = T)
  rownames(data)<-data$X
  all<-rbind(all,data)
}
drug_all<-all[-1,]
drug_all<-drug_all[,c(1,2,4,8,10,12,27,33)]


smiles<-read.csv("E:/T细胞障碍/肺癌/泛癌/药物敏感性/筛选过后的药物.csv", header = T)
smiles<-smiles[,c(1,3)]


sel_smiles<-merge(smiles,drug_all,by.x="Canonical.SMILES",by.y = "smiles")
rownames(sel_smiles)<-sel_smiles$X
sel_smiles<-sel_smiles[,-c(1,2)]
sel_smiles<-data.frame(t(sel_smiles))

write.csv(sel_smiles,"E:/T细胞障碍/肺癌/泛癌/药物敏感性/ADME_smiles.csv")#导出结果
#-----------#
library(pheatmap)
library(ggplot2)
pheatmap(sel_smiles, scale = "none", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
         
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",#表示行、列聚类使用的度量方法，默认为“euclidean”，也可为  "correlation"即按照 Pearson correlation方法进行聚类
         
         show_rownames = T,
         
         #cluster_rows = T,
         
         #main =cancer_name[i] ,
         
         cellwidth = 20, cellheight = 20,
         
         clustering_method = "complete", #表示聚类方法，包括：‘ward’, ‘ward.D’, ‘ward.D2’, ‘single’, ‘complete’, ‘average’, ‘mcquitty’, ‘median’, ‘centroid’
         
         cluster_rows = T,cluster_cols =F , #cluster_rows表示仅对行聚类，cluster_cols表示仅对列聚类，值为TRUE或FALSE
         
         cutree_rows = 2, #cutree_cols = , #若进行了行、列聚类，根据行、列聚类数量分隔热图行
         
         treeheight_row = 30, treeheight_col = 30,#、行、列聚类树高度调整
         
         number_format = "%.2f", #表示热图单元格上显示的数据格式，“%.2f” 表示两位小数,“%.1e”表示科学计数法
         
         number_color = "grey30", #表示热图单元格上显示的数据字体颜色
         display_numbers = T, #表示是否在单元格上显示原始数值或按照特殊条件进行区分标记
         
         fontsize_number = 6 #表示热图上显示数字的字体大小
         
         #color = colorRampPalette(c('#266EAF',"white",'#B52230'))(100), #表示热图颜色,(100)表示100个等级
         #color = c(colorRampPalette(colors = c("#266EAF","white"))(length(bk)/2),colorRampPalette(colors = c("white","#B52230"))(length(bk)/2)),
         #breaks=bk
) 
#--------#
#------------------#
setwd("F:/fsh/bulk/敏感癌型")
first_file_name<- list.files() 
downl<-NA
for (i in 1:length(first_file_name)) {
  data<-read.csv(first_file_name[i], header =T)
  rownames(data)<-data$X
  data<-data[,-1]
  UP<-data[data$Group=="Down",]
  downl<-rbind(downl,UP)
}
downl<-downl[-1,]

a<-intersect(rownames(downl),colnames(sel_smiles))
drug_downl<-downl[a,]
drug_downl$drug<-rownames(drug_downl)

b<-c( "Ulixertinib","Gefitinib","AMG.319","OSI.027", "Paclitaxel","Taselisib","Tanespimycin")
drug_downl_1<-drug_downl[b,]

#-----------#
ggplot(drug_downl_1, aes(drug,AveExpr))+
  geom_col(fill="#94C1E0")+ coord_flip()+
  labs( x="drug", y="IC50 Score")+
  theme_bw()+theme(panel.grid=element_blank())  

