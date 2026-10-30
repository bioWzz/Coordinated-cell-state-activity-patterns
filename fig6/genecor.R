library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(clustree)
library(patchwork)
library(corrplot)

#PAAD
load("D:/细胞状态2/论文/空间/st2/seurat_data.RData")
load("D:/细胞状态2/论文/空间/st/EST.RData")
load("/wangzhenzenlab/fsh/xff/st2/seurat_data.RData")
load("/wangzhenzenlab/fsh/xff/st2/EST.RData")


#合并表达数据与model1的值
head(seurat_data)
EST <- as.data.frame(EST)
EST[1:3,]
EST$cell <- rownames(EST)
expression_matrix <- GetAssayData(seurat_data, assay = "Spatial", slot = "data")
dim(expression_matrix)
expression_matrix[1:3,1:3]
expression_matrix2 <- as.data.frame(t(expression_matrix))
expression_matrix2[1:3,1:3]
expression_matrix2$cell <- rownames(expression_matrix2)


expression_matrix3 <- merge(expression_matrix2,EST,by="cell")
expression_matrix3$model1 <- rowMeans(expression_matrix3[,c("Angiogenesis","EMT","Invasion","Metastasis")])
expression_matrix3[1:3,1:3]
expression_matrix3$cell <- NULL
singlemt <- expression_matrix3
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
as.data.frame(cor_pearson[,"model1"])
write.csv(as.data.frame(cor_pearson[,c("model1","Angiogenesis","EMT","Invasion","Metastasis")]),file="D:/细胞状态2/论文/空间/paad_gene_cor.csv")
my_colors <- colorRampPalette(c("steelblue", "white", "#c11e21"))(100)
corrplot(cor_pearson[c("FN1","COL1A1","COL1A2","COL3A1","COL6A3","SPARC","LGALS1","BGN","model1","Angiogenesis","EMT","Invasion","Metastasis"),
                     c("FN1","COL1A1","COL1A2","COL3A1","COL6A3","SPARC","LGALS1","BGN","model1","Angiogenesis","EMT","Invasion","Metastasis")], 
         method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)

cor_pearson[c("FN1","COL1A1","COL1A2","COL3A1","COL6A3","SPARC","LGALS1","BGN","model1","Angiogenesis","EMT","Invasion","Metastasis"),
            c("FN1","COL1A1","COL1A2","COL3A1","COL6A3","SPARC","LGALS1","BGN","model1","Angiogenesis","EMT","Invasion","Metastasis")]




#STAD
load("D:/细胞状态2/论文/空间/stad/stadseurat_data.RData")
load("D:/细胞状态2/论文/空间/stad/stadEST.RData")
load("/wangzhenzenlab/fsh/xff/stad/stadseurat_data.RData")
load("/wangzhenzenlab/fsh/xff/stad/stadEST.RData")
#合并表达数据与model1的值
head(stadseurat_data)
stadEST <- as.data.frame(stadEST)
stadEST[1:3,]
stadEST$cell <- rownames(stadEST)
stadexpression_matrix <- GetAssayData(stadseurat_data, assay = "Spatial", slot = "data")
dim(stadexpression_matrix)
stadexpression_matrix[1:3,1:3]
stadexpression_matrix2 <- as.data.frame(t(stadexpression_matrix))
stadexpression_matrix2[1:3,1:3]
stadexpression_matrix2$cell <- rownames(stadexpression_matrix2)


stadexpression_matrix3 <- merge(stadexpression_matrix2,stadEST,by="cell")
stadexpression_matrix3$model1 <- rowMeans(stadexpression_matrix3[,c("Angiogenesis","EMT","Invasion","Metastasis")])
stadexpression_matrix3[1:3,1:3]
stadexpression_matrix3$cell <- NULL
singlemt <- stadexpression_matrix3
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson


write.csv(as.data.frame(cor_pearson[c("FN1","COL1A1","COL1A2","COL3A1","COL6A3","SPARC","LGALS1","BGN"),
                                    c("model1","Angiogenesis","EMT","Invasion","Metastasis")]),
          file="D:/细胞状态2/论文/空间/stad_gene_cor.csv")
write.csv(cor_pearson,file="/wangzhenzenlab/fsh/xff/stad/stad_gene_cor.csv")

cor_pearson[,c("model1","Angiogenesis","EMT","Invasion","Metastasis")]
my_colors <- colorRampPalette(c("steelblue", "white", "#c11e21"))(100)
corrplot(cor_pearson[c("FN1","COL1A1","COL1A2","COL3A1","COL6A3","SPARC","LGALS1","BGN","model1","Angiogenesis","EMT","Invasion","Metastasis"),
                     c("FN1","COL1A1","COL1A2","COL3A1","COL6A3","SPARC","LGALS1","BGN","model1","Angiogenesis","EMT","Invasion","Metastasis")], 
         method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)








paadtop10 <- expression_matrix2[,c("FN1","COL1A1","COL1A2","COL3A1","COL6A3","TIMP3","POSTN","VIM","SPARC","LGALS1")]
head(paadtop10)
write.csv(paadtop10,file="D:/细胞状态2/论文/空间/paadtop10.csv")

stadtop10 <- stadexpression_matrix2[,c("SPARC","COL1A2","COL1A1","IGFBP7","COL3A1","COL4A1","COL6A3","LUM","COL4A2","TGFBI")]
head(stadtop10)
write.csv(stadtop10,file="D:/细胞状态2/论文/空间/stadtop10.csv")


paadtop <- expression_matrix2[,c("FN1","COL1A1","COL1A2","COL3A1","COL6A3",
                                 "TIMP3","POSTN","VIM","SPARC","LGALS1",
                                 "COL11A1","GREM1","BGN","SULF1","INHBA",
                                 "MMP14","THBS2")]
head(paadtop)
write.csv(paadtop,file="D:/细胞状态2/论文/空间/paadtop.csv")

stadtop <- stadexpression_matrix2[,c("SPARC","COL1A2","COL1A1","IGFBP7","COL3A1",
                                     "COL4A1","COL6A3","LUM","COL4A2","TGFBI",
                                     "COL12A1","COL5A2","FTL","LGALS1","FSTL1",
                                     "PMEPA1","TGM2","NNMT","IFITM3","FN1",
                                     "A2M","TIMP1","THY1","BGN","FLNA",
                                     "PLAU","PDGFRB","MMP2","CTSB","SERPINH1",
                                     "FCER1G","CEBPB","DCN","MGP")]
head(stadtop)
write.csv(stadtop,file="D:/细胞状态2/论文/空间/stadtop.csv")





#单cox
setwd("D:/细胞状态2/论文/分化/fenhua/cox")
library(stringr)

library(survival)
library(survminer)

#Ensembl_ID转基因数据
genedata <- read.table("mart_export.txt", sep=",", header=TRUE)
head(genedata)

#临床数据
clin_all <- read.table("D:/细胞状态2/bulk/clinical/cli_all.txt",header = T)
head(clin_all)

#PAAD
paadbulk <- read.table("D:/细胞状态2/bulk/expression-counts/TCGA-PAAD.htseq_counts.tsv", sep="\t",header=TRUE)
paadbulk[1:3,1:3]

paadbulk$Ensembl_ID <- str_split_fixed(paadbulk$Ensembl_ID,'[.]',2)[,1]
new_data <- merge(paadbulk,genedata,by.x = "Ensembl_ID",by.y = "Gene.stable.ID")
col<- c(1, 184, 185, 186) 
new_df <- new_data[, -col]
unique_df <- new_df[!duplicated(new_df$Gene.name),]
rownames(unique_df) <- unique_df[,183]
unique_df <- unique_df[, -183]

unique_df[1:3,1:3]

# 进行log转化，避免零值影响
log_counts <- log2(unique_df + 1)
log_counts[1:3,1:3]
min(log_counts)


log_countsT <- t(log_counts)
log_countsT[1:3,1:3]
log_countsT <- as.data.frame(log_countsT)
log_countsT$sample <- rownames(log_countsT)
log_countsT$patient <- sub("\\.[^\\.]*$", "", log_countsT$sample)
head(log_countsT$patient)


clin_paad <- merge(log_countsT,clin_all,by.x="patient",by.y="PATIENT")
clin_paad[1:3,1:3]


#("FN1","COL1A1","COL1A2","COL3A1","COL6A3","TIMP3","POSTN","VIM","SPARC","LGALS1")
cox_model <- coxph(Surv(OS, EVENT) ~ FN1+COL1A1+COL1A2+COL3A1+COL6A3+TIMP3+POSTN+VIM+SPARC+LGALS1, data = clin_paad)
cox_model <- coxph(Surv(OS, EVENT) ~ FN1, data = clin_paad)
cox_model <- coxph(Surv(OS, EVENT) ~ COL1A1, data = clin_paad)
cox_model <- coxph(Surv(OS, EVENT) ~ COL1A2, data = clin_paad)
cox_model <- coxph(Surv(OS, EVENT) ~ COL3A1, data = clin_paad)
cox_model <- coxph(Surv(OS, EVENT) ~ COL6A3, data = clin_paad)
cox_model <- coxph(Surv(OS, EVENT) ~ TIMP3, data = clin_paad)
cox_model <- coxph(Surv(OS, EVENT) ~ POSTN, data = clin_paad)
cox_model <- coxph(Surv(OS, EVENT) ~ VIM, data = clin_paad)
cox_model <- coxph(Surv(OS, EVENT) ~ SPARC, data = clin_paad)
cox_model <- coxph(Surv(OS, EVENT) ~ LGALS1, data = clin_paad)
summary(cox_model)
# 绘制森林图
ggforest(cox_model, data = clin_paad)












#STAD
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

#c("SPARC","COL1A2","COL1A1","IGFBP7","COL3A1","COL4A1","COL6A3","LUM","COL4A2","TGFBI")

stadcox_model1 <- coxph(Surv(OS, EVENT) ~ SPARC+COL1A2+COL1A1+IGFBP7+COL3A1+COL4A1+COL6A3+LUM+COL4A2+TGFBI, data = clin_stad)
stadcox_model1 <- coxph(Surv(OS, EVENT) ~ SPARC, data = clin_stad)
stadcox_model1 <- coxph(Surv(OS, EVENT) ~ COL1A2, data = clin_stad)
stadcox_model1 <- coxph(Surv(OS, EVENT) ~ COL1A1, data = clin_stad)
stadcox_model1 <- coxph(Surv(OS, EVENT) ~ IGFBP7, data = clin_stad)
stadcox_model1 <- coxph(Surv(OS, EVENT) ~ COL3A1, data = clin_stad)
stadcox_model1 <- coxph(Surv(OS, EVENT) ~ COL4A1, data = clin_stad)
stadcox_model1 <- coxph(Surv(OS, EVENT) ~ COL6A3, data = clin_stad)
stadcox_model1 <- coxph(Surv(OS, EVENT) ~ LUM, data = clin_stad)
stadcox_model1 <- coxph(Surv(OS, EVENT) ~ COL4A2, data = clin_stad)
stadcox_model1 <- coxph(Surv(OS, EVENT) ~ TGFBI, data = clin_stad)
summary(stadcox_model1)

# 绘制森林图
ggforest(stadcox_model1, data = clin_stad)









#独立验证集的高度相关基因
load("/wangzhenzenlab/fsh/xff/stdata/ES.RData")
load("/wangzhenzenlab/fsh/xff/stdata/stadexp.RData")
expression_matrix0331 <- as.data.frame(t(expression_matrix0331))
expression_matrix0731 <- as.data.frame(t(expression_matrix0731))
expression_matrix0733 <- as.data.frame(t(expression_matrix0733))
expression_matrix0734 <- as.data.frame(t(expression_matrix0734))
expression_matrix1251 <- as.data.frame(t(expression_matrix1251))
expression_matrix1252 <- as.data.frame(t(expression_matrix1252))
expression_matrix1253 <- as.data.frame(t(expression_matrix1253))
expression_matrix1254 <- as.data.frame(t(expression_matrix1254))
expression_matrix1675 <- as.data.frame(t(expression_matrix1675))
expression_matrix0331[1:3,1:3]
expression_matrix0331$cell <- rownames(expression_matrix0331)
expression_matrix0731$cell <- rownames(expression_matrix0731)
expression_matrix0733$cell <- rownames(expression_matrix0733)
expression_matrix0734$cell <- rownames(expression_matrix0734)
expression_matrix1251$cell <- rownames(expression_matrix1251)
expression_matrix1252$cell <- rownames(expression_matrix1252)
expression_matrix1253$cell <- rownames(expression_matrix1253)
expression_matrix1254$cell <- rownames(expression_matrix1254)
expression_matrix1675$cell <- rownames(expression_matrix1675)

EST0331 <- as.data.frame(t(ES0331))
EST0731 <- as.data.frame(t(ES0731))
EST0733 <- as.data.frame(t(ES0733))
EST0734 <- as.data.frame(t(ES0734))
EST1251 <- as.data.frame(t(ES1251))
EST1252 <- as.data.frame(t(ES1252))
EST1253 <- as.data.frame(t(ES1253))
EST1254 <- as.data.frame(t(ES1254))
EST1675 <- as.data.frame(t(ES1675))

EST0331$cell <- rownames(EST0331)
EST0731$cell <- rownames(EST0731)
EST0733$cell <- rownames(EST0733)
EST0734$cell <- rownames(EST0734)
EST1251$cell <- rownames(EST1251)
EST1252$cell <- rownames(EST1252)
EST1253$cell <- rownames(EST1253)
EST1254$cell <- rownames(EST1254)
EST1675$cell <- rownames(EST1675)
EST0331[1:3,]

Mexpression_matrix0331 <- merge(expression_matrix0331, EST0331, by = "cell")
Mexpression_matrix0731 <- merge(expression_matrix0731, EST0731, by = "cell")
Mexpression_matrix0733 <- merge(expression_matrix0733, EST0733, by = "cell")
Mexpression_matrix0734 <- merge(expression_matrix0734, EST0734, by = "cell")
Mexpression_matrix1251 <- merge(expression_matrix1251, EST1251, by = "cell")
Mexpression_matrix1252 <- merge(expression_matrix1252, EST1252, by = "cell")
Mexpression_matrix1253 <- merge(expression_matrix1253, EST1253, by = "cell")
Mexpression_matrix1254 <- merge(expression_matrix1254, EST1254, by = "cell")
Mexpression_matrix1675 <- merge(expression_matrix1675, EST1675, by = "cell")
Mexpression_matrix0331[1:3,1:3]


Mexpression_matrix0331$model1 <- rowMeans(Mexpression_matrix0331[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
Mexpression_matrix0731$model1 <- rowMeans(Mexpression_matrix0731[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
Mexpression_matrix0733$model1 <- rowMeans(Mexpression_matrix0733[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
Mexpression_matrix0734$model1 <- rowMeans(Mexpression_matrix0734[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
Mexpression_matrix1251$model1 <- rowMeans(Mexpression_matrix1251[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
Mexpression_matrix1252$model1 <- rowMeans(Mexpression_matrix1252[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
Mexpression_matrix1253$model1 <- rowMeans(Mexpression_matrix1253[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
Mexpression_matrix1254$model1 <- rowMeans(Mexpression_matrix1254[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])
Mexpression_matrix1675$model1 <- rowMeans(Mexpression_matrix1675[, c("Angiogenesis", "EMT", "Invasion", "Metastasis")])



singlemt0331 <- scale(Mexpression_matrix0331[,-1])
singlemt0731 <- scale(Mexpression_matrix0731[,-1])
singlemt0733 <- scale(Mexpression_matrix0733[,-1])
singlemt0734 <- scale(Mexpression_matrix0734[,-1])
singlemt1251 <- scale(Mexpression_matrix1251[,-1])
singlemt1252 <- scale(Mexpression_matrix1252[,-1])
singlemt1253 <- scale(Mexpression_matrix1253[,-1])
singlemt1254 <- scale(Mexpression_matrix1254[,-1])
singlemt1675 <- scale(Mexpression_matrix1675[,-1])

cor_pearson0331 <- cor(singlemt0331, method = 'pearson')
cor_pearson0731 <- cor(singlemt0731, method = 'pearson')
cor_pearson0733 <- cor(singlemt0733, method = 'pearson')
cor_pearson0734 <- cor(singlemt0734, method = 'pearson')
cor_pearson1251 <- cor(singlemt1251, method = 'pearson')
cor_pearson1252 <- cor(singlemt1252, method = 'pearson')
cor_pearson1253 <- cor(singlemt1253, method = 'pearson')
cor_pearson1254 <- cor(singlemt1254, method = 'pearson')
cor_pearson1675 <- cor(singlemt1675, method = 'pearson')

write.csv(as.data.frame(cor_pearson0331[, c("model1", "Angiogenesis", "EMT", "Invasion", "Metastasis")]), file = "/wangzhenzenlab/fsh/xff/stdata/gene_cor0331.csv", row.names = TRUE)
write.csv(as.data.frame(cor_pearson0731[, c("model1", "Angiogenesis", "EMT", "Invasion", "Metastasis")]), file = "/wangzhenzenlab/fsh/xff/stdata/gene_cor0731.csv", row.names = TRUE)
write.csv(as.data.frame(cor_pearson0733[, c("model1", "Angiogenesis", "EMT", "Invasion", "Metastasis")]), file = "/wangzhenzenlab/fsh/xff/stdata/gene_cor0733.csv", row.names = TRUE)
write.csv(as.data.frame(cor_pearson0734[, c("model1", "Angiogenesis", "EMT", "Invasion", "Metastasis")]), file = "/wangzhenzenlab/fsh/xff/stdata/gene_cor0734.csv", row.names = TRUE)
write.csv(as.data.frame(cor_pearson1251[, c("model1", "Angiogenesis", "EMT", "Invasion", "Metastasis")]), file = "/wangzhenzenlab/fsh/xff/stdata/gene_cor1251.csv", row.names = TRUE)
write.csv(as.data.frame(cor_pearson1252[, c("model1", "Angiogenesis", "EMT", "Invasion", "Metastasis")]), file = "/wangzhenzenlab/fsh/xff/stdata/gene_cor1252.csv", row.names = TRUE)
write.csv(as.data.frame(cor_pearson1253[, c("model1", "Angiogenesis", "EMT", "Invasion", "Metastasis")]), file = "/wangzhenzenlab/fsh/xff/stdata/gene_cor1253.csv", row.names = TRUE)
write.csv(as.data.frame(cor_pearson1254[, c("model1", "Angiogenesis", "EMT", "Invasion", "Metastasis")]), file = "/wangzhenzenlab/fsh/xff/stdata/gene_cor1254.csv", row.names = TRUE)
write.csv(as.data.frame(cor_pearson1675[, c("model1", "Angiogenesis", "EMT", "Invasion", "Metastasis")]), file = "/wangzhenzenlab/fsh/xff/stdata/gene_cor1675.csv", row.names = TRUE)

# my_colors <- colorRampPalette(c("steelblue", "white", "#c11e21"))(100)
# corrplot(cor_pearson[c("FN1","COL1A1","COL1A2","COL3A1","COL6A3","SPARC","LGALS1","BGN","model1","Angiogenesis","EMT","Invasion","Metastasis"),
#                      c("FN1","COL1A1","COL1A2","COL3A1","COL6A3","SPARC","LGALS1","BGN","model1","Angiogenesis","EMT","Invasion","Metastasis")], 
#          method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)




singlemt0331 <- scale(Mexpression_matrix0331[,-1])
cor_pearson0331 <- cor(singlemt0331, method = 'pearson')



