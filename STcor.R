#PAAD

paaddealT2 <- read.csv("D:/细胞状态2/论文/空间/paad/paaddealT2 (2).csv")
paadstmetadata <- read.csv("D:/细胞状态2/论文/空间/paad/denggaoxian.csv")
paadstmetadata$result
paadstmetadata$X <- NULL
paadstmetadata$X.1 <- NULL
head(paaddealT2)
head(paadstmetadata)
paaddealT3 <- merge(paaddealT2,paadstmetadata,by.x="X",by.y="cell")
head(paaddealT3)

singlemt <- paaddealT3[,c(2:22,44)]
#singlemt <- paaddealT3[,c("MMP2","model1")]
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
library(corrplot)
my_colors <- colorRampPalette(c("steelblue", "white", "#c11e21"))(100)
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)

paadcor2<- as.data.frame(cor_pearson)
paadcor2$gene <- rownames(paadcor2)
paadcor3 <- paadcor2[,c("gene","model1")]
paadcor3 <- paadcor3[order(paadcor3$model1),]

write.csv(paadcor3,file = "D:/细胞状态2/论文/空间/paad/cor.csv")




#STAD
staddealT2 <- read.csv("D:/细胞状态2/论文/空间/stad/staddealT2 (2).csv")
stadstmeta <- read.csv("D:/细胞状态2/论文/空间/stad/stad/denggaoxian.csv")
head(staddealT2)
head(stadstmeta)
stadstmeta$X <- NULL
staddealT3 <- merge(staddealT2,stadstmeta,by.x="X",by.y="cell")
head(staddealT3)


singlemt <- staddealT3[,c(2:9,31)]
singlemt <- scale(singlemt)
dim(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
library(corrplot)
my_colors <- colorRampPalette(c("steelblue", "white", "#c11e21"))(100)
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)

stadcor2<- as.data.frame(cor_pearson)
stadcor2$gene <- rownames(stadcor2)
stadcor3 <- stadcor2[,c("gene","model1")]
stadcor3 <- stadcor3[order(stadcor3$model1),]

write.csv(stadcor3,file = "D:/细胞状态2/论文/空间/stad/stad/cor.csv")






# FN1+ COL1A1 +COL1A2  +COL3A1 + COL6A3+ SPARC +LGALS1 +BGN
