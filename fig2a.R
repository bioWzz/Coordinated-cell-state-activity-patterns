load("D:/细胞状态2/论文/fig1/TUMOR3.RData")
library(ggtree)
library(ape)
library(treeio)
library(ggplot2)

#合成model值
head(TUMOR3)
table(TUMOR3$celltype)
TUMOR3$model1 <- rowMeans(TUMOR3[,c("Angiogenesis","EMT","Invasion","Metastasis")])
TUMOR3$model2 <- TUMOR3$Apoptosis
TUMOR3$model3 <- TUMOR3$Hypoxia
TUMOR3$model4 <- rowMeans(TUMOR3[,c("Inflammation","Quiescence")])
TUMOR3$model5 <- rowMeans(TUMOR3[,c("Cell Cycle","DNA damage","DNA repair")])
TUMOR3$model6 <- TUMOR3$Proliferation
TUMOR3$model7 <- rowMeans(TUMOR3[,c("Differentiation","Stemness")])

TUMOR3$type <- paste(TUMOR3$cancertype, TUMOR3$celltype, sep = "_")

TUMOR3_means <- aggregate(. ~ type, data = TUMOR3[,32:39], mean)
head(TUMOR3_means)
class(TUMOR3_means)
rownames(TUMOR3_means) <- TUMOR3_means$type
treetumor <- TUMOR3_means[,2:8]

#分别按细胞和癌症染色
group1 <- TUMOR3[,c("celltype","type")]
group2 <- TUMOR3[,c("cancertype","type")]
colnames(group1) <- c("group","label")
colnames(group2) <- c("group","label")
group1 <- unique(group1)
group2 <- unique(group2)
head(group1)
table(group1$group)
table(group2$group)

hctumor <- hclust(dist(treetumor),method = "ave")
treeall = as.phylo(hctumor)
treealla1 <- full_join(treeall,group1,by="label")
treealla2 <- full_join(treeall,group2,by="label")


pdf("D:/细胞状态2/论文/fig2/fig2a/fig2a1.pdf",width = 30,height = 29)
ggtree(treealla1,layout = "circular",aes(color=group))+
  geom_tiplab2(offset=0.01,size=8)+
  scale_color_manual(values = c("#F8766C","#c11e21","#E48600","#E7B781","#44bd28",
                                "steelblue","#B883FF","#A3B5D6","#3C6882","#BFD3E6",
                                "#70C1B3","#73B3BE","#6D9EC3","#c6c322","#A7C66B"))
dev.off()

pdf("D:/细胞状态2/论文/fig2/fig2a/fig2a2.pdf",width = 30,height = 29)
ggtree(treealla2,layout = "circular",aes(color=group))+
  geom_tiplab2(offset=0.01,size=8)+
  scale_color_manual(values = c("#F8766C","#c11e21","#E48600","#E7B781","#44bd28",
                                "steelblue","#B883FF","#A3B5D6","#3C6882","#BFD3E6",
                                "#70C1B3","#73B3BE","#6D9EC3","#c6c322","#A7C66B"))
dev.off()

save(TUMOR3,file = "D:/细胞状态2/论文/fig1/TUMOR3.RData")
