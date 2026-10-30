load("D:/细胞状态2/论文/fig2/newTUMOR3.RData")
load("D:/细胞状态2/论文/fig3/fig3a/allgroup.RData")
library(dplyr)
library(ggplot2)

cescTP$value <- cescTP$TP
kircTP$value <- kircTP$TP
kircSR$value <- kircSR$SR
kircFP$value <- kircFP$FP
paadTP$value <- paadTP$TP
stadTP$value <- stadTP$TP
stadGS$value <- stadGS$GS
stadFP$value <- stadFP$FP
ucecTP$value <- ucecTP$TP
ucecFP$value <- ucecFP$FP

allgroup2 <- rbind(cescTP,kircTP,kircSR,kircFP,
                   paadTP,stadTP,stadGS,stadFP,
                   ucecTP,ucecFP)

allgroup2$modeltype <- as.factor(allgroup2$modeltype)
allgroup2$modeltype <- factor(allgroup2$modeltype,levels = c("cescTP","kircTP","kircSR","kircFP",
                                                           "paadTP","stadTP","stadGS","stadFP",
                                                           "ucecTP","ucecFP") )

head(allgroup2)
manyviolin <- allgroup2[,c("modeltype","celltype","value")]
head(manyviolin)

lines_data <- manyviolin %>%
  group_by(modeltype) %>%
  summarise(y_value = median(value)) %>%
  ungroup()


ggplot(manyviolin, aes(x = celltype, y = value, fill = celltype)) +
  geom_violin(trim = FALSE) +  # 添加小提琴图层，设置trim = FALSE以显示所有数据点
  scale_fill_brewer(palette = "Paired") +  # 使用颜色方案
  facet_wrap(~modeltype, ncol = 1)+
  scale_color_manual(values = c("#FB8072", "#c11e21", "#E48600", "#E7B781", "#c6c322", 
                                "#A7C66B", "#44bd28", "#BFD3E6", "#70C1B3", "#6D9EC3",
                                "#3C6882", "#70A3CC", "#B883FF", "#A3B5D6", "#73B3BE"))+
  scale_fill_manual(values = c("#FB8072", "#c11e21", "#E48600", "#E7B781", "#c6c322", 
                               "#A7C66B", "#44bd28", "#BFD3E6", "#70C1B3", "#6D9EC3",
                               "#3C6882", "#70A3CC", "#B883FF", "#A3B5D6", "#73B3BE"))+
  theme_bw() + 
  theme(panel.grid=element_blank())+
  geom_hline(data = lines_data, aes(yintercept = y_value), linetype = "dashed", color = "black") 



c("B"="#F8766C","CD4Tconv"="#c11e21","CD8T"="#E48600","CD8Tex"="#E7B781","Tprolif"="#c6c322",
  "Treg"="#A7C66B","DC"="#44bd28",,"Mast"="#BFD3E6","Mono/Macro"="#70C1B3","NK"="#6D9EC3",
  "Malignant"="#3C6882","Endothelial"="steelblue","Epithelial"="#B883FF","Fibroblasts"="#A3B5D6","Myofibroblasts"="#73B3BE")









prop.table(table(cescTP$celltype,cescTP$GROUP),margin = 1)

prop.table(table(kircTP$celltype,kircTP$GROUP),margin = 1)
prop.table(table(kircSR$celltype,kircSR$GROUP),margin = 1)
prop.table(table(kircFP$celltype,kircFP$GROUP),margin = 1)

prop.table(table(paadTP$celltype,paadTP$GROUP),margin = 1)

prop.table(table(stadTP$celltype,stadTP$GROUP),margin = 1)
prop.table(table(stadGS$celltype,stadGS$GROUP),margin = 1)
prop.table(table(stadFP$celltype,stadFP$GROUP),margin = 1)

prop.table(table(ucecTP$celltype,ucecTP$GROUP),margin = 1)
prop.table(table(ucecFP$celltype,ucecFP$GROUP),margin = 1)

