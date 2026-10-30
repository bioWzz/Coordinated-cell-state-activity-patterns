library(tidyr)
head(PAADfold)
PAADfoldloing1 <- PAADfold[,c("Angiogenesis","EMT","Invasion","Metastasis","celltype","model1")]
PAADfoldloing11 <- pivot_longer(
  data = PAADfoldloing1,
  cols = -celltype,  # 选择除了id之外的所有列进行转换
  names_to = "variable",
  values_to = "value"
)


ggplot(data = PAADfoldloing11, aes(x = celltype, y = value, fill = variable)) +
  geom_boxplot(width=0.5,
               position = position_dodge(0.8),
               outlier.shape = NA)+  #每组的小提琴图之间的距离
  scale_fill_manual(values = c( "steelblue","steelblue","steelblue","steelblue","#c11e21"))+
  theme_prism()+
  theme(axis.text.x = element_text(size = 20,angle = 30, hjust = 0.5, vjust = 0.5))


PAADfoldloing4 <- PAADfold[,c("Inflammation","Quiescence","celltype","model4")]
PAADfoldloing41 <- pivot_longer(
  data = PAADfoldloing4,
  cols = -celltype,  # 选择除了id之外的所有列进行转换
  names_to = "variable",
  values_to = "value"
)
PAADfoldloing41$variable <- factor(PAADfoldloing41$variable,levels=c("Inflammation","Quiescence","model4"))

ggplot(data = PAADfoldloing41, aes(x = celltype, y = value, fill = variable)) +
  geom_boxplot(width=0.5,
               position = position_dodge(0.8),
               outlier.shape = NA)+  #每组的小提琴图之间的距离
  scale_fill_manual(values = c( "steelblue","steelblue","#c11e21"))+
  theme_prism()+
  theme(axis.text.x = element_text(size = 20,angle = 30, hjust = 0.5, vjust = 0.5))



