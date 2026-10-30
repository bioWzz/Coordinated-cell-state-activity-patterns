setwd("D:/细胞状态2/论文/fig1")
load("D:/细胞状态2/论文/fig1/TUMOR3.RData")
head(TUMOR3)


# ggplot(allcellbar, aes(x = model, y = Count, fill = group)) +
#   geom_bar(stat = "identity", position = "stack") +
#   scale_fill_manual(values = c("up" = "#c11e21", "down" = "steelblue")) +
#   theme_minimal() +
#   labs(x = "model", y = "Count", fill = "Category", title = "Stacked Bar Chart")
# 
# dim(TUMOR3[TUMOR3$model1>0,])
# dim(TUMOR3[TUMOR3$model1<0,])
# dim(TUMOR3[TUMOR3$model2>0,])
# dim(TUMOR3[TUMOR3$model2<0,])
# dim(TUMOR3[TUMOR3$model3>0,])
# dim(TUMOR3[TUMOR3$model3<0,])
# dim(TUMOR3[TUMOR3$model4>0,])
# dim(TUMOR3[TUMOR3$model4<0,])
# dim(TUMOR3[TUMOR3$model5>0,])
# dim(TUMOR3[TUMOR3$model5<0,])
# dim(TUMOR3[TUMOR3$model6>0,])
# dim(TUMOR3[TUMOR3$model6<0,])
# dim(TUMOR3[TUMOR3$model7>0,])
# dim(TUMOR3[TUMOR3$model7<0,])
# 
# TUMOR3$celltype <- as.factor(TUMOR3$celltype)
# TUMOR3$celltype <- factor(TUMOR3$celltype,levels=c("B", "CD4Tconv", "CD8T", "CD8Tex", "Tprolif", "Treg", "DC", "Mast", "Mono/Macro", "NK", "Malignant", "Endothelial", "Epithelial", "Fibroblasts", "Myofibroblasts"))
# 
# allcellbar <- data.frame(
#   model = rep(c("model1", "model2", "model3", "model4", "model5", "model6", "model7"), each = 2),
#   group = rep(c("up", "down"), times = 7),
#   Count = c(859674, 497932, 265551, 1092055, 464691, 892915,
#             545271,812335,357101,1000505,1128077,229529,1223469,134137)
# )
# 
# 
# ggplot(allcellbar, aes(x = model, y = Count, fill = group)) +
#   geom_bar(stat = "identity", position = "stack") +
#   scale_fill_manual(values = c("up" = "#c11e21", "down" = "steelblue")) +
#   theme_minimal() +
#   labs(x = "model", y = "Count", fill = "Category", title = "Stacked Bar Chart")
# 
# 
# TUMOR3$group1 <- ifelse(TUMOR3$model1>median(TUMOR3$model1),"up","down")
# TUMOR3$group2 <- ifelse(TUMOR3$model2>median(TUMOR3$model2),"up","down")
# TUMOR3$group3 <- ifelse(TUMOR3$model3>median(TUMOR3$model3),"up","down")
# TUMOR3$group4 <- ifelse(TUMOR3$model4>median(TUMOR3$model4),"up","down")
# TUMOR3$group5 <- ifelse(TUMOR3$model5>median(TUMOR3$model5),"up","down")
# TUMOR3$group6 <- ifelse(TUMOR3$model6>median(TUMOR3$model6),"up","down")
# TUMOR3$group7 <- ifelse(TUMOR3$model7>median(TUMOR3$model7),"up","down")
# 
# 
# library(dplyr)
# table(TUMOR3$celltype)
# celltypebar <- TUMOR3%>%
#   group_by(celltype, group3) %>%
#   summarise(Count = n())
# 
# celltypebar2 <- celltypebar %>%
#   group_by(celltype) %>%
#   mutate(Count = Count / sum(Count)) %>%
#   ungroup()
# 
# pdf("D:/细胞状态2/论文/fig2/fig2b/fig2b3model1.pdf")
# ggplot(celltypebar2, aes(x = cancertype, y = Count, fill = group1)) +
#   geom_bar(stat = "identity", position = "stack") +
#   scale_fill_manual(values = c("up" = "#c11e21", "down" = "steelblue")) +
#   theme_minimal() +
#   labs(x = "model", y = "Count", fill = "Category", title = "Stacked Bar Chart")+
#   theme_minimal()+ 
#   theme(panel.grid = element_blank())
# 
# dev.off()


TUMOR3$TP <- rowMeans(TUMOR3[,c("Angiogenesis","EMT","Invasion","Metastasis")])
TUMOR3$SR <- rowMeans(TUMOR3[,c("Apoptosis","Hypoxia","Inflammation","Quiescence")])
TUMOR3$GS <- rowMeans(TUMOR3[,c("Cell Cycle","DNA damage","DNA repair")])
TUMOR3$FP <- rowMeans(TUMOR3[,c("Differentiation","Proliferation","Stemness")])

TUMOR3$groupTP <- ifelse(TUMOR3$TP>median(TUMOR3$TP),"high","low")
TUMOR3$groupSR <- ifelse(TUMOR3$SR>median(TUMOR3$SR),"high","low")
TUMOR3$groupGS <- ifelse(TUMOR3$GS>median(TUMOR3$GS),"high","low")
TUMOR3$groupFP <- ifelse(TUMOR3$FP>median(TUMOR3$FP),"high","low")

save(TUMOR3,file = "D:/细胞状态2/论文/fig2/newTUMOR3.RData")

library(dplyr)
table(TUMOR3$celltype)

#TP
celltypebar <- TUMOR3%>%
  group_by(celltype, groupTP) %>%
  summarise(Count = n())

celltypebar2 <- celltypebar %>%
  group_by(celltype) %>%
  mutate(Count = Count / sum(Count)) %>%
  ungroup()

pdf("D:/细胞状态2/论文/fig2/fig2b/fig2bTPcelltype.pdf")
ggplot(celltypebar2, aes(x = celltype, y = Count, fill = groupTP)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("high" = "#c11e21", "low" = "steelblue")) +
  theme_minimal() +
  labs(x = "model", y = "Count", fill = "Category", title = "Stacked Bar Chart")+
  theme_minimal()+ 
  theme(panel.grid = element_blank())
dev.off()

cancertypebar <- TUMOR3%>%
  group_by(cancertype, groupTP) %>%
  summarise(Count = n())

cancertypebar2 <- cancertypebar %>%
  group_by(cancertype) %>%
  mutate(Count = Count / sum(Count)) %>%
  ungroup()

pdf("D:/细胞状态2/论文/fig2/fig2b/fig2bTPcancertype.pdf")
ggplot(cancertypebar2, aes(x = cancertype, y = Count, fill = groupTP)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("high" = "#c11e21", "low" = "steelblue")) +
  theme_minimal() +
  labs(x = "model", y = "Count", fill = "Category", title = "Stacked Bar Chart")+
  theme_minimal()+ 
  theme(panel.grid = element_blank())
dev.off()



#SR
celltypebar <- TUMOR3%>%
  group_by(celltype, groupSR) %>%
  summarise(Count = n())

celltypebar2 <- celltypebar %>%
  group_by(celltype) %>%
  mutate(Count = Count / sum(Count)) %>%
  ungroup()

pdf("D:/细胞状态2/论文/fig2/fig2b/fig2bSRcelltype.pdf")
ggplot(celltypebar2, aes(x = celltype, y = Count, fill = groupSR)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("high" = "#c11e21", "low" = "steelblue")) +
  theme_minimal() +
  labs(x = "model", y = "Count", fill = "Category", title = "Stacked Bar Chart")+
  theme_minimal()+ 
  theme(panel.grid = element_blank())
dev.off()

cancertypebar <- TUMOR3%>%
  group_by(cancertype, groupSR) %>%
  summarise(Count = n())

cancertypebar2 <- cancertypebar %>%
  group_by(cancertype) %>%
  mutate(Count = Count / sum(Count)) %>%
  ungroup()

pdf("D:/细胞状态2/论文/fig2/fig2b/fig2bSRcancertype.pdf")
ggplot(cancertypebar2, aes(x = cancertype, y = Count, fill = groupSR)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("high" = "#c11e21", "low" = "steelblue")) +
  theme_minimal() +
  labs(x = "model", y = "Count", fill = "Category", title = "Stacked Bar Chart")+
  theme_minimal()+ 
  theme(panel.grid = element_blank())
dev.off()

#GS
celltypebar <- TUMOR3%>%
  group_by(celltype, groupGS) %>%
  summarise(Count = n())

celltypebar2 <- celltypebar %>%
  group_by(celltype) %>%
  mutate(Count = Count / sum(Count)) %>%
  ungroup()

pdf("D:/细胞状态2/论文/fig2/fig2b/fig2bGScelltype.pdf")
ggplot(celltypebar2, aes(x = celltype, y = Count, fill = groupGS)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("high" = "#c11e21", "low" = "steelblue")) +
  theme_minimal() +
  labs(x = "model", y = "Count", fill = "Category", title = "Stacked Bar Chart")+
  theme_minimal()+ 
  theme(panel.grid = element_blank())
dev.off()

cancertypebar <- TUMOR3%>%
  group_by(cancertype, groupGS) %>%
  summarise(Count = n())

cancertypebar2 <- cancertypebar %>%
  group_by(cancertype) %>%
  mutate(Count = Count / sum(Count)) %>%
  ungroup()

pdf("D:/细胞状态2/论文/fig2/fig2b/fig2bGScancertype.pdf")
ggplot(cancertypebar2, aes(x = cancertype, y = Count, fill = groupGS)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("high" = "#c11e21", "low" = "steelblue")) +
  theme_minimal() +
  labs(x = "model", y = "Count", fill = "Category", title = "Stacked Bar Chart")+
  theme_minimal()+ 
  theme(panel.grid = element_blank())
dev.off()

#FP
celltypebar <- TUMOR3%>%
  group_by(celltype, groupFP) %>%
  summarise(Count = n())

celltypebar2 <- celltypebar %>%
  group_by(celltype) %>%
  mutate(Count = Count / sum(Count)) %>%
  ungroup()

pdf("D:/细胞状态2/论文/fig2/fig2b/fig2bFPcelltype.pdf")
ggplot(celltypebar2, aes(x = celltype, y = Count, fill = groupFP)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("high" = "#c11e21", "low" = "steelblue")) +
  theme_minimal() +
  labs(x = "model", y = "Count", fill = "Category", title = "Stacked Bar Chart")+
  theme_minimal()+ 
  theme(panel.grid = element_blank())
dev.off()

cancertypebar <- TUMOR3%>%
  group_by(cancertype, groupFP) %>%
  summarise(Count = n())

cancertypebar2 <- cancertypebar %>%
  group_by(cancertype) %>%
  mutate(Count = Count / sum(Count)) %>%
  ungroup()

pdf("D:/细胞状态2/论文/fig2/fig2b/fig2bFPcancertype.pdf")
ggplot(cancertypebar2, aes(x = cancertype, y = Count, fill = groupFP)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("high" = "#c11e21", "low" = "steelblue")) +
  theme_minimal() +
  labs(x = "model", y = "Count", fill = "Category", title = "Stacked Bar Chart")+
  theme_minimal()+ 
  theme(panel.grid = element_blank())
dev.off()







#fig2c
#查看样本的上调和下调
# TUMOR3$Patient
# samplestate <- TUMOR3 %>%
#   group_by(Patient,cancertype) %>%
#   summarise(Sample1 = ifelse(sum(model1 > median(TUMOR3$model1))>sum(model1 < median(TUMOR3$model1)),"high","low"),
#             Sample2 = ifelse(sum(model2 > median(TUMOR3$model2))>sum(model2 < median(TUMOR3$model2)),"high","low"),
#             Sample3 = ifelse(sum(model3 > median(TUMOR3$model3))>sum(model3 < median(TUMOR3$model3)),"high","low"),
#             Sample4 = ifelse(sum(model4 > median(TUMOR3$model4))>sum(model4 < median(TUMOR3$model4)),"high","low"),
#             Sample5 = ifelse(sum(model5 > median(TUMOR3$model5))>sum(model5 < median(TUMOR3$model5)),"high","low"),
#             Sample6 = ifelse(sum(model6 > median(TUMOR3$model6))>sum(model6 < median(TUMOR3$model6)),"high","low"),
#             Sample7 = ifelse(sum(model7 > median(TUMOR3$model7))>sum(model7 < median(TUMOR3$model7)),"high","low")
#             )
# samplestate <- as.data.frame(samplestate)
# head(samplestate)
# 
# samplestate_long <- tidyr::pivot_longer(samplestate, cols = c(3:9), names_to = "Condition", values_to = "Value")
# #samplestate_long$group <- paste(samplestate_long$Condition,samplestate_long$Value,sep = "")
# 
# 
# samplestate2 <- samplestate_long %>%
#   group_by(Condition,Value) %>%
#   summarise(Count = n())
# library(tidyr)
# install.packages("ggalt")
# library(ggalt)
# samplestate3 <- pivot_wider(samplestate2, 
#                          id_cols = Condition, 
#                          names_from = Value , 
#                          values_from = Count)
# samplestate3$low <- -samplestate3$low
# 
# ggplot(aes(x=high,xend=low,y=Condition),data=samplestate3)+
#   geom_dumbbell(colour_x = "#FFB6C1",colour_xend = "#4169E1",size_x = 2,size_xend = 2,size=0.5,color="gray")+
#   theme_light()+
#   theme(panel.grid.minor.x =element_blank(),
#   )+
#   xlab("Annual Salary (in thousands)")+
#   theme_minimal()+ 
#   theme(panel.grid = element_blank())

# head(TUMOR3)
# samplestate
# 
# samplemodel1up <- as.data.frame(table(samplestate[samplestate$Sample1=="high",]$cancertype))
# samplemodel1down <- as.data.frame(table(samplestate[samplestate$Sample1=="low",]$cancertype))
# samplemodel2up <- as.data.frame(table(samplestate[samplestate$Sample2=="high",]$cancertype))
# samplemodel2down <- as.data.frame(table(samplestate[samplestate$Sample2=="low",]$cancertype))
# samplemodel3up <- as.data.frame(table(samplestate[samplestate$Sample3=="high",]$cancertype))
# samplemodel3down <- as.data.frame(table(samplestate[samplestate$Sample3=="low",]$cancertype))
# samplemodel4up <- as.data.frame(table(samplestate[samplestate$Sample4=="high",]$cancertype))
# samplemodel4down <- as.data.frame(table(samplestate[samplestate$Sample4=="low",]$cancertype))
# samplemodel5up <- as.data.frame(table(samplestate[samplestate$Sample5=="high",]$cancertype))
# samplemodel5down <- as.data.frame(table(samplestate[samplestate$Sample5=="low",]$cancertype))
# samplemodel6up <- as.data.frame(table(samplestate[samplestate$Sample6=="high",]$cancertype))
# samplemodel6down <- as.data.frame(table(samplestate[samplestate$Sample6=="low",]$cancertype))
# samplemodel7up <- as.data.frame(table(samplestate[samplestate$Sample7=="high",]$cancertype))
# samplemodel7down <- as.data.frame(table(samplestate[samplestate$Sample7=="low",]$cancertype))
# 
# 
# ggplot(samplemodel7down, aes(x = "", y = Freq, fill = Var1)) +
#   geom_bar(stat = "identity", width = 1) +
#   coord_polar(theta = "y") +
#   theme_void() +  # 移除坐标轴 
#   labs(title = "Pie Chart from Data Frame", fill = "Var1") +
#   scale_fill_manual(values = c("BCC"="#F8766C","BRCA"="#c11e21","CESC"="#E48600",
#                                 "CRC"="#E7B781","ESCA"="#44bd28","KIRC"="steelblue",
#                                 "LIHC"="#B883FF","NSCLC"="#A3B5D6","OV"="#3C6882",
#                                 "PAAD"="#79cedd","PRAD"="#70C1B3","SKCM"="#73B3BE",
#                                 "STAD"="#6D9EC3","THCA"="#c6c322","UCEC"="#A7C66B"))
# 





#NEWfig2c
samplestate <- TUMOR3 %>%
  group_by(Patient,cancertype) %>%
  summarise(SampleTP = ifelse(sum(TP > median(TUMOR3$TP))>sum(TP < median(TUMOR3$TP)),"high","low"),
            SampleSR = ifelse(sum(SR > median(TUMOR3$SR))>sum(SR < median(TUMOR3$SR)),"high","low"),
            SampleGS = ifelse(sum(GS > median(TUMOR3$GS))>sum(GS < median(TUMOR3$GS)),"high","low"),
            SampleFP = ifelse(sum(FP > median(TUMOR3$FP))>sum(FP < median(TUMOR3$FP)),"high","low")
  )

samplestate <- as.data.frame(samplestate)
head(samplestate)

samplestate_long <- tidyr::pivot_longer(samplestate, cols = c(3:6), names_to = "Condition", values_to = "Value")
samplestate_long


samplestate2 <- samplestate_long %>%
  group_by(Condition,Value) %>%
  summarise(Count = n())
library(tidyr)
install.packages("ggalt")
library(ggalt)
samplestate3 <- pivot_wider(samplestate2, 
                            id_cols = Condition, 
                            names_from = Value , 
                            values_from = Count)
samplestate3$low <- -samplestate3$low

ggplot(aes(x=high,xend=low,y=Condition),data=samplestate3)+
  geom_dumbbell(colour_x = "#FFB6C1",colour_xend = "#4169E1",size_x = 2,size_xend = 2,size=0.5,color="gray")+
  theme_light()+
  theme(panel.grid.minor.x =element_blank(),
  )+
  xlab("Annual Salary (in thousands)")+
  theme_minimal()+ 
  theme(panel.grid = element_blank())


SampleFPhigh <- as.data.frame(table(samplestate[samplestate$SampleFP=="high",]$cancertype))
SampleFPlow <- as.data.frame(table(samplestate[samplestate$SampleFP=="low",]$cancertype))
SampleGShigh <- as.data.frame(table(samplestate[samplestate$SampleGS=="high",]$cancertype))
SampleGSlow <- as.data.frame(table(samplestate[samplestate$SampleGS=="low",]$cancertype))
SampleSRhigh <- as.data.frame(table(samplestate[samplestate$SampleSR=="high",]$cancertype))
SampleSRlow <- as.data.frame(table(samplestate[samplestate$SampleSR=="low",]$cancertype))
SampleTPhigh <- as.data.frame(table(samplestate[samplestate$SampleTP=="high",]$cancertype))
SampleTPlow <- as.data.frame(table(samplestate[samplestate$SampleTP=="low",]$cancertype))


ggplot(SampleTPlow, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +  # 移除坐标轴 
  labs(title = "Pie Chart from Data Frame", fill = "Var1") +
  scale_fill_manual(values = c("BCC"="#F8766C","BRCA"="#c11e21","CESC"="#E48600",
                               "CRC"="#E7B781","ESCA"="#44bd28","KIRC"="steelblue",
                               "LIHC"="#B883FF","NSCLC"="#A3B5D6","OV"="#3C6882",
                               "PAAD"="#79cedd","PRAD"="#70C1B3","SKCM"="#73B3BE",
                               "STAD"="#6D9EC3","THCA"="#c6c322","UCEC"="#A7C66B"))


