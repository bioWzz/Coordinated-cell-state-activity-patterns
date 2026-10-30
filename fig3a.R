load("D:/细胞状态2/论文/fig2/newTUMOR3.RData")
head(TUMOR3)
table(TUMOR3$cancertype)
table(TUMOR3$celltype)

CESCsingle <- TUMOR3[TUMOR3$cancertype=="CESC",]
#ESCAsingle <- TUMOR3[TUMOR3$cancertype=="ESCA",]
KIRCsingle <- TUMOR3[TUMOR3$cancertype=="KIRC",]
PAADsingle <- TUMOR3[TUMOR3$cancertype=="PAAD",]
#SKCMsingle <- TUMOR3[TUMOR3$cancertype=="SKCM",]
STADsingle <- TUMOR3[TUMOR3$cancertype=="STAD",]
UCECsingle <- TUMOR3[TUMOR3$cancertype=="UCEC",]

#save(PAADsingle,SKCMsingle,STADsingle,file="D:/细胞状态2/论文/fig3/threecancersingle.RData")
save(PAADsingle,STADsingle,UCECsingle,file = "D:/细胞状态2/论文/fig3/newthreecancersingle.RData")

cescTP <- CESCsingle
kircTP <- KIRCsingle
kircSR <- KIRCsingle
kircFP <- KIRCsingle
paadTP <- PAADsingle
stadTP <- STADsingle
stadGS <- STADsingle
stadFP <- STADsingle
ucecTP <- UCECsingle
ucecFP <- UCECsingle

cescTP$modeltype <- "cescTP"
kircTP$modeltype <- "kircTP"
kircSR$modeltype <- "kircSR"
kircFP$modeltype <- "kircFP"
paadTP$modeltype <- "paadTP"
stadTP$modeltype <- "stadTP"
stadGS$modeltype <- "stadGS"
stadFP$modeltype <- "stadFP"
ucecTP$modeltype <- "ucecTP"
ucecFP$modeltype <- "ucecFP"


cescTP$GROUP <- ifelse(cescTP$TP>median(cescTP$TP),"high","low")
kircTP$GROUP <- ifelse(kircTP$TP>median(kircTP$TP),"high","low")
kircSR$GROUP <- ifelse(kircSR$SR>median(kircSR$SR),"high","low")
kircFP$GROUP <- ifelse(kircFP$FP>median(kircFP$FP),"high","low")
paadTP$GROUP <- ifelse(paadTP$TP>median(paadTP$TP),"high","low")
stadTP$GROUP <- ifelse(stadTP$TP>median(stadTP$TP),"high","low")
stadGS$GROUP <- ifelse(stadGS$GS>median(stadGS$GS),"high","low")
stadFP$GROUP <- ifelse(stadFP$FP>median(stadFP$FP),"high","low")
ucecTP$GROUP <- ifelse(ucecTP$TP>median(ucecTP$TP),"high","low")
ucecFP$GROUP <- ifelse(ucecFP$FP>median(ucecFP$FP),"high","low")

allgroup <- rbind(cescTP,
                  kircTP,kircSR,kircFP,
                  paadTP,
                  stadTP,stadGS,stadFP,
                  ucecTP,ucecFP)

allgroup$modeltype <- as.factor(allgroup$modeltype)
allgroup$modeltype <- factor(allgroup$modeltype,levels = c("cescTP","kircTP","kircSR","kircFP",
                                                           "paadTP","stadTP","stadGS","stadFP",
                                                           "ucecTP","ucecFP") )


highgroup <- allgroup[allgroup$GROUP=="high",]
lowgroup <- allgroup[allgroup$GROUP=="low",]


highbar <- highgroup%>%
  group_by(celltype, modeltype) %>%
  summarise(Count = n())

highbar2 <- highbar %>%
  group_by(modeltype) %>%
  mutate(Count = Count / sum(Count)) %>%
  ungroup()

ggplot(highbar2, aes(x = modeltype, y = Count, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("B"="#F8766C","CD4Tconv"="#c11e21","CD8T"="#E48600",
                               "CD8Tex"="#E7B781","DC"="#44bd28","Endothelial"="steelblue",
                               "Epithelial"="#B883FF","Fibroblasts"="#A3B5D6","Malignant"="#3C6882",
                               "Mast"="#BFD3E6","Mono/Macro"="#70C1B3","Myofibroblasts"="#73B3BE",
                               "NK"="#6D9EC3","Tprolif"="#c6c322","Treg"="#A7C66B")) +
  theme_minimal() +
  labs(x = "model", y = "Count", fill = "Category", title = "Stacked Bar Chart")+
  theme_minimal()+ 
  theme(panel.grid = element_blank())



lowbar <- lowgroup%>%
  group_by(celltype, modeltype) %>%
  summarise(Count = n())

lowbar2 <- lowbar %>%
  group_by(modeltype) %>%
  mutate(Count = Count / sum(Count)) %>%
  ungroup()

ggplot(lowbar2, aes(x = modeltype, y = Count, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("B"="#F8766C","CD4Tconv"="#c11e21","CD8T"="#E48600",
                               "CD8Tex"="#E7B781","DC"="#44bd28","Endothelial"="steelblue",
                               "Epithelial"="#B883FF","Fibroblasts"="#A3B5D6","Malignant"="#3C6882",
                               "Mast"="#BFD3E6","Mono/Macro"="#70C1B3","Myofibroblasts"="#73B3BE",
                               "NK"="#6D9EC3","Tprolif"="#c6c322","Treg"="#A7C66B")) +
  theme_minimal() +
  labs(x = "model", y = "Count", fill = "Category", title = "Stacked Bar Chart")+
  theme_minimal()+ 
  theme(panel.grid = element_blank())

save(allgroup,file="D:/细胞状态2/论文/fig3/fig3a/newallgroup.RData")

