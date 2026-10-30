estotal <- read.csv("D:/细胞状态2/bulk/total_ES.csv")
clin_all <- read.csv("D:/细胞状态2/bulk/clinical/clinical-all.csv")
estotal[1:3,1:3]
clin_all[1:3,1:3]
estotal$X.1 <- NULL
rownames(estotal) <- estotal$X
estotal$X <- NULL
estotalT <- t(estotal)
estotalT <- as.data.frame(estotalT)
head(estotalT)
estotalT$patient <- rownames(estotalT)
estotalT$patient <- sub("\\.$", "", estotalT$patient)
clin_all$patient <- gsub("-", ".", clin_all$X_PATIENT)

totalESTclin <- merge(clin_all,estotalT,by="patient")
totalESTclin$cancertype <- sub(".*\\((.*?)\\).*", "\\1", totalESTclin$X_cohort)
table(totalESTclin$cancertype)

intersectcancer<- intersect(names(table(totalESTclin$cancertype)),names(table(TUMOR3$cancertype)))
#"BRCA" "CESC" "ESCA" "KIRC" "LIHC" "OV"   "PAAD" "PRAD" "SKCM" "STAD" "THCA" "UCEC"
totalESTclin2 <- totalESTclin[totalESTclin$cancertype%in%intersectcancer,]

# totalESTclin2$model1 <- rowMeans(totalESTclin2[,c("Angiogenesis","EMT","Invasion","Metastasis")])
# totalESTclin2$model2 <- totalESTclin2$Apoptosis
# totalESTclin2$model3 <- totalESTclin2$Hypoxia
# totalESTclin2$model4 <- rowMeans(totalESTclin2[,c("Inflammation","Quiescence")])
# totalESTclin2$model5 <- rowMeans(totalESTclin2[,c("Cell Cycle","DNA damage","DNA repair")])
# totalESTclin2$model6 <- totalESTclin2$Proliferation
# totalESTclin2$model7 <- rowMeans(totalESTclin2[,c("Differentiation","Stemness")])

totalESTclin2$TP <- rowMeans(totalESTclin2[,c("Angiogenesis","EMT","Invasion","Metastasis")])
totalESTclin2$SR <- rowMeans(totalESTclin2[,c("Apoptosis","Hypoxia","Inflammation","Quiescence")])
totalESTclin2$GS <- rowMeans(totalESTclin2[,c("Cell Cycle","DNA damage","DNA repair")])
totalESTclin2$FP <- rowMeans(totalESTclin2[,c("Proliferation","Differentiation","Stemness")])

# totalESTclin2$group1 <- ifelse(totalESTclin2$model1> median(totalESTclin2$model1),"high" ,"low" )
# totalESTclin2$group2 <- ifelse(totalESTclin2$model2> median(totalESTclin2$model2),"high" ,"low" )
# totalESTclin2$group3 <- ifelse(totalESTclin2$model3> median(totalESTclin2$model3),"high" ,"low" )
# totalESTclin2$group4 <- ifelse(totalESTclin2$model4> median(totalESTclin2$model4),"high" ,"low" )
# totalESTclin2$group5 <- ifelse(totalESTclin2$model5> median(totalESTclin2$model5),"high" ,"low" )
# totalESTclin2$group6 <- ifelse(totalESTclin2$model6> median(totalESTclin2$model6),"high" ,"low" )
# totalESTclin2$group7 <- ifelse(totalESTclin2$model7> median(totalESTclin2$model7),"high" ,"low" )

totalESTclin2$groupTP <- ifelse(totalESTclin2$TP> median(totalESTclin2$TP),"high" ,"low" )
totalESTclin2$groupSR <- ifelse(totalESTclin2$SR> median(totalESTclin2$SR),"high" ,"low" )
totalESTclin2$groupGS <- ifelse(totalESTclin2$GS> median(totalESTclin2$GS),"high" ,"low" )
totalESTclin2$groupFP <- ifelse(totalESTclin2$FP> median(totalESTclin2$FP),"high" ,"low" )

totalESTclin2 <- totalESTclin2[!is.na(totalESTclin2$X_OS),]

library(survival)
library("survminer")
surv_object <- survfit(Surv(X_OS, X_EVENT) ~ groupTP, data=totalESTclin2)

ggsurvplot(surv_object,
           pval = TRUE, conf.int = F,
           risk.table = F,
           palette=c("#4b6b86","#ab2126" ))

for (i in intersectcancer) {
  #i <- "CESC"
  print(i)
  singlebulk2 <- totalESTclin2[totalESTclin2$cancertype==i,]
  singlebulk2$groupTP <- ifelse(singlebulk2$TP> median(singlebulk2$TP),"high" ,"low" )
  singlebulk2$groupSR <- ifelse(singlebulk2$SR> median(singlebulk2$SR),"high" ,"low" )
  singlebulk2$groupGS <- ifelse(singlebulk2$GS> median(singlebulk2$GS),"high" ,"low" )
  singlebulk2$groupFP <- ifelse(singlebulk2$FP> median(singlebulk2$FP),"high" ,"low" )
  
  for (m in c("groupTP","groupSR","groupGS","groupFP")) {
    #m <- "group5"
    myfilename <- paste("D:\\细胞状态2\\论文\\fig2\\fig2d\\gsvasur\\",i,"\\",m,".pdf",sep="")
    pdf(myfilename,onefile = F)
    surv_object <- survfit(Surv(X_OS, X_EVENT) ~ singlebulk2[,m], data=singlebulk2)
    p <- ggsurvplot(surv_object,
                    pval = TRUE, conf.int = F,
                    risk.table = F,
                    palette=c("#ab2126","#4b6b86" ))
    print(p)
    dev.off()
  }
}



#气泡图
# gsvamodeltype <- c("model1", "model2", "model3", "model4", "model5", "model6", "model7")
# gsvacancer_type <- names(table(totalESTclin2$cancertype))
# gsvap_value <- matrix(c(1,1,1,1,1,1,1,                                  #BRCA  
#                         0.004,1,0.0001,1,1,1,1,                         #CESC
#                         1,1,1,0.028,1,1,1,                              #ESCA
#                         0.0073,0.0036,1,0.0001,1,0.0001,1,              #KIRC
#                         1,1,1,1,1,1,1,                                  #LIHC
#                         1,1,1,1,1,1,1,                                  #OV
#                         0.026,1,0.0093,1,1,1,1,                         #PAAD
#                         1,1,1,1,1,1,1,                                  #PRAD
#                         1,1,1,0.0044,1,0.0034,1,                        #SKCM
#                         0.011,1,0.013,1,0.044,1,0.0016,            #STAD
#                         1,1,1,1,1,1,1,                                  #THCA
#                         0.035,1,1,1,1,1,1                        #UCEC
# ), nrow = 7, ncol = 12)  # 随机生成的p值矩阵，7行（模式）
# gsvagroup <- matrix(c(NA,NA,NA,NA,NA,NA,NA,                         #BRCA  
#                       "up",NA,"up",NA,NA,NA,NA,                     #CESC
#                       NA,NA,NA,"up",NA,NA,NA,                       #ESCA
#                       "up","up",NA,"up",NA,"up",NA,                 #KIRC
#                       NA,NA,NA,NA,NA,NA,NA,                         #LIHC
#                       NA,NA,NA,NA,NA,NA,NA,                         #OV
#                       "up",NA,"up",NA,NA,NA,NA,                     #PAAD
#                       NA,NA,NA,NA,NA,NA,NA,                         #PRAD
#                       NA,NA,NA,"down",NA,"down",NA,                 #SKCM
#                       "up",NA,"up",NA,"down",NA,"up",               #STAD
#                       NA,NA,NA,NA,NA,NA,NA,                         #THCA
#                       "up",NA,NA,NA,NA,NA,NA                        #UCEC
# ), nrow = 7, ncol = 12)
# 
# gsvaqipaotu <- expand.grid(model = gsvamodeltype, cancertype = gsvacancer_type)
# gsvaqipaotu$p_value <- as.vector(gsvap_value)
# gsvaqipaotu$group <- as.vector(gsvagroup)
# gsvaqipaotu <- gsvaqipaotu[!gsvaqipaotu$cancertype %in% c("BRCA","LIHC","OV","PRAD","THCA"),]
# gsvaqipaotu$log10pvalue <- -log10(gsvaqipaotu$p_value)
# 
# ggplot(gsvaqipaotu, aes(x = model, y = cancertype)) +
#   geom_point(aes(size = log10pvalue, color = group,fill=group), shape = 21) +
#   scale_size_continuous(range = c(1, 10)) +  # 调整气泡大小范围
#   scale_color_manual(values = c("down" = "steelblue", "up" = "#c11e21"),na.value = "white") +  # 设置分组颜色
#   scale_fill_manual(values = c("down" = "steelblue", "up" = "#c11e21"), na.value = "white") + 
#   labs(title = "Bubble Heatmap of p-values", x = "model", y = "CancerTypes", size = "-log10 pvalue", color = "Group") +
#   theme_minimal()+
#   theme_bw() + 
#   theme(panel.grid=element_blank())

#newfig2d
gsvamodeltype <- c("TP", "SR", "GS", "FP")
gsvacancer_type <- names(table(totalESTclin2$cancertype))
gsvap_value <- matrix(c(1,1,1,1,                           #BRCA  
                        0.004,1,1,1,                       #CESC
                        1,1,1,1,                           #ESCA
                        0.0073,0.0012,1,0.049,             #KIRC
                        1,1,1,1,                           #LIHC
                        1,1,1,1,                           #OV
                        0.026,1,1,1,                       #PAAD
                        1,1,1,1,                           #PRAD
                        1,1,1,1,                           #SKCM
                        0.011,1,0.044,0.0018,              #STAD
                        1,1,1,1,                           #THCA
                        0.035,1,1,0.0049                   #UCEC
), nrow = 4, ncol = 12)  # 生成的p值矩阵

gsvagroup <- matrix(c(NA,NA,NA,NA,                         #BRCA  
                      "high",NA,NA,NA,                     #CESC
                      NA,NA,NA,NA,                         #ESCA
                      "high","high",NA,"high",             #KIRC
                      NA,NA,NA,NA,                         #LIHC
                      NA,NA,NA,NA,                         #OV
                      "high",NA,NA,NA,                     #PAAD
                      NA,NA,NA,NA,                         #PRAD
                      NA,NA,NA,NA,                         #SKCM
                      "high",NA,"low","high",             #STAD
                      NA,NA,NA,NA,                         #THCA
                      "high",NA,NA,"high"                  #UCEC
), nrow = 4, ncol = 12)

gsvaqipaotu <- expand.grid(model = gsvamodeltype, cancertype = gsvacancer_type)
gsvaqipaotu$p_value <- as.vector(gsvap_value)
gsvaqipaotu$group <- as.vector(gsvagroup)
gsvaqipaotu <- gsvaqipaotu[!gsvaqipaotu$cancertype %in% c("BRCA","LIHC","OV","PRAD","THCA","SKCM","ESCA"),]
gsvaqipaotu$log10pvalue <- -log10(gsvaqipaotu$p_value)

ggplot(gsvaqipaotu, aes(x = model, y = cancertype)) +
  geom_point(aes(size = log10pvalue, color = group,fill=group), shape = 21) +
  scale_size_continuous(range = c(1, 10)) +  # 调整气泡大小范围
  scale_color_manual(values = c("low" = "steelblue", "high" = "#c11e21"),na.value = "white") +  # 设置分组颜色
  scale_fill_manual(values = c("low" = "steelblue", "high" = "#c11e21"), na.value = "white") + 
  labs(title = "Bubble Heatmap of p-values", x = "model", y = "CancerTypes", size = "-log10 pvalue", color = "Group") +
  theme_minimal()+
  theme_bw() + 
  theme(panel.grid=element_blank())


#集合图
library(UpSetR)
library(RColorBrewer)

# set1 <- c("CESC","KIRC","PAAD","STAD","UCEC")
# set2 <- c("KIRC")
# set3 <- c("CESC","PAAD","STAD")
# set4 <- c("ESCA","KIRC","SKCM")
# set5 <- c("STAD")
# set6 <- c("KIRC","SKCM")
# set7 <- c("STAD")
# 
# all_sets <- list(
#   model1 = set1,
#   model2 = set2,
#   model3 = set3,
#   model4 = set4,
#   model5 = set5,
#   model6 = set6,
#   model7 = set7
# )
# 
# pdf("D:/细胞状态2/论文/fig2/fig2e/fig2eupset.pdf")
# upset(fromList(all_sets),nsets=10, nintersects=100,
#       #sets=c("model1","model2","model3","model4","model5","model6","model7"),
#       number.angles = 0, #交互集合柱状图的柱标倾角
#       point.size=4, #图中点的大小
#       line.size=1, #图中连接线粗细
#       mainbar.y.label="Intersection size", #y轴的标签
#       main.bar.color = '#9D627E', #y轴柱状图颜色
#       matrix.color="black", #x轴点的颜色
#       sets.x.label="Set size",   #x轴的标签
#       sets.bar.color=c('#c11e21','#44bd28','#bab743','#5A84A7','#9280bb','#72c4d1',"#E48600"),
#       order.by = c("degree"),
#       decreasing = c(F),
#       keep.order = TRUE,
#       shade.color="#D2868C"
# )
# dev.off()
# 
# save(totalESTclin2,file="D:/细胞状态2/论文/fig2/totalESTclin2.RData")


set1 <- c("CESC","KIRC","PAAD","STAD","UCEC")
set2 <- c("KIRC")
set3 <- c("STAD")
set4 <- c("STAD","KIRC","UCEC")

all_sets <- list(
  TP = set1,
  SR = set2,
  GS = set3,
  FP = set4
)
"#ab2126""#009944""#009ccb""#9a0583"
pdf("D:/细胞状态2/论文/fig2/fig2e/newfig2eupset.pdf")
upset(fromList(all_sets),nsets=10, nintersects=100,
      #sets=c("model1","model2","model3","model4","model5","model6","model7"),
      number.angles = 0, #交互集合柱状图的柱标倾角
      point.size=4, #图中点的大小
      line.size=1, #图中连接线粗细
      mainbar.y.label="Intersection size", #y轴的标签
      main.bar.color = '#9D627E', #y轴柱状图颜色
      matrix.color="black", #x轴点的颜色
      sets.x.label="Set size",   #x轴的标签
      sets.bar.color=c("#ab2126","#9a0583","#009ccb","#009944"),
      order.by = c("degree"),
      decreasing = c(F),
      keep.order = TRUE,
      shade.color="#D2868C"
)
dev.off()

save(totalESTclin2,file="D:/细胞状态2/论文/fig2/newtotalESTclin2.RData")
