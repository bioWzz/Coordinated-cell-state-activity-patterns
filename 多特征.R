load("D:/细胞状态2/论文/挖掘/独立验证集/BRCA/newbrcametadata.RData")
load("D:/细胞状态2/论文/挖掘/独立验证集/KIRC/newkircmetadata.RData")
load("D:/细胞状态2/论文/挖掘/独立验证集/PAAD/newpaadmetadata.RData")
load("D:/细胞状态2/论文/挖掘/独立验证集/SKCM/newskcmmetadata.RData")
load("D:/细胞状态2/论文/挖掘/独立验证集/炎症/newstadmetadata.RData")

head(stmetadata3_0733)

stmetadata3_0733$cancertype <- "STAD"
paadstmetadata2_2005$cancertype <- "PAAD"
brcastmetadata7599_2$cancertype <- "BRCA"
skcmstmetadata7625_2$cancertype <- "SKCM"
kircstmetadata4030_2$cancertype <- "KIRC"


figBstmetadata <- rbind(stmetadata3_0733[,c("level","Inflammation","cancertype")],paadstmetadata2_2005[,c("level","Inflammation","cancertype")],
                        brcastmetadata7599_2[,c("level","Inflammation","cancertype")],skcmstmetadata7625_2[,c("level","Inflammation","cancertype")],
                        kircstmetadata4030_2[,c("level","Inflammation","cancertype")])
head(figBstmetadata)

library(ggplot2)
library(dplyr)

lines_data <- figBstmetadata %>%
  group_by(cancertype) %>%
  summarise(y_value = median(Inflammation)) %>%
  ungroup()


ggplot(figBstmetadata, aes(x = level, y = Inflammation, fill = level)) +
  geom_violin(trim = FALSE) +  # 添加小提琴图层，设置trim = FALSE以显示所有数据点
  scale_fill_brewer(palette = "Paired") +  # 使用颜色方案
  facet_wrap(~cancertype, ncol = 1)+
  scale_color_manual(values = c("steelblue", "#99bca0", "#e99593"))+
  scale_fill_manual(values = c("steelblue", "#99bca0", "#e99593"))+
  theme_bw() + 
  theme(panel.grid=element_blank())+
  geom_hline(data = lines_data, aes(yintercept = y_value), linetype = "dashed", color = "black") 







load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/st/stad/stmetadata5.RData")
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/挖掘/特征/paad/paadalldata.RData")
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/挖掘/特征/BRCA/brcaalldata.RData")
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/挖掘/特征/KIRC/kircalldata.RData")
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/挖掘/特征/SKCM/skcmalldata.RData")


daixieplotfun2 <- function(alldata,numb){
  df <- alldata[c("level", names(alldata)[which(names(alldata) == "Glycolysis / Gluconeogenesis"):ncol(alldata)])]
  colnames(df) <- gsub("/", "-", colnames(df))
  for (col in colnames(df)[2:ncol(df)]) {
    # 创建 PDF 文件
    pdf(paste0(numb,"/",col, ".pdf"))
    
    # 绘制小提琴图
    p <- ggplot(data = df, aes(x = level, y = get(col), fill = level)) +
      geom_violin(position = position_dodge(0.8)) +
      geom_boxplot(width = 0.5, position = position_dodge(0.8)) +
      stat_compare_means(aes(group = level)) +
      scale_fill_manual(values = c("steelblue", "#99bca0", "#e99593")) +
      theme(axis.text.x = element_text(size = 20, angle = 30, hjust = 0.5, vjust = 0.5))
    # 关闭 PDF 设备
    print(p)
    dev.off()
  }
}

scoreplotfun <- function(alldata,numb){
  df <- alldata[,c("level","TumorPurity","ImmuneScore")]
  pdf(paste0(numb,"TumorPurity.pdf"))
  p1 <- ggplot(data = df, aes(x = level, y = TumorPurity, fill = level)) +
    geom_violin(position = position_dodge(0.8)) +
    geom_boxplot(width = 0.5, position = position_dodge(0.8)) +
    stat_compare_means(aes(group = level)) +
    scale_fill_manual(values = c("steelblue", "#99bca0", "#e99593")) +
    theme(axis.text.x = element_text(size = 20, angle = 30, hjust = 0.5, vjust = 0.5))
  print(p1)
  dev.off()
  
  pdf(paste0(numb,"ImmuneScore.pdf"))
  p2 <- ggplot(data = df, aes(x = level, y = ImmuneScore, fill = level)) +
    geom_violin(position = position_dodge(0.8)) +
    geom_boxplot(width = 0.5, position = position_dodge(0.8)) +
    stat_compare_means(aes(group = level)) +
    scale_fill_manual(values = c("steelblue", "#99bca0", "#e99593")) +
    theme(axis.text.x = element_text(size = 20, angle = 30, hjust = 0.5, vjust = 0.5))
  print(p2)
  dev.off()
}

stmetadata5_0331[1:5,1:40]

library(ggplot2)
library(ggpubr)

#BRCA
setwd("/wangzhenzenlab/xff/cellstate/xff/cellstate2/挖掘/特征/BRCA/代谢/")

daixieplotfun2(alldata7599, "7599")
daixieplotfun2(alldata7601, "7601")
daixieplotfun2(alldata7603, "7603")

setwd("/wangzhenzenlab/xff/cellstate/xff/cellstate2/挖掘/特征/BRCA/estimate/")
scoreplotfun(alldata7599, "7599")
scoreplotfun(alldata7601, "7601")
scoreplotfun(alldata7603, "7603")


#PAAD
setwd("/wangzhenzenlab/xff/cellstate/xff/cellstate2/挖掘/特征/paad/代谢/")

daixieplotfun2(alldata1905, "1905")
daixieplotfun2(alldata1923, "1923")
daixieplotfun2(alldata1945, "1945")
daixieplotfun2(alldata1960, "1960")
daixieplotfun2(alldata2002, "2002")
daixieplotfun2(alldata2005, "2005")
daixieplotfun2(alldata2021, "2021")

setwd("/wangzhenzenlab/xff/cellstate/xff/cellstate2/挖掘/特征/paad/estimate/")
scoreplotfun(alldata1905, "1905")
scoreplotfun(alldata1923, "1923")
scoreplotfun(alldata1945, "1945")
scoreplotfun(alldata1960, "1960")
scoreplotfun(alldata2002, "2002")
scoreplotfun(alldata2005, "2005")
scoreplotfun(alldata2021, "2021")

#KIRC
setwd("/wangzhenzenlab/xff/cellstate/xff/cellstate2/挖掘/特征/KIRC/代谢/")

daixieplotfun2(alldata4030, "4030")
daixieplotfun2(alldata4031, "4031")
daixieplotfun2(alldata4032, "4032")

setwd("/wangzhenzenlab/xff/cellstate/xff/cellstate2/挖掘/特征/KIRC/estimate/")
scoreplotfun(alldata4030, "4030")
scoreplotfun(alldata4031, "4031")
scoreplotfun(alldata4032, "4032")

#SKCM
setwd("/wangzhenzenlab/xff/cellstate/xff/cellstate2/挖掘/特征/SKCM/代谢/")

daixieplotfun2(alldata7624, "7624")
daixieplotfun2(alldata7625, "7625")

setwd("/wangzhenzenlab/xff/cellstate/xff/cellstate2/挖掘/特征/SKCM/estimate/")
scoreplotfun(alldata7624, "7624")
scoreplotfun(alldata7625, "7625")


stmetadata5_0733$cancertype <- "STAD"
alldata2005$cancertype <- "PAAD"
alldata7599$cancertype <- "BRCA"
alldata7625$cancertype <- "SKCM"
alldata4030$cancertype <- "KIRC"
stmetadata5_0733$TumorPurity
stmetadata5_0733$'Glycolysis / Gluconeogenesis'

load("/wangzhenzenlab/xff/cellstate/xff/炎症/estimate/scores.RData")
load("/wangzhenzenlab/xff/cellstate/xff/炎症/estimate/scores0732.RData")
dim(scores0733T)
dim(alldata0733)
scores0733T$CELL <- rownames(scores0733T)
scores0733T$cell <- gsub("\\.", "-", scores0733T$CELL)
alldata0733 <- merge(stmetadata5_0733,scores0733T,by="cell")
head(alldata0733)

figBstmetadata2 <- rbind(alldata0733[,c("level","Inflammation","cancertype","TumorPurity","ImmuneScore",
                                             "exclude","exhaustion","senility","Tcell.dysfunction",
                                             "Glycolysis / Gluconeogenesis","model1")],
                         alldata2005[,c("level","Inflammation","cancertype","TumorPurity","ImmuneScore",
                                        "exclude","exhaustion","senility","Tcell.dysfunction",
                                        "Glycolysis / Gluconeogenesis","model1")],
                        alldata7599[,c("level","Inflammation","cancertype","TumorPurity","ImmuneScore",
                                       "exclude","exhaustion","senility","Tcell.dysfunction",
                                       "Glycolysis / Gluconeogenesis","model1")],
                        alldata7625[,c("level","Inflammation","cancertype","TumorPurity","ImmuneScore",
                                       "exclude","exhaustion","senility","Tcell.dysfunction",
                                       "Glycolysis / Gluconeogenesis","model1")],
                        alldata4030[,c("level","Inflammation","cancertype","TumorPurity","ImmuneScore",
                                       "exclude","exhaustion","senility","Tcell.dysfunction",
                                       "Glycolysis / Gluconeogenesis","model1")])

head(figBstmetadata2)

lines_data <- figBstmetadata2 %>%
  group_by(cancertype) %>%
  summarise(y_value = median(Inflammation)) %>%
  ungroup()

pdf("/wangzhenzenlab/xff/cellstate/xff/cellstate2/挖掘/特征/Inflammation.pdf",width = 3,height = 7)
ggplot(figBstmetadata2, aes(x = level, y = Inflammation, fill = level)) +
  geom_violin(trim = FALSE) +  # 添加小提琴图层，设置trim = FALSE以显示所有数据点
  scale_fill_brewer(palette = "Paired") +  # 使用颜色方案
  facet_wrap(~cancertype, ncol = 1)+
  scale_color_manual(values = c("steelblue", "#99bca0", "#e99593"))+
  scale_fill_manual(values = c("steelblue", "#99bca0", "#e99593"))+
  theme_bw() + 
  theme(panel.grid=element_blank())+
  geom_hline(data = lines_data, aes(yintercept = y_value), linetype = "dashed", color = "black") 
dev.off()


lines_data <- figBstmetadata2 %>%
  group_by(cancertype) %>%
  summarise(y_value = median(TumorPurity)) %>%
  ungroup()
pdf("/wangzhenzenlab/xff/cellstate/xff/cellstate2/挖掘/特征/TumorPurity.pdf",width = 3,height = 7)
ggplot(figBstmetadata2, aes(x = level, y = TumorPurity, fill = level)) +
  geom_violin(trim = FALSE) +  # 添加小提琴图层，设置trim = FALSE以显示所有数据点
  scale_fill_brewer(palette = "Paired") +  # 使用颜色方案
  facet_wrap(~cancertype, ncol = 1)+
  scale_color_manual(values = c("steelblue", "#99bca0", "#e99593"))+
  scale_fill_manual(values = c("steelblue", "#99bca0", "#e99593"))+
  theme_bw() + 
  theme(panel.grid=element_blank())+
  geom_hline(data = lines_data, aes(yintercept = y_value), linetype = "dashed", color = "black") 
dev.off()


lines_data <- figBstmetadata2 %>%
  group_by(cancertype) %>%
  summarise(y_value = median(ImmuneScore)) %>%
  ungroup()
pdf("/wangzhenzenlab/xff/cellstate/xff/cellstate2/挖掘/特征/ImmuneScore.pdf",width = 3,height = 7)
ggplot(figBstmetadata2, aes(x = level, y = ImmuneScore, fill = level)) +
  geom_violin(trim = FALSE) +  # 添加小提琴图层，设置trim = FALSE以显示所有数据点
  scale_fill_brewer(palette = "Paired") +  # 使用颜色方案
  facet_wrap(~cancertype, ncol = 1)+
  scale_color_manual(values = c("steelblue", "#99bca0", "#e99593"))+
  scale_fill_manual(values = c("steelblue", "#99bca0", "#e99593"))+
  theme_bw() + 
  theme(panel.grid=element_blank())+
  geom_hline(data = lines_data, aes(yintercept = y_value), linetype = "dashed", color = "black") 
dev.off()



lines_data <- figBstmetadata2 %>%
  group_by(cancertype) %>%
  summarise(y_value = median(exclude)) %>%
  ungroup()
pdf("/wangzhenzenlab/xff/cellstate/xff/cellstate2/挖掘/特征/exclude.pdf",width = 3,height = 7)
ggplot(figBstmetadata2, aes(x = level, y = exclude, fill = level)) +
  geom_violin(trim = FALSE) +  # 添加小提琴图层，设置trim = FALSE以显示所有数据点
  scale_fill_brewer(palette = "Paired") +  # 使用颜色方案
  facet_wrap(~cancertype, ncol = 1)+
  scale_color_manual(values = c("steelblue", "#99bca0", "#e99593"))+
  scale_fill_manual(values = c("steelblue", "#99bca0", "#e99593"))+
  theme_bw() + 
  theme(panel.grid=element_blank())+
  geom_hline(data = lines_data, aes(yintercept = y_value), linetype = "dashed", color = "black") 
dev.off()


lines_data <- figBstmetadata2 %>%
  group_by(cancertype) %>%
  summarise(y_value = median(exhaustion)) %>%
  ungroup()
pdf("/wangzhenzenlab/xff/cellstate/xff/cellstate2/挖掘/特征/exhaustion.pdf",width = 3,height = 7)
ggplot(figBstmetadata2, aes(x = level, y = exhaustion, fill = level)) +
  geom_violin(trim = FALSE) +  # 添加小提琴图层，设置trim = FALSE以显示所有数据点
  scale_fill_brewer(palette = "Paired") +  # 使用颜色方案
  facet_wrap(~cancertype, ncol = 1)+
  scale_color_manual(values = c("steelblue", "#99bca0", "#e99593"))+
  scale_fill_manual(values = c("steelblue", "#99bca0", "#e99593"))+
  theme_bw() + 
  theme(panel.grid=element_blank())+
  geom_hline(data = lines_data, aes(yintercept = y_value), linetype = "dashed", color = "black") 
dev.off()



lines_data <- figBstmetadata2 %>%
  group_by(cancertype) %>%
  summarise(y_value = median(senility)) %>%
  ungroup()
pdf("/wangzhenzenlab/xff/cellstate/xff/cellstate2/挖掘/特征/senility.pdf",width = 3,height = 7)
ggplot(figBstmetadata2, aes(x = level, y = senility, fill = level)) +
  geom_violin(trim = FALSE) +  # 添加小提琴图层，设置trim = FALSE以显示所有数据点
  scale_fill_brewer(palette = "Paired") +  # 使用颜色方案
  facet_wrap(~cancertype, ncol = 1)+
  scale_color_manual(values = c("steelblue", "#99bca0", "#e99593"))+
  scale_fill_manual(values = c("steelblue", "#99bca0", "#e99593"))+
  theme_bw() + 
  theme(panel.grid=element_blank())+
  geom_hline(data = lines_data, aes(yintercept = y_value), linetype = "dashed", color = "black") 
dev.off()


lines_data <- figBstmetadata2 %>%
  group_by(cancertype) %>%
  summarise(y_value = median(Tcell.dysfunction)) %>%
  ungroup()
pdf("/wangzhenzenlab/xff/cellstate/xff/cellstate2/挖掘/特征/Tcell.dysfunction.pdf",width = 3,height = 7)
ggplot(figBstmetadata2, aes(x = level, y = Tcell.dysfunction, fill = level)) +
  geom_violin(trim = FALSE) +  # 添加小提琴图层，设置trim = FALSE以显示所有数据点
  scale_fill_brewer(palette = "Paired") +  # 使用颜色方案
  facet_wrap(~cancertype, ncol = 1)+
  scale_color_manual(values = c("steelblue", "#99bca0", "#e99593"))+
  scale_fill_manual(values = c("steelblue", "#99bca0", "#e99593"))+
  theme_bw() + 
  theme(panel.grid=element_blank())+
  geom_hline(data = lines_data, aes(yintercept = y_value), linetype = "dashed", color = "black") 
dev.off()



# load("D:/细胞状态2/论文/挖掘/独立验证集/BRCA/代谢/BRCAdaixie.RData")
# metabolism.matrixT7599$`Glycolysis / Gluconeogenesis`

lines_data <- figBstmetadata2 %>%
  group_by(cancertype) %>%
  summarise(y_value = median(`Glycolysis / Gluconeogenesis`)) %>%
  ungroup()
pdf("/wangzhenzenlab/xff/cellstate/xff/cellstate2/挖掘/特征/Glycolysis-Gluconeogenesis.pdf",width = 3,height = 7)
ggplot(figBstmetadata2, aes(x = level, y = `Glycolysis / Gluconeogenesis`, fill = level)) +
  geom_violin(trim = FALSE) +  # 添加小提琴图层，设置trim = FALSE以显示所有数据点
  scale_fill_brewer(palette = "Paired") +  # 使用颜色方案
  facet_wrap(~cancertype, ncol = 1)+
  scale_color_manual(values = c("steelblue", "#99bca0", "#e99593"))+
  scale_fill_manual(values = c("steelblue", "#99bca0", "#e99593"))+
  theme_bw() + 
  theme(panel.grid=element_blank())+
  geom_hline(data = lines_data, aes(yintercept = y_value), linetype = "dashed", color = "black") 
dev.off()





figBstmetadata3 <- figBstmetadata2[,c("Inflammation","TumorPurity","ImmuneScore",
                                      "exclude","exhaustion","senility","Tcell.dysfunction",
                                      "Glycolysis / Gluconeogenesis","model1")]

head(figBstmetadata3)


singlemt <- figBstmetadata3
singlemt <- scale(singlemt)
cor_pearson <- cor(singlemt, method = 'pearson')
cor_pearson
cordataframe <- as.data.frame(cor_pearson)
write.csv(cordataframe,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/挖掘/特征/cor.csv")

library(corrplot)
my_colors <- colorRampPalette(c("steelblue", "white", "#c11e21"))(100)

pdf("/wangzhenzenlab/xff/cellstate/xff/cellstate2/挖掘/特征/cor.pdf")
corrplot(cor_pearson, method = 'pie', type = 'full', order = 'hclust', add = F, diag = FALSE,col = my_colors)
dev.off()



