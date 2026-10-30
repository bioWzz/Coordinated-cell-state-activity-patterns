library(ggvenn)
library(VennDiagram)
library(survival)
library("survminer")


load("D:/细胞状态2/论文/fig2/newtotalESTclin2.RData")
head(totalESTclin2)

#cescfold <- totalESTclin2[totalESTclin2$cancertype=="CESC",]
kircfold <- totalESTclin2[totalESTclin2$cancertype=="KIRC",]
#paadfold <- totalESTclin2[totalESTclin2$cancertype=="PAAD",]
stadfold <- totalESTclin2[totalESTclin2$cancertype=="STAD",]
ucecfold <- totalESTclin2[totalESTclin2$cancertype=="UCEC",]

c('#c11e21',"#E48600",'#bab743','#44bd28','#72c4d1','#5A84A7','#9280bb')

#KIRC
kircfold$groupTP <- ifelse(kircfold$TP> median(kircfold$TP),"high" ,"low" )
kircfold$groupSR <- ifelse(kircfold$SR> median(kircfold$SR),"high" ,"low" )
kircfold$groupFP <- ifelse(kircfold$FP> median(kircfold$FP),"high" ,"low" )

kircsetTP <- kircfold[kircfold$groupTP =="high","patient"]
kircsetSR <- kircfold[kircfold$groupSR =="high","patient"]
kircsetFP <- kircfold[kircfold$groupFP =="high","patient"]

kircset <- list(
  TP_Up=kircsetTP,
  SR_Up=kircsetSR,
  FP_Up=kircsetFP
)
ggvenn(kircset,
       fill_color=c('#ab2126',"#009944",'#9a0583'),
       fill_alpha = 0.4,
       stroke_color=NA
)


kircfold$GROUP <- ifelse(kircfold$groupTP=="high"& kircfold$groupSR=="high"& kircfold$groupFP=="high","T",
                         ifelse(kircfold$groupTP=="low"& kircfold$groupSR=="low"& kircfold$groupFP=="low","Z","O"))
surv_object <- survfit(Surv(X_OS, X_EVENT) ~ GROUP, data=kircfold)

ggsurvplot(surv_object,
           pval = TRUE, conf.int = F,
           risk.table = F,
           palette=c('#547c5d',"#ab2126","#4b6b86" )) 




#STAD
stadfold$groupTP <- ifelse(stadfold$TP> median(stadfold$TP),"high" ,"low" )
stadfold$groupGS <- ifelse(stadfold$GS> median(stadfold$GS),"high" ,"low" )
stadfold$groupFP <- ifelse(stadfold$FP> median(stadfold$FP),"high" ,"low" )

stadsetTP <- stadfold[stadfold$groupTP =="high","patient"]
stadsetGS <- stadfold[stadfold$groupGS =="low","patient"]
stadsetFP <- stadfold[stadfold$groupFP =="high","patient"]

stadset <- list(
  TP_Up=stadsetTP,
  GS_Down=stadsetGS,
  FP_Up=stadsetFP
)

ggvenn(stadset,
       fill_color=c('#ab2126','#009ccb','#9a0583'),
       fill_alpha = 0.4,
       stroke_color=NA
)

stadfold$GROUP <- ifelse(stadfold$groupTP=="high"&  stadfold$groupGS=="low"& stadfold$groupFP=="high","T",
                         ifelse(stadfold$groupTP=="low"& stadfold$groupGS=="high"& stadfold$groupFP=="low","Z","O"))
surv_object <- survfit(Surv(X_OS, X_EVENT) ~ GROUP, data=stadfold)

ggsurvplot(surv_object,
           pval = TRUE, conf.int = F,
           risk.table = F,
           palette=c('#547c5d',"#ab2126","#4b6b86"))



#UCEC
ucecfold$groupTP <- ifelse(ucecfold$TP> median(ucecfold$TP),"high" ,"low" )
ucecfold$groupFP <- ifelse(ucecfold$FP> median(ucecfold$FP),"high" ,"low" )
table(ucecfold$groupFP)

ucecsetTP <- ucecfold[ucecfold$groupTP =="high","patient"]
ucecsetFP <- ucecfold[ucecfold$groupFP =="high","patient"]

ucecset <- list(
  TP_Up=ucecsetTP,
  FP_Up=ucecsetFP
)

ggvenn(ucecset,
       fill_color=c('#ab2126','#9a0583'),
       fill_alpha = 0.4,
       stroke_color=NA
)

ucecfold$GROUP <- ifelse(ucecfold$groupTP=="high"&  ucecfold$groupFP=="high","T",
                         ifelse(ucecfold$groupTP=="low"& ucecfold$groupFP=="low","Z","O"))
surv_object <- survfit(Surv(X_OS, X_EVENT) ~ GROUP, data=ucecfold)

ggsurvplot(surv_object,
           pval = TRUE, conf.int = F,
           risk.table = F,
           palette=c('#547c5d',"#ab2126","#4b6b86"))

