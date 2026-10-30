setwd("/wangzhenzenlab/bio_wzz_data/multiomics-data_bulk/expression-fpkm")

paadcount<-read.table("/wangzhenzenlab/bio_wzz_data/multiomics-data_bulk/expression-counts/TCGA-PAAD.htseq_counts.tsv",header=T,sep="\t")
paadcount[1:3,1:3]

BiocManager::install("AnnotationDbi")
BiocManager::install("DBI")
library(AnnotationDbi)
library(org.Hs.eg.db)
library(stringr)

paadcount<-read.table("D:/细胞状态2/bulk/expression-counts/TCGA-PAAD.htseq_counts.tsv",header=T,sep="\t")
class(paadcount)
paadcount$ensemble_id=unlist(str_split(paadcount$Ensembl_ID,"[.]",simplify=T))[,1]
entrez_ids  <-  select(org.Hs.eg.db, 
                       keys = paadcount$ensemble_id, 
                       columns = "SYMBOL", 
                       keytype = "ENSEMBL")
entrez_ids$ENSEMBL
which(duplicated(entrez_ids$ENSEMBL) == T)

paadcount2 <- merge(paadcount,entrez_ids,by.x="ensemble_id",by.y="ENSEMBL")

paadcount2[1:3,1:3]
paadcount2$ensemble_id
paadcount2$SYMBOL
paadcount3 <- paadcount2[!is.na(paadcount2$SYMBOL),]

paadcount3sign <- paadcount3[paadcount3$SYMBOL%in%paadsig_EMT_genes,]
rownames(paadcount3sign) <- paadcount3sign$SYMBOL
paadcount3sign$ensemble_id <- NULL
paadcount3sign$Ensembl_ID <- NULL
paadcount3sign$SYMBOL <- NULL
paadcount3signT <- t(paadcount3sign)
paadcount3signT <- as.data.frame(paadcount3signT)
head(paadcount3signT)
paadcount3signT$sample <- rownames(paadcount3signT)
paadcount3signT$patient <- sub("\\.[^.]*$", "", paadcount3signT$sample)

clin_all <- read.csv("D:/细胞状态2/bulk/clinical/clinical-all.csv")
clin_all$patient <- gsub("-", ".", clin_all$X_PATIENT)
clin_all$cancertype <- sub(".*\\((.*?)\\).*", "\\1", clin_all$X_cohort)
dim(clin_all[clin_all$cancertype=="PAAD",])

paadsur_sign <- merge(paadcount3signT,clin_all,by="patient")
paadsur_sign$X_OS
paadsur_sign$X_OS_IND

library(survival)
library(survminer)
paadsig_EMT_genes
paadsur_sign$COL1A1group <- ifelse(paadsur_sign$COL1A1>median(paadsur_sign$COL1A1),"high","low")
paadsur_sign$BGNgroup <- ifelse(paadsur_sign$BGN>median(paadsur_sign$BGN),"high","low")
paadsur_sign$TAGLNgroup <- ifelse(paadsur_sign$TAGLN>median(paadsur_sign$TAGLN),"high","low")
paadsur_sign$MYL9group <- ifelse(paadsur_sign$MYL9>median(paadsur_sign$MYL9),"high","low")
paadsur_sign$FN1group <- ifelse(paadsur_sign$FN1>median(paadsur_sign$FN1),"high","low")
paadsur_sign$ACTA2group <- ifelse(paadsur_sign$ACTA2>median(paadsur_sign$ACTA2),"high","low")
paadsur_sign$FSTL1group <- ifelse(paadsur_sign$FSTL1>median(paadsur_sign$FSTL1),"high","low")
paadsur_sign$TIMP1group <- ifelse(paadsur_sign$TIMP1>median(paadsur_sign$TIMP1),"high","low")
paadsur_sign$IGFBP4group <- ifelse(paadsur_sign$IGFBP4>median(paadsur_sign$IGFBP4),"high","low")
paadsur_sign$VIMgroup <- ifelse(paadsur_sign$VIM>median(paadsur_sign$VIM),"high","low")

surv_object <- survfit(Surv(X_OS, X_EVENT) ~ COL1A1group, data=paadsur_sign)
ggsurvplot(surv_object,
           pval = TRUE, conf.int = F,
           risk.table = F,
           palette=c("steelblue","#c11e21" ))

cox_model <- coxph(Surv(X_OS, X_OS_IND) ~ COL1A1 + BGN + TAGLN+MYL9+FN1+ACTA2+FSTL1+TIMP1+IGFBP4+VIM, data = paadsur_sign)
summary(cox_model)
ggforest(cox_model,
         data=paadsur_sign,
         main = "Hazard ratio",        # 标题
         #cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
         fontsize = 0.8, # 字体大小
         refLabel = "reference", #显示因子的参考水平
         noDigits = 3)
