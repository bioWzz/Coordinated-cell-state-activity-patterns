library(Seurat)
library(hdf5r)


load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/newthreecancersingle.RData")


#读取合并癌型
#CESC
cescexp1 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/CESC_GSE168652/CESC_GSE168652_expression.h5")

CESCexp <- cescexp1
#save(CESCexp,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3c/CESCexp.RData")

#ESCA
escaexp1 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/ESCA_GSE154763/ESCA_GSE154763_expression.h5")
escaexp2 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/ESCA_GSE160269/ESCA_GSE160269_expression.h5")
escaexp3 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/ESCA_GSE173950/ESCA_GSE173950_expression.h5")

escainter <- Reduce(intersect,list(rownames(escaexp1),rownames(escaexp2),rownames(escaexp3)))
escainter <- unique(escainter)
ESCAexp <- cbind(escaexp1[escainter,],escaexp2[escainter,],escaexp3[escainter,])

#KIRC
kircexp1 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/KIRC_GSE111360/KIRC_GSE111360_expression.h5")
kircexp2 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/KIRC_GSE121636/KIRC_GSE121636_expression.h5")
kircexp3 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/KIRC_GSE139555/KIRC_GSE139555_expression.h5")
kircexp4 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/KIRC_GSE145281_aPDL1/KIRC_GSE145281_aPDL1_expression.h5")
kircexp5 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/KIRC_GSE159115/KIRC_GSE159115_expression.h5")
kircexp6 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/KIRC_GSE171306/KIRC_GSE171306_expression.h5")


kircinter <- Reduce(intersect,list(rownames(kircexp1),rownames(kircexp2),rownames(kircexp3),rownames(kircexp4),rownames(kircexp5),rownames(kircexp6)))
kircinter <- unique(kircinter)
KIRCexp <- cbind(kircexp1[kircinter,],kircexp2[kircinter,],kircexp3[kircinter,],kircexp4[kircinter,],kircexp5[kircinter,],kircexp6[kircinter,])

#PAAD
paadexp1 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/PAAD_CRA001160/PAAD_CRA001160_expression.h5")
paadexp2 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/PAAD_GSE111672/PAAD_GSE111672_expression.h5")
paadexp3 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/PAAD_GSE141017/PAAD_GSE141017_expression.h5")
paadexp4 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/PAAD_GSE148673/PAAD_GSE148673_expression.h5")
paadexp5 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/PAAD_GSE154763/PAAD_GSE154763_expression.h5")
paadexp6 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/PAAD_GSE154778/PAAD_GSE154778_expression.h5")
paadexp7 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/PAAD_GSE158356/PAAD_GSE158356_expression.h5")
paadexp8 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/PAAD_GSE162708/PAAD_GSE162708_expression.h5")
paadexp9 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/PAAD_GSE165399/PAAD_GSE165399_expression.h5")


paadinter <- Reduce(intersect,list(rownames(paadexp1),rownames(paadexp2),rownames(paadexp3),rownames(paadexp4),rownames(paadexp5),rownames(paadexp6),rownames(paadexp7),rownames(paadexp8),rownames(paadexp9)))
paadinter <- unique(paadinter)
PAADexp <- cbind(paadexp1[paadinter,],paadexp2[paadinter,],paadexp3[paadinter,],paadexp4[paadinter,],paadexp5[paadinter,],paadexp6[paadinter,],paadexp7[paadinter,],paadexp8[paadinter,],paadexp9[paadinter,])

#SKCM
skcmexp1 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/SKCM_GSE72056/SKCM_GSE72056_expression.h5")
skcmexp2 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/SKCM_GSE115978_aPD1/SKCM_GSE115978_aPD1_expression.h5")
skcmexp3 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/SKCM_GSE120575_aPD1aCTLA4/SKCM_GSE120575_aPD1aCTLA4_expression.h5")
skcmexp4 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/SKCM_GSE123139/SKCM_GSE123139_expression.h5")
skcmexp5 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/SKCM_GSE134388_aPD1/SKCM_GSE134388_aPD1_expression.h5")
skcmexp6 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/SKCM_GSE139249/SKCM_GSE139249_expression.h5")
skcmexp7 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/SKCM_GSE148190/SKCM_GSE148190_expression.h5")
skcmexp8 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/SKCM_GSE159251/SKCM_GSE159251_expression.h5")
skcmexp9 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/SKCM_GSE166181_aPD1/SKCM_GSE166181_aPD1_expression.h5")
skcmexp10 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/SKCM_GSE179373/SKCM_GSE179373_expression.h5")


skcminter <- Reduce(intersect,list(rownames(skcmexp1),rownames(skcmexp2),rownames(skcmexp3),rownames(skcmexp4),rownames(skcmexp5),rownames(skcmexp6),rownames(skcmexp7),rownames(skcmexp8),rownames(skcmexp9),rownames(skcmexp10)))
skcminter <- unique(skcminter)
SKCMexp <- cbind(skcmexp1[skcminter,],skcmexp2[skcminter,],skcmexp3[skcminter,],skcmexp4[skcminter,],skcmexp5[skcminter,],skcmexp6[skcminter,],skcmexp7[skcminter,],skcmexp8[skcminter,],skcmexp9[skcminter,],skcmexp10[skcminter,])

#STAD
stadexp1 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/STAD_GSE134520/STAD_GSE134520_expression.h5")
stadexp2 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/STAD_GSE167297/STAD_GSE167297_expression.h5")


stadinter <- Reduce(intersect,list(rownames(stadexp1),rownames(stadexp2)))
stadinter <- unique(stadinter)
STADexp <- cbind(stadexp1[stadinter,],stadexp2[stadinter,])

#UCEC
ucecexp1 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/UCEC_GSE139555/UCEC_GSE139555_expression.h5")
ucecexp2 <- Read10X_h5("/wangzhenzenlab/bio_wzz_data/pancancer_scRNA-seq/gdh_185pancancer_scRNA-seq/TISCH2_185pancancer_scRNA-seq/UCEC_GSE154763/UCEC_GSE154763_expression.h5")


ucecinter <- Reduce(intersect,list(rownames(ucecexp1),rownames(ucecexp2)))
ucecinter <- unique(ucecinter)
UCECexp <- cbind(ucecexp1[ucecinter,],ucecexp2[ucecinter,])

#提出不同癌症不同细胞类型
#cesc
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/threecancersingle.RData")
head(CESCsingle)       #来自fig3a.R
table(CESCsingle$celltype)
cescCD8T <- CESCsingle[CESCsingle$celltype=="CD8T",]
cescMono <- CESCsingle[CESCsingle$celltype=="Mono/Macro",]
cescMali <- CESCsingle[CESCsingle$celltype=="Malignant",]

cescCD8Tcell <- cescCD8T$Cell
cescMonocell <- cescMono$Cell
cescMalicell <- cescMali$Cell
save(cescCD8Tcell,cescMonocell,cescMalicell,file = "D:/细胞状态2/论文/fig3/fig3c/cesccell.RData")
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3c/cesccell.RData")

CESCexpCD8T <- CESCexp[,cescCD8Tcell]
CESCexpMono <- CESCexp[,cescMonocell]
CESCexpMali <- CESCexp[,cescMalicell]

pbmc <- CreateSeuratObject(counts = CESCexpMali, project = "pbmc3k", min.cells = 2, min.features = 100)


cesccd8tpbmc <- pbmc
cescmonopbmc <- pbmc
cescmalipbmc <- pbmc

save(cesccd8tpbmc,cescmonopbmc,cescmalipbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/cesc/cescpbmc.RData")

#esca
head(ESCAsingle)
table(ESCAsingle$celltype)
escaB <- ESCAsingle[ESCAsingle$celltype=="B",]
escaCD4T <- ESCAsingle[ESCAsingle$celltype%in%c("CD4Tconv","Treg"),]
escaCD8T <- ESCAsingle[ESCAsingle$celltype%in%c("CD8T","CD8Tex","Tprolif"),]
escaDC <- ESCAsingle[ESCAsingle$celltype=="DC",]
escaMast <- ESCAsingle[ESCAsingle$celltype=="Mast",]
escaMono <- ESCAsingle[ESCAsingle$celltype=="Mono/Macro",]
escaEnd <- ESCAsingle[ESCAsingle$celltype=="Endothelial",]
escaEpi <- ESCAsingle[ESCAsingle$celltype=="Epithelial",]
escaFib <- ESCAsingle[ESCAsingle$celltype%in%c("Fibroblasts","Myofibroblasts"),]

escaBcell <- escaB$Cell
escaCD4Tcell <- escaCD4T$Cell
escaCD8Tcell <- escaCD8T$Cell
escaDCcell <- escaDC$Cell
escaMastcell <- escaMast$Cell
escaMonocell <- escaMono$Cell
escaEndcell <- escaEnd$Cell
escaEpicell <- escaEpi$Cell
escaFibcell <- escaFib$Cell


save(escaBcell,escaCD4Tcell,escaCD8Tcell,escaDCcell,escaMastcell,escaMonocell,
     escaEndcell,escaEpicell,escaFibcell,file = "D:/细胞状态2/论文/fig3/fig3c/escacell.RData")
load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3c/escacell.RData")

#标注longmap的需要去重
ESCAB <- ESCAexp[,escaBcell]
ESCACD4T <- ESCAexp[,escaCD4Tcell]
ESCACD8T <- ESCAexp[,escaCD8Tcell]
ESCADC <- ESCAexp[,escaDCcell]
ESCAMast <- ESCAexp[,escaMastcell]#logmap
ESCAMono <- ESCAexp[,escaMonocell]#logmap
ESCAEnd <- ESCAexp[,escaEndcell]
ESCAEpi <- ESCAexp[,escaEpicell]
ESCAFib <- ESCAexp[,escaFibcell]


if (any(duplicated(colnames(ESCAMono)))) {
  colnames(ESCAMono) <- make.unique(colnames(ESCAMono))
}

pbmc <- CreateSeuratObject(counts = ESCACD8T, project = "pbmc3k", min.cells = 2, min.features = 100)

escaBpbmc <- pbmc
escaCD4Tpbmc <- pbmc
escaCD8Tpbmc <- pbmc
escaDCpbmc <- pbmc
escaMastpbmc <- pbmc
escaMonopbmc <- pbmc
escaEndpbmc <- pbmc
escaEpipbmc <- pbmc
escaFibpbmc <- pbmc

save(escaBpbmc,escaCD4Tpbmc,escaCD8Tpbmc,escaDCpbmc,
     escaEndpbmc,escaEpipbmc,escaFibpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/esca/escapbmc.RData")
save(escaCD8Tpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/esca/escaCD8Tpbmc.RData")
save(escaMastpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/esca/escaMastpbmc.RData")
save(escaMonopbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/esca/escaMonopbmc.RData")
save(escaEndpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/esca/escaEndpbmc.RData")
save(escaEpipbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/esca/escaEpipbmc.RData")
save(escaFibpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/esca/escaFibpbmc.RData")



#kirc
head(KIRCsingle)
table(KIRCsingle$celltype)

kircB <- KIRCsingle[KIRCsingle$celltype=="B",]
kircCD4T <- KIRCsingle[KIRCsingle$celltype%in%c("CD4Tconv","Treg"),]
kircCD8T <- KIRCsingle[KIRCsingle$celltype%in%c("CD8T","CD8Tex","Tprolif"),]
kircDC <- KIRCsingle[KIRCsingle$celltype=="DC",]
kircMast <- KIRCsingle[KIRCsingle$celltype=="Mast",]
kircMono <- KIRCsingle[KIRCsingle$celltype=="Mono/Macro",]
kircNK <- KIRCsingle[KIRCsingle$celltype=="NK",]
kircEnd <- KIRCsingle[KIRCsingle$celltype=="Endothelial",]
kircEpi <- KIRCsingle[KIRCsingle$celltype%in%c("Epithelial","Malignant"),]
kircFib <- KIRCsingle[KIRCsingle$celltype%in%c("Fibroblasts","Myofibroblasts"),]


kircBcell <- kircB$Cell
kircCD4Tcell <- kircCD4T$Cell
kircCD8Tcell <- kircCD8T$Cell
kircDCcell <- kircDC$Cell
kircMastcell <- kircMast$Cell
kircMonocell <- kircMono$Cell
kircNKcell <- kircNK$Cell
kircEndcell <- kircEnd$Cell
kircEpicell <- kircEpi$Cell
kircFibcell <- kircFib$Cell


save(kircBcell,kircCD4Tcell,kircCD8Tcell,kircDCcell,kircMastcell,kircMonocell,
     kircNKcell,kircEndcell,kircEpicell,kircFibcell,file = "D:/细胞状态2/论文/fig3/fig3c/kirccell.RData")

load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3c/kirccell.RData")

KIRCB <- KIRCexp[,kircBcell]
KIRCCD4T <- KIRCexp[,kircCD4Tcell]
KIRCCD8T <- KIRCexp[,kircCD8Tcell]
KIRCDC <- KIRCexp[,kircDCcell]
KIRCMast <- KIRCexp[,kircMastcell]
KIRCMono <- KIRCexp[,kircMonocell]
KIRCNK <- KIRCexp[,kircNKcell]
KIRCEnd <- KIRCexp[,kircEndcell]
KIRCEpi <- KIRCexp[,kircEpicell]
KIRCFib <- KIRCexp[,kircFibcell]

pbmc <- CreateSeuratObject(counts = KIRCEpi, project = "pbmc3k", min.cells = 2, min.features = 100)

kircBpbmc <- pbmc
kircCD4Tpbmc <- pbmc
kircCD8Tpbmc <- pbmc
kircDCpbmc <- pbmc
kircMastpbmc <- pbmc
kircMonopbmc <- pbmc
kircNKpbmc <- pbmc
kircEndpbmc <- pbmc
kircEpipbmc <- pbmc
kircFibpbmc <- pbmc

save(kircBpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/kirc/kircBpbmc.RData")
save(kircCD4Tpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/kirc/kircCD4Tpbmc.RData")
save(kircCD8Tpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/kirc/kircCD8Tpbmc.RData")
save(kircDCpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/kirc/kircDCpbmc.RData")
save(kircMastpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/kirc/kircMastpbmc.RData")
save(kircMonopbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/kirc/kircMonopbmc.RData")
save(kircNKpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/kirc/kircNKpbmc.RData")
save(kircEndpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/kirc/kircEndpbmc.RData")
save(kircEpipbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/kirc/kircEpipbmc.RData")
save(kircFibpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/kirc/kircFibpbmc.RData")

#paad
head(PAADsingle)
table(PAADsingle$celltype)

paadB <- PAADsingle[PAADsingle$celltype=="B",]
#paadCD4T <- PAADsingle[PAADsingle$celltype%in%c("CD4Tconv","Treg"),]
paadCD8T <- PAADsingle[PAADsingle$celltype%in%c("CD8T","CD8Tex","Tprolif"),]
paadDC <- PAADsingle[PAADsingle$celltype=="DC",]
paadMast <- PAADsingle[PAADsingle$celltype=="Mast",]
paadMono <- PAADsingle[PAADsingle$celltype=="Mono/Macro",]
paadNK <- PAADsingle[PAADsingle$celltype=="NK",]
paadEnd <- PAADsingle[PAADsingle$celltype=="Endothelial",]
paadEpi <- PAADsingle[PAADsingle$celltype%in%c("Epithelial","Malignant"),]
paadFib <- PAADsingle[PAADsingle$celltype=="Fibroblasts",]
paadmyoFib <- PAADsingle[PAADsingle$celltype=="Myofibroblasts",]

paadBcell <- paadB$Cell
#paadCD4Tcell <- paadCD4T$Cell
paadCD8Tcell <- paadCD8T$Cell
paadDCcell <- paadDC$Cell
paadMastcell <- paadMast$Cell
paadMonocell <- paadMono$Cell
paadNKcell <- paadNK$Cell
paadEndcell <- paadEnd$Cell
paadEpicell <- paadEpi$Cell
paadFibcell <- paadFib$Cell
paadmyoFibcell <- paadmyoFib$Cell

save(paadBcell,paadCD8Tcell,paadDCcell,paadMastcell,paadMonocell,
     paadNKcell,paadEndcell,paadEpicell,paadFibcell,file = "D:/细胞状态2/论文/fig3/fig3c/paadcell.RData")

load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3c/paadcell.RData")

PAADB <- PAADexp[,paadBcell]
PAADCD8T <- PAADexp[,paadCD8Tcell]
PAADDC <- PAADexp[,paadDCcell]
PAADMast <- PAADexp[,paadMastcell]
PAADMono <- PAADexp[,paadMonocell]
PAADNK <- PAADexp[,paadNKcell]
PAADEnd <- PAADexp[,paadEndcell]
PAADEpi <- PAADexp[,paadEpicell]
PAADFib <- PAADexp[,paadFibcell]
PAADmyoFib <- PAADexp[,paadmyoFibcell]

pbmc <- CreateSeuratObject(counts = PAADmyoFib, project = "pbmc3k", min.cells = 2, min.features = 100)

paadBpbmc <- pbmc
paadCD8Tpbmc <- pbmc
paadDCpbmc <- pbmc
paadMastpbmc <- pbmc
paadMonopbmc <- pbmc
paadNKpbmc <- pbmc
paadEndpbmc <- pbmc
paadEpipbmc <- pbmc
paadFibpbmc <- pbmc
paadmyoFibpbmc <- pbmc

save(paadBpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/paad/paadBpbmc.RData")
save(paadCD8Tpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/paad/paadCD8Tpbmc.RData")
save(paadDCpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/paad/paadDCpbmc.RData")
save(paadMastpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/paad/paadMastpbmc.RData")
save(paadMonopbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/paad/paadMonopbmc.RData")
save(paadNKpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/paad/paadNKpbmc.RData")
save(paadEndpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/paad/paadEndpbmc.RData")
save(paadEpipbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/paad/paadEpipbmc.RData")
save(paadFibpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/paad/paadFibpbmc.RData")
save(paadmyoFibpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/paad/paadmyoFibpbmc.RData")
#save(paadFibpbmc,file = "/wangzhenzenlab/fsh/xff/fenhua/paad/paadFibpbmc.RData")
#save(paadmyoFibpbmc,file = "/wangzhenzenlab/fsh/xff/fenhua/paad/paadmyoFibpbmc.RData")

#skcm
head(SKCMsingle)
table(SKCMsingle$celltype)

skcmB <- SKCMsingle[SKCMsingle$celltype=="B",]
skcmCD4T <- SKCMsingle[SKCMsingle$celltype%in%c("CD4Tconv","Treg"),]
skcmCD8T <- SKCMsingle[SKCMsingle$celltype%in%c("CD8T","CD8Tex","Tprolif"),]
skcmDC <- SKCMsingle[SKCMsingle$celltype=="DC",]
skcmMast <- SKCMsingle[SKCMsingle$celltype=="Mast",]
skcmMono <- SKCMsingle[SKCMsingle$celltype=="Mono/Macro",]
skcmNK <- SKCMsingle[SKCMsingle$celltype=="NK",]
skcmMali <- SKCMsingle[SKCMsingle$celltype=="Malignant",]
skcmEnd <- SKCMsingle[SKCMsingle$celltype=="Endothelial",]
#skcmEpi <- SKCMsingle[SKCMsingle$celltype=="Epithelial",]
skcmFib <- SKCMsingle[SKCMsingle$celltype%in%c("Fibroblasts","Myofibroblasts"),]
skcmTprolif <- SKCMsingle[SKCMsingle$celltype%in%c("Tprolif"),]

skcmBcell <- skcmB$Cell
skcmCD4Tcell <- skcmCD4T$Cell
skcmCD8Tcell <- skcmCD8T$Cell
skcmDCcell <- skcmDC$Cell
skcmMastcell <- skcmMast$Cell
skcmMonocell <- skcmMono$Cell
skcmNKcell <- skcmNK$Cell
skcmMalicell <- skcmMali$Cell
skcmEndcell <- skcmEnd$Cell
#skcmEpicell <- skcmEpi$Cell
skcmFibcell <- skcmFib$Cell
skcmTprolifcell <- skcmTprolif$Cell

save(skcmBcell,skcmCD4Tcell,skcmCD8Tcell,skcmDCcell,skcmMastcell,skcmMonocell,
     skcmNKcell,skcmMalicell,skcmEndcell,skcmFibcell,file = "D:/细胞状态2/论文/fig3/fig3c/skcmcell.RData")

load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3c/skcmcell.RData")

SKCMB <- SKCMexp[,skcmBcell]#logmap
SKCMCD4T <- SKCMexp[,skcmCD4Tcell]#logmap
SKCMCD8T <- SKCMexp[,skcmCD8Tcell]#logmap
SKCMDC <- SKCMexp[,skcmDCcell]
SKCMMast <- SKCMexp[,skcmMastcell]
SKCMMono <- SKCMexp[,skcmMonocell]#logmap
SKCMNK <- SKCMexp[,skcmNKcell]#logmap
SKCMMali <- SKCMexp[,skcmMalicell]#logmap
SKCMEnd <- SKCMexp[,skcmEndcell]#logmap
SKCMFib <- SKCMexp[,skcmFibcell]#logmap
SKCMRprolif <- SKCMexp[,skcmTprolifcell]

#logmap的项进行如下操作
if (any(duplicated(colnames(SKCMCD4T)))) {
  colnames(SKCMCD4T) <- make.unique(colnames(SKCMCD4T))
}


pbmc <- CreateSeuratObject(counts = SKCMRprolif, project = "pbmc3k", min.cells = 2, min.features = 100)

skcmBpbmc <- pbmc
skcmCD4Tpbmc <- pbmc
skcmCD8Tpbmc <- pbmc
skcmDCpbmc <- pbmc
skcmMastpbmc <- pbmc
skcmMonopbmc <- pbmc
skcmNKpbmc <- pbmc
skcmMalipbmc <- pbmc
skcmEndpbmc <- pbmc
skcmFibpbmc <- pbmc
skcmTprolifpbmc <- pbmc

save(skcmBpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/skcm/skcmBpbmc.RData")
save(skcmCD4Tpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/skcm/skcmCD4Tpbmc.RData")
save(skcmCD8Tpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/skcm/skcmCD8Tpbmc.RData")
save(skcmDCpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/skcm/skcmDCpbmc.RData")
save(skcmMastpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/skcm/skcmMastpbmc.RData")
save(skcmMonopbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/skcm/skcmMonopbmc.RData")
save(skcmNKpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/skcm/skcmNKpbmc.RData")
save(skcmMalipbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/skcm/skcmMalipbmc.RData")
save(skcmEndpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/skcm/skcmEndpbmc.RData")
save(skcmFibpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/skcm/skcmFibpbmc.RData")
save(skcmTprolifpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig4/fig4d/skcm/skcmTprolifpbmc.RData")
save(skcmTprolifpbmc,file = "/wangzhenzenlab/fsh/xff/fenhua/zhushi/skcm/Tprolif/skcmTprolifpbmc.RData")

#stad
head(STADsingle)
table(STADsingle$celltype)

stadB <- STADsingle[STADsingle$celltype=="B",]
#stadCD4T <- STADsingle[STADsingle$celltype%in%c("CD4Tconv","Treg"),]
stadCD8T <- STADsingle[STADsingle$celltype%in%c("CD8T","CD8Tex","Tprolif"),]
stadDC <- STADsingle[STADsingle$celltype=="DC",]
stadMast <- STADsingle[STADsingle$celltype=="Mast",]
stadMono <- STADsingle[STADsingle$celltype=="Mono/Macro",]
#stadNK <- STADsingle[STADsingle$celltype=="NK",]
#stadMali <- STADsingle[STADsingle$celltype=="Malignant",]
stadEnd <- STADsingle[STADsingle$celltype=="Endothelial",]
stadEpi <- STADsingle[STADsingle$celltype=="Epithelial",]
stadFib <- STADsingle[STADsingle$celltype%in%c("Fibroblasts","Myofibroblasts"),]


stadBcell <- stadB$Cell
#stadCD4Tcell <- stadCD4T$Cell
stadCD8Tcell <- stadCD8T$Cell
stadDCcell <- stadDC$Cell
stadMastcell <- stadMast$Cell
stadMonocell <- stadMono$Cell
#stadNKcell <- stadNK$Cell
#stadMalicell <- stadMali$Cell
stadEndcell <- stadEnd$Cell
stadEpicell <- stadEpi$Cell
stadFibcell <- stadFib$Cell


save(stadBcell,stadCD8Tcell,stadDCcell,stadMastcell,stadMonocell,
     stadEndcell,stadEpicell,stadFibcell,file = "D:/细胞状态2/论文/fig3/fig3c/stadcell.RData")

load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3c/stadcell.RData")

STADB <- STADexp[,stadBcell]
STADCD8T <- STADexp[,stadCD8Tcell]
STADDC <- STADexp[,stadDCcell]
STADMast <- STADexp[,stadMastcell]
STADMono <- STADexp[,stadMonocell]
STADEnd <- STADexp[,stadEndcell]
STADEpi <- STADexp[,stadEpicell]
STADFib <- STADexp[,stadFibcell]


pbmc <- CreateSeuratObject(counts = STADEnd, project = "pbmc3k", min.cells = 2, min.features = 100)

stadBpbmc <- pbmc
stadCD8Tpbmc <- pbmc
stadDCpbmc <- pbmc
stadMastpbmc <- pbmc
stadMonopbmc <- pbmc
stadEndpbmc <- pbmc
stadEpipbmc <- pbmc
stadFibpbmc <- pbmc

save(stadBpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/stad/stadBpbmc.RData")
save(stadCD8Tpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/stad/stadCD8Tpbmc.RData")
save(stadDCpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/stad/stadDCpbmc.RData")
save(stadMastpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/stad/stadMastpbmc.RData")
save(stadMonopbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/stad/stadMonopbmc.RData")
save(stadEndpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/stad/stadEndpbmc.RData")
save(stadEpipbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/stad/stadEpipbmc.RData")
save(stadFibpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/stad/stadFibpbmc.RData")

#ucec
head(UCECsingle)
table(UCECsingle$celltype)

#ucecB <- UCECsingle[UCECsingle$celltype=="B",]
ucecCD4T <- UCECsingle[UCECsingle$celltype%in%c("CD4Tconv","Treg"),]
ucecCD8T <- UCECsingle[UCECsingle$celltype%in%c("CD8T","CD8Tex","Tprolif"),]
ucecDC <- UCECsingle[UCECsingle$celltype=="DC",]
ucecMast <- UCECsingle[UCECsingle$celltype=="Mast",]
ucecMono <- UCECsingle[UCECsingle$celltype=="Mono/Macro",]
#ucecNK <- UCECsingle[UCECsingle$celltype=="NK",]
#ucecMali <- UCECsingle[UCECsingle$celltype=="Malignant",]
#ucecEnd <- UCECsingle[UCECsingle$celltype=="Endothelial",]
#ucecEpi <- UCECsingle[UCECsingle$celltype=="Epithelial",]
ucecFib <- UCECsingle[UCECsingle$celltype%in%c("Fibroblasts","Myofibroblasts"),]


#ucecBcell <- ucecB$Cell
ucecCD4Tcell <- ucecCD4T$Cell
ucecCD8Tcell <- ucecCD8T$Cell
ucecDCcell <- ucecDC$Cell
ucecMastcell <- ucecMast$Cell
ucecMonocell <- ucecMono$Cell
#ucecNKcell <- ucecNK$Cell
#ucecMalicell <- ucecMali$Cell
#ucecEndcell <- ucecEnd$Cell
#ucecEpicell <- ucecEpi$Cell
ucecFibcell <- ucecFib$Cell


save(ucecCD4Tcell,ucecCD8Tcell,ucecDCcell,ucecMastcell,ucecMonocell,ucecFibcell,
     file = "D:/细胞状态2/论文/fig3/fig3c/uceccell.RData")

load("/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3c/uceccell.RData")


UCECCD4T <- UCECexp[,ucecCD4Tcell]
UCECCD8T <- UCECexp[,ucecCD8Tcell]
UCECDC <- UCECexp[,ucecDCcell]
UCECMast <- UCECexp[,ucecMastcell]
UCECMono <- UCECexp[,ucecMonocell]
UCECFib <- UCECexp[,ucecFibcell]

pbmc <- CreateSeuratObject(counts = UCECFib, project = "pbmc3k", min.cells = 2, min.features = 100)

ucecCD4Tpbmc <- pbmc
ucecCD8Tpbmc <- pbmc
ucecDCpbmc <- pbmc
ucecMastpbmc <- pbmc
ucecMonopbmc <- pbmc
ucecFibpbmc <- pbmc

save(ucecCD4Tpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/ucec/ucecCD4Tpbmc.RData")
save(ucecCD8Tpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/ucec/ucecCD8Tpbmc.RData")
save(ucecDCpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/ucec/ucecDCpbmc.RData")
save(ucecMastpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/ucec/ucecMastpbmc.RData")
save(ucecMonopbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/ucec/ucecMonopbmc.RData")
save(ucecFibpbmc,file = "/wangzhenzenlab/xff/cellstate/xff/cellstate2/fig3/fig3d/ucec/ucecFibpbmc.RData")

#seurat
pbmc <- seurat_object
#QC质量控制来筛选高质量的细胞
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#QC指标过滤细胞
pbmc <- subset(pbmc, subset = nFeature_RNA > 100 & nFeature_RNA < 2000 & percent.mt < 5)
#数据标准化
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
#识别高变基因
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#数据归一化
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
#线性降维PCA主成分分析
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

#用KNN方法对细胞聚类
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
#非线性降维(UMAP/tSNE)
#UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

UMAPPlot(pbmc)


