library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
library(ggplot2)

load("D:/细胞状态2/论文/fig3/threecancersingle.RData")
head(SKCMsingle)
head(PAADsingle)
head(STADsingle)

cellstategene <- read.csv("D:/细胞状态2/cell_state.csv",header = T)
head(cellstategene)
table(cellstategene$Signature)

model1gene <- cellstategene[cellstategene$Signature%in%c("Angiogenesis","EMT","Invasion","Metastasis"),]
model5gene <- cellstategene[cellstategene$Signature%in%c("Cell Cycle","DNA damage","DNA repair"),]
model6gene <- cellstategene[cellstategene$Signature%in%c("Proliferation"),]

#stadFib
load("D:/细胞状态2/论文/分化/fenhua/zhushi/stad/Fib/cds.RData")
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=10)
max(Track_genes$morans_I)
min(Track_genes$morans_I)
median(Track_genes$morans_I)

Track_genes_sig <- Track_genes %>% top_n(n=100, morans_I) %>%  pull(gene_short_name) %>% as.character()

#plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="pseudotime", 
#                         min_expr=0.05, ncol = 2)

#DNA damage
model5int <- intersect(Track_genes_sig,model5gene$GeneName)
plot_genes_in_pseudotime(cds[model5int,], color_cells_by="pseudotime", 
                         min_expr=0.2, ncol = 2)

model5gene[model5gene$GeneName%in%model5int,]
Track_genes[model5int,]



#stadEnd
load("D:/细胞状态2/论文/分化/fenhua/zhushi/stad/End/cds.RData")
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=10)
max(Track_genes$morans_I)
min(Track_genes$morans_I)
median(Track_genes$morans_I)

Track_genes_sig <- Track_genes %>% top_n(n=100, morans_I) %>%  pull(gene_short_name) %>% as.character()

#EMT   FBLN5 COL6A2 CD44
#Invasion  COL6A2 LGALS1
#Metastasis CD44 CXCL12 EGR1 PDPN SELE
stadEndmodel1int <- intersect(Track_genes_sig,model1gene$GeneName)
plot_genes_in_pseudotime(cds[stadEndmodel1int,], color_cells_by="pseudotime", 
                         min_expr=0.1, ncol = 2)
model1gene[model1gene$GeneName%in%stadEndmodel1int,]
Track_genes[stadEndmodel1int,]


#DNA damage
stadEndmodel5int <- intersect(Track_genes_sig,model5gene$GeneName)
plot_genes_in_pseudotime(cds[stadEndmodel5int,], color_cells_by="pseudotime", 
                         min_expr=0.5, ncol = 2)
model5gene[model5gene$GeneName%in%stadEndmodel5int,]
Track_genes[stadEndmodel5int,]




#paadFib
load("D:/细胞状态2/论文/分化/fenhua/zhushi/paad/Fib/cds.RData")
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=10)
max(Track_genes$morans_I)
min(Track_genes$morans_I)
median(Track_genes$morans_I)

Track_genes_sig <- Track_genes %>% top_n(n=100, morans_I) %>%  pull(gene_short_name) %>% as.character()

#Angiogenesis     EGR3 FN1 NR4A1 MMP14 PDGFA VCAN SPP1
#EMT              BGN MMP2 ACTA2 COL1A1 VCAN FN1 MMP14 HTRA1 S100A4
#Invasion         AEBP1 BGN CALD1 COL1A1 VCAN CTSK FN1 INHBA MMP2 MMP11
#Metastasis       ANXA2 ATF3 CTSK EGR1 FN1 MMP11 MMP14 MMP2 PTP4A3 S100A4 SPP1 VCAN


paadFibmodel1int <- intersect(Track_genes_sig,model1gene$GeneName)
plot_genes_in_pseudotime(cds[paadFibmodel1int,], color_cells_by="pseudotime", 
                         min_expr=0.5, ncol = 2)
model1gene[model1gene$GeneName%in%paadFibmodel1int,]
Track_genes[paadFibmodel1int,]




#skcmFib
load("D:/细胞状态2/论文/分化/fenhua/zhushi/skcm/Fib/cds.RData")
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=10)
max(Track_genes$morans_I)
min(Track_genes$morans_I)
median(Track_genes$morans_I)

Track_genes_sig <- Track_genes %>% top_n(n=100, morans_I) %>%  pull(gene_short_name) %>% as.character()

#FGFR1 PTPRC MIF CCL5
skcmFibmodel6int <- intersect(Track_genes_sig,model6gene$GeneName)
class(colData(cds))
head(colData(cds))    # 查看 cds 的 colData
head(SKCMsingle)
table(colData(cds)$cell %in% SKCMsingle$Cell)
colData(cds) <- merge(colData(cds),SKCMsingle,by.x="cell",by.y="Cell")
plot_genes_in_pseudotime(cds[skcmFibmodel6int,], color_cells_by="pseudotime", 
                         min_expr=0.05, ncol = 2)
model6gene[model6gene$GeneName%in%skcmFibmodel6int,]
Track_genes[skcmFibmodel6int,]



