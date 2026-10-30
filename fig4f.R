install.packages("WGCNA")
install.packages("C:\\Users\\hy\\AppData\\Local\\Temp\\RtmpmuZWa7\\downloaded_packages\\magick_2.8.4.tar.gz", type = "source")
devtools::install_github("junjunlab/jjAnno")
devtools::install_github("junjunlab/ClusterGVis")
library(WGCNA)
library(monocle3)
library(dplyr)
library(ClusterGVis)

modulated_genes<-graph_test(cds,neighbor_graph="principal_graph",cores=4)
genes<-row.names(subset(modulated_genes,q_value==0&morans_I>0.25))
mat<-pre_pseudotime_matrix(cds_obj=cds,
                           gene_list=genes)
head(mat[1:5,1:5])
getwd()
setwd("D:/细胞状态2/论文/分化/fenhua/zhushi/test")
paadck<-clusterData(exp=mat,
                cluster.method="kmeans",
                cluster.num=2)
visCluster(object=paadck,
           plot.type="both",
           add.sampleanno=F,
           markGenes=sample(rownames(mat),replace=F))








stadmodulated_genes<-graph_test(stadcds,neighbor_graph="principal_graph",cores=4)
stadgenes<-row.names(subset(stadmodulated_genes,q_value==0&morans_I>0.25))
stadmat<-pre_pseudotime_matrix(cds_obj=stadcds,
                           gene_list=stadgenes)
head(stadmat[1:5,1:5])
stadck<-clusterData(exp=stadmat,
                    cluster.method="kmeans",
                    cluster.num=5)
pdf("stadvisCluster.pdf",onefile=F)
visCluster(object=stadck,
           plot.type="both",
           add.sampleanno=F,
           markGenes=sample(rownames(stadmat),30,replace=F))

dev.off()

save(modulated_genes,skcmmodulated_genes,stadmodulated_genes,file="D:/细胞状态2/论文/fig3/fig3f/modulated_genes.RData")







plot_cells(cds,
           color_cells_by = "group1",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=6)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)]%>%filter(q_value<1e-3)
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()

library(monocle3)
library(dplyr)
library(Seurat)
library(patchwork)
library(CytoTRACE)
library(stringr)


load()
colData(cds) <- merge(colData(cds),)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
ciliated_cds_pr_test_res <- graph_test(cds, 
                                       neighbor_graph="principal_graph", 
                                       cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))

gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(10^seq(-6,-1)))

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$group5)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")

plot_cells(cds,
           genes=gene_module_df %>% filter(module %in% c(15, 11, 12, 10,16,8,18)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)


modulegenes=gene_module_df %>% filter(module %in% c(15, 11, 12, 10,16,8,18))



Track_genes <- graph_test(cds, neighbor_graph="principal_graph")
# 提取前100个差异基因
Track_genes_sig <- Track_genes %>%
  top_n(n = 100, wt = morans_I) %>%
  pull(gene_short_name) %>%
  as.character()
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by = "group5", min_expr = 0.5)




gene_module_df[gene_module_df$id%in%genes,]
AFD_genes <- c("PDPN","GJA5","ACKR1")
genesAFD_genes <- c("PDPN","GJA5","ACKR1","PTPRC","SEMA3G" ,"CDKN1A" ,"RELN","CLU","FBLN5","TFF3" )
AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes,
                       colData(cds)$group1 %in% c("AFD")]
AFD_lineage_cds <- order_cells(AFD_lineage_cds)

plot_genes_in_pseudotime(cds["GJA5",],
                         color_cells_by="group5",
                         min_expr=0.5)

modulegenesid<- modulegenes$id
for (i in modulegenesid) {
  print(i)
  myfilename <- paste("D:/细胞状态2/论文/fig3/fig3f/stadendgene/",i,".pdf",sep = "")
  pdf(myfilename)
  p <- plot_genes_in_pseudotime(cds[i,],
                           color_cells_by="group1",
                           min_expr=0.5)
  print(p)
  dev.off()
}




#gene与model
load("D:/细胞状态2/论文/fig3/fig3d/stad/Endothelial/cds.RData")
cdsmodel1 <- cds@colData
class(cdsmodel1)
cdsmodel1 <- as.data.frame(cdsmodel1)
cdsmodel1[,1:8] <- NULL
head(cdsmodel1)

expression_data <-  exprs(cds)
expression_data[1:3,1:3]
expression_data <- t(expression_data)
class(expression_data)
expression_data <- as.data.frame(expression_data)


expression_data2 <- cbind(expression_data, cdsmodel1)

head(expression_data2[,14430:14435])
class(expression_data2[,1:14435])
model1correlations <- apply(expression_data2[,1:14429], 2, function(x) cor(x, expression_data2[,14430], use = "complete.obs"))
EMTcorrelations <- apply(expression_data2[,1:14429], 2, function(x) cor(x, expression_data2[,14431], use = "complete.obs"))
Angiogenesiscorrelations <- apply(expression_data2[,1:14429], 2, function(x) cor(x, expression_data2[,14432], use = "complete.obs"))
Invasioncorrelations <- apply(expression_data2[,1:14429], 2, function(x) cor(x, expression_data2[,14433], use = "complete.obs"))
Metastasiscorrelations <- apply(expression_data2[,1:14429], 2, function(x) cor(x, expression_data2[,14434], use = "complete.obs"))


model1high_cor_genes <- names(which(model1correlations > 0.2))
EMThigh_cor_genes <- names(which(EMTcorrelations > 0.2))
Angiogenesishigh_cor_genes <- names(which(Angiogenesiscorrelations > 0.2))
Invasionhigh_cor_genes <- names(which(Invasioncorrelations > 0.2))
Metastasishigh_cor_genes <- names(which(Metastasiscorrelations > 0.2))

model1highgene <- model1correlations[model1high_cor_genes] 
EMThighgene <-EMTcorrelations[EMThigh_cor_genes] 
Angiogenesishighgene <- Angiogenesiscorrelations[Angiogenesishigh_cor_genes] 
Invasionhighgene <- Invasioncorrelations[Invasionhigh_cor_genes] 
Metastasishighgene <- Metastasiscorrelations[Metastasishigh_cor_genes] 


Track_genes <- graph_test(cds, neighbor_graph="principal_graph")
# 提取前100个差异基因
Track_genes_sig <- Track_genes %>%
  top_n(n = 100, wt = morans_I) %>%
  pull(gene_short_name) %>%
  as.character()


model1testgene <- intersect(Track_genes_sig,names(model1highgene))
EMTtestgene <- intersect(Track_genes_sig,names(EMThighgene))
Angiogenesistestgene <- intersect(Track_genes_sig,names(Angiogenesishighgene))
Invasiontestgene <- intersect(Track_genes_sig,names(Invasionhighgene))
Metastasistestgene <- intersect(Track_genes_sig,names(Metastasishighgene))


model1sign <- model1highgene[model1testgene]
EMTsign <- EMThighgene[EMTtestgene]
Angiogenesissign <- Angiogenesishighgene[Angiogenesistestgene]
Invasionsign <- Invasionhighgene[Invasiontestgene]
Metastasissign <- Metastasishighgene[Metastasistestgene]

save(model1sign,EMTsign,Angiogenesissign,Invasionsign,Metastasissign,file="D:/细胞状态2/论文/fig3/fig3f/stadendgene/sign.RData")


colData(cds)
plot_genes_in_pseudotime(cds[model1testgene,],
                         color_cells_by ="group1",
                         min_expr=0.3)

a

load("D:/细胞状态2/论文/分化/fenhua/zhushi/paad/Fib/cds.RData")




load("D:/细胞状态2/论文/分化/fenhua/zhushi/stad/Fib/cds.RData")
cdsmodel5 <- cds@colData
class(cdsmodel5)
cdsmodel5 <- as.data.frame(cdsmodel5)
cdsmodel5[,1:8] <- NULL
head(cdsmodel5)

expression_data <-  exprs(cds)
expression_data[1:3,1:3]
expression_data <- t(expression_data)
class(expression_data)
expression_data <- as.data.frame(as.matrix(expression_data))


expression_data2 <- cbind(expression_data, cdsmodel5)
head(expression_data2[,14785:14789])
class(expression_data2[,1:14784])
model5correlations <- apply(expression_data2[,1:14784], 2, function(x) cor(x, expression_data2[,14785], use = "complete.obs"))
CellCyclecorrelations <- apply(expression_data2[,1:14784], 2, function(x) cor(x, expression_data2[,14786], use = "complete.obs"))
DNAdamagecorrelations <- apply(expression_data2[,1:14784], 2, function(x) cor(x, expression_data2[,14787], use = "complete.obs"))
DNArepaircorrelations <- apply(expression_data2[,1:14784], 2, function(x) cor(x, expression_data2[,14788], use = "complete.obs"))

model5high_cor_genes <- names(which(model5correlations > 0.2))
CellCyclehigh_cor_genes <- names(which(CellCyclecorrelations > 0.2))
DNAdamagehigh_cor_genes <- names(which(DNAdamagecorrelations > 0.2))
DNArepairhigh_cor_genes <- names(which(DNArepaircorrelations > 0.2))


model5highgene <- model5correlations[model5high_cor_genes] 
CellCyclehighgene <-CellCyclecorrelations[CellCyclehigh_cor_genes] 
DNAdamagehighgene <- DNAdamagecorrelations[DNAdamagehigh_cor_genes] 
DNArepairhighgene <- DNArepaircorrelations[DNArepairhigh_cor_genes] 



expression_data2$cell <- rownames(expression_data2)

#pseudotime <- cds@principal_graph_aux$UMAP$pseudotime
#pseudotime <- as.data.frame(pseudotime)
#pseudotime$cell <- rownames(pseudotime)
#pseudotime <- pseudotime[!is.infinite(pseudotime$pseudotime),]

#pseudotime2 <- merge(pseudotime,expression_data2,by="cell")
#rownames(pseudotime2) <- pseudotime2$cell
#pseudotime2[1:3,1:3]
#class(pseudotime2)
#pseudotime2$cell <- NULL
#pseudotime_column <- as.numeric(pseudotime2[,1])
#gene_expressions <- as.numeric(pseudotime2[,2:14771])
#correlations <- apply(pseudotime2[,2:14771], 2, function(x) cor(x, pseudotime_column, use = "complete.obs"))
#high_cor_genes <- names(which(correlations > 0.2))

#pseudotime2highgene <- correlations[high_cor_genes] 

#intersect(names(pseudotime2highgene),names(model5highgene))
#intersect(names(pseudotime2highgene),names(CellCyclehighgene))
#intersect(names(pseudotime2highgene),names(DNAdamagehighgene))
#intersect(names(pseudotime2highgene),names(DNArepairhighgene))

Track_genes <- graph_test(cds, neighbor_graph="principal_graph")
# 提取前100个差异基因
Track_genes_sig <- Track_genes %>%
  top_n(n = 100, wt = morans_I) %>%
  pull(gene_short_name) %>%
  as.character()

testgene <- intersect(Track_genes_sig,names(model5highgene))
testgene2 <- intersect(Track_genes_sig,names(CellCyclehighgene))
testgene3 <- intersect(Track_genes_sig,names(DNAdamagehighgene))
testgene4 <- intersect(Track_genes_sig,names(DNArepairhighgene))


model5sign <- model5highgene[testgene]
CellCyclesign <- CellCyclehighgene[testgene2]
DNAdamagesign <- DNAdamagehighgene[testgene3]
DNArepairsign <- DNArepairhighgene[testgene4]

save(model5sign,CellCyclesign,DNAdamagesign,DNArepairsign,file="D:/细胞状态2/论文/fig3/fig3f/stadfibgene/sign.RData")


colData(cds)
plot_genes_in_pseudotime(cds[testgene,],
                         color_cells_by ="group5",
                         min_expr=0.3)


?plot_genes_in_pseudotime
