library(ggplot2)
library(jjPlot)
library(reshape2)
library(ggnewscale)

head(PAADsingle)
head(STADsingle)
head(UCECsingle)
table(PAADsingle$celltype)

#PAAD model1 Fibroblasts Myofibroblasts
paadmodel1 <- PAADsingle[PAADsingle$celltype%in%c("Fibroblasts","Myofibroblasts"),]

result <- paadmodel1 %>%
  group_by(celltype) %>%
  summarise(Angiogenesis = median(Angiogenesis),
            EMT = median(EMT),
            Invasion = median(Invasion),
            Metastasis = median(Metastasis)
  )

ggplot(result,
       aes(x = "",y = celltype)) +
  geom_jjtriangle(aes(fill = Angiogenesis),type = 'mu') +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "#c11e21",
                       midpoint = 0, limit = c(-0.5,0.5)) +
  # new legend
  new_scale_fill() +
  geom_jjtriangle(aes(fill = EMT),type = 'mb') +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "#c11e21",
                       midpoint = 0, limit = c(-0.5,0.5)) +
  # new legend
  new_scale_fill() +
  geom_jjtriangle(aes(fill = Invasion),type = 'ml') +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "#c11e21",
                       midpoint = 0, limit = c(-0.5,0.5)) +
  # new legend
  new_scale_fill() +
  geom_jjtriangle(aes(fill = Metastasis),type = 'mr') +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "#c11e21",
                       midpoint = 0, limit = c(-0.5,0.5)) +
  theme(legend.position = 'top') +
  coord_fixed()+
  theme_bw() + 
  theme(panel.grid=element_blank())


#stad model1 Mast Mono/Macro Endothelial Fibroblasts
stadmodel1 <- STADsingle[STADsingle$celltype%in%c("Mast","NK","Endothelial","Fibroblasts"),]
result4 <- stadmodel1 %>%
  group_by(celltype) %>%
  summarise(Angiogenesis = median(Angiogenesis),
            EMT = median(EMT),
            Invasion = median(Invasion),
            Metastasis = median(Metastasis)
  )
ggplot(result4,
       aes(x = "",y = celltype)) +
  geom_jjtriangle(aes(fill = Angiogenesis),type = 'mu') +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "#c11e21",
                       midpoint = 0, limit = c(-0.5,0.5)) +
  # new legend
  new_scale_fill() +
  geom_jjtriangle(aes(fill = EMT),type = 'mb') +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "#c11e21",
                       midpoint = 0, limit = c(-0.5,0.5)) +
  # new legend
  new_scale_fill() +
  geom_jjtriangle(aes(fill = Invasion),type = 'ml') +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "#c11e21",
                       midpoint = 0, limit = c(-0.5,0.5)) +
  # new legend
  new_scale_fill() +
  geom_jjtriangle(aes(fill = Metastasis),type = 'mr') +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "#c11e21",
                       midpoint = 0, limit = c(-0.5,0.5)) +
  theme(legend.position = 'top') +
  coord_fixed()+
  theme_bw() + 
  theme(panel.grid=element_blank())

#stad model5 DC Mast Mono/Macro Endothelial Fibroblasts
table(STADsingle$celltype)
stadmodel5 <- STADsingle[STADsingle$celltype%in%c("Mono/Macro","Endothelial","Fibroblasts"),]
result5 <- stadmodel5 %>%
  group_by(celltype) %>%
  summarise(`Cell Cycle` = median(`Cell Cycle`),
            `DNA damage` = median(`DNA damage`),
            `DNA repair` = median(`DNA repair`)
  )
ggplot(result5,
       aes(x = "",y = celltype)) +
  geom_jjtriangle(aes(fill = `Cell Cycle`),type = 'mu') +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "#c11e21",
                       midpoint = 0, limit = c(-0.5,0.5)) +
  # new legend
  new_scale_fill() +
  geom_jjtriangle(aes(fill = `DNA damage`),type = 'ml') +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "#c11e21",
                       midpoint = 0, limit = c(-0.5,0.5)) +
  # new legend
  new_scale_fill() +
  geom_jjtriangle(aes(fill = `DNA repair`),type = 'mr') +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "#c11e21",
                       midpoint = 0, limit = c(-0.5,0.5)) +
  theme(legend.position = 'top') +
  coord_fixed()+
  theme_bw() + 
  theme(panel.grid=element_blank())


#UCEC TP DC Mast Fibroblasts
ucecmodel1 <- UCECsingle[UCECsingle$celltype%in%c("DC","Mast","Fibroblasts"),]
result3 <- ucecmodel1 %>%
  group_by(celltype) %>%
  summarise(Angiogenesis = median(Angiogenesis),
            EMT = median(EMT),
            Invasion = median(Invasion),
            Metastasis = median(Metastasis)
  )
ggplot(result3,
       aes(x = "",y = celltype)) +
  geom_jjtriangle(aes(fill = Angiogenesis),type = 'mu') +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "#c11e21",
                       midpoint = 0, limit = c(-0.5,0.5)) +
  # new legend
  new_scale_fill() +
  geom_jjtriangle(aes(fill = EMT),type = 'mb') +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "#c11e21",
                       midpoint = 0, limit = c(-0.5,0.5)) +
  # new legend
  new_scale_fill() +
  geom_jjtriangle(aes(fill = Invasion),type = 'ml') +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "#c11e21",
                       midpoint = 0, limit = c(-0.5,0.5)) +
  # new legend
  new_scale_fill() +
  geom_jjtriangle(aes(fill = Metastasis),type = 'mr') +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "#c11e21",
                       midpoint = 0, limit = c(-0.5,0.5)) +
  theme(legend.position = 'top') +
  coord_fixed()+
  theme_bw() + 
  theme(panel.grid=element_blank())
