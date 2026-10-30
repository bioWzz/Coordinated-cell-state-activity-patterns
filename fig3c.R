#devtools::install_github('Asa12138/pcutils',dependencies=T)
library(pcutils)

sunbursttest <- read.csv("D:/细胞状态2/论文/fig3/fig3c/sunburst.csv")
#my_sunburst(sunbursttest)

gghuan2(sunbursttest)

