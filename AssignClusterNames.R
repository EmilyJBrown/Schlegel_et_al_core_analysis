library(Seurat)
library(dplyr)

mar.merge <- readRDS("../martin.merge.CCA.rds")
summary(mar.merge@ident)

oldidents <- as.character(unique(mar.merge@ident))
newidents <- c("DCs CD4+", "B cells", "Developing T cells", "T cells 1", "T regs", 
               "CD8+ T cells/ NKT/ NK", "APOE Mø", "ILC3 and yd T cells", "T cells 2",
               "Granulocytes", "Monocytes", "Immature T cells", "LPL Mø", "DCs CD8+")
mar.merge@ident <- plyr::mapvalues(x=mar.merge@ident, from = oldidents, to=newidents)

