library(Seurat)
library(dplyr)
library(SingleR)
setSeurat <- function(dirpath, name, projectName){
  name.data <- Read10X(data.dir = dirpath)
  name <- CreateSeuratObject(raw.data = name.data, project = projectName, min.cells = 3, min.genes = 200)
  name@meta.data$expt <- projectName
  name.mito <- grep(pattern = "mt-", x=rownames(x = name@data), value=TRUE)
  name.percent.mito <- Matrix::colSums(name@raw.data[name.mito, ])/Matrix::colSums(name@raw.data)
  name <- AddMetaData(object = name, metadata = name.percent.mito, col.name="percent.mito")
  return(name)
}

plotToFilter <- function(seurat){
  VlnPlot(object = seurat, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
}

filterMitoGenes <- function(seurat, maxmito, mingenes, maxgenes){
  seurat <- FilterCells(object = seurat, subset.names = c("percent.mito", "nGene"),
                        low.thresholds = c(-Inf, mingenes), 
                        high.thresholds = c(maxmito, maxgenes))
}

topVarGenes <- function(seurat){
  return(head(rownames(seurat@hvg.info),1000))
}

wt.bl <- setSeurat("../10x_Files/Martin.WT.BL/", wt.bl, "WT.BL")
plotToFilter(wt.bl)
wt.bl <- filterMitoGenes(wt.bl, 0.05, 500, 3000)
wt.bl <- NormalizeData(wt.bl, display.progress = F)
wt.bl <- ScaleData(wt.bl, vars.to.regress = c("nUMI", "percent.mito"))
wt.bl <- FindVariableGenes(wt.bl, do.plot = F)
saveRDS(wt.bl, file = "martin.wt.bl.norm.scale.regressnUMImito.rds")

ko.bl <- setSeurat("../10x_Files/Martin.KO.BL/", ko.bl, "KO.BL")
plotToFilter(ko.bl)
ko.bl <- filterMitoGenes(ko.bl, 0.05, 500, 2750)
ko.bl <- NormalizeData(ko.bl, display.progress = F)
ko.bl <- ScaleData(ko.bl, vars.to.regress = c("nUMI", "percent.mito"))
ko.bl <- FindVariableGenes(ko.bl, do.plot = F)
saveRDS(ko.bl, file = "martin.ko.bl.norm.scale.regressnUMImito.rds")

wt.reg <- setSeurat("../10x_Files/Martin.WT.Reg/", wt.reg, "WT.Reg")
plotToFilter(wt.reg)
wt.reg <- filterMitoGenes(wt.reg, 0.075, 0, 2000)
wt.reg <- NormalizeData(wt.reg, display.progress = F)
wt.reg <- ScaleData(wt.reg, vars.to.regress = c("nUMI", "percent.mito"))
wt.reg <- FindVariableGenes(wt.reg, do.plot = F)
saveRDS(wt.reg, file = "martin.wt.reg.norm.scale.regressnUMImito.rds")

ko.reg <- setSeurat("../10x_Files/Martin.KO.Reg/", ko.reg, "KO.Reg")
plotToFilter(ko.reg)
ko.reg <- filterMitoGenes(ko.reg, 0.05, 0, 3000)
ko.reg <- NormalizeData(ko.reg, display.progress = F)
ko.reg <- ScaleData(ko.reg, vars.to.regress = c("nUMI", "percent.mito"))
ko.reg <- FindVariableGenes(ko.reg, do.plot = F)
saveRDS(ko.reg, file = "martin.ko.reg.norm.scale.regressnUMImito.rds")

g.wt.bl <- topVarGenes(wt.bl)
g.ko.bl <- topVarGenes(ko.bl)
g.wt.reg <- topVarGenes(wt.reg)
g.ko.reg <- topVarGenes(ko.reg)
genes.use <- unique(c(g.wt.bl, g.ko.bl, g.wt.reg, g.ko.reg))
length(genes.use)
genes.use <- intersect(genes.use, rownames(wt.bl@scale.data))
genes.use <- intersect(genes.use, rownames(ko.bl@scale.data))
genes.use <- intersect(genes.use, rownames(wt.reg@scale.data))
genes.use <- intersect(genes.use, rownames(ko.reg@scale.data))

mar.list <- c(wt.bl, ko.bl, wt.reg, ko.reg)
cellids <- c("WT.BL", "KO.BL", "WT.Reg", "KO.Reg")

mar.merge.cca <- RunMultiCCA(object.list = mar.list, genes.use = genes.use, 
                             add.cell.ids = cellids, num.ccs = 25)
p1 <- DimPlot(object = mar.merge.cca, reduction.use = "cca", group.by = "expt",
              pt.size = 0.5, do.return = T)
p2 <- VlnPlot(mar.merge.cca, features.plot = "CC1", group.by = "expt", do.return = T)
p3 <- VlnPlot(mar.merge.cca, features.plot = "CC2", group.by = "expt", do.return = T)
p1
plot_grid(p2, p3)

mar.merge.cca <- AlignSubspace(mar.merge.cca, reduction.type = "cca", 
                               grouping.var = "expt", dims.align = 1:25)
VlnPlot(mar.merge.cca, features.plot = c("ACC1", "ACC2"), group.by = "expt", nCol=2)

mar.merge.cca <- RunTSNE(mar.merge.cca, reduction.use = "cca.aligned", dims.use = 1:25,
                         do.fast=T)
mar.merge.cca <- FindClusters(mar.merge.cca, reduction.type = "cca.aligned", 
                              resolution = 0.6, dims.use = 1:25, print.output = F)
saveRDS(mar.merge.cca, file = "martin.merge.CCA.rds")