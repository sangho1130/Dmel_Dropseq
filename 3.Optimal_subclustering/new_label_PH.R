### Prohemocyte ###
library(Seurat)
library(ggplot2)
library(patchwork)
library(plyr)
library(dplyr)

filtergenes <- readRDS('tmp/filtergenes.Rds')
ph <- readRDS('tmp/ph.Rds')

head(ph@meta.data)

ph@meta.data$supergroup <- mapvalues(ph@meta.data$supergroup,
                                     from = unique(ph@meta.data$supergroup),
                                     to = c('PH 1', 'PH 2', 'PH 6', 'PH 3', 'PH 4', 'PH 5'))
ph@meta.data$supergroup <- factor(ph@meta.data$supergroup, levels = c('PH 1', 'PH 2', 'PH 3', 'PH 4', 'PH 5', 'PH 6'))
DimPlot(ph, reduction = 'tsne', label = T, group.by = 'supergroup')

savelabel <- data.frame(ph@meta.data, check.rows = F, check.names = F)
savelabel <- savelabel[, c(2:8, 11)]
head(savelabel)
saveRDS(savelabel, 'tmp/ph.newlabel.Rds')

