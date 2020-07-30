### Plasmatocytes ###
library(Seurat)
library(ggplot2)
library(patchwork)
library(plyr)
library(dplyr)

filtergenes <- readRDS('tmp/filtergenes.Rds')
pm <- readRDS('tmp/pm.Rds')

head(pm@meta.data)
Idents(pm) <- factor(levels(Idents(pm)), levels = rev(levels(Idents(pm))))

DotPlot(pm, features = rev(c('srp', 'Tep4', 'Ance', 'Hml', 'Pxn', 
                             'Cdk1', 'stg', 'CycB', 'polo', 'aurB', 'Det', 'dUTPase', 'Pen', 'sle',
                             'Nplp2', 'vkg', 'NimC1', 'eater', 'Pvr', 'Ama', 'fax', 'vir-1', 'crq', 'msn'))) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

pm@meta.data$supergroup <- mapvalues(pm@meta.data$supergroup,
                                     from = unique(pm@meta.data$supergroup),
                                     to = c('PM 1', 'PM 2', 'PM 3', 'PM 4'))
DimPlot(pm, reduction = 'tsne', label = T, group.by = 'supergroup')

savelabel <- data.frame(pm@meta.data, check.rows = F, check.names = F)
savelabel <- savelabel[, c(2:8, 11)]
head(savelabel)
saveRDS(savelabel, 'tmp/pm.newlabel.Rds')

