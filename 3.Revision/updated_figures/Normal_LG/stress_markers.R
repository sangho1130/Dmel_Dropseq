library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)


ph <- readRDS('../../PH/total ex_scbkgenes/tmp/ph.Rds')
ph_newlabel <- readRDS('../../PH/total ex_scbkgenes/tmp/ph.newlabel.Rds')
ph@meta.data <- ph_newlabel
Idents(ph) <- 'supergroup'

DimPlot(ph, reduction = 'tsne', label = T, label.size = 1) +
  scale_color_manual(values = c('#dee6d4', '#b9dcb1', '#78c0ac', '#359bb7', '#1e659b', '#1c306d')) + 
  theme_void() +
  theme(legend.position = 'None')

Idents(ph) <- factor(Idents(ph), levels = rev(levels(Idents(ph))))

DotPlot(ph, features = rev(c('CecA1', 'CecA2', 'CecB', 'CecC', 'PGRP-LB', 'AttA', 'Dro', ))) +
  labs(x = '', y = '') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

