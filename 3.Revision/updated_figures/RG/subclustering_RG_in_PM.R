### Plasmatocytes ###
library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

lymphgland.combined <- readRDS('../../../Drop-seq_alignment/1.seurat3_alignment_withMuscle_regress-Library-nUMI/lymphgland.combined.flt_allscaled.Rds')
head(lymphgland.combined@meta.data)

Idents(lymphgland.combined) <- 'anno_simple'
pm <- subset(lymphgland.combined, idents = c('PM', 'RG'))
remove(lymphgland.combined)

head(pm@meta.data)
unique(pm@meta.data$anno_simple)

VlnPlot(pm, features = c('Hml', 'phm'), group.by = 'anno_simple', cols = c('red2', 'purple'), pt.size = .01)
ggsave('rg_markers.pdf', units = 'cm', width = 6, height = 4.5)

