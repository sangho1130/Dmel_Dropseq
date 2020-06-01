library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

lymphgland <- readRDS('tmp/lymphgland.combined_transferred_filtered.Rds')
head(lymphgland@meta.data)
DimPlot(lymphgland, split.by = 'timepoint')

### marker list 20200429 ###
DefaultAssay(lymphgland) <- 'RNA'
lymphgland@meta.data$labelTransfer_s5d <- lymphgland@meta.data$labelTransfer
lymphgland@meta.data$labelTransfer_s5d <- factor(lymphgland@meta.data$labelTransfer_s5d, levels = rev(levels(lymphgland@meta.data$labelTransfer)))

DotPlot(lymphgland, group.by = 'labelTransfer_s5d', 
        split.by = 'timepoint', cols = c('dodgerblue', 'red2'), 
        features = rev(c('srp', 'Tep4', 'Ance', 'Hml', 'Pxn', 
                         'Antp', 'mthl7', 'tau', 'kn', 'chrb', 'Ilp6',
                         'Dl', 'E(spl)m4-BFM', 'E(spl)m3-HLH', 'dome', 'shg', 'Ppa', 'ci',
                         'NimB3', 'IM18', 'CecA2', 'yellow-f',
                         'Nplp2', 'vkg', 'NimC1', 'eater', 'Pvr', 'Ama', 'vir-1', 'crq',
                         'msn', 'Pvf2', 'mthl4', 'atilla', 'PPO3',
                         'lz', 'peb', 'PPO1', 'PPO2', 'fok', 'Men', 
                         'CG18547', 'CG3397', 'dys')), 
        dot.scale = 2) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = '', y = 'Subclusters')
#ggsave('FigureS5D_20200429.pdf', units = 'cm', width = 20, height = 15)



DotPlot(subset(lymphgland, cells = rownames(subset(lymphgland@meta.data, timepoint == 'Infested'))), 
        col.max = 2, col.min = -1.5,
        group.by = 'labelTransfer_s5d', 
        cols = c('grey90', 'red2'), 
        features = rev(c('srp', 'Tep4', 'Ance', 'Hml', 'Pxn', 
                         'Antp', 'mthl7', 'tau', 'kn', 'chrb', 'Ilp6',
                         'Dl', 'E(spl)m4-BFM', 'E(spl)m3-HLH', 'dome', 'shg', 'Ppa', 'ci',
                         'NimB3', 'IM18', 'CecA2', 'yellow-f',
                         'Nplp2', 'vkg', 'NimC1', 'eater', 'Pvr', 'Ama', 'vir-1', 'crq',
                         'msn', 'Pvf2', 'mthl4', 'atilla', 'PPO3',
                         'lz', 'peb', 'PPO1', 'PPO2', 'fok', 'Men', 
                         'CG18547', 'CG3397', 'dys')), 
        dot.scale = 3) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = '', y = 'Subclusters')
#ggsave('FigureS5D_Infested_20200429.pdf', units = 'cm', width = 20, height = 15)
DotPlot(subset(lymphgland, cells = rownames(subset(lymphgland@meta.data, timepoint == 'Normal'))), 
        col.max = 2, col.min = -1.5,
        group.by = 'labelTransfer_s5d', 
        cols = c('grey90', 'dodgerblue'), 
        features = rev(c('srp', 'Tep4', 'Ance', 'Hml', 'Pxn', 
                         'Antp', 'mthl7', 'tau', 'kn', 'chrb', 'Ilp6',
                         'Dl', 'E(spl)m4-BFM', 'E(spl)m3-HLH', 'dome', 'shg', 'Ppa', 'ci',
                         'NimB3', 'IM18', 'CecA2', 'yellow-f',
                         'Nplp2', 'vkg', 'NimC1', 'eater', 'Pvr', 'Ama', 'vir-1', 'crq',
                         'msn', 'Pvf2', 'mthl4', 'atilla', 'PPO3',
                         'lz', 'peb', 'PPO1', 'PPO2', 'fok', 'Men', 
                         'CG18547', 'CG3397', 'dys')), 
        dot.scale = 3) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = '', y = 'Subclusters')
#ggsave('FigureS5D_Normal_20200429.pdf', units = 'cm', width = 20, height = 15)

