library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

lymphgland <- readRDS('rdata/lymphgland.Rds')

### Subclustering ###
### PH
DimPlot(lymphgland, reduction = 'tsne', pt.size = .5, label = T, label.size = 1) + 
  scale_color_manual(values = c('grey90', #PSC
                                '#dee6d4', '#b9dcb1', '#78c0ac', '#359bb7', '#1e659b', '#1c306d', #PH
                                'grey90', 'grey90', 'grey90', 'grey90', #PM
                                'grey90', 'grey90', 'grey90', 'grey90', #LM, CC
                                'grey90', 'grey90', #GST, Adipo
                                'grey90', 'grey90', 'grey90')) + #DV, RG, Neu
  theme_void() +
  theme(legend.position = 'None')
#ggsave('subclustering_ph.pdf', units = 'cm', height = 4, width = 4)
plt <- DimPlot(lymphgland, reduction = 'tsne', pt.size = .5) + 
  scale_color_manual(values = c('grey90', #PSC
                                '#dee6d4', '#b9dcb1', '#78c0ac', '#359bb7', '#1e659b', '#1c306d', #PH
                                'grey90', 'grey90', 'grey90', 'grey90', #PM
                                'grey90', 'grey90', 'grey90', 'grey90', #LM, CC
                                'grey90', 'grey90', #GST, Adipo
                                'grey90', 'grey90', 'grey90')) + #DV, RG, Neu
  theme_void() +
  theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('subclustering_ph_augment.pdf', units = 'cm', height = 4, width = 4)

tsne <- data.frame(Embeddings(lymphgland, reduction = 'tsne'), check.rows = F, check.names = F)
tsne$celltype <- lymphgland@meta.data$anno_simple
tsne$subclustering <- lymphgland@meta.data$new_subclustering
head(tsne)
ggplot() +
  geom_point(data = subset(tsne, celltype != 'PH'), aes(tSNE_1, tSNE_2), col = 'grey90', size =.1) +
  geom_point(data = subset(tsne, celltype == 'PH'), aes(tSNE_1, tSNE_2, col = subclustering), size =.1) +
  scale_color_manual(values = c('#dee6d4', '#b9dcb1', '#78c0ac', '#359bb7', '#1e659b', '#1c306d')) +
  theme_void() +
  theme(legend.position = 'None')
#ggsave('subclustering_ph.png', units = 'cm', height = 4, width = 4)



### PM 
DimPlot(lymphgland, reduction = 'tsne', pt.size = .5, label = T, label.size = 1) + 
  scale_color_manual(values = c('grey90', #PSC
                                'grey90', 'grey90', 'grey90', 'grey90', 'grey90', 'grey90', 'grey90', #PH
                                '#ecd5a5', '#f6ae72', '#ea6740', '#cf231c', #PM
                                'grey90', 'grey90', 'grey90', 'grey90', #LM, CC
                                'grey90', 'grey90', #GST, Adipo
                                'grey90', 'grey90', 'grey90')) + #DV, RG, Neu
  theme_void() +
  theme(legend.position = 'None')
#ggsave('subclustering_pm.pdf', units = 'cm', height = 4, width = 4)
plt <- DimPlot(lymphgland, reduction = 'tsne', pt.size = .5) + 
  scale_color_manual(values = c('grey90', #PSC
                                'grey90', 'grey90', 'grey90', 'grey90', 'grey90', 'grey90', 'grey90', #PH
                                '#ecd5a5', '#f6ae72', '#ea6740', '#cf231c', #PM
                                'grey90', 'grey90', 'grey90', 'grey90', #LM, CC
                                'grey90', 'grey90', #GST, Adipo
                                'grey90', 'grey90', 'grey90')) + #DV, RG, Neu
  theme_void() +
  theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('subclustering_pm_augment.pdf', units = 'cm', height = 4, width = 4)

tsne <- data.frame(Embeddings(lymphgland, reduction = 'tsne'), check.rows = F, check.names = F)
tsne$celltype <- lymphgland@meta.data$anno_simple
tsne$subclustering <- lymphgland@meta.data$new_subclustering
head(tsne)
ggplot() +
  geom_point(data = subset(tsne, celltype != 'PM'), aes(tSNE_1, tSNE_2), col = 'grey90', size =.1) +
  geom_point(data = subset(tsne, celltype == 'PM'), aes(tSNE_1, tSNE_2, col = subclustering), size =.1) +
  scale_color_manual(values = c('#ecd5a5', '#f6ae72', '#ea6740', '#cf231c')) +
  theme_void() +
  theme(legend.position = 'None')
#ggsave('subclustering_pm.png', units = 'cm', height = 4, width = 4)


### Dot plot ###
lymphgland@meta.data$new_subclustering <- factor(lymphgland@meta.data$new_subclustering, levels = rev(levels(Idents(lymphgland))))
DotPlot(lymphgland, group.by = 'new_subclustering', features = rev(c('stg', 'polo', 'Cdk1', 'aurB', 'Det', 'CycB', 'dUTPase', 'Pen', 'sle')), 
        cols = c('grey90', 'red2'), dot.scale = 4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = '', y = 'Subclusters')
#ggsave('Figure2B_cellcycle.pdf', units = 'cm', width = 13, height = 17)
DotPlot(lymphgland, group.by = 'new_subclustering', 
        features = rev(c('CG15550', 'CG6023', 'mthl7', 'ImpL2', 'Pepck', #PSC
                         
                         'pncr002:3R', 'CG43164', 'E(spl)m4-BFM', 'E(spl)m2-BFM', 'E(spl)malpha-BFM', # PH1
                         'NimB3', 'loh', 'Mal-A5', 'Hand', 'egr', # PH2
                         'CG15115', 'rdgA', 'CG30054', 'dlp', 'IM18', # PH3
                         'CG13160', 'NK7.1', 'Swim', 'CG11236', 'Sr-CIII', # PH4
                         'Lip4', 'CG10433', 'CG31174', 'CG33225', 'CR46093', # PH5
                         'CG34296', 'Dtg', 'Obp99b', 'sn', 'PGRP-SC2', # PH6
                         
                         'CG34054', 'lectin-24A', 'Tig', 'CG34437', 'Ance-5', # PM1
                         'Drs', 'NimC1', 'CG14291', 'CG31673', 'CG8501', # PM2
                         'Cpr49Ac', 'fax', 'Ama', 'betaTub60D', 'Hsp23', # PM3
                         'CG7778', 'Gs1', 'NLaz', 'CG13321', 'CG4950', # PM4

                         'CG42792', 'CG18557', 'CG33458', 'mthl4', 'CR44316', # LM1
                         'atilla', 'CG1208', 'CG14610', 'CR44948', 'CG2556', # LM2
                         
                         'PPO2', 'PPO1', 'peb', 'lz', 'klu', # CC1
                         'fok', 'CG15343', 'CG9119', 'CG17109', 'CG5828', # CC2
                         
                         'CG15784', 'CR44430', 'CG18547', 'dys', 'CG3397', # GST
                         'E23', 'CG11899', 'hid', 'Gal', 'GstT4', # Adipo
                         
                         'Mlc2', 'Mf', 'Mlp60A', 'Mlc1', 'up',# DV
                         'CG15506', 'Timp', 'CG4408', 'CG42680', 'if', # RG
                         'CR43283', 'Appl', 'bru-3', 'Argk', 'CG13928'# Neuron
                         )), 
        cols = c('grey90', 'red2'), dot.scale = 4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = '', y = 'Subclusters')
#ggsave('Figure2B_markers.pdf', units = 'cm', width = 45, height = 17)

### Figure 2B marker list 20200429 ###
DotPlot(lymphgland, group.by = 'new_subclustering', 
        features = rev(c('srp', 'Tep4', 'Ance', 'Hml', 'Pxn', 
                         'Antp', 'mthl7', 'tau', 'kn', 'chrb', 'Ilp6',
                         'Dl', 'E(spl)m4-BFM', 'E(spl)m3-HLH', 'dome', 'shg', 'Ppa', 'ci',
                         'NimB3', 'IM18', 'CecA2', 'yellow-f',
                         'Nplp2', 'vkg', 'NimC1', 'eater', 'Pvr', 'Ama', 'vir-1', 'crq',
                         'msn', 'Pvf2', 'mthl4', 'atilla', 'PPO3',
                         'lz', 'peb', 'PPO1', 'PPO2', 'fok', 'Men', 
                         'CG18547', 'CG3397', 'dys', 'Sirup', 'Jabba', 'Lsd-2', 
                         'Mlc2', 'Hand', 'phm', 'nSyb')), 
        cols = c('grey90', 'red2'), dot.scale = 4) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = '', y = 'Subclusters')
#ggsave('Figure2B_20200429.pdf', units = 'cm', width = 35, height = 17)



### Main Figure 2A
ph <- readRDS('../../PH/total ex_scbkgenes/tmp/ph.Rds')
ph_newlabel <- readRDS('../../PH/total ex_scbkgenes/tmp/ph.newlabel.Rds')
ph@meta.data <- ph_newlabel
Idents(ph) <- 'supergroup'

DimPlot(ph, reduction = 'tsne', label = T, label.size = 1) +
  scale_color_manual(values = c('#dee6d4', '#b9dcb1', '#78c0ac', '#359bb7', '#1e659b', '#1c306d')) + 
  theme_void() +
  theme(legend.position = 'None')
#ggsave('Figure2A_ph.pdf', units = 'cm', height = 4, width = 4)
DimPlot(ph, reduction = 'tsne', label = T, label.size = 1) +
  scale_color_manual(values = c('#dee6d4', '#b9dcb1', '#78c0ac', '#359bb7', '#1e659b', '#1c306d')) + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(colour = 'black'), axis.ticks = element_line(colour = 'black'), legend.position = 'None')
#ggsave('Figure2A_ph_2.pdf', units = 'cm', height = 8, width = 8.5)
plt <- DimPlot(ph, reduction = 'tsne', pt.size = 1) +
  scale_color_manual(values = c('#dee6d4', '#b9dcb1', '#78c0ac', '#359bb7', '#1e659b', '#1c306d')) + 
  theme_void() +
  theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('Figure2A_ph_augment.pdf', units = 'cm', height = 4, width = 4)
DimPlot(ph, reduction = 'tsne', pt.size = .1) +
  scale_color_manual(values = c('#dee6d4', '#b9dcb1', '#78c0ac', '#359bb7', '#1e659b', '#1c306d')) + 
  theme_void() +
  theme(legend.position = 'None')
#ggsave('Figure2A_ph.png', units = 'cm', height = 4, width = 4)


pm <- readRDS('../../PM/total ex_scbkgenes/tmp/pm.Rds')
pm_newlabel <- readRDS('../../PM/total ex_scbkgenes/tmp/pm.newlabel.Rds')
pm@meta.data <- pm_newlabel
Idents(pm) <- 'supergroup'

DimPlot(pm, reduction = 'tsne', label = T, label.size = 1) +
  scale_color_manual(values = c('#FCB17B', '#F16C4B', '#C81C12', '#7F0000')) + 
  theme_void() +
  theme(legend.position = 'None')
#ggsave('Figure2A_pm.pdf', units = 'cm', height = 4, width = 4)
DimPlot(pm, reduction = 'tsne', label = T, label.size = 1) +
  scale_color_manual(values = c('#FCB17B', '#F16C4B', '#C81C12', '#7F0000')) + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(colour = 'black'), axis.ticks = element_line(colour = 'black'), legend.position = 'None')
#ggsave('Figure2A_pm_2.pdf', units = 'cm', height = 8, width = 8.5)
plt <- DimPlot(pm, reduction = 'tsne', pt.size = 1) +
  scale_color_manual(values = c('#FCB17B', '#F16C4B', '#C81C12', '#7F0000')) + 
  theme_void() +
  theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('Figure2A_pm_augment.pdf', units = 'cm', height = 4, width = 4)
DimPlot(pm, reduction = 'tsne', pt.size = .1) +
  scale_color_manual(values = c('#FCB17B', '#F16C4B', '#C81C12', '#7F0000')) + 
  theme_void() +
  theme(legend.position = 'None')
#ggsave('Figure2A_pm.png', units = 'cm', height = 4, width = 4)
FeaturePlot(pm, reduction = 'tsne', features = c('Hml', 'NimC1', 'AttA', 'CecA1', 'fax', 'stg'))
FeaturePlot(pm, reduction = 'tsne', features = c('Ance', 'Nplp2', 'Hml', 'NimC1', 'stg'), cols = c('grey90', 'red2'), pt.size = 1, max.cutoff = 2.5)
FeaturePlot(pm, reduction = 'tsne', features = 'NimC1', pt.size = 1, max.cutoff = 2.5)

head(lymphgland@meta.data)
lymphgland@meta.data$new_anno_simple <- mapvalues(lymphgland@meta.data$new_subclustering,
                                                  from = levels(lymphgland@meta.data$new_subclustering),
                                                  to = c('PSC', 'PH', 'PH', 'PH', 'PH', 'PH', 'PH', 
                                                         'PM', 'PM', 'PM', 'PM', 'LM', 'LM', 'CC', 'CC',
                                                         'GST-rich', 'Adipohemocyte', 'DV', 'RG', 'Neurons'))

lymphgland@meta.data <- lymphgland@meta.data[, c(2:6, 9, 10)]
writeTable <- data.frame(Barcode = rownames(lymphgland@meta.data), lymphgland@meta.data, check.rows = F, check.names = F)
write.table(writeTable, 'new_clustering_label.txt', quote = F, sep = '\t', row.names = F, col.names = T)



### Cell cycle ###
head(lymphgland@meta.data)
VlnPlot(lymphgland, features = c('stg', 'polo', 'Cdk1', 'aurB', 'Det', 'CycB', 'dUTPase', 'Pen', 'sle'), ncol = 3)
#ggsave('subclustering_cellcycling.pdf', units = 'cm', height = 30, width = 40)

Idents(lymphgland) <- 'new_subclustering'
lymphgland@meta.data$new_anno_simple <- mapvalues(lymphgland@meta.data$new_subclustering,
                                                  from = levels(lymphgland@meta.data$new_subclustering),
                                                  to = c('PSC', 'PH', 'PH', 'PH', 'PH', 'PH', 'PH', 'PM', 'PM', 'PM', 'PM', 
                                                         'LM', 'LM', 'CC', 'CC', 'GST-rich', 'Adipohemocyte', 'DV', 'RG', 'Neurons'))
Idents(lymphgland) <- 'new_anno_simple'
lymphgland_phpm <- subset(lymphgland, idents = c('PH', 'PM'))
Idents(lymphgland_phpm) <- 'new_subclustering'
lymphgland_phpm@meta.data <- droplevels(lymphgland_phpm@meta.data)
levels(lymphgland_phpm@meta.data$new_subclustering)

lymphgland_phpm_v2 <- subset(lymphgland_phpm, cells = rownames(subset(lymphgland_phpm@meta.data, new_subclustering != 'PH 1' & new_subclustering != 'PH 2')))
DotPlot(lymphgland_phpm_v2, features = rev(c('Cdk1', 'stg', 'CycB', 'polo', 'aurB', 'Det', 'dUTPase', 'Pen', 'sle'))) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

lymphgland_phpm_v3 <- subset(lymphgland_phpm, cells = rownames(subset(lymphgland_phpm@meta.data, new_subclustering != 'PH 1' & new_subclustering != 'PH 2' & new_subclustering != 'PH 3')))
DotPlot(lymphgland_phpm_v3, features = rev(c('Cdk1', 'stg', 'CycB', 'polo', 'aurB', 'Det', 'dUTPase', 'Pen', 'sle'))) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))


empty <- data.frame(matrix(nrow = 9, ncol = 10))
rownames(empty) <- c('Cdk1', 'stg', 'CycB', 'polo', 'aurB', 'Det', 'dUTPase', 'Pen', 'sle')
colnames(empty) <- levels(lymphgland_phpm@meta.data$new_subclustering)
empty
Idents(lymphgland_phpm) <- 'new_subclustering'
for (subcluster in levels(lymphgland_phpm@meta.data$new_subclustering)) {
  tmpobj <- subset(lymphgland_phpm, idents = subcluster)
  for (ccgene in c('Cdk1', 'stg', 'CycB', 'polo', 'aurB', 'Det', 'dUTPase', 'Pen', 'sle')) {
    tmpexpr <- GetAssayData(tmpobj, slot = 'data')[ccgene,]
    tmpexpr_mean <- mean(tmpexpr)
    tmpexpr_median <- median(tmpexpr)
    empty[ccgene, subcluster] <- tmpexpr_mean
  }
}

empty



FeaturePlot(lymphgland, features = c('Ance', 'Nplp2', 'Hml', 'NimC1'), cols = c('grey90', 'red2'))
