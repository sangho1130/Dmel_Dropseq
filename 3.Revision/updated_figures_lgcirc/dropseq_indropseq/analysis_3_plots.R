library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

blood.combined_flt <- readRDS('tmp/blood.combined_newsubclustering_flt.Rds')
DefaultAssay(blood.combined_flt)
head(blood.combined_flt@meta.data)
levels(blood.combined_flt@meta.data$timepoint)

DimPlot(blood.combined_flt, dims = c(1,2), group.by = 'origin', cols = c('#ffa500', '#7ac5cd'))
#ggsave('umap/__flt__umap.1_2.byOrigin.seed320158.mindist_0.4_.pdf', units = 'cm', width = 14, height = 12)
plt <- DimPlot(blood.combined_flt, dims = c(1,2), group.by = 'origin', cells = rownames(subset(blood.combined_flt@meta.data, origin == 'LG')), cols = c('#ffa500'), pt.size = .5) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('umap/__flt__umap.1_2.byOrigin-LG.seed320158.mindist_0.4_.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined_flt, dims = c(1,2), group.by = 'origin', cells = rownames(subset(blood.combined_flt@meta.data, origin == 'Circ')), cols = c('#7ac5cd'), pt.size = .5) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('umap/__flt__umap.1_2.byOrigin-Circ.seed320158.mindist_0.4_.pdf', units = 'cm', width = 4, height = 4)


### timepoint v1 
plt <- DimPlot(blood.combined_flt, dims = c(1,2), group.by = 'origin', cells = rownames(subset(blood.combined_flt@meta.data, origin == 'LG' & timepoint == 'LG_drop_96')), cols = c('#ffda92'), pt.size = .5) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('umap/__flt__umap.1_2.byOrigin-LG_96.seed320158.mindist_0.4_.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined_flt, dims = c(1,2), group.by = 'origin', cells = rownames(subset(blood.combined_flt@meta.data, origin == 'LG' & timepoint == 'LG_drop_120')), cols = c('#ffa500'), pt.size = .5) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('umap/__flt__umap.1_2.byOrigin-LG_120.seed320158.mindist_0.4_.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined_flt, dims = c(1,2), group.by = 'origin', cells = rownames(subset(blood.combined_flt@meta.data, origin == 'Circ' & timepoint == 'Circ_drop_96')), cols = c('#a8e2d3'), pt.size = .5) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('umap/__flt__umap.1_2.byOrigin-Circ_96.seed320158.mindist_0.4_.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined_flt, dims = c(1,2), group.by = 'origin', cells = rownames(subset(blood.combined_flt@meta.data, origin == 'Circ' & timepoint == 'Circ_drop_120' | timepoint == 'Circ_indrop_120')), cols = c('#7ac5cd'), pt.size = .5) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('umap/__flt__umap.1_2.byOrigin-Circ_120.seed320158.mindist_0.4_.pdf', units = 'cm', width = 4, height = 4)


### timepoint v2 
DimPlot(blood.combined_flt, dims = c(1,2), group.by = 'timepoint', cols = c('#ffda92', '#ffa500', '#a8e2d3', '#7ac5cd', '#7ac5cd'))
#ggsave('umap/__flt__umap.1_2.bytimepoint.seed320158.mindist_0.4_.pdf', units = 'cm', width = 15, height = 12)
DimPlot(object = blood.combined_flt, dims = c(1,2), reduction = "umap", group.by = "timepoint", cols = c('#ffda92', '#ffa500', '#a8e2d3', '#7ac5cd', '#7ac5cd')) + facet_wrap(~timepoint, ncol = 2)
#ggsave('umap/__flt__umap.1_2.bytimepoint.seed320158.mindist_0.4_sep.pdf', units = 'cm', width = 16, height = 14)


### Annotation - simple
DimPlot(blood.combined_flt, dims = c(1,2), group.by = 'anno_simple', cols = c('#f15fa6', '#207eb3', '#a80d0c', '#f0a142', '#25a9b0', '#a4a4a4', '#1a1a1a'))
#ggsave('umap/__flt__umap.1_2.byanno_simple.seed320158.mindist_0.4_.pdf', units = 'cm', width = 16, height = 12)
plt <- DimPlot(blood.combined_flt, dims = c(1,2), group.by = 'anno_simple', cols = c('#f15fa6', '#207eb3', '#a80d0c', '#f0a142', '#25a9b0', '#a4a4a4', '#1a1a1a')) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('umap/__flt__umap.1_2.byanno_simple.seed320158.mindist_0.4_.augment.pdf', units = 'cm', width = 4, height = 4)


### Ubx ###
dir.create('markers')
plt <- FeaturePlot(blood.combined_flt, dims = c(1,2), cells = rownames(subset(blood.combined_flt@meta.data, origin == 'LG')), features = 'Ubx', pt.size = 1, cols = c('grey90', 'red2')) +
  theme_void() + theme(legend.position = 'None', plot.title = element_blank()); plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('markers/__flt__umap.1_2.Ubx_LG.pdf', units = 'cm', width = 4, height = 4)

plt <- FeaturePlot(blood.combined_flt, dims = c(1,2), cells = rownames(subset(blood.combined_flt@meta.data, origin == 'Circ')), features = 'Ubx', pt.size = 1, cols = c('grey90', 'red2')) +
  theme_void() + theme(legend.position = 'None', plot.title = element_blank()); plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('markers/__flt__umap.1_2.Ubx_Circ.pdf', units = 'cm', width = 4, height = 4)


### cell types
dir.create('umap/highlight')
# PH
plt <- DimPlot(blood.combined_flt, dims = c(1,2), cells = rownames(subset(blood.combined_flt@meta.data, origin == 'LG')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, anno_simple == 'PH')), cols.highlight = '#207eb3') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.PH-LG.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined_flt, dims = c(1,2), cells = rownames(subset(blood.combined_flt@meta.data, origin == 'Circ')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, anno_simple == 'PH')), cols.highlight = '#207eb3') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.PH-Circ.pdf', units = 'cm', width = 4, height = 4)


# PH 1
plt <- DimPlot(blood.combined_flt, dims = c(1,2), cells = rownames(subset(blood.combined_flt@meta.data, origin == 'LG')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, Subclustering == 'PH 1')), cols.highlight = '#207eb3') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.PH1-LG.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined_flt, dims = c(1,2), cells = rownames(subset(blood.combined_flt@meta.data, origin == 'Circ')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, Subclustering == 'PH 1')), cols.highlight = '#207eb3') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.PH1-Circ.pdf', units = 'cm', width = 4, height = 4)


# PH 4
plt <- DimPlot(blood.combined_flt, dims = c(1,2), cells = rownames(subset(blood.combined_flt@meta.data, origin == 'LG')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, Subclustering == 'PH 4')), cols.highlight = '#207eb3') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.PH4-LG.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined_flt, dims = c(1,2), cells = rownames(subset(blood.combined_flt@meta.data, origin == 'Circ')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, Subclustering == 'PH 4')), cols.highlight = '#207eb3') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.PH4-Circ.pdf', units = 'cm', width = 4, height = 4)


# PM
plt <- DimPlot(blood.combined_flt, dims = c(1,2), cells = rownames(subset(blood.combined_flt@meta.data, origin == 'LG')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, anno_simple == 'PM')), cols.highlight = '#a80d0c') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.PM-LG.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined_flt, dims = c(1,2), cells = rownames(subset(blood.combined_flt@meta.data, origin == 'Circ')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, anno_simple == 'PM')), cols.highlight = '#a80d0c') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.PM-Circ.pdf', units = 'cm', width = 4, height = 4)


# CC
plt <- DimPlot(blood.combined_flt, dims = c(1,2), cells = rownames(subset(blood.combined_flt@meta.data, origin == 'LG')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, anno_simple == 'CC')), cols.highlight = '#25a9b0') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.CC-LG.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined_flt, dims = c(1,2), cells = rownames(subset(blood.combined_flt@meta.data, origin == 'Circ')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, anno_simple == 'CC')), cols.highlight = '#25a9b0') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.CC-Circ.pdf', units = 'cm', width = 4, height = 4)


# CC 1
plt <- DimPlot(blood.combined_flt, dims = c(1,2), cells = rownames(subset(blood.combined_flt@meta.data, origin == 'LG')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, Subclustering == 'CC 1')), cols.highlight = '#72c4b2') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.CC1-LG.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined_flt, dims = c(1,2), cells = rownames(subset(blood.combined_flt@meta.data, origin == 'Circ')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, Subclustering == 'CC 1')), cols.highlight = '#72c4b2') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.CC1-Circ.pdf', units = 'cm', width = 4, height = 4)


# CC 2
plt <- DimPlot(blood.combined_flt, dims = c(1,2), cells = rownames(subset(blood.combined_flt@meta.data, origin == 'LG')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, Subclustering == 'CC 2')), cols.highlight = '#1b7e7d') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.CC2-LG.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined_flt, dims = c(1,2), cells = rownames(subset(blood.combined_flt@meta.data, origin == 'Circ')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, Subclustering == 'CC 2')), cols.highlight = '#1b7e7d') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.CC2-Circ.pdf', units = 'cm', width = 4, height = 4)

