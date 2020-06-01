library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

dir.create('tmp')

### LG ###
lymphgland <- readRDS('../rdata/lymphgland.Rds')
levels(Idents(lymphgland))
lymphgland.markers <- FindAllMarkers(object = lymphgland, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
top10 <- lymphgland.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

DoHeatmap(object = lymphgland, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('lymphgland_subclustering.pdf', units = 'cm', width = 40, height = 30)
write.table(lymphgland.markers, 'lymphgland_subclustering.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
saveRDS(lymphgland.markers, 'tmp/lymphgland_subclustering.markers.Rds')

head(lymphgland@meta.data)
lymphgland@meta.data$new_anno_simple <- mapvalues(lymphgland@meta.data$new_subclustering,
                                                  from = levels(lymphgland@meta.data$new_subclustering),
                                                  to = c('PSC', 'PH', 'PH', 'PH', 'PH', 'PH', 'PH', 
                                                         'PM', 'PM', 'PM', 'PM', 'LM', 'LM', 'CC', 'CC',
                                                         'GST-rich', 'Adipohemocyte', 'DV', 'RG', 'Neurons'))
Idents(lymphgland) <- 'new_anno_simple'
levels(Idents(lymphgland))
lymphgland.markers <- FindAllMarkers(object = lymphgland, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
top10 <- lymphgland.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

DoHeatmap(object = lymphgland, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('lymphgland_celltype.pdf', units = 'cm', width = 40, height = 30)
write.table(lymphgland.markers, 'lymphgland_celltype.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
saveRDS(lymphgland.markers, 'tmp/lymphgland_celltype.markers.Rds')


### PH ###
ph <- readRDS('../../../PH/total ex_scbkgenes/tmp/ph.Rds')
ph_newlabel <- readRDS('../../../PH/total ex_scbkgenes/tmp/ph.newlabel.Rds')
ph@meta.data <- ph_newlabel
Idents(ph) <- 'supergroup'

ph.markers <- FindAllMarkers(object = ph, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
top10 <- ph.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

DoHeatmap(object = ph, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('ph.pdf', units = 'cm', width = 30, height = 20)
write.table(ph.markers, 'ph.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
saveRDS(ph.markers, 'tmp/ph.markers.Rds')

### PM ###
pm <- readRDS('../../..//PM/total ex_scbkgenes/tmp/pm.Rds')
pm_newlabel <- readRDS('../../../PM/total ex_scbkgenes/tmp/pm.newlabel.Rds')
pm@meta.data <- pm_newlabel
Idents(pm) <- 'supergroup'

pm.markers <- FindAllMarkers(object = pm, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
top10 <- pm.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

DoHeatmap(object = pm, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('pm.pdf', units = 'cm', width = 30, height = 20)
write.table(pm.markers, 'pm.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
saveRDS(pm.markers, 'tmp/pm.markers.Rds')
