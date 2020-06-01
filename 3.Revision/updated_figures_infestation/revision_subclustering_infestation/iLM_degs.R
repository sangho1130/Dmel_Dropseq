library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)


#dir.create('degs_lm')
lymphgland.combined <- readRDS('tmp/lymphgland.combined_transferred_filtered.Rds')

### LM 1 vs LM 2 ###
DefaultAssay(lymphgland.combined) <- 'RNA'
head(lymphgland.combined@meta.data)
lymphgland.combined.LM <- subset(lymphgland.combined, idents = c('LM 1', 'LM 2'))
lymphgland.combined.LM <- subset(lymphgland.combined.LM, cells = rownames(subset(lymphgland.combined.LM@meta.data, timepoint == 'Infested')))
lymphgland.combined.LM <- ScaleData(lymphgland.combined.LM, vars.to.regress = c('Library', 'nCount_RNA'), features = rownames(lymphgland.combined.LM))


LM.markers <- FindAllMarkers(lymphgland.combined.LM, min.pct = .25, logfc.threshold = .25, only.pos = T)
#write.table(LM.markers, 'degs_lm/findAllMarkers.LM.labelTransfer.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- LM.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

dehm <- DoHeatmap(object = lymphgland.combined.LM, features = top10$gene, angle = 90, size = 3, raster = F, draw.lines = F); dehm
#ggsave('degs_lm/findAllMarkers.LM.labelTransfer.pdf', units = 'cm', width = 12, height = 10)
