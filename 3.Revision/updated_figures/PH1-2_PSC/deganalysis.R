library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

dir.create('tmp')
lymphgland <- readRDS('../Normal_LG/rdata/lymphgland.Rds')
levels(Idents(lymphgland))

earlycells <- subset(lymphgland, idents = c('PH 1', 'PH 2', 'PSC'))
remove(lymphgland)
earlycells@meta.data <- droplevels(earlycells@meta.data)
earlycells@meta.data$new_subclustering <- factor(earlycells@meta.data$new_subclustering, levels = c('PH 1', 'PH 2', 'PSC'))
Idents(earlycells) <- 'new_subclustering'
earlycells <- ScaleData(earlycells, vars.to.regress = c('Library', 'nCount_RNA'))
#saveRDS(earlycells, 'tmp/earlycells.Rds')
earlycells <- readRDS('tmp/earlycells.Rds')

earlycells.markers <- FindAllMarkers(object = earlycells, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
top10 <- earlycells.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

DoHeatmap(object = earlycells, features = c('N', 'Dl', 'shg', 'Socs36E', 'kn', top10$gene), angle = 90, size = 3, draw.lines = F, group.colors = c('#B9DCB1', '#78C0AC', '#EB76A7'), raster = F)
#ggsave('earlycells.pdf', units = 'cm', width = 15, height = 15)
#write.table(earlycells.markers, 'earlycells.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
#saveRDS(earlycells.markers, 'tmp/earlycells.markers.Rds')


VlnPlot(earlycells, features = c('kn'), group.by = 'new_subclustering', pt.size = 0, col = c('#B9DCB1', '#78C0AC', '#EB76A7')) + theme(legend.position = 'None')
ggsave('earlycells.kn.pdf', units = 'cm', width = 4, height = 6)

### Dl N in PH 1-2 ###
PH1PH2 <- subset(earlycells, idents = c('PH 1', 'PH 2'))
DefaultAssay(PH1PH2)
PH1PH2@meta.data <- droplevels(PH1PH2@meta.data)
PH1PH2 <- ScaleData(object = PH1PH2, vars.to.regress = c('Library', 'nCount_RNA'), features = rownames(PH1PH2))
#saveRDS(PH1PH2, 'tmp/PH1PH2.Rds')

Dl_N <- c()
smallset <- data.frame(GetAssayData(PH1PH2, slot = 'counts', assay = 'RNA'), check.rows = F, check.names = F)[c('Dl', 'N'),]
for (onecell in rownames(PH1PH2@meta.data)) {
  if (smallset['Dl', onecell] == 0) {
    if (smallset['N', onecell] == 0) {
      Dl_N <- append(Dl_N, 'Dl-N-')
    } else if (smallset['N', onecell] != 0) {
      Dl_N <- append(Dl_N, 'Dl-N+')
    }
  } else {
    if (smallset['N', onecell] == 0) {
      Dl_N <- append(Dl_N, 'Dl+N-')
    } else if (smallset['N', onecell] != 0) {
      Dl_N <- append(Dl_N, 'Dl+N+')
    }
  }
}
PH1PH2@meta.data$Dl_N <- Dl_N
Idents(PH1PH2) <- 'Dl_N'
Idents(PH1PH2) <- factor(Idents(PH1PH2), levels = c('Dl+N-', 'Dl+N+', 'Dl-N+', 'Dl-N-'))
#saveRDS(PH1PH2, 'tmp/PH1PH2.Rds')


counts <- data.frame(t(data.matrix(GetAssayData(PH1PH2, slot = 'data', assay = 'RNA'))[c('Dl', 'N'), ]), check.rows = F, check.names = F)
identical(rownames(counts), rownames(PH1PH2@meta.data))
counts$Subcluster <- PH1PH2@meta.data$new_subclustering
counts$Group <- PH1PH2@meta.data$Dl_N
counts$Group <- factor(counts$Group, levels = c('Dl+N-', 'Dl+N+', 'Dl-N+', 'Dl-N-'))
head(counts)
cor.test(counts$Dl, counts$N, method = 'spearman') # S = 310770, p-value = 3.275e-13, rho = 0.5361099

ggplot(counts, aes(Dl, N, col = Group, shape = Subcluster)) +
  geom_point(size = 3) + 
  scale_color_manual(values = c('#F79694', '#F4180D', '#931207', '#5B2A27')) +
  theme_bw() + 
  theme(axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        panel.grid = element_blank())
#ggsave('run_seurat3.Dl~N.pdf', units = 'cm', width = 10.5, height = 8)
ggplot(counts, aes(Dl, N, col = Group, shape = Subcluster)) +
  geom_point(size = 3) + 
  scale_color_manual(values = c('#F79694', '#F4180D', '#931207', '#5B2A27')) +
  theme_void() +
  theme(legend.position = 'None')
#ggsave('run_seurat3.Dl~N.points.pdf', units = 'cm', width = 6, height = 6)

PH1PH2 <- readRDS('tmp/PH1PH2.Rds')
PH1PH2.markers <- readRDS('tmp/PH1PH2.markers.Rds')
PH1PH2.markers <- FindAllMarkers(object = PH1PH2, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
top10 <- PH1PH2.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

DoHeatmap(object = PH1PH2, features = c('Dl', 'N', 'dome', 'shg', 'hop', 'Stat92E', 'Socs36E', top10$gene), angle = 90, size = 3, draw.lines = F, group.colors = c('#F79694', '#F4180D', '#931207', '#5B2A27'), raster = F)
#ggsave('PH1PH2_Dl-N.pdf', units = 'cm', width = 15, height = 15)
#write.table(PH1PH2.markers, 'PH1PH2_Dl-N.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
#saveRDS(PH1PH2.markers, 'tmp/PH1PH2.markers.Rds')


head(PH1PH2@meta.data)
VlnPlot(PH1PH2, features = c('Dl', 'N', 'shg', 'dome', 'Stat92E', 'Socs36E', 'Ance', 'Tep4', 'IM18'), split.by = 'new_subclustering', ncol = 3, pt.size = 0, cols = c('#B9DCB1', '#78C0AC'))
ggsave('run_seurat3.Dl~N.markersViolin.pdf', units = 'cm', width = 15, height = 15)


