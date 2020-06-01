library(Seurat)
library(ggplot2)
library(ggrepel)
library(plyr)
library(dplyr)

# biased genes
scGenes <- read.delim('../dropseq/__filter_scGenes__/bulk_sc_pseudo_pt3-scGenes.scGenes_v2.txt', header = F)
scGenes <- scGenes$V1
head(scGenes); length(scGenes)

bkGenes <- read.delim('../dropseq/__filter_scGenes__/bulk_sc_pseudo_pt3-scGenes.bkGenes_v2.txt', header = F)
bkGenes <- bkGenes$V1
head(bkGenes); length(bkGenes)

# Dataset 
blood.combined_flt <- readRDS('tmp/blood.combined_newsubclustering_flt.Rds')
#blood.combined_flt <- ScaleData(object = blood.combined_flt, vars.to.regress = c('Library', 'nCount_RNA', 'percent.mt'), features = rownames(blood.combined_flt))
#saveRDS(blood.combined_flt, 'tmp/blood.combined_newsubclustering_flt.Rds')

DimPlot(blood.combined_flt, dims = c(1, 2))

dir.create('compare')
##############
### Origin ###
##############
dir.create('compare/origin')
head(blood.combined_flt@meta.data)
Idents(blood.combined_flt) <- "origin"

blood.combined_flt_markers <- FindAllMarkers(blood.combined_flt, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
blood.combined_flt_markers <- subset(blood.combined_flt_markers, p_val_adj <= 0.05 & avg_logFC >= 1)
blood.combined_flt_markers <- blood.combined_flt_markers[setdiff(rownames(blood.combined_flt_markers), scGenes), ]
nrow(blood.combined_flt_markers)
blood.combined_flt_markers <- blood.combined_flt_markers[setdiff(rownames(blood.combined_flt_markers), bkGenes), ]##
nrow(blood.combined_flt_markers)

top10 <- blood.combined_flt_markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = blood.combined_flt, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('compare/origin/blood.combined_flt_MAST_degs.pdf', units = 'cm', width = 20, height = 20)
write.table(blood.combined_flt_markers, 'compare/origin/blood.combined_flt_MAST_degs.origin.txt', quote = F, sep = '\t', row.names = F, col.names = T)

avg_blood.combined_flt <- log1p(AverageExpression(blood.combined_flt, verbose = FALSE)$RNA)
avg_blood.combined_flt$gene <- rownames(avg_blood.combined_flt)

avg_blood.combined_flt$labelgenes <- NA
avg_blood.combined_flt[c(top10$gene), 'labelgenes'] <- c(top10$gene)

avg_blood.combined_flt$colorgenes <- 'Others'
avg_blood.combined_flt[rownames(subset(blood.combined_flt_markers, avg_logFC >= 1 & cluster == 'LG')), 'colorgenes'] <- 'L'
avg_blood.combined_flt[rownames(subset(blood.combined_flt_markers, avg_logFC >= 1 & cluster == 'Circ')), 'colorgenes'] <- 'C'

plt <- ggplot(avg_blood.combined_flt, aes(LG, Circ, label = labelgenes, col = colorgenes)) + 
  geom_point() + scale_color_manual(values = c('#69b9c2', '#fd9409', 'grey90')) + 
  ggtitle("blood.combined_flt") + xlim(0, 8) + ylim(0, 8) + 
  geom_text_repel(data = avg_blood.combined_flt[c(top10$gene[1:10]), ], color = '#fd9409',
                  nudge_x = 8 - avg_blood.combined_flt[c(top10$gene[1:10]), ]$LG,
                  segment.size = 0.2, segment.color = "grey50", direction = "y", hjust = 1) +
  geom_text_repel(data = avg_blood.combined_flt[top10$gene[11:20], ], 
                  nudge_y = 8 - avg_blood.combined_flt[top10$gene[11:20], ]$Circ,
                  segment.size = 0.2, segment.color = "grey50", direction = "x", vjust = 0, angle = 90) +
  theme_bw() +
  theme(legend.position = 'None', panel.grid = element_blank(), plot.title = element_blank());plt
ggsave('compare/origin/blood.combined_flt_scatter.origin.pdf', units = 'cm', width = 10, height = 10)

plt <- ggplot() +
  geom_point(data = subset(avg_blood.combined_flt, colorgenes == 'Others'), aes(LG, Circ, col = colorgenes), size = .5) + 
  geom_point(data = subset(avg_blood.combined_flt, colorgenes == 'C'), aes(LG, Circ, col = colorgenes), size = .5) + 
  geom_point(data = subset(avg_blood.combined_flt, colorgenes == 'L'), aes(LG, Circ, col = colorgenes), size = .5) + 
  scale_color_manual(values = c('#69b9c2', '#fd9409', 'grey90')) + 
  theme_void() + theme(legend.position = 'None') + xlim(0, 8) + ylim(0, 8);plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('compare/origin/blood.combined_flt_scatter.origin.augment.pdf', units = 'cm', width = 4, height = 4)
ggsave('compare/origin/blood.combined_flt_scatter.origin.augment.png', units = 'cm', width = 4, height = 4)


##########
### PH ###
##########
dir.create('compare/ph')
head(blood.combined_flt@meta.data)
Idents(blood.combined_flt) <- "anno_simple"
ph <- subset(blood.combined_flt, idents = "PH")
Idents(ph) <- "origin"

ph_markers <- FindAllMarkers(ph, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
ph_markers <- subset(ph_markers, p_val_adj <= 0.05 & avg_logFC >= 1)
ph_markers <- ph_markers[setdiff(rownames(ph_markers), scGenes), ]
nrow(ph_markers)
ph_markers <- ph_markers[setdiff(rownames(ph_markers), bkGenes), ]
nrow(ph_markers)


top10 <- ph_markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = ph, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('compare/ph/blood.combined_flt_MAST_degs.PH.pdf', units = 'cm', width = 30, height = 20)
write.table(ph_markers, 'compare/ph/blood.combined_flt_MAST_degs.PH.txt', quote = F, sep = '\t', row.names = F, col.names = T)

avg_ph <- log1p(AverageExpression(ph, verbose = FALSE)$RNA)
avg_ph$gene <- rownames(avg_ph)

avg_ph$labelgenes <- NA
avg_ph[c(top10$gene), 'labelgenes'] <- c(top10$gene)

avg_ph$colorgenes <- 'Others'
avg_ph[rownames(subset(ph_markers, avg_logFC >= 1 & cluster == 'LG')), 'colorgenes'] <- 'L'
avg_ph[rownames(subset(ph_markers, avg_logFC >= 1 & cluster == 'Circ')), 'colorgenes'] <- 'C'

plt <- ggplot(avg_ph, aes(LG, Circ, label = labelgenes, col = colorgenes)) + 
  geom_point() + scale_color_manual(values = c('#69b9c2', '#fd9409', 'grey90')) + 
  ggtitle("blood.combined_flt") + xlim(0, 8) + ylim(0, 8) + 
  geom_text_repel(data = avg_ph[c(top10$gene[1:10]), ], color = '#fd9409',
                  nudge_x = 8 - avg_ph[c(top10$gene[1:10]), ]$LG,
                  segment.size = 0.2, segment.color = "grey50", direction = "y", hjust = 1) +
  geom_text_repel(data = avg_ph[top10$gene[11:20], ], 
                  nudge_y = 8 - avg_ph[top10$gene[11:20], ]$Circ,
                  segment.size = 0.2, segment.color = "grey50", direction = "x", vjust = 0, angle = 90) +
  theme_bw() +
  theme(legend.position = 'None', panel.grid = element_blank(), plot.title = element_blank());plt
ggsave('compare/ph/blood.combined_flt_scatter.PH.pdf', units = 'cm', width = 10, height = 10)


plt <- ggplot() +
  geom_point(data = subset(avg_ph, colorgenes == 'Others'), aes(LG, Circ, col = colorgenes), size = .5) + 
  geom_point(data = subset(avg_ph, colorgenes == 'C'), aes(LG, Circ, col = colorgenes), size = .5) + 
  geom_point(data = subset(avg_ph, colorgenes == 'L'), aes(LG, Circ, col = colorgenes), size = .5) + 
  scale_color_manual(values = c('#69b9c2', '#fd9409', 'grey90')) + 
  theme_void() + theme(legend.position = 'None') + xlim(0, 8) + ylim(0, 8);plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('compare/ph/blood.combined_flt_scatter.PH.augment.pdf', units = 'cm', width = 4, height = 4)
ggsave('compare/ph/blood.combined_flt_scatter.PH.augment.png', units = 'cm', width = 4, height = 4)



############
### PH 1 ###
############
dir.create('compare/ph1')
Idents(blood.combined_flt) <- 'Subclustering'
levels(Idents(blood.combined_flt))

conserved_markers_ph1 <- FindConservedMarkers(blood.combined_flt, ident.1 = 'PH 1', grouping.var = 'origin', min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
head(conserved_markers_ph1)
conserved_markers_ph1 <- subset(conserved_markers_ph1, conserved_markers_ph1[,5] <= 0.05 &  conserved_markers_ph1[,10] <= 0.05 &
                                  conserved_markers_ph1[,2] >= 1 & conserved_markers_ph1[,7] >= 1)
conserved_markers_ph1 <- conserved_markers_ph1[setdiff(rownames(conserved_markers_ph1), scGenes), ]
nrow(conserved_markers_ph1)
conserved_markers_ph1_write <- data.frame(Symbol = rownames(conserved_markers_ph1), conserved_markers_ph1)
write.table(conserved_markers_ph1_write, 'compare/ph1/conserved_markers_ph1.txt', quote = F, sep = '\t', row.names = F, col.names = T)
DoHeatmap(object = blood.combined_flt, features = rownames(conserved_markers_ph1), angle = 90, size = 3, raster = T, draw.lines = F)



###
Idents(blood.combined_flt) <- "Subclustering"
ph <- subset(blood.combined_flt, idents = "PH 1")
Idents(ph) <- "origin"

ph_markers <- FindAllMarkers(ph, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
ph_markers <- subset(ph_markers, p_val_adj <= 0.05 & avg_logFC >= 1)
ph_markers <- ph_markers[setdiff(rownames(ph_markers), scGenes), ]
nrow(ph_markers)
ph_markers <- ph_markers[setdiff(rownames(ph_markers), bkGenes), ]##
nrow(ph_markers)


top10 <- ph_markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = ph, features = top10$gene, angle = 90, size = 3, draw.lines = F)
ggsave('compare/ph1/blood.combined_MAST_degs.pdf', units = 'cm', width = 15, height = 10)
write.table(ph_markers, 'compare/ph1/blood.combined_MAST_degs.PH_1.txt', quote = F, sep = '\t', row.names = F, col.names = T)

avg_ph <- log1p(AverageExpression(ph, verbose = FALSE)$RNA)
avg_ph$gene <- rownames(avg_ph)

avg_ph$labelgenes <- NA
avg_ph[c(top10$gene), 'labelgenes'] <- c(top10$gene)
avg_ph[rownames(conserved_markers_ph1), 'labelgenes'] <- rownames(conserved_markers_ph1)


avg_ph$colorgenes <- 'Others'
avg_ph[rownames(subset(ph_markers, avg_logFC >= 1 & cluster == 'LG')), 'colorgenes'] <- 'L'
avg_ph[rownames(subset(ph_markers, avg_logFC >= 1 & cluster == 'Circ')), 'colorgenes'] <- 'C'
avg_ph[rownames(conserved_markers_ph1), 'colorgenes'] <- 'Common'

plt <- ggplot(avg_ph, aes(LG, Circ, label = labelgenes, col = colorgenes)) + 
  geom_point() + scale_color_manual(values = c('#69b9c2', 'black', '#fd9409', 'grey90')) + 
  ggtitle("blood.combined") + xlim(0, 8) + ylim(0, 8) + 
  geom_text_repel(data = avg_ph[c(top10$gene[1:10]), ], color = '#fd9409',
                  nudge_x = 8 - avg_ph[c(top10$gene[1:10]), ]$LG,
                  segment.size = 0.2, segment.color = "grey50", direction = "y", hjust = 1) +
  geom_text_repel(data = avg_ph[top10$gene[11:20], ], 
                  nudge_y = 8 - avg_ph[top10$gene[11:20], ]$Circ,
                  segment.size = 0.2, segment.color = "grey50", direction = "x", vjust = 0, angle = 90) +
  geom_text_repel(data = avg_ph[rownames(conserved_markers_ph1), ],
                  nudge_x = 6 - avg_ph[rownames(conserved_markers_ph1), ]$LG,
                  nudge_y = 6 - avg_ph[rownames(conserved_markers_ph1), ]$Circ,
                  segment.size = 0.2, segment.color = "grey50", direction = "y", hjust = 1) + 
  theme_bw() +
  theme(legend.position = 'None', panel.grid = element_blank(), plot.title = element_blank());plt
ggsave('compare/ph1/blood.combined_scatter.PH_1_v2.pdf', units = 'cm', width = 10, height = 10)

plt <- ggplot() +
  geom_point(data = subset(avg_ph, colorgenes == 'Others'), aes(LG, Circ, col = colorgenes), size = .5) + 
  geom_point(data = subset(avg_ph, colorgenes == 'C'), aes(LG, Circ, col = colorgenes), size = .5) + 
  geom_point(data = subset(avg_ph, colorgenes == 'L'), aes(LG, Circ, col = colorgenes), size = .5) + 
  geom_point(data = subset(avg_ph, colorgenes == 'Common'), aes(LG, Circ, col = colorgenes), size = .5) + 
  scale_color_manual(values = c('#69b9c2', 'black', '#fd9409', 'grey90')) + 
  theme_void() + theme(legend.position = 'None') + xlim(0, 8) + ylim(0, 8);plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('compare/ph1/blood.combined_scatter.PH_1_v2.augment.pdf', units = 'cm', width = 4, height = 4)
ggsave('compare/ph1/blood.combined_scatter.PH_1_v2.augment.png', units = 'cm', width = 4, height = 4)
###


############
### PH 4 ###
############
dir.create('compare/ph4')
head(blood.combined_flt@meta.data)
Idents(blood.combined_flt) <- "Subclustering"
ph <- subset(blood.combined_flt, idents = "PH 4")
Idents(ph) <- "origin"

ph_markers <- FindAllMarkers(ph, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
ph_markers <- subset(ph_markers, p_val_adj <= 0.05 & avg_logFC >= 1)
ph_markers <- ph_markers[setdiff(rownames(ph_markers), scGenes), ]
nrow(ph_markers)
ph_markers <- ph_markers[setdiff(rownames(ph_markers), bkGenes), ]
nrow(ph_markers)


top10 <- ph_markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = ph, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('compare/ph4/blood.combined_flt_MAST_degs.PH4.pdf', units = 'cm', width = 30, height = 20)
write.table(ph_markers, 'compare/ph4/blood.combined_flt_MAST_degs.PH4.txt', quote = F, sep = '\t', row.names = F, col.names = T)

avg_ph <- log1p(AverageExpression(ph, verbose = FALSE)$RNA)
avg_ph$gene <- rownames(avg_ph)

avg_ph$labelgenes <- NA
avg_ph[c(top10$gene), 'labelgenes'] <- c(top10$gene)

avg_ph$colorgenes <- 'Others'
avg_ph[rownames(subset(ph_markers, avg_logFC >= 1 & cluster == 'LG')), 'colorgenes'] <- 'L'
avg_ph[rownames(subset(ph_markers, avg_logFC >= 1 & cluster == 'Circ')), 'colorgenes'] <- 'C'

plt <- ggplot(avg_ph, aes(LG, Circ, label = labelgenes, col = colorgenes)) + 
  geom_point() + scale_color_manual(values = c('#69b9c2', '#fd9409', 'grey90')) + 
  ggtitle("blood.combined_flt") + xlim(0, 8) + ylim(0, 8) + 
  geom_text_repel(data = avg_ph[c(top10$gene[1:7]), ], color = '#fd9409',
                  nudge_x = 8 - avg_ph[c(top10$gene[1:7]), ]$LG,
                  segment.size = 0.2, segment.color = "grey50", direction = "y", hjust = 1) +
  geom_text_repel(data = avg_ph[top10$gene[8:17], ], 
                  nudge_y = 8 - avg_ph[top10$gene[8:17], ]$Circ,
                  segment.size = 0.2, segment.color = "grey50", direction = "x", vjust = 0, angle = 90) +
  theme_bw() +
  theme(legend.position = 'None', panel.grid = element_blank(), plot.title = element_blank());plt
ggsave('compare/ph4/blood.combined_flt_scatter.PH4.pdf', units = 'cm', width = 10, height = 10)


plt <- ggplot() +
  geom_point(data = subset(avg_ph, colorgenes == 'Others'), aes(LG, Circ, col = colorgenes), size = .5) + 
  geom_point(data = subset(avg_ph, colorgenes == 'C'), aes(LG, Circ, col = colorgenes), size = .5) + 
  geom_point(data = subset(avg_ph, colorgenes == 'L'), aes(LG, Circ, col = colorgenes), size = .5) + 
  scale_color_manual(values = c('#69b9c2', '#fd9409', 'grey90')) + 
  theme_void() + theme(legend.position = 'None') + xlim(0, 8) + ylim(0, 8);plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('compare/ph4/blood.combined_flt_scatter.PH4.augment.pdf', units = 'cm', width = 4, height = 4)



##########
### PM ###
##########
dir.create('compare/pm')
head(blood.combined_flt@meta.data)
Idents(blood.combined_flt) <- "anno_simple"
pm <- subset(blood.combined_flt, idents = "PM")
Idents(pm) <- "origin"
summary(subset(pm@meta.data, origin == 'LG')$Subclustering)
summary(subset(pm@meta.data, origin == 'Circ')$Subclustering)

pm_markers <- FindAllMarkers(pm, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
pm_markers <- subset(pm_markers, p_val_adj <= 0.05 & avg_logFC >= 1)
pm_markers <- pm_markers[setdiff(rownames(pm_markers), scGenes), ]
nrow(pm_markers)
pm_markers <- pm_markers[setdiff(rownames(pm_markers), bkGenes), ]###
nrow(pm_markers)


top10 <- pm_markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = pm, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('compare/pm/blood.combined_flt_MAST_degs.PM.pdf', units = 'cm', width = 30, height = 20)
write.table(pm_markers, 'compare/pm/blood.combined_flt_MAST_degs.PM.txt', quote = F, sep = '\t', row.names = F, col.names = T)

avg_pm <- log1p(AverageExpression(pm, verbose = FALSE)$RNA)
avg_pm$gene <- rownames(avg_pm)

avg_pm$labelgenes <- NA
avg_pm[c(top10$gene), 'labelgenes'] <- c(top10$gene)

avg_pm$colorgenes <- 'Others'
avg_pm[rownames(subset(pm_markers, avg_logFC >= 1 & cluster == 'LG')), 'colorgenes'] <- 'L'
avg_pm[rownames(subset(pm_markers, avg_logFC >= 1 & cluster == 'Circ')), 'colorgenes'] <- 'C'

plt <- ggplot(avg_pm, aes(LG, Circ, label = labelgenes, col = colorgenes)) + 
  geom_point() + scale_color_manual(values = c('#69b9c2', '#fd9409', 'grey90')) + 
  ggtitle("blood.combined_flt") + xlim(0, 8) + ylim(0, 8) + 
  geom_text_repel(data = avg_pm[c(top10$gene[1:10]), ], color = '#fd9409',
                  nudge_x = 8 - avg_pm[c(top10$gene[1:10]), ]$LG,
                  segment.size = 0.2, segment.color = "grey50", direction = "y", hjust = 1) +
  geom_text_repel(data = avg_pm[top10$gene[11:20], ], 
                  nudge_y = 8 - avg_pm[top10$gene[11:20], ]$Circ,
                  segment.size = 0.2, segment.color = "grey50", direction = "x", vjust = 0, angle = 90) +
  theme_bw() +
  theme(legend.position = 'None', panel.grid = element_blank(), plot.title = element_blank());plt
ggsave('compare/pm/blood.combined_flt_scatter.PM.pdf', units = 'cm', width = 10, height = 10)

plt <- ggplot() +
  geom_point(data = subset(avg_pm, colorgenes == 'Others'), aes(LG, Circ, col = colorgenes), size = .5) + 
  geom_point(data = subset(avg_pm, colorgenes == 'C'), aes(LG, Circ, col = colorgenes), size = .5) + 
  geom_point(data = subset(avg_pm, colorgenes == 'L'), aes(LG, Circ, col = colorgenes), size = .5) + 
  scale_color_manual(values = c('#69b9c2', '#fd9409', 'grey90')) + 
  theme_void() + theme(legend.position = 'None') + xlim(0, 8) + ylim(0, 8);plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('compare/pm/blood.combined_flt_scatter.PM.augment.pdf', units = 'cm', width = 4, height = 4)
ggsave('compare/pm/blood.combined_flt_scatter.PM.augment.png', units = 'cm', width = 4, height = 4)



############
### PM 1 ###
############
dir.create('compare/pm1')
head(blood.combined_flt@meta.data)
Idents(blood.combined_flt) <- "anno_simple"

pm <- subset(blood.combined_flt, idents = "PM")
pm <- subset(pm, cells = rownames(subset(pm@meta.data, Subclustering == 'PM 1')))
summary(subset(pm@meta.data, origin == 'LG')$Subclustering)
summary(subset(pm@meta.data, origin == 'Circ')$Subclustering)
Idents(pm) <- "origin"

pm_markers <- FindAllMarkers(pm, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
pm_markers <- subset(pm_markers, p_val_adj <= 0.05 & avg_logFC >= 1)
pm_markers <- pm_markers[setdiff(rownames(pm_markers), scGenes), ]
nrow(pm_markers)
pm_markers <- pm_markers[setdiff(rownames(pm_markers), bkGenes), ]###
nrow(pm_markers)


top10 <- pm_markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = pm, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('compare/pm1/blood.combined_flt_MAST_degs.PM_1.pdf', units = 'cm', width = 30, height = 20)
write.table(pm_markers, 'compare/pm1/blood.combined_flt_MAST_degs.PM_1.txt', quote = F, sep = '\t', row.names = F, col.names = T)

avg_pm <- log1p(AverageExpression(pm, verbose = FALSE)$RNA)
avg_pm$gene <- rownames(avg_pm)

avg_pm$labelgenes <- NA
avg_pm[c(top10$gene), 'labelgenes'] <- c(top10$gene)

avg_pm$colorgenes <- 'Others'
avg_pm[rownames(subset(pm_markers, avg_logFC >= 1 & cluster == 'LG')), 'colorgenes'] <- 'L'
avg_pm[rownames(subset(pm_markers, avg_logFC >= 1 & cluster == 'Circ')), 'colorgenes'] <- 'C'

plt <- ggplot(avg_pm, aes(LG, Circ, label = labelgenes, col = colorgenes)) + 
  geom_point() + scale_color_manual(values = c('#69b9c2', '#fd9409', 'grey90')) + 
  ggtitle("blood.combined_flt") + xlim(0, 8) + ylim(0, 8) + 
  geom_text_repel(data = avg_pm[c(top10$gene[1:10]), ], color = '#fd9409',
                  nudge_x = 8 - avg_pm[c(top10$gene[1:10]), ]$LG,
                  segment.size = 0.2, segment.color = "grey50", direction = "y", hjust = 1) +
  geom_text_repel(data = avg_pm[top10$gene[11:20], ], 
                  nudge_y = 8 - avg_pm[top10$gene[11:20], ]$Circ,
                  segment.size = 0.2, segment.color = "grey50", direction = "x", vjust = 0, angle = 90) +
  theme_bw() +
  theme(legend.position = 'None', panel.grid = element_blank(), plot.title = element_blank());plt
ggsave('compare/pm1/blood.combined_flt_scatter.PM_1.pdf', units = 'cm', width = 10, height = 10)

plt <- ggplot() +
  geom_point(data = subset(avg_pm, colorgenes == 'Others'), aes(LG, Circ, col = colorgenes), size = .5) + 
  geom_point(data = subset(avg_pm, colorgenes == 'C'), aes(LG, Circ, col = colorgenes), size = .5) + 
  geom_point(data = subset(avg_pm, colorgenes == 'L'), aes(LG, Circ, col = colorgenes), size = .5) + 
  scale_color_manual(values = c('#69b9c2', '#fd9409', 'grey90')) + 
  theme_void() + theme(legend.position = 'None') + xlim(0, 8) + ylim(0, 8);plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('compare/pm1/blood.combined_flt_scatter.PM_1.augment.pdf', units = 'cm', width = 4, height = 4)
ggsave('compare/pm1/blood.combined_flt_scatter.PM_1.augment.png', units = 'cm', width = 4, height = 4)



##########
### CC ###
##########
dir.create('compare/cc')
Idents(blood.combined_flt) <- "anno_simple"
cc <- subset(blood.combined_flt, idents = "CC")
Idents(cc) <- "origin"

cc_markers <- FindAllMarkers(cc, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
cc_markers <- subset(cc_markers, p_val_adj <= 0.05 & avg_logFC >= 1)
cc_markers <- cc_markers[setdiff(rownames(cc_markers), scGenes), ]
nrow(cc_markers)
cc_markers <- cc_markers[setdiff(rownames(cc_markers), bkGenes), ]
nrow(cc_markers)

top10 <- cc_markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = cc, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('compare/cc/blood.combined_flt_MAST_degs.CC.pdf', units = 'cm', width = 30, height = 20)
write.table(cc_markers, 'compare/cc/blood.combined_flt_MAST_degs.CC.txt', quote = F, sep = '\t', row.names = F, col.names = T)

avg_cc <- log1p(AverageExpression(cc, verbose = FALSE)$RNA)
avg_cc$gene <- rownames(avg_cc)

avg_cc$labelgenes <- NA
avg_cc[c(top10$gene), 'labelgenes'] <- c(top10$gene)

avg_cc$colorgenes <- 'Others'
avg_cc[rownames(subset(cc_markers, avg_logFC >= 1 & cluster == 'LG')), 'colorgenes'] <- 'L'
avg_cc[rownames(subset(cc_markers, avg_logFC >= 1 & cluster == 'Circ')), 'colorgenes'] <- 'C'

plt <- ggplot(avg_cc, aes(LG, Circ, label = labelgenes, col = colorgenes)) + 
  geom_point() + scale_color_manual(values = c('#69b9c2', '#fd9409', 'grey90')) + 
  ggtitle("blood.combined_flt") + xlim(0, 8) + ylim(0, 8) + 
  geom_text_repel(data = avg_cc[c(top10$gene[1:2]), ], color = '#fd9409',
                  nudge_x = 8 - avg_cc[c(top10$gene[1:2]), ]$LG,
                  segment.size = 0.2, segment.color = "grey50", direction = "y", hjust = 1) +
  geom_text_repel(data = avg_cc[top10$gene[3:12], ], 
                  nudge_y = 8 - avg_cc[top10$gene[3:12], ]$Circ,
                  segment.size = 0.2, segment.color = "grey50", direction = "x", vjust = 0, angle = 90) +
  theme_bw() +
  theme(legend.position = 'None', panel.grid = element_blank(), plot.title = element_blank());plt
ggsave('compare/cc/blood.combined_flt_scatter.CC.pdf', units = 'cm', width = 10, height = 10)

plt <- ggplot() +
  geom_point(data = subset(avg_cc, colorgenes == 'Others'), aes(LG, Circ, col = colorgenes), size = .5) + 
  geom_point(data = subset(avg_cc, colorgenes == 'C'), aes(LG, Circ, col = colorgenes), size = .5) + 
  geom_point(data = subset(avg_cc, colorgenes == 'L'), aes(LG, Circ, col = colorgenes), size = .5) + 
  scale_color_manual(values = c('#69b9c2', '#fd9409', 'grey90')) + 
  theme_void() + theme(legend.position = 'None') + xlim(0, 8) + ylim(0, 8);plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('compare/cc/blood.combined_flt_scatter.CC.augment.pdf', units = 'cm', width = 4, height = 4)
ggsave('compare/cc/blood.combined_flt_scatter.CC.augment.png', units = 'cm', width = 4, height = 4)



############
### CC 1 ###
############
dir.create('compare/cc1')
Idents(blood.combined_flt) <- "Subclustering"
cc <- subset(blood.combined_flt, idents = "CC 1")
Idents(cc) <- "origin"

cc_markers <- FindAllMarkers(cc, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
cc_markers <- subset(cc_markers, p_val_adj <= 0.05 & avg_logFC >= 1)
cc_markers <- cc_markers[setdiff(rownames(cc_markers), scGenes), ]
nrow(cc_markers)
cc_markers <- cc_markers[setdiff(rownames(cc_markers), bkGenes), ]
nrow(cc_markers)

top10 <- cc_markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = cc, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('compare/cc1/blood.combined_flt_MAST_degs.CC_1.pdf', units = 'cm', width = 30, height = 20)
write.table(cc_markers, 'compare/cc1/blood.combined_flt_MAST_degs.CC_1.txt', quote = F, sep = '\t', row.names = F, col.names = T)

avg_cc <- log1p(AverageExpression(cc, verbose = FALSE)$RNA)
avg_cc$gene <- rownames(avg_cc)

avg_cc$labelgenes <- NA
avg_cc[c(top10$gene), 'labelgenes'] <- c(top10$gene)

avg_cc$colorgenes <- 'Others'
avg_cc[rownames(subset(cc_markers, avg_logFC >= 1 & cluster == 'LG')), 'colorgenes'] <- 'L'
avg_cc[rownames(subset(cc_markers, avg_logFC >= 1 & cluster == 'Circ')), 'colorgenes'] <- 'C'

plt <- ggplot(avg_cc, aes(LG, Circ, label = labelgenes, col = colorgenes)) + 
  geom_point() + scale_color_manual(values = c('#69b9c2', '#fd9409', 'grey90')) + 
  ggtitle("blood.combined_flt") + xlim(0, 8) + ylim(0, 8) + 
  geom_text_repel(data = avg_cc[c(top10$gene[1:4]), ], color = '#fd9409',
                  nudge_x = 8 - avg_cc[c(top10$gene[1:4]), ]$LG,
                  segment.size = 0.2, segment.color = "grey50", direction = "y", hjust = 1) +
  geom_text_repel(data = avg_cc[top10$gene[5:14], ], 
                  nudge_y = 8 - avg_cc[top10$gene[5:14], ]$Circ,
                  segment.size = 0.2, segment.color = "grey50", direction = "x", vjust = 0, angle = 90) +
  theme_bw() +
  theme(legend.position = 'None', panel.grid = element_blank(), plot.title = element_blank());plt
ggsave('compare/cc1/blood.combined_flt_scatter.CC_1.pdf', units = 'cm', width = 10, height = 10)

plt <- ggplot() +
  geom_point(data = subset(avg_cc, colorgenes == 'Others'), aes(LG, Circ, col = colorgenes), size = .5) + 
  geom_point(data = subset(avg_cc, colorgenes == 'C'), aes(LG, Circ, col = colorgenes), size = .5) + 
  geom_point(data = subset(avg_cc, colorgenes == 'L'), aes(LG, Circ, col = colorgenes), size = .5) + 
  scale_color_manual(values = c('#69b9c2', '#fd9409', 'grey90')) + 
  theme_void() + theme(legend.position = 'None') + xlim(0, 8) + ylim(0, 8);plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('compare/cc1/blood.combined_flt_scatter.CC_1.augment.pdf', units = 'cm', width = 4, height = 4)
ggsave('compare/cc1/blood.combined_flt_scatter.CC_1.augment.png', units = 'cm', width = 4, height = 4)



############
### CC 2 ###
############
dir.create('compare/cc2')
Idents(blood.combined_flt) <- "Subclustering"
cc <- subset(blood.combined_flt, idents = "CC 2")
Idents(cc) <- "origin"

cc_markers <- FindAllMarkers(cc, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
cc_markers <- subset(cc_markers, p_val_adj <= 0.05 & avg_logFC >= 1)
cc_markers <- cc_markers[setdiff(rownames(cc_markers), scGenes), ]
nrow(cc_markers)
cc_markers <- cc_markers[setdiff(rownames(cc_markers), bkGenes), ]
nrow(cc_markers)

top10 <- cc_markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = cc, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('compare/cc2/blood.combined_flt_MAST_degs.CC_2.pdf', units = 'cm', width = 30, height = 20)
write.table(cc_markers, 'compare/cc2/blood.combined_flt_MAST_degs.CC_2.txt', quote = F, sep = '\t', row.names = F, col.names = T)

avg_cc <- log1p(AverageExpression(cc, verbose = FALSE)$RNA)
avg_cc$gene <- rownames(avg_cc)

avg_cc$labelgenes <- NA
avg_cc[c(top10$gene), 'labelgenes'] <- c(top10$gene)

avg_cc$colorgenes <- 'Others'
avg_cc[rownames(subset(cc_markers, avg_logFC >= 1 & cluster == 'LG')), 'colorgenes'] <- 'L'
avg_cc[rownames(subset(cc_markers, avg_logFC >= 1 & cluster == 'Circ')), 'colorgenes'] <- 'C'

plt <- ggplot(avg_cc, aes(LG, Circ, label = labelgenes, col = colorgenes)) + 
  geom_point() + scale_color_manual(values = c('#69b9c2', '#fd9409', 'grey90')) + 
  ggtitle("blood.combined_flt") + xlim(0, 8) + ylim(0, 8) + 
  geom_text_repel(data = avg_cc[c(top10$gene[1:5]), ], color = '#fd9409',
                  nudge_x = 8 - avg_cc[c(top10$gene[1:5]), ]$LG,
                  segment.size = 0.2, segment.color = "grey50", direction = "y", hjust = 1) +
  geom_text_repel(data = avg_cc[top10$gene[6:15], ], 
                  nudge_y = 8 - avg_cc[top10$gene[6:15], ]$Circ,
                  segment.size = 0.2, segment.color = "grey50", direction = "x", vjust = 0, angle = 90) +
  theme_bw() +
  theme(legend.position = 'None', panel.grid = element_blank(), plot.title = element_blank());plt
ggsave('compare/cc2/blood.combined_flt_scatter.CC_2.pdf', units = 'cm', width = 10, height = 10)

plt <- ggplot() +
  geom_point(data = subset(avg_cc, colorgenes == 'Others'), aes(LG, Circ, col = colorgenes), size = .5) + 
  geom_point(data = subset(avg_cc, colorgenes == 'C'), aes(LG, Circ, col = colorgenes), size = .5) + 
  geom_point(data = subset(avg_cc, colorgenes == 'L'), aes(LG, Circ, col = colorgenes), size = .5) + 
  scale_color_manual(values = c('#69b9c2', '#fd9409', 'grey90')) + 
  theme_void() + theme(legend.position = 'None') + xlim(0, 8) + ylim(0, 8);plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('compare/cc2/blood.combined_flt_scatter.CC_2.augment.pdf', units = 'cm', width = 4, height = 4)
ggsave('compare/cc2/blood.combined_flt_scatter.CC_2.augment.png', units = 'cm', width = 4, height = 4)


