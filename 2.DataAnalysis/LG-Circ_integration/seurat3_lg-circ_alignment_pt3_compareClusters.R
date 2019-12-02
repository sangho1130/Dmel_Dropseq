library(Seurat)
library(ggplot2)
library(ggrepel)
library(Matrix)
library(cowplot)
library(plyr)
library(dplyr)

blood.combined_flt <- readRDS('blood.combined_flt.Rds')
DefaultAssay(blood.combined_flt)


scGenes <- read.delim('__filter_scGenes__/bulk_sc_pseudo_pt3-scGenes.scGenes_v2.Sym2ID.ID2Sym.txt', header = F)
scGenes <- scGenes$V1
head(scGenes); length(scGenes)

bkGenes <- read.delim('__filter_scGenes__/bulk_sc_pseudo_pt3-scGenes.bkGenes_v2.Sym2ID.ID2Sym.txt', header = F)
bkGenes <- bkGenes$V1
head(bkGenes); length(bkGenes)


##############
### Origin ###
##############
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
#ggsave('Compare_flt/origin/blood.combined_flt_MAST_degs.pdf', units = 'cm', width = 20, height = 20)
#write.table(blood.combined_flt_markers, 'Compare_flt/origin/blood.combined_flt_MAST_degs.origin.txt', quote = F, sep = '\t', row.names = F, col.names = T)

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
#ggsave('Compare_flt/origin/blood.combined_flt_scatter.origin.pdf', units = 'cm', width = 10, height = 10)

plt <- ggplot(avg_blood.combined_flt, aes(LG, Circ, col = colorgenes)) + 
  geom_point() + scale_color_manual(values = c('#69b9c2', '#fd9409', 'grey90')) + 
  theme_void() + theme(legend.position = 'None') + xlim(0, 8) + ylim(0, 8);plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('Compare_flt/origin/blood.combined_flt_scatter.origin.augment.pdf', units = 'cm', width = 4, height = 4)


##########
### PH ###
##########
head(blood.combined_flt@meta.data)
Idents(blood.combined_flt) <- "anno_simple"
ph <- subset(blood.combined_flt, idents = "PH")
Idents(ph) <- "origin"

ph_markers <- FindAllMarkers(ph, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
ph_markers <- subset(ph_markers, p_val_adj <= 0.05 & avg_logFC >= 1)
ph_markers <- ph_markers[setdiff(rownames(ph_markers), scGenes), ]
nrow(ph_markers)
ph_markers <- ph_markers[setdiff(rownames(ph_markers), c), ]##
nrow(ph_markers)

top10 <- ph_markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = ph, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
#ggsave('Compare_flt/PH/blood.combined_flt_MAST_degs.PH.pdf', units = 'cm', width = 30, height = 20)
#write.table(ph_markers, 'Compare_flt/PH/blood.combined_flt_MAST_degs.PH.txt', quote = F, sep = '\t', row.names = F, col.names = T)

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
#ggsave('Compare_flt/PH/blood.combined_flt_scatter.PH.pdf', units = 'cm', width = 10, height = 10)

plt <- ggplot(avg_ph, aes(LG, Circ, col = colorgenes)) + 
  geom_point() + scale_color_manual(values = c('#69b9c2', '#fd9409', 'grey90')) + 
  theme_void() + theme(legend.position = 'None') + xlim(0, 8) + ylim(0, 8);plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('Compare_flt/PH/blood.combined_flt_scatter.PH.augment.pdf', units = 'cm', width = 4, height = 4)


############
### PH 1 ###
############
levels(Idents(blood.combined_flt))

conserved_markers_ph1 <- FindConservedMarkers(blood.combined_flt, ident.1 = 'PH 1', grouping.var = 'origin', min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
conserved_markers_ph1
conserved_markers_ph1 <- subset(conserved_markers_ph1, conserved_markers_ph1[,5] <= 0.05 &  conserved_markers_ph1[,10] <= 0.05 &
                                  conserved_markers_ph1[,2] >= 1 & conserved_markers_ph1[,7] >= 1)
conserved_markers_ph1 <- conserved_markers_ph1[setdiff(rownames(conserved_markers_ph1), scGenes), ]
nrow(conserved_markers_ph1)
#saveRDS(conserved_markers_ph1, 'degs/conserved_markers_ph1_vs_PHs.Rds')
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
#ggsave('Compare_flt/PH_1/blood.combined_MAST_degs.pdf', units = 'cm', width = 15, height = 10)
#write.table(ph_markers, 'Compare_flt/PH_1/blood.combined_MAST_degs.PH_1.txt', quote = F, sep = '\t', row.names = F, col.names = T)

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
#ggsave('Compare_flt/PH_1/blood.combined_scatter.PH_1_v2.pdf', units = 'cm', width = 10, height = 10)

plt <- ggplot(avg_ph, aes(LG, Circ, col = colorgenes)) + 
  geom_point() + scale_color_manual(values = c('#69b9c2', 'black', '#fd9409', 'grey90')) + 
  theme_void() + theme(legend.position = 'None') + xlim(0, 8) + ylim(0, 8);plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('Compare_flt/PH_1/blood.combined_scatter.PH_1_v2.augment.pdf', units = 'cm', width = 4, height = 4)
###



##########
### PM ###
##########
head(blood.combined_flt@meta.data)
Idents(blood.combined_flt) <- "anno_simple"
pm <- subset(blood.combined_flt, idents = "PM")
Idents(pm) <- "origin"

pm_markers <- FindAllMarkers(pm, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
pm_markers <- subset(pm_markers, p_val_adj <= 0.05 & avg_logFC >= 1)
pm_markers <- pm_markers[setdiff(rownames(pm_markers), scGenes), ]
nrow(pm_markers)
pm_markers <- pm_markers[setdiff(rownames(pm_markers), bkGenes), ]###
nrow(pm_markers)


top10 <- pm_markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(object = pm, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
#ggsave('Compare_flt/PM/blood.combined_flt_MAST_degs.PM.pdf', units = 'cm', width = 30, height = 20)
#write.table(pm_markers, 'Compare_flt/PM/blood.combined_flt_MAST_degs.PM.txt', quote = F, sep = '\t', row.names = F, col.names = T)

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
#ggsave('Compare_flt/PM/blood.combined_flt_scatter.PM.pdf', units = 'cm', width = 10, height = 10)

plt <- ggplot(avg_pm, aes(LG, Circ, col = colorgenes)) + 
  geom_point() + scale_color_manual(values = c('#69b9c2', '#fd9409', 'grey90')) + 
  theme_void() + theme(legend.position = 'None') + xlim(0, 8) + ylim(0, 8);plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('Compare_flt/PM/blood.combined_flt_scatter.PM.augment.pdf', units = 'cm', width = 4, height = 4)


##############
### PM 4-7 ###
##############
head(blood.combined_flt@meta.data)
Idents(blood.combined_flt) <- "anno_simple"
pm <- subset(blood.combined_flt, idents = "PM")
pm <- subset(pm, cells = rownames(subset(pm@meta.data, Subclustering == 'PM 4'|Subclustering == 'PM 5'|Subclustering == 'PM 6'|Subclustering == 'PM 7')))
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
#ggsave('Compare_flt/PM_4-7/blood.combined_flt_MAST_degs.PM_4-7.pdf', units = 'cm', width = 30, height = 20)
#write.table(pm_markers, 'Compare_flt/PM_4-7/blood.combined_flt_MAST_degs.PM_4-7.txt', quote = F, sep = '\t', row.names = F, col.names = T)

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
  geom_text_repel(data = avg_pm[c(top10$gene[1:6]), ], color = '#fd9409',
                  nudge_x = 8 - avg_pm[c(top10$gene[1:6]), ]$LG,
                  segment.size = 0.2, segment.color = "grey50", direction = "y", hjust = 1) +
  geom_text_repel(data = avg_pm[top10$gene[7:16], ], 
                  nudge_y = 8 - avg_pm[top10$gene[7:16], ]$Circ,
                  segment.size = 0.2, segment.color = "grey50", direction = "x", vjust = 0, angle = 90) +
  theme_bw() +
  theme(legend.position = 'None', panel.grid = element_blank(), plot.title = element_blank());plt
#ggsave('Compare_flt/PM_4-7/blood.combined_flt_scatter.PM_4-7.pdf', units = 'cm', width = 10, height = 10)

plt <- ggplot(avg_pm, aes(LG, Circ, col = colorgenes)) + 
  geom_point() + scale_color_manual(values = c('#69b9c2', '#fd9409', 'grey90')) + 
  theme_void() + theme(legend.position = 'None') + xlim(0, 8) + ylim(0, 8);plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('Compare_flt/PM_4-7/blood.combined_flt_scatter.PM_4-7.augment.pdf', units = 'cm', width = 4, height = 4)



summary(as.factor(subset(blood.combined_flt@meta.data, origin == 'LG')$Subclustering))
summary(as.factor(subset(blood.combined_flt@meta.data, origin == 'Circ')$Subclustering))
summary(as.factor(subset(blood.combined_flt@meta.data, origin == 'LG')$anno_simple))
summary(as.factor(subset(blood.combined_flt@meta.data, origin == 'Circ')$anno_simple))


VlnPlot(blood.combined_flt, features = c('Hml', 'Pxn', 'NimC1', 'vkg', 'vir-1', 'stg'), group.by = 'Subclustering', pt.size = 0, ncol = 2)


###
head(blood.combined_flt@meta.data)
Idents(blood.combined_flt) <- "anno_simple"
cc <- subset(blood.combined_flt, idents = "CC")
Idents(cc) <- "origin"

cc_markers <- FindAllMarkers(cc, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
cc_markers <- subset(cc_markers, p_val_adj <= 0.05 & avg_logFC >= 1)
cc_markers <- cc_markers[setdiff(rownames(cc_markers), scGenes), ]
nrow(cc_markers)
cc_markers <- cc_markers[setdiff(rownames(cc_markers), bkGenes), ]###
nrow(cc_markers)

#
head(blood.combined_flt@meta.data)
Idents(blood.combined_flt) <- "Subclustering"
cc <- subset(blood.combined_flt, idents = "CC 1")
Idents(cc) <- "origin"

cc_markers <- FindAllMarkers(cc, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
cc_markers <- subset(cc_markers, p_val_adj <= 0.05 & avg_logFC >= 1)
cc_markers <- cc_markers[setdiff(rownames(cc_markers), scGenes), ]
nrow(cc_markers)
cc_markers <- cc_markers[setdiff(rownames(cc_markers), bkGenes), ]###
nrow(cc_markers)

#
head(blood.combined_flt@meta.data)
Idents(blood.combined_flt) <- "Subclustering"
cc <- subset(blood.combined_flt, idents = "CC 2")
Idents(cc) <- "origin"

cc_markers <- FindAllMarkers(cc, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
cc_markers <- subset(cc_markers, p_val_adj <= 0.05 & avg_logFC >= 1)
cc_markers <- cc_markers[setdiff(rownames(cc_markers), scGenes), ]
nrow(cc_markers)
cc_markers <- cc_markers[setdiff(rownames(cc_markers), bkGenes), ]###
nrow(cc_markers)


DimPlot(blood.combined_flt, group.by = 'origin', cols = c('#ffa500', '#7ac5cd'))
#ggsave('umap/__flt__umap.1_2.byOrigin.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 14, height = 12)
plt <- DimPlot(blood.combined_flt, group.by = 'origin', cells = rownames(subset(blood.combined_flt@meta.data, origin == 'LG')), cols = c('#ffa500')) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('umap/__flt__umap.1_2.byOrigin-LG.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined_flt, group.by = 'origin', cells = rownames(subset(blood.combined_flt@meta.data, origin == 'Circ')), cols = c('#7ac5cd')) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('umap/__flt__umap.1_2.byOrigin-Circ.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 4, height = 4)
DimPlot(blood.combined_flt, group.by = 'origin', cols = c('#ffa500', '#7ac5cd'), dims = c(1, 3))
#ggsave('umap/__flt__umap.1_3.byOrigin.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 14, height = 12)
DimPlot(blood.combined_flt, group.by = 'origin', cols = c('#ffa500', '#7ac5cd'), dims = c(2, 3))
#ggsave('umap/__flt__umap.2_3.byOrigin.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 14, height = 12)

DimPlot(blood.combined_flt, group.by = 'timepoint', cols = c('#ffda92', '#ffa500', '#a8e2d3', '#7ac5cd', '#4f9bab'))
#ggsave('umap/__flt__umap.1_2.bytimepoint.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 16, height = 12)
DimPlot(object = blood.combined_flt, reduction = "umap", group.by = "timepoint", cols = c('#ffda92', '#ffa500', '#a8e2d3', '#7ac5cd', '#4f9bab')) + facet_wrap(~timepoint, ncol = 2)
#ggsave('umap/__flt__umap.1_2.bytimepointSep.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 16, height = 18)
DimPlot(blood.combined_flt, group.by = 'timepoint', cols = c('#ffda92', '#ffa500', '#a8e2d3', '#7ac5cd', '#4f9bab'), dims = c(1, 3))
#ggsave('umap/__flt__umap.1_3.bytimepoint.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 16, height = 12)
DimPlot(blood.combined_flt, group.by = 'timepoint', cols = c('#ffda92', '#ffa500', '#a8e2d3', '#7ac5cd', '#4f9bab'), dims = c(2, 3))
#ggsave('umap/__flt__umap.2_3.bytimepoint.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 16, height = 12)

DimPlot(blood.combined_flt, group.by = 'anno_simple', cols = c('#f15fa6', '#207eb3', '#a80d0c', '#f0a142', '#25a9b0', '#a4a4a4', '#1a1a1a'))
#ggsave('umap/__flt__umap.1_2.byanno_simple.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 16, height = 12)
plt <- DimPlot(blood.combined_flt, group.by = 'anno_simple', cols = c('#f15fa6', '#207eb3', '#a80d0c', '#f0a142', '#25a9b0', '#a4a4a4', '#1a1a1a')) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('umap/__flt__umap.1_2.byanno_simple.seed1021367.mindist_0.4_.augment.pdf', units = 'cm', width = 4, height = 4)
DimPlot(blood.combined_flt, group.by = 'anno_simple', cols = c('#f15fa6', '#207eb3', '#a80d0c', '#f0a142', '#25a9b0', '#a4a4a4', '#1a1a1a'), dims = c(1, 3))
#ggsave('umap/__flt__umap.1_3.byanno_simple.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 16, height = 12)
DimPlot(blood.combined_flt, group.by = 'anno_simple', cols = c('#f15fa6', '#207eb3', '#a80d0c', '#f0a142', '#25a9b0', '#a4a4a4', '#1a1a1a'), dims = c(2, 3))
#ggsave('umap/__flt__umap.2_3.byanno_simple.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 16, height = 12)


### cell types
head(blood.combined_flt@meta.data)
# PH
plt <- DimPlot(blood.combined_flt, cells = rownames(subset(blood.combined_flt@meta.data, origin == 'LG')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, anno_simple == 'PH')), cols.highlight = '#207eb3') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.PH-LG.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined_flt, cells = rownames(subset(blood.combined_flt@meta.data, origin == 'Circ')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, anno_simple == 'PH')), cols.highlight = '#207eb3') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.PH-Circ.pdf', units = 'cm', width = 4, height = 4)

head(blood.combined_flt@meta.data)
# PH 1
plt <- DimPlot(blood.combined_flt, cells = rownames(subset(blood.combined_flt@meta.data, origin == 'LG')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, Subclustering == 'PH 1')), cols.highlight = '#207eb3') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.PH1-LG.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined_flt, cells = rownames(subset(blood.combined_flt@meta.data, origin == 'Circ')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, Subclustering == 'PH 1')), cols.highlight = '#207eb3') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.PH1-Circ.pdf', units = 'cm', width = 4, height = 4)

# PM
plt <- DimPlot(blood.combined_flt, cells = rownames(subset(blood.combined_flt@meta.data, origin == 'LG')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, anno_simple == 'PM')), cols.highlight = '#a80d0c') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.PM-LG.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined_flt, cells = rownames(subset(blood.combined_flt@meta.data, origin == 'Circ')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, anno_simple == 'PM')), cols.highlight = '#a80d0c') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.PM-Circ.pdf', units = 'cm', width = 4, height = 4)

# CC
plt <- DimPlot(blood.combined_flt, cells = rownames(subset(blood.combined_flt@meta.data, origin == 'LG')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, anno_simple == 'CC')), cols.highlight = '#25a9b0') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.CC-LG.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined_flt, cells = rownames(subset(blood.combined_flt@meta.data, origin == 'Circ')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, anno_simple == 'CC')), cols.highlight = '#25a9b0') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.CC-Circ.pdf', units = 'cm', width = 4, height = 4)

# CC 1
plt <- DimPlot(blood.combined_flt, cells = rownames(subset(blood.combined_flt@meta.data, origin == 'LG')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, Subclustering == 'CC 1')), cols.highlight = '#72c4b2') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.CC1-LG.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined_flt, cells = rownames(subset(blood.combined_flt@meta.data, origin == 'Circ')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, Subclustering == 'CC 1')), cols.highlight = '#72c4b2') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.CC1-Circ.pdf', units = 'cm', width = 4, height = 4)

# CC 2
plt <- DimPlot(blood.combined_flt, cells = rownames(subset(blood.combined_flt@meta.data, origin == 'LG')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, Subclustering == 'CC 2')), cols.highlight = '#1b7e7d') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.CC2-LG.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined_flt, cells = rownames(subset(blood.combined_flt@meta.data, origin == 'Circ')), 
               cells.highlight = rownames(subset(blood.combined_flt@meta.data, Subclustering == 'CC 2')), cols.highlight = '#1b7e7d') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/__flt__umap.1_2.CC2-Circ.pdf', units = 'cm', width = 4, height = 4)


#########
#########
Idents(blood.combined_flt) <- factor(Idents(blood.combined_flt), levels = rev(levels(Idents(blood.combined_flt))))
DotPlot(blood.combined_flt, features = rev(c('Ance', 'Tep4', 'Dl', 'N', 'Stat92E', 'dome', 'E(spl)m3-HLH', 'nw', 'lncRNA:cherub', 'scrt')), split.by = 'origin', cols = c('#ffa500', '#7ac5cd')) + 
  labs(x = '', y = '') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave('Compare_flt/PH_1/knownmarkers.pdf', units = 'cm', width = 14, height = 18)

