#!/home/sangho/miniconda3/envs/seurat3/bin/R
library(Seurat)
library(ggplot2)
library(Matrix)
library(cowplot)
library(plyr)
library(dplyr)

### Circulation inDrop-seq ###
circulation_expr <- read.delim('sudhir.merged.uninjured.expr.symbol.txt', header = T, sep = '\t', check.names = F, row.names = 1)
circulation_label <- read.delim('sudhir.merged.indropseq.mtcut.label.txt', header = T, sep = '\t', check.names = F, row.names = 1)
circulation_label <- circulation_label[, c(4, 1)]
colnames(circulation_label) <- c('Library', 'timepoint')
# uninjured_1
uninjured_1_label <- subset(circulation_label, Library == 'uninjured_1')
uninjured_1_label$timepoint <- 'Circ_indrop_120'
uninjured_1_label$anno_simple <- 'Circ_indrop_120'
uninjured_1_label$Subclustering <- 'Circ_indrop_120'
uninjured_1_expr <- circulation_expr[, rownames(uninjured_1_label)]
uninjured_1_expr <- as(as.matrix(uninjured_1_expr), "dgCMatrix")
# uninjured_2
uninjured_2_label <- subset(circulation_label, Library == 'uninjured_2')
uninjured_2_label$timepoint <- 'Circ_indrop_120'
uninjured_2_label$anno_simple <- 'Circ_indrop_120'
uninjured_2_label$Subclustering <- 'Circ_indrop_120'
uninjured_2_expr <- circulation_expr[, rownames(uninjured_2_label)]
uninjured_2_expr <- as(as.matrix(uninjured_2_expr), "dgCMatrix")


### Circulation 96 Drop-seq ###
circdrop96_expr <- read.delim('merged.exprs.circ.96AEL.IDs.symbol.txt', header = T, sep = '\t', check.names = F, row.names = 1)
circdrop96_label <- read.delim('dropseq.mtcut.label.96AEL.txt', header = T, sep = '\t', check.names = F, row.names = 1)
circdrop96_label <- circdrop96_label[, c(4, 1)]
colnames(circdrop96_label) <- c('Library', 'timepoint')
circdrop96_label$timepoint <- 'Circ_drop_96'
circdrop96_label$anno_simple <- 'Circ_drop_96'
circdrop96_label$Subclustering <- 'Circ_drop_96'
circdrop96_expr <- circdrop96_expr[, rownames(circdrop96_label)]
circdrop96_expr <- as(as.matrix(circdrop96_expr), "dgCMatrix")


### Circulation 120 Drop-seq ###
circdrop120_expr <- read.delim('merged.exprs.circ.120AEL.IDs.symbol.txt', header = T, sep = '\t', check.names = F, row.names = 1)
circdrop120_label <- read.delim('dropseq.mtcut.label.120AEL.txt', header = T, sep = '\t', check.names = F, row.names = 1)
circdrop120_label <- circdrop120_label[, c(4, 1)]
colnames(circdrop120_label) <- c('Library', 'timepoint')
circdrop120_label$timepoint <- 'Circ_drop_120'
circdrop120_label$anno_simple <- 'Circ_drop_120'
circdrop120_label$Subclustering <- 'Circ_drop_120'
circdrop120_expr <- circdrop120_expr[, rownames(circdrop120_label)]
circdrop120_expr <- as(as.matrix(circdrop120_expr), "dgCMatrix")


### lymph gland ###
lymphgland_expr <- read.delim('merged.3tps.expr.allGenes.IDs.symbol.txt', header = T, sep = '\t', check.names = F, row.names = 1)
lymphgland_label <- read.delim('labels.txt', header = T, sep = '\t', check.names = F, row.names = 1)
lymphgland_label <- subset(lymphgland_label, Subclustering != 'DV' & Subclustering != 'RG' & Subclustering != 'Neurons'); lymphgland_label <- droplevels(lymphgland_label)
lymphgland_label <- lymphgland_label[, c(4, 6,  7, 8)] # Library, timepoint, anno_simple, Subclustering
# LG 96
lymphgland96_label <- subset(lymphgland_label, timepoint == 'AEL96hr')
lymphgland96_label$timepoint <- 'LG_drop_96'
lymphgland96_expr <- lymphgland_expr[, rownames(lymphgland96_label)]
lymphgland96_expr <- as(as.matrix(lymphgland96_expr), "dgCMatrix")
# LG 120
lymphgland120_label <- subset(lymphgland_label, timepoint == 'AEL120hr')
lymphgland120_label$timepoint <- 'LG_drop_120'
lymphgland120_expr <- lymphgland_expr[, rownames(lymphgland120_label)]
lymphgland120_expr <- as(as.matrix(lymphgland120_expr), "dgCMatrix")


### Seurat objects for each dataset ###
uninjured_1 <- CreateSeuratObject(counts = uninjured_1_expr, project = "lg2circ")
uninjured_1 <- AddMetaData(object = uninjured_1, metadata = uninjured_1_label)
uninjured_1 <- NormalizeData(uninjured_1, normalization.method = "LogNormalize", scale.factor = 10000)
uninjured_1 <- FindVariableFeatures(object = uninjured_1, selection.method = "vst", nfeatures = 2000)
head(uninjured_1@meta.data)

uninjured_2 <- CreateSeuratObject(counts = uninjured_2_expr, project = "lg2circ")
uninjured_2 <- AddMetaData(object = uninjured_2, metadata = uninjured_2_label)
uninjured_2 <- NormalizeData(uninjured_2, normalization.method = "LogNormalize", scale.factor = 10000)
uninjured_2 <- FindVariableFeatures(object = uninjured_2, selection.method = "vst", nfeatures = 2000)
head(uninjured_2@meta.data)

circdrop96 <- CreateSeuratObject(counts = circdrop96_expr, project = "lg2circ")
circdrop96 <- AddMetaData(object = circdrop96, metadata = circdrop96_label)
circdrop96 <- NormalizeData(circdrop96, normalization.method = "LogNormalize", scale.factor = 10000)
circdrop96 <- FindVariableFeatures(object = circdrop96, selection.method = "vst", nfeatures = 2000)
head(circdrop96@meta.data)

circdrop120 <- CreateSeuratObject(counts = circdrop120_expr, project = "lg2circ")
circdrop120 <- AddMetaData(object = circdrop120, metadata = circdrop120_label)
circdrop120 <- NormalizeData(circdrop120, normalization.method = "LogNormalize", scale.factor = 10000)
circdrop120 <- FindVariableFeatures(object = circdrop120, selection.method = "vst", nfeatures = 2000)
head(circdrop120@meta.data)

lymphgland96 <- CreateSeuratObject(counts = lymphgland96_expr, project = "lg2circ")
lymphgland96 <- AddMetaData(object = lymphgland96, metadata = lymphgland96_label)
lymphgland96 <- NormalizeData(lymphgland96, normalization.method = "LogNormalize", scale.factor = 10000)
lymphgland96 <- FindVariableFeatures(object = lymphgland96, selection.method = "vst", nfeatures = 2000)
head(lymphgland96@meta.data)

lymphgland120 <- CreateSeuratObject(counts = lymphgland120_expr, project = "lg2circ")
lymphgland120 <- AddMetaData(object = lymphgland120, metadata = lymphgland120_label)
lymphgland120 <- NormalizeData(lymphgland120, normalization.method = "LogNormalize", scale.factor = 10000)
lymphgland120 <- FindVariableFeatures(object = lymphgland120, selection.method = "vst", nfeatures = 2000)
head(lymphgland120@meta.data)
save.image('tmp1.Rdata')


### Transfer labels ###
# LG_drop_120 -> Circ_indrop_120-1
transfer.anchors <- FindTransferAnchors(reference = lymphgland120, query = uninjured_1, dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = lymphgland120$Subclustering, dims = 1:30)
uninjured_1@meta.data$Subclustering <- predictions$predicted.id
# LG_drop_120 -> Circ_indrop_120-2
transfer.anchors <- FindTransferAnchors(reference = lymphgland120, query = uninjured_2, dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = lymphgland120$Subclustering, dims = 1:30)
uninjured_2@meta.data$Subclustering <- predictions$predicted.id
# LG_drop_120 -> Circ_drop_120
transfer.anchors <- FindTransferAnchors(reference = lymphgland120, query = circdrop120, dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = lymphgland120$Subclustering, dims = 1:30)
circdrop120@meta.data$Subclustering <- predictions$predicted.id
# LG_drop_96 -> Circ_drop_96
transfer.anchors <- FindTransferAnchors(reference = lymphgland96, query = circdrop96, dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = lymphgland96$Subclustering, dims = 1:30)
circdrop96@meta.data$Subclustering <- predictions$predicted.id
save.image('tmp2.Rdata')


### Anchor datasets ###
blood.anchors <- FindIntegrationAnchors(object.list = list(lymphgland96, lymphgland120, circdrop96, circdrop120, uninjured_1, uninjured_2), dims = 1:30)
blood.combined <- IntegrateData(anchorset = blood.anchors, dims = 1:30)
DefaultAssay(blood.combined) <- 'RNA'
blood.combined[["percent.mt"]] <- PercentageFeatureSet(object = blood.combined, pattern = "^mt:")

plot1 <- FeatureScatter(object = blood.combined, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'Library') + geom_abline(intercept = 10, col = 'red2', slope = 0)
plot2 <- FeatureScatter(object = blood.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'Library')
CombinePlots(plots = list(plot1, plot2))
ggsave('stats/2-2.stats.percent.mt_nFeature_bynCount_RNA.byLibrary.pdf', units = 'cm', width = 26, height = 10)

blood.combined@meta.data$timepoint <- factor(blood.combined@meta.data$timepoint, c('LG_drop_96', 'LG_drop_120', 'Circ_drop_96', 'Circ_drop_120', 'Circ_indrop_120'))
DefaultAssay(object = blood.combined) <- "integrated"
head(blood.combined@meta.data)


### The standard workflow ###
blood.combined <- ScaleData(object = blood.combined, vars.to.regress = c('Library', 'nCount_RNA', 'percent.mt'))
blood.combined <- RunPCA(object = blood.combined, npcs = 80)
blood.combined <- JackStraw(object = blood.combined, num.replicate = 100, dims = 80)
blood.combined <- ScoreJackStraw(object = blood.combined, dims = 1:80)

JackStrawPlot(blood.combined, dims = 1:80) # PCs
ggsave('stats/2-4.pca.JackStrawPlot.pdf', units = 'cm', width = 35, height = 14)
ElbowPlot(blood.combined, ndims = 80) # PCs
ggsave('stats/2-4.pca.ElbowPlot.pdf', units = 'cm', width = 15, height = 10)
save.image('tmp3.Rdata')


### UMAP ###
for (tmpseed in c(1021351:1021399)){
  blood.combined <- RunUMAP(blood.combined, reduction = "pca", dims = 1:54, min.dist = 0.4, n.components = 3, seed.use = tmpseed, umap.method = 'umap-learn', metric = 'correlation')
  
  DimPlot(object = blood.combined, reduction = "umap", group.by = "Subclustering", pt.size = .5, label = T, label.size = 2)
  ggsave(paste(c('umap/compare/umap.1_2.bySubclustering.seed', tmpseed, '.mindist_0.4_.pdf'), collapse = ''), units = 'cm', width = 18, height = 12)
  DimPlot(object = blood.combined, reduction = "umap", dims = c(1,3), group.by = "Subclustering", pt.size = .5, label = T, label.size = 2)
  ggsave(paste(c('umap/compare/umap.1_3.bySubclustering.seed', tmpseed, '.mindist_0.4_.pdf'), collapse = ''), units = 'cm', width = 18, height = 12)
  DimPlot(object = blood.combined, reduction = "umap", dims = c(2,3), group.by = "Subclustering", pt.size = .5, label = T, label.size = 2)
  ggsave(paste(c('umap/compare/umap.2_3.bySubclustering.seed', tmpseed, '.mindist_0.4_.pdf'), collapse = ''), units = 'cm', width = 18, height = 12)
}

DefaultAssay(blood.combined) <- 'RNA'
# 1021290 0.3; 1021367 0.4
blood.combined <- RunUMAP(blood.combined, reduction = "pca", dims = 1:54, min.dist = 0.4, n.components = 3, seed.use = 1021367, umap.method = 'umap-learn', metric = 'correlation')

FeaturePlot(blood.combined, features = c('Antp', 'Dl', 'Ance', 'IM18', 'Hml', 'Ama', 'atilla', 'PPO1', 'dysf'), cols = c("grey","red"))
ggsave('umap/umap.1_2.markers.seed1021367.mindist_0.4.pdf', units = 'cm', width = 27, height = 22)
FeaturePlot(blood.combined, features = c('Antp', 'Dl', 'Ance', 'IM18', 'Hml', 'Ama', 'atilla', 'PPO1', 'dysf'), cols = c("grey","red"), dims = c(1,3))
ggsave('umap/umap.1_3.markers.seed1021367.mindist_0.4.pdf', units = 'cm', width = 27, height = 22)
FeaturePlot(blood.combined, features = c('Antp', 'Dl', 'Ance', 'IM18', 'Hml', 'Ama', 'atilla', 'PPO1', 'dysf'), cols = c("grey","red"), dims = c(2,3))
ggsave('umap/umap.2_3.markers.seed1021367.mindist_0.4.pdf', units = 'cm', width = 27, height = 22)

DimPlot(object = blood.combined, reduction = "umap", group.by = "timepoint", pt.size = .5)
ggsave('umap/umap.1_2.bytimepoint.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 16, height = 12)
DimPlot(object = blood.combined, reduction = "umap", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
ggsave('umap/umap.1_2.bytimepointSep.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 20, height = 12)

DimPlot(object = blood.combined, reduction = "umap", group.by = "Library", pt.size = .5)
ggsave('umap/umap.1_2.byLibrary.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 16, height = 12)
DimPlot(object = blood.combined, reduction = "umap", group.by = "Library", pt.size = .5) + facet_wrap(~Library)
ggsave('umap/umap.1_2.byLibrarySep.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 28, height = 22)
DimPlot(object = blood.combined, reduction = "umap", dims = c(1,3), group.by = "Library", pt.size = .5) + facet_wrap(~Library)
ggsave('umap/umap.1_3.byLibrarySep.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 28, height = 22)
DimPlot(object = blood.combined, reduction = "umap", dims = c(2,3), group.by = "Library", pt.size = .5) + facet_wrap(~Library)
ggsave('umap/umap.2_3.byLibrarySep.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 28, height = 22)

DimPlot(object = blood.combined, reduction = "umap", group.by = "anno_simple", pt.size = .5, label = T, label.size = 2)
ggsave('umap/umap.1_2.byanno_simple.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 16, height = 12)
DimPlot(object = blood.combined, reduction = "umap", dims = c(1,3), group.by = "anno_simple", pt.size = .5, label = T, label.size = 2)
ggsave('umap/umap.1_3.byanno_simple.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 16, height = 12)
DimPlot(object = blood.combined, reduction = "umap", dims = c(2,3), group.by = "anno_simple", pt.size = .5, label = T, label.size = 2)
ggsave('umap/umap.2_3.byanno_simple.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 16, height = 12)

DimPlot(object = blood.combined, reduction = "umap", group.by = "Subclustering", pt.size = .5, label = T, label.size = 2)
ggsave('umap/umap.1_2.bySubclustering.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 18, height = 12)
DimPlot(object = blood.combined, reduction = "umap", dims = c(1,3), group.by = "Subclustering", pt.size = .5, label = T, label.size = 2)
ggsave('umap/umap.1_3.bySubclustering.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 18, height = 12)
DimPlot(object = blood.combined, reduction = "umap", dims = c(2,3), group.by = "Subclustering", pt.size = .5, label = T, label.size = 2)
ggsave('umap/umap.2_3.bySubclustering.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 18, height = 12)

save.image('tmp4.Rdata')
saveRDS(blood.combined, 'blood.combined.Rds')



#####################
head(blood.combined@meta.data)

DimPlot(blood.combined, group.by = 'origin', cols = c('#ffa500', '#7ac5cd'))
#ggsave('umap/umap.1_2.byOrigin.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 14, height = 12)
plt <- DimPlot(blood.combined, group.by = 'origin', cells = rownames(subset(blood.combined@meta.data, origin == 'LG')), cols = c('#ffa500')) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('umap/umap.1_2.byOrigin-LG.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined, group.by = 'origin', cells = rownames(subset(blood.combined@meta.data, origin == 'Circ')), cols = c('#7ac5cd')) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('umap/umap.1_2.byOrigin-Circ.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 4, height = 4)
DimPlot(blood.combined, group.by = 'origin', cols = c('#ffa500', '#7ac5cd'), dims = c(1, 3))
#ggsave('umap/umap.1_3.byOrigin.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 14, height = 12)
DimPlot(blood.combined, group.by = 'origin', cols = c('#ffa500', '#7ac5cd'), dims = c(2, 3))
#ggsave('umap/umap.2_3.byOrigin.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 14, height = 12)

DimPlot(blood.combined, group.by = 'timepoint', cols = c('#ffda92', '#ffa500', '#a8e2d3', '#7ac5cd', '#4f9bab'))
#ggsave('umap/umap.1_2.bytimepoint.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 16, height = 12)
DimPlot(object = blood.combined, reduction = "umap", group.by = "timepoint", cols = c('#ffda92', '#ffa500', '#a8e2d3', '#7ac5cd', '#4f9bab')) + facet_wrap(~timepoint, ncol = 2)
#ggsave('umap/umap.1_2.bytimepointSep.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 16, height = 18)
DimPlot(blood.combined, group.by = 'timepoint', cols = c('#ffda92', '#ffa500', '#a8e2d3', '#7ac5cd', '#4f9bab'), dims = c(1, 3))
#ggsave('umap/umap.1_3.bytimepoint.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 16, height = 12)
DimPlot(blood.combined, group.by = 'timepoint', cols = c('#ffda92', '#ffa500', '#a8e2d3', '#7ac5cd', '#4f9bab'), dims = c(2, 3))
#ggsave('umap/umap.2_3.bytimepoint.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 16, height = 12)

DimPlot(blood.combined, group.by = 'anno_simple', cols = c('#f15fa6', '#207eb3', '#a80d0c', '#f0a142', '#25a9b0', '#a4a4a4', '#1a1a1a'))
#ggsave('umap/umap.1_2.byanno_simple.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 16, height = 12)
plt <- DimPlot(blood.combined, group.by = 'anno_simple', cols = c('#f15fa6', '#207eb3', '#a80d0c', '#f0a142', '#25a9b0', '#a4a4a4', '#1a1a1a')) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('umap/umap.1_2.byanno_simple.seed1021367.mindist_0.4_.augment.pdf', units = 'cm', width = 4, height = 4)
DimPlot(blood.combined, group.by = 'anno_simple', cols = c('#f15fa6', '#207eb3', '#a80d0c', '#f0a142', '#25a9b0', '#a4a4a4', '#1a1a1a'), dims = c(1, 3))
#ggsave('umap/umap.1_3.byanno_simple.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 16, height = 12)
DimPlot(blood.combined, group.by = 'anno_simple', cols = c('#f15fa6', '#207eb3', '#a80d0c', '#f0a142', '#25a9b0', '#a4a4a4', '#1a1a1a'), dims = c(2, 3))
#ggsave('umap/umap.2_3.byanno_simple.seed1021367.mindist_0.4_.pdf', units = 'cm', width = 16, height = 12)



### cell types
head(blood.combined@meta.data)
# PH
plt <- DimPlot(blood.combined, cells = rownames(subset(blood.combined@meta.data, origin == 'LG')), 
               cells.highlight = rownames(subset(blood.combined@meta.data, anno_simple == 'PH')), cols.highlight = '#207eb3') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/umap.1_2.PH-LG.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined, cells = rownames(subset(blood.combined@meta.data, origin == 'Circ')), 
               cells.highlight = rownames(subset(blood.combined@meta.data, anno_simple == 'PH')), cols.highlight = '#207eb3') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/umap.1_2.PH-Circ.pdf', units = 'cm', width = 4, height = 4)

# PM
plt <- DimPlot(blood.combined, cells = rownames(subset(blood.combined@meta.data, origin == 'LG')), 
               cells.highlight = rownames(subset(blood.combined@meta.data, anno_simple == 'PM')), cols.highlight = '#a80d0c') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/umap.1_2.PM-LG.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined, cells = rownames(subset(blood.combined@meta.data, origin == 'Circ')), 
               cells.highlight = rownames(subset(blood.combined@meta.data, anno_simple == 'PM')), cols.highlight = '#a80d0c') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/umap.1_2.PM-Circ.pdf', units = 'cm', width = 4, height = 4)

# LM
plt <- DimPlot(blood.combined, cells = rownames(subset(blood.combined@meta.data, origin == 'LG')), 
               cells.highlight = rownames(subset(blood.combined@meta.data, anno_simple == 'LM')), cols.highlight = '#f0a142') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/umap.1_2.LM-LG.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined, cells = rownames(subset(blood.combined@meta.data, origin == 'Circ')), 
               cells.highlight = rownames(subset(blood.combined@meta.data, anno_simple == 'LM')), cols.highlight = '#f0a142') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/umap.1_2.LM-Circ.pdf', units = 'cm', width = 4, height = 4)

# CC
plt <- DimPlot(blood.combined, cells = rownames(subset(blood.combined@meta.data, origin == 'LG')), 
               cells.highlight = rownames(subset(blood.combined@meta.data, anno_simple == 'CC')), cols.highlight = '#25a9b0') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/umap.1_2.CC-LG.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined, cells = rownames(subset(blood.combined@meta.data, origin == 'Circ')), 
               cells.highlight = rownames(subset(blood.combined@meta.data, anno_simple == 'CC')), cols.highlight = '#25a9b0') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/umap.1_2.CC-Circ.pdf', units = 'cm', width = 4, height = 4)

# GST-rich
plt <- DimPlot(blood.combined, cells = rownames(subset(blood.combined@meta.data, origin == 'LG')), 
               cells.highlight = rownames(subset(blood.combined@meta.data, anno_simple == 'GST-rich')), cols.highlight = '#787878') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/umap.1_2.GST-LG.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined, cells = rownames(subset(blood.combined@meta.data, origin == 'Circ')), 
               cells.highlight = rownames(subset(blood.combined@meta.data, anno_simple == 'GST-rich')), cols.highlight = '#787878') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/umap.1_2.GST-Circ.pdf', units = 'cm', width = 4, height = 4)

# Adipohemocyte
plt <- DimPlot(blood.combined, cells = rownames(subset(blood.combined@meta.data, origin == 'LG')),  
               cells.highlight = rownames(subset(blood.combined@meta.data, anno_simple == 'Adipohemocyte')), cols.highlight = '#1a1a1a') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/umap.1_2.Adipohemocyte-LG.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(blood.combined, cells = rownames(subset(blood.combined@meta.data, origin == 'Circ')), 
               cells.highlight = rownames(subset(blood.combined@meta.data, anno_simple == 'Adipohemocyte')), cols.highlight = '#1a1a1a') +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('umap/highlight/umap.1_2.Adipohemocyte-Circ.pdf', units = 'cm', width = 4, height = 4)

