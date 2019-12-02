#!/home/sangho/miniconda3/envs/seurat3/bin/R
library(Seurat)
library(ggplot2)
library(Matrix)
library(cowplot)


### 72hr ###
lymphgland_72_expr <- read.delim('72hAEL_merged.expr.allGenes.txt', header = T, sep = '\t', check.names = F, row.names = 1)
lymphgland_72_label <- read.delim('72hAEL_merged.label.allGenes.txt', header = T, sep = '\t', check.names = F, row.names = 1)

lymphgland_72_label <- lymphgland_72_label[, c(1:5, ncol(lymphgland_72_label))]
colnames(lymphgland_72_label)[6] <- 'orig.anno'
lymphgland_72_expr <- lymphgland_72_expr[, rownames(lymphgland_72_label)]


### 96hr ###
lymphgland_96_expr <- read.delim('96hAEL_merged.expr.allGenes.txt', header = T, sep = '\t', check.names = F, row.names = 1)
lymphgland_96_label <- read.delim('96hAEL_merged.label.allGenes.txt', header = T, sep = '\t', check.names = F, row.names = 1)

lymphgland_96_label <- lymphgland_96_label[, c(1:5, ncol(lymphgland_96_label))]
colnames(lymphgland_96_label)[6] <- 'orig.anno'
lymphgland_96_expr <- lymphgland_96_expr[, rownames(lymphgland_96_label)]


### 120hr ###
lymphgland_120_expr <- read.delim('120hAEL_merged.expr.allGenes.txt', header = T, sep = '\t', check.names = F, row.names = 1)
lymphgland_120_label <- read.delim('120hAEL_merged.label.allGenes.txt', header = T, sep = '\t', check.names = F, row.names = 1)

lymphgland_120_label <- lymphgland_120_label[, c(1:5, ncol(lymphgland_120_label))]
colnames(lymphgland_120_label)[6] <- 'orig.anno'
lymphgland_120_expr <- lymphgland_120_expr[, rownames(lymphgland_120_label)]


lymphgland_72_expr <- as(as.matrix(lymphgland_72_expr), "dgCMatrix")
lymphgland_96_expr <- as(as.matrix(lymphgland_96_expr), "dgCMatrix")
lymphgland_120_expr <- as(as.matrix(lymphgland_120_expr), "dgCMatrix")


### Seurat objects for each dataset ###
lymphgland_72 <- CreateSeuratObject(counts = lymphgland_72_expr, project = "AEL72hr")
lymphgland_72 <- AddMetaData(object = lymphgland_72, metadata = lymphgland_72_label)
lymphgland_72 <- NormalizeData(lymphgland_72, normalization.method = "LogNormalize", scale.factor = 10000)
lymphgland_72 <- FindVariableFeatures(object = lymphgland_72, selection.method = "vst", nfeatures = 2000)
lymphgland_72@meta.data$timepoint <- "AEL72hr"
head(lymphgland_72@meta.data)

lymphgland_96 <- CreateSeuratObject(counts = lymphgland_96_expr, project = "AEL72hr")
lymphgland_96 <- AddMetaData(object = lymphgland_96, metadata = lymphgland_96_label)
lymphgland_96 <- NormalizeData(lymphgland_96, normalization.method = "LogNormalize", scale.factor = 10000)
lymphgland_96 <- FindVariableFeatures(object = lymphgland_96, selection.method = "vst", nfeatures = 2000)
lymphgland_96@meta.data$timepoint <- "AEL96hr"
head(lymphgland_96@meta.data)

lymphgland_120 <- CreateSeuratObject(counts = lymphgland_120_expr, project = "AEL72hr")
lymphgland_120 <- AddMetaData(object = lymphgland_120, metadata = lymphgland_120_label)
lymphgland_120 <- NormalizeData(lymphgland_120, normalization.method = "LogNormalize", scale.factor = 10000)
lymphgland_120 <- FindVariableFeatures(object = lymphgland_120, selection.method = "vst", nfeatures = 2000)
lymphgland_120@meta.data$timepoint <- "AEL120hr"
head(lymphgland_120@meta.data)


lymphgland.anchors <- FindIntegrationAnchors(object.list = list(lymphgland_72, lymphgland_96, lymphgland_120), dims = 1:30)
lymphgland.combined <- IntegrateData(anchorset = lymphgland.anchors, dims = 1:30)
head(lymphgland.combined@meta.data)
DefaultAssay(object = lymphgland.combined) <- "integrated"
lymphgland.combined@meta.data$timepoint <- factor(lymphgland.combined@meta.data$timepoint, c('AEL72hr', 'AEL96hr', 'AEL120hr'))


### The standard workflow ###
lymphgland.combined <- ScaleData(object = lymphgland.combined, vars.to.regress = c('Library','nCount_RNA'))
lymphgland.combined <- RunPCA(object = lymphgland.combined, npcs = 80)

lymphgland.combined <- JackStraw(object = lymphgland.combined, num.replicate = 100, dims = 80)
lymphgland.combined <- ScoreJackStraw(object = lymphgland.combined, dims = 1:80)

dir.create('stats')
JackStrawPlot(lymphgland.combined, dims = 1:80) #PC52
ggsave('stats/2-4.pca.JackStrawPlot.mitoCut.pdf', units = 'cm', width = 30, height = 10)
ElbowPlot(lymphgland.combined, ndims = 80) #PC52
ggsave('stats/2-4.pca.ElbowPlot.mitoCut.pdf', units = 'cm', width = 15, height = 10)


### UMAP ###
for (tmpseed in c(618201:618250)){
  lymphgland.combined <- RunUMAP(lymphgland.combined, reduction = "pca", dims = 1:52, min.dist = 0.3, n.components = 3, seed.use = tmpseed)
  FeaturePlot(lymphgland.combined, features = c("Antp", "Dl", "Ance", "IM18", "Hml", "Ama", "mthl4", "PPO1", "Mlc1"), cols = c("grey","red"))
  ggsave(paste(c('umap/compare/umap.mitoCut.markers.mindist_0.3.seed', tmpseed, '.pdf'), collapse = ''), units = 'cm', width = 27, height = 22)
  DimPlot(object = lymphgland.combined, reduction = "umap", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
  ggsave(paste(c('umap/compare/umap.mitoCut.bySampleSep.seed', tmpseed, '.mindist_0.3_.pdf'), collapse = ''), units = 'cm', width = 36, height = 12)
}
save.image('tmp3.Rdata')


load('tmp3.Rdata')
lymphgland.combined <- RunUMAP(lymphgland.combined, reduction = "pca", dims = 1:52, min.dist = 0.3, n.components = 3, seed.use = 618213)
FeaturePlot(lymphgland.combined, features = c("Antp", "Dl", "IM18", "Ance", "Hml", "Ama", "mthl4", "PPO1", "Mlc1"), cols = c("grey","red"))
ggsave('umap/umap.mitoCut.markers.mindist_0.3.seed618213.pdf', units = 'cm', width = 27, height = 22)
DimPlot(object = lymphgland.combined, reduction = "umap", group.by = "timepoint", pt.size = .5)
ggsave('umap/umap.mitoCut.bySample.seed618213.mindist_0.3_.pdf', units = 'cm', width = 14, height = 12)
DimPlot(object = lymphgland.combined, reduction = "umap", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
ggsave('umap/umap.mitoCut.bySampleSep.seed618213.mindist0.3_.pdf', units = 'cm', width = 36, height = 12)
save.image('tmp3.Rdata')


### t-SNE
for (tmpseed in c(618151:618199)){
  lymphgland.combined <- RunTSNE(object = lymphgland.combined, dims = 1:52, reduction.key = 'tSNE', dim.embed = 3, seed.use = tmpseed)
  FeaturePlot(lymphgland.combined, features = c("Antp", "Dl", "Ance", "IM18", "Hml", "Ama", "mthl4", "PPO1", "Mlc1"), cols = c("grey","red"), reduction = "tsne")
  ggsave(paste(c('tsne/compare/tsne.mitoCut.1_2.markers.seed', tmpseed, '.pdf'), collapse = ''), units = 'cm', width = 27, height = 22)
  DimPlot(object = lymphgland.combined, reduction = "tsne", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
  ggsave(paste(c('tsne/compare/tsne.mitoCut.1_2.bySampleSep.seed', tmpseed, '.pdf'), collapse = ''), units = 'cm', width = 36, height = 12)
}
save.image('tmp3.Rdata')

load('tmp3.Rdata')
lymphgland.combined <- RunTSNE(object = lymphgland.combined, dims = 1:52, reduction.key = 'tSNE', dim.embed = 3, seed.use = 618139)
FeaturePlot(lymphgland.combined, features = c("Antp", "Dl", "IM18", "Ance", "Hml", "Ama", "mthl4", "PPO1", "Mlc1"), cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.mitoCut.1_2.markers.seed618139.pdf', units = 'cm', width = 27, height = 22)
DimPlot(object = lymphgland.combined, reduction = "tsne", group.by = "timepoint", pt.size = .5)
ggsave('tsne/tsne.mitoCut.1_2.bySample.seed618139.pdf', units = 'cm', units = 'cm', width = 14, height = 12)
DimPlot(object = lymphgland.combined, reduction = "tsne", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
ggsave('tsne/tsne.mitoCut.1_2.bySampleSep.seed618139.pdf', units = 'cm', width = 36, height = 12)
save.image('tmp3.Rdata')


### res ###
lymphgland.combined <- FindNeighbors(object = lymphgland.combined, dims = 1:52)
for (res in c(3:20)){
  res <- res/10
  lymphgland.combined <- FindClusters(object = lymphgland.combined, resolution = res)
  DimPlot(lymphgland.combined, label = T, pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed618139.mitoCut.1_2.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 13.5, height = 10)
  DimPlot(lymphgland.combined, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed618139.mitoCut.1_3.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 13.5, height = 10)
  DimPlot(lymphgland.combined, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed618139.mitoCut.2_3.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 13.5, height = 10)
}

lymphgland.combined <- FindClusters(object = lymphgland.combined, resolution = 0.8)
DimPlot(lymphgland.combined, label = T, pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed618139.mitoCut.1_2.res_0.8.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(lymphgland.combined, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed618139.mitoCut.1_3.res_0.8.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(lymphgland.combined, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed618139.mitoCut.2_3.res_0.8.pdf', units = 'cm', width = 12.5, height = 10)

levels(lymphgland.combined@meta.data$seurat_clusters)
seurat_clusters_anno_marker <- c('PH', 'PM', 'PM', 'PM', 'PM', 'PM', 'PH', 'PH', 'PH', 'PM', 'PH', 'LM', 'Adipohemocyte', 'CC', 'PM', 'GST-rich', 'PSC', 'PH', 'DV')
names(x = seurat_clusters_anno_marker) <- levels(x = lymphgland.combined)
lymphgland.combined <- RenameIdents(lymphgland.combined, seurat_clusters_anno_marker)

saveRDS(lymphgland.combined, 'lymphgland.combined.Rds')
save.image('tmp4.Rdata')


