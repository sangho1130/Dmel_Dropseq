### Subclustering - PSC, PH 1, PH 2 ###
library(Seurat)
library(ggplot2)
library(Matrix)
library(cowplot)
library(plyr)

lymphgland.combined.flt <- readRDS('lymphgland.combined.flt.Rds')

PSC_PHs <- subset(lymphgland.combined.flt, cells = rownames(subset(lymphgland.combined.flt@meta.data, Subclustering == 'PSC' | Subclustering == 'PH 1' | Subclustering == 'PH 2')))
PSC_PHs@meta.data <- droplevels(PSC_PHs@meta.data)

PSC_PHs <- ScaleData(object = PSC_PHs, vars.to.regress = c('Library', 'nCount_RNA'))
PSC_PHs <- RunPCA(object = PSC_PHs, npcs = 30)
PSC_PHs <- JackStraw(object = PSC_PHs, num.replicate = 100, dims = 30)
PSC_PHs <- ScoreJackStraw(object = PSC_PHs, dims = 1:30)

dir.create('stats')
JackStrawPlot(PSC_PHs, dims = 1:30) #PC 21
ggsave('stats/2-4.pca.JackStrawPlot.mitoCut.pdf', units = 'cm', width = 20, height = 10)
ElbowPlot(PSC_PHs, ndims = 30) #PC 21
ggsave('stats/2-4.pca.ElbowPlot.mitoCut.pdf', units = 'cm', width = 15, height = 10)


### t-SNE ###
for (tmpseed in c(704201:704230)){
  PSC_PHs <- RunTSNE(object = PSC_PHs, dims = 1:21, reduction.key = 'tSNE', dim.embed = 3, seed.use = tmpseed)
  FeaturePlot(PSC_PHs, features = c("Antp", "N", "Dl", "fne", "NimB3", "IM18", "Ance", "phm", "Mlc1"), cols = c("grey","red"), reduction = "tsne")
  ggsave(paste(c('tsne/compare/tsne.mitoCut.1_2.markers.seed', tmpseed, '.pdf'), collapse = ''), units = 'cm', width = 27, height = 22)
  DimPlot(object = PSC_PHs, reduction = "tsne", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
  ggsave(paste(c('tsne/compare/tsne.mitoCut.1_2.bySampleSep.seed', tmpseed, '.pdf'), collapse = ''), units = 'cm', width = 36, height = 12)
}

PSC_PHs <- RunTSNE(object = PSC_PHs, dims = 1:21, reduction.key = 'tSNE', dim.embed = 3, seed.use = 704228)
FeaturePlot(PSC_PHs, features = c("Antp", "N", "Dl", "fne", "NimB3", "IM18", "Ance", "phm", "Mlc1"), cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.mitoCut.1_2.markers.seed704228.pdf', units = 'cm', width = 28, height = 22)
DimPlot(PSC_PHs, reduction = "tsne", group.by = "timepoint", pt.size = 1) + facet_wrap(~timepoint)
ggsave('tsne/tsne.mitoCut.1_2.bySampleSep.seed704228.pdf', units = 'cm', width = 24, height = 8)


### UMAP ###
for (tmpseed in c(704251:704270)){
  PSC_PHs <- RunUMAP(PSC_PHs, reduction = "pca", dims = 1:21, min.dist = 0.4, n.components = 3, seed.use = tmpseed)
  FeaturePlot(PSC_PHs, features = c("Antp", "N", "Dl", "fne", "NimB3", "IM18", "Ance", "phm", "Mlc1"), cols = c("grey","red"), pt.size = .5)
  ggsave(paste(c('umap/compare/umap.mitoCut.markers.mindist_0.4.seed', tmpseed, '.pdf'), collapse = ''), units = 'cm', width = 27, height = 22)
  DimPlot(object = PSC_PHs, reduction = "umap", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
  ggsave(paste(c('umap/compare/umap.mitoCut.bySampleSep.mindist_0.4.seed', tmpseed, '.pdf'), collapse = ''), units = 'cm', width = 36, height = 12)
}

PSC_PHs <- RunUMAP(PSC_PHs, reduction = "pca", dims = 1:21, min.dist = 0.3, n.components = 3, seed.use = 704261)
FeaturePlot(PSC_PHs, features = c("Antp", "N", "Dl", "fne", "NimB3", "IM18", "Ance", "phm", "Mlc1"), cols = c("grey","red"))
ggsave('umap/umap.mitoCut.1_2.markers.mindist_0.3.seed704261.pdf', units = 'cm', width = 27, height = 22)
DimPlot(PSC_PHs, reduction = "umap", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
ggsave('umap/umap.mitoCut.1_2.bySampleSep.mindist_0.3.seed704261.pdf', units = 'cm', width = 36, height = 12)

DimPlot(PSC_PHs, reduction = "umap", group.by = "Subclustering", pt.size = .5)
ggsave('umap/umap.mitoCut.1_2.bySubclustering.mindist_0.3.seed704261.pdf', units = 'cm', width = 13, height = 10)
DimPlot(PSC_PHs, reduction = "umap", group.by = "Subclustering", dims = c(1,3), pt.size = .5)
ggsave('umap/umap.mitoCut.1_3.bySubclustering.mindist_0.3.seed704261.pdf', units = 'cm', width = 13, height = 10)
DimPlot(PSC_PHs, reduction = "umap", group.by = "Subclustering", dims = c(2,3),pt.size = .5)
ggsave('umap/umap.mitoCut.2_3.bySubclustering.mindist_0.3.seed704261.pdf', units = 'cm', width = 13, height = 10)


PSC_PHs <- FindNeighbors(PSC_PHs, dims = 1:21)
for (res in c(1:15)){
  res <- res/10
  PSC_PHs <- FindClusters(PSC_PHs, resolution = res)
  
  DimPlot(PSC_PHs, label = T, pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed704228.mitoCut.1_2.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 12, height = 10)
  DimPlot(PSC_PHs, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed704228.mitoCut.1_3.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 12, height = 10)
  DimPlot(PSC_PHs, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed704228.mitoCut.2_3.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 12, height = 10)
  
  DimPlot(PSC_PHs, label = T, pt.size = 0.5, reduction = 'umap')
  ggsave(paste(c('res/compare/umap.seed704261.mitoCut.1_2.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 12, height = 10)
  DimPlot(PSC_PHs, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'umap')
  ggsave(paste(c('res/compare/umap.seed704261.mitoCut.1_3.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 12, height = 10)
  DimPlot(PSC_PHs, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'umap')
  ggsave(paste(c('res/compare/umap.seed704261.mitoCut.2_3.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 12, height = 10)
}

head(PSC_PHs@meta.data)
PSC_PHs <- FindClusters(PSC_PHs, resolution = 0.7)
DimPlot(PSC_PHs, label = T, pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed704228.mitoCut.1_2.res_0.7.pdf', units = 'cm', width = 12, height = 10)
DimPlot(PSC_PHs, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed704228.mitoCut.1_3.res_0.7.pdf', units = 'cm', width = 12, height = 10)
DimPlot(PSC_PHs, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed704228.mitoCut.2_3.res_0.7.pdf', units = 'cm', width = 12, height = 10)


DimPlot(PSC_PHs, label = T, pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed704261.mitoCut.1_2.min_dist_0.3.res_0.7.pdf', units = 'cm', width = 12, height = 10)
DimPlot(PSC_PHs, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed704261.mitoCut.1_3.min_dist_0.3.res_0.7.pdf', units = 'cm', width = 12, height = 10)
DimPlot(PSC_PHs, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed704261.mitoCut.2_3.min_dist_0.3.res_0.7.pdf', units = 'cm', width = 12, height = 10)


### DEG ###

levels(PSC_PHs@meta.data$seurat_clusters)
# 0 1 2 3 4 5
seurat_clusters_anno <- c('PSC_1', 'PH2', 'PH1_1', 'PSC_2', 'PH1_2', 'PSC_3')
names(x = seurat_clusters_anno) <- levels(x = PSC_PHs)
PSC_PHs <- RenameIdents(PSC_PHs, seurat_clusters_anno)
PSC_PHs@active.ident <- factor(PSC_PHs@active.ident,levels = c('PG1_1', 'PH1_2', 'PH2', 'PSC_1', 'PSC_2', 'PSC_3'))
levels(PSC_PHs@active.ident)

PSC_PHs.markers <- FindAllMarkers(object = PSC_PHs, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
PSC_PHs.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(PSC_PHs.markers, 'degs/findAllMarkers.mitoCut.pct25.thre0.25.onlyPos.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
length(rownames(PSC_PHs.markers))
top10 <- PSC_PHs.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
head(top10)

DoHeatmap(object = PSC_PHs, features = c('N', 'Dl', 'shg', top10$gene), angle = 90, size = 3, raster = F, draw.lines = F)
ggsave('degs/degs.mitoCut.DEgenes.top10_ver2.pdf', units = 'cm', width = 30, height = 20)

DimPlot(PSC_PHs, reduction = "umap", group.by = 'Subclustering', pt.size = 1, label = T, label.size = 4)
ggsave('res/umap.seed704261.mitoCut.1_2.min_dist_0.3.Subclustering.pdf', units = 'cm', width = 14, height = 11)
DimPlot(PSC_PHs, reduction = "umap", group.by = 'Subclustering', pt.size = 1, dims = c(1, 3), label = T, label.size = 4)
ggsave('res/umap.seed704261.mitoCut.1_3.min_dist_0.3.Subclustering.pdf', units = 'cm', width = 14, height = 11)
DimPlot(PSC_PHs, reduction = "umap", group.by = 'Subclustering', pt.size = 1, dims = c(2, 3), label = T, label.size = 4)
ggsave('res/umap.seed704261.mitoCut.2_3.min_dist_0.3.Subclustering.pdf', units = 'cm', width = 14, height = 11)

DimPlot(PSC_PHs, reduction = "umap", group.by = 'RNA_snn_res.0.8_anno', pt.size = 1, label = T, label.size = 4)
ggsave('res/umap.seed704261.mitoCut.1_2.min_dist_0.3.RNA_snn_res.0.8_anno.pdf', units = 'cm', width = 14, height = 11)
DimPlot(PSC_PHs, reduction = "umap", group.by = 'RNA_snn_res.0.8_anno', pt.size = 1, dims = c(1, 3), label = T, label.size = 4)
ggsave('res/umap.seed704261.mitoCut.1_3.min_dist_0.3.RNA_snn_res.0.8_anno.pdf', units = 'cm', width = 14, height = 11)
DimPlot(PSC_PHs, reduction = "umap", group.by = 'RNA_snn_res.0.8_anno', pt.size = 1, dims = c(2, 3), label = T, label.size = 4)
ggsave('res/umap.seed704261.mitoCut.2_3.min_dist_0.3.RNA_snn_res.0.8_anno.pdf', units = 'cm', width = 14, height = 11)
save.image('PSC+SCs.Rdata')


#################################
### Filtering Ambiguous cells ###
#################################


### 1. PSC -> PH1, or vice versa ###
rmcells <- c()
for (simpleType in unique(PSC_PHs@meta.data$RNA_snn_res.0.8_anno_simple)){
  tmpSimpleTypeSub <- subset(PSC_PHs@meta.data, RNA_snn_res.0.8_anno_simple == simpleType)
  if (simpleType == 'PH'){
    for (bc in rownames(tmpSimpleTypeSub)){
      if (is.element(tmpSimpleTypeSub[bc, 'PSC_PHs_res_0.7'], c('PSC_1', 'PSC_2', 'PSC_3'))){
        rmcells <- append(values = bc, rmcells)
      }
    } 
  }else { 
    for (bc in rownames(tmpSimpleTypeSub)){
      if (is.element(tmpSimpleTypeSub[bc, 'PSC_PHs_res_0.7'], c('PH1_1', 'PH1_2', 'PH2'))){
        rmcells <- append(values = bc, rmcells)
      }
    }
  }
}

PSC_PHs.rm <- subset(PSC_PHs, cells = setdiff(rownames(PSC_PHs@meta.data), rmcells))
nrow(PSC_PHs.rm@meta.data) #373, originally 386 cells

### mixed in cluster dimension; ambiguous ###
ambigcells <- data.frame(Embeddings(PSC_PHs.rm, reduction = 'umap'))
ambigcells <- rownames(subset(ambigcells, ambigcells$UMAP_1 < -2))
ambigcells <- rownames(subset(PSC_PHs.rm@meta.data[ambigcells,], RNA_snn_res.0.8_anno == 'PSC'))

DimPlot(PSC_PHs.rm, reduction = "umap", pt.size = 1, cells.highlight = ambigcells) + theme(legend.position = 'None')
ggsave('res/__removed__.umap.seed704261.mitoCut.1_2.min_dist_0.3.Ambiguous.pdf', units = 'cm', width = 11.5, height = 11)
DimPlot(PSC_PHs.rm, reduction = "umap", dims = c(1, 3), pt.size = 1, cells.highlight = ambigcells) + theme(legend.position = 'None')
ggsave('res/__removed__.umap.seed704261.mitoCut.1_3.min_dist_0.3.Ambiguous.pdf', units = 'cm', width = 11.5, height = 11)
DimPlot(PSC_PHs.rm, reduction = "umap", dims = c(2, 3), pt.size = 1, cells.highlight = ambigcells) + theme(legend.position = 'None')
ggsave('res/__removed__.umap.seed704261.mitoCut.2_3.min_dist_0.3.Ambiguous.pdf', units = 'cm', width = 11.5, height = 11)

PSC_PHs.rm <- subset(PSC_PHs.rm, cells = setdiff(rownames(PSC_PHs.rm@meta.data), ambigcells))
nrow(PSC_PHs.rm@meta.data) #370, previously 373 cells

removedtotal <- data.frame(barcode = c(rmcells, ambigcells), note = c('doub_ident', 'doub_ident', 'doub_ident', 'doub_ident', 'doub_ident', 'doub_ident',
                                                                      'doub_ident', 'doub_ident', 'doub_ident', 'doub_ident', 'doub_ident', 'doub_ident', 'doub_ident',
                                                                      'ambiguous', 'ambiguous', 'ambiguous'))
write.table(removedtotal, '__removed_cells__.txt', quote = F, row.names = F, col.names = T, sep = '\t')



### new results ###
DimPlot(PSC_PHs.rm, label = T, pt.size = 1, reduction = 'tsne', label.size = 2)
ggsave('res/__removed__.tsne.seed704228.mitoCut.1_2.res_0.7.pdf', units = 'cm', width = 14, height = 11)
DimPlot(PSC_PHs.rm, label = T, dims = c(1,3), pt.size = 1, reduction = 'tsne', label.size = 2)
ggsave('res/__removed__.tsne.seed704228.mitoCut.1_3.res_0.7.pdf', units = 'cm', width = 14, height = 11)
DimPlot(PSC_PHs.rm, label = T, dims = c(2,3), pt.size = 1, reduction = 'tsne', label.size = 2)
ggsave('res/__removed__.tsne.seed704228.mitoCut.2_3.res_0.7.pdf', units = 'cm', width = 14, height = 11)

DimPlot(PSC_PHs.rm, reduction = "umap", do.return = TRUE, pt.size = 1, label = T, label.size = 4)
ggsave('res/__removed__.umap.seed704261.mitoCut.1_2.min_dist_0.3.res_0.7.pdf', units = 'cm', width = 14, height = 11)
DimPlot(PSC_PHs.rm, reduction = "umap", do.return = TRUE, pt.size = 1, dims = c(1, 3), label = T, label.size = 4)
ggsave('res/__removed__.umap.seed704261.mitoCut.1_3.min_dist_0.3.res_0.7.pdf', units = 'cm', width = 14, height = 11)
DimPlot(PSC_PHs.rm, reduction = "umap", do.return = TRUE, pt.size = 1, dims = c(2, 3), label = T, label.size = 4)
ggsave('res/__removed__.umap.seed704261.mitoCut.2_3.min_dist_0.3.res_0.7.pdf', units = 'cm', width = 14, height = 11)

DimPlot(PSC_PHs.rm, reduction = "umap", group.by = 'Subclustering', pt.size = 1, label = T, label.size = 4)
ggsave('res/__removed__.umap.seed704261.mitoCut.1_2.min_dist_0.3.Subclustering.pdf', units = 'cm', width = 14, height = 11)
DimPlot(PSC_PHs.rm, reduction = "umap", group.by = 'Subclustering', pt.size = 1, dims = c(1, 3), label = T, label.size = 4)
ggsave('res/__removed__.umap.seed704261.mitoCut.1_3.min_dist_0.3.Subclustering.pdf', units = 'cm', width = 14, height = 11)
DimPlot(PSC_PHs.rm, reduction = "umap", group.by = 'Subclustering', pt.size = 1, dims = c(2, 3), label = T, label.size = 4)
ggsave('res/__removed__.umap.seed704261.mitoCut.2_3.min_dist_0.3.Subclustering.pdf', units = 'cm', width = 14, height = 11)

DimPlot(PSC_PHs.rm, reduction = "umap", group.by = 'RNA_snn_res.0.8_anno', pt.size = 1, label = T, label.size = 4)
ggsave('res/__removed__.umap.seed704261.mitoCut.1_2.min_dist_0.3.RNA_snn_res.0.8_anno.pdf', units = 'cm', width = 14, height = 11)
DimPlot(PSC_PHs.rm, reduction = "umap", group.by = 'RNA_snn_res.0.8_anno', pt.size = 1, dims = c(1, 3), label = T, label.size = 4)
ggsave('res/__removed__.umap.seed704261.mitoCut.1_3.min_dist_0.3.RNA_snn_res.0.8_anno.pdf', units = 'cm', width = 14, height = 11)
DimPlot(PSC_PHs.rm, reduction = "umap", group.by = 'RNA_snn_res.0.8_anno', pt.size = 1, dims = c(2, 3), label = T, label.size = 4)
ggsave('res/__removed__.umap.seed704261.mitoCut.2_3.min_dist_0.3.RNA_snn_res.0.8_anno.pdf', units = 'cm', width = 14, height = 11)

levels(PSC_PHs.rm@active.ident)
PSC_PHs.rm.markers <- FindAllMarkers(object = PSC_PHs.rm, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
PSC_PHs.rm.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(PSC_PHs.rm.markers, 'degs/__removed__.findAllMarkers.mitoCut.pct25.thre0.25.onlyPos.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
length(rownames(PSC_PHs.rm.markers))
top10 <- PSC_PHs.rm.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

DoHeatmap(object = PSC_PHs.rm, features = c(top10$gene), angle = 90, size = 3, raster = F, draw.lines = F)
ggsave('degs/__removed__.degs.mitoCut.DEgenes.top10.pdf', units = 'cm', width = 30, height = 20)
DoHeatmap(object = PSC_PHs.rm, features = c('N', 'Dl', 'shg', 'Socs36E', top10$gene), angle = 90, size = 3, raster = F, draw.lines = F)
ggsave('degs/__removed__.degs.mitoCut.DEgenes.top10_ver2.pdf', units = 'cm', width = 30, height = 20)


PSC_PHs.rm@meta.data <- PSC_PHs.rm@meta.data[, -c(11:26)]
head(PSC_PHs.rm@meta.data)
writeLable <-as.matrix(PSC_PHs.rm@meta.data)
writeLable <- data.frame(Barcode = rownames(writeLable), writeLable, check.names = F)
write.table(writeLable, 'PSC+SCs.__removed__.txt', sep = '\t', quote = F, row.names = F, col.names = T)
#save.image('PSC+SCs.__removed__.Rdata')

###
# PH1_2 expresses neuronal marker genes
#load('PSC+SCs.__removed__.Rdata')

PSC_PHs.rm@meta.data$PSC_PHs_res_0.7 <- mapvalues(PSC_PHs.rm@meta.data$PSC_PHs_res_0.7, 
                                                  from = levels(PSC_PHs.rm@meta.data$PSC_PHs_res_0.7),
                                                  to = c('PH 1', 'Neurons', 'PH 2', 'PSC', 'PSC', 'PSC'))
Idents(PSC_PHs.rm) <- 'PSC_PHs_res_0.7'
VlnPlot(PSC_PHs.rm, features = c('nSyb', 'brp', 'Syt1', 'Syt4'), assay = 'RNA', pt.size = 0.5, ncol = 2)
#ggsave('degs/__removed__.markers.neurons_20190913.pdf', units = 'cm', width = 6, height = 10)


