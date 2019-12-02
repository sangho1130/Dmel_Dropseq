### Subclusteirng CC ###
library(Seurat)
library(ggplot2)
library(Matrix)
library(cowplot)
library(plyr)

lymphgland.combined <- readRDS('lymphgland.combined.Rds')

lymphgland.combined.CC <- subset(lymphgland.combined, idents = 'CC')
lymphgland.combined.CC <- ScaleData(object = lymphgland.combined.CC, vars.to.regress = c('Library', 'nCount_RNA'))
lymphgland.combined.CC <- RunPCA(object = lymphgland.combined.CC, npcs = 60)
lymphgland.combined.CC <- JackStraw(object = lymphgland.combined.CC, num.replicate = 100, dims = 60)
lymphgland.combined.CC <- ScoreJackStraw(object = lymphgland.combined.CC, dims = 1:60)

dir.create('stats')
JackStrawPlot(lymphgland.combined.CC, dims = 1:60)
ggsave('stats/2-4.pca.JackStrawPlot.mitoCut.pdf', units = 'cm', width = 30, height = 10)
ElbowPlot(lymphgland.combined.CC, ndims = 60)
ggsave('stats/2-4.pca.ElbowPlot.mitoCut.pdf', units = 'cm', width = 15, height = 10)



### t-SNE ###
for (tmpseed in c(701101:701120)){
  lymphgland.combined.CC <- RunTSNE(object = lymphgland.combined.CC, dims = 1:39, reduction.key = 'tSNE', dim.embed = 3, seed.use = tmpseed)
  FeaturePlot(lymphgland.combined.CC, features = c('mthl7', "Tep4", "Ance", "Hml", "Pxn", "atilla", "Ama", "PPO2", "PPO1"), cols = c("grey","red"), reduction = "tsne")
  ggsave(paste(c('tsne/compare/tsne.mitoCut.1_2.markers.seed', tmpseed, '.pdf'), collapse = ''), units = 'cm', width = 27, height = 22)
  DimPlot(object = lymphgland.combined.CC, reduction = "tsne", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
  ggsave(paste(c('tsne/compare/tsne.mitoCut.1_2.bySampleSep.seed', tmpseed, '.pdf'), collapse = ''), units = 'cm', width = 36, height = 12)
}
lymphgland.combined.CC <- RunTSNE(object = lymphgland.combined.CC, dims = 1:39, reduction.key = 'tSNE', dim.embed = 3, seed.use = 701106)
FeaturePlot(lymphgland.combined.CC, features = c('mthl7', "Tep4", "Ance", "Hml", "Pxn", "atilla", "Ama", "PPO2", "PPO1"), cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.mitoCut.1_2.markers.seed701106.pdf', units = 'cm', width = 27, height = 22)
DimPlot(lymphgland.combined.CC, reduction = "tsne", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
ggsave('tsne/tsne.mitoCut.1_2.bySampleSep.seed701106.pdf', units = 'cm', width = 36, height = 12)



### UMAP ###
for (tmpseed in c(701151:701170)){
  lymphgland.combined.CC <- RunUMAP(lymphgland.combined.CC, reduction = "pca", dims = 1:39, min.dist = 0.3, n.components = 3, seed.use = tmpseed)
  FeaturePlot(lymphgland.combined.CC, features = c('mthl7', "Tep4", "Ance", "Hml", "Pxn", "atilla", "Ama", "PPO2", "PPO1"), cols = c("grey","red"))
  ggsave(paste(c('umap/compare/umap.mitoCut.markers.mindist_0.3.seed', tmpseed, '.pdf'), collapse = ''), units = 'cm', width = 27, height = 22)
  DimPlot(object = lymphgland.combined.CC, reduction = "umap", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
  ggsave(paste(c('umap/compare/umap.mitoCut.bySampleSep.seed', tmpseed, '.mindist_0.3_.pdf'), collapse = ''), units = 'cm', width = 36, height = 12)
}

lymphgland.combined.CC <- RunUMAP(lymphgland.combined.CC, reduction = "pca", dims = 1:39, min.dist = 0.3, n.components = 3, seed.use = 701164)
FeaturePlot(lymphgland.combined.CC, features = c('mthl7', "Tep4", "Ance", "Hml", "Pxn", "atilla", "Ama", "PPO2", "PPO1"), cols = c("grey","red"))
ggsave('umap/umap.mitoCut.markers.mindist_0.3.seed701164.pdf', units = 'cm', width = 27, height = 22)
DimPlot(lymphgland.combined.CC, reduction = "umap", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
ggsave('umap/umap.mitoCut.bySampleSep.mindist_0.3.seed701164.pdf', units = 'cm', width = 36, height = 12)
save.image('tmp1_subclustering_CC.Rdata')



### Clustering ###
lymphgland.combined.CC <- FindNeighbors(lymphgland.combined.CC, dims = 1:5)
for (res in c(3:20)){
  res <- res/10
  lymphgland.combined.CC <- FindClusters(lymphgland.combined.CC, resolution = res)
  DimPlot(lymphgland.combined.CC, label = T, pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed701106.mitoCut.1_2.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 12, height = 10)
  DimPlot(lymphgland.combined.CC, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed701106.mitoCut.1_3.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 12, height = 10)
  DimPlot(lymphgland.combined.CC, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed701106.mitoCut.2_3.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 12, height = 10)
}

lymphgland.combined.CC <- FindClusters(lymphgland.combined.CC, resolution = 0.1)
DimPlot(lymphgland.combined.CC, label = T, pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed701106.mitoCut.1_2.res_0.1.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(lymphgland.combined.CC, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed701106.mitoCut.1_3.res_0.1.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(lymphgland.combined.CC, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed701106.mitoCut.2_3.res_0.1.pdf', units = 'cm', width = 12.5, height = 10)

DimPlot(lymphgland.combined.CC, label = T, pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed701164.mitoCut.1_2.min_dist_0.3.res_0.1.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(lymphgland.combined.CC, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed701164.mitoCut.1_3.min_dist_0.3.res_0.1.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(lymphgland.combined.CC, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed701164.mitoCut.2_3.min_dist_0.3.res_0.1.pdf', units = 'cm', width = 12.5, height = 10)
save.image('tmp2_subclustering_CC.Rdata')


levels(lymphgland.combined.CC@meta.data$seurat_clusters)
seurat_clusters_anno <- c('1_CC_PPO1-', '3_CC-PPO1+', '2_CC-PPO1+')

names(x = seurat_clusters_anno) <- levels(x = lymphgland.combined.CC)
lymphgland.combined.CC <- RenameIdents(lymphgland.combined.CC, seurat_clusters_anno)
lymphgland.combined.CC@active.ident <- factor(lymphgland.combined.CC@active.ident, levels = c('1_CC_PPO1-', '2_CC-PPO1+', '3_CC-PPO1+'))
levels(lymphgland.combined.CC@active.ident)
lymphgland.combined.CC@meta.data$RNA_snn_res.0.1_anno <- lymphgland.combined.CC@active.ident
saveRDS(lymphgland.combined.CC, 'lymphgland.combined.CC.Rds')



### Frequencies ###
Piechart <- function(columnNames, Values, PlotTitle, outputPdf){
  library(ggplot2)
  library(scales)
  library(RColorBrewer)
  
  data <- data.frame( group = columnNames, value = Values )
  data$group <- factor(data$group, columnNames)
  colourCount <- length(Values)
  getPalette <- colorRampPalette(brewer.pal(9, "Spectral"))
  pie <- ggplot(data, aes(x="", y=value, fill=factor(group))) +
    geom_bar(width = 1, stat = "identity") + 
    coord_polar("y", start=0) + 
    scale_fill_manual(values = getPalette(colourCount)) + 
    labs(title = PlotTitle) + 
    coord_polar(theta = "y", direction = -1) +
    theme_void() + 
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank()) 
  
  ggsave(outputPdf, units = 'cm', height = 10, width = 16)
}


library(RColorBrewer)
library(frequency)
getPalette <- colorRampPalette(brewer.pal(9, "Spectral"))
dir.create('frequencies')

for (subcluster in sort(unique(lymphgland.combined.CC@meta.data$RNA_snn_res.0.1_anno))){
  subcluster_data <- subset(lymphgland.combined.CC@meta.data[,c('Library', 'RNA_snn_res.0.1_anno')], RNA_snn_res.0.1_anno == subcluster)
  numTable <- data.frame(table(subcluster_data$Library))
  freqTable <- freq(subcluster_data$Library)$x[c(1:nrow(numTable)),]
  freqTable <- data.frame(row.names = freqTable$x, perc = freqTable$Percent)
  numTable <- cbind(numTable, percent = freqTable[as.character(numTable$Var1), ])
  write.table(numTable, paste(c('frequencies/number_percent_', subcluster, '.txt'), collapse = ''), quote = F, col.names = T, row.names = F, sep = '\t')
  
  Piechart(numTable$Var1, numTable$Freq, subcluster, paste(c('frequencies/number_percent_', subcluster, '.pdf'), collapse = ''))
}
# 2_CC-PPO1+ mainly originated from a single library



### Correlation ###
dir.create('correlation')
fltcells <- rownames(subset(lymphgland.combined.CC@meta.data, RNA_snn_res.0.1_anno != '2_CC-PPO1+'))
lymphgland.combined.CC_flt <- subset(lymphgland.combined.CC, cells = fltcells)
lymphgland.combined.CC_flt@meta.data <- droplevels(lymphgland.combined.CC_flt@meta.data)

DefaultAssay(lymphgland.combined.CC_flt) <- 'RNA'
pseudo_exprs <- data.frame(row.names = rownames(lymphgland.combined.CC_flt))
for (subcluster in unique(lymphgland.combined.CC_flt@meta.data$RNA_snn_res.0.1_anno)){
  tmp_obj <- subset(lymphgland.combined.PG, idents = subcluster)
  tmp_expr <- as.matrix(GetAssayData(tmp_obj, slot = 'data', assay = 'RNA'))
  psuedo_expr <- data.frame(row.names = rownames(tmp_expr), sum = rowSums(tmp_expr))
  if (identical(rownames(psuedo_expr), rownames(pseudo_exprs))){
    pseudo_exprs <- cbind(pseudo_exprs, psuedo_expr)
  }
}
colnames(pseudo_exprs) <- unique(lymphgland.combined.CC_flt@meta.data$RNA_snn_res.0.1_anno)

corMat <- cor(pseudo_exprs, method = 'spearman')

callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
pheatmap::pheatmap(corMat, color = colorRampPalette(c("steelblue", 'grey60', "red2"))(n = 10),
                   breaks = c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0),
                   fontsize_row = 7, fontsize_col = 9, border_color = NA, 
                   cellwidth = 10, cellheight = 10,
                   clustering_distance_rows = 'euclidean', clustering_method = 'complete', clustering_callback = callback, 
                   main = 'CC subclustering', filename = 'correlation/SpearmanCorr.CCsub.allGenes.pdf')

pheatmap::pheatmap(corMat, color = colorRampPalette(c("grey60", 'grey60', 'grey60', 'grey60', 'grey60', 'grey60', 'grey60', 'grey60', 'grey60', "red2"))(n = 10), 
                   breaks = c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0),
                   fontsize_row = 7, fontsize_col = 9, border_color = NA, 
                   cellwidth = 10, cellheight = 10,
                   clustering_distance_rows = 'euclidean', clustering_method = 'complete', clustering_callback = callback, 
                   main = 'CC subclustering', filename = 'correlation/SpearmanCorr.CCsub.allGenes.contrast.pdf')

corMat <- data.frame(Subcluster = rownames(corMat), corMat, check.names = F)
write.table(corMat, 'correlation/SpearmanCorr.CCsub.allGenes.txt', quote = F, sep = '\t', row.names = F, col.names = T)
saveRDS(lymphgland.combined.CC_flt, 'lymphgland.combined.CC_flt.Rds')
