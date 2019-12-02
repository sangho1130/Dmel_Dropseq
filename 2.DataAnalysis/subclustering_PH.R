### Subclusteirng PH ###
library(Seurat)
library(ggplot2)
library(Matrix)
library(cowplot)
library(plyr)

lymphgland.combined <- readRDS('lymphgland.combined.Rds')

lymphgland.combined.PH <- subset(lymphgland.combined, idents = 'PH')
lymphgland.combined.PH <- ScaleData(object = lymphgland.combined.PH, vars.to.regress = c('Library', 'nCount_RNA'))
lymphgland.combined.PH <- RunPCA(object = lymphgland.combined.PH, npcs = 60)
lymphgland.combined.PH <- JackStraw(object = lymphgland.combined.PH, num.replicate = 100, dims = 60)
lymphgland.combined.PH <- ScoreJackStraw(object = lymphgland.combined.PH, dims = 1:60)

dir.create('stats')
JackStrawPlot(lymphgland.combined.PH, dims = 1:60)
ggsave('stats/2-4.pca.JackStrawPlot.mitoCut.pdf', units = 'cm', width = 30, height = 10)
ElbowPlot(lymphgland.combined.PH, ndims = 60)
ggsave('stats/2-4.pca.ElbowPlot.mitoCut.pdf', units = 'cm', width = 15, height = 10)



### t-SNE ###
for (tmpseed in c(619201:619250)){
  lymphgland.combined.PH <- RunTSNE(object = lymphgland.combined.PH, dims = 1:39, reduction.key = 'tSNE', dim.embed = 3, seed.use = tmpseed)
  FeaturePlot(lymphgland.combined.PH, features = c("N", "Dl", "rdgA", "IM18", "Tep4", "Ance", "CecA2", "PGRP-SC2", "PCNA"), cols = c("grey","red"), reduction = "tsne")
  ggsave(paste(c('tsne/compare/tsne.mitoCut.1_2.markers.seed', tmpseed, '.pdf'), collapse = ''), units = 'cm', width = 27, height = 22)
  DimPlot(object = lymphgland.combined.PH, reduction = "tsne", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
  ggsave(paste(c('tsne/compare/tsne.mitoCut.1_2.bySampleSep.seed', tmpseed, '.pdf'), collapse = ''), units = 'cm', width = 36, height = 12)
}
lymphgland.combined.PH <- RunTSNE(object = lymphgland.combined.PH, dims = 1:39, reduction.key = 'tSNE', dim.embed = 3, seed.use = 619216)
FeaturePlot(lymphgland.combined.PH, features = c("N", "Dl", "rdgA", "IM18", "Tep4", "Ance", "CecA2", "PGRP-SC2", "PCNA"), cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.mitoCut.1_2.markers.seed619216.pdf', units = 'cm', width = 27, height = 22)
DimPlot(lymphgland.combined.PH, reduction = "tsne", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
ggsave('tsne/tsne.mitoCut.1_2.bySampleSep.seed619216.pdf', units = 'cm', width = 36, height = 12)



### UMAP ###
for (tmpseed in c(619251:619270)){
  lymphgland.combined.PH <- RunUMAP(lymphgland.combined.PH, reduction = "pca", dims = 1:39, min.dist = 0.3, n.components = 3, seed.use = tmpseed)
  FeaturePlot(lymphgland.combined.PH, features = c("N", "Dl", "rdgA", "IM18", "Tep4", "Ance", "CecA2", "PGRP-SC2", "PCNA"), cols = c("grey","red"))
  ggsave(paste(c('umap/compare/umap.mitoCut.markers.mindist_0.3.seed', tmpseed, '.pdf'), collapse = ''), units = 'cm', width = 27, height = 22)
  DimPlot(object = lymphgland.combined.PH, reduction = "umap", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
  ggsave(paste(c('umap/compare/umap.mitoCut.bySampleSep.seed', tmpseed, '.mindist_0.3_.pdf'), collapse = ''), units = 'cm', width = 36, height = 12)
}

lymphgland.combined.PH <- RunUMAP(lymphgland.combined.PH, reduction = "pca", dims = 1:39, min.dist = 0.3, n.components = 3, seed.use = 619252)
FeaturePlot(lymphgland.combined.PH, features = c("N", "Dl", "rdgA", "IM18", "Tep4", "Ance", "CecA2", "PGRP-SC2", "PCNA"), cols = c("grey","red"))
ggsave('umap/umap.mitoCut.markers.mindist_0.3.seed619252.pdf', units = 'cm', width = 27, height = 22)
DimPlot(lymphgland.combined.PH, reduction = "umap", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
ggsave('umap/umap.mitoCut.bySampleSep.mindist_0.3.seed619252.pdf', units = 'cm', width = 36, height = 12)
save.image('tmp1_subclustering_PH.Rdata')



### Clustering ###
lymphgland.combined.PH <- FindNeighbors(lymphgland.combined.PH, dims = 1:39)
for (res in c(3:20)){
  res <- res/10
  lymphgland.combined.PH <- FindClusters(lymphgland.combined.PH, resolution = res)
  DimPlot(lymphgland.combined.PH, label = T, pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed619216.mitoCut.1_2.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 12, height = 10)
  DimPlot(lymphgland.combined.PH, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed619216.mitoCut.1_3.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 12, height = 10)
  DimPlot(lymphgland.combined.PH, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed619216.mitoCut.2_3.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 12, height = 10)
}

lymphgland.combined.PH <- FindClusters(lymphgland.combined.PH, resolution = 0.8)
DimPlot(lymphgland.combined.PH, label = T, pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed619216.mitoCut.1_2.res_0.8.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(lymphgland.combined.PH, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed619216.mitoCut.1_3.res_0.8.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(lymphgland.combined.PH, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed619216.mitoCut.2_3.res_0.8.pdf', units = 'cm', width = 12.5, height = 10)

DimPlot(lymphgland.combined.PH, label = T, pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed619252.mitoCut.1_2.min_dist_0.3.res_0.8.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(lymphgland.combined.PH, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed619252.mitoCut.1_3.min_dist_0.3.res_0.8.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(lymphgland.combined.PH, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed619252.mitoCut.2_3.min_dist_0.3.res_0.8.pdf', units = 'cm', width = 12.5, height = 10)
save.image('tmp2_subclustering_PG.Rdata')


levels(lymphgland.combined.PH@meta.data$seurat_clusters)
seurat_clusters_anno <- c('11_PGsub-CecA2', '8_PGsub-HmgD', '5_PGsub-dUTPase', '9_PGsub-Hml', '10_PGsub-CecA2', '4_PGsub-IM18', '12_PGsub-sea',
                          '13_PGsub-fax', '6_PGsub-Pen', '7_PGsub-GstD10', '1_PGsub-E(spl)m4-BFM', '3_PGsub-CR45018', '2_PGsub-NimB3')

names(x = seurat_clusters_anno) <- levels(x = lymphgland.combined.PH)
lymphgland.combined.PH <- RenameIdents(lymphgland.combined.PH, seurat_clusters_anno)
lymphgland.combined.PH@active.ident <- factor(lymphgland.combined.PH@active.ident, 
                                              levels = c('1_PGsub-E(spl)m4-BFM', '2_PGsub-NimB3', '3_PGsub-CR45018', '4_PGsub-IM18', '5_PGsub-dUTPase', '6_PGsub-Pen',
                                                         '7_PGsub-GstD10', '8_PGsub-HmgD', '9_PGsub-Hml', '10_PGsub-CecA2', '11_PGsub-CecA2', '12_PGsub-sea', '13_PGsub-fax'))
levels(lymphgland.combined.PH@active.ident)
lymphgland.combined.PH@meta.data$RNA_snn_res.0.8_anno <- lymphgland.combined.PH@active.ident
saveRDS(lymphgland.combined.PH, 'lymphgland.combined.PG.Rds')



### Frequencies ###
Piechart <- function(columnNames, Values, PlotTitle, outputPdf){
  library(ggplot2)
  library(scales)
  library(RColorBrewer)
  
  data <- data.frame(
    group = columnNames,
    value = Values
  )
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
    theme(plot.title = element_text(hjust = 0.5),
          legend.title = element_blank()) 
  
  ggsave(outputPdf, units = 'cm', height = 10, width = 16)
}


library(RColorBrewer)
library(frequency)
getPalette <- colorRampPalette(brewer.pal(9, "Spectral"))
dir.create('frequencies')

for (subcluster in sort(unique(lymphgland.combined.PH@meta.data$RNA_snn_res.0.8_anno))){
  subcluster_data <- subset(lymphgland.combined.PH@meta.data[,c('Library', 'RNA_snn_res.0.8_anno')], RNA_snn_res.0.8_anno == subcluster)
  numTable <- data.frame(table(subcluster_data$Library))
  freqTable <- freq(subcluster_data$Library)$x[c(1:nrow(numTable)),]
  freqTable <- data.frame(row.names = freqTable$x, perc = freqTable$Percent)
  numTable <- cbind(numTable, percent = freqTable[as.character(numTable$Var1), ])
  write.table(numTable, paste(c('frequencies/number_percent_', subcluster, '.txt'), collapse = ''), quote = F, col.names = T, row.names = F, sep = '\t')
  
  Piechart(numTable$Var1, numTable$Freq, subcluster, paste(c('frequencies/number_percent_', subcluster, '.pdf'), collapse = ''))
}
# 3_PGsub-CR45018 and 12_PGsub-sea mainly originated from a single library


### Correlation ###
dir.create('correlation')
fltcells <- rownames(subset(lymphgland.combined.PH@meta.data, RNA_snn_res.0.8_anno != '3_PGsub-CR45018' & RNA_snn_res.0.8_anno != '12_PGsub-sea'))
lymphgland.combined.PH_flt <- subset(lymphgland.combined.PH, cells = fltcells)
lymphgland.combined.PH_flt@meta.data <- droplevels(lymphgland.combined.PH_flt@meta.data)

DefaultAssay(lymphgland.combined.PH_flt) <- 'RNA'
pseudo_exprs <- data.frame(row.names = rownames(lymphgland.combined.PH_flt))
for (subcluster in unique(lymphgland.combined.PH_flt@meta.data$RNA_snn_res.0.8_anno)){
  tmp_obj <- subset(lymphgland.combined.PG, idents = subcluster)
  tmp_expr <- as.matrix(GetAssayData(tmp_obj, slot = 'data', assay = 'RNA'))
  tmp_expr[1:6, 1:6]
  psuedo_expr <- data.frame(row.names = rownames(tmp_expr), sum = rowSums(tmp_expr))
  if (identical(rownames(psuedo_expr), rownames(pseudo_exprs))){
    pseudo_exprs <- cbind(pseudo_exprs, psuedo_expr)
  }
}
colnames(pseudo_exprs) <- unique(lymphgland.combined.PH_flt@meta.data$RNA_snn_res.0.8_anno)

corMat <- cor(pseudo_exprs, method = 'spearman')

callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
pheatmap::pheatmap(corMat, color = colorRampPalette(c("steelblue", 'grey60', "red2"))(n = 14),
                   breaks = c(0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0),
                   fontsize_row = 7, fontsize_col = 9, border_color = NA, 
                   cellwidth = 10, cellheight = 10,
                   clustering_distance_rows = 'euclidean', clustering_method = 'complete', clustering_callback = callback, 
                   main = 'PH subclustering', filename = 'correlation/SpearmanCorr.PHsub.allGenes.pdf')

pheatmap::pheatmap(corMat, color = colorRampPalette(c("grey60", 'grey60', 'grey60', 'grey60', "grey60", 'grey60', 'grey60',
                                                      'grey60', 'grey60', 'grey60', 'grey60', 'grey60', 'grey60', "red2"))(n = 14), 
                   breaks = c(0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0),
                   fontsize_row = 7, fontsize_col = 9, border_color = NA, 
                   cellwidth = 10, cellheight = 10,
                   clustering_distance_rows = 'euclidean', clustering_method = 'complete', clustering_callback = callback, 
                   main = 'PH subclustering', filename = 'correlation/SpearmanCorr.PHsub.allGenes.contrast.pdf')

corMat <- data.frame(Subcluster = rownames(corMat), corMat, check.names = F)
write.table(corMat, 'correlation/SpearmanCorr.PHsub.allGenes.txt', quote = F, sep = '\t', row.names = F, col.names = T)
saveRDS(lymphgland.combined.PH_flt, 'lymphgland.combined.PH_flt.Rds')
