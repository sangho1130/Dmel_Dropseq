### Subclusteirng PM ###
library(Seurat)
library(ggplot2)
library(Matrix)
library(cowplot)
library(plyr)

lymphgland.combined <- readRDS('lymphgland.combined.Rds')

lymphgland.combined.PM <- subset(lymphgland.combined, idents = 'PM')
lymphgland.combined.PM <- ScaleData(object = lymphgland.combined.PM, vars.to.regress = c('Library', 'nCount_RNA'))
lymphgland.combined.PM <- RunPCA(object = lymphgland.combined.PM, npcs = 60)
lymphgland.combined.PM <- JackStraw(object = lymphgland.combined.PM, num.replicate = 100, dims = 60)
lymphgland.combined.PM <- ScoreJackStraw(object = lymphgland.combined.PM, dims = 1:60)

dir.create('stats')
JackStrawPlot(lymphgland.combined.PM, dims = 1:60)
ggsave('stats/2-4.pca.JackStrawPlot.mitoCut.pdf', units = 'cm', width = 30, height = 10)
ElbowPlot(lymphgland.combined.PM, ndims = 60)
ggsave('stats/2-4.pca.ElbowPlot.mitoCut.pdf', units = 'cm', width = 15, height = 10)



### t-SNE ###
for (tmpseed in c(619101:619150)){
  lymphgland.combined.PM <- RunTSNE(object = lymphgland.combined.PM, dims = 1:39, reduction.key = 'tSNE', dim.embed = 3, seed.use = tmpseed)
  FeaturePlot(lymphgland.combined.PM, features = c("Tep4", "Ance", "stg", "Hml", "Pxn", "Nplp2", "Ppn", "NimC1", "Ama"), cols = c("grey","red"), reduction = "tsne")
  ggsave(paste(c('tsne/compare/tsne.mitoCut.1_2.markers.seed', tmpseed, '.pdf'), collapse = ''), units = 'cm', width = 27, height = 22)
  DimPlot(object = lymphgland.combined.PM, reduction = "tsne", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
  ggsave(paste(c('tsne/compare/tsne.mitoCut.1_2.bySampleSep.seed', tmpseed, '.pdf'), collapse = ''), units = 'cm', width = 36, height = 12)
}
lymphgland.combined.PM <- RunTSNE(object = lymphgland.combined.PM, dims = 1:39, reduction.key = 'tSNE', dim.embed = 3, seed.use = 619150)
FeaturePlot(lymphgland.combined.PM, features = c("Tep4", "Ance", "stg", "Hml", "Pxn", "Nplp2", "Ppn", "NimC1", "Ama"), cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.mitoCut.1_2.markers.seed619150.pdf', units = 'cm', width = 27, height = 22)
DimPlot(lymphgland.combined.PM, reduction = "tsne", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
ggsave('tsne/tsne.mitoCut.1_2.bySampleSep.seed619150.pdf', units = 'cm', width = 36, height = 12)



### UMAP ###
for (tmpseed in c(619151:619170)){
  lymphgland.combined.PM <- RunUMAP(lymphgland.combined.PM, reduction = "pca", dims = 1:39, min.dist = 0.3, n.components = 3, seed.use = tmpseed)
  FeaturePlot(lymphgland.combined.PM, features = c("Tep4", "Ance", "stg", "Hml", "Pxn", "Nplp2", "Ppn", "NimC1", "Ama"), cols = c("grey","red"))
  ggsave(paste(c('umap/compare/umap.mitoCut.markers.mindist_0.3.seed', tmpseed, '.pdf'), collapse = ''), units = 'cm', width = 27, height = 22)
  DimPlot(object = lymphgland.combined.PM, reduction = "umap", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
  ggsave(paste(c('umap/compare/umap.mitoCut.bySampleSep.seed', tmpseed, '.mindist_0.3_.pdf'), collapse = ''), units = 'cm', width = 36, height = 12)
}

lymphgland.combined.PM <- RunUMAP(lymphgland.combined.PM, reduction = "pca", dims = 1:39, min.dist = 0.3, n.components = 3, seed.use = 619154)
FeaturePlot(lymphgland.combined.PM, features = c("Tep4", "Ance", "stg", "Hml", "Pxn", "Nplp2", "Ppn", "NimC1", "Ama"), cols = c("grey","red"))
ggsave('umap/umap.mitoCut.markers.mindist_0.3.seed619154.pdf', units = 'cm', width = 27, height = 22)
DimPlot(lymphgland.combined.PM, reduction = "umap", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
ggsave('umap/umap.mitoCut.bySampleSep.mindist_0.3.seed619154.pdf', units = 'cm', width = 36, height = 12)
save.image('tmp1_subclustering_PM.Rdata')



### Clustering ###
lymphgland.combined.PM <- FindNeighbors(lymphgland.combined.PM, dims = 1:42)
for (res in c(3:20)){
  res <- res/10
  lymphgland.combined.PM <- FindClusters(lymphgland.combined.PM, resolution = res)
  DimPlot(lymphgland.combined.PM, label = T, pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed619150.mitoCut.1_2.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 12, height = 10)
  DimPlot(lymphgland.combined.PM, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed619150.mitoCut.1_3.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 12, height = 10)
  DimPlot(lymphgland.combined.PM, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed619150.mitoCut.2_3.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 12, height = 10)
}

lymphgland.combined.PM <- FindClusters(lymphgland.combined.PM, resolution = 0.9)
DimPlot(lymphgland.combined.PM, label = T, pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed619150.mitoCut.1_2.res_0.9.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(lymphgland.combined.PM, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed619150.mitoCut.1_3.res_0.9.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(lymphgland.combined.PM, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed619150.mitoCut.2_3.res_0.9.pdf', units = 'cm', width = 12.5, height = 10)

DimPlot(lymphgland.combined.PM, label = T, pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed619154.mitoCut.1_2.min_dist_0.3.res_0.9.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(lymphgland.combined.PM, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed619154.mitoCut.1_3.min_dist_0.3.res_0.9.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(lymphgland.combined.PM, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed619154.mitoCut.2_3.min_dist_0.3.res_0.9.pdf', units = 'cm', width = 12.5, height = 10)
save.image('tmp2_subclustering_PM.Rdata')


levels(lymphgland.combined.PM@meta.data$seurat_clusters)
seurat_clusters_anno <- c('2_PMsub-Nplp2', '5_PMsub-dUTPase', '1_PMsub-CecA2', '6_PMsub-NimC1', '9_PMsub-Cpr49Ac', '7_PMsub-Cyt-c-p', '10_PMsub-roX1', 
                          '8_PMsub-lectin-24A', '12_PMsub-betaTub60D', '11_PMsub-fax', '3_PMsub-sle', '4_PMsub-stg', '14_PMsub-phm', '13_PMsub-PPO1')

names(x = seurat_clusters_anno) <- levels(x = lymphgland.combined.PM)
lymphgland.combined.PM <- RenameIdents(lymphgland.combined.PM, seurat_clusters_anno)
lymphgland.combined.PM@active.ident <- factor(lymphgland.combined.PM@active.ident, 
                                              levels = c('1_PMsub-CecA2', '2_PMsub-Nplp2', '3_PMsub-sle', '4_PMsub-stg', '5_PMsub-dUTPase', '6_PMsub-NimC1', '7_PMsub-Cyt-c-p', 
                                                         '8_PMsub-lectin-24A', '9_PMsub-Cpr49Ac', '10_PMsub-roX1', '11_PMsub-fax', '12_PMsub-betaTub60D', '13_PMsub-PPO1', '14_PMsub-phm'))
levels(lymphgland.combined.PM@active.ident)
lymphgland.combined.PM@meta.data$RNA_snn_res.0.9_anno <- lymphgland.combined.PM@active.ident
saveRDS(lymphgland.combined.PM, 'lymphgland.combined.PM.Rds')



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

for (subcluster in sort(unique(lymphgland.combined.PM@meta.data$RNA_snn_res.0.9_anno))){
  subcluster_data <- subset(lymphgland.combined.PM@meta.data[,c('Library', 'RNA_snn_res.0.9_anno')], RNA_snn_res.0.9_anno == subcluster)
  numTable <- data.frame(table(subcluster_data$Library))
  freqTable <- freq(subcluster_data$Library)$x[c(1:nrow(numTable)),]
  freqTable <- data.frame(row.names = freqTable$x, perc = freqTable$Percent)
  numTable <- cbind(numTable, percent = freqTable[as.character(numTable$Var1), ])
  write.table(numTable, paste(c('frequencies/number_percent_', subcluster, '.txt'), collapse = ''), quote = F, col.names = T, row.names = F, sep = '\t')
  
  Piechart(numTable$Var1, numTable$Freq, subcluster, paste(c('frequencies/number_percent_', subcluster, '.pdf'), collapse = ''))
}
# 13_PMsub-PPO1, 10_PMsub-roX1, 9_PMsub-Cpr49Ac mainly originated from a single library


### Correlation ###
dir.create('correlation')
fltcells <- rownames(subset(lymphgland.combined.PM@meta.data, RNA_snn_res.0.9_anno != '13_PMsub-PPO1' & RNA_snn_res.0.9_anno != '10_PMsub-roX1' & RNA_snn_res.0.9_anno != '9_PMsub-Cpr49Ac'))
lymphgland.combined.PM_flt <- subset(lymphgland.combined.PM, cells = fltcells)
lymphgland.combined.PM_flt@meta.data <- droplevels(lymphgland.combined.PM_flt@meta.data)

DefaultAssay(lymphgland.combined.PM_flt) <- 'RNA'
pseudo_exprs <- data.frame(row.names = rownames(lymphgland.combined.PM_flt))
for (subcluster in unique(lymphgland.combined.PM_flt@meta.data$RNA_snn_res.0.9_anno)){
  tmp_obj <- subset(lymphgland.combined.PG, idents = subcluster)
  tmp_expr <- as.matrix(GetAssayData(tmp_obj, slot = 'data', assay = 'RNA'))
  tmp_expr[1:6, 1:6]
  psuedo_expr <- data.frame(row.names = rownames(tmp_expr), sum = rowSums(tmp_expr))
  if (identical(rownames(psuedo_expr), rownames(pseudo_exprs))){
    pseudo_exprs <- cbind(pseudo_exprs, psuedo_expr)
  }
}
colnames(pseudo_exprs) <- unique(lymphgland.combined.PM_flt@meta.data$RNA_snn_res.0.9_anno)

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
                   main = 'PM subclustering', filename = 'correlation/SpearmanCorr.PMsub.allGenes.pdf')

pheatmap::pheatmap(corMat, color = colorRampPalette(c("grey60", 'grey60', 'grey60', 'grey60', "grey60", 'grey60', 'grey60',
                                                      'grey60', 'grey60', 'grey60', 'grey60', 'grey60', 'grey60', "red2"))(n = 14), 
                   breaks = c(0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0),
                   fontsize_row = 7, fontsize_col = 9, border_color = NA, 
                   cellwidth = 10, cellheight = 10,
                   clustering_distance_rows = 'euclidean', clustering_method = 'complete', clustering_callback = callback, 
                   main = 'PM subclustering', filename = 'correlation/SpearmanCorr.PMsub.allGenes.contrast.pdf')

corMat <- data.frame(Subcluster = rownames(corMat), corMat, check.names = F)
write.table(corMat, 'correlation/SpearmanCorr.PMsub.allGenes.txt', quote = F, sep = '\t', row.names = F, col.names = T)
saveRDS(lymphgland.combined.PM_flt, 'lymphgland.combined.PM_flt.Rds')
