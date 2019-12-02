### Subclusteirng PH ###
library(Seurat)
library(ggplot2)
library(Matrix)
library(cowplot)
library(plyr)

lymphgland.combined <- readRDS('lymphgland.combined.Rds')

lymphgland.combined.adipo <- subset(lymphgland.combined, idents = 'Adipohemocyte')
lymphgland.combined.adipo <- ScaleData(object = lymphgland.combined.adipo, vars.to.regress = c('Library', 'nCount_RNA'))
lymphgland.combined.adipo <- RunPCA(object = lymphgland.combined.adipo, npcs = 60)
lymphgland.combined.adipo <- JackStraw(object = lymphgland.combined.adipo, num.replicate = 100, dims = 60)
lymphgland.combined.adipo <- ScoreJackStraw(object = lymphgland.combined.adipo, dims = 1:60)

dir.create('stats')
JackStrawPlot(lymphgland.combined.adipo, dims = 1:60)
ggsave('stats/2-4.pca.JackStrawPlot.mitoCut.pdf', units = 'cm', width = 30, height = 10)
ElbowPlot(lymphgland.combined.adipo, ndims = 60)
ggsave('stats/2-4.pca.ElbowPlot.mitoCut.pdf', units = 'cm', width = 15, height = 10)



### t-SNE ###
for (tmpseed in c(620101:620110)){
  lymphgland.combined.adipo <- RunTSNE(object = lymphgland.combined.adipo, dims = 1:39, reduction.key = 'tSNE', dim.embed = 3, seed.use = tmpseed)
  FeaturePlot(lymphgland.combined.adipo, features = c("Tep4", "Ance", "Hml", "Pxn", "atilla", "mthl4", "Ama", "PPO2", "hid"), cols = c("grey","red"), reduction = "tsne")
  ggsave(paste(c('tsne/compare/tsne.mitoCut.1_2.markers.seed', tmpseed, '.pdf'), collapse = ''), units = 'cm', width = 27, height = 22)
  DimPlot(object = lymphgland.combined.adipo, reduction = "tsne", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
  ggsave(paste(c('tsne/compare/tsne.mitoCut.1_2.bySampleSep.seed', tmpseed, '.pdf'), collapse = ''), units = 'cm', width = 36, height = 12)
}
lymphgland.combined.adipo <- RunTSNE(object = lymphgland.combined.adipo, dims = 1:39, reduction.key = 'tSNE', dim.embed = 3, seed.use = 620110)
FeaturePlot(lymphgland.combined.adipo, features = c("Tep4", "Ance", "Hml", "Pxn", "atilla", "mthl4", "Ama", "PPO2", "hid"), cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.mitoCut.1_2.markers.seed620110.pdf', units = 'cm', width = 27, height = 22)
DimPlot(lymphgland.combined.adipo, reduction = "tsne", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
ggsave('tsne/tsne.mitoCut.1_2.bySampleSep.seed620110.pdf', units = 'cm', width = 36, height = 12)



### UMAP ###
for (tmpseed in c(620151:620160)){
  lymphgland.combined.adipo <- RunUMAP(lymphgland.combined.adipo, reduction = "pca", dims = 1:39, min.dist = 0.3, n.components = 3, seed.use = tmpseed)
  FeaturePlot(lymphgland.combined.adipo, features = c("Tep4", "Ance", "Hml", "Pxn", "atilla", "mthl4", "Ama", "PPO2", "hid"), cols = c("grey","red"))
  ggsave(paste(c('umap/compare/umap.mitoCut.markers.mindist_0.3.seed', tmpseed, '.pdf'), collapse = ''), units = 'cm', width = 27, height = 22)
  DimPlot(object = lymphgland.combined.adipo, reduction = "umap", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
  ggsave(paste(c('umap/compare/umap.mitoCut.bySampleSep.seed', tmpseed, '.mindist_0.3_.pdf'), collapse = ''), units = 'cm', width = 36, height = 12)
}

lymphgland.combined.adipo <- RunUMAP(lymphgland.combined.adipo, reduction = "pca", dims = 1:39, min.dist = 0.3, n.components = 3, seed.use = 620158)
FeaturePlot(lymphgland.combined.adipo, features = c("Tep4", "Ance", "Hml", "Pxn", "atilla", "mthl4", "Ama", "PPO2", "hid"), cols = c("grey","red"))
ggsave('umap/umap.mitoCut.markers.mindist_0.3.seed620158.pdf', units = 'cm', width = 27, height = 22)
DimPlot(lymphgland.combined.adipo, reduction = "umap", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
ggsave('umap/umap.mitoCut.bySampleSep.mindist_0.3.seed620158.pdf', units = 'cm', width = 36, height = 12)
save.image('tmp1_subclustering_adipo.Rdata')



### Clustering ###
lymphgland.combined.adipo <- FindNeighbors(lymphgland.combined.adipo, dims = 1:16)
for (res in c(3:20)){
  res <- res/10
  lymphgland.combined.adipo <- FindClusters(lymphgland.combined.adipo, resolution = res)
  DimPlot(lymphgland.combined.adipo, label = T, pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed620110.mitoCut.1_2.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 12, height = 10)
  DimPlot(lymphgland.combined.adipo, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed620110.mitoCut.1_3.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 12, height = 10)
  DimPlot(lymphgland.combined.adipo, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'tsne')
  ggsave(paste(c('res/compare/tsne.seed620110.mitoCut.2_3.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 12, height = 10)
}

lymphgland.combined.adipo <- FindClusters(lymphgland.combined.adipo, resolution = 0.3)
DimPlot(lymphgland.combined.adipo, label = T, pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed620110.mitoCut.1_2.res_0.3.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(lymphgland.combined.adipo, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed620110.mitoCut.1_3.res_0.3.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(lymphgland.combined.adipo, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'tsne')
ggsave('res/tsne.seed620110.mitoCut.2_3.res_0.3.pdf', units = 'cm', width = 12.5, height = 10)

DimPlot(lymphgland.combined.adipo, label = T, pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed620158.mitoCut.1_2.min_dist_0.3.res_0.3.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(lymphgland.combined.adipo, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed620158.mitoCut.1_3.min_dist_0.3.res_0.3.pdf', units = 'cm', width = 12.5, height = 10)
DimPlot(lymphgland.combined.adipo, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'umap')
ggsave('res/umap.seed620158.mitoCut.2_3.min_dist_0.3.res_0.3.pdf', units = 'cm', width = 12.5, height = 10)
save.image('tmp2_subclustering_adipo.Rdata')


levels(lymphgland.combined.adipo@meta.data$seurat_clusters)
seurat_clusters_anno <- c('adipo1', 'adipo2', 'adipo3')

names(x = seurat_clusters_anno) <- levels(x = lymphgland.combined.adipo)
lymphgland.combined.adipo <- RenameIdents(lymphgland.combined.adipo, seurat_clusters_anno)
lymphgland.combined.adipo@active.ident <- factor(lymphgland.combined.adipo@active.ident, levels = c('adipo1', 'adipo2', 'adipo3'))
levels(lymphgland.combined.adipo@active.ident)
lymphgland.combined.adipo@meta.data$RNA_snn_res.0.3_anno <- lymphgland.combined.adipo@active.ident
saveRDS(lymphgland.combined.adipo, 'lymphgland.combined.adipo.Rds')



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

for (subcluster in sort(unique(lymphgland.combined.adipo@meta.data$RNA_snn_res.0.3_anno))){
  subcluster_data <- subset(lymphgland.combined.adipo@meta.data[,c('Library', 'RNA_snn_res.0.3_anno')], RNA_snn_res.0.3_anno == subcluster)
  numTable <- data.frame(table(subcluster_data$Library))
  freqTable <- freq(subcluster_data$Library)$x[c(1:nrow(numTable)),]
  freqTable <- data.frame(row.names = freqTable$x, perc = freqTable$Percent)
  numTable <- cbind(numTable, percent = freqTable[as.character(numTable$Var1), ])
  write.table(numTable, paste(c('frequencies/number_percent_', subcluster, '.txt'), collapse = ''), quote = F, col.names = T, row.names = F, sep = '\t')
  
  Piechart(numTable$Var1, numTable$Freq, subcluster, paste(c('frequencies/number_percent_', subcluster, '.pdf'), collapse = ''))
}
# adipo2 & adipo3 mainly originated from a single library


### Correlation ###
dir.create('correlation')
fltcells <- rownames(subset(lymphgland.combined.adipo@meta.data, RNA_snn_res.0.3_anno != '13_PMsub-PPO1' & RNA_snn_res.0.3_anno != '10_PMsub-roX1' & RNA_snn_res.0.3_anno != '9_PMsub-Cpr49Ac'))
lymphgland.combined.adipo_flt <- subset(lymphgland.combined.adipo, cells = fltcells)
lymphgland.combined.adipo_flt@meta.data <- droplevels(lymphgland.combined.adipo_flt@meta.data)

DefaultAssay(lymphgland.combined.adipo_flt) <- 'RNA'
pseudo_exprs <- data.frame(row.names = rownames(lymphgland.combined.adipo_flt))
for (subcluster in unique(lymphgland.combined.adipo_flt@meta.data$RNA_snn_res.0.3_anno)){
  tmp_obj <- subset(lymphgland.combined.PG, idents = subcluster)
  tmp_expr <- as.matrix(GetAssayData(tmp_obj, slot = 'data', assay = 'RNA'))
  tmp_expr[1:6, 1:6]
  psuedo_expr <- data.frame(row.names = rownames(tmp_expr), sum = rowSums(tmp_expr))
  if (identical(rownames(psuedo_expr), rownames(pseudo_exprs))){
    pseudo_exprs <- cbind(pseudo_exprs, psuedo_expr)
  }
}
colnames(pseudo_exprs) <- unique(lymphgland.combined.adipo_flt@meta.data$RNA_snn_res.0.3_anno)

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
                   main = 'Adipo subclustering', filename = 'correlation/SpearmanCorr.adiposub.allGenes.pdf')

pheatmap::pheatmap(corMat, color = colorRampPalette(c("grey60", 'grey60', 'grey60', 'grey60', "grey60", 'grey60', 'grey60',
                                                      'grey60', 'grey60', 'grey60', 'grey60', 'grey60', 'grey60', "red2"))(n = 14), 
                   breaks = c(0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0),
                   fontsize_row = 7, fontsize_col = 9, border_color = NA, 
                   cellwidth = 10, cellheight = 10,
                   clustering_distance_rows = 'euclidean', clustering_method = 'complete', clustering_callback = callback, 
                   main = 'Adipo subclustering', filename = 'correlation/SpearmanCorr.adiposub.allGenes.contrast.pdf')

corMat <- data.frame(Subcluster = rownames(corMat), corMat, check.names = F)
write.table(corMat, 'correlation/SpearmanCorr.adiposub.allGenes.txt', quote = F, sep = '\t', row.names = F, col.names = T)
saveRDS(lymphgland.combined.adipo_flt, 'lymphgland.combined.adipo_flt.Rds')
