library(Seurat)
library(ggplot2)
library(patchwork)
library(plyr)
library(dplyr)

setwd('/home/egpark/2020_new/projects/tumor_scRNA-seq/analysis_hg19/seurat3/subclustering/T\ cell/')

### 3. Fix seed values ###
lymphgland <- readRDS('tmp/lymphgland.Rds')
mygene <- c("Antp", "Dl", "NimB3", "Ance", "Hml", "NimC1", "PPO1", "mthl4", "dys")
# tsne fix
lymphgland <- RunTSNE(lymphgland, dims = 1:29, reduction.key = 'tSNE', dim.embed = 3, seed.use = 200513122)
FeaturePlot(lymphgland, features = mygene, cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.mitoCut.1_2.markers.seed200513122.png', units = 'cm', width = 27, height = 22)
# umap fix
lymphgland <- RunUMAP(lymphgland, reduction = "pca", dims = 1:29, min.dist = 0.3, n.components = 3, seed.use = 200513198)#, umap.method='umap-learn', metric='correlation')
FeaturePlot(lymphgland, features = mygene, cols = c("grey","red"))
ggsave('umap/umap.mitoCut.markers.mindist_0.3.seed200513198.png', units = 'cm', width = 27, height = 22)
saveRDS(lymphgland, 'tmp/lymphgland.Rds')
