#################
### Functions ###
#################

# 2. tSNE & UMAP
calc_lat <- function(inputObj, outputPath, pcNum, inputseed, iter = 20, markergenes) {
  dir.create(paste0(c(outputPath, 'tsne'), collapse = '/'))
  dir.create(paste0(c(outputPath, 'tsne/compare'), collapse = '/'))
  dir.create(paste0(c(outputPath, 'umap'), collapse = '/'))
  dir.create(paste0(c(outputPath, 'umap/compare'), collapse = '/'))
  
  seedstart <- inputseed
  seedend <- inputseed + iter
  for (tmpseed in c(seedstart:seedend)){
    inputObj <- RunTSNE(object = inputObj, dims = 1:pcNum, reduction.key = 'tSNE', dim.embed = 3, seed.use = tmpseed)
    FeaturePlot(inputObj, features = markergenes, cols = c("grey","red"), reduction = "tsne")
    ggsave(paste(c(outputPath, '/tsne/compare/tsne.mitoCut.1_2.markers.seed', tmpseed, '.png'), collapse = ''), units = 'cm', width = 27, height = 22)
  }
  seedstart <- inputseed +50
  seedend <- inputseed +50 + iter
  for (tmpseed in c(seedstart:seedend)){
    inputObj <- RunUMAP(inputObj, reduction = "pca", dims = 1:pcNum, min.dist = 0.2, n.components = 3, seed.use = tmpseed)
    FeaturePlot(inputObj, features = markergenes, cols = c("grey","red"))
    ggsave(paste(c(outputPath, '/umap/compare/umap.mitoCut.markers.mindist_0.3.seed', tmpseed, '.png'), collapse = ''), units = 'cm', width = 27, height = 22)
  }
  return(inputObj)
}

######################################################################################################################################
library(Seurat)
library(ggplot2)
library(patchwork)
library(plyr)
library(dplyr)

### 2. Get seeds ###
lymphgland <- readRDS('tmp/lymphgland.Rds')
mygene <- c("Antp", "Dl", "NimB3", "Ance", "Hml", "NimC1", "PPO1", "mthl4", "dys")
lymphgland <- calc_lat(lymphgland, './', 32, 200514200, iter=20, markergenes=mygene) # tsne 200514204, umap 200514250
saveRDS(lymphgland, 'tmp/lymphgland.Rds')
