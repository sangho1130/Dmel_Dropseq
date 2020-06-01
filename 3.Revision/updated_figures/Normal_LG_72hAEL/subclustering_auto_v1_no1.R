#################
### Functions ###
#################
# 1. PCA
calc_pcs <- function(inputObj, outputPath) {
  inputObj <- ScaleData(object = inputObj, vars.to.regress = c('Library', 'nCount_RNA'))
  inputObj <- RunPCA(object = inputObj, npcs = 80)
  inputObj <- JackStraw(object = inputObj, num.replicate = 100, dims = 80)
  inputObj <- ScoreJackStraw(object = inputObj, dims = 1:80)
  
  dir.create(paste0(c(outputPath, 'stats'), collapse = '/'))
  JackStrawPlot(inputObj, dims = 1:80) #
  ggsave(paste0(c(outputPath, 'stats/2-4.pca.JackStrawPlot.mitoCut.pdf'), collapse = '/'), units = 'cm', width = 30, height = 12)
  ElbowPlot(inputObj, ndims = 80) #
  ggsave(paste0(c(outputPath, 'stats/2-4.pca.ElbowPlot.mitoCut.pdf'), collapse = '/'), units = 'cm', width = 15, height = 10)
  return(inputObj)
}

######################################################################################################################################
library(Seurat)
library(ggplot2)
library(patchwork)
library(plyr)
library(dplyr)


lymphgland <- readRDS('../../../../Project1_LymphGland/Drop-seq_alignment/1.seurat3_alignment_withMuscle_regress-Library-nUMI/lymphgland.combined.flt_allscaled.Rds')
Idents(lymphgland) <- 'timepoint'
lymphgland <- subset(lymphgland, idents = 'AEL72hr')
Idents(lymphgland) <- 'Subclustering'
lymphgland <- subset(lymphgland, idents = setdiff(levels(Idents(lymphgland)), c("DV", "RG", "Neurons")))
lymphgland@meta.data <- droplevels(lymphgland@meta.data)
head(lymphgland@meta.data)
levels(lymphgland@meta.data$anno_simple)
levels(lymphgland@meta.data$Subclustering)

dir.create('tmp')
dir.create('tsne')
dir.create('tsne/compare')
dir.create('umap')
dir.create('umap/compare')

### Filter genes ###
filtergenes <- readRDS('../../PH/total ex_scbkgenes/tmp/filtergenes.Rds')
lymphgland <- NormalizeData(object = lymphgland, normalization.method = "LogNormalize", scale.factor = 10000)
lymphgland <- FindVariableFeatures(lymphgland, selection.method = "vst", nfeatures = 2000)
#write.table(intersect(VariableFeatures(lymphgland), filtergenes), 'filteredgenes_72hAEL.txt', quote = F, sep = '\n', row.names = F, col.names = F)

DefaultAssay(lymphgland)
head(lymphgland@assays$RNA@meta.features)
lymphgland@assays$RNA@meta.features$vst.variable <- FALSE
testit <- lymphgland@assays$RNA@meta.features[order(lymphgland@assays$RNA@meta.features$vst.variance.standardized, decreasing = T),]
newvargenes <- setdiff(rownames(testit), filtergenes)[1:2000]
lymphgland@assays$RNA@var.features <- newvargenes
lymphgland@assays$RNA@meta.features[newvargenes, 'vst.variable'] <- TRUE

top10 <- head(VariableFeatures(lymphgland), 10)
plot1 <- VariableFeaturePlot(lymphgland)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave('stats/2-3.variableFeatures.pdf', units = 'cm', width = 17, height = 10)

### 1. Determine PCs ###
lymphgland <- calc_pcs(lymphgland, './') # 29 PCs
saveRDS(lymphgland, 'tmp/lymphgland.Rds')
