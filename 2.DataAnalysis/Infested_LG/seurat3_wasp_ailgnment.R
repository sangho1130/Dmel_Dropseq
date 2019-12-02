#!/home/sangho/miniconda3/envs/seurat3/bin/R
library(Seurat)
library(ggplot2)
library(Matrix)
library(cowplot)
library(plyr)
library(dplyr)

### Normal lymph gland ###
lg96_normal_expr <- read.delim('merged.3tps.expr.allGenes.txt', header = T, sep = '\t', check.names = F, row.names = 1)
lg96_normal_label <- read.delim('__filtered__SCrm__seurat3.3tps.newlabels.txt', header = T, sep = '\t', check.names = F, row.names = 1)
lg96_normal_label <- subset(lg96_normal_label, timepoint == 'AEL96hr')
head(lg96_normal_label); nrow(lg96_normal_label) #9383

lg96_normal_label <- lg96_normal_label[, c(2,3,4,5,8,9,11)]
colnames(lg96_normal_label)[5:7] <- c('cluster_anno', 'simple_anno', 'subcluster_anno')
head(lg96_normal_label)
lg96_normal_expr <- lg96_normal_expr[, rownames(lg96_normal_label)]


### Infestes lymph gland ###
lg96_infested_expr <- read.delim('merged.pi24lg.expr.allGenes.txt', header = T, sep = '\t', check.names = F, row.names = 1)
lg96_infested_label <- read.delim('merged.pi24lg.label.allGenes.txt', header = T, sep = '\t', check.names = F, row.names = 1)
head(lg96_infested_label); nrow(lg96_infested_label) #10179

lg96_infested_label <- lg96_infested_label[, c(2:5, ncol(lg96_infested_label), ncol(lg96_infested_label), ncol(lg96_infested_label))]
colnames(lg96_infested_label)[5:7] <- c('cluster_anno', 'simple_anno', 'subcluster_anno')
head(lg96_infested_label)
lg96_infested_expr <- lg96_infested_expr[, rownames(lg96_infested_label)]

identical(colnames(lg96_normal_label), colnames(lg96_infested_label))
lg96_normal_expr <- as(as.matrix(lg96_normal_expr), "dgCMatrix")
lg96_infested_expr <- as(as.matrix(lg96_infested_expr), "dgCMatrix")

### Seurat objects for each dataset ###
lg96_normal <- CreateSeuratObject(counts = lg96_normal_expr, project = "LG_wasp")
lg96_normal <- AddMetaData(object = lg96_normal, metadata = lg96_normal_label)
lg96_normal <- NormalizeData(lg96_normal, normalization.method = "LogNormalize", scale.factor = 10000)
lg96_normal <- FindVariableFeatures(object = lg96_normal, selection.method = "vst", nfeatures = 2000)
lg96_normal@meta.data$timepoint <- "96AEL_Normal"
head(lg96_normal@meta.data)

lg96_infested <- CreateSeuratObject(counts = lg96_infested_expr, project = "LG_wasp")
lg96_infested <- AddMetaData(object = lg96_infested, metadata = lg96_infested_label)
lg96_infested <- NormalizeData(lg96_infested, normalization.method = "LogNormalize", scale.factor = 10000)
lg96_infested <- FindVariableFeatures(object = lg96_infested, selection.method = "vst", nfeatures = 2000)
lg96_infested@meta.data$timepoint <- "96AEL_Infested"
head(lg96_infested@meta.data)


#load('tmp1.Rdata')
lymphgland.anchors <- FindIntegrationAnchors(object.list = list(lg96_normal, lg96_infested), dims = 1:30)
lymphgland.combined <- IntegrateData(anchorset = lymphgland.anchors, dims = 1:30)
head(lymphgland.combined@meta.data)
DefaultAssay(object = lymphgland.combined) <- "integrated"
lymphgland.combined@meta.data$timepoint <- factor(lymphgland.combined@meta.data$timepoint, c('96AEL_Normal', '96AEL_Infested'))

### The standard workflow ###
lymphgland.combined <- ScaleData(object = lymphgland.combined, vars.to.regress = c('Library', 'nCount_RNA'))
lymphgland.combined <- RunPCA(object = lymphgland.combined, npcs = 80)

lymphgland.combined <- JackStraw(object = lymphgland.combined, num.replicate = 100, dims = 80)
lymphgland.combined <- ScoreJackStraw(object = lymphgland.combined, dims = 1:80)
JackStrawPlot(lymphgland.combined, dims = 1:80) #PC 46
ggsave('stats/2-4.pca.JackStrawPlot.mitoCut.pdf', units = 'cm', width = 30, height = 12)
ElbowPlot(lymphgland.combined, ndims = 80) #PC 46
ggsave('stats/2-4.pca.ElbowPlot.mitoCut.pdf', units = 'cm', width = 15, height = 10)
save.image('tmp2.Rdata')


#load('tmp2.Rdata')
for (tmpseed in c(96711151:96711160)){
  lymphgland.combined <- RunUMAP(lymphgland.combined, reduction = "pca", dims = 1:46, min.dist = 0.4, n.components = 3, seed.use = tmpseed)
  FeaturePlot(lymphgland.combined, features = c("Antp", "Dl", "Ance", "IM18", "Hml", "Ama", "mthl4", "PPO1", "Mlc1"), cols = c("grey","red"))
  ggsave(paste(c('umap/compare/umap.mitoCut.markers.mindist_0.4.seed', tmpseed, '.pdf'), collapse = ''), units = 'cm', width = 27, height = 22)
}

for (tmpseed in c(96710101:96710130)){
  lymphgland.combined <- RunTSNE(object = lymphgland.combined, dims = 1:46, reduction.key = 'tSNE', dim.embed = 3, seed.use = tmpseed)
  FeaturePlot(lymphgland.combined, features = c("Antp", "Dl", "Ance", "IM18", "Hml", "Ama", "mthl4", "PPO1", "Mlc1"), cols = c("grey","red"), reduction = "tsne")
  ggsave(paste(c('tsne/compare/tsne.mitoCut.1_2.markers.seed', tmpseed, '.pdf'), collapse = ''), units = 'cm', width = 27, height = 22)
}
save.image('tmp3.Rdata')

#'cluster_anno', 'simple_anno', 'subcluster_anno'
lymphgland.combined <- RunUMAP(lymphgland.combined, reduction = "pca", dims = 1:46, min.dist = 0.4, n.components = 3, seed.use = 96711161)
FeaturePlot(lymphgland.combined, features = c("Antp", "Dl", "IM18", "Ance", "Hml", "Ama", "mthl4", "PPO1", "Mlc1"), cols = c("grey","red"))
ggsave('umap/umap.mitoCut.markers.mindist_0.4.seed96711161.pdf', units = 'cm', width = 27, height = 22)
DimPlot(object = lymphgland.combined, reduction = "umap", group.by = "timepoint", pt.size = .5)
ggsave('umap/umap.mitoCut.bySample.seed96711161.mindist_0.4_.pdf', units = 'cm', width = 17, height = 12)
DimPlot(object = lymphgland.combined, reduction = "umap", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
ggsave('umap/umap.mitoCut.bySampleSep.seed96711161.mindist_0.4_.pdf', units = 'cm', width = 26, height = 12)


lymphgland.combined <- RunTSNE(object = lymphgland.combined, dims = 1:46, reduction.key = 'tSNE', dim.embed = 3, seed.use = 96710126)
FeaturePlot(lymphgland.combined, features = c("Antp", "Dl", "IM18", "Ance", "Hml", "Ama", "mthl4", "PPO1", "Mlc1"), cols = c("grey","red"), reduction = "tsne")
ggsave('tsne/tsne.mitoCut.1_2.markers.seed96710126.pdf', units = 'cm', width = 27, height = 22)
DimPlot(object = lymphgland.combined, reduction = "tsne", group.by = "timepoint", pt.size = .5)
ggsave('tsne/tsne.mitoCut.1_2.bySample.seed96710126.pdf', units = 'cm', width = 17, height = 12)
DimPlot(object = lymphgland.combined, reduction = "tsne", group.by = "timepoint", pt.size = .5) + facet_wrap(~timepoint)
ggsave('tsne/tsne.mitoCut.1_2.bySampleSep.seed96710126.pdf', units = 'cm', width = 26, height = 12)

save.image('tmp3.Rdata')


### bring subclustering labels ###
newlabel <- read.delim('seurat3.96lgs.mitoCut.label.allGenes.normalized.ver4.txt', header = T, sep = '\t', row.names = 1)
identical(rownames(lymphgland.combined@meta.data), rownames(newlabel))
lymphgland.combined@meta.data <- newlabel
head(lymphgland.combined@meta.data)
lymphgland.combined@meta.data$LG_aligned_res.0.8_anno_v2 <- factor(lymphgland.combined@meta.data$LG_aligned_res.0.8_anno_v2,
                                                                   levels = rev(c('PSC', 'SC-like', 'PH 1', 'PH 2', 'PH 3', 'PH 4',
                                                                                  'PM 1', 'PM 2', 'PM 3', 'PM 4', 'PM 5', 'PM 6',
                                                                                  'LM', 'CC', 'GST-rich', 'DV', 'RG', 'Neurons')))
lymphgland.combined@meta.data$LG_aligned_res.0.8_simple_v2 <- factor(lymphgland.combined@meta.data$LG_aligned_res.0.8_simple_v2,
                                                                     levels = rev(c('PSC', 'SC-like', 'PH', 'PM', 'LM', 'CC', 'GST-rich', 'DV', 'RG', 'Neurons')))
head(lymphgland.combined@meta.data)
unique(lymphgland.combined@meta.data$subcluster_anno)

DimPlot(lymphgland.combined, group.by = 'subcluster_anno', cells = rownames(subset(lymphgland.combined@meta.data, timepoint == 'Normal')), label = T) + facet_wrap(~subcluster_anno)

refanno <- subset(lymphgland.combined@meta.data, timepoint == 'Normal')
refanno <- droplevels(refanno)
components <- data.frame(row.names = levels(ref_comp$subcluster_anno))
for (anno in unique(refanno$LG_aligned_res.0.8_anno_v2)){
  ref_comp <- subset(refanno, LG_aligned_res.0.8_anno_v2 == anno)
  components <- cbind(components, summary(ref_comp$subcluster_anno))
}
colnames(components) <- unique(refanno$LG_aligned_res.0.8_anno_v2)
head(components)


### Transfer ###
lg96_normal
head(lg96_normal@meta.data)
lg96_infested
head(lg96_infested@meta.data)

transfer.anchors <- FindTransferAnchors(reference = lg96_normal, query = lg96_infested, dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = lg96_normal$subcluster_anno, dims = 1:30)
lg96_infested <- AddMetaData(lg96_infested, metadata = predictions)
head(lg96_infested@meta.data)
identical(rownames(lg96_infested@meta.data), rownames(subset(lymphgland.combined@meta.data, timepoint == 'Infested')))

lymphgland.combined@meta.data$labelTransfer <- as.character(lymphgland.combined@meta.data$subcluster_anno)
lymphgland.combined@meta.data[rownames(lg96_infested@meta.data), 'labelTransfer'] <- lg96_infested@meta.data$predicted.id
unique(lymphgland.combined@meta.data$labelTransfer)

lymphgland.combined@meta.data$labelTransfer <- mapvalues(lymphgland.combined@meta.data$labelTransfer,
                                                         from = unique(lymphgland.combined@meta.data$labelTransfer),
                                                         to = c("PSC", "PH 2", "PH 1", "Neurons", "PH 10", "PH 11", "PH 3", "PH 7", "PH 4", "PH 8", "PH 5", "PH 6", 
                                                                "PM 4", "PM 3", "PM 5", "PM 2", "PM 1", "PM 6", "PM 10", "RG", "LM 1", "LM 2", "CC 2", "CC 1", "GST-rich", "Adipohemocyte", "DV"))
lymphgland.combined@meta.data$labelTransfer <- factor(lymphgland.combined@meta.data$labelTransfer, 
                                                      levels = c("PSC", "PH 1", "PH 2", "PH 3", "PH 4", "PH 5", "PH 6", "PH 7", "PH 8", "PH 10", "PH 11", 
                                                                 "PM 1", "PM 2", "PM 3", "PM 4", "PM 5", "PM 6", "PM 10",  "LM 1", "LM 2", "CC 1", "CC 2", "GST-rich", "Adipohemocyte",
                                                                 "DV", "RG", "Neurons"))


lymphgland.combined@meta.data$labelTransfer_simple <- mapvalues(lymphgland.combined@meta.data$labelTransfer,
                                                                from = unique(lymphgland.combined@meta.data$labelTransfer),
                                                                to = c("PSC", "PH", "PH", "Neurons", "PH", "PH", "PH", "PH", "PH", "PH", "PH", "PH", 
                                                                       "PM", "PM", "PM", "PM", "PM", "PM", "PM", "RG", "LM", "LM", "CC", "CC", "GST-rich", "Adipohemocyte", "DV"))
lymphgland.combined@meta.data$labelTransfer_simple <- factor(lymphgland.combined@meta.data$labelTransfer_simple, 
                                                             levels = c("PSC", "PH", "PM", "LM", "CC", "GST-rich", "Adipohemocyte", "DV", "RG", "Neurons"))


DimPlot(lymphgland.combined, group.by = 'labelTransfer', label = T, dims = c(1,2))# + facet_wrap(~labelTransfer)
summary(lymphgland.combined@meta.data$labelTransfer)
DimPlot(lymphgland.combined, group.by = 'labelTransfer_simple', label = T, dims = c(1,2))# + facet_wrap(~labelTransfer)
summary(lymphgland.combined@meta.data$labelTransfer_simple)

lymphgland.combined@meta.data$LG_aligned_res.0.8_anno <- NULL
lymphgland.combined@meta.data$LG_aligned_res.0.8_simple <- NULL
head(lymphgland.combined@meta.data)

writeLable <- as.matrix(lymphgland.combined@meta.data)
writeLable <- data.frame(Barcode = rownames(writeLable), writeLable, check.names = F)
#write.table(writeLable, 'seurat3.96lgs.mitoCut.label.allGenes.normalized.ver5.labelTransferred.txt', sep = '\t', quote = F, row.names = F, col.names = T)
#save.image('Seurat3.alignment.96LGs.mitoCut.labelTransferred.Rdata')


##############################################################
##############################################################
##############################################################
library(Seurat)
library(ggplot2)
library(Matrix)
library(cowplot)
library(plyr)
library(dplyr)
library(RColorBrewer)

load('Seurat3.alignment.96LGs.mitoCut.labelTransferred.Rdata')

head(lymphgland.combined@meta.data)
Idents(lymphgland.combined) <- 'labelTransfer'

### Simple clustering ###
DimPlot(lymphgland.combined, group.by = 'labelTransfer_simple', label = T, label.size = 2, cols = c('#f18eb6', '#fc5c1b', '#207eb3', '#a80d0c', '#f0a142', '#5fccbc', '#a4a4a4', '#1a1a1a', '#c259aa', '#54040b', '#0d611b'))
#ggsave('labelTransferred/umap.seed96711161.mitoCut.1_2.min_dist_0.4_.labelTransfer_simple.pdf', units = 'cm', width = 16, height = 11)
plt <- DimPlot(lymphgland.combined, group.by = 'labelTransfer_simple', cols = c('#f18eb6', '#fc5c1b', '#207eb3', '#a80d0c', '#f0a142', '#5fccbc', '#a4a4a4', '#1a1a1a', '#c259aa', '#54040b', '#0d611b')) +
  theme_void(base_size = 0) + theme(legend.position = 'None');plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('labelTransferred/umap.seed96711161.mitoCut.1_2.min_dist_0.4_.labelTransfer_simple.augment.pdf', units = 'cm', width = 4, height = 4)

plt <- DimPlot(lymphgland.combined, group.by = 'labelTransfer_simple',
               cells = rownames(subset(lymphgland.combined@meta.data, timepoint == 'Normal')),
               cols = c('#f18eb6', '#fc5c1b', '#207eb3', '#a80d0c', '#f0a142', '#5fccbc', '#a4a4a4', '#1a1a1a', '#c259aa', '#54040b', '#0d611b')) +
  theme_void(base_size = 0) + theme(legend.position = 'None');plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('labelTransferred/umap.seed96711161.mitoCut.1_2.min_dist_0.4_.labelTransfer_simple.normal.augment.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(lymphgland.combined, group.by = 'labelTransfer_simple',
               cells = rownames(subset(lymphgland.combined@meta.data, timepoint == 'Infested')),
               cols = c('#f18eb6', '#fc5c1b', '#207eb3', '#a80d0c', '#f0a142', '#5fccbc', '#a4a4a4', '#c259aa', '#54040b', '#0d611b')) +
  theme_void(base_size = 0) + theme(legend.position = 'None');plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('labelTransferred/umap.seed96711161.mitoCut.1_2.min_dist_0.4_.labelTransfer_simple.infested.augment.pdf', units = 'cm', width = 4, height = 4)


### Subclustering ###
DimPlot(lymphgland.combined, group.by = 'labelTransfer', label = T, label.size = 2,
        cols = c('#f18eb6', '#fc5c1b', '#cbe0c3', '#b9dcb1', '#9bcba3', '#78c0ac', '#54b0ba', '#359bb7', '#257fb2', '#224a8d', '#1c306d', 
                 '#f5e7c9', '#ecd5a5', '#eac487', '#f6ae72', '#f38956', '#ea6740', '#901407',
                 '#e7c693', '#e18d30', '#71c0b0', '#177e7d', 
                 '#a4a4a4', '#1a1a1a', '#c259aa', '#54040b', '#0d611b'))
#ggsave('labelTransferred/umap.seed96711161.mitoCut.1_2.min_dist_0.4_.labelTransfer.pdf', units = 'cm', width = 19, height = 11)

### PH ###
lymphgland.combined@meta.data$labelTransfer_pg <- as.character(lymphgland.combined@meta.data$labelTransfer)
lymphgland.combined@meta.data[rownames(subset(lymphgland.combined@meta.data, labelTransfer_simple != 'PSC' & labelTransfer_simple != 'SC-like' & labelTransfer_simple != 'PH')), 'labelTransfer_pg'] <- 'Others'
unique(lymphgland.combined@meta.data$labelTransfer_pg)
lymphgland.combined@meta.data$labelTransfer_pg <- factor(lymphgland.combined@meta.data$labelTransfer_pg, 
                                                         levels = c("PSC", "PH 1", "PH 2", "PH 3", "PH 4", "PH 5", "PH 6", "PH 7", "PH 8", "PH 10", "PH 11", 'Others'))

DimPlot(lymphgland.combined, group.by = 'labelTransfer_pg', label = T, label.size = 2,
        cells = rownames(subset(lymphgland.combined@meta.data, timepoint == 'Normal')),
        cols = c('#f18eb6', '#fc5c1b', '#cbe0c3', '#b9dcb1', '#9bcba3', '#78c0ac', '#54b0ba', '#359bb7', '#257fb2', '#224a8d', '#1c306d', 'grey90'))
#ggsave('labelTransferred/umap.seed96711161.mitoCut.1_2.min_dist_0.4_.labelTransfer_pg.normal.pdf', units = 'cm', width = 14, height = 11)
plt <- DimPlot(lymphgland.combined, group.by = 'labelTransfer_pg',
               cells = rownames(subset(lymphgland.combined@meta.data, timepoint == 'Normal')),
               cols = c('#f18eb6', '#fc5c1b', '#cbe0c3', '#b9dcb1', '#9bcba3', '#78c0ac', '#54b0ba', '#359bb7', '#257fb2', '#224a8d', '#1c306d', 'grey90')) +
  theme_void(base_size = 0) + theme(legend.position = 'None');plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('labelTransferred/umap.seed96711161.mitoCut.1_2.min_dist_0.4_.labelTransfer_pg.normal.augment.pdf', units = 'cm', width = 4, height = 4)

DimPlot(lymphgland.combined, group.by = 'labelTransfer_pg', label = T, label.size = 2,
        cells = rownames(subset(lymphgland.combined@meta.data, timepoint == 'Infested')),
        cols = c('#f18eb6', '#fc5c1b', '#cbe0c3', '#b9dcb1', '#9bcba3', '#78c0ac', '#54b0ba', '#359bb7', '#257fb2', '#224a8d', 'grey90'))
#ggsave('labelTransferred/umap.seed96711161.mitoCut.1_2.min_dist_0.4_.labelTransfer_pg.infested.pdf', units = 'cm', width = 14, height = 11)
plt <- DimPlot(lymphgland.combined, group.by = 'labelTransfer_pg',
               cells = rownames(subset(lymphgland.combined@meta.data, timepoint == 'Infested')),
               cols = c('#f18eb6', '#fc5c1b', '#cbe0c3', '#b9dcb1', '#9bcba3', '#78c0ac', '#54b0ba', '#359bb7', '#257fb2', '#224a8d', 'grey90')) +
  theme_void(base_size = 0) + theme(legend.position = 'None');plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('labelTransferred/umap.seed96711161.mitoCut.1_2.min_dist_0.4_.labelTransfer_pg.infested.augment.pdf', units = 'cm', width = 4, height = 4)



### PM ###
lymphgland.combined@meta.data$labelTransfer_pm <- as.character(lymphgland.combined@meta.data$labelTransfer)
lymphgland.combined@meta.data[rownames(subset(lymphgland.combined@meta.data, labelTransfer_simple != 'PM')), 'labelTransfer_pm'] <- 'Others'
unique(lymphgland.combined@meta.data$labelTransfer_pm)
lymphgland.combined@meta.data$labelTransfer_pm <- factor(lymphgland.combined@meta.data$labelTransfer_pm, 
                                                         levels = c("PM 1", "PM 2", "PM 3", "PM 4", "PM 5", "PM 6", "PM 10", 'Others'))

DimPlot(lymphgland.combined, group.by = 'labelTransfer_pm', label = T, label.size = 2,
        cells = rownames(subset(lymphgland.combined@meta.data, timepoint == 'Normal')),
        cols = c('#f5e7c9', '#ecd5a5', '#eac487', '#f6ae72', '#f38956', '#ea6740', '#901407', 'grey90'))
#ggsave('labelTransferred/umap.seed96711161.mitoCut.1_2.min_dist_0.4_.labelTransfer_pm.normal.pdf', units = 'cm', width = 14, height = 11)
plt <- DimPlot(lymphgland.combined, group.by = 'labelTransfer_pm',
               cells = rownames(subset(lymphgland.combined@meta.data, timepoint == 'Normal')),
               cols = c('#f5e7c9', '#ecd5a5', '#eac487', '#f6ae72', '#f38956', '#ea6740', '#901407', 'grey90')) +
  theme_void(base_size = 0) + theme(legend.position = 'None');plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('labelTransferred/umap.seed96711161.mitoCut.1_2.min_dist_0.4_.labelTransfer_pm.normal.augment.pdf', units = 'cm', width = 4, height = 4)

DimPlot(lymphgland.combined, group.by = 'labelTransfer_pm', label = T, label.size = 2,
        cells = rownames(subset(lymphgland.combined@meta.data, timepoint == 'Infested')),
        cols = c('#f5e7c9', '#ecd5a5', '#eac487', '#f6ae72', '#f38956', '#ea6740', 'grey90'))
#ggsave('labelTransferred/umap.seed96711161.mitoCut.1_2.min_dist_0.4_.labelTransfer_pm.infested.pdf', units = 'cm', width = 14, height = 11)
plt <- DimPlot(lymphgland.combined, group.by = 'labelTransfer_pm',
               cells = rownames(subset(lymphgland.combined@meta.data, timepoint == 'Infested')),
               cols = c('#f5e7c9', '#ecd5a5', '#eac487', '#f6ae72', '#f38956', '#ea6740', 'grey90')) +
  theme_void(base_size = 0) + theme(legend.position = 'None');plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('labelTransferred/umap.seed96711161.mitoCut.1_2.min_dist_0.4_.labelTransfer_pm.infested.augment.pdf', units = 'cm', width = 4, height = 4)



### Others ###
lymphgland.combined@meta.data$labelTransfer_others <- as.character(lymphgland.combined@meta.data$labelTransfer)
lymphgland.combined@meta.data[rownames(subset(lymphgland.combined@meta.data, labelTransfer_simple != 'LM' & labelTransfer_simple != 'CC' & labelTransfer_simple != 'GST-rich' & 
                                                labelTransfer_simple != 'DV' & labelTransfer_simple != 'RG' & labelTransfer_simple != 'Neurons' & labelTransfer_simple != 'Adipohemocyte')), 'labelTransfer_others'] <- 'Others'
unique(lymphgland.combined@meta.data$labelTransfer_others)
lymphgland.combined@meta.data$labelTransfer_others <- factor(lymphgland.combined@meta.data$labelTransfer_others, 
                                                             levels = c("LM 1", "LM 2", "CC 1", "CC 2", "GST-rich", "Adipohemocyte", "DV", "RG", "Neurons", 'Others'))

DimPlot(lymphgland.combined, group.by = 'labelTransfer_others', label = T, label.size = 2,
        cells = rownames(subset(lymphgland.combined@meta.data, timepoint == 'Normal')),
        cols = c('#e7c693', '#e18d30', '#71c0b0', '#177e7d', '#a4a4a4', '#1a1a1a', '#c259aa', '#54040b', '#0d611b', 'grey90'))
#ggsave('labelTransferred/umap.seed96711161.mitoCut.1_2.min_dist_0.4_.labelTransfer_others.normal.pdf', units = 'cm', width = 14, height = 11)
plt <- DimPlot(lymphgland.combined, group.by = 'labelTransfer_others',
               cells = rownames(subset(lymphgland.combined@meta.data, timepoint == 'Normal')),
               cols = c('#e7c693', '#e18d30', '#71c0b0', '#177e7d', '#a4a4a4', '#1a1a1a', '#c259aa', '#54040b', '#0d611b', 'grey90')) +
  theme_void(base_size = 0) + theme(legend.position = 'None');plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('labelTransferred/umap.seed96711161.mitoCut.1_2.min_dist_0.4_.labelTransfer_others.normal.augment.pdf', units = 'cm', width = 4, height = 4)

DimPlot(lymphgland.combined, group.by = 'labelTransfer_others', label = T, label.size = 2,
        cells = rownames(subset(lymphgland.combined@meta.data, timepoint == 'Infested')),
        cols = c('#e7c693', '#e18d30', '#71c0b0', '#a4a4a4', '#c259aa', '#54040b', '#0d611b', 'grey90'))
#ggsave('labelTransferred/umap.seed96711161.mitoCut.1_2.min_dist_0.4_.labelTransfer_others.infested.pdf', units = 'cm', width = 14, height = 11)
plt <- DimPlot(lymphgland.combined, group.by = 'labelTransfer_others',
               cells = rownames(subset(lymphgland.combined@meta.data, timepoint == 'Infested')),
               cols = c('#e7c693', '#e18d30', '#71c0b0', '#a4a4a4', '#c259aa', '#54040b', '#0d611b', 'grey90')) +
  theme_void(base_size = 0) + theme(legend.position = 'None');plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('labelTransferred/umap.seed96711161.mitoCut.1_2.min_dist_0.4_.labelTransfer_others.infested.augment.pdf', units = 'cm', width = 4, height = 4)



### DEGs ###
#Idents(lymphgland.combined)
lymphgland.combined.markers <- FindAllMarkers(object = lymphgland.combined, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
lymphgland.combined.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(lymphgland.combined.markers, 'degs/findAllMarkers.labelTransfer.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
length(rownames(lymphgland.combined.markers))
top10 <- lymphgland.combined.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
head(top10)

dehm <- DoHeatmap(object = lymphgland.combined, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('degs/findAllMarkers.labelTransfer.pdf', units = 'cm', width = 60, height = 40)


Idents(lymphgland.combined) <- 'labelTransfer_simple'
lymphgland.combined.markers <- FindAllMarkers(object = lymphgland.combined, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
lymphgland.combined.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(lymphgland.combined.markers, 'degs/findAllMarkers.labelTransfer_simple.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
length(rownames(lymphgland.combined.markers))
top10 <- lymphgland.combined.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
head(top10)

dehm <- DoHeatmap(object = lymphgland.combined, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('degs/findAllMarkers.labelTransfer_simple.pdf', units = 'cm', width = 40, height = 30)
