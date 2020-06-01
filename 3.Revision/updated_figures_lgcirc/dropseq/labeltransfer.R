library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

dir.create('stats')
dir.create('tmp')
dir.create('tsne')
dir.create('umap')

### Circulation 96 Drop-seq ###
circdrop96_expr <- read.delim('../../dropseq_circulation_20190917/tables/merged.exprs.circ.96AEL.txt', header = T, sep = '\t', check.names = F, row.names = 1)
circdrop96_label <- read.delim('../../dropseq_circulation_20190917/96hAEL/seurat3/dropseq.mtcut.label.txt', header = T, sep = '\t', check.names = F, row.names = 1)
head(circdrop96_label)
circdrop96_label <- circdrop96_label[, c(4, 1)]
colnames(circdrop96_label) <- c('Library', 'timepoint')
circdrop96_label$timepoint <- 'Circ_96'
circdrop96_label$anno_simple <- 'Circ_96'
circdrop96_label$Subclustering <- 'Circ_96'
circdrop96_expr <- circdrop96_expr[, rownames(circdrop96_label)]
circdrop96_expr <- as(as.matrix(circdrop96_expr), "dgCMatrix")
### Circulation 120 Drop-seq ###
circdrop120_expr <- read.delim('../../dropseq_circulation_20190917/tables/merged.exprs.circ.120AEL.txt', header = T, sep = '\t', check.names = F, row.names = 1)
circdrop120_label <- read.delim('../../dropseq_circulation_20190917/120hAEL/seurat3/dropseq.mtcut.label.txt', header = T, sep = '\t', check.names = F, row.names = 1)
circdrop120_label <- circdrop120_label[, c(4, 1)]
colnames(circdrop120_label) <- c('Library', 'timepoint')
circdrop120_label$timepoint <- 'Circ_120'
circdrop120_label$anno_simple <- 'Circ_120'
circdrop120_label$Subclustering <- 'Circ_120'
circdrop120_expr <- circdrop120_expr[, rownames(circdrop120_label)]
circdrop120_expr <- as(as.matrix(circdrop120_expr), "dgCMatrix")


### lymph gland ###
lymphgland_expr <- read.delim('../../../Project1_LymphGland/Drop-seq_alignment/tables/merged.3tps.expr.allGenes.txt', header = T, sep = '\t', check.names = F, row.names = 1)
lymphgland_label <- readRDS('../../../Project1_LymphGland/revision_subclustering/updated_subclustering/Normal_LG/rdata/label.Rds')
lymphgland_label <- subset(lymphgland_label, new_subclustering != 'DV' & new_subclustering != 'RG' & new_subclustering != 'Neurons'); lymphgland_label <- droplevels(lymphgland_label)
lymphgland_label <- lymphgland_label[, c(4, 6, 7, 9)] # Library, timepoint, anno_simple, new_subclustering
colnames(lymphgland_label)[4] <- 'Subclustering'
head(lymphgland_label)
# LG 96
lymphgland96_label <- subset(lymphgland_label, timepoint == 'AEL96hr')
lymphgland96_label$timepoint <- 'LG_96'
lymphgland96_expr <- lymphgland_expr[, rownames(lymphgland96_label)]
lymphgland96_expr <- as(as.matrix(lymphgland96_expr), "dgCMatrix")
# LG 120
lymphgland120_label <- subset(lymphgland_label, timepoint == 'AEL120hr')
lymphgland120_label$timepoint <- 'LG_120'
lymphgland120_expr <- lymphgland_expr[, rownames(lymphgland120_label)]
lymphgland120_expr <- as(as.matrix(lymphgland120_expr), "dgCMatrix")
remove(lymphgland_expr)
remove(lymphgland_label)


### Seurat objects for each dataset ###
circdrop96 <- CreateSeuratObject(counts = circdrop96_expr, project = "lg2circ")
circdrop96 <- AddMetaData(object = circdrop96, metadata = circdrop96_label)
circdrop96 <- NormalizeData(circdrop96, normalization.method = "LogNormalize", scale.factor = 10000)
circdrop96 <- FindVariableFeatures(object = circdrop96, selection.method = "vst", nfeatures = 2000)
saveRDS(circdrop96, 'tmp/circdrop96.Rds')
head(circdrop96@meta.data)

circdrop120 <- CreateSeuratObject(counts = circdrop120_expr, project = "lg2circ")
circdrop120 <- AddMetaData(object = circdrop120, metadata = circdrop120_label)
circdrop120 <- NormalizeData(circdrop120, normalization.method = "LogNormalize", scale.factor = 10000)
circdrop120 <- FindVariableFeatures(object = circdrop120, selection.method = "vst", nfeatures = 2000)
saveRDS(circdrop120, 'tmp/circdrop120.Rds')
head(circdrop120@meta.data)

lymphgland96 <- CreateSeuratObject(counts = lymphgland96_expr, project = "lg2circ")
lymphgland96 <- AddMetaData(object = lymphgland96, metadata = lymphgland96_label)
lymphgland96 <- NormalizeData(lymphgland96, normalization.method = "LogNormalize", scale.factor = 10000)
lymphgland96 <- FindVariableFeatures(object = lymphgland96, selection.method = "vst", nfeatures = 2000)
saveRDS(lymphgland96, 'tmp/lymphgland96.Rds')
head(lymphgland96@meta.data)

lymphgland120 <- CreateSeuratObject(counts = lymphgland120_expr, project = "lg2circ")
lymphgland120 <- AddMetaData(object = lymphgland120, metadata = lymphgland120_label)
lymphgland120 <- NormalizeData(lymphgland120, normalization.method = "LogNormalize", scale.factor = 10000)
lymphgland120 <- FindVariableFeatures(object = lymphgland120, selection.method = "vst", nfeatures = 2000)
saveRDS(lymphgland120, 'tmp/lymphgland120.Rds')
head(lymphgland120@meta.data)
save.image('tmp/tmp1.Rdata')


### Transfer labels ###
# LG_drop_120 -> Circ_drop_120
transfer.anchors <- FindTransferAnchors(reference = lymphgland120, query = circdrop120, dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = lymphgland120$Subclustering, dims = 1:30)
saveRDS(predictions, 'tmp/predictions_120.Rds')
circdrop120@meta.data$Subclustering <- predictions$predicted.id
summary(as.factor(circdrop120@meta.data$Subclustering))
# LG_drop_96 -> Circ_drop_96
transfer.anchors <- FindTransferAnchors(reference = lymphgland96, query = circdrop96, dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = lymphgland96$Subclustering, dims = 1:30)
saveRDS(predictions, 'tmp/predictions_96.Rds')
circdrop96@meta.data$Subclustering <- predictions$predicted.id
summary(as.factor(circdrop96@meta.data$Subclustering))
save.image('tmp/tmp2.Rdata')


###
blood.combined <- readRDS('../../dropseq_circulation_20190917/seurat_alignment/alignment_dropseq/tmp/blood.combined.Rds')
head(blood.combined@meta.data)

blood.combined@meta.data[rownames(lymphgland96@meta.data), 'Subclustering'] <- lymphgland96@meta.data$Subclustering
blood.combined@meta.data[rownames(lymphgland120@meta.data), 'Subclustering'] <- lymphgland120@meta.data$Subclustering
blood.combined@meta.data[rownames(circdrop96@meta.data), 'Subclustering'] <- circdrop96@meta.data$Subclustering
blood.combined@meta.data[rownames(circdrop120@meta.data), 'Subclustering'] <- circdrop120@meta.data$Subclustering
blood.combined@meta.data <- droplevels(blood.combined@meta.data)
summary(blood.combined@meta.data$Subclustering)
Idents(blood.combined) <- 'anno_simple'

saveRDS(blood.combined, 'tmp/blood.combined_newsubclustering.Rds')
save.image('tmp/tmp4.Rdata')


blood.combined <- readRDS('tmp/blood.combined_newsubclustering.Rds')
nrow(blood.combined@meta.data)
writeLabel <- data.frame(Barcode = rownames(blood.combined@meta.data), blood.combined@meta.data, check.rows = F, check.names = F)
write.table(writeLabel, 'label.txt', quote = F, sep = '\t', row.names = F, col.names = T)
