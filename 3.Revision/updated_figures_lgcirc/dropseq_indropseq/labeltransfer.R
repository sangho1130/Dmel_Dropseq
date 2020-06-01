library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

dir.create('stats')
dir.create('tmp')
dir.create('tsne')
dir.create('umap')

### Circulation inDrop-seq ###
circulation_expr <- read.delim('../../dropseq_circulation_20190917/sudhir_table/merged.uninjured.expr.symbol.txt', header = T, sep = '\t', check.names = F, row.names = 1)
circulation_label <- read.delim('../../dropseq_circulation_20190917/120hAEL_inDrop/seurat3/indropseq.mtcut.label.txt', header = T, sep = '\t', check.names = F, row.names = 1)
circulation_label <- circulation_label[, c(4, 1)]
colnames(circulation_label) <- c('Library', 'timepoint')
# uninjured_1
uninjured_1_label <- subset(circulation_label, Library == 'uninjured_1')
uninjured_1_label$timepoint <- 'Circ_indrop_120'
uninjured_1_label$anno_simple <- 'Circ_indrop_120'
uninjured_1_label$Subclustering <- 'Circ_indrop_120'
uninjured_1_expr <- circulation_expr[, rownames(uninjured_1_label)]
uninjured_1_expr <- as(as.matrix(uninjured_1_expr), "dgCMatrix")
# uninjured_2
uninjured_2_label <- subset(circulation_label, Library == 'uninjured_2')
uninjured_2_label$timepoint <- 'Circ_indrop_120'
uninjured_2_label$anno_simple <- 'Circ_indrop_120'
uninjured_2_label$Subclustering <- 'Circ_indrop_120'
uninjured_2_expr <- circulation_expr[, rownames(uninjured_2_label)]
uninjured_2_expr <- as(as.matrix(uninjured_2_expr), "dgCMatrix")
remove(circulation_expr)


### Circulation 96 Drop-seq ###
circdrop96_expr <- read.delim('../../dropseq_circulation_20190917/tables/merged.exprs.circ.96AEL.IDs.symbol.txt', header = T, sep = '\t', check.names = F, row.names = 1)
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
circdrop120_expr <- read.delim('../../dropseq_circulation_20190917/tables/merged.exprs.circ.120AEL.IDs.symbol.txt', header = T, sep = '\t', check.names = F, row.names = 1)
circdrop120_label <- read.delim('../../dropseq_circulation_20190917/120hAEL/seurat3/dropseq.mtcut.label.txt', header = T, sep = '\t', check.names = F, row.names = 1)
circdrop120_label <- circdrop120_label[, c(4, 1)]
colnames(circdrop120_label) <- c('Library', 'timepoint')
circdrop120_label$timepoint <- 'Circ_120'
circdrop120_label$anno_simple <- 'Circ_120'
circdrop120_label$Subclustering <- 'Circ_120'
circdrop120_expr <- circdrop120_expr[, rownames(circdrop120_label)]
circdrop120_expr <- as(as.matrix(circdrop120_expr), "dgCMatrix")


### lymph gland ###
lymphgland_expr <- read.delim('../../dropseq_circulation_20190917/tables/merged.3tps.expr.allGenes.IDs.symbol.txt', header = T, sep = '\t', check.names = F, row.names = 1)
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
uninjured_1 <- CreateSeuratObject(counts = uninjured_1_expr, project = "lg2circ")
uninjured_1 <- AddMetaData(object = uninjured_1, metadata = uninjured_1_label)
uninjured_1 <- NormalizeData(uninjured_1, normalization.method = "LogNormalize", scale.factor = 10000)
uninjured_1 <- FindVariableFeatures(object = uninjured_1, selection.method = "vst", nfeatures = 2000)
saveRDS(uninjured_1, 'tmp/uninjured_1.Rds')
head(uninjured_1@meta.data)

uninjured_2 <- CreateSeuratObject(counts = uninjured_2_expr, project = "lg2circ")
uninjured_2 <- AddMetaData(object = uninjured_2, metadata = uninjured_2_label)
uninjured_2 <- NormalizeData(uninjured_2, normalization.method = "LogNormalize", scale.factor = 10000)
uninjured_2 <- FindVariableFeatures(object = uninjured_2, selection.method = "vst", nfeatures = 2000)
saveRDS(uninjured_2, 'tmp/uninjured_2.Rds')
head(uninjured_2@meta.data)

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
# LG_drop_120 -> Circ_indrop_120-1
transfer.anchors <- FindTransferAnchors(reference = lymphgland120, query = uninjured_1, dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = lymphgland120$Subclustering, dims = 1:30)
saveRDS(predictions, 'tmp/predictions_indrop120_1.Rds')
uninjured_1@meta.data$Subclustering <- predictions$predicted.id
# LG_drop_120 -> Circ_indrop_120-2
transfer.anchors <- FindTransferAnchors(reference = lymphgland120, query = uninjured_2, dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = lymphgland120$Subclustering, dims = 1:30)
saveRDS(predictions, 'tmp/predictions_indrop120_2.Rds')
uninjured_2@meta.data$Subclustering <- predictions$predicted.id
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
blood.combined <- readRDS('../../dropseq_circulation_20190917/seurat_alignment/alignment/blood.combined.Rds')
head(blood.combined@meta.data)
identical(sort(rownames(blood.combined@meta.data)),
          sort(c(rownames(uninjured_1@meta.data), rownames(uninjured_2@meta.data),
                 rownames(circdrop120@meta.data), rownames(circdrop96@meta.data),
                 rownames(lymphgland96@meta.data), rownames(lymphgland120@meta.data))))

blood.combined@meta.data$anno_simple <- NA
blood.combined@meta.data$Subclustering <- NA
blood.combined@meta.data[rownames(uninjured_1@meta.data), 'Subclustering'] <- as.character(uninjured_1@meta.data$Subclustering)
blood.combined@meta.data[rownames(uninjured_2@meta.data), 'Subclustering'] <- as.character(uninjured_2@meta.data$Subclustering)
blood.combined@meta.data[rownames(lymphgland96@meta.data), 'Subclustering'] <- as.character(lymphgland96@meta.data$Subclustering)
blood.combined@meta.data[rownames(lymphgland120@meta.data), 'Subclustering'] <- as.character(lymphgland120@meta.data$Subclustering)
blood.combined@meta.data[rownames(circdrop96@meta.data), 'Subclustering'] <- as.character(circdrop96@meta.data$Subclustering)
blood.combined@meta.data[rownames(circdrop120@meta.data), 'Subclustering'] <- as.character(circdrop120@meta.data$Subclustering)
blood.combined@meta.data$Subclustering <- factor(blood.combined@meta.data$Subclustering, 
                                                 levels = levels(lymphgland120@meta.data$Subclustering))
summary(blood.combined@meta.data$Subclustering)

blood.combined@meta.data$anno_simple <- mapvalues(blood.combined@meta.data$Subclustering,
                                                  from = levels(blood.combined@meta.data$Subclustering),
                                                  to = c('PSC', 'PH', 'PH', 'PH', 'PH', 'PH', 'PH', 
                                                         'PM', 'PM', 'PM', 'PM', 'LM', 'LM', 'CC', 'CC', 'GST-rich', 'Adipohemocyte'))
summary(blood.combined@meta.data$anno_simple)

blood.combined@meta.data <- droplevels(blood.combined@meta.data)
Idents(blood.combined) <- 'anno_simple'

#saveRDS(blood.combined, 'tmp/blood.combined_newsubclustering.Rds')
#save.image('tmp/tmp4.Rdata')

blood.combined <- readRDS('tmp/blood.combined_newsubclustering.Rds')
nrow(blood.combined@meta.data)
writeLabel <- data.frame(Barcode = rownames(blood.combined@meta.data), blood.combined@meta.data, check.rows = F, check.names = F)
write.table(writeLabel, 'label.txt', quote = F, sep = '\t', row.names = F, col.names = T)


