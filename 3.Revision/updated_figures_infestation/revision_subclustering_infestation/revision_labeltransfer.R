library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)
library(RColorBrewer)

#dir.create('tmp')
#load('../Drop-seq_96hr_LymphGland_alignment/Seurat3.alignment.96LGs.mitoCut.Rdata')
#saveRDS(lg96_normal, 'tmp/lg96_normal.Rds')
#saveRDS(lg96_infested, 'tmp/lg96_infested.Rds')
#saveRDS(lymphgland.combined, 'tmp/lymphgland.combined.Rds')

lymphgland.combined <- readRDS('tmp/lymphgland.combined.Rds')
lg96_normal <- readRDS('tmp/lg96_normal.Rds')
lg96_infested <- readRDS('tmp/lg96_infested.Rds')
newlabel <- readRDS('../../Project1_LymphGland/revision_subclustering/updated_subclustering/Normal_LG/rdata/label.Rds')
newlabel <- newlabel[rownames(lg96_normal@meta.data), ]
head(newlabel)

### Transfer ###
lg96_normal
lg96_normal@meta.data <- lg96_normal@meta.data[, c(1:5, 9)]
lg96_normal@meta.data$newsubclustering <- newlabel$new_subclustering
head(lg96_normal@meta.data)

lg96_infested
lg96_infested@meta.data <- lg96_infested@meta.data[, c(1:5, 9)]
lg96_infested@meta.data$newsubclustering <- 'ns'
head(lg96_infested@meta.data)


transfer.anchors <- FindTransferAnchors(reference = lg96_normal, query = lg96_infested, dims = 1:30)
saveRDS(transfer.anchors, 'tmp/transfer.anchors.Rds')
predictions <- TransferData(anchorset = transfer.anchors, refdata = lg96_normal@meta.data$newsubclustering, dims = 1:30)
saveRDS(predictions, 'tmp/predictions.Rds')
head(predictions)
lg96_infested@meta.data$newsubclustering <- as.character(predictions$predicted.id)
#saveRDS(lg96_infested, 'tmp/lg96_infested_transferred.Rds')


head(lymphgland.combined@meta.data)
lymphgland.combined@meta.data <- lymphgland.combined@meta.data[, c(1:5, 9)]

lymphgland.combined@meta.data$labelTransfer <- NA
lymphgland.combined@meta.data[rownames(lg96_normal@meta.data), 'labelTransfer'] <- as.character(lg96_normal@meta.data$newsubclustering)
lymphgland.combined@meta.data[rownames(lg96_infested@meta.data), 'labelTransfer'] <- lg96_infested@meta.data$newsubclustering


lymphgland.combined@meta.data$labelTransfer <- factor(lymphgland.combined@meta.data$labelTransfer,
                                                      levels = c('PSC', 'PH 1', 'PH 2', 'PH 3', 'PH 4',  'PH 6', 
                                                                 "PM 1", "PM 2", "PM 4", "LM 1", "LM 2", "CC 1", "CC 2",
                                                                 'GST-rich', 'Adipohemocyte', 'DV', 'RG', 'Neurons'))
levels(lymphgland.combined@meta.data$labelTransfer)

lymphgland.combined@meta.data$labelTransfer_simple <- mapvalues(lymphgland.combined@meta.data$labelTransfer,
                                                                from = levels(lymphgland.combined@meta.data$labelTransfer),
                                                                to = c('PSC', 'PH', 'PH', 'PH', 'PH', 'PH', 'PM', 'PM', 'PM', 'LM', 'LM', 
                                                                       'CC', 'CC', 'GST-rich', 'Adipohemocyte', 'DV', 'RG', 'Neurons'))
levels(lymphgland.combined@meta.data$labelTransfer_simple)



DimPlot(lymphgland.combined, group.by = 'labelTransfer', label = T, dims = c(1,2))# + facet_wrap(~labelTransfer)
summary(lymphgland.combined@meta.data$labelTransfer)
DimPlot(lymphgland.combined, group.by = 'labelTransfer_simple', label = T, dims = c(1,2))# + facet_wrap(~labelTransfer)
summary(lymphgland.combined@meta.data$labelTransfer_simple)

head(lymphgland.combined@meta.data)
writeLable <- as.matrix(lymphgland.combined@meta.data)
writeLable <- data.frame(Barcode = rownames(writeLable), writeLable, check.names = F)
#write.table(writeLable, 'lymphgland.combined.label.txt', sep = '\t', quote = F, row.names = F, col.names = T)
#saveRDS(lymphgland.combined, 'tmp/lymphgland.combined_transferred.Rds')



#############################
### minor cell filtration ###
#############################
nrow(lymphgland.combined@meta.data)*0.001 #19.562
summary(lymphgland.combined@meta.data$labelTransfer)

#nrow(subset(lymphgland.combined@meta.data, timepoint == 'Normal'))*0.001
#summary(subset(lymphgland.combined@meta.data, timepoint == 'Normal')$labelTransfer)
#nrow(subset(lymphgland.combined@meta.data, timepoint == 'Infested'))*0.001
#summary(subset(lymphgland.combined@meta.data, timepoint == 'Infested')$labelTransfer)

Idents(lymphgland.combined) <- 'labelTransfer'
head(lymphgland.combined@meta.data)
lymphgland.combined <- subset(lymphgland.combined, idents = c('PSC', 'PH 1', 'PH 2', 'PH 3', 'PH 4', 'PM 1', 'LM 1', 'LM 2', 'CC 1', 'CC 2', 'GST-rich'))
#saveRDS(lymphgland.combined, 'tmp/lymphgland.combined_transferred_filtered.Rds')
#saveRDS(lymphgland.combined@meta.data, 'tmp/infested_label.Rds')


lymphgland.combined <- readRDS('tmp/lymphgland.combined_transferred_filtered.Rds')


#dir.create('labelTransferred')
plt <- DimPlot(lymphgland.combined, group.by = 'labelTransfer_simple',
               cells = rownames(subset(lymphgland.combined@meta.data, timepoint == 'Normal')),
               cols = c('#EB76A7', '#1C6AA4', '#95000C', '#EA9034', '#51C3AE', '#A4A4A4')) +
  theme_void(base_size = 0) + theme(legend.position = 'None');plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('labelTransferred/labelTransfer_simple.normal.augment.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(lymphgland.combined, group.by = 'labelTransfer_simple',
               cells = rownames(subset(lymphgland.combined@meta.data, timepoint == 'Infested')),
               cols = c('#EB76A7', '#1C6AA4', '#95000C', '#EA9034', '#51C3AE', '#A4A4A4')) +
  theme_void(base_size = 0) + theme(legend.position = 'None');plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('labelTransferred/labelTransfer_simple.infested.augment.pdf', units = 'cm', width = 4, height = 4)


plt <- DimPlot(lymphgland.combined, group.by = 'labelTransfer',
               cells = rownames(subset(lymphgland.combined@meta.data, timepoint == 'Normal')),
               cols = c('#EB76A7', '#dee6d4', '#b9dcb1', '#78c0ac', '#359bb7', '#95000C',
                        '#E7C693', '#E18D30', '#71C0B0', '#177E7D', '#A4A4A4')) +
  theme_void(base_size = 0) + theme(legend.position = 'None');plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('labelTransferred/labelTransfer.normal.augment.pdf', units = 'cm', width = 4, height = 4)
plt <- DimPlot(lymphgland.combined, group.by = 'labelTransfer',
               cells = rownames(subset(lymphgland.combined@meta.data, timepoint == 'Infested')),
               cols = c('#EB76A7', '#dee6d4', '#b9dcb1', '#78c0ac', '#359bb7', '#95000C',
                        '#E7C693', '#E18D30', '#71C0B0', '#A4A4A4')) +
  theme_void(base_size = 0) + theme(legend.position = 'None');plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('labelTransferred/labelTransfer.infested.augment.pdf', units = 'cm', width = 4, height = 4)







