library(Seurat)
library(plyr)
library(ggplot2)
library(reshape2)

lg96 <- read.delim('96hAEL_label_0.7.txt')
head(lg96)


### timepoint clustering ###
lg96$supergroup_new_subclustering <- mapvalues(lg96$supergroup, from = unique(lg96$supergroup), to = c('PSC', 'PH1/2', 'PH', 'PH3', 'PH_cellcycle', 'PM', 'GST-rich', 'CC'))
lg96$supergroup_new_subclustering <- factor(lg96$supergroup_new_subclustering, levels = c('PSC', 'PH1/2', 'PH3', 'PH_cellcycle', 'PH', 'PM', 'CC', 'GST-rich'))
lg96$supergroup_new_celltype <- mapvalues(lg96$supergroup_new_subclustering, from = levels(lg96$supergroup_new_subclustering),
                                          to = c('PSC', 'PH', 'PH', 'PH', 'PH', 'PM', 'CC', 'GST-rich'))

### new subcluster and cell type ###
lg96$new_subclustering <- factor(lg96$new_subclustering,
                                 levels = c('PSC', 'PH 1', 'PH 2', 'PH 3', 'PH 4', 'PH 6', 'PM 1', 'PM 2', 'PM 4', 'LM 1', 'LM 2', 'CC 1', 'CC 2', 'GST-rich', 'Adipohemocyte'))
lg96$new_celltype <- mapvalues(lg96$new_subclustering,
                               from = levels(lg96$new_subclustering),
                               to = c('PSC', 'PH', 'PH', 'PH', 'PH', 'PH', 'PM', 'PM', 'PM', 'LM', 'LM', 'CC', 'CC', 'GST-rich', 'Adipohemocyte'))


lg96 <- arrange(lg96, lg96$supergroup_new_celltype)
#lg96$supergroup_order <- c(1:nrow(lg96))
lg96$supergroup_order <- NA
count <- 0
for (i in c(1:length(levels(lg96$supergroup_new_celltype)))) {
  bcs <- rownames(subset(lg96, supergroup_new_celltype==levels(lg96$supergroup_new_celltype)[i]))
  lg96[bcs, 'supergroup_order'] <- as.numeric(bcs) + count
  count <- count + 250
}


head(lg96)
lg96 <- arrange(lg96, lg96$new_celltype)
#lg96$new_subclustering_order <- c(1:nrow(lg96))
lg96$new_subclustering_order <- NA
count <- 0
for (i in c(1:length(levels(lg96$new_celltype)))) {
  bcs <- rownames(subset(lg96, new_celltype==levels(lg96$new_celltype)[i]))
  lg96[bcs, 'new_subclustering_order'] <- as.numeric(bcs) + count
  count <- count + 250
}
head(lg96)

lg96_m <- melt(lg96[, c(1, 15:18)])
head(lg96_m)
summary(lg96_m)
lg96_m$variable <- mapvalues(lg96_m$variable, from = levels(lg96_m$variable), to = c('96h AEL', 'Original'))

ggplot(lg96_m, aes(variable, value, group = Barcode)) + 
  geom_point(data = subset(lg96_m, variable == 'Original'), aes(variable, value, col = new_celltype)) +
  geom_point(data = subset(lg96_m, variable == '96h AEL'), aes(variable, value, col = supergroup_new_celltype)) +
  scale_color_manual(values = c('#F15FA6', '#207EB3', '#A80D0C', '#F0A142', '#25A9B0', '#A4A4A4', '#1A1A1A')) +
  geom_line(size=.5, alpha=0.01) +
  theme_void() +
  theme(axis.ticks = element_blank(), 
        axis.text = element_text(colour = 'black'), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), 
        axis.text.y = element_blank())
ggsave('96hAEL_label_0.7_matchCelltype.pdf', units = 'cm', width = 8, height = 20)



library(Seurat)
library(plyr)
library(ggplot2)

lg96 <- read.delim('label_0.7.txt')
lg96$supergroup_new_subclustering <- mapvalues(lg96$supergroup, from = unique(lg96$supergroup), to = c('PSC', 'PH1/2', 'PH', 'PH3', 'PH_cellcycle', 'PM', 'GST-rich', 'CC'))
lg96$supergroup_new_subclustering <- factor(lg96$supergroup_new_subclustering, levels = c('PSC', 'PH1/2', 'PH3', 'PH_cellcycle', 'PH', 'PM', 'CC', 'GST-rich'))
lg96$supergroup_new_celltype <- mapvalues(lg96$supergroup_new_subclustering, from = levels(lg96$supergroup_new_subclustering),
                                          to = c('PSC', 'PH', 'PH', 'PH', 'PH', 'PM', 'CC', 'GST-rich'))
head(lg96)

lymphgland <- readRDS('../../Normal_LG/rdata/lymphgland.Rds')

lymphgland_96 <- subset(lymphgland, timepoint == 'AEL96hr')
lymphgland_96@meta.data <- droplevels(lymphgland_96@meta.data)
head(lymphgland_96@meta.data)

rownames(lg96) <- lg96$Barcode
lymphgland_96@meta.data$tp_cluster <- lymphgland_96@meta.data$anno_simple
lymphgland_96@meta.data[rownames(lg96), 'tp_cluster'] <- lg96$supergroup_new_celltype
lymphgland_96@meta.data <- droplevels(lymphgland_96@meta.data)
#writetable <- data.frame(Barcode = rownames(lymphgland_96@meta.data), lymphgland_96@meta.data, check.rows = F, check.names = F)
#write.table(writetable, 'label_0.7_96hAELmeta.txt', quote = F, sep = '\t', row.names = F, col.names = T)


# '#f179a8', '#1d91c5', '#b9060a', '#f5b053', '#a4a4a4', '#cf77e0', '#975a1a', '#bc2aff' 
plt <- DimPlot(lymphgland_96, reduction = 'tsne', group.by = 'tp_cluster') +
  scale_color_manual(values = c('#f179a8', '#1d91c5', '#b9060a', '#23c7bc', '#a4a4a4', '#cf77e0', '#975a1a', '#bc2aff' )); plt
ggsave('Fig1D.tsne.tpcluster.pdf', units = 'cm', width = 11, height = 8)
AugmentPlot(plt + theme_void() + theme(legend.position = 'None'), width = 3, height = 3, dpi = 300)
ggsave('Fig1D.tsne.tpcluster.augment.pdf', units = 'cm', width = 3, height = 3)


summary(lymphgland_96@meta.data$tp_cluster)