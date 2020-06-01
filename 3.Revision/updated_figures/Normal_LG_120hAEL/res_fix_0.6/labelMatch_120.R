library(Seurat)
library(plyr)
library(ggplot2)
library(reshape2)

lg120 <- read.delim('120hAEL_label_0.6.txt')
head(lg120)

### timepoint clustering ###
lg120$supergroup_new_subclustering <- mapvalues(lg120$supergroup, from = unique(lg120$supergroup),
                                                to = c('PSC', 'PH1-3', 'PH4/5', 'PM3/4', 'PM1/2', 'LM', 'CC', 'GST-rich', 'PM2', 'Adipohemocyte', 'PH6'))
lg120$supergroup_new_subclustering <- factor(lg120$supergroup_new_subclustering, 
                                             levels = c('PSC', 'PH1-3', 'PH4/5', 'PH6', 'PM1/2', 'PM2', 'PM3/4', 'LM', 'CC', 'GST-rich', 'Adipohemocyte'))
lg120$supergroup_new_celltype <- mapvalues(lg120$supergroup_new_subclustering, from = levels(lg120$supergroup_new_subclustering),
                                          to = c('PSC', 'PH', 'PH', 'PH', 'PM', 'PM', 'PM', 'LM', 'CC', 'GST-rich', 'Adipohemocyte'))

### new subcluster and cell type ###
lg120$new_subclustering <- factor(lg120$new_subclustering,
                                 levels = c('PSC', 'PH 1', 'PH 2', 'PH 3', 'PH 4', 'PH 5', 'PH 6', 'PM 1', 'PM 2', 'PM 3', 'PM 4',
                                            'LM 1', 'LM 2', 'CC 1', 'CC 2', 'GST-rich', 'Adipohemocyte'))
lg120$new_celltype <- mapvalues(lg120$new_subclustering,
                                from = levels(lg120$new_subclustering),
                                to = c('PSC', 'PH', 'PH', 'PH', 'PH', 'PH', 'PH', 'PM', 'PM', 'PM', 'PM', 'LM', 'LM', 'CC', 'CC', 'GST-rich', 'Adipohemocyte'))


lg120 <- arrange(lg120, lg120$supergroup_new_celltype)
#lg120$supergroup_order <- c(1:nrow(lg120))
lg120$supergroup_order <- NA
count <- 0
for (i in c(1:length(levels(lg120$supergroup_new_celltype)))) {
  bcs <- rownames(subset(lg120, supergroup_new_celltype==levels(lg120$supergroup_new_celltype)[i]))
  lg120[bcs, 'supergroup_order'] <- as.numeric(bcs) + count
  count <- count + 250
}


head(lg120)
lg120 <- arrange(lg120, lg120$new_celltype)
#lg120$new_subclustering_order <- c(1:nrow(lg120))
lg120$new_subclustering_order <- NA
count <- 0
for (i in c(1:length(levels(lg120$new_celltype)))) {
  bcs <- rownames(subset(lg120, new_celltype==levels(lg120$new_celltype)[i]))
  lg120[bcs, 'new_subclustering_order'] <- as.numeric(bcs) + count
  count <- count + 250
}
head(lg120)

lg120_m <- melt(lg120[, c(1, 16:19)])
head(lg120_m)
summary(lg120_m)
lg120_m$variable <- mapvalues(lg120_m$variable, from = levels(lg120_m$variable), to = c('120h AEL', 'Original'))
ggplot(lg120_m, aes(variable, value, group = Barcode)) + 
  geom_point(data = subset(lg120_m, variable == 'Original'), aes(variable, value, col = new_celltype)) +
  geom_point(data = subset(lg120_m, variable == '120h AEL'), aes(variable, value, col = supergroup_new_celltype)) +
  scale_color_manual(values = c('#F15FA6', '#207EB3', '#A80D0C', '#F0A142', '#25A9B0', '#A4A4A4', '#1A1A1A')) +
  geom_line(size=.5, alpha=0.01) +
  theme_void() +
  theme(axis.ticks = element_blank(), 
        axis.text = element_text(colour = 'black'), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), 
        axis.text.y = element_blank())
ggsave('120hAEL_label_0.6_matchCelltype.pdf', units = 'cm', width = 8, height = 20)




library(Seurat)
library(plyr)
library(ggplot2)

lg120 <- read.delim('label_0.6.txt')
lg120$supergroup_new_subclustering <- mapvalues(lg120$supergroup, from = unique(lg120$supergroup),
                                                to = c('PSC', 'PH1-3', 'PH4/5', 'PM3/4', 'PM1/2', 'LM', 'CC', 'GST-rich', 'PM2', 'Adipohemocyte', 'PH6'))
lg120$supergroup_new_subclustering <- factor(lg120$supergroup_new_subclustering, 
                                             levels = c('PSC', 'PH1-3', 'PH4/5', 'PH6', 'PM1/2', 'PM2', 'PM3/4', 'LM', 'CC', 'GST-rich', 'Adipohemocyte'))
lg120$supergroup_new_celltype <- mapvalues(lg120$supergroup_new_subclustering, from = levels(lg120$supergroup_new_subclustering),
                                           to = c('PSC', 'PH', 'PH', 'PH', 'PM', 'PM', 'PM', 'LM', 'CC', 'GST-rich', 'Adipohemocyte'))
head(lg120)

lymphgland <- readRDS('../../Normal_LG/rdata/lymphgland.Rds')

lymphgland_120 <- subset(lymphgland, timepoint == 'AEL120hr')
lymphgland_120@meta.data <- droplevels(lymphgland_120@meta.data)
head(lymphgland_120@meta.data)

rownames(lg120) <- lg120$Barcode
lymphgland_120@meta.data$tp_cluster <- lymphgland_120@meta.data$anno_simple
lymphgland_120@meta.data[rownames(lg120), 'tp_cluster'] <- lg120$supergroup_new_celltype
lymphgland_120@meta.data <- droplevels(lymphgland_120@meta.data)
#writetable <- data.frame(Barcode = rownames(lymphgland_120@meta.data), lymphgland_120@meta.data, check.rows = F, check.names = F)
#write.table(writetable, 'label_0.6_120hAELmeta.txt', quote = F, sep = '\t', row.names = F, col.names = T)

plt <- DimPlot(lymphgland_120, reduction = 'tsne', group.by = 'tp_cluster') +
  scale_color_manual(values = c('#f179a8', '#1d91c5', '#b9060a', '#f5b053', '#23c7bc', '#a4a4a4', '#1a1a1a', '#cf77e0', '#975a1a')); plt
ggsave('Fig1D.tsne.tpcluster.pdf', units = 'cm', width = 12, height = 8)
AugmentPlot(plt + theme_void() + theme(legend.position = 'None'), width = 3, height = 3, dpi = 300)
ggsave('Fig1D.tsne.tpcluster.augment.pdf', units = 'cm', width = 3, height = 3)


summary(lymphgland_120@meta.data$tp_cluster)