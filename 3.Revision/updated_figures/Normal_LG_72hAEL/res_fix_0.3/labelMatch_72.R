library(Seurat)
library(plyr)
library(ggplot2)
library(reshape2)

lg72 <- read.delim('label_0.3.txt')
head(lg72)


### timepoint clustering ###
lg72$supergroup_new_subclustering <- mapvalues(lg72$supergroup, from = unique(lg72$supergroup), to = c('PSC', 'PH1-3', 'PH', 'PM', 'PM-NimC1'))
lg72$supergroup_new_subclustering <- factor(lg72$supergroup_new_subclustering, levels = c('PSC', 'PH1-3', 'PH', 'PM', 'PM-NimC1'))
lg72$supergroup_new_celltype <- mapvalues(lg72$supergroup_new_subclustering, from = levels(lg72$supergroup_new_subclustering),
                                          to = c('PSC', 'PH', 'PH', 'PM', 'PM'))

### new subcluster and cell type ###
unique(lg72$new_subclustering)
lg72$new_subclustering <- factor(lg72$new_subclustering,
                                 levels = c('PSC', 'PH 1', 'PH 2', 'PH 3', 'PH 4', 'PH 6', 'PM 1', 'PM 4', 'LM 1', 'CC 1', 'GST-rich', 'Adipohemocyte'))
lg72$new_celltype <- mapvalues(lg72$new_subclustering,
                               from = levels(lg72$new_subclustering),
                               to = c('PSC', 'PH', 'PH', 'PH', 'PH', 'PH', 'PM', 'PM', 'LM', 'CC', 'GST-rich', 'Adipohemocyte'))


lg72 <- arrange(lg72, lg72$supergroup_new_celltype)
#lg72$supergroup_order <- c(1:nrow(lg72))
lg72$supergroup_order <- NA
count <- 0
for (i in c(1:length(levels(lg72$supergroup_new_celltype)))) {
  bcs <- rownames(subset(lg72, supergroup_new_celltype==levels(lg72$supergroup_new_celltype)[i]))
  lg72[bcs, 'supergroup_order'] <- as.numeric(bcs) + count
  count <- count + 50
}


head(lg72)
lg72 <- arrange(lg72, lg72$new_celltype)
#lg72$new_subclustering_order <- c(1:nrow(lg72))
lg72$new_subclustering_order <- NA
count <- 0
for (i in c(1:length(levels(lg72$new_celltype)))) {
  bcs <- rownames(subset(lg72, new_celltype==levels(lg72$new_celltype)[i]))
  lg72[bcs, 'new_subclustering_order'] <- as.numeric(bcs) + count
  count <- count + 50
}
head(lg72)

lg72_m <- melt(lg72[, c(1, 15:18)])
head(lg72_m)
summary(lg72_m)
lg72_m$variable <- mapvalues(lg72_m$variable, from = levels(lg72_m$variable), to = c('72h AEL', 'Original'))

ggplot(lg72_m, aes(variable, value, group = Barcode)) + 
  geom_point(data = subset(lg72_m, variable == 'Original'), aes(variable, value, col = new_celltype)) +
  geom_point(data = subset(lg72_m, variable == '72h AEL'), aes(variable, value, col = supergroup_new_celltype)) +
  scale_color_manual(values = c('#F15FA6', '#207EB3', '#A80D0C', '#F0A142', '#25A9B0', '#A4A4A4', '#1A1A1A')) +
  geom_line(size=.5, alpha=0.1) +
  theme_void() +
  theme(axis.ticks = element_blank(), 
        axis.text = element_text(colour = 'black'), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5), 
        axis.text.y = element_blank())
#ggsave('72hAEL_label_0.3_matchCelltype.pdf', units = 'cm', width = 8, height = 15)



library(Seurat)
library(plyr)
library(ggplot2)

lg72 <- read.delim('label_0.3.txt')
lg72$supergroup_new_subclustering <- mapvalues(lg72$supergroup, from = unique(lg72$supergroup), to = c('PSC', 'PH1-3', 'PH', 'PM', 'PM-NimC1'))
lg72$supergroup_new_subclustering <- factor(lg72$supergroup_new_subclustering, levels = c('PSC', 'PH1-3', 'PH', 'PM', 'PM-NimC1'))
lg72$supergroup_new_celltype <- mapvalues(lg72$supergroup_new_subclustering, from = levels(lg72$supergroup_new_subclustering),
                                          to = c('PSC', 'PH', 'PH', 'PM', 'PM'))
head(lg72)

lymphgland <- readRDS('../../Normal_LG/rdata/lymphgland.Rds')

lymphgland_72 <- subset(lymphgland, timepoint == 'AEL72hr')
lymphgland_72@meta.data <- droplevels(lymphgland_72@meta.data)
head(lymphgland_72@meta.data)

rownames(lg72) <- lg72$Barcode
lymphgland_72@meta.data$tp_cluster <- lymphgland_72@meta.data$anno_simple
lymphgland_72@meta.data[rownames(lg72), 'tp_cluster'] <- lg72$supergroup_new_celltype
lymphgland_72@meta.data <- droplevels(lymphgland_72@meta.data)
#writetable <- data.frame(Barcode = rownames(lymphgland_72@meta.data), lymphgland_72@meta.data, check.rows = F, check.names = F)
#write.table(writetable, 'label_0.3_72hAELmeta.txt', quote = F, sep = '\t', row.names = F, col.names = T)


plt <- DimPlot(lymphgland_72, reduction = 'tsne', group.by = 'tp_cluster') +
  scale_color_manual(values = c('#f179a8', '#1d91c5', '#b9060a', '#f5b053', '#975a1a', '#bc2aff'))
ggsave('Fig1D.tsne.tpcluster.pdf', units = 'cm', width = 11, height = 8)
AugmentPlot(plt + theme_void() + theme(legend.position = 'None'), width = 3, height = 3, dpi = 300)
ggsave('Fig1D.tsne.tpcluster.augment.pdf', units = 'cm', width = 3, height = 3)


