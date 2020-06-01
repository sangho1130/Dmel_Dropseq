library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)

lymphgland <- readRDS('tmp/lymphgland.combined_transferred_filtered.Rds')
head(lymphgland@meta.data)
Idents(lymphgland) <- 'labelTransfer'
lymphgland@meta.data <- droplevels(lymphgland@meta.data)
summary(lymphgland@meta.data$labelTransfer)

lymphgland@meta.data$labelTransfer_2 <- paste(lymphgland@meta.data$labelTransfer, lymphgland@meta.data$timepoint, sep='_')
Idents(lymphgland) <- 'labelTransfer_2'
lymphgland@meta.data$labelTransfer_2 <- factor(lymphgland@meta.data$labelTransfer_2,
                                               levels = c("PSC_Normal", "PSC_Infested", "PH 1_Normal", "PH 1_Infested", "PH 2_Normal", "PH 2_Infested",
                                                          "PH 3_Normal", "PH 3_Infested", "PH 4_Normal", "PH 4_Infested", "PM 1_Normal", "PM 1_Infested",
                                                          "LM 1_Normal", "LM 1_Infested", "LM 2_Normal", "LM 2_Infested", "CC 1_Normal", "CC 1_Infested", "CC 2_Normal",
                                                          "GST-rich_Normal", "GST-rich_Infested"))


cc_cells_orig <- data.frame(t(as.matrix(GetAssayData(lymphgland, slot = 'data', assay = 'RNA'))), check.rows = F, check.names = F)
cc_cells_orig <- cc_cells_orig[, c('Cdk1', 'CycD', 'CycE', 'stg', 'CycB', 'CycA', 'polo', 'aurB', 'Det')]
head(cc_cells_orig)

# cell cycle stage score
cc_cells_orig$G1 <- rowMeans(cc_cells_orig[, c('Cdk1', 'CycD', 'CycE')])
cc_cells_orig$G2 <- rowMeans(cc_cells_orig[, c('stg', 'CycB', 'CycA')])
cc_cells_orig$M <- rowMeans(cc_cells_orig[, c('polo', 'aurB', 'Det')])
cc_cells_orig <- cc_cells_orig[,c('G1', 'G2', 'M')]
head(cc_cells_orig)

# filter low cell cycling cells part-1
cc_cells_orig$expsum <- rowSums(cc_cells_orig)
nrow(cc_cells_orig)
summary(cc_cells_orig$expsum)
cc_cells <- subset(cc_cells_orig, expsum >= 0.9952)
nrow(cc_cells)
head(cc_cells)
cc_cells$expsum <- NULL

cc_cells$labelTransfer_2 <- lymphgland@meta.data[rownames(cc_cells), 'labelTransfer_2']
cc_cells <- cc_cells[, c(4,1:3)]
head(cc_cells)

# filter low cell cycling cell part-2
summary(as.factor(lymphgland@meta.data$labelTransfer_2))
summary(as.factor(cc_cells$labelTransfer_2))
summary(as.factor(cc_cells$labelTransfer_2))/summary(as.factor(lymphgland@meta.data$labelTransfer_2))[names(summary(as.factor(cc_cells$labelTransfer_2)))]*100 #> 25
cc_cells <- droplevels(cc_cells)

# summary cell cycle scores by subclusters
head(cc_cells)
cc_cell_median <- data.frame(matrix(nrow = length(unique(cc_cells$labelTransfer_2)), ncol = 5)) ###
colnames(cc_cell_median) <- c(colnames(cc_cells), 'Prop')
cc_cell_median <- cc_cell_median[, c(1,5, 2:4)]
cc_cell_median$labelTransfer_2 <- levels(cc_cells$labelTransfer_2)
rownames(cc_cell_median) <- cc_cell_median$labelTransfer_2
calc_prop <- summary(as.factor(cc_cells$labelTransfer_2))/summary(as.factor(lymphgland@meta.data$labelTransfer_2))[names(summary(as.factor(cc_cells$labelTransfer_2)))]*100
calc_prop[is.nan(calc_prop)] <- 0
cc_cell_median$Prop <- calc_prop


for (subclust in levels(cc_cells$labelTransfer_2)) {
  tmp <- subset(cc_cells, labelTransfer_2 == subclust)
  for (stage in c('G1', 'G2', 'M')) {
    cc_cell_median[subclust, stage] <- median(tmp[, stage])
  }
}
cc_cell_median
cc_cell_median[cc_cell_median$Prop <= 25, ][, c(2)] <- 10
cc_cell_median$Prop <- as.character(cc_cell_median$Prop)
cc_cell_median$G1 <- scale(cc_cell_median$G1)###
cc_cell_median$G2 <- scale(cc_cell_median$G2)###
cc_cell_median$M <- scale(cc_cell_median$M)###

cc_cell_median_m <- melt(cc_cell_median)
cc_cell_median_m$Prop <- as.numeric(cc_cell_median_m$Prop)
colnames(cc_cell_median_m) <- c('Subcluster', 'Proportion', 'Stage', 'Expression')
head(cc_cell_median_m)
cc_cell_median_m$Subcluster <- factor(cc_cell_median_m$Subcluster, levels = rev(levels(cc_cells$labelTransfer_2)))

plt <- ggplot(cc_cell_median_m, aes(Subcluster, Stage, size = Proportion, colour = Expression)) + 
  geom_point(alpha = .8) + 
  coord_flip() + 
  scale_color_gradient2(low = 'grey90', mid = 'royalblue', high = 'red2', midpoint = median(cc_cell_median_m$Expression)) + 
  scale_size_continuous(range = c(1, 5), breaks = c(10, 25, 50, 75)) + 
  theme_bw(base_size = 9) + 
  theme(text = element_text(colour = 'black', size = 9),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        panel.grid = element_blank()) + 
  labs(title = 'Cell cycle score', x = '', y = ''); plt
ggsave('analysis_cellcycling_v2_scaled.pdf', units = 'cm', width = 7, height = 12)

