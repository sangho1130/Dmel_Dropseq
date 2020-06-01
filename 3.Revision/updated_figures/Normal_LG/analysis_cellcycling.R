library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)

lymphgland <- readRDS('rdata/lymphgland.Rds')
Idents(lymphgland) <- 'anno_simple'
levels(Idents(lymphgland))
lymphgland <- subset(lymphgland, idents = c('PSC', 'PH', 'PM', 'LM', 'CC', 'GST-rich', 'Adipohemocyte'))
Idents(lymphgland) <- 'new_subclustering'
summary(lymphgland@meta.data$new_subclustering)

  
cc_cells_orig <- data.frame(t(as.matrix(GetAssayData(lymphgland, slot = 'data'))), check.rows = F, check.names = F)
cc_cells_orig <- cc_cells_orig[, c('Cdk1', 'stg', 'CycB', 'polo', 'aurB', 'Det')]
head(cc_cells_orig)

# cell cycle stage score
cc_cells_orig$G1 <- cc_cells_orig$Cdk1
cc_cells_orig$G2 <- rowMeans(cc_cells_orig[, c('stg', 'CycB')])
cc_cells_orig$M <- rowMeans(cc_cells_orig[, c('polo', 'aurB', 'Det')])
cc_cells_orig <- cc_cells_orig[,-c(1:6)]
head(cc_cells_orig)

# filter low cell cycling cells part-1
cc_cells_orig$expsum <- rowSums(cc_cells_orig)
nrow(cc_cells_orig)
summary(cc_cells_orig$expsum)
cc_cells <- subset(cc_cells_orig, expsum >= 1)
nrow(cc_cells)
head(cc_cells)
cc_cells$expsum <- NULL

cc_cells$new_subclustering <- lymphgland@meta.data[rownames(cc_cells), 'new_subclustering']
cc_cells <- cc_cells[, c(4,1:3)]
head(cc_cells)

# filter low cell cycling cell part-2
summary(lymphgland@meta.data$new_subclustering)
summary(cc_cells$new_subclustering)
summary(cc_cells$new_subclustering)/summary(lymphgland@meta.data$new_subclustering)*100 > 33
median(na.omit(summary(cc_cells$new_subclustering)/summary(lymphgland@meta.data$new_subclustering)*100))

check_cc_cells <- mapply(function(X) { if (X[1] > 33) return(TRUE) else return(FALSE)}, na.omit(summary(cc_cells$new_subclustering)/summary(lymphgland@meta.data$new_subclustering)*100))
cc_cell_subclusters <- names(check_cc_cells[check_cc_cells == TRUE])
cc_cells <- subset(cc_cells, new_subclustering %in% cc_cell_subclusters)
cc_cells <- droplevels(cc_cells)
unique(cc_cells$new_subclustering)

# summary cell cycle scores by subclusters
head(cc_cells)
cc_cell_median <- data.frame(matrix(nrow = length(unique(cc_cells$new_subclustering)), ncol = 5)) ###
colnames(cc_cell_median) <- c(colnames(cc_cells), 'Prop')
cc_cell_median <- cc_cell_median[, c(1,5, 2:4)]
cc_cell_median$new_subclustering <- unique(cc_cells$new_subclustering)
rownames(cc_cell_median) <- cc_cell_median$new_subclustering
cc_cell_median$Prop <- as.character(summary(cc_cells$new_subclustering)/summary(droplevels(subset(lymphgland@meta.data, new_subclustering %in% unique(cc_cells$new_subclustering)))$new_subclustering) * 100)

for (subclust in unique(cc_cells$new_subclustering)) {
  tmp <- subset(cc_cells, new_subclustering == subclust)
  for (stage in c('G1', 'G2', 'M')) {
    cc_cell_median[subclust, stage] <- mean(tmp[, stage])
  }
}

cc_cell_median_m <- melt(cc_cell_median)
cc_cell_median_m$Prop <- as.numeric(cc_cell_median_m$Prop)
colnames(cc_cell_median_m) <- c('Subcluster', 'Proportion', 'Stage', 'Expression')
head(cc_cell_median_m)

plt <- ggplot(cc_cell_median_m, aes(Subcluster, Stage, size = Expression, colour = Proportion)) + 
  geom_point(alpha = .8) + 
  coord_flip() + 
  scale_color_gradient2(low = 'grey90', high = 'dodgerblue2') + 
  scale_size_continuous(range = c(1, 6)) + 
  theme_bw(base_size = 9) + 
  theme(text = element_text(colour = 'black', size = 9),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        panel.grid = element_blank()) + 
  labs(title = 'Cell cycle score', x = '', y = ''); plt
 

