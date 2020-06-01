library(monocle3)
library(ggplot2)
library(plyr)
library(ggridges)
library(viridis)
library(tidyr)
options(scipen = 100)
load('../../../Drop-seq_alignment/2.monocle3/filter__PSC+DV+Neurons+PM11__20190829_v3/run_monocle3_major_v3.Rdata')
newlabel <- readRDS('../Normal_LG/rdata/label.Rds')
head(newlabel)

head(pData(cds))
pData(cds)$new_subclustering <- subset(newlabel, anno_simple != 'PSC' & anno_simple != 'DV' & anno_simple != 'Neurons' & anno_simple != 'RG')$new_subclustering
pData(cds)$anno_simple <- subset(newlabel, anno_simple != 'PSC' & anno_simple != 'DV' & anno_simple != 'Neurons' & anno_simple != 'RG')$anno_simple

cds_pseudotime <- data.frame(pseudotime = pseudotime(cds, reduction_method = 'UMAP'))
head(cds_pseudotime)
identical(rownames(cds_pseudotime), rownames(pData(cds)))

cds_pseudotime$timepoint <- pData(cds)$timepoint
cds_pseudotime$new_subclustering <- pData(cds)$new_subclustering
cds_pseudotime$anno_simple <- pData(cds)$anno_simple
cds_pseudotime <- subset(cds_pseudotime, pseudotime != Inf)
cds_pseudotime <- droplevels(cds_pseudotime)
head(cds_pseudotime)


ggplot(cds_pseudotime, aes(x = pseudotime, y = new_subclustering, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 2.5, rel_min_height = 0.01, gradient_lwd = 1, bandwidth = 2) +
  scale_x_continuous(expand = c(0.001, 0)) +
  scale_y_discrete(expand = c(0.001, 0)) +
  scale_fill_viridis(name = "Pseudotime", option = "C") +
  labs(title = 'Monocle 3', x = 'Pseudotime', y = '') + 
  theme_ridges(font_size = 9, grid = TRUE) + 
  theme(axis.title.y = element_blank(),
        axis.text = element_text(colour = 'black'),
        axis.title.x = element_text(hjust = .5),
        plot.title = element_text(hjust = .5))
#ggsave('density/densityplot.new_subclustering.pdf', units = 'cm', width = 12, height = 10)
#

head(pData(cds))
cds_pseudotime$anno_simple_v2 <- mapvalues(cds_pseudotime$new_subclustering, 
                                           from = levels(cds_pseudotime$new_subclustering),
                                           to = c('PH 1', 'PH 2', 'PH', 'PH', 'PH', 'PH',
                                                  'PM', 'PM', 'PM', 'PM', 
                                                  'LM', 'LM', 'CC', 'CC', 'GST-rich', 'Adipohemocyte'))

ggplot(cds_pseudotime, aes(x = pseudotime, y = anno_simple_v2, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 2.5, rel_min_height = 0.01, gradient_lwd = 1, bandwidth = 1.5) +
  scale_x_continuous(expand = c(0.001, 0)) +
  scale_y_discrete(expand = c(0.001, 0)) +
  scale_fill_viridis(name = "Pseudotime", option = "C") +
  labs(title = 'Monocle 3', x = 'Pseudotime', y = '') + 
  theme_ridges(font_size = 9, grid = TRUE) + 
  theme(axis.title.y = element_blank(), axis.text = element_text(colour = 'black'), axis.title.x = element_text(hjust = .5), plot.title = element_text(hjust = .5))
#ggsave('density/densityplot.anno_simple_v2.pdf', units = 'cm', width = 10, height = 7)

