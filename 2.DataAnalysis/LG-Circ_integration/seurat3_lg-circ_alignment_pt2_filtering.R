# PSC = C-120-2
# PH 2 = C-96-3
# PH 3 = C-96-1
# PH 4 = C-120-1 (24/45 = 53.33%) 
# PH 5 = C-120-3 (2/3 = 66.67%)
# PH 6 = C-120-1 (5/6 = 83.33%)
# PH 11 = C-96-3 (2/4 = 50.00%)
# PM 2 = C-96-3 (18/30 = 60.00%)
# PM 3 = C-120-2
# PM 9 = C-120-2 (31/38 = 81.58%)
# PM 10 = C-120-2
# Adipohemocyte = C-120-2 (8/15 = 53.33%)
# GST-rich = uninjured_2 (37/55 = 67.27%)
# LM 2 = C-120-2 (9/16 = 56.25%)

library(Seurat)
library(ggplot2)
library(ggrepel)
library(Matrix)
library(cowplot)
library(plyr)
library(dplyr)

blood.combined <- readRDS('blood.combined.Rds')
head(blood.combined@meta.data)
notusecells <- rownames(subset(subset(blood.combined@meta.data, origin == 'Circ'),
                            Subclustering == 'PSC' | Subclustering == 'PH 2' | Subclustering == 'PH 3' | Subclustering == 'PH 4' |
                              Subclustering == 'PH 5' | Subclustering == 'PH 6' | Subclustering == 'PH 11' | Subclustering == 'PM 2' |
                              Subclustering == 'PM 3' | Subclustering == 'PM 9' | Subclustering == 'PM 10' | Subclustering == 'Adipohemocyte' | 
                              Subclustering == 'GST-rich' | Subclustering == 'LM 2'))
usecells <- setdiff(rownames(blood.combined@meta.data), notusecells)
blood.combined_flt <- subset(blood.combined, cells = usecells)
blood.combined_flt@meta.data <- droplevels(blood.combined_flt@meta.data)

summary(as.factor(subset(blood.combined@meta.data, origin == 'LG')$Subclustering))
summary(as.factor(subset(blood.combined_flt@meta.data, origin == 'LG')$Subclustering))
summary(as.factor(subset(blood.combined_flt@meta.data, origin == 'Circ')$Subclustering))
nrow(subset(blood.combined_flt@meta.data, origin == 'Circ'))*0.01 # = 37.18

# PH 7 (5); PM 1 (5), PM 8 (2)
filterCells <- rownames(subset(subset(blood.combined@meta.data, origin == 'Circ'),
                               Subclustering == 'PH 7' | Subclustering == 'PM 1' | Subclustering == 'PM 8'))
usecells_2 <- setdiff(rownames(blood.combined_flt@meta.data), filterCells)
blood.combined_flt <- subset(blood.combined_flt, cells = usecells_2)
blood.combined_flt@meta.data <- droplevels(blood.combined_flt@meta.data)

write.table(data.frame(Barcode = rownames(blood.combined_flt@meta.data), blood.combined_flt@meta.data, check.rows = F, check.names = F),
            'tmp7_flt.label.txt', col.names = T, row.names = F, sep = '\t', quote = F)
DefaultAssay(blood.combined_flt)
save.image('tmp7_flt.Rdata')
saveRDS(blood.combined_flt, 'blood.combined_flt.Rds')

count_data <- as.matrix(GetAssayData(blood.combined_flt, slot = 'counts'))
count_data <- data.frame(Symbol = rownames(count_data), count_data, check.rows = F, check.names = F)
#write.table(count_data, 'blood.combined_flt.count_data.txt', quote = F, sep = '\t', row.names = F, col.names = T)

