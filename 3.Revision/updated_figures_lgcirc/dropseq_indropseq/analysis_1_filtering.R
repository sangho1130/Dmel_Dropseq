library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

blood.combined <- readRDS('tmp/blood.combined_newsubclustering.Rds')
head(blood.combined@meta.data)

levels(blood.combined@meta.data$timepoint)
circ_only <- subset(blood.combined, cells = rownames(subset(blood.combined@meta.data, timepoint == "Circ_drop_96" | timepoint == "Circ_drop_120" | timepoint == "Circ_indrop_120")))
circ_only@meta.data <- droplevels(circ_only@meta.data); nrow(circ_only@meta.data)

countdf <- data.frame(row.names = unique(circ_only@meta.data$Library))
for (subcluster in levels(circ_only@meta.data$Subclustering)) {
  tmpdf <- subset(circ_only@meta.data, Subclustering == subcluster)
  tmpcount <- data.frame(count = summary(as.factor(tmpdf$Library)))
  tmpcount <- data.frame(row.names = unique(circ_only@meta.data$Library), count = tmpcount[unique(circ_only@meta.data$Library), ])
  countdf <- cbind(countdf, tmpcount)
}
colnames(countdf) <- levels(circ_only@meta.data$Subclustering)
countdf[is.na(countdf)] <- 0
head(countdf)

propdf <- countdf
for (subcluster in colnames(propdf)){
  propdf[,subcluster] <- propdf[,subcluster]/sum(propdf[,subcluster])*100
}

countdf <- data.frame(Library = rownames(countdf), countdf, check.rows = F, check.names = F)
write.table(countdf, 'stats/subcluster_count_unfiltered.txt', quote = F, sep = '\t', row.names = F, col.names = T)
propdf <- data.frame(Library = rownames(propdf), propdf, check.rows = F, check.names = F)
write.table(propdf, 'stats/subcluster_proportion_unfiltered.txt', quote = F, sep = '\t', row.names = F, col.names = T)


### Filtering library-biased and low count cell subclusters ###
head(blood.combined@meta.data)
blood.combined@meta.data$origin <- mapvalues(blood.combined@meta.data$timepoint, 
                                             from = levels(blood.combined@meta.data$timepoint),
                                             to = c('LG', 'LG', 'Circ', 'Circ', 'Circ'))
summary(subset(blood.combined@meta.data, origin == 'LG' & timepoint == 'LG_drop_96')$Subclustering)
summary(subset(blood.combined@meta.data, origin == 'LG' & timepoint == 'LG_drop_120')$Subclustering)


filtercells <- rownames(subset(subset(blood.combined@meta.data, origin == 'Circ'),
                               Subclustering == 'PSC' | Subclustering == 'PH 2' | Subclustering == 'PH 3' | Subclustering == 'PH 6' |
                                 Subclustering == 'PM 2' | Subclustering == 'PM 3' | Subclustering == 'PM 4' | 
                                 Subclustering == 'Adipohemocyte' | Subclustering == 'GST-rich' | Subclustering == 'LM 2'))
usecells <- setdiff(rownames(blood.combined@meta.data), filtercells)
blood.combined_flt <- subset(blood.combined, cells = usecells)
Idents(blood.combined_flt) <- 'anno_simple'
#DimPlot(blood.combined_flt, dims = c(1,2), split.by = 'origin')

#blood.combined_flt <- readRDS('tmp/blood.combined_newsubclustering_flt.Rds')
blood.combined_flt@meta.data$anno_simple <- mapvalues(blood.combined_flt@meta.data$Subclustering, 
                                                      from = levels(blood.combined_flt@meta.data$Subclustering),
                                                      to = c('PSC', 'PH', 'PH', 'PH', 'PH', 'PH', 'PH', 
                                                             'PM', 'PM', 'PM', 'PM', 'LM', 'LM', 'CC', 'CC', 'GST-rich', 'Adipohemocyte'))

writeLable <- data.frame(Barcode = rownames(blood.combined_flt@meta.data), blood.combined_flt@meta.data, check.rows = F, check.names = F)
write.table(writeLable, 'label_filtered.txt', quote = F, sep = '\t', row.names = F, col.names = T)

#saveRDS(blood.combined_flt, 'tmp/blood.combined_newsubclustering_flt.Rds')
