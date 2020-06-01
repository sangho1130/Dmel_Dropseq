library(Seurat)
library(plyr)

lymphgland <- readRDS('tmp/lymphglandt.fix.res0.3.Rds')
newlabel <- readRDS('../Normal_LG/rdata/label.Rds')

lymphgland@meta.data$new_subclustering <- newlabel[rownames(lymphgland@meta.data), ]$new_subclustering
head(lymphgland@meta.data)
writeTable <- data.frame(Barcode = rownames(lymphgland@meta.data), lymphgland@meta.data, check.rows = F, check.names = F)
write.table(writeTable, 'res_fix_0.3//label_0.3.txt', quote = F, sep = '\t', row.names = F, col.names = T)


lymphgland <- readRDS('tmp/lymphglandt.fix.res0.6.Rds')
newlabel <- readRDS('../Normal_LG/rdata/label.Rds')

lymphgland@meta.data$new_subclustering <- newlabel[rownames(lymphgland@meta.data), ]$new_subclustering
head(lymphgland@meta.data)
writeTable <- data.frame(Barcode = rownames(lymphgland@meta.data), lymphgland@meta.data, check.rows = F, check.names = F)
write.table(writeTable, 'res_fix_0.6/label_0.6.txt', quote = F, sep = '\t', row.names = F, col.names = T)
