library(Seurat)
library(plyr)

lymphgland <- readRDS('tmp/lymphglandt.fix.Rds')
newlabel <- readRDS('../Normal_LG/rdata/label.Rds')

lymphgland@meta.data$new_subclustering <- newlabel[rownames(lymphgland@meta.data), ]$new_subclustering
head(lymphgland@meta.data)
writeTable <- data.frame(Barcode = rownames(lymphgland@meta.data), lymphgland@meta.data, check.rows = F, check.names = F)
write.table(writeTable, 'res_fix_0.7/label_0.7.txt', quote = F, sep = '\t', row.names = F, col.names = T)
