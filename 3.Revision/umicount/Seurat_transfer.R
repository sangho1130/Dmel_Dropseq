library(Seurat)
library(ggplot2)

### Objects ###
refObj <- readRDS('../Drop-seq_alignment/1.seurat3_alignment_withMuscle_regress-Library-nUMI/lymphgland.combined.flt_allscaled.Rds')
head(refObj@meta.data)
#DefaultAssay(refObj) <- 'RNA'

ggplot(refObj@meta.data, aes(anno_simple, nCount_RNA)) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, width = .1) +
  labs(x = '', y = 'UMI conut') +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave('celltypes_umicount.pdf', units = 'cm', width = 18, height = 10)
ggplot(refObj@meta.data, aes(anno_simple, nFeature_RNA)) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, width = .1) +
  labs(x = '', y = 'Gene conut') +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave('celltypes_genecount.pdf', units = 'cm', width = 18, height = 10)


