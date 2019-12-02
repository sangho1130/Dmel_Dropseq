library(Seurat)
library(ggplot2)
library(ggrepel)

tpsmeans <- read.delim('bulk_sc_pseudo.tpm_vs_count.txt', sep = '\t', row.names = 1)
head(tpsmeans)
nrow(tpsmeans)

# 72
tpsmeans_72 <- read.delim('bulk_sc72_pseudo.tpm_vs_count.txt', sep = '\t', row.names = 1, check.names = F)
head(tpsmeans_72)
# 96
tpsmeans_96 <- read.delim('bulk_sc96_pseudo.tpm_vs_count.txt', sep = '\t', row.names = 1, check.names = F)
head(tpsmeans_96)
# 120
tpsmeans_120 <- read.delim('bulk_sc120_pseudo.tpm_vs_count.txt', sep = '\t', row.names = 1, check.names = F)
head(tpsmeans_120)

nrow(tpsmeans_72)
nrow(tpsmeans_96)
nrow(tpsmeans_120)


tpsmeans <- data.frame(bulk_72 = numeric(nrow(tpsmeans)),
                       bulk_96 = numeric(nrow(tpsmeans)),
                       bulk_120 = numeric(nrow(tpsmeans)),
                       pseudo_72 = numeric(nrow(tpsmeans)),
                       pseudo_96 = numeric(nrow(tpsmeans)),
                       pseudo_120 = numeric(nrow(tpsmeans)),
                       row.names = rownames(tpsmeans),
                       check.rows = F, 
                       check.names = F)

tpsmeans[rownames(tpsmeans_72), 'bulk_72'] <- tpsmeans_72$`72AEL`
tpsmeans[rownames(tpsmeans_96), 'bulk_96'] <- tpsmeans_96$`96AEL`
tpsmeans[rownames(tpsmeans_120), 'bulk_120'] <- tpsmeans_120$`120AEL`
tpsmeans[rownames(tpsmeans_72), 'pseudo_72'] <- tpsmeans_72$`pseudo`
tpsmeans[rownames(tpsmeans_96), 'pseudo_96'] <- tpsmeans_96$`pseudo`
tpsmeans[rownames(tpsmeans_120), 'pseudo_120'] <- tpsmeans_120$`pseudo`
head(tpsmeans)


tpsmeans$pseudo_72_norm <- tpsmeans$pseudo_72/(sum(tpsmeans$pseudo_72)/1000000)
tpsmeans$pseudo_96_norm <- tpsmeans$pseudo_96/(sum(tpsmeans$pseudo_96)/1000000)
tpsmeans$pseudo_120_norm <- tpsmeans$pseudo_120/(sum(tpsmeans$pseudo_120)/1000000)

tpsmeans$bulk_mean <- (tpsmeans$bulk_72 + tpsmeans$bulk_96 + tpsmeans$bulk_120)/3
tpsmeans$pseudo_mean <- (tpsmeans$pseudo_72 + tpsmeans$pseudo_96 + tpsmeans$pseudo_120)/3
tpsmeans$pseudo_mean_norm <- tpsmeans$pseudo_mean/(sum(tpsmeans$pseudo_mean)/1000000)
head(tpsmeans)


### GST genes ###
gstgenes <- c('CG15784', 'CR44430', 'CG18547', 'dys', 'CG3397', 'Dh31', 'AOX1', 'JhI-26', 'CG10638', 'CG2909')
tpsmeans$GST <- 'None'
tpsmeans[gstgenes, 'GST'] <- 'GST'
tpsmeans$GST <- factor(tpsmeans$GST, levels = c('None', 'GST'))

tpsmeans$genelabel <- NA
tpsmeans[gstgenes, 'genelabel'] <- gstgenes

ggplot(tpsmeans, aes(log10(bulk_mean+1), log10(pseudo_mean_norm+1), col = GST)) + 
  scale_color_manual(values = c('grey95', 'red2')) + 
  geom_point(size = .5, alpha = 0.5) + 
  geom_text_repel(aes(label = genelabel), nudge_x = 4, size = 2) + 
  theme_bw() + theme(axis.text = element_text(colour = 'black'), panel.grid = element_blank(), legend.position = 'None') + 
  xlim(0, 6) + ylim(0, 6)
#ggsave('compare_gstrich.pdf', units = 'cm', width = 8, height = 8)


###########
### vis ###
###########
tpsmeans$GST <- NA
tpsmeans[gstgenes, 'GST'] <- 'GST'
tpsmeans$GST <- factor(tpsmeans$GST, levels = c('None', 'GST'))

ggplot(tpsmeans, aes(log10(bulk_mean+1), log10(pseudo_mean_norm+1), col = GST)) + 
  scale_color_manual(values = c('red2')) + 
  geom_point(size = .5) + 
  theme_bw() + theme(axis.text = element_text(colour = 'black'), panel.grid = element_blank(), legend.position = 'None') + 
  xlim(0, 6) + ylim(0, 6)
#ggsave('compare_gstrich_gstrich_only.pdf', units = 'cm', width = 8, height = 8)

tpsmeans$GST <- 'None'
tpsmeans[gstgenes, 'GST'] <- NA
tpsmeans$GST <- factor(tpsmeans$GST, levels = c('None', 'GST'))

plt <- ggplot(tpsmeans, aes(log10(bulk_mean+1), log10(pseudo_mean_norm+1), col = GST)) + 
  scale_color_manual(values = c('grey90')) + 
  geom_point(size = 1) + 
  theme_void() + theme(legend.position = 'None') + 
  xlim(0, 6) + ylim(0, 6);plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('compare_gstrich_others_only.pdf', units = 'cm', width = 4, height = 4)



load('../../Drop-seq_alignment/1.seurat3_alignment_withMuscle_regress-Library-nUMI/__filtered__SCrm__labelModi__Seurat3.alignment.3tps.mitoCut.Rdata')
levels(lymphgland.combined.flt@meta.data$anno_simple)

DefaultAssay(lymphgland.combined.flt) <- 'RNA'
VlnPlot(lymphgland.combined.flt, group.by = 'anno_simple', pt.size = 0, ncol = 2,
        features = gstgenes)
#ggsave('compare_gstrich_vlnplot.pdf', units = 'cm', width = 20, height = 30)



### Adipohemocyte genes ###
adipogenes <- c('E23', 'CG11899', 'hid', 'Gal', 'Gs1', 'GstT4', 'Grp', 'Sirup', 'CG7860', 'CG9989')
tpsmeans$Adipo <- 'None'
tpsmeans[adipogenes, 'Adipo'] <- 'Adipo'
tpsmeans$Adipo <- factor(tpsmeans$Adipo, levels = c('None', 'Adipo'))

tpsmeans$genelabel <- NA
tpsmeans[adipogenes, 'genelabel'] <- adipogenes


ggplot(tpsmeans, aes(log10(bulk_mean+1), log10(pseudo_mean_norm+1), col = Adipo)) + 
  scale_color_manual(values = c('grey95', 'red2')) + 
  geom_point(size = .5, alpha = 0.5) + 
  geom_text_repel(aes(label = genelabel), nudge_x = 4, size = 2) + 
  theme_bw() + theme(axis.text = element_text(colour = 'black'), panel.grid = element_blank(), legend.position = 'None') + 
  xlim(0, 6) + ylim(0, 6)
#ggsave('adipohemocyte.pdf', units = 'cm', width = 8, height = 8)


###########
### vis ###
###########
tpsmeans$Adipo <- NA
tpsmeans[adipogenes, 'Adipo'] <- 'Adipo'
tpsmeans$Adipo <- factor(tpsmeans$Adipo, levels = c('None', 'Adipo'))

ggplot(tpsmeans, aes(log10(bulk_mean+1), log10(pseudo_mean_norm+1), col = Adipo)) + 
  scale_color_manual(values = c('red2')) + 
  geom_point(size = .5) + 
  theme_bw() + theme(axis.text = element_text(colour = 'black'), panel.grid = element_blank(), legend.position = 'None') + 
  xlim(0, 6) + ylim(0, 6)
#ggsave('adipohemocyte_gstrich_only.pdf', units = 'cm', width = 8, height = 8)

tpsmeans$Adipo <- 'None'
tpsmeans[adipogenes, 'Adipo'] <- NA
tpsmeans$Adipo <- factor(tpsmeans$Adipo, levels = c('None', 'Adipo'))

plt <- ggplot(tpsmeans, aes(log10(bulk_mean+1), log10(pseudo_mean_norm+1), col = Adipo)) + 
  scale_color_manual(values = c('grey90')) + 
  geom_point(size = 1) + 
  theme_void() + theme(legend.position = 'None') + 
  xlim(0, 6) + ylim(0, 6);plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
#ggsave('adipohemocyte_others_only.pdf', units = 'cm', width = 4, height = 4)




VlnPlot(lymphgland.combined.flt, group.by = 'anno_simple', pt.size = 0, ncol = 2,
        features = adipogenes)
#ggsave('compare_adipohemocyte_vlnplot.pdf', units = 'cm', width = 20, height = 30)
