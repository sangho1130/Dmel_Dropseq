library(ggplot2)

bulk_sc_pseudo <- read.delim('bulk_sc_pseudo.tpm_vs_count.txt', sep = '\t', row.names = 1)
head(bulk_sc_pseudo)
nrow(bulk_sc_pseudo)

# 72
bulk_sc_pseudo_72 <- read.delim('bulk_sc72_pseudo.tpm_vs_count.txt', sep = '\t', row.names = 1, check.names = F)
head(bulk_sc_pseudo_72)
# 96
bulk_sc_pseudo_96 <- read.delim('bulk_sc96_pseudo.tpm_vs_count.txt', sep = '\t', row.names = 1, check.names = F)
head(bulk_sc_pseudo_96)
# 120
bulk_sc_pseudo_120 <- read.delim('bulk_sc120_pseudo.tpm_vs_count.txt', sep = '\t', row.names = 1, check.names = F)
head(bulk_sc_pseudo_120)

nrow(bulk_sc_pseudo_72)
nrow(bulk_sc_pseudo_96)
nrow(bulk_sc_pseudo_120)


tpsmeans <- data.frame(bulk_72 = numeric(nrow(bulk_sc_pseudo)),
                       bulk_96 = numeric(nrow(bulk_sc_pseudo)),
                       bulk_120 = numeric(nrow(bulk_sc_pseudo)),
                       pseudo_72 = numeric(nrow(bulk_sc_pseudo)),
                       pseudo_96 = numeric(nrow(bulk_sc_pseudo)),
                       pseudo_120 = numeric(nrow(bulk_sc_pseudo)),
                       row.names = rownames(bulk_sc_pseudo),
                       check.rows = F, 
                       check.names = F)

tpsmeans[rownames(bulk_sc_pseudo_72), 'bulk_72'] <- bulk_sc_pseudo_72$`72AEL`
tpsmeans[rownames(bulk_sc_pseudo_96), 'bulk_96'] <- bulk_sc_pseudo_96$`96AEL`
tpsmeans[rownames(bulk_sc_pseudo_120), 'bulk_120'] <- bulk_sc_pseudo_120$`120AEL`
tpsmeans[rownames(bulk_sc_pseudo_72), 'pseudo_72'] <- bulk_sc_pseudo_72$`pseudo`
tpsmeans[rownames(bulk_sc_pseudo_96), 'pseudo_96'] <- bulk_sc_pseudo_96$`pseudo`
tpsmeans[rownames(bulk_sc_pseudo_120), 'pseudo_120'] <- bulk_sc_pseudo_120$`pseudo`

head(tpsmeans)

tpsmeans$bulk_mean <- rowSums(tpsmeans[, c(1:3)])/3
tpsmeans$bulk_sum <- rowSums(tpsmeans[, c(1:3)])
tpsmeans$pseudo_mean <- rowSums(tpsmeans[, c(4:6)])/3
tpsmeans$pseudo_sum <- rowSums(tpsmeans[, c(4:6)])
head(tpsmeans)

scGenes <- rownames(subset(tpsmeans, bulk_mean == 0 | bulk_mean != 0 & pseudo_mean/bulk_mean >= 500))
#write.table(scGenes, 'bulk_sc_pseudo_pt3-scGenes.scGenes.txt', row.names = F, col.names = F, quote = F)
length(scGenes)

tpsmeans$sc_spc <- 'No'
tpsmeans[scGenes, 'sc_spc'] <- 'Yes'


summary(tpsmeans$bulk_mean)
summary(tpsmeans$pseudo_mean)

ggplot(tpsmeans, aes(log10(bulk_mean+1), log10(pseudo_mean+1), col = sc_spc)) + 
  geom_point(size = .5) + 
  scale_color_manual(values = c('grey90', 'red2')) +
  theme_bw() + theme(axis.text = element_text(colour = 'black'), panel.grid = element_blank()) + xlim(0, 7) + ylim(0, 7)

ggplot(tpsmeans, aes(log10(bulk_sum+1), log10(pseudo_sum+1), col = sc_spc)) + 
  geom_point(size = .5) + 
  scale_color_manual(values = c('grey90', 'red2')) +
  theme_bw() + 
  theme(axis.text = element_text(colour = 'black'), 
        panel.grid = element_blank(),
        legend.position = 'None') + 
  xlim(0, 7) + ylim(0, 7)
#ggsave('bulk_sc_pseudo_pt3-scGenes.scGenes.pdf', units = 'cm', width = 8, height = 8)


###
head(tpsmeans)
tpsmeans$pseudo_72_norm <- tpsmeans$pseudo_72/(sum(tpsmeans$pseudo_72)/1000000)
tpsmeans$pseudo_96_norm <- tpsmeans$pseudo_96/(sum(tpsmeans$pseudo_96)/1000000)
tpsmeans$pseudo_120_norm <- tpsmeans$pseudo_120/(sum(tpsmeans$pseudo_120)/1000000)
tpsmeans$pseudo_norm_mean <- rowSums(tpsmeans[, c(12:14)])/3
tpsmeans$pseudo_norm_sum <- rowSums(tpsmeans[, c(12:14)])

scGenes_v2 <- rownames(subset(tpsmeans, bulk_mean == 0 | bulk_mean != 0 & pseudo_norm_mean/bulk_mean >= 10))
#write.table(scGenes_v2, 'bulk_sc_pseudo_pt3-scGenes.scGenes_v2.txt', row.names = F, col.names = F, quote = F)
length(scGenes_v2)

tpsmeans$sc_spc_v2 <- 'No'
tpsmeans[scGenes_v2, 'sc_spc_v2'] <- 'Yes'


bkGenes_v2 <- rownames(subset(tpsmeans, pseudo_norm_mean == 0 | pseudo_norm_mean != 0 & bulk_mean/pseudo_norm_mean >= 10))
#write.table(bkGenes_v2, 'bulk_sc_pseudo_pt3-scGenes.bkGenes_v2.txt', row.names = F, col.names = F, quote = F)
length(bkGenes_v2)

tpsmeans$bk_spc_v2 <- 'No'
tpsmeans[bkGenes_v2, 'bk_spc_v2'] <- 'Yes'

tpsmeans$sc_bk_spc_v2 <- 'No'
tpsmeans[scGenes_v2, 'sc_bk_spc_v2'] <- 'Single-cell'
tpsmeans[bkGenes_v2, 'sc_bk_spc_v2'] <- 'Bulk'
tpsmeans$sc_bk_spc_v2 <- factor(tpsmeans$sc_bk_spc_v2, levels = c('Single-cell', 'Bulk', 'No'))

ggplot(tpsmeans, aes(log10(bulk_sum+1), log10(pseudo_norm_sum+1), col = sc_bk_spc_v2)) + 
  geom_point(size = .5) + 
  scale_color_manual(values = c('red2', 'steelblue2', 'grey90')) +
  theme_bw() + 
  theme(axis.text = element_text(colour = 'black'), 
        panel.grid = element_blank(),
        legend.position = 'None') + 
  xlim(0, 6) + ylim(0, 6)
#ggsave('bulk_sc_pseudo_pt3-scGenes.scGenes_v2.pdf', units = 'cm', width = 8, height = 8)
#ggsave('bulk_sc_pseudo_pt3-scGenes.sc_bk_spc_v2.pdf', units = 'cm', width = 8, height = 8)

plt <- ggplot(tpsmeans, aes(log10(bulk_sum+1), log10(pseudo_norm_sum+1), col = sc_bk_spc_v2)) + 
  geom_point(size = 1) + 
  scale_color_manual(values = c('red2', 'steelblue2', 'grey90')) +
  theme_void() + theme(legend.position = 'None') + xlim(0, 6) + ylim(0, 6)
Seurat::AugmentPlot(plt, dpi = 300, width = 4, height = 4) 
#ggsave('bulk_sc_pseudo_pt3-scGenes.sc_bk_spc_v2.augment.pdf', units = 'cm', width = 4, height = 4)
