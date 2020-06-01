origin <- read.delim('origin/blood.combined_flt_MAST_degs.origin.txt')
origin <- subset(origin, p_val_adj <= 0.05 & pct.2 <= 0.5)
head(origin)

ph <- read.delim('ph/blood.combined_flt_MAST_degs.PH.txt')
ph <- subset(ph, p_val_adj <= 0.05 & pct.2 <= 0.5)
head(ph)

ph1 <- read.delim('ph1/blood.combined_MAST_degs.PH_1.txt')
ph1 <- subset(ph1, p_val_adj <= 0.05 & pct.2 <= 0.5)
head(ph1)

pm <- read.delim('pm/blood.combined_flt_MAST_degs.PM.txt')
pm <- subset(pm, p_val_adj <= 0.05 & pct.2 <= 0.5)
head(pm)

cc <- read.delim('cc/blood.combined_flt_MAST_degs.CC.txt')
cc <- subset(cc, p_val_adj <= 0.05 & pct.2 <= 0.5)
head(cc)
cc1 <- read.delim('cc1/blood.combined_flt_MAST_degs.CC_1.txt')
cc1 <- subset(cc1, p_val_adj <= 0.05 & pct.2 <= 0.5)
head(cc1)
cc2 <- read.delim('cc2/blood.combined_flt_MAST_degs.CC_2.txt')
cc2 <- subset(cc2, p_val_adj <= 0.05 & pct.2 <= 0.5)
head(cc2)

### Origin - LG ###
origin_lg <-
  intersect(
    intersect(subset(origin, cluster == 'LG')$gene, 
              ph$gene), 
    pm$gene)
### Origin - Circ ###
origin_circ <-
  intersect(
    intersect(subset(origin, cluster == 'Circ')$gene, 
              ph$gene), 
    pm$gene)

### PH-LG ###
ph_lg <- 
  setdiff(
    setdiff(
      setdiff(
        setdiff(subset(ph, cluster == 'LG')$gene, pm$gene),
                ph1$gene),
      cc1$gene),
    cc2$gene)


### PH 1-LG ###
ph1_lg <- 
  setdiff(
    setdiff(
      setdiff(
        setdiff(subset(ph1, cluster == 'LG')$gene, pm$gene),
                ph$gene),
      cc1$gene),
    cc2$gene)
### PH 1-Circ ###
ph1_circ <- 
  setdiff(
    setdiff(
      setdiff(
        setdiff(subset(ph1, cluster == 'Circ')$gene, pm$gene),
                subset(ph, cluster == 'LG')$gene),
      cc1$gene),
    cc2$gene)


### PM-LG ###
pm_lg <- 
  setdiff(
    setdiff(
      setdiff(
        setdiff(subset(pm, cluster == 'LG')$gene, ph$gene),
                ph1$gene),
      cc1$gene),
    cc2$gene)
### PM-Circ ###
pm_circ <-
  setdiff(
    setdiff(
      setdiff(
        setdiff(subset(pm, cluster == 'Circ')$gene, ph$gene),
                ph1$gene),
      cc1$gene),
    cc2$gene)


### CC-Circ ###
cc_circ <- 
  setdiff(
    setdiff(
      setdiff(subset(cc, cluster == 'Circ')$gene, ph$gene),
      pm$gene),
    ph1$gene)


### CC 1-LG ###
cc1_lg <- 
  setdiff(
    setdiff(
      setdiff(
        setdiff(subset(cc1, cluster == 'LG')$gene, ph$gene),
        pm$gene),
      ph1$gene),
    cc2$gene)
### CC 1-Circ ###
cc1_circ <- 
  setdiff(
    setdiff(
      setdiff(
        setdiff(subset(cc1, cluster == 'Circ')$gene, ph$gene),
        pm$gene),
      ph1$gene),
    cc2$gene)


### CC 2-LG ###
cc2_lg <- 
  setdiff(
    setdiff(
      setdiff(
        setdiff(subset(cc2, cluster == 'LG')$gene, ph$gene),
        pm$gene),
      ph1$gene),
    cc1$gene)

### CC 2-Circ ###
cc2_circ <- 
  setdiff(
    setdiff(
      setdiff(
        setdiff(subset(cc2, cluster == 'Circ')$gene, ph$gene),
        pm$gene),
      ph1$gene),
    cc1$gene)


ph_lg
ph1_lg
ph1_circ
pm_lg
pm_circ
cc_circ
cc1_lg
cc1_circ
cc2_lg
cc2_circ
unique(c(ph_lg, ph1_lg, ph1_circ, pm_lg, pm_circ, cc_circ, cc1_lg, cc1_circ, cc2_lg, cc2_circ))


library(Seurat)
library(ggplot2)
library(plyr)
blood.combined_flt <- readRDS('../tmp/blood.combined_newsubclustering_flt.Rds')
DefaultAssay(blood.combined_flt)
head(blood.combined_flt@meta.data)

blood.combined_flt@meta.data$anno_simple_v2 <- mapvalues(blood.combined_flt@meta.data$Subclustering, 
                                                         from = levels(blood.combined_flt@meta.data$Subclustering),
                                                         to = c('PSC', 'PH 1', 'PH', 'PH', 'PH', 'PH', 'PH',
                                                                'PM', 'PM', 'PM', 'PM',
                                                                'LM', 'LM', 'CC 1', 'CC 2', 'GST-rich', 'Adipohemocyte'))
blood.combined_flt@meta.data$anno_simple_v2 <- factor(blood.combined_flt@meta.data$anno_simple_v2, levels = rev(levels(blood.combined_flt@meta.data$anno_simple_v2)))
blood.combined_flt@meta.data$origin <- factor(blood.combined_flt@meta.data$origin, levels = rev(levels(blood.combined_flt@meta.data$origin)))

DotPlot(blood.combined_flt, group.by = 'anno_simple_v2', split.by = 'origin', 
        features = unique(rev(
          c(origin_lg, origin_circ,
            ph1_lg, ph1_circ,
            pm_lg, pm_circ,
            cc1_lg, cc1_circ,
            cc2_lg, cc2_circ, 'Df31')
        )), 
        dot.scale = 3, cols = c('#fd9409', '#69b9c2')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('intersect_genes.pdf', units = 'cm', width = 25, height = 12)


Idents(blood.combined_flt) <- 'origin'
DotPlot(subset(blood.combined_flt, idents = 'LG'), 
        group.by = 'anno_simple_v2', 
        features = unique(rev(
          c(origin_lg, origin_circ,
            ph1_lg, ph1_circ,
            pm_lg, pm_circ,
            cc1_lg, cc1_circ,
            cc2_lg, cc2_circ, 'Df31')
        )), 
        dot.scale = 3, cols = c('grey90', '#fd9409')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('intersect_genes_LG.pdf', units = 'cm', width = 25, height = 12)

DotPlot(subset(blood.combined_flt, idents = 'Circ'), 
        group.by = 'anno_simple_v2', 
        features = unique(rev(
          c(origin_lg, origin_circ,
            ph1_lg, ph1_circ,
            pm_lg, pm_circ,
            cc1_lg, cc1_circ,
            cc2_lg, cc2_circ, 'Df31')
        )), 
        dot.scale = 3, cols = c('grey90', '#69b9c2')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('intersect_genes_Circ.pdf', units = 'cm', width = 25, height = 12)










