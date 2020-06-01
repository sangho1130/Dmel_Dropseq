origin <- read.delim('origin/blood.combined_flt_MAST_degs.origin.txt')
origin <- subset(origin, p_val_adj <= 0.05 & pct.2 <= 0.5)
head(origin)

ph <- read.delim('ph/blood.combined_flt_MAST_degs.PH.txt')
ph <- subset(ph, p_val_adj <= 0.05 & pct.2 <= 0.5)
head(ph)

ph1 <- read.delim('ph1/blood.combined_MAST_degs.PH_1.txt')
ph1 <- subset(ph1, p_val_adj <= 0.05 & pct.2 <= 0.5)
head(ph1)

ph4 <- read.delim('ph4/blood.combined_flt_MAST_degs.PH4.txt')
ph4 <- subset(ph4, p_val_adj <= 0.05 & pct.2 <= 0.5)
head(ph4)

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
    intersect(
      intersect(subset(origin, cluster == 'LG')$gene, 
                ph$gene),
      ph4$gene),
    pm$gene); origin_lg


### Origin - Circ ###
origin_circ <-
  intersect(
    intersect(
      intersect(subset(origin, cluster == 'Circ')$gene, 
                ph$gene),
      ph4$gene),
    pm$gene); origin_circ


### PH 1-LG ###
ph1_lg <- 
  setdiff(
    setdiff(
      setdiff(
        setdiff(
          setdiff(subset(ph1, cluster == 'LG')$gene, 
                  subset(ph, cluster == 'LG')$gene),
          subset(ph4, cluster == 'LG')$gene),
        pm$gene),
      cc1$gene),
    cc2$gene); ph1_lg
### PH 1-Circ ###
ph1_circ <- 
  setdiff(
    setdiff(
      setdiff(
        setdiff(
          setdiff(subset(ph1, cluster == 'Circ')$gene, 
                  subset(ph, cluster == 'Circ')$gene),
          subset(ph4, cluster == 'Circ')$gene),
        pm$gene),
      cc1$gene),
    cc2$gene); ph1_circ


### PH 4-LG ###
ph4_lg <- 
  setdiff(
    setdiff(
      setdiff(
        setdiff(
          setdiff(subset(ph4, cluster == 'LG')$gene, 
                  subset(ph, cluster == 'LG')$gene),
          subset(ph1, cluster == 'LG')$gene),
        pm$gene),
      cc1$gene),
    cc2$gene); ph4_lg
### PH 4-Circ ###
ph4_circ <- 
  setdiff(
    setdiff(
      setdiff(
        setdiff(
          setdiff(subset(ph4, cluster == 'Circ')$gene, 
                  subset(ph, cluster == 'Circ')$gene),
          subset(ph1, cluster == 'Circ')$gene),
        pm$gene),
      cc1$gene),
    cc2$gene); ph4_circ


### PM-LG ###
pm_lg <- 
  setdiff(
    setdiff(
      setdiff(
        setdiff(
          setdiff(subset(pm, cluster == 'LG')$gene, ph$gene),
          ph4$gene), 
        ph1$gene),
      cc1$gene),
    cc2$gene); pm_lg
### PM-Circ ###
pm_circ <-
  setdiff(
    setdiff(
      setdiff(
        setdiff(
          setdiff(subset(pm, cluster == 'Circ')$gene, ph$gene),
          ph4$gene),
        ph1$gene),
      cc1$gene),
    cc2$gene); pm_circ


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
        setdiff(
          setdiff(subset(cc1, cluster == 'LG')$gene, ph$gene),
          ph4$gene),
        ph1$gene),
      pm$gene),
    cc2$gene); cc1_lg
### CC 1-Circ ###
cc1_circ <- 
  setdiff(
    setdiff(
      setdiff(
        setdiff(
          setdiff(subset(cc1, cluster == 'Circ')$gene, ph$gene),
          ph4$gene),
        ph1$gene),
      pm$gene),
    cc2$gene); cc1_circ


### CC 2-LG ###
cc2_lg <- 
  setdiff(
    setdiff(
      setdiff(
        setdiff(
          setdiff(subset(cc2, cluster == 'LG')$gene, ph$gene),
          ph4$gene),
        ph1$gene),
      pm$gene),
    cc1$gene); cc2_lg

### CC 2-Circ ###
cc2_circ <- 
  setdiff(
    setdiff(
      setdiff(
        setdiff(
          setdiff(subset(cc2, cluster == 'Circ')$gene, ph$gene),
          ph4$gene),
        ph1$gene),
      pm$gene),
    cc1$gene); cc2_circ


ph_lg
ph1_lg
ph1_circ
ph4_lg
ph4_circ
pm_lg
pm_circ
cc_circ
cc1_lg
cc1_circ
cc2_lg
cc2_circ
unique(c(ph_lg, ph1_lg, ph1_circ, pm_lg, pm_circ, cc_circ, cc1_lg, cc1_circ, cc2_lg, cc2_circ))
unique(c(ph1_lg, ph1_circ, ph4_lg, ph4_circ, pm_lg, pm_circ, cc_circ, cc1_lg, cc1_circ, cc2_lg, cc2_circ))


library(Seurat)
library(ggplot2)
library(plyr)
blood.combined_flt <- readRDS('../tmp/blood.combined_newsubclustering_flt.Rds')
DefaultAssay(blood.combined_flt)
head(blood.combined_flt@meta.data)

blood.combined_flt@meta.data$anno_simple_v2 <- mapvalues(blood.combined_flt@meta.data$Subclustering, 
                                                         from = levels(blood.combined_flt@meta.data$Subclustering),
                                                         to = c('PSC', 'PH 1', 'PH', 'PH', 'PH 4', 'PH', 'PH',
                                                                'PM', 'PM', 'PM', 'PM',
                                                                'LM', 'LM', 'CC 1', 'CC 2', 'GST-rich', 'Adipohemocyte'))
blood.combined_flt@meta.data$anno_simple_v2 <- factor(blood.combined_flt@meta.data$anno_simple_v2, 
                                                      levels = rev(c("PSC", "PH", "PH 1", "PH 4", "PM", "LM", "CC 1", "CC 2", "GST-rich", "Adipohemocyte")))
blood.combined_flt@meta.data$origin <- factor(blood.combined_flt@meta.data$origin, levels = rev(levels(blood.combined_flt@meta.data$origin)))

DotPlot(blood.combined_flt, group.by = 'anno_simple_v2', split.by = 'origin', 
        features = unique(rev(
          c(origin_lg, origin_circ,
            ph1_lg[1:5], ph1_circ[1:5],
            ph4_lg, ph4_circ,
            pm_lg, pm_circ[1:5],
            cc1_lg, cc1_circ,
            cc2_lg, cc2_circ[1:5], '14-3-3zeta')
        )), 
        dot.scale = 3, cols = c('#fd9409', '#69b9c2')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('intersect_genes_20200504_v2.pdf', units = 'cm', width = 25, height = 12)


Idents(blood.combined_flt) <- 'origin'
DotPlot(subset(blood.combined_flt, idents = 'LG'), 
        group.by = 'anno_simple_v2', 
        features = unique(rev(
          c(origin_lg, origin_circ,
            ph1_lg[1:5], ph1_circ[1:5],
            ph4_lg, ph4_circ,
            pm_lg, pm_circ[1:5],
            cc1_lg, cc1_circ,
            cc2_lg, cc2_circ[1:5], '14-3-3zeta')
        )), 
        dot.scale = 3, cols = c('grey90', '#fd9409')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('intersect_genes_LG_20200504_v2.pdf', units = 'cm', width = 25, height = 12)

DotPlot(subset(blood.combined_flt, idents = 'Circ'), 
        group.by = 'anno_simple_v2', 
        features = unique(rev(
          c(origin_lg, origin_circ,
            ph1_lg[1:5], ph1_circ[1:5],
            ph4_lg, ph4_circ,
            pm_lg, pm_circ[1:5],
            cc1_lg, cc1_circ,
            cc2_lg, cc2_circ[1:5], '14-3-3zeta')
        )), 
        dot.scale = 3, cols = c('grey90', '#69b9c2')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
#ggsave('intersect_genes_Circ_20200504_v2.pdf', units = 'cm', width = 25, height = 12)










