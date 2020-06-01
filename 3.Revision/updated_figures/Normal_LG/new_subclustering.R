library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)

#dir.create('rdata')
lymphgland <- readRDS('../../../../Project1_LymphGland/Drop-seq_alignment/1.seurat3_alignment_withMuscle_regress-Library-nUMI/lymphgland.combined.flt_allscaled.Rds')
ph_newlabel <- readRDS('../../PH/total ex_scbkgenes/tmp/ph.newlabel.Rds')
pm_newlabel <- readRDS('../../PM/total ex_scbkgenes/tmp/pm.newlabel.Rds')

head(lymphgland@meta.data)
lymphgland@meta.data$new_subclustering <- as.character(lymphgland@meta.data$Subclustering)

lymphgland@meta.data[rownames(ph_newlabel), 'new_subclustering'] <- as.character(ph_newlabel$supergroup)
lymphgland@meta.data[rownames(pm_newlabel), 'new_subclustering'] <- as.character(pm_newlabel$supergroup)

unique(lymphgland@meta.data$new_subclustering)
levels(lymphgland@meta.data$Subclustering)

lymphgland@meta.data$new_subclustering <- factor(lymphgland@meta.data$new_subclustering, 
                                                 levels = c("PSC", "PH 1", "PH 2", "PH 3", "PH 4", "PH 5", "PH 6",
                                                            "PM 1", "PM 2", "PM 3", "PM 4", 
                                                            "LM 1", "LM 2", "CC 1", "CC 2", 
                                                            "GST-rich", "Adipohemocyte", "DV", "RG", "Neurons"))
levels(lymphgland@meta.data$new_subclustering)
Idents(lymphgland) <- 'new_subclustering'
#saveRDS(lymphgland, 'rdata/lymphgland.Rds')
#saveRDS(lymphgland@meta.data, 'rdata/label.Rds')

