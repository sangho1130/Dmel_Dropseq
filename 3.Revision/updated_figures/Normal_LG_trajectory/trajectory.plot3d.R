library(monocle3)
library(Matrix)
library(ggplot2)
library(Seurat)
library(pheatmap)
library(rgl)
library(plyr)
library(dplyr)
options(scipen = 100)
load('../../../Drop-seq_alignment/2.monocle3/filter__PSC+DV+Neurons+PM11__20190829_v3/run_monocle3_major_v3.Rdata')
newlabel <- readRDS('../Normal_LG/rdata/label.Rds')
head(newlabel)

head(pData(cds))
pData(cds)$new_subclustering <- subset(newlabel, anno_simple != 'PSC' & anno_simple != 'DV' & anno_simple != 'Neurons' & anno_simple != 'RG')$new_subclustering
pData(cds)$anno_simple <- subset(newlabel, anno_simple != 'PSC' & anno_simple != 'DV' & anno_simple != 'Neurons' & anno_simple != 'RG')$anno_simple
pData(cds) <- droplevels(pData(cds))

umap <- data.frame(cds@reducedDims$UMAP)
colnames(umap) <- c('UMAP1', 'UMAP2', 'UMAP3')

head(pData(cds))
umap$timepoint <- pData(cds)$timepoint
umap$anno_simple <- pData(cds)$anno_simple
umap$new_subclustering <- pData(cds)$new_subclustering
head(umap)


### subclustering ###
subColors <- mapvalues(umap$new_subclustering, from = levels(umap$new_subclustering), 
                       c('#dee6d4', '#b9dcb1', '#78c0ac', '#359bb7', '#1e659b', '#1c306d',
                         '#FCB17B', '#F16C4B', '#C81C12', '#7F0000',
                         '#e7c693', '#e18d30', '#71c0b0', '#177e7d', '#a4a4a4', '#1a1a1a'))

plot3d(umap[,c(1,3,2)], col = subColors, xlab = '', ylab = '', zlab = '')
decorate3d(box = F, xlab = 'UMAP 1', ylab = 'UMAP 2', zlab = 'UMAP 3')
view3d(theta = 140, phi = 10, zoom = .9)
#rgl.postscript(filename = "trajectory/traj.merged.3d_new_subclustering.pdf", fmt = "pdf", drawText = TRUE)



### PH ###
subColors <- mapvalues(umap$new_subclustering, from = levels(umap$new_subclustering), 
                       c('#dee6d4', '#b9dcb1', '#78c0ac', '#359bb7', '#1e659b', '#1c306d',
                         "grey90", "grey90", "grey90", "grey90",
                         "grey90", "grey90", "grey90", "grey90", "grey90", "grey90"))

plot3d(umap[,c(1,3,2)], col = subColors, xlab = '', ylab = '', zlab = '')
decorate3d(box = F, xlab = 'UMAP 1', ylab = 'UMAP 2', zlab = 'UMAP 3')
view3d(theta = 140, phi = 10, zoom = .9)
#rgl.postscript(filename = "trajectory/traj.merged.3d_new_subclustering_PH.pdf", fmt = "pdf", drawText = TRUE)

### PM ###
subColors <- mapvalues(umap$new_subclustering, from = levels(umap$new_subclustering), 
                       c("grey90", "grey90", "grey90", "grey90", "grey90", "grey90", 
                         '#FCB17B', '#F16C4B', '#C81C12', '#7F0000',
                         "grey90", "grey90", "grey90", "grey90", "grey90", "grey90"))

plot3d(umap[,c(1,3,2)], col = subColors, xlab = '', ylab = '', zlab = '')
decorate3d(box = F, xlab = 'UMAP 1', ylab = 'UMAP 2', zlab = 'UMAP 3')
view3d(theta = 140, phi = 10, zoom = .9)
#rgl.postscript(filename = "trajectory/traj.merged.3d_new_subclustering_PM.pdf", fmt = "pdf", drawText = TRUE)


