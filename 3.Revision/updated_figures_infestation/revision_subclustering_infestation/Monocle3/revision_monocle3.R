library(monocle3)
library(ggplot2)
library(plyr)
library(dplyr)
library(pheatmap)
library(rgl)

#dir.create('tmp')
#load('../../Drop-seq_96hr_LymphGland_alignment/Monocle3/allCells/run_monocle3_allCells.Rdata')
#newlabel <- readRDS('../tmp/infested_label.Rds')
#newlabel <- newlabel[rownames(pData(cds)), ]
#newlabel <- droplevels(newlabel)

#head(pData(cds))
#pData(cds)$labelTransfer <- newlabel$labelTransfer
#pData(cds)$labelTransfer_simple <- newlabel$labelTransfer_simple
#saveRDS(cds, 'tmp/cds_new_subclustering.Rds')


#dir.create('trajectory')
cds <- readRDS('tmp/cds_new_subclustering.Rds')


### principal graph ###
x <- 2
y <- 3
z <- 1

ica_space_df <- t(cds@principal_graph_aux[['UMAP']]$dp_mst) %>% 
  as.data.frame() %>% 
  dplyr::select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y, prin_graph_dim_3 = z) %>% 
  dplyr::mutate(sample_name = rownames(.), sample_state = rownames(.))

dp_mst <- cds@principal_graph[['UMAP']]
edge_df <- dp_mst %>% igraph::as_data_frame() %>% 
  dplyr::select_(source = "from", target = "to") %>% 
  dplyr::left_join(ica_space_df %>% 
                     dplyr::select_(source = "sample_name", source_prin_graph_dim_1 = "prin_graph_dim_1", source_prin_graph_dim_2 = "prin_graph_dim_2", source_prin_graph_dim_3 = "prin_graph_dim_3"), by = "source") %>% 
  dplyr::left_join(ica_space_df %>% 
                     dplyr::select_(target = "sample_name", target_prin_graph_dim_1 = "prin_graph_dim_1", target_prin_graph_dim_2 = "prin_graph_dim_2", target_prin_graph_dim_3 = "prin_graph_dim_3"), by = "target")
head(edge_df)


### umap coord ###
umap <- data.frame(cds@reducedDims$UMAP)
colnames(umap) <- c('UMAP1', 'UMAP2', 'UMAP3')

head(pData(cds))
umap$labelTransfer <- pData(cds)$labelTransfer
umap$labelTransfer_simple <- pData(cds)$labelTransfer_simple
head(umap)

simpleColors <- mapvalues(umap$labelTransfer_simple, from = levels(umap$labelTransfer_simple), c('#1C6AA4', '#95000C', '#EA9034', '#51C3AE', '#A4A4A4'))
subColors <- mapvalues(umap$labelTransfer, from = levels(umap$labelTransfer), c('#dee6d4', '#b9dcb1', '#78c0ac', '#359bb7', '#95000C', '#E7C693', '#E18D30', '#04635C', '#A4A4A4'))

### subColors
plot3d(umap[, c(2,3,1)], col = subColors)
segments3d(x=as.vector(t(edge_df[,c(3,6)])), y=as.vector(t(edge_df[,c(4,7)])), z=as.vector(t(edge_df[,c(5,8)])))
decorate3d(box = F)
#view3d(theta = 180, phi = 110, zoom = .9) #rotate

view3d(theta = 90, phi = 10, zoom = .9) # need to be reflected (horizontally)
#play3d(spin3d(axis = c(0, 0, 1), rpm = 2.5))
#rgl.postscript(filename = "trajectory/3d_labelTransfer.pdf", fmt = "pdf", drawText = TRUE)


### simpleColors
plot3d(umap[, c(2,3,1)], col = simpleColors)
segments3d(x=as.vector(t(edge_df[,c(3,6)])), y=as.vector(t(edge_df[,c(4,7)])), z=as.vector(t(edge_df[,c(5,8)])))
decorate3d(box = F)
#view3d(theta = 180, phi = 110, zoom = .9) #rotate

view3d(theta = 90, phi = 10, zoom = .9) # need to be reflected (horizontally)
#play3d(spin3d(axis = c(0, 0, 1), rpm = 2.5))
#rgl.postscript(filename = "trajectory/3d_labelTransfer_simple.pdf", fmt = "pdf", drawText = TRUE)


### pseudotime ###
pseudotime(cds, reduction_method = 'UMAP')
identical(rownames(umap), names(pseudotime(cds, reduction_method = 'UMAP')))
umap$pseudotime <- pseudotime(cds, reduction_method = 'UMAP')
head(umap)

myColorRamp <- function(colors, values) {
  v <- (values - min(values))/diff(range(values))
  x <- colorRamp(colors)(v)
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}
cols <- myColorRamp(c("#0d0887", '#8305a7', '#ef7e50', "#f0f921"), umap$pseudotime) 

plot3d(umap[, c(2,3,1)], col = cols, xlab = '', ylab = '', zlab = '')
segments3d(x=as.vector(t(edge_df[,c(3,6)])), y=as.vector(t(edge_df[,c(4,7)])), z=as.vector(t(edge_df[,c(5,8)])))

#view3d(theta = 180, phi = 110, zoom = .9) #
view3d(theta = 90, phi = 10, zoom = .9) # need to be reflected (horizontally)
#rgl.postscript(filename = "trajectory/3d_pseudotime.pdf", fmt = "pdf", drawText = TRUE)

rgl.close()
