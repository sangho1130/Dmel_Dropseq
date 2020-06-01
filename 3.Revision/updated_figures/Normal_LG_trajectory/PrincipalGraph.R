library(monocle3)
library(Matrix)
library(ggplot2)
library(plyr)
library(rgl)

load('../../../Drop-seq_alignment/2.monocle3/filter__PSC+DV+Neurons+PM11__20190829_v3/run_monocle3_major_v3.Rdata')
newlabel <- readRDS('../Normal_LG/rdata/label.Rds')
head(newlabel)

head(pData(cds))
pData(cds)$new_subclustering <- subset(newlabel, anno_simple != 'PSC' & anno_simple != 'DV' & anno_simple != 'Neurons' & anno_simple != 'RG')$new_subclustering
pData(cds)$anno_simple <- subset(newlabel, anno_simple != 'PSC' & anno_simple != 'DV' & anno_simple != 'Neurons' & anno_simple != 'RG')$anno_simple
pData(cds) <- droplevels(pData(cds))
cds

x <- 1
y <- 3
z <- 2

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



umap <- data.frame(cds@reducedDims$UMAP)
colnames(umap) <- c('UMAP1', 'UMAP2', 'UMAP3')
umap$anno_simple <- pData(cds)$anno_simple
umap$new_subclustering <- pData(cds)$new_subclustering
umap$timepoint <- pData(cds)$timepoint
umap$pseudotime <- pseudotime(cds, reduction_method = 'UMAP')
head(umap)

# new_subcluster
useCols <- mapvalues(umap$new_subclustering, from = levels(umap$new_subclustering), 
                     c('#dee6d4', '#b9dcb1', '#78c0ac', '#359bb7', '#1e659b', '#1c306d',
                       '#FCB17B', '#F16C4B', '#C81C12', '#7F0000',
                       '#e7c693', '#e18d30', '#71c0b0', '#177e7d', '#a4a4a4', '#1a1a1a'))

plot3d(umap[,c(x, y, z)], col = useCols, xlab = '', ylab = '', zlab = '')
decorate3d(box = F, xlab = 'UMAP 1', ylab = 'UMAP 2', zlab = 'UMAP 3')
view3d(theta = 140, phi = 10, zoom = .9)

segments3d(x=as.vector(t(edge_df[,c(3,6)])), y=as.vector(t(edge_df[,c(4,7)])), z=as.vector(t(edge_df[,c(5,8)])))

#rgl.postscript(filename = "trajectory/traj.merged.3d_new_subclustering+principalGraph.pdf", fmt = "pdf", drawText = TRUE)
rgl.close()
