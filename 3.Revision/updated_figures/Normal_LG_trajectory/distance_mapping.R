library(monocle3)
library(ggplot2)
library(plyr)
library(ggridges)
library(viridis)
library(tidyr)
library(geosphere)
library(rgl)
options(scipen = 100)
load('../../../Drop-seq_alignment/2.monocle3/filter__PSC+DV+Neurons+PM11__20190829_v3/run_monocle3_major_v3.Rdata')
newlabel <- readRDS('../Normal_LG/rdata/label.Rds')
head(newlabel)

head(pData(cds))
pData(cds)$new_subclustering <- subset(newlabel, anno_simple != 'PSC' & anno_simple != 'DV' & anno_simple != 'Neurons' & anno_simple != 'RG')$new_subclustering
pData(cds)$anno_simple <- subset(newlabel, anno_simple != 'PSC' & anno_simple != 'DV' & anno_simple != 'Neurons' & anno_simple != 'RG')$anno_simple
pData(cds) <- droplevels(pData(cds))

head(pData(cds))

cds_pseudotime <- data.frame(pseudotime = pseudotime(cds, reduction_method = 'UMAP'))
cds_pseudotime$timepoint <- pData(cds)$timepoint
cds_pseudotime$new_subclustering <- pData(cds)$new_subclustering
cds_pseudotime$anno_simple <- pData(cds)$anno_simple
cds_pseudotime <- subset(cds_pseudotime, pseudotime != Inf)
head(cds_pseudotime)


###
umap <- data.frame(cds@reducedDims$UMAP)
colnames(umap) <- c('UMAP1', 'UMAP2', 'UMAP3')
umap$timepoint <- cds_pseudotime$timepoint
umap$new_subclustering <- cds_pseudotime$new_subclustering
umap$anno_simple <- cds_pseudotime$anno_simple
head(umap)

avg_umap <- data.frame(row.names = levels(umap$new_subclustering), 
                       UMAP1 = numeric(16), UMAP2 = numeric(16), UMAP3 = numeric(16),
                       timepoint = character(16),
                       new_subclustering = levels(umap$new_subclustering),
                       anno_simple = levels(umap$new_subclustering))
avg_umap$timepoint <- 'merged'
avg_umap$anno_simple <- mapvalues(avg_umap$anno_simple,
                                  from = levels(avg_umap$anno_simple),
                                  to = c('Adipohemocyte', 'CC', 'CC', 'GST-rich', 'LM', 'LM', 
                                         'PH', 'PH', 'PH', 'PH', 'PH', 'PH',
                                         'PM', 'PM', 'PM', 'PM'))
head(avg_umap)

for (subtype in levels(umap$new_subclustering)){
  tmp <- subset(umap, new_subclustering == subtype)
  avg_umap[subtype, 'UMAP1'] <- mean(tmp$UMAP1)
  avg_umap[subtype, 'UMAP2'] <- mean(tmp$UMAP2)
  avg_umap[subtype, 'UMAP3'] <- mean(tmp$UMAP3)
}

umap$mean <- 'Original'
avg_umap$mean <- rownames(avg_umap)

umap$class <- 'Original'
avg_umap$class <- 'Mean'

umap$label <- NA
avg_umap$label <- rownames(avg_umap)

umap_wAvg <- rbind(umap, avg_umap)
umap_wAvg$mean <- factor(umap_wAvg$mean, levels = c('Original', levels(umap$new_subclustering)))
umap_wAvg$class <- factor(umap_wAvg$class, levels = c('Original', 'Mean'))
head(umap_wAvg)

ggplot(umap_wAvg, aes(UMAP1, UMAP2, col = mean)) + 
  geom_point(aes(shape = class)) + 
  geom_text(aes(label = label), hjust = 0, vjust = 0, col = 'black') +
  scale_color_manual(values = c('grey90',
                                '#dee6d4', '#b9dcb1', '#78c0ac', '#359bb7', '#1e659b', '#1c306d',
                                '#FCB17B', '#F16C4B', '#C81C12', '#7F0000',
                                '#e7c693', '#e18d30', '#71c0b0', '#177e7d', '#a4a4a4', '#1a1a1a')) + 
  scale_shape_manual(values = c(16, 8)) +
  theme_bw() + 
  theme(panel.grid = element_blank())
#ggsave('distance/traj.1_2.pdf', units = 'cm', width = 18, height = 12)


distmat <- data.matrix(dist(avg_umap[,1:3], diag = T, upper = T))


library(pheatmap)

callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

pheatmap(as.matrix(distmat), cluster_rows = T, cluster_cols = T, scale = 'none', treeheight_row = 20, treeheight_col = 20,  
         color = colorRampPalette(c("red", "grey60", 'white'))(n = 101),
         cellwidth = 5, cellheight = 5, fontsize_row = 6, fontsize_col = 6, border_color = NA,
         clustering_distance_rows = 'euclidean', clustering_method = 'complete', clustering_callback = callback)#,
         #filename = 'distance/distance_mapping.pdf')


distmat <- data.frame(subclusters = rownames(distmat), distmat, check.names = F, check.rows = F)
#write.table(distmat, 'distance/distance_mapping.txt', quote = F, sep = '\t', row.names = F, col.names = T)



###
plot3d(umap[,c(1,3,2)], col = 'grey90',  xlab = '', ylab = '', zlab = '')
decorate3d(box = F, xlab = 'UMAP 1', ylab = 'UMAP 2', zlab = 'UMAP 3')

umap_wAvg_label <- subset(umap_wAvg, mean != 'Original')
umap_wAvg_label <- droplevels(umap_wAvg_label)
head(umap_wAvg_label)
colors <- mapvalues(umap_wAvg_label$mean, from = levels(umap_wAvg_label$mean), to = c('#dee6d4', '#b9dcb1', '#78c0ac', '#359bb7', '#1e659b', '#1c306d',
                                                                                      '#FCB17B', '#F16C4B', '#C81C12', '#7F0000',
                                                                                      '#e7c693', '#e18d30', '#71c0b0', '#177e7d', '#a4a4a4', '#1a1a1a'))
plot3d(umap_wAvg_label[,c(1,3,2)], size = 10, col = colors,  xlab = '', ylab = '', zlab = '', add = T)
text3d(umap_wAvg_label[,c(1,3,2)], texts = umap_wAvg_label$label,  xlab = '', ylab = '', zlab = '')
view3d(theta = 140, phi = 10, zoom = .9)
#rgl.postscript(filename = "distance/traj.merged.3d_subcluster_mean_v2.pdf", fmt = "pdf", drawText = TRUE)

### add principal graph
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

segments3d(x=as.vector(t(edge_df[,c(3,6)])), y=as.vector(t(edge_df[,c(4,7)])), z=as.vector(t(edge_df[,c(5,8)])))
#rgl.postscript(filename = "distance/traj.merged.3d_subcluster_mean_v2+principalGraph.pdf", fmt = "pdf", drawText = TRUE)

rgl.close()


