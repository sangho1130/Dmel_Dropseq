library(monocle3)
library(ggplot2)
library(plyr)
library(dplyr)
library(pheatmap)
library(rgl)

#dir.create('trajectory_lm')
cds <- readRDS('tmp/cds_new_subclustering.Rds')

### Choose LM region ###
plot_cells(cds, color_cells_by = "labelTransfer", x = 1, y = 3, label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE, group_label_size = 3, cell_size = 1) +
  scale_color_manual(values = c('#dee6d4', '#b9dcb1', '#78c0ac', '#359bb7', '#95000C', '#E7C693', '#E18D30', '#04635C', '#A4A4A4')) + 
  geom_abline(slope = 8, intercept = 16, color = 'red2') + 
  geom_vline(xintercept = -1.5, color = 'blue')# + facet_wrap(~labelTransfer)

head(pData(cds))
umap <- data.frame(cds@reducedDims$UMAP); colnames(umap) <- c('UMAP1', 'UMAP2', 'UMAP3')
ggplot(umap, aes(UMAP1, UMAP3)) + geom_point(color = 'grey60') + geom_abline(slope = 8, intercept = 16, color = 'red2') + geom_vline(xintercept = -1.5, color = 'blue')
#umap <- subset(umap, UMAP1 <= -1.5)
umap <- subset(umap, UMAP3 >= 8*UMAP1 + 16) # v2
usecells <- rownames(umap)

cds_subset <- cds[, intersect(usecells, rownames(subset(pData(cds), labelTransfer_simple != 'GST-rich' & labelTransfer != 'PH 3' )))]
#saveRDS(cds_subset, 'tmp/cds_subset.Rds')
#cds_subset <- readRDS('tmp/cds_subset.Rds')
pData(cds_subset) <- droplevels(pData(cds_subset))
summary(pData(cds_subset)$labelTransfer)
plot_cells(cds_subset, x = 1, y = 3, color_cells_by = "labelTransfer", label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE, show_trajectory_graph=F, group_label_size = 7)# + geom_vline(xintercept = -1.5)

plot_cells(cds_subset, color_cells_by = "labelTransfer", x = 1, y = 3, label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE, group_label_size = 3, cell_size = 1) +
  scale_color_manual(values = c('#359bb7', '#95000C', '#E7C693', '#E18D30'))


subset_pr_test_res = graph_test(cds_subset, neighbor_graph="principal_graph", cores=3)
#saveRDS(subset_pr_test_res, 'tmp/lmsubset_graphtest.Rds')
pr_deg_ids = row.names(subset(subset_pr_test_res, q_value < 0.05))

gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=1e-2)
#saveRDS(gene_module_df, 'tmp/lmsubset_gene_module_df.Rds')

cell_group_df <- tibble::tibble(cell=row.names(colData(cds_subset)), cell_group=colData(cds_subset)$labelTransfer)
agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df, cell_group_df)
module_dendro <- hclust(dist(agg_mat))
gene_module_df$module <- factor(gene_module_df$module, levels = row.names(agg_mat)[module_dendro$order])

plot_cells(cds_subset, x = 1, y = 3, genes=gene_module_df, label_cell_groups=FALSE, show_trajectory_graph=FALSE)
pheatmap(agg_mat, scale="column", clustering_method="ward.D2", color = colorRampPalette(c("blue", 'grey60', "red"))(n = 100))#, filename = 'trajectory_lm/inset_gene_module_df_labelTransfer.pdf'); dev.off()
#write.table(gene_module_df, 'trajectory_lm/inset_gene_module_df_labelTransfer.txt', row.names = F, col.names = T, sep = '\t', quote = F)


### remove gene modules 
gene_module_df <- readRDS('tmp/lmsubset_gene_module_df.Rds')
gene_module_df <- subset(gene_module_df, module != 12 & module != 16 & module != 23)
gene_module_df <- droplevels(gene_module_df)

cell_group_df <- tibble::tibble(cell=row.names(colData(cds_subset)), cell_group=colData(cds_subset)$labelTransfer)
agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df, cell_group_df)
module_dendro <- hclust(dist(agg_mat))
gene_module_df$module <- factor(gene_module_df$module, levels = row.names(agg_mat)[module_dendro$order])

plot_cells(cds_subset, x = 1, y = 3, genes=gene_module_df, label_cell_groups=FALSE, show_trajectory_graph=FALSE)
ggsave('trajectory_lm/inset_gene_module_df_labelTransfer_flt_expression.pdf', units = 'cm', width = 25, height = 20)

agg_mat <- agg_mat[, c(3,4,1,2)]
pheatmap(agg_mat, scale="column", cluster_cols = F, clustering_method="ward.D2", color = colorRampPalette(c("blue", 'grey60', "red"))(n = 100))#, filename = 'trajectory_lm/inset_gene_module_df_labelTransfer_fltv2.pdf'); dev.off()



library(rgl)
umap <- umap[intersect(usecells, rownames(subset(pData(cds), labelTransfer_simple != 'GST-rich' & labelTransfer != 'PH 3' ))), ]
umap$labelTransfer <- pData(cds_subset)$labelTransfer
umap <- droplevels(umap)
subColors <- mapvalues(umap$labelTransfer, from = levels(umap$labelTransfer), c('#359bb7', '#95000C', '#E7C693', '#E18D30'))

### subColors

### principal graph ###
x <- 2
y <- 3
z <- 1

ica_space_df <- t(cds_subset@principal_graph_aux[['UMAP']]$dp_mst) %>% 
  as.data.frame() %>% 
  dplyr::select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y, prin_graph_dim_3 = z) %>% 
  dplyr::mutate(sample_name = rownames(.), sample_state = rownames(.))

dp_mst <- cds_subset@principal_graph[['UMAP']]
edge_df <- dp_mst %>% igraph::as_data_frame() %>% 
  dplyr::select_(source = "from", target = "to") %>% 
  dplyr::left_join(ica_space_df %>% 
                     dplyr::select_(source = "sample_name", source_prin_graph_dim_1 = "prin_graph_dim_1", source_prin_graph_dim_2 = "prin_graph_dim_2", source_prin_graph_dim_3 = "prin_graph_dim_3"), by = "source") %>% 
  dplyr::left_join(ica_space_df %>% 
                     dplyr::select_(target = "sample_name", target_prin_graph_dim_1 = "prin_graph_dim_1", target_prin_graph_dim_2 = "prin_graph_dim_2", target_prin_graph_dim_3 = "prin_graph_dim_3"), by = "target")
head(edge_df)

plot3d(umap[, c(2,3,1)], col = subColors, xlab = '', ylab = '', zlab = '')
segments3d(x=as.vector(t(edge_df[,c(3,6)])), y=as.vector(t(edge_df[,c(4,7)])), z=as.vector(t(edge_df[,c(5,8)])))
decorate3d(box = F)
view3d(theta = 90, phi = 10, zoom = .9) # need to be reflected (horizontally)
#rgl.postscript(filename = "trajectory_lm/3d_inset_lamellocytes.pdf", fmt = "pdf", drawText = TRUE)

plot_cells(cds_subset,
           color_cells_by = "labelTransfer", x = 1, y = 3, label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE, group_label_size = 3, cell_size = 1) +
  scale_color_manual(values = c('#359bb7', '#95000C', '#E7C693', '#E18D30')) + 
  geom_abline(slope = 8, intercept = 16, size = 2)
#ggsave('trajectory_lm/2d_inset_lamellocytes.pdf', units = 'cm', height = 10, width = 11)

rgl.close()
