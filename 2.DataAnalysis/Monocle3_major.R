library(monocle3)
library(Matrix)
library(ggplot2)
library(Seurat)
library(pheatmap)
library(plyr)
options(scipen = 100)

labelData <- read.delim('../../1.seurat3_alignment_withMuscle_regress-Library-nUMI/__filtered__SCrm__seurat3.3tps.newlabels.txt', header = T, sep = '\t', check.names = F, row.names = 1)
head(labelData); nrow(labelData) #19439
unique(labelData$Subclustering)

labelData <- subset(labelData, Subclustering != 'RG' & Subclustering != 'PSC' & Subclustering != 'DV' & Subclustering != 'Neurons')
nrow(labelData) #19143
labelData$orig.ident <- NULL
labelData$orig.anno <- NULL
labelData$RNA_snn_res.0.8_anno <- NULL
labelData$Subclustering <- NULL
labelData <- droplevels(labelData)
head(labelData)

labelData$anno_simple <- as.character(mapvalues(labelData$RNA_snn_res.0.8_anno_simple, from = levels(labelData$RNA_snn_res.0.8_anno_simple), to = c('CC', 'LM', 'PH', 'PM', 'GST-rich', 'Adipohemocyte')))
labelData$RNA_snn_res.0.8_anno_simple <- NULL
labelData[rownames(subset(labelData, Subclustering_v2 == 'SC-like')), 'anno_simple'] <- 'SC-like'
labelData$anno_simple <- factor(labelData$anno_simple, levels = c('SC-like', 'PH', 'PM', 'LM', 'CC', 'GST-rich', 'Adipohemocyte'))
levels(labelData$anno_simple)

labelData$Subclustering_v2 <- mapvalues(labelData$Subclustering_v2, from = levels(labelData$Subclustering_v2),
                                        to = c('CC 1', 'CC 2', 'LM 1', 'LM 2', 'PH 10', 'PH 11', 'PH 2', 'PH 3', 'PH 4', 'PH 5', 'PH 6', 'PH 7', 'PH 8', 'PH 9',
                                               'PM 1', 'PM 10', 'PM 2', 'PM 3', 'PM 4', 'PM 5', 'PM 6', 'PM 7', 'PM 8', 'PM 9', 'SC-like', 'GST-rich', 'Adipohemocyte'))
labelData$Subclustering_v2 <- factor(labelData$Subclustering_v2,
                                     levels = c('SC-like', 'PH 2', 'PH 3', 'PH 4', 'PH 5', 'PH 6', 'PH 7', 'PH 8', 'PH 9', 'PH 10', 'PH 11', 
                                                'PM 1', 'PM 2', 'PM 3', 'PM 4', 'PM 5', 'PM 6', 'PM 7', 'PM 8', 'PM 9', 'PM 10', 'LM 1', 'LM 2', 'CC 1', 'CC 2', 'GST-rich', 'Adipohemocyte'))
levels(labelData$Subclustering_v2)


exprs <- read.delim('merged.3tps.expr.allGenes.txt', row.names = 1, header = T, check.names = F, sep = '\t')
exprs_flt <- exprs[, rownames(labelData)]; identical(rownames(labelData), colnames(exprs_flt))

gene_metadata <- data.frame(row.names = rownames(exprs_flt), gene_short_name = rownames(exprs_flt))
cds <- new_cell_data_set(as(as.matrix(exprs_flt), "sparseMatrix"), cell_metadata = labelData, gene_metadata = gene_metadata)
head(exprs(cds))
head(pData(cds))
pData(cds)$timepoint <- factor(pData(cds)$timepoint, levels = c('AEL72hr', 'AEL96hr', 'AEL120hr'))

### Pre-process the data
cds <- preprocess_cds(cds, num_dim = 75, residual_model_formula_str = "~Library + nCount_RNA + percent.mt")
plot_pc_variance_explained(cds)
ggsave('stats/plot_pc_variance_explained.pdf', units = 'cm', height = 6, width = 12)


cds <- reduce_dimension(cds, max_components = 3, umap.min_dist = 0.4)
plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "anno_simple", x = 1, y = 2, group_label_size = 3, show_trajectory_graph = F, cell_size = 1) +
  scale_color_manual(values = c('#f84417', '#1c6aa4', '#95000c', '#ea9034', '#25a9b0', '#a4a4a4', 'black'))
plot_cells_3d(cds, color_cells_by="anno_simple", show_trajectory_graph = F)
plot_cells_3d(cds, color_cells_by = "Subclustering_v2", show_trajectory_graph = F)


plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "percent.mt", x = 1, y = 2, group_label_size = 3, show_trajectory_graph = F, cell_size = 1)
ggsave('trajectory/traj.merged.dim1_2.by__percent.mt__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 13)
plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "percent.mt", x = 1, y = 3, group_label_size = 3, show_trajectory_graph = F, cell_size = 1)
ggsave('trajectory/traj.merged.dim1_3.by__percent.mt__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 13)
plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "percent.mt", x = 2, y = 3, group_label_size = 3, show_trajectory_graph = F, cell_size = 1)
ggsave('trajectory/traj.merged.dim2_3.by__percent.mt__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 13)

plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "anno_simple", x = 1, y = 2, group_label_size = 3, show_trajectory_graph = F, cell_size = 1) +
  scale_color_manual(values = c('#f84417', '#1c6aa4', '#95000c', '#ea9034', '#25a9b0', '#a4a4a4', 'black'))
ggsave('trajectory/traj.merged.dim1_2.by__anno_simple__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 11)
plt <- plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "anno_simple", x = 1, y = 2, group_label_size = 3, show_trajectory_graph = F, label_cell_groups = F, cell_size = 1) +
  scale_color_manual(values = c('#f84417', '#1c6aa4', '#95000c', '#ea9034', '#25a9b0', '#a4a4a4', 'black')) + theme_void() + theme(legend.position = 'None');plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('trajectory/traj.merged.dim1_2.by__anno_simple__.min_dist__0.4__.augment.pdf', units = 'cm', height = 4, width = 4)

plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "anno_simple", x = 1, y = 3, group_label_size = 3, show_trajectory_graph = F, cell_size = 1) +
  scale_color_manual(values = c('#1c6aa4', '#1c6aa4', '#95000c', '#ea9034', '#25a9b0', '#a4a4a4', 'black'))
ggsave('trajectory/traj.merged.dim1_3.by__anno_simple__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 11)
plt <- plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "anno_simple", x = 1, y = 3, group_label_size = 3, show_trajectory_graph = F, label_cell_groups = F, cell_size = 1) +
  scale_color_manual(values = c('#1c6aa4', '#1c6aa4', '#95000c', '#ea9034', '#25a9b0', '#a4a4a4', 'black')) + theme_void() + theme(legend.position = 'None');plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('trajectory/traj.merged.dim1_3.by__anno_simple__.min_dist__0.4__.augment.pdf', units = 'cm', height = 4, width = 4)

plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "anno_simple", x = 2, y = 3, group_label_size = 3, show_trajectory_graph = F, cell_size = 1) +
  scale_color_manual(values = c('#f84417', '#1c6aa4', '#95000c', '#ea9034', '#25a9b0', '#a4a4a4', 'black'))
ggsave('trajectory/traj.merged.dim2_3.by__anno_simple__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 11)
plt <- plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "anno_simple", x = 2, y = 3, group_label_size = 3, show_trajectory_graph = F, label_cell_groups = F, cell_size = 1) +
  scale_color_manual(values = c('#f84417', '#1c6aa4', '#95000c', '#ea9034', '#25a9b0', '#a4a4a4', 'black')) + theme_void() + theme(legend.position = 'None');plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('trajectory/traj.merged.dim2_3.by__anno_simple__.min_dist__0.4__.augment.pdf', units = 'cm', height = 4, width = 4)



plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_v2", x = 1, y = 2, group_label_size = 3, show_trajectory_graph = F, cell_size = 1) +
  scale_color_manual(values = c('#f84417', '#cbe0c3', '#b9dcb1', '#9bcba3', '#78c0ac', '#54b0ba', '#359bb7', '#257fb2', '#1e659b', '#224a8d', '#1c306d',
                                '#f5e7c9', '#ecd5a5', '#eac487', '#f6ae72', '#f38956', '#ea6740', '#e14634', '#cf231c', '#b30a11', '#901407',
                                '#e7c693', '#e18d30', '#71c0b0', '#177e7d', '#a4a4a4', 'black'))
ggsave('trajectory/traj.merged.dim1_2.by__Subclustering_v2__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 11)
plt <- plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_v2", x = 1, y = 2, group_label_size = 3, show_trajectory_graph = F, label_cell_groups = F, cell_size = 1) +
  scale_color_manual(values = c('#f84417', '#cbe0c3', '#b9dcb1', '#9bcba3', '#78c0ac', '#54b0ba', '#359bb7', '#257fb2', '#1e659b', '#224a8d', '#1c306d',
                                '#f5e7c9', '#ecd5a5', '#eac487', '#f6ae72', '#f38956', '#ea6740', '#e14634', '#cf231c', '#b30a11', '#901407',
                                '#e7c693', '#e18d30', '#71c0b0', '#177e7d', '#a4a4a4', 'black')) + theme_void() + theme(legend.position = 'None');plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('trajectory/traj.merged.dim1_2.by__Subclustering_v2__.min_dist__0.4__.augment.pdf', units = 'cm', height = 4, width = 4)
plt <- plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_v2", x = 1, y = 2, group_label_size = 3, show_trajectory_graph = F, label_cell_groups = F, cell_size = 1) +
  scale_color_manual(values = c('#f84417', '#cbe0c3', '#b9dcb1', '#9bcba3', '#78c0ac', '#54b0ba', '#359bb7', '#257fb2', '#1e659b', '#224a8d', '#1c306d',
                                '#f5e7c9', '#ecd5a5', '#eac487', '#f6ae72', '#f38956', '#ea6740', '#e14634', '#cf231c', '#b30a11', '#901407',
                                '#e7c693', '#e18d30', '#71c0b0', '#177e7d', '#a4a4a4', 'black')) + theme_void() + theme(legend.position = 'None') + facet_wrap(~Subclustering_v2, ncol = 4);plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('trajectory/traj.merged.dim1_2.by__Subclustering_v2__.min_dist__0.4__.sep.pdf', units = 'cm', height = 20, width = 20)

plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_v2", x = 1, y = 3, group_label_size = 3, show_trajectory_graph = F, cell_size = 1) +
  scale_color_manual(values = c('#f84417', '#cbe0c3', '#b9dcb1', '#9bcba3', '#78c0ac', '#54b0ba', '#359bb7', '#257fb2', '#1e659b', '#224a8d', '#1c306d',
                                '#f5e7c9', '#ecd5a5', '#eac487', '#f6ae72', '#f38956', '#ea6740', '#e14634', '#cf231c', '#b30a11', '#901407',
                                '#e7c693', '#e18d30', '#71c0b0', '#177e7d', '#a4a4a4', 'black'))
ggsave('trajectory/traj.merged.dim1_3.by__Subclustering_v2__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 11)
plt <- plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_v2", x = 1, y = 3, group_label_size = 3, show_trajectory_graph = F, label_cell_groups = F, cell_size = 1) +
  scale_color_manual(values = c('#f84417', '#cbe0c3', '#b9dcb1', '#9bcba3', '#78c0ac', '#54b0ba', '#359bb7', '#257fb2', '#1e659b', '#224a8d', '#1c306d',
                                '#f5e7c9', '#ecd5a5', '#eac487', '#f6ae72', '#f38956', '#ea6740', '#e14634', '#cf231c', '#b30a11', '#901407',
                                '#e7c693', '#e18d30', '#71c0b0', '#177e7d', '#a4a4a4', 'black')) + theme_void() + theme(legend.position = 'None');plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('trajectory/traj.merged.dim1_3.by__Subclustering_v2__.min_dist__0.4__.augment.pdf', units = 'cm', height = 4, width = 4)
plt <- plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_v2", x = 1, y = 3, group_label_size = 3, show_trajectory_graph = F, label_cell_groups = F, cell_size = 1) +
  scale_color_manual(values = c('#f84417', '#cbe0c3', '#b9dcb1', '#9bcba3', '#78c0ac', '#54b0ba', '#359bb7', '#257fb2', '#1e659b', '#224a8d', '#1c306d',
                                '#f5e7c9', '#ecd5a5', '#eac487', '#f6ae72', '#f38956', '#ea6740', '#e14634', '#cf231c', '#b30a11', '#901407',
                                '#e7c693', '#e18d30', '#71c0b0', '#177e7d', '#a4a4a4', 'black')) + theme_void() + theme(legend.position = 'None') + facet_wrap(~Subclustering_v2, ncol = 4);plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('trajectory/traj.merged.dim1_3.by__Subclustering_v2__.min_dist__0.4__.sep.pdf', units = 'cm', height = 20, width = 20)

plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_v2", x = 2, y = 3, group_label_size = 3, show_trajectory_graph = F, cell_size = 1) +
  scale_color_manual(values = c('#f84417', '#cbe0c3', '#b9dcb1', '#9bcba3', '#78c0ac', '#54b0ba', '#359bb7', '#257fb2', '#1e659b', '#224a8d', '#1c306d',
                                '#f5e7c9', '#ecd5a5', '#eac487', '#f6ae72', '#f38956', '#ea6740', '#e14634', '#cf231c', '#b30a11', '#901407',
                                '#e7c693', '#e18d30', '#71c0b0', '#177e7d', '#a4a4a4', 'black'))
ggsave('trajectory/traj.merged.dim2_3.by__Subclustering_v2__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 11)
plt <- plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_v2", x = 2, y = 3, group_label_size = 3, show_trajectory_graph = F, label_cell_groups = F, cell_size = 1) +
  scale_color_manual(values = c('#f84417', '#cbe0c3', '#b9dcb1', '#9bcba3', '#78c0ac', '#54b0ba', '#359bb7', '#257fb2', '#1e659b', '#224a8d', '#1c306d',
                                '#f5e7c9', '#ecd5a5', '#eac487', '#f6ae72', '#f38956', '#ea6740', '#e14634', '#cf231c', '#b30a11', '#901407',
                                '#e7c693', '#e18d30', '#71c0b0', '#177e7d', '#a4a4a4', 'black')) + theme_void() + theme(legend.position = 'None');plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('trajectory/traj.merged.dim2_3.by__Subclustering_v2__.min_dist__0.4__.augment.pdf', units = 'cm', height = 4, width = 4)
plt <- plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_v2", x = 2, y = 3, group_label_size = 3, show_trajectory_graph = F, label_cell_groups = F, cell_size = 1) +
  scale_color_manual(values = c('#f84417', '#cbe0c3', '#b9dcb1', '#9bcba3', '#78c0ac', '#54b0ba', '#359bb7', '#257fb2', '#1e659b', '#224a8d', '#1c306d',
                                '#f5e7c9', '#ecd5a5', '#eac487', '#f6ae72', '#f38956', '#ea6740', '#e14634', '#cf231c', '#b30a11', '#901407',
                                '#e7c693', '#e18d30', '#71c0b0', '#177e7d', '#a4a4a4', 'black')) + theme_void() + theme(legend.position = 'None') + facet_wrap(~Subclustering_v2, ncol = 4);plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('trajectory/traj.merged.dim2_3.by__Subclustering_v2__.min_dist__0.4__.sep.pdf', units = 'cm', height = 20, width = 20)


cds <- cluster_cells(cds, resolution = 0.001); plot_cells(cds, color_cells_by = "partition", x = 1, y = 2) # or 0.4/0.0009
cds <- learn_graph(cds); plot_cells(cds, color_cells_by = "partition", x = 1, y = 2)


plot_cells(cds, color_cells_by = "anno_simple", label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE, group_label_size = 3, cell_size = 1) +
  scale_color_manual(values = c('#f84417', '#1c6aa4', '#95000c', '#ea9034', '#25a9b0', '#a4a4a4', 'black'))
ggsave('trajectory/traj.merged.dim1_2.by__partition__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 11)
plot_cells(cds, color_cells_by = "anno_simple", x = 1, y = 3, label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE, group_label_size = 3, cell_size = 1) +
  scale_color_manual(values = c('#f84417', '#1c6aa4', '#95000c', '#ea9034', '#25a9b0', '#a4a4a4', 'black'))
ggsave('trajectory/traj.merged.dim1_3.by__partition__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 11)
plot_cells(cds, color_cells_by = "anno_simple", x = 2, y = 3, label_groups_by_cluster=FALSE, label_leaves=FALSE, label_branch_points=FALSE, group_label_size = 3, cell_size = 1) +
  scale_color_manual(values = c('#f84417', '#1c6aa4', '#95000c', '#ea9034', '#25a9b0', '#a4a4a4', 'black'))
ggsave('trajectory/traj.merged.dim2_3.by__partition__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 11)


cds = order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, show_trajectory_graph = FALSE, graph_label_size = 3, cell_size = 1)
ggsave('trajectory/traj.merged.dim1_2.by__pseudotime__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 13)
plt <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, show_trajectory_graph = FALSE, cell_size = 1) + theme_void() + theme(legend.position = 'None');plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('trajectory/traj.merged.dim1_2.by__pseudotime__.min_dist__0.4__.augment.pdf', units = 'cm', width = 4, height = 4)
plot_cells(cds, color_cells_by = "pseudotime", x = 1, y = 3, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, show_trajectory_graph = FALSE, graph_label_size = 3, cell_size = 1)
ggsave('trajectory/traj.merged.dim1_3.by__pseudotime__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 13)
plt <- plot_cells(cds, color_cells_by = "pseudotime", x = 1, y = 3, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, show_trajectory_graph = FALSE, cell_size = 1) + theme_void() + theme(legend.position = 'None');plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('trajectory/traj.merged.dim1_3.by__pseudotime__.min_dist__0.4__.augment.pdf', units = 'cm', width = 4, height = 4)
plot_cells(cds, color_cells_by = "pseudotime", x = 2, y = 3, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, show_trajectory_graph = FALSE, graph_label_size = 3, cell_size = 1)
ggsave('trajectory/traj.merged.dim2_3.by__pseudotime__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 13)
plt <- plot_cells(cds, color_cells_by = "pseudotime", x = 2, y = 3, label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, show_trajectory_graph = FALSE, cell_size = 1) + theme_void() + theme(legend.position = 'None');plt
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('trajectory/traj.merged.dim2_3.by__pseudotime__.min_dist__0.4__.augment.pdf', units = 'cm', width = 4, height = 4)



###
via_pseudotime = graph_test(cds, neighbor_graph="principal_graph", cores=3)
head(via_pseudotime)
deg_ids_via_pseudotime = row.names(subset(via_pseudotime, q_value < 0.05))

gene_module_df <- find_gene_modules(cds[deg_ids_via_pseudotime,], resolution=c(0,10^seq(-6,-1)))
write.table(gene_module_df, 'genemodules/gene_module_df.txt', row.names = F, col.names = T, sep = '\t', quote = F)

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$Subclustering_v2)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
pheatmap(agg_mat, scale="column", clustering_method="ward.D2", color = colorRampPalette(c("blue", 'grey70', "red"))(n = 100), filename = 'genemodules/heatmap.by__Subclustering_v2__.pdf')

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$anno_simple)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
pheatmap(agg_mat, scale="column", clustering_method="ward.D2", color = colorRampPalette(c("blue", 'grey70', "red"))(n = 100), filename = 'genemodules/heatmap.by__anno_simple__.pdf')

#save.image('run_monocle3_major_v3.Rdata')




##########################################################
##########################################################
library(monocle3)
library(Matrix)
library(ggplot2)
library(Seurat)
library(pheatmap)
options(scipen = 100)
load('run_monocle3_major_v3.Rdata')

head(pData(cds))

### PG ### 
pData(cds)$Subclustering_pg <- as.character(pData(cds)$Subclustering_v2)
pData(cds)$anno_simple
pData(cds)[rownames(subset(pData(cds), anno_simple != 'SC-like' & anno_simple != 'PH')), 'Subclustering_pg'] <- 'Others'
unique(pData(cds)$Subclustering_pg)
pData(cds)$Subclustering_pg <- factor(pData(cds)$Subclustering_pg, levels = c("SC-like", "PH 2", "PH 3", "PH 4", "PH 5", "PH 6", "PH 7", "PH 8", "PH 9", "PH 10", "PH 11", 'Others'))

plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_pg", x = 1, y = 2, group_label_size = 3, cell_size = .5, show_trajectory_graph = F, alpha = .5) + 
  scale_color_manual(values = c('#f84417', '#cbe0c3', '#b9dcb1', '#9bcba3', '#78c0ac', '#54b0ba', '#359bb7', '#257fb2', '#1e659b', '#224a8d', '#1c306d', 'grey95'))
ggsave('trajectory/traj.merged.dim1_2.__subclustering_PH__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 11)
plt <- plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_pg", x = 1, y = 2, labels_per_group = F, cell_size = 1, show_trajectory_graph = F, alpha = .5) + 
  scale_color_manual(values = c('#f84417', '#cbe0c3', '#b9dcb1', '#9bcba3', '#78c0ac', '#54b0ba', '#359bb7', '#257fb2', '#1e659b', '#224a8d', '#1c306d', 'grey95')) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('trajectory/traj.merged.dim1_2.__subclustering_PH__.min_dist__0.4__.augment.pdf', units = 'cm', height = 4, width = 4)

plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_pg", x = 1, y = 3, group_label_size = 3, cell_size = .5, show_trajectory_graph = F, alpha = .5) + 
  scale_color_manual(values = c('#f84417', '#cbe0c3', '#b9dcb1', '#9bcba3', '#78c0ac', '#54b0ba', '#359bb7', '#257fb2', '#1e659b', '#224a8d', '#1c306d', 'grey95'))
ggsave('trajectory/traj.merged.dim1_3.__subclustering_PH__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 11)
plt <- plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_pg", x = 1, y = 3, labels_per_group = F, cell_size = 1, show_trajectory_graph = F, alpha = .5) + 
  scale_color_manual(values = c('#f84417', '#cbe0c3', '#b9dcb1', '#9bcba3', '#78c0ac', '#54b0ba', '#359bb7', '#257fb2', '#1e659b', '#224a8d', '#1c306d', 'grey95')) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('trajectory/traj.merged.dim1_3.__subclustering_PH__.min_dist__0.4__.augment.pdf', units = 'cm', height = 4, width = 4)

plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_pg", x = 2, y = 3, group_label_size = 3, cell_size = .5, show_trajectory_graph = F, alpha = .5) + 
  scale_color_manual(values = c('#f84417', '#cbe0c3', '#b9dcb1', '#9bcba3', '#78c0ac', '#54b0ba', '#359bb7', '#257fb2', '#1e659b', '#224a8d', '#1c306d', 'grey95'))
ggsave('trajectory/traj.merged.dim2_3.__subclustering_PH__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 11)
plt <- plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_pg", x = 2, y = 3, labels_per_group = F, cell_size = 1, show_trajectory_graph = F, alpha = .5) + 
  scale_color_manual(values = c('#f84417', '#cbe0c3', '#b9dcb1', '#9bcba3', '#78c0ac', '#54b0ba', '#359bb7', '#257fb2', '#1e659b', '#224a8d', '#1c306d', 'grey95')) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('trajectory/traj.merged.dim2_3.__subclustering_PH__.min_dist__0.4__.augment.pdf', units = 'cm', height = 4, width = 4)

### PM ### 
pData(cds)$Subclustering_pm <- as.character(pData(cds)$Subclustering_v2)
pData(cds)[rownames(subset(pData(cds), anno_simple != 'PM')), 'Subclustering_pm'] <- 'Others'
unique(pData(cds)$Subclustering_pm)
pData(cds)$Subclustering_pm <- factor(pData(cds)$Subclustering_pm, levels = c("PM 1", "PM 2", "PM 3", "PM 4", "PM 5", "PM 6", "PM 7", "PM 8", "PM 9", "PM 10", "Others"))

plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_pm", x = 1, y = 2, group_label_size = 3, cell_size = .5, show_trajectory_graph = F, alpha = .5) + 
  scale_color_manual(values = c('#f5e7c9', '#ecd5a5', '#eac487', '#f6ae72', '#f38956', '#ea6740', '#e14634', '#cf231c', '#b30a11', '#901407', 'grey95'))
ggsave('trajectory/traj.merged.dim1_2.__subclustering_PM__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 11)
plt <- plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_pm", x = 1, y = 2, labels_per_group = F, cell_size = 1, show_trajectory_graph = F, alpha = .5) + 
  scale_color_manual(values = c('#f5e7c9', '#ecd5a5', '#eac487', '#f6ae72', '#f38956', '#ea6740', '#e14634', '#cf231c', '#b30a11', '#901407', 'grey95')) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('trajectory/traj.merged.dim1_2.__subclustering_PM__.min_dist__0.4__.augment.pdf', units = 'cm', height = 4, width = 4)

plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_pm", x = 1, y = 3, group_label_size = 3, cell_size = .5, show_trajectory_graph = F, alpha = .5) + 
  scale_color_manual(values = c('#f5e7c9', '#ecd5a5', '#eac487', '#f6ae72', '#f38956', '#ea6740', '#e14634', '#cf231c', '#b30a11', '#901407', 'grey95'))
ggsave('trajectory/traj.merged.dim1_3.__subclustering_PM__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 11)
plt <- plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_pm", x = 1, y = 3, labels_per_group = F, cell_size = 1, show_trajectory_graph = F, alpha = .5) + 
  scale_color_manual(values = c('#f5e7c9', '#ecd5a5', '#eac487', '#f6ae72', '#f38956', '#ea6740', '#e14634', '#cf231c', '#b30a11', '#901407', 'grey95')) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('trajectory/traj.merged.dim1_3.__subclustering_PM__.min_dist__0.4__.augment.pdf', units = 'cm', height = 4, width = 4)

plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_pm", x = 2, y = 3, group_label_size = 3, cell_size = .5, show_trajectory_graph = F, alpha = .5) + 
  scale_color_manual(values = c('#f5e7c9', '#ecd5a5', '#eac487', '#f6ae72', '#f38956', '#ea6740', '#e14634', '#cf231c', '#b30a11', '#901407', 'grey95'))
ggsave('trajectory/traj.merged.dim2_3.__subclustering_PM__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 11)
plt <- plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_pm", x = 2, y = 3, labels_per_group = F, cell_size = 1, show_trajectory_graph = F, alpha = .5) + 
  scale_color_manual(values = c('#f5e7c9', '#ecd5a5', '#eac487', '#f6ae72', '#f38956', '#ea6740', '#e14634', '#cf231c', '#b30a11', '#901407', 'grey95')) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('trajectory/traj.merged.dim2_3.__subclustering_PM__.min_dist__0.4__.augment.pdf', units = 'cm', height = 4, width = 4)

### CC ### 
pData(cds)$Subclustering_cc <- as.character(pData(cds)$Subclustering_v2)
pData(cds)[rownames(subset(pData(cds), anno_simple != 'CC')), 'Subclustering_cc'] <- 'Others'
unique(pData(cds)$Subclustering_cc)
pData(cds)$Subclustering_cc <- factor(pData(cds)$Subclustering_cc, levels = c("CC 1", "CC 2", 'Others'))

plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_cc", x = 1, y = 2, group_label_size = 3, cell_size = .5, show_trajectory_graph = F, alpha = .5) + 
  scale_color_manual(values = c('#71c0b0', '#177e7d', 'grey95'))
ggsave('trajectory/traj.merged.dim1_2.__subclustering_CC__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 11)
plt <- plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_cc", x = 1, y = 2, labels_per_group = F, cell_size = 1, show_trajectory_graph = F, alpha = .5) + 
  scale_color_manual(values = c('#71c0b0', '#177e7d', 'grey95')) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('trajectory/traj.merged.dim1_2.__subclustering_CC__.min_dist__0.4__.augment.pdf', units = 'cm', height = 4, width = 4)

plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_cc", x = 1, y = 3, group_label_size = 3, cell_size = .5, show_trajectory_graph = F, alpha = .5) + 
  scale_color_manual(values = c('#71c0b0', '#177e7d', 'grey95'))
ggsave('trajectory/traj.merged.dim1_3.__subclustering_CC__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 11)
plt <- plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_cc", x = 1, y = 3, labels_per_group = F, cell_size = 1, show_trajectory_graph = F, alpha = .5) + 
  scale_color_manual(values = c('#71c0b0', '#177e7d', 'grey95')) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('trajectory/traj.merged.dim1_3.__subclustering_CC__.min_dist__0.4__.augment.pdf', units = 'cm', height = 4, width = 4)

plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_cc", x = 2, y = 3, group_label_size = 3, cell_size = .5, show_trajectory_graph = F, alpha = .5) + 
  scale_color_manual(values = c('#71c0b0', '#177e7d', 'grey95'))
ggsave('trajectory/traj.merged.dim2_3.__subclustering_CC__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 11)
plt <- plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_cc", x = 2, y = 3, labels_per_group = F, cell_size = 1, show_trajectory_graph = F, alpha = .5) + 
  scale_color_manual(values = c('#71c0b0', '#177e7d', 'grey95')) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('trajectory/traj.merged.dim2_3.__subclustering_CC__.min_dist__0.4__.augment.pdf', units = 'cm', height = 4, width = 4)


### LM ### 
pData(cds)$Subclustering_lm <- as.character(pData(cds)$Subclustering_v2)
pData(cds)[rownames(subset(pData(cds), anno_simple != 'LM')), 'Subclustering_lm'] <- 'Others'
unique(pData(cds)$Subclustering_lm)
pData(cds)$Subclustering_lm <- factor(pData(cds)$Subclustering_lm, levels = c("LM 1", "LM 2", "Others"))

plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_lm", x = 1, y = 2, group_label_size = 3, cell_size = .5, show_trajectory_graph = F, alpha = .5) + 
  scale_color_manual(values = c('#e7c693', '#e18d30', 'grey95'))
ggsave('trajectory/traj.merged.dim1_2.__subclustering_LM__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 11)
plt <- plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_lm", x = 1, y = 2, labels_per_group = F, cell_size = 1, show_trajectory_graph = F, alpha = .5) + 
  scale_color_manual(values = c('#e7c693', '#e18d30', 'grey95')) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('trajectory/traj.merged.dim1_2.__subclustering_LM__.min_dist__0.4__.augment.pdf', units = 'cm', height = 4, width = 4)

plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_lm", x = 1, y = 3, group_label_size = 3, cell_size = .5, show_trajectory_graph = F, alpha = .5) + 
  scale_color_manual(values = c('#e7c693', '#e18d30', 'grey95'))
ggsave('trajectory/traj.merged.dim1_3.__subclustering_LM__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 11)
plt <- plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_lm", x = 1, y = 3, labels_per_group = F, cell_size = 1, show_trajectory_graph = F, alpha = .5) + 
  scale_color_manual(values = c('#e7c693', '#e18d30', 'grey95')) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('trajectory/traj.merged.dim1_3.__subclustering_LM__.min_dist__0.4__.augment.pdf', units = 'cm', height = 4, width = 4)

plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_lm", x = 2, y = 3, group_label_size = 3, cell_size = .5, show_trajectory_graph = F, alpha = .5) + 
  scale_color_manual(values = c('#e7c693', '#e18d30', 'grey95'))
ggsave('trajectory/traj.merged.dim2_3.__subclustering_LM__.min_dist__0.4__.pdf', units = 'cm', height = 10, width = 11)
plt <- plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Subclustering_lm", x = 2, y = 3, labels_per_group = F, cell_size = 1, show_trajectory_graph = F, alpha = .5) + 
  scale_color_manual(values = c('#e7c693', '#e18d30', 'grey95')) +
  theme_void() + theme(legend.position = 'None')
AugmentPlot(plt, dpi = 300, width = 4, height = 4)
ggsave('trajectory/traj.merged.dim2_3.__subclustering_LM__.min_dist__0.4__.augment.pdf', units = 'cm', height = 4, width = 4)




#########################################################
library(monocle3)
library(Matrix)
library(ggplot2)
library(Seurat)
library(pheatmap)
options(scipen = 100)
load('run_monocle3_major_v3.Rdata')

head(pData(cds))
plot_cells_3d(cds, color_cells_by = "Subclustering_v2", show_trajectory_graph = F)
plot_cells_3d(cds, color_cells_by = "anno_simple", show_trajectory_graph = F)
plot_cells_3d(cds, color_cells_by = "timepoint", show_trajectory_graph = F)

module_int <- subset(gene_module_df, module == 31)
module_int$id


library(rgl)
library(plyr)
library(dplyr)

umap <- data.frame(cds@reducedDims$UMAP)
colnames(umap) <- c('UMAP1', 'UMAP2', 'UMAP3')

head(pData(cds))
umap$timepoint <- pData(cds)$timepoint
umap$simple_v2 <- pData(cds)$anno_simple
umap$simple_v2 <- mapvalues(umap$simple_v2, from = levels(umap$simple_v2), to = c('PH', levels(umap$simple_v2)[2:7]))
umap$Subclustering_v2 <- pData(cds)$Subclustering_v2
head(umap)
simpleColors <- mapvalues(umap$simple_v2, from = levels(umap$simple_v2), c('#1c6aa4', '#95000c', '#ea9034', '#25a9b0', '#a4a4a4', '#1a1a1a'))

plot3d(umap[,c(1,3,2)], col = simpleColors, xlab = '', ylab = '', zlab = '')
decorate3d(box = F, xlab = 'UMAP 1', ylab = 'UMAP 2', zlab = 'UMAP 3')
legend3d('right', legend = levels(umap$simple_v2), col = c('#1c6aa4', '#95000c', '#ea9034', '#25a9b0', '#a4a4a4', '#1a1a1a'), pch = 16, cex=1.5, inset = c(0.01))

view3d(theta = 0, phi = 10, zoom = .9)
view3d(theta = 140, phi = 10, zoom = .9)
view3d(theta = 180, phi = 10, zoom = .9)
play3d(spin3d(axis = c(0, 1, 0), rpm = 2.5))

#rgl.postscript(filename = "trajectory/traj.merged.3d.pdf", fmt = "pdf", drawText = TRUE)
rgl.close()


#########################################################
#########################################################
#########################################################
library(monocle3)
library(Matrix)
library(ggplot2)
library(Seurat)
library(pheatmap)
options(scipen = 100)
load('run_monocle3_major_v3.Rdata')

head(pData(cds))
# remove 43, 48, 50
gene_module_df <- subset(gene_module_df, module != 43 & module != 48 & module != 50)
gene_module_df <- droplevels(gene_module_df)

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$Subclustering_v2)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
pheatmap(agg_mat, scale="column", clustering_method="ward.D2",
         color = colorRampPalette(c("blue", 'grey70', "red"))(n = 100),
         filename = 'genemodules/heatmap.by__Subclustering_v2__rmModules.pdf')

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$anno_simple)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
pheatmap(agg_mat, scale="column", clustering_method="ward.D2", 
         color = colorRampPalette(c("blue", 'grey70', "red"))(n = 100), 
         filename = 'genemodules/heatmap.by__anno_simple__rmModules.pdf')



