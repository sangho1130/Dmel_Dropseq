library(monocle3)
library(Matrix)
library(ggplot2)
library(Seurat)
library(pheatmap)
library(plyr)
#options(scipen = 100)

load('../../../Drop-seq_alignment/2.monocle3/filter__PSC+DV+Neurons+PM11__20190829_v3/run_monocle3_major_v3.Rdata')
newlabel <- readRDS('../Normal_LG/rdata/label.Rds')
head(newlabel)

head(pData(cds))
pData(cds)$new_subclustering <- subset(newlabel, anno_simple != 'PSC' & anno_simple != 'DV' & anno_simple != 'Neurons' & anno_simple != 'RG')$new_subclustering
pData(cds)$anno_simple <- subset(newlabel, anno_simple != 'PSC' & anno_simple != 'DV' & anno_simple != 'Neurons' & anno_simple != 'RG')$anno_simple
pData(cds) <- droplevels(pData(cds))


###
via_pseudotime = graph_test(cds, neighbor_graph="principal_graph", cores=3)
head(via_pseudotime)
deg_ids_via_pseudotime = row.names(subset(via_pseudotime, q_value < 0.05))

gene_module_df <- find_gene_modules(cds[deg_ids_via_pseudotime,], resolution=c(0,10^seq(-6,-1)))
#write.table(gene_module_df, 'genemodules/gene_module_df.txt', row.names = F, col.names = T, sep = '\t', quote = F)

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$new_subclustering)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
pheatmap(agg_mat, scale="column", clustering_method="ward.D2", color = colorRampPalette(c("blue", 'grey70', "red"))(n = 100), filename = 'genemodules/heatmap.by__new_subclustering__.pdf')

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$anno_simple)
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
pheatmap(agg_mat, scale="column", clustering_method="ward.D2", color = colorRampPalette(c("blue", 'grey70', "red"))(n = 100), filename = 'genemodules/heatmap.by__anno_simple__.pdf')


head(pData(cds))
# remove 30, 41, 42, 44, 46, 48, 49, 50, 51, 52
gene_module_df_sub <- subset(gene_module_df, module != 30 & module != 41 & module != 42 & module != 44 & module != 46 & 
                           module != 48 & module != 49 & module != 50 & module != 51 & module != 52)
gene_module_df_sub <- droplevels(gene_module_df_sub)

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$new_subclustering)
agg_mat <- aggregate_gene_expression(cds, gene_module_df_sub, cell_group_df)
row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
pheatmap(agg_mat, scale="column", clustering_method="ward.D2", border_color = NA,
         color = colorRampPalette(c("blue", 'grey70', "red"))(n = 100),
         filename = 'genemodules/heatmap.by__new_subclustering__rmModules.pdf')

cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$anno_simple)
agg_mat <- aggregate_gene_expression(cds, gene_module_df_sub, cell_group_df)
row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
pheatmap(agg_mat, scale="column", clustering_method="ward.D2", border_color = NA,
         color = colorRampPalette(c("blue", 'grey70', "red"))(n = 100), 
         filename = 'genemodules/heatmap.by__anno_simple__rmModules.pdf')

