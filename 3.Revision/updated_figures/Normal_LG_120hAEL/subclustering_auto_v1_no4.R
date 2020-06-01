#################
### Functions ###
#################

# 3. res v3 - excluding least correlated cells 
calc_res <- function(inputObj, outputPath, pcNum, scaleFactor, rangeMax = 2, cormethod = 'spearman') {
  dir.create(paste0(c(outputPath, 'res'), collapse = '/'))
  dir.create(paste0(c(outputPath, 'res/tsne'), collapse = '/'))
  dir.create(paste0(c(outputPath, 'res/umap'), collapse = '/'))
  dir.create(paste0(c(outputPath, 'res/heatmap'), collapse = '/'))

  inputObj <- FindNeighbors(inputObj, dims = 1:pcNum)
  clusterCount <- data.frame(matrix(ncol = 8, nrow = 0))
  
  endpoint <- rangeMax*scaleFactor
  for (res in c(1:endpoint)){
    res <- res/scaleFactor
    inputObj <- FindClusters(inputObj, resolution = res)
    clusternum <- length(unique(inputObj@meta.data$seurat_clusters))
    #####################
    if (clusternum != 1) {
      vardf <- data.frame(as.matrix(GetAssayData(inputObj, slot = 'data')), check.rows = F, check.names = F)
      vardf <- vardf[VariableFeatures(inputObj), ]

      emptydf <- data.frame(row.names = VariableFeatures(inputObj), matrix(nrow = 2000, ncol = length(unique(inputObj@meta.data$seurat_clusters))))
      for (tmpcluster in levels(inputObj@meta.data$seurat_clusters)) {
        tmp_obj <- subset(inputObj, idents = tmpcluster)
        tmp_cells <- rownames(tmp_obj@meta.data)
        tmp_pseudo <- rowSums(vardf[, tmp_cells])
        emptydf[,as.numeric(tmpcluster) + 1] <- tmp_pseudo
      }
      colnames(emptydf) <- levels(inputObj@meta.data$seurat_clusters)
      
      corMat <- data.frame(cor(emptydf, method = cormethod), check.names = F) # pearson spearman
      pseudo_exprs <- emptydf; head(pseudo_exprs)
      inputObj@meta.data$supergroup <- as.character(inputObj@meta.data$seurat_clusters)
      
      if (clusternum > 2) {
        round <- 1
        if (min(corMat) < 0.9 & ncol(corMat) >= 2){
          while (min(corMat) < 0.9 & ncol(corMat) >= 2) {
            leastcor <- names(which.min((rowSums(corMat)-1)/(ncol(corMat)-1)))
            
            if (round == 1) {
              inputObj_sub <- subset(inputObj, idents = setdiff(levels(Idents(inputObj)), leastcor))
            } else {
              inputObj_sub <- subset(inputObj_sub, idents = setdiff(levels(Idents(inputObj_sub)), leastcor))
            }
            inputObj_sub@meta.data <- droplevels(inputObj_sub@meta.data)
            inputObj_sub <- FindVariableFeatures(inputObj_sub, selection.method = "vst", nfeatures = 2000)
            
            inputObj_sub_pseudo <- data.frame(row.names = VariableFeatures(inputObj_sub), matrix(nrow = 2000, ncol = length(unique(inputObj_sub@meta.data$seurat_clusters))))

            for (tmpcluster_idx in c(1:ncol(inputObj_sub_pseudo))) {
              tmp_obj <- subset(inputObj_sub, idents = levels(inputObj_sub@meta.data$seurat_clusters)[tmpcluster_idx])
              tmp_cells <- rownames(tmp_obj@meta.data)
              tmp_pseudo <- rowSums(vardf[, tmp_cells])
              inputObj_sub_pseudo[,tmpcluster_idx] <- tmp_pseudo
            }
            
            colnames(inputObj_sub_pseudo) <- levels(inputObj_sub@meta.data$seurat_clusters)

            corMat <- data.frame(cor(inputObj_sub_pseudo, method = cormethod), check.names = F)
            round <- round + 1
          }
        } else {
          inputObj_sub <- inputObj
          inputObj_sub_pseudo <- pseudo_exprs
        }
      } else {
        inputObj_sub <- inputObj
        inputObj_sub_pseudo <- pseudo_exprs
      }
      
      diag(corMat) <- 0
      while (max(corMat) >= 0.95) {
        results <- data.frame(v1=character(0), v2=character(0), ccor=numeric(0), stringsAsFactors=FALSE)
        while (sum(corMat>0)>1) {
          maxval <- max(corMat)
          max <- which(corMat==maxval, arr.ind=TRUE)[1,]
          results <- rbind(results, data.frame(v1=rownames(corMat)[max[1]], v2=colnames(corMat)[max[2]], cor=maxval))
          corMat[max[1],] <- 0
          corMat[,max[1]] <- 0
          corMat[max[2],] <- 0
          corMat[,max[2]] <- 0
        }
        maxpair <- data.frame(results[1, c(1,2)])

        inputObj_sub@meta.data[rownames(subset(inputObj_sub@meta.data, supergroup == as.character(maxpair$v1))), 'supergroup'] <- paste0(c(as.character(maxpair$v1), as.character(maxpair$v2)), collapse = '_')
        inputObj_sub@meta.data[rownames(subset(inputObj_sub@meta.data, supergroup == as.character(maxpair$v2))), 'supergroup'] <- paste0(c(as.character(maxpair$v1), as.character(maxpair$v2)), collapse = '_')
        
        inputObj_sub_pseudo[, paste0(c(as.character(maxpair$v1), as.character(maxpair$v2)), collapse = '_')] <- inputObj_sub_pseudo[, as.character(maxpair$v1)] + inputObj_sub_pseudo[, as.character(maxpair$v2)]
        inputObj_sub_pseudo[, as.character(maxpair$v1)] <- NULL
        inputObj_sub_pseudo[, as.character(maxpair$v2)] <- NULL
        
        corMat <- cor(inputObj_sub_pseudo, method = cormethod)
        diag(corMat) <- 0
      }
      head(inputObj@meta.data)
      head(inputObj_sub@meta.data)
      Idents(inputObj) <- 'supergroup'
      Idents(inputObj_sub) <- 'supergroup'
      
    } else {
      inputObj_sub <- inputObj
      inputObj@meta.data$supergroup <- inputObj@meta.data$seurat_clusters
      inputObj_sub@meta.data$supergroup <- inputObj@meta.data$seurat_clusters
    }
    #####################
    
    inputObj@meta.data[rownames(inputObj_sub@meta.data), 'supergroup'] <- inputObj_sub@meta.data$supergroup####
    Idents(inputObj) <- 'supergroup'
    print (levels(Idents(inputObj)))
    DimPlot(inputObj, reduction = 'tsne', pt.size = .5, group.by = 'supergroup', label = T, label.size = 2) + theme(legend.position = 'None')
    ggsave(paste(c(outputPath, '/res/tsne/tsne', res, '_.pdf'), collapse = ''), units = 'cm', width = 12, height = 12)
    DimPlot(inputObj, reduction = 'umap', pt.size = .5, group.by = 'supergroup', label = T, label.size = 2) + theme(legend.position = 'None')
    ggsave(paste(c(outputPath, '/res/umap/umap', res, '_.pdf'), collapse = ''), units = 'cm', width = 12, height = 12)
    clusternum <- length(unique(inputObj@meta.data$supergroup)) 
    
    if (clusternum != 1) {
      supergroupnum <- clusternum#length(unique(inputObj@meta.data$supergroup)) 
      
      clusterMarkers <- FindAllMarkers(object = inputObj, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
      heatmapgenes_all <- clusterMarkers %>% group_by(cluster) %>% top_n(10, avg_logFC)
      dehm <- DoHeatmap(object = inputObj, features = heatmapgenes_all$gene, angle = 90, size = 3, raster = T, draw.lines = F)
      ggsave(paste(c(outputPath, '/res/heatmap/heatmap_all', res, '_.pdf'), collapse = ''), units = 'cm', width = 30, height = 20)
      
      clusterMarkers_top <- subset(clusterMarkers, clusterMarkers$p_val_adj <= 0.05)
      clusterMarkers_top <- length(unique(clusterMarkers$gene)) 
      
      clusterMarkers_sig <- subset(clusterMarkers, clusterMarkers$p_val_adj <= 0.05 & avg_logFC >= 1)
      clusterMarkers_sig <- clusterMarkers_sig %>% group_by(cluster) %>% top_n(10, avg_logFC) ###
      clusterMarkers_sig <- length(unique(clusterMarkers_sig$gene)) 
      
    } else {
      clusternum <- 0
      supergroupnum <- 0
      clusterMarkers_top <- 0
      clusterMarkers_sig <- 0
    }
    
    ### most cell degs ### v3
    most_cell_clusters_num <- length(unique(Idents(inputObj_sub)))
    if (most_cell_clusters_num > 1) {
      mostMarkers <- FindAllMarkers(object = inputObj_sub, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
      heatmapgenes <- mostMarkers %>% group_by(cluster) %>% top_n(10, avg_logFC)
      dehm <- DoHeatmap(object = inputObj_sub, features = heatmapgenes$gene, angle = 90, size = 3, raster = T, draw.lines = F)
      ggsave(paste(c(outputPath, '/res/heatmap/heatmap_mostcell', res, '_.pdf'), collapse = ''), units = 'cm', width = 30, height = 20)
      
      most_cell_clusters_degs <- subset(mostMarkers, mostMarkers$p_val_adj <= 0.05)
      most_cell_clusters_degs <- length(unique(mostMarkers$gene)) 
      
      most_cell_clusters_degs_sig <- subset(mostMarkers, mostMarkers$p_val_adj <= 0.05 & avg_logFC >= 1)
      most_cell_clusters_degs_sig <- most_cell_clusters_degs_sig %>% group_by(cluster) %>% top_n(10, avg_logFC) ###
      most_cell_clusters_degs_sig <- length(unique(most_cell_clusters_degs_sig$gene)) 
      
    } else {
      most_cell_clusters_degs <- 0
      most_cell_clusters_degs_sig <- 0
    }
    ###
    
    if (nrow(clusterCount) == 0) {
      clusterCount[1, ] <- c(res, clusternum, supergroupnum, clusterMarkers_top, clusterMarkers_sig, most_cell_clusters_num, most_cell_clusters_degs, most_cell_clusters_degs_sig)
    } else {
      clusterCount <- rbind(clusterCount, res = c(res, clusternum, supergroupnum, clusterMarkers_top, clusterMarkers_sig, most_cell_clusters_num, most_cell_clusters_degs, most_cell_clusters_degs_sig))
    }
  }
  colnames(clusterCount) <- c('res', 'cluster', 'supergroup', 'deg', 'deg_sig_top10', 'mostcell', 'mostcell_deg', 'mostcell_deg_sig_top10')
  
  write.table(clusterCount, paste(c(outputPath, '/res/clusterCount_', 1/scaleFactor, '-', nrow(clusterCount)/scaleFactor, '_.txt'), collapse = ''), quote = F, sep = '\t', row.names = F, col.names = T)
  ggplot(clusterCount, aes(cluster)) + 
    geom_histogram(binwidth = 1) +
    theme_bw() + theme(panel.grid = element_blank())
  ggsave(paste(c(outputPath, '/res/clusterCount_', 1/scaleFactor, '-', nrow(clusterCount)/scaleFactor, '_.pdf'), collapse = ''), units = 'cm', width = 8, height = 8)
  ggplot(clusterCount, aes(cluster, deg)) + 
    geom_point() +
    theme_bw() + theme(panel.grid = element_blank())
  ggsave(paste(c(outputPath, '/res/deg_top10_', 1/scaleFactor, '-', nrow(clusterCount)/scaleFactor, '_.pdf'), collapse = ''), units = 'cm', width = 12, height = 6)
  ggplot(clusterCount, aes(cluster, deg_sig_top10)) + 
    geom_point() +
    theme_bw() + theme(panel.grid = element_blank())
  ggsave(paste(c(outputPath, '/res/deg_sig_top10_', 1/scaleFactor, '-', nrow(clusterCount)/scaleFactor, '_.pdf'), collapse = ''), units = 'cm', width = 12, height = 6)
  
}

######################################################################################################################################
library(Seurat)
library(ggplot2)
library(patchwork)
library(plyr)
library(dplyr)

### 4. Iterate resolution ###
lymphgland <- readRDS('tmp/lymphgland.Rds')
calc_res(lymphgland, './', pcNum = 32, scaleFactor = 10, rangeMax = 3, cormethod = 'pearson')
