#################
### Functions ###
#################

# 4. resolution fix
fix_res <- function(inputObj, outputPath, pcNum, res, cormethod = 'spearman') {
  dir.create(paste0(c(outputPath, '/res_fix_', res), collapse = ''))
  
  inputObj <- FindNeighbors(inputObj, dims = 1:pcNum)
  inputObj <- FindClusters(inputObj, resolution = res)
  
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
  head(emptydf)
  
  
  corMat <- data.frame(cor(emptydf, method = cormethod), check.names = F) # pearson spearman
  pseudo_exprs <- emptydf; head(pseudo_exprs)
  inputObj@meta.data$supergroup <- as.character(inputObj@meta.data$seurat_clusters)
  
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
      head(inputObj_sub_pseudo)
      
      for (tmpcluster_idx in c(1:ncol(inputObj_sub_pseudo))) {
        tmp_obj <- subset(inputObj_sub, idents = levels(inputObj_sub@meta.data$seurat_clusters)[tmpcluster_idx])
        tmp_cells <- rownames(tmp_obj@meta.data)
        tmp_pseudo <- rowSums(vardf[, tmp_cells])
        inputObj_sub_pseudo[,tmpcluster_idx] <- tmp_pseudo
      }
      
      colnames(inputObj_sub_pseudo) <- levels(inputObj_sub@meta.data$seurat_clusters)
      head(inputObj_sub_pseudo)
      
      corMat <- data.frame(cor(inputObj_sub_pseudo, method = cormethod), check.names = F)
      round <- round + 1
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
    print (maxpair)
    
    inputObj_sub@meta.data[rownames(subset(inputObj_sub@meta.data, supergroup == as.character(maxpair$v1))), 'supergroup'] <- paste0(c(as.character(maxpair$v1), as.character(maxpair$v2)), collapse = '_')
    inputObj_sub@meta.data[rownames(subset(inputObj_sub@meta.data, supergroup == as.character(maxpair$v2))), 'supergroup'] <- paste0(c(as.character(maxpair$v1), as.character(maxpair$v2)), collapse = '_')
    
    inputObj_sub_pseudo[, paste0(c(as.character(maxpair$v1), as.character(maxpair$v2)), collapse = '_')] <- inputObj_sub_pseudo[, as.character(maxpair$v1)] + inputObj_sub_pseudo[, as.character(maxpair$v2)]
    inputObj_sub_pseudo[, as.character(maxpair$v1)] <- NULL
    inputObj_sub_pseudo[, as.character(maxpair$v2)] <- NULL
    
    corMat <- cor(inputObj_sub_pseudo, method = cormethod)
    diag(corMat) <- 0
    print (head(corMat))
  }
  
  Idents(inputObj_sub) <- 'supergroup'
  
  inputObj@meta.data[rownames(inputObj_sub@meta.data), 'supergroup'] <- inputObj_sub@meta.data$supergroup
  Idents(inputObj) <- 'supergroup'
  print (levels(Idents(inputObj)))
  DimPlot(inputObj, reduction = 'tsne', pt.size = .5, group.by = 'supergroup', label = T, label.size = 2) + theme(legend.position = 'None')
  ggsave(paste(c(outputPath, '/res_fix_', res, '/tsne_', res, '.pdf'), collapse = ''), units = 'cm', width = 12, height = 12)
  DimPlot(inputObj, reduction = 'umap', pt.size = .5, group.by = 'supergroup', label = T, label.size = 2) + theme(legend.position = 'None')
  ggsave(paste(c(outputPath, '/res_fix_', res, '/umap_', res, '.pdf'), collapse = ''), units = 'cm', width = 12, height = 12)
  
  ###
  clusterMarkers <- FindAllMarkers(object = inputObj, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
  heatmapgenes_all <- clusterMarkers %>% group_by(cluster) %>% top_n(10, avg_logFC)
  dehm <- DoHeatmap(object = inputObj, features = heatmapgenes_all$gene, angle = 90, size = 3, raster = T, draw.lines = F)
  ggsave(paste(c(outputPath, '/res_fix_', res, '/heatmap_all_', res, '.pdf'), collapse = ''), units = 'cm', width = 30, height = 20)
  write.table(clusterMarkers, paste(c(outputPath, '/res_fix_', res, '/heatmap_all_', res, '.txt'), collapse = ''), quote = F, sep = '\t', row.names = F, col.names = T)
  
  if (length(unique(Idents(inputObj_sub))) > 1) {
    mostMarkers <- FindAllMarkers(object = inputObj_sub, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
    heatmapgenes <- mostMarkers %>% group_by(cluster) %>% top_n(10, avg_logFC)
    dehm <- DoHeatmap(object = inputObj_sub, features = heatmapgenes$gene, angle = 90, size = 3, raster = T, draw.lines = F)
    ggsave(paste(c(outputPath, '/res_fix_', res, '/heatmap_mostcell_', res, '.pdf'), collapse = ''), units = 'cm', width = 30, height = 20)
    write.table(mostMarkers, paste(c(outputPath, '/res_fix_', res, '/heatmap_mostcell_', res, '.txt'), collapse = ''), quote = F, sep = '\t', row.names = F, col.names = T)
  }
  writetable <- data.frame(inputObj@meta.data, check.rows = F, check.names = F)
  writetable <- data.frame(Barcode = rownames(writetable), writetable, check.rows = F, check.names = F)
  write.table(writetable, paste(c(outputPath, '/res_fix_', res, '/label_', res, '.txt'), collapse = ''), quote = F, sep = '\t', row.names = F, col.names = T)
  
  return(inputObj)
}

######################################################################################################################################
library(Seurat)
library(ggplot2)
library(patchwork)
library(plyr)
library(dplyr)

### 5. Fix resolution value ###
lymphgland <- readRDS('tmp/lymphgland.Rds')
lymphgland <- ScaleData(object = lymphgland, features = rownames(lymphgland), vars.to.regress = c('Library', 'nCount_RNA'))
lymphgland <- fix_res(lymphgland, './', pcNum = 32, res = 0.3, cormethod = 'pearson') # new
saveRDS(lymphgland, 'tmp/lymphglandt.fix.res0.3.Rds')


### 5. Fix resolution value - res 0.6 ###
lymphgland <- readRDS('tmp/lymphglandt.fix.res0.3.Rds')
lymphgland <- fix_res(lymphgland, './', pcNum = 32, res = 0.6, cormethod = 'pearson') # new
saveRDS(lymphgland, 'tmp/lymphglandt.fix.res0.6.Rds')
