### Plasmatocytes ###
library(Seurat)
library(ggplot2)
library(patchwork)
library(plyr)
library(dplyr)

#dir.create('tmp')
scGenes <- read.delim('__filter_scGenes__/bulk_sc_pseudo_pt3-scGenes.scGenes_v2.txt', header = F)
scGenes <- scGenes$V1
bkGenes <- read.delim('__filter_scGenes__/bulk_sc_pseudo_pt3-scGenes.bkGenes_v2.txt', header = F)
bkGenes <- bkGenes$V1
filtergenes <- c(as.character(scGenes), as.character(bkGenes))
length(filtergenes)
#saveRDS(filtergenes, 'tmp/filtergenes.Rds')
#filtergenes <- readRDS('tmp/filtergenes.Rds')

lymphgland.combined <- readRDS('../../../Drop-seq_alignment/1.seurat3_alignment_withMuscle_regress-Library-nUMI/lymphgland.combined.flt_allscaled.Rds')
head(lymphgland.combined@meta.data)

Idents(lymphgland.combined) <- 'anno_simple'
pm <- subset(lymphgland.combined, idents = 'PM')
remove(lymphgland.combined)


### New variable genes ###
pm <- FindVariableFeatures(pm, selection.method = "vst", nfeatures = 2000)
#write.table(intersect(VariableFeatures(pm), filtergenes), 'filteredgenes_PM.txt', quote = F, sep = '\n', row.names = F, col.names = F)

DefaultAssay(pm)
head(pm@assays$RNA@meta.features)
pm@assays$RNA@meta.features$vst.variable <- FALSE
testit <- pm@assays$RNA@meta.features[order(pm@assays$RNA@meta.features$vst.variance.standardized, decreasing = T),]
newvargenes <- setdiff(rownames(testit), filtergenes)[1:2000]
pm@assays$RNA@var.features <- newvargenes
pm@assays$RNA@meta.features[newvargenes, 'vst.variable'] <- TRUE

top10 <- head(VariableFeatures(pm), 10)
plot1 <- VariableFeaturePlot(pm)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
dir.create('stats')
ggsave('stats/2-3.new_variableFeatures.pdf', units = 'cm', width = 17, height = 10)


###
#saveRDS(pm, 'tmp/pm.Rds')
#pm <- readRDS('tmp/pm.Rds')

calc_pcs(pm, './') # 31 PCs
calc_lat(pm, './', 31, 200330101) # tsne 200330115 umap 200330160
pm <- RunTSNE(pm, dims = 1:31, reduction.key = 'tSNE', dim.embed = 3, seed.use = 200330115)
pm <- RunUMAP(pm, reduction = "pca", dims = 1:31, min.dist = 0.3, n.components = 3, seed.use = 200330160)

calc_res(pm, './', pcNum = 31, scaleFactor = 10, rangeMax = 3, cormethod = 'pearson')

pm <- fix_res(pm, './', pcNum = 31, res = 0.7, cormethod = 'pearson') # old
pm <- fix_res(pm, './', pcNum = 31, res = 0.6, cormethod = 'pearson') # new
saveRDS(pm, 'tmp/pm.Rds')
summary(pm@meta.data)


#################
### Functions ###
#################
# 1. PCA
calc_pcs <- function(inputObj, outputPath) {
  inputObj <- ScaleData(object = inputObj, vars.to.regress = c('Library', 'nCount_RNA'))
  inputObj <- RunPCA(object = inputObj, npcs = 40)
  inputObj <- JackStraw(object = inputObj, num.replicate = 100, dims = 40)
  inputObj <- ScoreJackStraw(object = inputObj, dims = 1:40)
  
  dir.create(paste0(c(outputPath, 'stats'), collapse = '/'))
  JackStrawPlot(inputObj, dims = 1:40) #
  ggsave(paste0(c(outputPath, 'stats/2-4.pca.JackStrawPlot.mitoCut.pdf'), collapse = '/'), units = 'cm', width = 20, height = 12)
  ElbowPlot(inputObj, ndims = 40) #
  ggsave(paste0(c(outputPath, 'stats/2-4.pca.ElbowPlot.mitoCut.pdf'), collapse = '/'), units = 'cm', width = 15, height = 10)
}
# 2. tSNE & UMAP
calc_lat <- function(inputObj, outputPath, pcNum, inputseed, iter = 20) {
  dir.create(paste0(c(outputPath, 'tsne'), collapse = '/'))
  dir.create(paste0(c(outputPath, 'tsne/compare'), collapse = '/'))
  dir.create(paste0(c(outputPath, 'umap'), collapse = '/'))
  dir.create(paste0(c(outputPath, 'umap/compare'), collapse = '/'))
  
  seedstart <- inputseed
  seedend <- inputseed + iter
  for (tmpseed in c(seedstart:seedend)){
    inputObj <- RunTSNE(object = inputObj, dims = 1:pcNum, reduction.key = 'tSNE', dim.embed = 3, seed.use = tmpseed)
    FeaturePlot(inputObj, features = c("Nplp2", "Hml", "Pxn", "Ama", "NimC1", "eater", "Jabba", "stg", "AttA"), cols = c("grey","red"), reduction = "tsne")
    ggsave(paste(c(outputPath, '/tsne/compare/tsne.mitoCut.1_2.markers.seed', tmpseed, '.pdf'), collapse = ''), units = 'cm', width = 27, height = 22)
    }
  seedstart <- inputseed +50
  seedend <- inputseed +50 + iter
  for (tmpseed in c(seedstart:seedend)){
    inputObj <- RunUMAP(inputObj, reduction = "pca", dims = 1:pcNum, min.dist = 0.3, n.components = 3, seed.use = tmpseed)
    FeaturePlot(inputObj, features = c("Nplp2", "Hml", "Pxn", "Ama", "NimC1", "eater", "Jabba", "stg", "AttA"), cols = c("grey","red"))
    ggsave(paste(c(outputPath, '/umap/compare/umap.mitoCut.markers.mindist_0.3.seed', tmpseed, '.pdf'), collapse = ''), units = 'cm', width = 27, height = 22)
    }
}
# 3. res v3 - excluding least correlated cells 
calc_res <- function(inputObj, outputPath, pcNum, scaleFactor, rangeMax = 2, cormethod = 'spearman') {
  dir.create(paste0(c(outputPath, 'res'), collapse = '/'))
  dir.create(paste0(c(outputPath, 'res/tsne'), collapse = '/'))
  dir.create(paste0(c(outputPath, 'res/umap'), collapse = '/'))
  dir.create(paste0(c(outputPath, 'res/heatmap'), collapse = '/'))
  
  ### testor ###
  #inputObj <- ph
  #outputPath <- './'
  #pcNum <- 17
  #scaleFactor <- 10
  #rangeMax = 3
  #cormethod = 'pearson'
  #res <- 3
  
  ### testor ###
  
  inputObj <- FindNeighbors(inputObj, dims = 1:pcNum)
  clusterCount <- data.frame(matrix(ncol = 8, nrow = 0))
  
  endpoint <- rangeMax*scaleFactor
  for (res in c(1:endpoint)){
    res <- res/scaleFactor
    inputObj <- FindClusters(inputObj, resolution = res)
    clusternum <- length(unique(inputObj@meta.data$seurat_clusters))
    print (res)
    #####################
    if (clusternum != 1) {
      vardf <- data.frame(as.matrix(GetAssayData(inputObj, slot = 'data')), check.rows = F, check.names = F)
      vardf <- vardf[VariableFeatures(inputObj), ]
      dim(vardf)
      
      emptydf <- data.frame(row.names = VariableFeatures(inputObj), matrix(nrow = 2000, ncol = length(unique(inputObj@meta.data$seurat_clusters))))
      head(emptydf)
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
# 4. resolution fix
fix_res <- function(inputObj, outputPath, pcNum, res, cormethod = 'spearman') {
  dir.create(paste0(c(outputPath, '/res_fix_', res), collapse = ''))
  
  ### test ###
  inputObj <- ph
  outputPath <- './'
  pcNum = 17
  res = 1.3
  cormethod = 'pearson'
  ###
  
  
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



######################
### Filtered genes ###
######################
filtergenes <- readRDS('tmp/filtergenes.Rds')

pm_oldvari <- FindVariableFeatures(pm, selection.method = "vst", nfeatures = 2000)
Idents(pm_oldvari) <- 'Subclustering'

foundfiltergenes <- intersect(as.character(VariableFeatures(pm_oldvari)), as.character(filtergenes))[1:100]; length(foundfiltergenes)
DoHeatmap(object = pm_oldvari, features = foundfiltergenes, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('filteredgenes_PM_vari_top100.pdf', units = 'cm', width = 30, height = 20)

testdata <- data.frame(as.matrix(GetAssayData(pm_oldvari, slot = 'data')), check.rows = F, check.names = F)
testdata$pseudo <- rowSums(testdata)
testdata <- subset(testdata, pseudo > 1)
testdata <- testdata[order(testdata$pseudo, decreasing = T), ]
foundfiltergenes <- intersect(rownames(testdata), as.character(filtergenes))[1:100]; length(foundfiltergenes)
DoHeatmap(object = pm_oldvari, features = foundfiltergenes, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('filteredgenes_PM_top100.pdf', units = 'cm', width = 30, height = 20)



#############
### Stats ###
#############
label <- read.delim('res_fix_0.7_dump/label_0.7.txt', row.names = 1)
label <- read.delim('res_fix_0.6/label_0.6.txt', row.names = 1)
head(label)
label$Subclustering <- factor(label$Subclustering, levels = c("PM 1", "PM 2", "PM 3", "PM 4", "PM 5", "PM 6", "PM 7", "PM 8", "PM 9", "PM 10"))
levels(label$Subclustering)

count_subcluster <- data.frame(row.names = levels(label$Subclustering), matrix(nrow = length(levels(label$Subclustering)), ncol = length(unique(label$supergroup))))
colnames(count_subcluster) <- unique(label$supergroup)
for (newcluster in unique(label$supergroup)) {
  count_subcluster[, newcluster] <- summary(subset(label, supergroup == newcluster)$Subclustering)
}
count_subcluster <- data.frame(Subcluster = rownames(count_subcluster), orig_count = summary(label$Subclustering), count_subcluster, check.names = F, check.rows = F)
write.table(count_subcluster, 'res_fix_0.6/label_0.6_to_subclustering.txt', quote = F, sep = '\t', row.names = F, col.names = T)


