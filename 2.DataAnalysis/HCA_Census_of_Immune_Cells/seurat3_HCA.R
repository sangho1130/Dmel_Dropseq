#!/home/sangho/miniconda3/envs/seurat3/bin/R

library(Seurat)
library(ggplot2)
library(Matrix)
library(cowplot)
library(plyr)
library(dplyr)
setwd('/home/sangho/YSH/Project/drosophila-bloodcells/other_dataset/HCA/Census_of_Immune_Cells/seurat3_v2-2')

ica_bone_marrow <- Read10X_h5('../ica_bone_marrow_h5.h5')
ica_bm <- CreateSeuratObject(counts = ica_bone_marrow)

label <- read.delim('/home/sangho/YSH/Project/drosophila-bloodcells/other_dataset/HCA/Census_of_Immune_Cells/barcode_label_v2.txt', row.names = 1, sep = '\t')
ica_bm <- AddMetaData(ica_bm, label)
head(ica_bm@meta.data)

plt <- ggplot(ica_bm@meta.data, aes(Patient, nCount_RNA)) + geom_jitter(size = 0.25) + 
  geom_hline(yintercept = c(70000, 80000), col = c('red2', 'blue')) + 
  theme_bw() + theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank());plt
ggsave('stats/2-1.stats.nCount_RNA.byPatient.pdf', units = 'cm', width = 8, height = 6)
plt <- ggplot(ica_bm@meta.data, aes(Patient, nCount_RNA)) + geom_jitter(size = 0.25) + 
  ylim(0, 5000) + theme_bw() + theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank());plt
ggsave('stats/2-1.stats.nCount_RNA_low.byPatient.pdf', units = 'cm', width = 8, height = 6)
plt <- ggplot(ica_bm@meta.data, aes(Patient, nFeature_RNA)) + geom_jitter(size = 0.25) + 
  theme_bw() + theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank());plt
ggsave(paste(c('stats/2-1.stats.nFeature_RNA.byPatient.pdf'), collapse = ''), units = 'cm', width = 8, height = 6)


ica_bm_list <- SplitObject(ica_bm, split.by = "Patient")
ica_bm_list
for (tmpident in c('BM1', 'BM2', 'BM3', 'BM4', 'BM5', 'BM6', 'BM7', 'BM8')){
  tmplib <- eval(parse(text=paste("ica_bm_list$`", tmpident, '`', sep = "")))
  if (tmpident == 'BM4'){
    tmplib <- subset(tmplib, subset = nCount_RNA < 70000 & nFeature_RNA >= 500)
  } else{
    tmplib <- subset(tmplib, subset = nCount_RNA < 80000 & nFeature_RNA >= 500)
  }
  umihighcut <- as.integer(mean(tmplib@meta.data$nCount_RNA) + 2*sd(tmplib@meta.data$nCount_RNA))
  tmplib <- subset(tmplib, subset = nCount_RNA < umihighcut)
  ica_bm_list[[tmpident]] <- tmplib
}
ncol(ica_bm_list[['BM1']]@assays$RNA) + ncol(ica_bm_list[['BM2']]@assays$RNA) + ncol(ica_bm_list[['BM3']]@assays$RNA) + ncol(ica_bm_list[['BM4']]@assays$RNA) + 
  ncol(ica_bm_list[['BM5']]@assays$RNA) + ncol(ica_bm_list[['BM6']]@assays$RNA) + ncol(ica_bm_list[['BM7']]@assays$RNA) + ncol(ica_bm_list[['BM8']]@assays$RNA) #267598


objList <- list()
for (tmpident in c('BM1', 'BM2', 'BM3', 'BM4', 'BM5', 'BM6', 'BM7', 'BM8')){
  tmplib <- eval(parse(text=paste("ica_bm_list$`", tmpident, '`', sep = "")))
  tmpExpr <- as(as.matrix(GetAssayData(tmplib)), "dgCMatrix")
  tmpLabel <- tmplib@meta.data
  
  tmpObj <- CreateSeuratObject(counts = tmpExpr, project = tmpident)
  tmpObj <- AddMetaData(object = tmpObj, metadata = tmpLabel)
  tmpObj <- NormalizeData(tmpObj, normalization.method = "LogNormalize", scale.factor = 10000)
  tmpObj <- FindVariableFeatures(object = tmpObj, selection.method = "vst", nfeatures = 2000)
  tmpObj@meta.data$Patient <- tmpident
  tmpObj[["percent.mt"]] <- PercentageFeatureSet(object = tmpObj, pattern = "^MT-")
  
  objList[[tmpident]] <- tmpObj
}

ica_bm.anchors <- FindIntegrationAnchors(object.list = objList, dims = 1:30)
ica_bm.combined <- IntegrateData(anchorset = ica_bm.anchors, dims = 1:30)
DefaultAssay(object = ica_bm.combined) <- "integrated"
ica_bm.combined@meta.data$Patient <- factor(ica_bm.combined@meta.data$Patient, c('BM1', 'BM2', 'BM3', 'BM4', 'BM5', 'BM6', 'BM7', 'BM8'))

### Filter high MT ###
ica_bm.combined <- subset(ica_bm.combined, subset = percent.mt <= 10)
print (nrow(ica_bm.combined@meta.data))
print ('cells detected')

plt <- ggplot(ica_bm.combined@meta.data, aes(Patient, nCount_RNA)) + geom_jitter(size = 0.25) + 
  theme_bw() + theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank());plt
ggsave('stats/2-2.stats.nCount_RNA.byPatient.pdf', units = 'cm', width = 8, height = 6)
plt <- ggplot(ica_bm.combined@meta.data, aes(Patient, nFeature_RNA)) + geom_jitter(size = 0.25) + 
  theme_bw() + theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank());plt
ggsave(paste(c('stats/2-2.stats.nFeature_RNA.byPatient.pdf'), collapse = ''), units = 'cm', width = 8, height = 6)

plt <- ggplot(ica_bm.combined@meta.data, aes(Patient, percent.mt)) + geom_jitter(size = 0.25) + 
  theme_bw() + theme(text = element_text(size = 7), axis.text = element_text(size = 7), panel.grid = element_blank());plt
ggsave('stats/2-2.stats.percent.mt.byPatient.pdf', units = 'cm', width = 8, height = 6)
plot1 <- FeatureScatter(object = ica_bm.combined, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'Patient') + geom_abline(intercept = 10, col = 'red2', slope = 0)
plot2 <- FeatureScatter(object = ica_bm.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'Patient')
CombinePlots(plots = list(plot1, plot2))
ggsave('stats/2-2.stats.percent.mt_nFeature_bynCount_RNA.byPatient.pdf', units = 'cm', width = 26, height = 10)

writeLabel <- as.matrix(ica_bm.combined@meta.data); writeLabel <- data.frame(Barcode = rownames(writeLabel), writeLabel, check.names = F)
write.table(writeLabel, 'tmp1.label.txt', sep = '\t', quote = F, row.names = F, col.names = T)


### Standard ###
ica_bm.combined <- ScaleData(object = ica_bm.combined, vars.to.regress = c('Patient', 'nCount_RNA'))
ica_bm.combined <- RunPCA(object = ica_bm.combined, npcs = 80)
ica_bm.combined <- JackStraw(object = ica_bm.combined, num.replicate = 100, dims = 80)
ica_bm.combined <- ScoreJackStraw(object = ica_bm.combined, dims = 1:80)

JackStrawPlot(ica_bm.combined, dims = 1:80) # 61 PCs
ggsave('stats/2-4.pca.JackStrawPlot.pdf', units = 'cm', width = 30, height = 12)
ElbowPlot(ica_bm.combined, ndims = 80) # 61 PCs
ggsave('stats/2-4.pca.ElbowPlot.pdf', units = 'cm', width = 15, height = 10)
save.image('tmp2.Rdata')

for (tmpseed in c(1017151:1017170)){
  ica_bm.combined <- RunUMAP(ica_bm.combined, reduction = "pca", dims = 1:61, min.dist = 0.3, n.components = 3, seed.use = tmpseed, umap.method = 'umap-learn', metric = 'correlation')
  for (testgene in c("CD34", "CD8A", "NKG7", "CD14")) {
    plt <- FeaturePlot(ica_bm.combined, features = testgene, cols = c("grey", "red"), pt.size = 0.5)
    AugmentPlot(plt, dpi = 300, width = 12, height = 10)
    ggsave(paste(c('umap/compare/umap.', testgene, '.mindist_0.3.seed', tmpseed, '.pdf'), collapse = ''), units = 'cm', width = 12, height = 10)
  }
} # 

save.image('tmp3.Rdata')

# 1017164
ica_bm.combined <- RunUMAP(ica_bm.combined, reduction = "pca", dims = 1:61, min.dist = 0.3, n.components = 3, seed.use = 1017164, umap.method = 'umap-learn', metric = 'correlation')


for (tsetgene in c("CD34", "CD4", "CD8A", "NKG7", "CD79B", "CD14", "B3GAT1", "CD69", "LYZ", "HLA-DRA", "FCN1", "NCAM1", "FCGR3A")){
  plt <- FeaturePlot(ica_bm.combined, features = tsetgene, cols = c("grey", "red"), pt.size = 0.5)
  AugmentPlot(plt, dpi = 300, width = 12, height = 10)
  ggsave(paste(c('umap/umap.markers-', tsetgene, '.mindist_0.3.seed1017164.pdf'), collapse = ''), units = 'cm', width = 10.5, height = 10)
}


### res ###
ica_bm.combined <- FindNeighbors(object = ica_bm.combined, dims = 1:61)
for (res in c(5:12)){
  res <- res/10
  ica_bm.combined <- FindClusters(object = ica_bm.combined, resolution = res)
  plt <- DimPlot(ica_bm.combined, label = T, pt.size = 0.5, reduction = 'umap')
  AugmentPlot(plt, dpi = 300, width = 12, height = 12)
  ggsave(paste(c('res/compare/umap.seed.1_2.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 12, height = 12)
  plt <- DimPlot(ica_bm.combined, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'umap')
  AugmentPlot(plt, dpi = 300, width = 12, height = 12)
  ggsave(paste(c('res/compare/umap.seed.1_3.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 12, height = 12)
  plt <- DimPlot(ica_bm.combined, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'umap')
  AugmentPlot(plt, dpi = 300, width = 12, height = 12)
  ggsave(paste(c('res/compare/umap.seed.2_3.res_',res,'.pdf'), collapse = ''), units = 'cm', width = 12, height = 12)
}
save.image('tmp4.Rdata')


ica_bm.combined <- FindClusters(object = ica_bm.combined, resolution = 0.8)
plt <- DimPlot(ica_bm.combined, label = T, pt.size = 0.5, reduction = 'umap')
AugmentPlot(plt, dpi = 300, width = 12, height = 12)
ggsave('res/umap.seed1017164.1_2.res_0.6.pdf', units = 'cm', width = 12, height = 12)
plt <- DimPlot(ica_bm.combined, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'umap')
AugmentPlot(plt, dpi = 300, width = 12, height = 12)
ggsave('res/umap.seed1017164.1_3.res_0.6.pdf', units = 'cm', width = 12, height = 12)
plt <- DimPlot(ica_bm.combined, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'umap')
AugmentPlot(plt, dpi = 300, width = 12, height = 12)
ggsave('res/umap.seed1017164.2_3.res_0.6.pdf', units = 'cm', width = 12, height = 12)

ica_bm.combined@meta.data <- ica_bm.combined@meta.data[, c(1:5, 7, 10)]
ica_bm.combined@meta.data$use_anno <- ica_bm.combined@active.ident
writeLable <- as.matrix(ica_bm.combined@meta.data)
writeLable <- data.frame(Barcode = rownames(writeLable), writeLable, check.names = F)
write.table(writeLable, 'tmp4.label.txt', sep = '\t', quote = F, row.names = F, col.names = T)


ica_bm.combined.markers <- FindAllMarkers(object = ica_bm.combined, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
ica_bm.combined.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(ica_bm.combined.markers, 'degs/findAllMarkers.pct25.thre0.25.onlyPos.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- ica_bm.combined.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

save.image('tmp5.Rdata')

DoHeatmap(object = ica_bm.combined, features = top10$gene, angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('degs/degs.DEgenes.top10.pdf', units = 'cm', width = 50, height = 40, limitsize=F)


seurat_clusters_anno_marker <- c('Naive T 1', 'Monocyte 1', 'NK 1', 'Naive T 5', 'B cell 1', 'Naive T 3', 'NK 2',
                                 'Naive CD8 1', 'B cell 2', 'Eos/Neu', 'Naive T 6', 'Cytotoxic 1', 'Cytotoxic 2', 'Pro-B 2',
                                 'Naive CD8 2', 'Granulocyte progenitors', 'Naive T 2', 'B cell 3', 'Erythroblast 2', 'CD8 1', 'cDC',
                                 'B precursor', 'HSC~MPP', 'NK 3', 'Monocyte 2', 'pDC', 'Naive T 4', 'CD8 2',
                                 'Erythroblast 1', 'Plasma cell', 'Pro-B 1', 'T small', 'Stromal', 'Platelet')
names(x = seurat_clusters_anno_marker) <- levels(x = ica_bm.combined)
ica_bm.combined <- RenameIdents(ica_bm.combined, seurat_clusters_anno_marker)

ica_bm.combined@active.ident <- factor(ica_bm.combined@active.ident, 
                                       levels = c('HSC~MPP', 'Erythroblast 1', 'Erythroblast 2', 'Platelet',
                                                  'Granulocyte progenitors', 'Eos/Neu', 'Monocyte 1', 'Monocyte 2', 'cDC', 'pDC',
                                                  'B precursor', 'Pro-B 1', 'Pro-B 2', 'B cell 1', 'B cell 2', 'B cell 3', 'Plasma cell',
                                                  'T small', 'Naive T 1', 'Naive T 2', 'Naive T 3', 'Naive T 4', 'Naive T 5', 'Naive T 6',
                                                  'Naive CD8 1', 'Naive CD8 2', 'CD8 1', 'CD8 2', 'Cytotoxic 1', 'Cytotoxic 2', 
                                                  'NK 1', 'NK 2', 'NK 3', 
                                                  'Stromal'))


ica_bm.combined.markers <- FindAllMarkers(object = ica_bm.combined, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
ica_bm.combined.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(ica_bm.combined.markers, 'degs/findAllMarkers.pct25.thre0.25.onlyPos.txt', sep = '\t', quote = F, col.names = T, row.names = F) 
top10 <- ica_bm.combined.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

ica_bm.combined@meta.data$use_anno <- ica_bm.combined@active.ident
writeLable <- as.matrix(ica_bm.combined@meta.data)
writeLable <- data.frame(Barcode = rownames(writeLable), writeLable, check.names = F)
write.table(writeLable, 'tmp6.label.txt', sep = '\t', quote = F, row.names = F, col.names = T)
save.image('tmp6.Rdata')

### new figures ###
DoHeatmap(object = ica_bm.combined, features = unique(top10$gene), angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('degs/degs.DEgenes.top10.pdf', units = 'cm', width = 50, height = 40, limitsize = F)

DotPlot(object = ica_bm.combined, features = unique(top10$gene), dot.scale = 3)
ggsave('degs/dotplt.DEgenes.top10.pdf', units = 'cm', width = 50, height = 40, limitsize = F)
# FCGR3A for NK and CD16+ monocyte; NCAM1 for NK; B3GAT1 for effecter; CD69 for memory; 
DotPlot(object = ica_bm.combined, features = c('CD34', 'FCGR3A', 'NCAM1', 'B3GAT1', 'CD69', 'CD8A', 'CD4'), dot.scale = 6, assay='RNA') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('degs/dotplt.markers1.pdf', units = 'cm', width = 20, height = 30)


plt <- DimPlot(ica_bm.combined, label = T, pt.size = 0.5, reduction = 'umap')
AugmentPlot(plt, dpi = 300, width = 12, height = 12)
ggsave('res/umap.seed1017164.1_2.res_0.6.pdf', units = 'cm', width = 12, height = 12)
plt <- DimPlot(ica_bm.combined, label = T, dims = c(1,3), pt.size = 0.5, reduction = 'umap')
AugmentPlot(plt, dpi = 300, width = 12, height = 12)
ggsave('res/umap.seed1017164.1_3.res_0.6.pdf', units = 'cm', width = 12, height = 12)
plt <- DimPlot(ica_bm.combined, label = T, dims = c(2,3), pt.size = 0.5, reduction = 'umap')
AugmentPlot(plt, dpi = 300, width = 12, height = 12)
ggsave('res/umap.seed1017164.2_3.res_0.6.pdf', units = 'cm', width = 12, height = 12)


### save temp files ###
dir.create('tmp')
saveRDS(ica_bone_marrow, 'tmp/ica_bone_marrow.Rds')
saveRDS(label, 'tmp/label.Rds')
saveRDS(ica_bm, 'tmp/ica_bm.Rds')
saveRDS(ica_bm.anchors, 'tmp/ica_bm.anchors.Rds')
saveRDS(ica_bm.combined, 'tmp/ica_bm.combined.Rds')
saveRDS(ica_bm.combined.markers, 'tmp/ica_bm.combined.markers.Rds')



DefaultAssay(ica_bm.combined) <- 'RNA'

### similarity ###
pseudo_exprs <- data.frame(row.names = rownames(ica_bm.combined))
for (celltype in unique(ica_bm.combined@meta.data$use_anno)) {
  tmp_obj <- subset(ica_bm.combined, idents = celltype)
  tmp_expr <- as.matrix(GetAssayData(tmp_obj, slot = 'data', assay = 'RNA'))
  
  psuedo_expr <- data.frame(row.names = rownames(tmp_expr), sum = rowSums(tmp_expr))
  if (identical(rownames(psuedo_expr), rownames(pseudo_exprs))){
    pseudo_exprs <- cbind(pseudo_exprs, psuedo_expr)
  }
}

colnames(pseudo_exprs) <- unique(ica_bm.combined@meta.data$use_anno)
saveRDS(pseudo_exprs, 'tmp/pseudo_exprs.Rds')
head(pseudo_exprs)

dir.create('correlation')
print (dim(pseudo_exprs)) ##
pseudo_exprs <- pseudo_exprs[rownames(subset(pseudo_exprs, rowSums(pseudo_exprs) != 0)), ] ##
print (dim(pseudo_exprs)) ##
corMat <- cor(pseudo_exprs, method = 'spearman')
saveRDS(corMat, 'tmp/corMat.Rds')
head(corMat); min(corMat)
write.table(data.frame(Celltype = rownames(corMat), corMat, check.rows = F, check.names = F), 'correlation/SpearmanCorr.allGenes.txt', quote = F, sep = '\t', row.names = F, col.names = T)

###
library(pheatmap)

callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

pheatmap::pheatmap(corMat, color = colorRampPalette(c("steelblue", 'orange', "red2"))(n = 30),
                   breaks = c(70:100)/100,
                   fontsize_row = 9, fontsize_col = 9, border_color = NA, cellwidth = 8, cellheight = 8,
                   clustering_distance_rows = 'euclidean', clustering_method = 'complete', clustering_callback = callback, 
                   main = 'HCA Census', filename = 'correlation/SpearmanCorr.allGenes.pdf')

pheatmap::pheatmap(corMat, color = colorRampPalette(c("steelblue", 'orange', "red2", 'red4'))(n = 4),
                   breaks = c(0.8, 0.85, 0.9, 0.95, 1.0),
                   fontsize_row = 9, fontsize_col = 9, border_color = NA, cellwidth = 8, cellheight = 8,
                   clustering_distance_rows = 'euclidean', clustering_method = 'complete', clustering_callback = callback,
                   main = 'HCA Census', filename = 'correlation/SpearmanCorr.allGenes.categorized.pdf')



###
ica_bm.combined <- readRDS('tmp/ica_bm.combined.Rds')
DefaultAssay(ica_bm.combined) <- 'RNA'
ica_bm.combined@meta.data$supergroup <- as.character(ica_bm.combined@meta.data$use_anno)

corMat <- readRDS('tmp/corMat.Rds')
diag(corMat) <- 0

pseudo_exprs <- readRDS('tmp/pseudo_exprs.Rds')
pseudo_exprs <- pseudo_exprs[rownames(subset(pseudo_exprs, rowSums(pseudo_exprs) != 0)), ] ##
print (dim(pseudo_exprs)) ##

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
  
  ica_bm.combined@meta.data[rownames(subset(ica_bm.combined@meta.data, supergroup == as.character(maxpair$v1))), 'supergroup'] <- paste0(c(as.character(maxpair$v1), as.character(maxpair$v2)), collapse = '_')
  ica_bm.combined@meta.data[rownames(subset(ica_bm.combined@meta.data, supergroup == as.character(maxpair$v2))), 'supergroup'] <- paste0(c(as.character(maxpair$v1), as.character(maxpair$v2)), collapse = '_')
  
  pseudo_exprs[, paste0(c(as.character(maxpair$v1), as.character(maxpair$v2)), collapse = '_')] <- pseudo_exprs[, as.character(maxpair$v1)] + pseudo_exprs[, as.character(maxpair$v2)]
  pseudo_exprs[, as.character(maxpair$v1)] <- NULL
  pseudo_exprs[, as.character(maxpair$v2)] <- NULL
  
  corMat <- cor(pseudo_exprs, method = 'spearman')
  diag(corMat) <- 0
  print (head(corMat))
}

print (unique(ica_bm.combined@meta.data$supergroup))
saveRDS(ica_bm.combined, 'tmp/ica_bm.combined.supergroup.Rds')
saveRDS(corMat, 'tmp/corMat.supergroup.Rds')
write.table(data.frame(Celltype = rownames(corMat), corMat, check.rows = F, check.names = F), 'correlation/SpearmanCorr.allGenes.supergroup.txt', quote = F, sep = '\t', row.names = F, col.names = T)

###
library(pheatmap)

callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

pheatmap::pheatmap(corMat, color = colorRampPalette(c("steelblue", 'grey60', "red2"))(n = 20),
                   breaks = c(0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0),
                   fontsize_row = 9, fontsize_col = 9, border_color = NA, cellwidth = 8, cellheight = 8,
                   clustering_distance_rows = 'euclidean', clustering_method = 'complete', clustering_callback = callback, 
                   main = 'HCA Census', filename = 'correlation/SpearmanCorr.allGenes.supergroup.pdf')
pheatmap::pheatmap(corMat, color = colorRampPalette(c("steelblue", 'orange', "red2", 'red4'))(n = 4),
                   breaks = c(0.8, 0.85, 0.9, 0.95, 1.0),
                   fontsize_row = 9, fontsize_col = 9, border_color = NA, cellwidth = 8, cellheight = 8,
                   clustering_distance_rows = 'euclidean', clustering_method = 'complete', clustering_callback = callback,
                   main = 'HCA Census', filename = 'correlation/SpearmanCorr.allGenes.supergroup.categorized.pdf')


###
ica_bm.combined.supergroup@meta.data$supergroup <- mapvalues(ica_bm.combined.supergroup@meta.data$supergroup,
                                                             from = unique(ica_bm.combined.supergroup@meta.data$supergroup),
                                                             to = c('Naive T 1', 'Monocyte 1', 'Pro-B 1', 'NK', 'B cell', 'HSC~MPP',
                                                                    'cDC', 'Naive T 2', 'Pro-B 2', 'CD8', 'B precursor', 'Erythroblast 2', 
                                                                    'Stromal', 'Granulocyte progenitors', 'pDC', 'Erythroblast 1', 'T small',
                                                                    'Plasma cell', 'Monocyte 2', 'Platelet'))

ica_bm.combined.supergroup@meta.data$supergroup <- factor(ica_bm.combined.supergroup@meta.data$supergroup,
                                                          levels = c('HSC~MPP', 'Erythroblast 1', 'Erythroblast 2', 'Platelet', 'Granulocyte progenitors',
                                                                     'Monocyte 1', 'Monocyte 2', 'cDC', 'pDC', 'B precursor', 'Pro-B 1', 'Pro-B 2', 'B cell', 
                                                                     'Plasma cell', 'T small', 'Naive T 1', 'Naive T 2', 'CD8', 'NK', 'Stromal'))
Idents(ica_bm.combined.supergroup) <- 'supergroup'

ica_bm.combined.supergroup.markers <- FindAllMarkers(object = ica_bm.combined.supergroup, min.pct = 0.25, logfc.threshold = 0.25, only.pos = T)
ica_bm.combined.supergroup.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
write.table(ica_bm.combined.supergroup.markers, 'degs/findAllMarkers.pct25.thre0.25.onlyPos.txt', sep = '\t', quote = F, col.names = T, row.names = F)
top10 <- ica_bm.combined.supergroup.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

ica_bm.combined.supergroup@meta.data$use_anno <- ica_bm.combined.supergroup@active.ident
writeLable <- as.matrix(ica_bm.combined.supergroup@meta.data)
writeLable <- data.frame(Barcode = rownames(writeLable), writeLable, check.names = F)
write.table(writeLable, 'tmp7.label.txt', sep = '\t', quote = F, row.names = F, col.names = T)
save.image('tmp7.Rdata')

###
#ica_bm.combined.supergroup <- ScaleData(object = ica_bm.combined.supergroup, vars.to.regress = c('Patient', 'nCount_RNA')) ### takes too long time
DefaultAssay(ica_bm.combined.supergroup) <- 'integrated'
saveRDS(ica_bm.combined.supergroup, 'tmp/ica_bm.combined.supergroup.Rds')

DoHeatmap(object = ica_bm.combined.supergroup, features = unique(top10$gene), angle = 90, size = 3, raster = T, draw.lines = F)
ggsave('degs/degs.DEgenes.top10.pdf', units = 'cm', width = 50, height = 40, limitsize = F)

DotPlot(object = ica_bm.combined.supergroup, features = unique(top10$gene), dot.scale = 3)
ggsave('degs/dotplt.DEgenes.top10.pdf', units = 'cm', width = 50, height = 40, limitsize = F)
# FCGR3A for NK and CD16+ monocyte; NCAM1 for NK; B3GAT1 for effecter; CD69 for memory;
DotPlot(object = ica_bm.combined.supergroup, features = c('CD34', 'FCGR3A', 'NCAM1', 'B3GAT1', 'CD69', 'CD8A', 'CD4'), dot.scale = 6, assay='RNA') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('degs/dotplt.markers1.pdf', units = 'cm', width = 20, height = 30)

