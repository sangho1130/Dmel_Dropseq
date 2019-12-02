library('GSEABase')
library('GSVA')

data <- read.delim('hca.pseudo_exprs.supergroup.dmelSymbols.txt', header = T, row.names = 1, check.names = F)
data$Stromal <- NULL

pData <- AnnotatedDataFrame(data=data.frame(row.names = colnames(data), cell = colnames(data)))
dataSet <- ExpressionSet(as.matrix(data), phenoData = pData, annotation = 'Symbol')

##
geneSet <- read.delim('testDEGs.anno_simple.all.txt', header = F, row.names = 1)
geneSet <- geneSet[-1,] # remove PSC
celltype <- data.frame(row.names = rownames(geneSet), rownames(geneSet))
colnames(celltype) <- 'classification'
geneSet$V2 <- NULL
geneSet <- data.frame(t(geneSet))

###
geneSet_unbiased <- data.frame(matrix(ncol = ncol(geneSet), nrow = 30))
colnames(geneSet_unbiased) <- colnames(geneSet)
for (targetcelltype in colnames(geneSet)) {
  foundGenes <- intersect(geneSet[,targetcelltype], rownames(data))
  
  print (targetcelltype)
  print (length(foundGenes))
  geneSet_unbiased[, targetcelltype] <- foundGenes[1:30]
}
dim(geneSet_unbiased)


identical(rownames(celltype), colnames(geneSet))
phenocelltype <- AnnotatedDataFrame(data=celltype)
signatureSet <- ExpressionSet(as.matrix(geneSet_unbiased), phenoData=phenocelltype)
exprs(signatureSet)

ph <- GeneSet(exprs(signatureSet)[,1][exprs(signatureSet)[,1]!=''], setName = as.character(signatureSet$classification[1]))
pm <- GeneSet(exprs(signatureSet)[,2][exprs(signatureSet)[,2]!=''], setName = as.character(signatureSet$classification[2]))
lm <- GeneSet(exprs(signatureSet)[,3][exprs(signatureSet)[,3]!=''], setName = as.character(signatureSet$classification[3]))
cc <- GeneSet(exprs(signatureSet)[,4][exprs(signatureSet)[,4]!=''], setName = as.character(signatureSet$classification[4]))
gst <- GeneSet(exprs(signatureSet)[,5][exprs(signatureSet)[,5]!=''], setName = as.character(signatureSet$classification[5]))
adipo <- GeneSet(exprs(signatureSet)[,6][exprs(signatureSet)[,6]!=''], setName = as.character(signatureSet$classification[6]))
mySignatures <- GeneSetCollection(object = ph, pm, lm, cc, gst, adipo)


### GSVA ###
es.gsva <- gsva(expr = dataSet, gset.idx.list = mySignatures, method = 'gsva', verbose=T, kcdf = 'Gaussian', parallel.sz = 1)

gsvaScores <- as.data.frame(t(exprs(es.gsva)))
gsvaScores <- cbind(anno_simple = es.gsva$cell, gsvaScores)

callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

pheatmap(as.matrix(t(scoredata)), cluster_rows = F, cluster_cols = T, scale = 'none', treeheight_row = 15, treeheight_col = 15, cellwidth = 10, cellheight = 10, fontsize_row = 9, fontsize_col = 9,  
         color = colorRampPalette(c("#3a5fcd", '#fefebd', "#ee0000"))(n = 100), border_color = NA, filename = 'testDEGs.anno_simple.all.pdf')
