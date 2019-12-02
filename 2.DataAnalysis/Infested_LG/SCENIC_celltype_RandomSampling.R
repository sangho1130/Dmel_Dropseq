library('SCENIC')#; packageVersion('SCENIC')
library('RcisTarget')#; packageVersion('RcisTarget')
library('plyr')

expr <- read.delim('merged.3tps.expr.allGenes.txt', header = T, sep = '\t', check.names = F, row.names = 1)
label <- read.delim('merged.3tps.clustering_labels.txt', header = T, sep = '\t', check.names = F, row.names = 1)

colnames(label) <- c('Timepoint', 'Celltype', 'Subcluster')
label <- subset(label, Celltype != 'DV' & Celltype != 'Neurons' & Celltype != 'RG'); label <- droplevels(label)

label$Celltype <- factor(label$Celltype, levels = c('PSC', 'PH', 'PM', 'LM', 'CC', 'GST-rich', 'Adipohemocyte'))
nrow(label) # 19332 cells

expr <- expr[, rownames(label)]
identical(colnames(expr), rownames(label)); dim(expr)

summary(label$Celltype)
summary(label$Subcluster)

for (i in c(1:33)){
  if (i < 10){
      tmpDir <- paste0(c('trial_00', i), collapse = '')
  } else {
      tmpDir <- paste0(c('trial_0', i), collapse = '')
  }
  dir.create(tmpDir)
  setwd(tmpDir)

  rand_barcodes <- c()
  for (celltype in levels(label$Subcluster)){
    celltmp <- subset(label, Subcluster == celltype)
    targetRows <- sort(sample(x = 1:nrow(celltmp), size = 42))
    targetBcs <- rownames(celltmp[targetRows, ])
    rand_barcodes <- append(rand_barcodes, targetBcs)
  }
  rand_expr <- as.matrix(expr[,rand_barcodes])
  rand_label <- label[rand_barcodes, ]
  #identical(rownames(rand_label), colnames(rand_expr))
  
  rand_label$Subcluster <- NULL
  rand_label$Celltype <- factor(rand_label$Celltype, levels = c('PSC', 'PH', 'PM', 'LM', 'CC', 'GST-rich', 'Adipohemocyte'))
  colVars <- list(CellType = c("PSC"="#f15fa6", "PH"="#207eb3", "PM"="#a80d0c", "LM"="#f0a142", "CC"="#25a9b0", "GST-rich"="#a4a4a4", "Adipohemocyte"="#1a1a1a"))
  
  org <- 'dmel'
  dbDir <- '../../databases' 
  myDatasetTitle='Major'
  scenicOptions <- initializeScenic(org=org, dbDir=dbDir, datasetTitle=myDatasetTitle, nCores=30)
  
  dir.create("int")
  saveRDS(label, file="int/cellInfo.Rds")
  saveRDS(colVars, file="int/colVars.Rds")
  
  scenicOptions@inputDatasetInfo$label <- "int/cellInfo.Rds"
  scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
  saveRDS(scenicOptions, file="int/scenicOptions.Rds")
  
  ### 
  genesKept <- geneFiltering(rand_expr, scenicOptions=scenicOptions,
                             minCountsPerGene=3*.01*ncol(rand_expr),
                             minSamples=ncol(rand_expr)*.01)

  saveRDS(genesKept, file=getIntName(scenicOptions, "genesKept"))
  expr_filtered <- rand_expr[genesKept, ]
  expr_filtered <- as.matrix(expr_filtered)
  
  ### Correlation ###
  runCorrelation(expr_filtered, scenicOptions)
  expr_filtered <- log2(expr_filtered + 1)
  ### GENIE3 ###
  runGenie3(as.matrix(expr_filtered), scenicOptions)
  
  scenicOptions@settings$verbose <- TRUE
  scenicOptions@settings$seed <- i
  saveRDS(scenicOptions, file="int/scenicOptions.Rds")
  
  runSCENIC_1_coexNetwork2modules(scenicOptions)
  runSCENIC_2_createRegulons(scenicOptions)
  runSCENIC_3_scoreCells(scenicOptions, log2(rand_expr+1))
  runSCENIC_4_aucell_binarize(scenicOptions)

  setwd('../')
}
