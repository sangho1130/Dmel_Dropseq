
#library(Seurat)
#lymphgland <- readRDS('../../Normal_LG/rdata/lymphgland.Rds')
#levels(Idents(lymphgland))

#earlycells <- subset(lymphgland, idents = c('PH 1', 'PH 2', 'PSC'))
#earlycells@meta.data <- droplevels(earlycells@meta.data)
#saveRDS(earlycells@meta.data, 'earlycells_labels.Rds')

#exprs <- as.matrix(GetAssayData(earlycells, slot = 'counts', assay = 'RNA'))
#saveRDS(exprs, 'earlycells_counts.Rds')
#dim(exprs)

###
library('SCENIC'); packageVersion('SCENIC')
library('RcisTarget'); packageVersion('RcisTarget')

expr <- readRDS('earlycells_counts.Rds')
label <- readRDS('earlycells_labels.Rds')
label <- label[, c(2, 3, 5, 6, 9)]
head(label); nrow(label) # 348

org <- 'dmel'
dbDir <- 'databases' # RcisTarget databases location
myDatasetTitle='PSC_PH1-2'
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, datasetTitle=myDatasetTitle, nCores=20)

label$nCount_RNA <- as.numeric(label$nCount_RNA)
label$nFeature_RNA <- as.numeric(label$nFeature_RNA)
label$percent.mt <- as.numeric(label$percent.mt)

dir.create("int")
saveRDS(label, file="int/cellInfo.Rds")

scenicOptions@inputDatasetInfo$label <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

### Co-expression network ###
nCellsPerGene <- apply(expr, 1, function(x) sum(x>0))
nCountsPerGene <- apply(expr, 1, sum)
summary(nCellsPerGene)
summary(nCountsPerGene)
#max(expr) # 1001
#sum(expr>0) / sum(expr==0) # 0.08391485

# 348*0.025 = 8.7
# 348*0.05 = 17.4
# First filter:
minReads <- 3*0.03*ncol(expr);minReads # 31.32
genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > minReads)]
length(genesLeft_minReads) # 4209 genes
# Second filter:
minSamples <- ncol(expr)*0.03;minSamples #10.44
nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > minSamples)]
length(genesLeft_minCells) # 4125 genes

motifRankings <- importRankings(getDatabases(scenicOptions)[[1]]) # either one, they should have the same genes
genesInDatabase <- colnames(getRanking(motifRankings))

genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]
length(genesLeft_minCells_inDatabases) # 3844
interestingGenes <- c("Dl", 'Stat92E', 'kn', 'Hand') ### Stat92 and kn
interestingGenes <- c('Atg5', 'CG13197', 'dome', 'E(spl)m4-BFM', 'Eaf6', 'jigr1', 'KrT95D', 'melt', 'Snoo', 'DOR', 'LIMK1', 'dpp', 'Mmp2', 'Adar', 'bnl', 'CR45102', 'Dl', 'DnaJ-H', 'pgant6', 'qsm', 'PGRP-LC', 'vimar', 'CG6040') ### Stat92 targets
interestingGenes[which(!interestingGenes %in% genesLeft_minCells_inDatabases)]


genesKept <- genesLeft_minCells_inDatabases
saveRDS(genesKept, file=getIntName(scenicOptions, "genesKept"))
expr_filtered <- expr[genesKept, ]
expr_filtered <- as.matrix(expr_filtered)

### Correlation ###
corrMat <- cor(t(expr_filtered), method="spearman")
saveRDS(corrMat, file=getIntName(scenicOptions, "corrMat"))
save.image('tmp1.Rdata')


### GENIE3 ###
#load('tmp1.Rdata')
expr_filtered <- log2(expr_filtered + 1)
runGenie3(as.matrix(expr_filtered), scenicOptions)
save.image('tmp2.Rdata')

### runSCENIC ###
#load('tmp2.Rdata')
expr <- as.matrix(expr)
expr <- log2(expr+1)
dim(expr) # 
#scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$seed <- 2004093
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, expr)
save.image('tmp3.Rdata')


#load('tmp3.Rdata')
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=c(15, 20, 25, 30, 35, 40, 45, 50), perpl=c(10, 20, 30, 40, 50))
pdf('int/tSNE_compare.pdf')
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=T, varName="new_subclustering", cex=.5) # saving in pdf files
dev.off()
pdf('int/tSNE_compare.oneshot.pdf', width = 30, height = 22)
par(mfcol=c(5, 8))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=T, varName="new_subclustering", cex=1)
dev.off()

# Using only "high-confidence" regulons (normally similar)
fileNames.highconf <- tsneAUC(scenicOptions, aucType="AUC", nPcs=c(15, 20, 25, 30, 35, 40, 45, 50), perpl=c(10, 20, 30, 40, 50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
pdf('int/tSNE_compare.oHC.pdf')
fileNames.highconf <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames.highconf, scenicOptions, showLegend=T, varName="new_subclustering", cex=.5)
dev.off()
pdf('int/tSNE_compare.oHC.oneshot.pdf', width = 30, height = 22)
par(mfcol=c(5, 8))
fileNames.highconf <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames.highconf, scenicOptions, showLegend=T, varName="new_subclustering", cex=1)
dev.off()
save.image('tmp4.post_tSNE.Rdata')

### 20200409 trial5 ###
scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 45
scenicOptions@settings$defaultTsne$perpl <- 10
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

runSCENIC_3_scoreCells(scenicOptions, expr) 
save.image('tmp4.post_tSNE.scoreCells.Rdata')


######################
### Local computer ###
######################
library('SCENIC'); packageVersion('SCENIC')
library('RcisTarget'); packageVersion('RcisTarget')

load('tmp4.post_tSNE.scoreCells.Rdata')
scenicOptions@settings$nCores <- 3
exprMat <- expr
aucellApp <- plotTsne_AUCellApp(scenicOptions, expr)
savedSelections <- shiny::runApp(aucellApp)
savedSelections <- shiny::runApp(aucellApp)

# Save the modified thresholds:
#scenicOptions@settings$seed <- 123
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
#save.image('tmp4.post_tSNE.scoreCells.newthresholds2.Rdata')

#load('tmp4.post_tSNE.scoreCells.newthresholds.Rdata')
runSCENIC_4_aucell_binarize(scenicOptions)
#save.image('tmp5.threshold_updates.ver1.Rdata')
#load('tmp5.threshold_updates.ver1.Rdata')

